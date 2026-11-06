import os
import glob
import subprocess
import pickle
import time
import numpy as np
import cv2
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import astropy.units as u
from scipy.ndimage import binary_dilation
import vtk
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# ============================================================
# DisPerSE 実行 + FITS 切り出し 保存用スクリプト(最終版)
#
# FilFinder版パイプラインと対になる、DisPerSEベースの暫定ラベル生成
# ============================================================

# ------------------------------------------------------------
# 設定
# ------------------------------------------------------------
root_dir = "/mnt/sp4t/fujiwara.r"
fits_name = "HGBS_aquilaM2_hires_column_density_map"
fits_path = os.path.join(root_dir, fits_name + ".fits")

DISPERSE_BIN_DIR = "/home1/fujiwara.r/PycharmProjects/DisPerSE/build/src/mse"

def bin_path(name):
    return os.path.join(DISPERSE_BIN_DIR, name) if DISPERSE_BIN_DIR else name

# 探索の結果、目視確認で妥当と判断した値
NSIG = 1

# マスク幅の方針:
#   "skeleton_only": 太らせない(線が含まれる領域そのまま)
#   "physical"      : 典型物理幅WIDTH_PCを仮定して太らせる
WIDTH_MODE = "skeleton_only"
WIDTH_PC = 0.1
distance_pc = 436  # FilFinderスクリプトと同じ(Aquila)

# タイル切り出し設定(FilFinder版と共通)
size = 256
step_rate = 0.0625  # ステップ = 16px
tiles_per_dir = 1000

# ------------------------------------------------------------
# 出力ディレクトリ
# ------------------------------------------------------------
disperse_root = os.path.join(root_dir, fits_name, "disperse")
raw_dir  = os.path.join(disperse_root, "raw")   # mse/skelconv 生出力
full_dir = os.path.join(disperse_root, "full")  # 正規化画像・マスク全面画像
vis_dir  = os.path.join(disperse_root, "vis")   # オーバーレイ可視化
in_dir   = os.path.join(disperse_root, "in")    # 入力パッチ(FilFinder版のin/と同じ形式)
out_dir  = os.path.join(disperse_root, "out")   # ラベルパッチ

for d in [raw_dir, full_dir, vis_dir, in_dir, out_dir]:
    os.makedirs(d, exist_ok=True)


def check_make_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)


def save_pkl(fpath, val):
    with open(fpath, 'wb') as f:
        pickle.dump(val, f)


# ============================================================
# ステップ1: FITS読み込み・正規化(FilFinder版のdata_bak_mlと同じロジック)
# ============================================================
hdu = fits.open(fits_path)[0]
w = WCS(hdu)

img_raw = hdu.data.astype(np.float64).copy()
img_raw[~np.isfinite(img_raw)] = np.nan
finite_mask = np.isfinite(img_raw)

img_log = np.where(finite_mask & (img_raw > 0), np.log10(img_raw), np.nan)
des_log = pd.DataFrame(pd.Series(img_log[np.isfinite(img_log)].ravel()).describe(
    percentiles=[0.005, 0.995]))
data_bak_ml = img_log.copy()
data_bak_ml -= des_log[0]["0.5%"]
data_bak_ml /= ((des_log[0]["99.5%"] - des_log[0]["0.5%"]) / 255)
data_bak_ml = np.clip(data_bak_ml, 0, 255)
data_bak_ml[~finite_mask] = np.nan

input_u8 = np.where(finite_mask, data_bak_ml, 0).astype("uint8")
cv2.imwrite(os.path.join(full_dir, fits_name + "_input_full_ml.png"), input_u8)
print("正規化済み入力画像を保存しました。")


# ============================================================
# ステップ2: mse 実行(nsig=1)
# ============================================================
param_suffix = f"_s{NSIG}"
expected_ndskl_pattern = os.path.join(raw_dir, f"*{param_suffix}.up.NDskl")
ndskl_candidates = sorted(glob.glob(expected_ndskl_pattern))

if not ndskl_candidates:
    cmd = [bin_path("mse"), fits_path, "-upSkl", "-nsig", str(NSIG), "-outDir", raw_dir]
    print("RUN:", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    print(res.stdout[-2000:])
    if res.returncode != 0:
        print(res.stderr[-2000:])
        raise RuntimeError("mse failed")
    ndskl_candidates = sorted(glob.glob(expected_ndskl_pattern))

if not ndskl_candidates:
    raise FileNotFoundError(
        f"'{expected_ndskl_pattern}' に一致するファイルがありません。"
        f"raw_dir の中身: {os.listdir(raw_dir)}"
    )
ndskl = ndskl_candidates[-1]
print("使用する NDskl:", ndskl)


# ============================================================
# ステップ3: skelconv → VTK
# ============================================================
expected_vtk_pattern = ndskl + "*.vtk"
existing_vtk = sorted(glob.glob(expected_vtk_pattern))

if not existing_vtk:
    cmd = [bin_path("skelconv"), ndskl, "-to", "vtk", "-outDir", raw_dir]
    print("RUN:", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    print(res.stdout[-2000:])
    if res.returncode != 0:
        print(res.stderr[-2000:])
        raise RuntimeError("skelconv failed")
    existing_vtk = sorted(glob.glob(expected_vtk_pattern))

if not existing_vtk:
    raise FileNotFoundError(f"VTKファイルが見つかりません: {expected_vtk_pattern}")
vtk_path = existing_vtk[0]
print("VTK:", vtk_path)


# ============================================================
# ステップ4: VTK読み込み・ラスタ化
# ============================================================
def read_vtk_skeleton(path):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    poly = reader.GetOutput()
    filaments = []
    lines = poly.GetLines()
    lines.InitTraversal()
    idList = vtk.vtkIdList()
    while lines.GetNextCell(idList):
        pts = [poly.GetPoint(idList.GetId(i))[:2] for i in range(idList.GetNumberOfIds())]
        filaments.append(np.array(pts))
    return filaments

filaments = read_vtk_skeleton(vtk_path)
print("検出フィラメント数:", len(filaments))

shape = hdu.data.shape
skeleton = np.zeros(shape, dtype=np.uint8)
for pts in filaments:
    xs, ys = pts[:, 0], pts[:, 1]
    for k in range(len(pts) - 1):
        p1 = (int(round(xs[k])), int(round(ys[k])))
        p2 = (int(round(xs[k + 1])), int(round(ys[k + 1])))
        cv2.line(skeleton, p1, p2, 1, 1)

cv2.imwrite(os.path.join(full_dir, fits_name + "_skeleton_01.png"), skeleton)
cv2.imwrite(os.path.join(full_dir, fits_name + "_skeleton_label01.png"), skeleton)
print("skeleton画素数:", skeleton.sum(), "/ 割合: {:.4f}%".format(100 * skeleton.sum() / skeleton.size))


# ============================================================
# ステップ5: マスク化(create_mask相当)
# ============================================================
if WIDTH_MODE == "skeleton_only":
    mask = skeleton.copy()
    print("マスクモード: skeleton_only(太らせなし)")

elif WIDTH_MODE == "physical":
    pixel_scale_deg = proj_plane_pixel_scales(w)[:2].mean()
    pixel_scale_rad = (pixel_scale_deg * u.deg).to(u.rad).value
    pc_per_pixel = distance_pc * pixel_scale_rad
    width_px = WIDTH_PC / pc_per_pixel
    radius = max(1, int(round(width_px / 2)))

    yy, xx = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    struct = (yy**2 + xx**2 <= radius**2)
    mask = binary_dilation(skeleton.astype(bool), structure=struct).astype(np.uint8)
    print(f"マスクモード: physical ({WIDTH_PC}pc → {width_px:.2f}px → radius={radius}px)")

else:
    raise ValueError("WIDTH_MODE は 'skeleton_only' か 'physical' を指定してください")

cv2.imwrite(os.path.join(full_dir, fits_name + "_mask_01.png"), mask)
cv2.imwrite(os.path.join(full_dir, fits_name + "_mask_label01.png"), mask)


# ============================================================
# ステップ6: 可視化(matplotlibベース、WCS回転なしのピクセル直接オーバーレイ)
# ============================================================
def save_overlay_png(background_u8, overlay_mask, out_path, color=(1, 0, 0), alpha=0.9):
    fig, ax = plt.subplots(figsize=(16, 16))
    ax.imshow(background_u8, cmap="gray", origin="upper", vmin=0, vmax=255)
    overlay_rgba = np.zeros((*overlay_mask.shape, 4), dtype=float)
    overlay_rgba[..., 0] = color[0]
    overlay_rgba[..., 1] = color[1]
    overlay_rgba[..., 2] = color[2]
    overlay_rgba[..., 3] = overlay_mask.astype(float) * alpha
    ax.imshow(overlay_rgba, origin="upper")
    ax.axis("off")
    fig.savefig(out_path, dpi=150, bbox_inches="tight", pad_inches=0)
    plt.close(fig)

save_overlay_png(
    input_u8, skeleton.astype(bool),
    os.path.join(vis_dir, f"{fits_name}_disperse_skeleton_overlay_nsig{NSIG}.png"),
)
save_overlay_png(
    input_u8, mask.astype(bool),
    os.path.join(vis_dir, f"{fits_name}_disperse_mask_overlay_nsig{NSIG}.png"),
)
print("オーバーレイ画像を保存しました。")


# ============================================================
# ステップ7: タイル切り出し(FilFinder版のcrop_data()と同じ仕様)
#   in/  : 入力パッチ(data_bak_ml、0-255正規化済み)
#   out/ : ラベルパッチ(mask、0/1)
#   NaN判定は data_bak_ml を基準にする(観測領域外を除外)
# ============================================================
def crop_data(nda, nda_check_nan, base_dir, suf):
    w_, h_ = nda.shape[1::-1]
    print(f'canvas w,h = {w_}, {h_} ({suf})')
    step = int(step_rate * size)

    tile_count = 0
    y = 0
    while (y + size <= h_):
        x = 0
        while (x + size <= w_):
            cropped = nda[y:y + size, x:x + size]
            cropped_check = nda_check_nan[y:y + size, x:x + size]

            o_name = fits_name + "_" + str(x).zfill(5) + "-" + str(y).zfill(5)
            x += step

            if np.isnan(cropped_check).any():
                continue

            sub_dir = os.path.join(base_dir, "tiles_{:03d}".format(tile_count // tiles_per_dir))
            check_make_dir(sub_dir)

            o_path = os.path.join(sub_dir, o_name + "_" + suf + ".pkl")
            save_pkl(o_path, cropped)
            o_path = os.path.join(sub_dir, o_name + "_" + suf + ".png")
            cv2.imwrite(o_path, cropped.astype("int"))

            tile_count += 1
        y += step

    print(f'Total tiles saved ({suf}): {tile_count}')

crop_data(data_bak_ml, data_bak_ml, in_dir, "in")     # 入力パッチ
crop_data(mask, data_bak_ml, out_dir, "out")          # ラベルパッチ(mask基準)

print("DONE")
