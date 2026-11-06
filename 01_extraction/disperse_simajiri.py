#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DisPerSE 実行 + FITS 処理・可視化スクリプト (最終統合版)
- Av=4未満をNaNにしてメモリ削減
- タイル切り出し(in, out)を完全廃止
- 端の赤枠ノイズを cv2.erode で完全除去
- 背景に全体のグレースケール画像を採用
- 可視化の向きを FilFinder と同じ origin="lower" に統一

[修正 2026-07-17]
- ステップ6に save_science_figure() を追加。
  inferno カラーマップ＋カラーバー＋タイトル(天体名/距離/beam/フィラメント数)＋
  x/y軸(pixel)＋5pcスケールバー付きの「科学図」を、
  (1) skeleton(赤い細線)を重ねたもの
  (2) mask(赤い太い半透明領域)を重ねたもの
  の2枚として保存する。
"""

import os
import glob
import subprocess
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
import matplotlib.colors

# ============================================================
# 設定 (CONFIG)
# ============================================================
root_dir = "/mnt/sp4t/fujiwara.r"
# 対象の天体FITS（必要に応じて変更してください）
fits_name = "HGBS_aquilaM2_hires_column_density_map"
fits_path = os.path.join(root_dir, fits_name + ".fits")

# science_figure のタイトル・スケールバー計算に使う情報
display_name = "aquilaM2"   # タイトルに表示する天体名（必要に応じて変更）
beam_arcsec = 18.2            # ビームサイズ["]（hires列密度マップの典型値。実測値がわかれば修正）

# DisPerSE の mse コマンドがあるバイナリディレクトリ
DISPERSE_BIN_DIR = "/home1/fujiwara.r/PycharmProjects/DisPerSE/build/src/mse"

def bin_path(name):
    return os.path.join(DISPERSE_BIN_DIR, name) if DISPERSE_BIN_DIR else namex

# DisPerSE 固有パラメータ
NSIG = 1

# マスク幅の方針
WIDTH_MODE = "skeleton_only"
WIDTH_PC = 0.1
distance_pc = 440.0  # TauN3/Orion B 等に合わせて調整してください

# Av = 4 相当の柱密度しきい値 (1.88e21 x 2 = 3.76e21 H2 cm-2)
AV4_THRESH = 7.52e+21

# ------------------------------------------------------------
# 出力ディレクトリの作成
# ------------------------------------------------------------
disperse_root = os.path.join(root_dir, fits_name, "disperse")
raw_dir  = os.path.join(disperse_root, "raw")   # mse/skelconv 生出力
full_dir = os.path.join(root_dir, fits_name, "full")  # 全面・FITS出力フォルダ
vis_dir  = os.path.join(root_dir, fits_name, "vis")   # 可視化画像フォルダ

for d in [raw_dir, full_dir, vis_dir]:
    os.makedirs(d, exist_ok=True)


# ============================================================
# ステップ1: FITS読み込み・しきい値(Av < 4)処理・正規化
# ============================================================
# DisPerSE用の一時的な「Av4マスク済みFITSファイル」のパス
masked_fits_path = os.path.join(raw_dir, fits_name + "_av4_masked.fits")

with fits.open(fits_path, memmap=True) as hdul:
    hdu_header = hdul[0].header.copy()
    # 2次元に平坦化
    img_raw = np.squeeze(np.asarray(hdul[0].data, dtype=np.float64))

# WCSの取得
w = WCS(hdu_header)
if w.naxis > 2:
    w = w.celestial

# --- 可視化用：元のデータ（Av < 4 の背景も含んだ全体）の正規化画像を作成 ---
finite_mask_all = np.isfinite(img_raw)
img_log_all = np.where(finite_mask_all & (img_raw > 0), np.log10(img_raw), np.nan)
des_log_all = pd.DataFrame(pd.Series(img_log_all[np.isfinite(img_log_all)].ravel()).describe(
    percentiles=[0.005, 0.995]))

data_bak_ml = img_log_all.copy()
data_bak_ml -= des_log_all[0]["0.5%"]
data_bak_ml /= ((des_log_all[0]["99.5%"] - des_log_all[0]["0.5%"]) / 255)
data_bak_ml = np.clip(data_bak_ml, 0, 255)
data_bak_ml[~finite_mask_all] = np.nan

# この input_u8 が、背景となる「全体グレースケール画像」になります
input_u8 = np.where(finite_mask_all, data_bak_ml, 0).astype("uint8")
cv2.imwrite(os.path.join(full_dir, fits_name + "_input_full_ml.png"), input_u8)
print("全体グレースケール（背景用）入力を保存しました。")


# --- DisPerSE用：Av = 4 未満を NaN にしたマスク済みFITSの作成 ---
img_raw_masked = img_raw.copy()
print(f"DisPerSE用：Av = 4 閾値 ({AV4_THRESH}) 未満の領域を NaN に置き換えています...")
img_raw_masked[img_raw_masked < AV4_THRESH] = np.nan

masked_hdu = fits.PrimaryHDU(data=img_raw_masked, header=hdu_header)
masked_hdu.writeto(masked_fits_path, overwrite=True)
print(f"DisPerSE用のマスク済みFITSを保存しました: {masked_fits_path}")


# ============================================================
# ステップ2: mse 実行 (マスク済みFITSを入力にする)
# ============================================================
param_suffix = f"_s{NSIG}"
expected_ndskl_pattern = os.path.join(raw_dir, f"*{param_suffix}.up.NDskl")
ndskl_candidates = sorted(glob.glob(expected_ndskl_pattern))

if not ndskl_candidates:
    cmd = [bin_path("mse"), masked_fits_path, "-upSkl", "-nsig", str(NSIG), "-outDir", raw_dir]
    print("RUN:", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    print(res.stdout[-2000:])
    if res.returncode != 0:
        print(res.stderr[-2000:])
        raise RuntimeError("mse failed")
    ndskl_candidates = sorted(glob.glob(expected_ndskl_pattern))

if not ndskl_candidates:
    raise FileNotFoundError(f"'{expected_ndskl_pattern}' に一致するファイルがありません。")
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
# ステップ4: VTK読み込み・ラスタ化 (境界のノイズ除去＆線幅を太く)
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

shape = img_raw.shape
skeleton = np.zeros(shape, dtype=np.uint8)

# 白いフィラメントの上に重ねる赤線を太く描くため、thicknessを 3 に設定
for pts in filaments:
    xs, ys = pts[:, 0], pts[:, 1]
    for k in range(len(pts) - 1):
        p1 = (int(round(xs[k])), int(round(ys[k])))
        p2 = (int(round(xs[k + 1])), int(round(ys[k + 1])))
        cv2.line(skeleton, p1, p2, 1, thickness=3)

# 境界の偽フィラメント（赤線ノイズ）を消すため、Av >= 4 マスクを作成
valid_mask = np.isfinite(img_raw_masked)
kernel = np.ones((5, 5), np.uint8)
valid_mask_eroded = cv2.erode(valid_mask.astype(np.uint8), kernel, iterations=2).astype(bool)

# マスク外（Av < 4）および境界ギリギリのスケルトンを強制除去
skeleton[~valid_mask_eroded] = 0

cv2.imwrite(os.path.join(full_dir, fits_name + "_skeleton_01.png"), skeleton)
cv2.imwrite(os.path.join(full_dir, fits_name + "_skeleton_label01.png"), skeleton)
print("skeleton画素数（境界ノイズ除去後）:", skeleton.sum())


# ============================================================
# ステップ5: マスク化
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
    mask[~valid_mask_eroded] = 0
    print(f"マスクモード: physical ({WIDTH_PC}pc → {width_px:.2f}px → radius={radius}px)")

else:
    raise ValueError("WIDTH_MODE は 'skeleton_only' か 'physical' を指定してください")

cv2.imwrite(os.path.join(full_dir, fits_name + "_mask_01.png"), mask)
cv2.imwrite(os.path.join(full_dir, fits_name + "_mask_label01.png"), mask)


# ============================================================
# ステップ6: 可視化 (FilFinderと同一の origin="lower" で全体背景重ね合わせ)
# ============================================================
def save_overlay_png(background_u8, overlay_mask, out_path, color=(1, 0, 0), alpha=1.0):
    fig, ax = plt.subplots(figsize=(16, 16))
    
    # FilFinderと同じ向きにするため origin="lower" を明示
    ax.imshow(background_u8, cmap="gray", origin="lower", vmin=0, vmax=255)
    
    overlay_rgba = np.zeros((*overlay_mask.shape, 4), dtype=float)
    overlay_rgba[..., 0] = color[0]  # R (赤)
    overlay_rgba[..., 1] = color[1]  # G
    overlay_rgba[..., 2] = color[2]  # B
    overlay_rgba[..., 3] = overlay_mask.astype(float) * alpha
    
    ax.imshow(overlay_rgba, origin="lower")
    ax.axis("off")
    fig.savefig(out_path, dpi=150, bbox_inches="tight", pad_inches=0)
    plt.close(fig)

# 全体グレースケール画像 (input_u8) の上に、太いスケルトンを重ね合わせます
save_overlay_png(
    input_u8, skeleton.astype(bool),
    os.path.join(vis_dir, f"{fits_name}_disperse_skeleton_overlay_nsig{NSIG}.png"),
    color=(1, 0, 0),
    alpha=1.0  # 赤線をはっきり見せるために 1.0 に設定
)

save_overlay_png(
    input_u8, mask.astype(bool),
    os.path.join(vis_dir, f"{fits_name}_disperse_mask_overlay_nsig{NSIG}.png"),
    color=(1, 0, 0),
    alpha=0.6
)
print("グレースケール背景の上に、FilFinderと同じ向きで赤線を重ねた画像を保存しました。")


# ------------------------------------------------------------
# 論文・FilFinder出力と同じ体裁の科学図（inferno + colorbar + タイトル + pixel軸 +
# 赤いフィラメント線/マスク + 右下の横棒スケールバー）を追加で保存する。
# ------------------------------------------------------------
def save_science_figure(data2d, overlay_mask, out_path, wcs, distance_pc,
                         display_name, beam_arcsec, n_filaments,
                         scalebar_pc=5.0, vmax=None, cmap="inferno",
                         filament_color="red", filament_alpha=1.0):
    """
    data2d        : 物理値(線形, cm^-2)の柱密度2次元配列。NaNは黒表示になる。
    overlay_mask  : data2dと同形状の bool配列。Trueの画素に色を重ねて描く。
                    skeleton(1px芯線)を渡せば細線オーバーレイ、mask(太らせ済み)を
                    渡せば太い領域オーバーレイになる。
    filament_color: matplotlibが解釈できる色名 (例: "red", "cyan")
    filament_alpha: 重ね描きの不透明度 (skeletonは1.0推奨、maskは0.6程度推奨)
    n_filaments   : タイトルに表示するフィラメント数。
                    注意: read_vtk_skeleton()で数えたVTKセグメント数は、
                    arc-assembly（隣接アークの接続）を経ていないため
                    「科学的に意味のあるフィラメント本数」より過大になりやすい。
                    正式な本数があれば置き換えること。
    """
    ny, nx = data2d.shape

    if vmax is None:
        vmax = np.nanpercentile(data2d[np.isfinite(data2d)], 99.5)

    fig, ax = plt.subplots(figsize=(11, 8.5))

    data_masked = np.ma.masked_invalid(data2d)
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad("black")

    im = ax.imshow(data_masked, cmap=cmap_obj, origin="lower",
                    vmin=0, vmax=vmax, aspect="equal")

    # フィラメント（skeletonまたはmask）を指定色で重ね描き
    color_rgb = matplotlib.colors.to_rgb(filament_color)
    fil_rgba = np.zeros((*overlay_mask.shape, 4), dtype=float)
    fil_rgba[..., 0] = color_rgb[0]
    fil_rgba[..., 1] = color_rgb[1]
    fil_rgba[..., 2] = color_rgb[2]
    fil_rgba[..., 3] = overlay_mask.astype(float) * filament_alpha
    ax.imshow(fil_rgba, origin="lower", aspect="equal")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"Column Density $N(H_2)$ [cm$^{-2}$]")

    ax.set_title(f"{display_name} (Herschel Space Observatory) "
                 f"d={distance_pc:.0f} pc, beam={beam_arcsec:.1f}\" — {n_filaments} Filaments")
    ax.set_xlabel("x [pixels]")
    ax.set_ylabel("y [pixels]")

    # -----------------------------------------------------------------
    # 5 pc スケールバーの自動計算と描画（右下、横棒のみ）
    # -----------------------------------------------------------------
    pixel_scale_deg = proj_plane_pixel_scales(wcs)[:2].mean()
    pixel_scale_arcsec = pixel_scale_deg * 3600.0
    pc_per_pix = (pixel_scale_arcsec / 206265.0) * distance_pc
    pixels_for_scalebar = scalebar_pc / pc_per_pix

    bar_x_end = nx - int(nx * 0.05)
    bar_x_start = bar_x_end - pixels_for_scalebar
    bar_y = int(ny * 0.08)

    ax.plot([bar_x_start, bar_x_end], [bar_y, bar_y], color="white", lw=3.0)
    ax.text((bar_x_start + bar_x_end) / 2.0, bar_y + (ny * 0.02), f"{scalebar_pc:g} pc",
            color="white", fontsize=12, fontweight="bold", ha="center", va="bottom")

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"科学図を保存しました: {out_path}")


# 科学図(1): 背景の上にskeleton(赤い細線)を重ねたもの
save_science_figure(
    data2d=img_raw,
    overlay_mask=skeleton.astype(bool),
    out_path=os.path.join(vis_dir, f"{fits_name}_science_figure_skeleton_nsig{NSIG}.png"),
    wcs=w,
    distance_pc=distance_pc,
    display_name=display_name,
    beam_arcsec=beam_arcsec,
    n_filaments=len(filaments),  # 上記docstring注記のとおり、要件に応じてassemble後の本数に置き換え可
    filament_color="red",
    filament_alpha=1.0,
)

# 科学図(2): 背景の上にmask(赤い太い半透明領域)を重ねたもの
save_science_figure(
    data2d=img_raw,
    overlay_mask=mask.astype(bool),
    out_path=os.path.join(vis_dir, f"{fits_name}_science_figure_mask_nsig{NSIG}.png"),
    wcs=w,
    distance_pc=distance_pc,
    display_name=display_name,
    beam_arcsec=beam_arcsec,
    n_filaments=len(filaments),
    filament_color="red",
    filament_alpha=0.6,
)

print("DONE")
