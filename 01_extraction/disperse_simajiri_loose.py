#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DisPerSE 実行 + FITS 処理・可視化スクリプト 【loose版】

strict版(元スクリプト)からの変更点:
  - NSIG を 1 -> 0.5 に緩和
  - AV4_THRESH を strict の半分に (3.76e+21 -> 1.88e+21)
  - 出力先を raw_loose / full_loose / vis_loose に変更(strict出力は上書きしない)
  - WIDTH_MODE, WIDTH_PC は不変(マスク幅の話であり、検出の緩さとは別軸のため)

他の領域(ic5146, orionB, polaris, tauN3)に展開する場合は、
AV4_THRESH(各領域のstrict値)を半分にした値に置き換えてください。
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
fits_name = "HGBS_orionB_hires_column_density_map"
fits_path = os.path.join(root_dir, fits_name + ".fits")

display_name = "orionB [loose]"
beam_arcsec = 18.2

DISPERSE_BIN_DIR = "/home1/fujiwara.r/PycharmProjects/DisPerSE/build/src/mse"

def bin_path(name):
    return os.path.join(DISPERSE_BIN_DIR, name) if DISPERSE_BIN_DIR else name

# DisPerSE 固有パラメータ 【変更】 1 -> 0.5 に緩和
NSIG = 0.5

# マスク幅の方針 — strictと同一(不変)
WIDTH_MODE = "skeleton_only"
WIDTH_PC = 0.1
distance_pc = 400.0  # TauN3/Orion B 等に合わせて調整してください

# 【変更】 AV4_THRESHをstrictの半分(3.76e+21 -> 1.88e+21)に
AV4_THRESH = 2.82e+21

# ------------------------------------------------------------
# 出力ディレクトリの作成 【変更】 loose用に分離
# ------------------------------------------------------------
disperse_root = os.path.join(root_dir, fits_name, "disperse")
raw_dir  = os.path.join(disperse_root, "raw_loose")
full_dir = os.path.join(root_dir, fits_name, "full_loose")
vis_dir  = os.path.join(root_dir, fits_name, "vis_loose")

for d in [raw_dir, full_dir, vis_dir]:
    os.makedirs(d, exist_ok=True)


# ============================================================
# ステップ1: FITS読み込み・しきい値(Av < 2相当)処理・正規化
# ============================================================
masked_fits_path = os.path.join(raw_dir, fits_name + "_av_masked_loose.fits")

with fits.open(fits_path, memmap=True) as hdul:
    hdu_header = hdul[0].header.copy()
    img_raw = np.squeeze(np.asarray(hdul[0].data, dtype=np.float64))

w = WCS(hdu_header)
if w.naxis > 2:
    w = w.celestial

finite_mask_all = np.isfinite(img_raw)
img_log_all = np.where(finite_mask_all & (img_raw > 0), np.log10(img_raw), np.nan)
des_log_all = pd.DataFrame(pd.Series(img_log_all[np.isfinite(img_log_all)].ravel()).describe(
    percentiles=[0.005, 0.995]))

data_bak_ml = img_log_all.copy()
data_bak_ml -= des_log_all[0]["0.5%"]
data_bak_ml /= ((des_log_all[0]["99.5%"] - des_log_all[0]["0.5%"]) / 255)
data_bak_ml = np.clip(data_bak_ml, 0, 255)
data_bak_ml[~finite_mask_all] = np.nan

input_u8 = np.where(finite_mask_all, data_bak_ml, 0).astype("uint8")
cv2.imwrite(os.path.join(full_dir, fits_name + "_input_full_ml.png"), input_u8)
print("全体グレースケール（背景用, loose）入力を保存しました。")


img_raw_masked = img_raw.copy()
print(f"DisPerSE用(loose)：閾値 ({AV4_THRESH}) 未満の領域を NaN に置き換えています...")
img_raw_masked[img_raw_masked < AV4_THRESH] = np.nan

masked_hdu = fits.PrimaryHDU(data=img_raw_masked, header=hdu_header)
masked_hdu.writeto(masked_fits_path, overwrite=True)
print(f"DisPerSE用(loose)のマスク済みFITSを保存しました: {masked_fits_path}")


# ============================================================
# ステップ2: mse 実行
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
print("検出フィラメント数(loose):", len(filaments))

shape = img_raw.shape
skeleton = np.zeros(shape, dtype=np.uint8)

for pts in filaments:
    xs, ys = pts[:, 0], pts[:, 1]
    for k in range(len(pts) - 1):
        p1 = (int(round(xs[k])), int(round(ys[k])))
        p2 = (int(round(xs[k + 1])), int(round(ys[k + 1])))
        cv2.line(skeleton, p1, p2, 1, thickness=3)

valid_mask = np.isfinite(img_raw_masked)
kernel = np.ones((5, 5), np.uint8)
valid_mask_eroded = cv2.erode(valid_mask.astype(np.uint8), kernel, iterations=2).astype(bool)

skeleton[~valid_mask_eroded] = 0

cv2.imwrite(os.path.join(full_dir, fits_name + "_skeleton_01.png"), skeleton)
cv2.imwrite(os.path.join(full_dir, fits_name + "_skeleton_label01.png"), skeleton)
print("skeleton画素数（境界ノイズ除去後, loose）:", skeleton.sum())


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
# ステップ6: 可視化
# ============================================================
def save_overlay_png(background_u8, overlay_mask, out_path, color=(1, 0, 0), alpha=1.0):
    fig, ax = plt.subplots(figsize=(16, 16))
    ax.imshow(background_u8, cmap="gray", origin="lower", vmin=0, vmax=255)

    overlay_rgba = np.zeros((*overlay_mask.shape, 4), dtype=float)
    overlay_rgba[..., 0] = color[0]
    overlay_rgba[..., 1] = color[1]
    overlay_rgba[..., 2] = color[2]
    overlay_rgba[..., 3] = overlay_mask.astype(float) * alpha

    ax.imshow(overlay_rgba, origin="lower")
    ax.axis("off")
    fig.savefig(out_path, dpi=150, bbox_inches="tight", pad_inches=0)
    plt.close(fig)

save_overlay_png(
    input_u8, skeleton.astype(bool),
    os.path.join(vis_dir, f"{fits_name}_disperse_skeleton_overlay_nsig{NSIG}.png"),
    color=(1, 0, 0),
    alpha=1.0
)

save_overlay_png(
    input_u8, mask.astype(bool),
    os.path.join(vis_dir, f"{fits_name}_disperse_mask_overlay_nsig{NSIG}.png"),
    color=(1, 0, 0),
    alpha=0.6
)
print("グレースケール背景の上に赤線を重ねた画像(loose)を保存しました。")


def save_science_figure(data2d, overlay_mask, out_path, wcs, distance_pc,
                         display_name, beam_arcsec, n_filaments,
                         scalebar_pc=5.0, vmax=None, cmap="inferno",
                         filament_color="red", filament_alpha=1.0):
    ny, nx = data2d.shape

    if vmax is None:
        vmax = np.nanpercentile(data2d[np.isfinite(data2d)], 99.5)

    fig, ax = plt.subplots(figsize=(11, 8.5))

    data_masked = np.ma.masked_invalid(data2d)
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad("black")

    im = ax.imshow(data_masked, cmap=cmap_obj, origin="lower",
                    vmin=0, vmax=vmax, aspect="equal")

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
    print(f"科学図を保存しました(loose): {out_path}")


save_science_figure(
    data2d=img_raw,
    overlay_mask=skeleton.astype(bool),
    out_path=os.path.join(vis_dir, f"{fits_name}_science_figure_skeleton_nsig{NSIG}.png"),
    wcs=w,
    distance_pc=distance_pc,
    display_name=display_name,
    beam_arcsec=beam_arcsec,
    n_filaments=len(filaments),
    filament_color="red",
    filament_alpha=1.0,
)

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

print("DONE (loose)")
