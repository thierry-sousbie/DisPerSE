#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pickle
import numpy as np
import cv2
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import binary_dilation, label, find_objects
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from reproject import reproject_interp

# ----------------------------------------------------------------------
# 設定
# ----------------------------------------------------------------------

ROOT = "/mnt/sp4t/fujiwara.r"
FITS_NAMES = [
    "HGBS_aquilaM2_hires_column_density_map",
    "HGBS_ic5146_hires_column_density_map",
    "HGBS_orionB_hires_column_density_map",
    "HGBS_polaris_hires_column_density_map",
    "HGBS_tauN3_hires_column_density_map",
]

DISTANCE_PC = {
    "aquilam2": 440,
    "ic5146": 600,
    "orionb": 400,
    "polaris": 150,
    "taun3": 140,
}

WIDTH_PC = 0.1
MIN_FILAMENT_LENGTH_PC = 0.3

def herschel_path_for(root_dir, fits_name):
    base = fits_name.replace("_hires_column_density_map", "")
    return os.path.join(root_dir, base + "_skeleton_map.fits")

def guess_distance_pc(fits_name, distance_table=DISTANCE_PC):
    name_lower = fits_name.lower()
    for key, dist in distance_table.items():
        if key in name_lower:
            return dist
    return None

def load_wcs_and_pixelscale(fits_path):
    with fits.open(fits_path) as hdul:
        header = hdul[0].header
        w = WCS(header)
        pixel_scale_deg = proj_plane_pixel_scales(w)[:2].mean()
        pixel_scale_arcsec = pixel_scale_deg * 3600.0
    return w, header, pixel_scale_arcsec

def load_column_density(fits_path):
    return fits.getdata(fits_path)

def dilate_to_width(skeleton_mask, pc_per_pixel, width_pc=WIDTH_PC):
    radius_px = max(1, int(round((width_pc / pc_per_pixel) / 2)))
    structure = np.ones((2 * radius_px + 1, 2 * radius_px + 1), dtype=bool)
    return binary_dilation(skeleton_mask.astype(bool), structure=structure)

def filter_short_segments(mask, pc_per_pixel, min_length_pc=MIN_FILAMENT_LENGTH_PC, label_name=""):
    structure_8conn = np.ones((3, 3), dtype=int)
    labeled, n_total = label(mask, structure=structure_8conn)
    if n_total == 0:
        return mask.copy()

    objects = find_objects(labeled)
    keep = np.zeros_like(mask, dtype=bool)
    n_removed = 0

    for i, slc in enumerate(objects, start=1):
        if slc is None:
            continue
        height_px = slc[0].stop - slc[0].start
        width_px = slc[1].stop - slc[1].start
        diag_px = np.hypot(height_px, width_px)
        length_pc = diag_px * pc_per_pixel
        if length_pc >= min_length_pc:
            keep |= (labeled == i)
        else:
            n_removed += 1

    print(f"  [{label_name}] {n_total}個中{n_removed}個を除外 (残り{n_total - n_removed}個)")
    return keep

def load_disperse_mask(full_dir, fits_name, pc_per_pixel):
    skel_path = os.path.join(full_dir, fits_name + "_skeleton_01.png")
    skeleton_raw = cv2.imread(skel_path, cv2.IMREAD_GRAYSCALE)
    if skeleton_raw is None:
        raise FileNotFoundError(f"DisPerSE画像が見つかりません: {skel_path}")
    skeleton = skeleton_raw.astype(bool)
    skeleton = filter_short_segments(skeleton, pc_per_pixel, label_name="DisPerSE skeleton")
    return dilate_to_width(skeleton, pc_per_pixel)

def load_and_reproject_herschel(herschel_path, ref_header):
    """HerschelのFITSを読み込み、ベースのFITSヘッダー(ref_header)のWCSに位置合わせする"""
    if not os.path.exists(herschel_path):
        raise FileNotFoundError(f"Herschelマップが見つかりません: {herschel_path}")

    with fits.open(herschel_path) as hdul:
        data = hdul[0].data
        header = hdul[0].header

    mask = (data > 0).astype(float)
    ref_wcs = WCS(ref_header)
    ref_shape = (ref_header["NAXIS2"], ref_header["NAXIS1"])

    # WCSをもとに座標を再投影
    reprojected, _ = reproject_interp(
        (mask, WCS(header)),
        ref_wcs,
        shape_out=ref_shape,
        order="nearest-neighbor"
    )
    return (np.nan_to_num(reprojected) > 0.5).astype(bool)

def draw_scalebar(ax, nx, ny, pixel_scale_arcsec, distance_pc, scalebar_pc=5.0, color="white"):
    if pixel_scale_arcsec is None or distance_pc is None:
        return
    pc_per_pix = (pixel_scale_arcsec / 206265.0) * distance_pc
    pixels_for_scalebar = scalebar_pc / pc_per_pix
    bar_x_end = nx - int(nx * 0.05)
    bar_x_start = bar_x_end - pixels_for_scalebar
    bar_y = int(ny * 0.08)

    ax.plot([bar_x_start, bar_x_end], [bar_y, bar_y], color=color, lw=3.0)
    ax.text((bar_x_start + bar_x_end) / 2.0, bar_y + (ny * 0.02),
            f"{scalebar_pc:g} pc", color=color, fontsize=12,
            fontweight="bold", ha="center", va="bottom")

def plot_two_way_overlay(method_mask, herschel_mask, method_name, method_color,
                         fits_name, out_path, background=None, bg_cmap="inferno",
                         method_vis_dilate_iter=0, herschel_vis_dilate_iter=0,
                         pixel_scale_arcsec=None, distance_pc=None, scalebar_pc=5.0):
    ny, nx = method_mask.shape
    fig, ax = plt.subplots(figsize=(10, 8))

    if background is not None:
        finite_vals = background[np.isfinite(background)]
        bg_vmax = np.nanpercentile(finite_vals, 99.5) if finite_vals.size > 0 else 1.0
        im_bg = ax.imshow(np.nan_to_num(background, nan=0.0) * 0.6,
                          cmap=bg_cmap, origin="lower", vmin=0, vmax=bg_vmax)
        cbar = fig.colorbar(im_bg, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(r"Column Density $N(H_2)$ [cm$^{-2}$]")

    method_vis = method_mask
    if method_vis_dilate_iter > 0:
        method_vis = binary_dilation(method_mask, iterations=method_vis_dilate_iter)

    herschel_vis = herschel_mask
    if herschel_vis_dilate_iter > 0:
        herschel_vis = binary_dilation(herschel_mask, iterations=herschel_vis_dilate_iter)

    both = method_vis & herschel_vis
    method_only = method_vis & ~herschel_vis
    herschel_only = ~method_vis & herschel_vis

    overlay = np.zeros((*method_mask.shape, 4), dtype=float)
    yellow = mcolors.to_rgb("yellow")
    method_rgb = mcolors.to_rgb(method_color)
    cyan = mcolors.to_rgb("cyan")

    for mask_, color in [(both, yellow), (method_only, method_rgb), (herschel_only, cyan)]:
        overlay[mask_, 0] = color[0]
        overlay[mask_, 1] = color[1]
        overlay[mask_, 2] = color[2]
        overlay[mask_, 3] = 0.9

    ax.imshow(overlay, origin="lower")
    ax.set_title(f"{fits_name}: {method_name} vs Herschel\n"
                 f"Yellow=Both / {method_name} only / Cyan=Herschel only")
    ax.set_xlabel("x [pixels]")
    ax.set_ylabel("y [pixels]")

    from matplotlib.patches import Patch
    legend_elems = [
        Patch(facecolor="yellow", label="Both"),
        Patch(facecolor=method_color, label=f"{method_name} only"),
        Patch(facecolor="cyan", label="Herschel only"),
    ]
    ax.legend(handles=legend_elems, loc="upper right", fontsize=8)
    draw_scalebar(ax, nx, ny, pixel_scale_arcsec, distance_pc, scalebar_pc=scalebar_pc)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> 保存完了: {out_path}")

# ----------------------------------------------------------------------
# メイン処理
# ----------------------------------------------------------------------

def run_comparison_for_fits_name(root_dir, fits_name):
    print(f"\n=== {fits_name} ===")

    region_dir = os.path.join(root_dir, fits_name)
    full_dir = os.path.join(region_dir, "full")
    out_dir = os.path.join(region_dir, "compare_herschel")
    os.makedirs(out_dir, exist_ok=True)

    original_fits_path = os.path.join(root_dir, fits_name + ".fits")
    herschel_path = herschel_path_for(root_dir, fits_name)

    distance_pc = guess_distance_pc(fits_name)
    if distance_pc is None:
        print(f"[スキップ] {fits_name}: 距離不明")
        return None

    wcs_obj, ref_header, pixel_scale_arcsec = load_wcs_and_pixelscale(original_fits_path)
    pc_per_pixel = (pixel_scale_arcsec / 206265.0) * distance_pc
    img_raw = load_column_density(original_fits_path)

    # 1. DisPerSEの読み込み
    disperse_mask = load_disperse_mask(full_dir, fits_name, pc_per_pixel)

    # 2. Herschelの読み込み ＆ 天体座標(WCS)ベースでの位置合わせ(Reproject)
    herschel_mask = load_and_reproject_herschel(herschel_path, ref_header)
    herschel_mask = filter_short_segments(herschel_mask, pc_per_pixel, label_name="Herschel skeleton")

    # DisPerSEと背景画像のサイズ確認（もしDisPerSEの画像サイズが画像FITSと異なれば警告）
    if disperse_mask.shape != img_raw.shape:
        print(f"[警告] DisPerSE画像のサイズ ({disperse_mask.shape}) と背景FITS ({img_raw.shape}) が一致していません！")
        # 必要に応じてリサイズ処理
        disperse_mask = cv2.resize(disperse_mask.astype(np.uint8), (img_raw.shape[1], img_raw.shape[0]), interpolation=cv2.INTER_NEAREST).astype(bool)

    print(f"  DisPerSE(physical) True画素数: {disperse_mask.sum()}")
    print(f"  Herschel(reprojected) True画素数: {herschel_mask.sum()}")

    vis_dilate_px = max(1, int(round((WIDTH_PC / pc_per_pixel) / 2)))

    plot_two_way_overlay(
        disperse_mask, herschel_mask, "DisPerSE", "red",
        fits_name, os.path.join(out_dir, f"{fits_name}_disperse_vs_herschel.png"),
        background=img_raw,
        method_vis_dilate_iter=0,
        herschel_vis_dilate_iter=vis_dilate_px,
        pixel_scale_arcsec=pixel_scale_arcsec, distance_pc=distance_pc, scalebar_pc=5.0,
    )

def main():
    for fits_name in FITS_NAMES:
        try:
            run_comparison_for_fits_name(ROOT, fits_name)
        except FileNotFoundError as e:
            print(f"[スキップ] {fits_name}: {e}")
            continue

if __name__ == "__main__":
    main()
