#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
compare_vs_herschel.py

FilFinder / DisPerSE それぞれと、Herschel(公開)スケルトンマップとの比較のみを行うスクリプト。
(FilFinder vs DisPerSE の手法間比較は行わない)
定量指標(IoU/Precision/Recall)は計算せず、色分けオーバーレイ画像の出力のみ行う。

比較の組み合わせ:
    1. FilFinder  vs Herschel   (外部基準との照合)
    2. DisPerSE   vs Herschel   (外部基準との照合)

パス構成は compare_disperse_filfinder_fixed.py と揃えてある:
    root_dir/{fits_name}.fits                              柱密度FITS(背景・WCS取得元)
    root_dir/{fits_name}/cache/{fits_name}_FilFinder2D-medskel.pkl      FilFinderスケルトン
    root_dir/{fits_name}/cache/{fits_name}_FilFinder2D-create_mask.pkl  FilFinderマスク
    root_dir/{fits_name}/full/{fits_name}_skeleton_01.png               DisPerSEスケルトン
    root_dir/{fits_name}/full/{fits_name}_mask_01.png                   DisPerSEマスク
    root_dir/{HGBS_xxx}_skeleton_map.fits                               Herschelスケルトン
        (例: HGBS_aquilaM2_hires_column_density_map -> HGBS_aquilaM2_skeleton_map.fits)
出力は root_dir/{fits_name}/compare_herschel/ 配下に保存する。

前提・要修正箇所:
    - ASCIIカタログ(座標リスト)の場合は load_herschel_mask() の
      書き換えが必要(コメント参照)。
    - FilFinder2Dオブジェクトの属性名(.skeleton, .mask)はバージョンにより
      異なる場合があるため要確認。

依存:
    numpy, astropy, scipy, matplotlib, opencv-python
    pip install numpy astropy scipy matplotlib opencv-python --break-system-packages
"""

import os
import pickle
import numpy as np
import cv2
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import binary_dilation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# ----------------------------------------------------------------------
# 設定(compare_disperse_filfinder_fixed.py と揃えてある)
# ----------------------------------------------------------------------

ROOT = "/mnt/sp4t/fujiwara.r"
FITS_NAMES = [
    "HGBS_aquilaM2_hires_column_density_map",
    "HGBS_ic5146_hires_column_density_map",
    "HGBS_orionB_hires_column_density_map",
    "HGBS_polaris_hires_column_density_map",
    "HGBS_tauN3_hires_column_density_map",
]

# fits_name(小文字化)に含まれる文字列 -> 距離[pc]。要確認・要更新。
DISTANCE_PC = {
    "aquilam2": 440,
    "ic5146": 600,
    "orionb": 400,
    "polaris": 150,
    "taun3": 140,
}

# フィラメント幅[pc]。DisPerSE側マスク生成スクリプトと同じ値に揃える。
WIDTH_PC = 0.1

# Herschelスケルトンマップのパス。root_dir直下に
# "{HGBS_xxx}_skeleton_map.fits" として置かれている想定。
# 例: HGBS_aquilaM2_hires_column_density_map -> HGBS_aquilaM2_skeleton_map.fits
def herschel_path_for(root_dir, fits_name):
    base = fits_name.replace("_hires_column_density_map", "")
    return os.path.join(root_dir, base + "_skeleton_map.fits")


def guess_distance_pc(fits_name, distance_table=DISTANCE_PC):
    name_lower = fits_name.lower()
    for key, dist in distance_table.items():
        if key in name_lower:
            return dist
    print(f"[警告] fits_name='{fits_name}' から距離を推定できませんでした。")
    return None


# ----------------------------------------------------------------------
# ユーティリティ
# ----------------------------------------------------------------------

def load_wcs_and_pixelscale(fits_path):
    """FITSヘッダからWCSとpixelスケール[arcsec/pix]を取得する。"""
    if not os.path.exists(fits_path):
        raise FileNotFoundError(f"FITSファイルが見つかりません: {fits_path}")
    header = fits.getheader(fits_path)
    w = WCS(header)
    pixel_scale_deg = proj_plane_pixel_scales(w)[:2].mean()
    pixel_scale_arcsec = pixel_scale_deg * 3600.0
    return w, pixel_scale_arcsec


def load_column_density(fits_path):
    """柱密度FITSのデータ本体(2次元配列, 物理値cm^-2)を読み込む。"""
    if not os.path.exists(fits_path):
        raise FileNotFoundError(f"柱密度FITSが見つかりません: {fits_path}")
    return fits.getdata(fits_path)


def dilate_to_width(skeleton_mask, pc_per_pixel, width_pc=WIDTH_PC):
    """1pxスケルトンを物理幅でdilationして面マスクにする(physicalモード)。"""
    radius_px = max(1, int(round((width_pc / pc_per_pixel) / 2)))
    structure = np.ones((2 * radius_px + 1, 2 * radius_px + 1), dtype=bool)
    return binary_dilation(skeleton_mask.astype(bool), structure=structure)


# ----------------------------------------------------------------------
# 読み込み関数(環境に合わせて要調整)
# ----------------------------------------------------------------------

def load_filfinder_mask(cache_dir, fits_name, mode, pc_per_pixel):
    """
    root_dir/{fits_name}/cache/{fits_name}_FilFinder2D-*.pkl から読み込む。
    mode: "skeleton_only" -> medskel.pkl の .skeleton (1px)
          "physical"      -> create_mask.pkl の .mask (既にphysical幅) を優先し、
                              無ければ .skeleton を0.1pcでdilationして代用
    """
    skel_path = os.path.join(cache_dir, fits_name + "_FilFinder2D-medskel.pkl")
    mask_path = os.path.join(cache_dir, fits_name + "_FilFinder2D-create_mask.pkl")

    with open(skel_path, "rb") as f:
        skel_obj = pickle.load(f)
    skeleton = skel_obj.skeleton.astype(bool)

    if mode == "skeleton_only":
        return skeleton
    elif mode == "physical":
        if os.path.exists(mask_path):
            with open(mask_path, "rb") as f:
                mask_obj = pickle.load(f)
            return mask_obj.mask.astype(bool)
        return dilate_to_width(skeleton, pc_per_pixel)
    else:
        raise ValueError(f"不明なmode: {mode}")


def load_disperse_mask(full_dir, fits_name, mode, pc_per_pixel):
    """
    root_dir/{fits_name}/full/{fits_name}_{skeleton,mask}_01.png から読み込む。
    """
    skel_path = os.path.join(full_dir, fits_name + "_skeleton_01.png")
    mask_path = os.path.join(full_dir, fits_name + "_mask_01.png")

    skeleton_raw = cv2.imread(skel_path, cv2.IMREAD_GRAYSCALE)
    if skeleton_raw is None:
        raise FileNotFoundError(f"見つかりません: {skel_path}")
    skeleton = skeleton_raw.astype(bool)

    if mode == "skeleton_only":
        return skeleton
    elif mode == "physical":
        mask_raw = cv2.imread(mask_path, cv2.IMREAD_GRAYSCALE)
        if mask_raw is not None:
            return mask_raw.astype(bool)
        return dilate_to_width(skeleton, pc_per_pixel)
    else:
        raise ValueError(f"不明なmode: {mode}")


def load_herschel_mask(path):
    """
    公開Herschelスケルトンマップを読み込む。

    NOTE: カタログがFITS形式のマスク/スケルトンでない場合
    （例: getsf由来のASCII座標リストなど）は、ここで
    座標 -> pixelマスクへの変換処理を実装してください。
    例:
        table = ascii.read(path)
        wcs = WCS(reference_header)
        x, y = wcs.world_to_pixel(SkyCoord(table['ra'], table['dec'], unit='deg'))
        mask = np.zeros(shape, dtype=bool)
        mask[y.astype(int), x.astype(int)] = True
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Herschelスケルトンマップが見つかりません: {path}\n"
            f"root_dir直下の{{HGBS_xxx}}_skeleton_map.fitsという命名で"
            f"置かれているか確認してください。"
        )
    if path.endswith(".fits") or path.endswith(".fit"):
        data = fits.getdata(path)
        return (data > 0).astype(bool)
    else:
        raise NotImplementedError(
            f"Herschelカタログの形式 {path} は未対応です。"
            "load_herschel_mask()内のコメントを参照して実装してください。"
        )


def reproject_to_reference(mask, mask_header, ref_header):
    """
    グリッドが異なる場合に、referenceのWCS/shapeに合わせてreprojectする。
    2値マスクなので order='nearest-neighbor' を使う。
    """
    from reproject import reproject_interp
    ref_wcs = WCS(ref_header)
    ref_shape = (ref_header["NAXIS2"], ref_header["NAXIS1"])
    reprojected, _ = reproject_interp(
        (mask.astype(float), WCS(mask_header)),
        ref_wcs, shape_out=ref_shape, order="nearest-neighbor"
    )
    return (reprojected > 0.5).astype(bool)



# ----------------------------------------------------------------------
# 可視化(2値比較: 黄=両方一致、method_color=手法のみ、水色=Herschelのみ)
# ----------------------------------------------------------------------

def plot_two_way_overlay(method_mask, herschel_mask, method_name, method_color,
                          fits_name, out_path, background=None, bg_cmap="inferno",
                          method_vis_dilate_iter=0, herschel_vis_dilate_iter=0):
    """
    method_mask(FilFinderまたはDisPerSE) と herschel_mask の2値オーバーレイを保存する。
    黄=両方一致、method_color=method_maskのみ、水色(cyan)=herschel_maskのみ。
    背景に柱密度の物理値(background)を渡せばinfernoカラーマップで表示する。

    method_vis_dilate_iter, herschel_vis_dilate_iter:
        可視化のためだけにmaskを膨張させる回数。1px幅のskeleton同士だと
        線が細すぎて見えなくなるため、表示用にのみ太らせる。
        指標(IoU/Precision/Recall)の計算には使わない(呼び出し元で
        別途、膨張前のmaskから計算しておくこと)。
    """
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

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> オーバーレイ画像を保存: {out_path}")
    print(f"     黄(Both)={both.sum()}px, {method_color}({method_name} only)="
          f"{method_only.sum()}px, cyan(Herschel only)={herschel_only.sum()}px")


# ----------------------------------------------------------------------
# メイン処理
# ----------------------------------------------------------------------

def run_comparison_for_fits_name(root_dir, fits_name):
    """
    FilFinderはskeleton_only(1px)、DisPerSEはphysical(0.1pc dilation)で固定して
    それぞれHerschelスケルトンマップと比較する。
    """
    print(f"\n=== {fits_name} ===")

    region_dir = os.path.join(root_dir, fits_name)
    cache_dir = os.path.join(region_dir, "cache")
    full_dir = os.path.join(region_dir, "full")
    out_dir = os.path.join(region_dir, "compare_herschel")
    os.makedirs(out_dir, exist_ok=True)

    original_fits_path = os.path.join(root_dir, fits_name + ".fits")
    herschel_path = herschel_path_for(root_dir, fits_name)

    distance_pc = guess_distance_pc(fits_name)
    if distance_pc is None:
        print(f"[スキップ] {fits_name}: 距離が不明なためスキップします。")
        return None

    wcs_obj, pixel_scale_arcsec = load_wcs_and_pixelscale(original_fits_path)
    pc_per_pixel = (pixel_scale_arcsec / 206265.0) * distance_pc
    img_raw = load_column_density(original_fits_path)

    fil_mask = load_filfinder_mask(cache_dir, fits_name, "skeleton_only", pc_per_pixel)
    disperse_mask = load_disperse_mask(full_dir, fits_name, "physical", pc_per_pixel)
    herschel_mask = load_herschel_mask(herschel_path)

    # shapeが一致しない場合は reproject_to_reference() で先に揃えること
    assert fil_mask.shape == disperse_mask.shape == herschel_mask.shape, (
        f"マスクのshapeが一致していません: "
        f"FilFinder={fil_mask.shape}, DisPerSE={disperse_mask.shape}, "
        f"Herschel={herschel_mask.shape}。reproject_to_reference()で揃えてください。"
    )

    print(f"  fil_mask(FilFinder skeleton)  の True画素数: {fil_mask.sum()}")
    print(f"  disperse_mask(DisPerSE physical) の True画素数: {disperse_mask.sum()}")
    print(f"  herschel_mask(Herschel skeleton) の True画素数: {herschel_mask.sum()}")

    # 可視化専用: 1px幅のskeleton(FilFinder, Herschel)はそのままだと細すぎて
    # 見えないため、WIDTH_PC相当の太さに膨張させて表示する。
    # DisPerSEは既にphysical(0.1pc)で太いため、そのまま(膨張なし)で使う。
    vis_dilate_px = max(1, int(round((WIDTH_PC / pc_per_pixel) / 2)))
    print(f"  可視化用の膨張半径: {vis_dilate_px}px (WIDTH_PC={WIDTH_PC}pc換算)")

    plot_two_way_overlay(
        fil_mask, herschel_mask, "FilFinder", "red",
        fits_name, os.path.join(out_dir, f"{fits_name}_filfinder_vs_herschel.png"),
        background=img_raw,
        method_vis_dilate_iter=vis_dilate_px,
        herschel_vis_dilate_iter=vis_dilate_px,
    )
    plot_two_way_overlay(
        disperse_mask, herschel_mask, "DisPerSE", "red",
        fits_name, os.path.join(out_dir, f"{fits_name}_disperse_vs_herschel.png"),
        background=img_raw,
        method_vis_dilate_iter=0,
        herschel_vis_dilate_iter=vis_dilate_px,
    )

    return {"fits_name": fits_name}



def main():
    for fits_name in FITS_NAMES:
        try:
            run_comparison_for_fits_name(ROOT, fits_name)
        except FileNotFoundError as e:
            print(f"[スキップ] {fits_name}: {e}")
            continue


if __name__ == "__main__":
    main()
