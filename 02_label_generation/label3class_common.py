#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
label3class_common.py
extract_filfinder_masks.py / extract_disperse_masks.py / generate_3class_labels.py /
generate_and_visualize_3class_map.py の各スクリプトで共通利用するユーティリティ。
同じディレクトリに置いて `from label3class_common import ...` で読み込むこと。

【今回の修正点】
    footprint(観測領域)の境界付近を侵食(erosion)して除外する
    compute_eroded_footprint() を追加。loose閾値のマスクは境界付近の淡い
    ノイズを拾いやすいため、この関数で境界をあらかじめ無効化してから
    使うことを想定している。
"""
import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import binary_dilation, binary_erosion, label, find_objects

# ----------------------------------------------------------------------
# 設定(共通)
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

# フィラメント幅[pc]。skeletonを面マスク化する際の膨張幅。
# 【要確認】0.1pcから0.03pcに変更されている。意図的な変更か要確認。
WIDTH_PC = 0.1
# footprint(観測領域)の境界から何pixel分を無効化するか。
# loose閾値マスクが境界付近の淡いノイズを拾う対策。
FOOTPRINT_ERODE_PX = 10

# フィラメントとみなす最小長さ[pc]。これ未満の連結成分はノイズとして除外する。
# 特にDisPerSEは希薄な領域(Polaris等)でごく短いノイズ的セグメントを
# 大量に生成しやすいため、物理幅への膨張前に必ずこのフィルタを通す。
MIN_FILAMENT_LENGTH_PC = 0.3

# 中間ファイル(strict/loose物理マスク)の保存先サブディレクトリ名
MASKS_SUBDIR = "masks_3class"


def guess_distance_pc(fits_name, distance_table=DISTANCE_PC):
    name_lower = fits_name.lower()
    for key, dist in distance_table.items():
        if key in name_lower:
            return dist
    print(f"[警告] fits_name='{fits_name}' から距離を推定できませんでした。")
    return None


def load_wcs_and_pixelscale(fits_path):
    if not os.path.exists(fits_path):
        raise FileNotFoundError(f"FITSファイルが見つかりません: {fits_path}")
    header = fits.getheader(fits_path)
    w = WCS(header)
    pixel_scale_deg = proj_plane_pixel_scales(w)[:2].mean()
    return w, pixel_scale_deg * 3600.0


def compute_pc_per_pixel(fits_path, distance_pc):
    _, pixel_scale_arcsec = load_wcs_and_pixelscale(fits_path)
    return (pixel_scale_arcsec / 206265.0) * distance_pc


def dilate_to_width(skeleton_mask, pc_per_pixel, width_pc=WIDTH_PC):
    """1pxスケルトンを物理幅width_pc[pc]でdilateして面マスクにする。"""
    radius_px = max(1, int(round((width_pc / pc_per_pixel) / 2)))
    structure = np.ones((2 * radius_px + 1, 2 * radius_px + 1), dtype=bool)
    return binary_dilation(skeleton_mask.astype(bool), structure=structure)


def filter_short_segments(mask, pc_per_pixel, min_length_pc=MIN_FILAMENT_LENGTH_PC,
                           label_name=""):
    """
    2値マスク(1pxスケルトン)を連結成分に分解し、各成分のバウンディングボックス
    対角長を物理長[pc]に換算して、min_length_pc未満の成分を除去する。

    特にDisPerSEは希薄・乱流的な領域(Polaris等)で大量の短いノイズ的セグメントを
    生成しやすいため、dilate_to_width()で物理幅に太らせる前に必ず適用すること。

    連結成分判定は8方向連結(対角も接続とみなす)を使う。1px幅のスケルトンは
    斜めにしかつながっていない箇所が多く、デフォルトの4方向連結だと本来ひと続きの
    フィラメントが細かく分断され、誤って大量除去されてしまうため。
    """
    structure_8conn = np.ones((3, 3), dtype=int)
    labeled, n_total = label(mask, structure=structure_8conn)
    if n_total == 0:
        print(f"  [{label_name}] 連結成分なし。フィルタ処理をスキップ。")
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

    print(f"  [{label_name}] {n_total}個中{n_removed}個を "
          f"{min_length_pc}pc未満として除外(残り{n_total - n_removed}個)")
    return keep


def compute_eroded_footprint(img_raw, erode_px=FOOTPRINT_ERODE_PX):
    """
    観測領域(NaNでない画素)を計算し、境界からerode_px分だけ侵食した
    「安全な内側領域」を返す。loose閾値のマスクをこの領域でANDすることで、
    footprint境界のギザギザに沿ったノイズ検出を除外できる。
    """
    footprint = np.isfinite(img_raw)
    if erode_px <= 0:
        return footprint
    structure = np.ones((2 * erode_px + 1, 2 * erode_px + 1), dtype=bool)
    return binary_erosion(footprint, structure=structure)


def masks_out_dir(root_dir, fits_name):
    d = os.path.join(root_dir, fits_name, MASKS_SUBDIR)
    os.makedirs(d, exist_ok=True)
    return d


def save_mask_npy(root_dir, fits_name, method, loose, mask):
    """
    method: "filfinder" or "disperse"
    loose : True/False
    保存先: root_dir/{fits_name}/masks_3class/{method}_{strict|loose}.npy
    """
    tag = "loose" if loose else "strict"
    out_path = os.path.join(masks_out_dir(root_dir, fits_name), f"{method}_{tag}.npy")
    np.save(out_path, mask)
    print(f"  -> 保存: {out_path} (True画素数={mask.sum()})")
    return out_path


def load_mask_npy(root_dir, fits_name, method, loose):
    tag = "loose" if loose else "strict"
    path = os.path.join(masks_out_dir(root_dir, fits_name), f"{method}_{tag}.npy")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{method}({tag})の中間マスクが見つかりません: {path}\n"
            f"先に extract_{method}_masks.py を実行してください。"
        )
    return np.load(path)
