#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
generate_and_visualize_3class_map.py

FilFinder と DisPerSE の検出結果（strict / loose）から3クラス合成ラベルを作成し、
確定コア（赤色: Class 1）と 候補（水色: Class 2）を重ね合わせた全体可視化マップ（PNG）および
3クラスラベルFITSを出力する専用スクリプト。
"""

import os
import glob
import pickle
import numpy as np
import cv2
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch

from label3class_common import (
    guess_distance_pc, load_wcs_and_pixelscale, dilate_to_width, WIDTH_PC,
    compute_eroded_footprint, FOOTPRINT_ERODE_PX,
    filter_short_segments, MIN_FILAMENT_LENGTH_PC,
)

# ============================================================
# 設定 (CONFIG)
# ============================================================
ROOT_DIR  = "/mnt/sp4t/fujiwara.r"
FITS_NAME = "HGBS_orionB_hires_column_density_map"  # 対象領域に合わせて変更

REGION_DIR = os.path.join(ROOT_DIR, FITS_NAME)
OUT_DIR    = os.path.join(REGION_DIR, "vis_3class_map")

# 可視化色定義 (RGBA)
COLOR_CORE_RED  = (1.0, 1.0, 0.5, 1.0)    # Class 1: 確定コア (チェリーピンク/赤)
COLOR_CAND_CYAN = (0.22, 0.74, 0.97, 0.7)  # Class 2: 候補/周縁 (水色)

# 誤って全域マップ画像（column_density等）をマスクとして読まないための除外キーワード
EXCLUDE_KEYWORDS = ["input_full", "av_masked", "full_ml"]

# ノイズ除去閾値（膨張後にこのピクセル数未満の孤立アイランドはノイズとして破棄）
MIN_BLOB_AREA_PX = 15

# Aquila以外の画像端をNaN扱いにするフレーム幅 [pixel]
BORDER_NAN_PX = 10


# ============================================================
# 補助関数
# ============================================================
def remove_small_blobs(mask, min_size_px=15):
    """膨張によってできた小さな孤立ゴミ（点ノイズや1pxヒゲの肥大化）を除去"""
    if mask.sum() == 0 or min_size_px <= 0:
        return mask
    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(mask.astype(np.uint8))
    clean_mask = np.zeros_like(mask, dtype=bool)
    for i in range(1, num_labels):
        if stats[i, cv2.CC_STAT_AREA] >= min_size_px:
            clean_mask[labels == i] = True
    return clean_mask


def draw_scalebar(ax, nx, ny, pixel_scale_arcsec, distance_pc, scalebar_pc=5.0, color="white"):
    """右下に物理長スケールバーを描画"""
    if pixel_scale_arcsec is None or distance_pc is None:
        return

    pc_per_pix = (pixel_scale_arcsec / 206265.0) * distance_pc
    pixels_for_scalebar = scalebar_pc / pc_per_pix

    bar_x_end = nx - int(nx * 0.05)
    bar_x_start = bar_x_end - pixels_for_scalebar
    bar_y = int(ny * 0.05)

    ax.plot([bar_x_start, bar_x_end], [bar_y, bar_y], color=color, lw=4.5)
    ax.text((bar_x_start + bar_x_end) / 2.0, bar_y + (ny * 0.02),
            f"{scalebar_pc:g} pc", color=color, fontsize=12,
            fontweight="bold", ha="center", va="bottom")


# ============================================================
# 1. マスク自動読み込み関数
# ============================================================
def load_mask_from_dir(search_dir, target_shape):
    """指定ディレクトリから骨格マスク (.fits, .npy, .png, .pkl) を探索してロード"""
    if not os.path.exists(search_dir):
        print(f"  [警告] ディレクトリが存在しません: {search_dir}")
        return np.zeros(target_shape, dtype=bool)

    patterns = ["*.fits", "*.npy", "*.png", "*.pkl"]
    found_files = []
    for pat in patterns:
        found_files.extend(glob.glob(os.path.join(search_dir, pat)))

    def is_skeleton_file(filepath):
        fname = os.path.basename(filepath).lower()
        return ("skeleton" in fname) or ("medskel" in fname) or ("skel" in fname)

    def file_priority(filepath):
        fname = os.path.basename(filepath).lower()
        if is_skeleton_file(filepath):
            return 1
        if any(ex in fname for ex in EXCLUDE_KEYWORDS):
            return 999
        if "mask" in fname:
            return 2
        return 10

    sorted_files = sorted(found_files, key=file_priority)

    for path in sorted_files:
        fname = os.path.basename(path).lower()

        if not is_skeleton_file(filepath=path) and any(ex in fname for ex in EXCLUDE_KEYWORDS):
            print(f"  [スキップ (全域画像除外)]: {os.path.basename(path)}")
            continue

        try:
            data = None
            if path.endswith(".fits"):
                data = fits.getdata(path)
            elif path.endswith(".npy"):
                data = np.load(path)
            elif path.endswith(".png"):
                data = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
            elif path.endswith(".pkl"):
                with open(path, "rb") as f:
                    obj = pickle.load(f)
                    if isinstance(obj, dict):
                        data = obj.get("skeleton", obj.get("medskel", None))
                    else:
                        data = getattr(obj, "skeleton", getattr(obj, "medskel", getattr(obj, "mask", None)))

            if data is None:
                continue

            data = np.squeeze(data)
            mask = (data > 0).astype(bool)

            valid_ratio = mask.sum() / mask.size
            if valid_ratio > 0.3:
                print(f"  [スキップ (全域塗潰し検知: {valid_ratio*100:.1f}%)]: {os.path.basename(path)}")
                continue

            if mask.shape != target_shape:
                print(f"  [調整] マスクリサイズ: {mask.shape} -> {target_shape}")
                mask_u8 = cv2.resize(mask.astype(np.uint8), (target_shape[1], target_shape[0]),
                                     interpolation=cv2.INTER_NEAREST)
                mask = (mask_u8 > 0)

            if mask.sum() > 0:
                print(f"  -> 読込成功: {os.path.basename(path)} (有効ピクセル: {mask.sum()})")
                return mask

        except Exception:
            continue

    print(f"  [エラー] 有効な骨格マスクが見つかりませんでした: {search_dir}")
    return np.zeros(target_shape, dtype=bool)


# ============================================================
# 2. 3クラスラベル合成処理
# ============================================================
def build_3class_label(img_raw, header, fits_path):
    """FilFinder と DisPerSE の検出結果から 3クラスラベルを合成"""
    shape = img_raw.shape

    # --- 距離・pixelスケールを取得 ---
    distance_pc = guess_distance_pc(FITS_NAME)
    if distance_pc is None:
        raise RuntimeError(
            f"距離が推定できませんでした。FITS_NAME='{FITS_NAME}' が"
            f"label3class_common.DISTANCE_PC のキーに対応しているか確認してください。"
        )
    _, pixel_scale_arcsec = load_wcs_and_pixelscale(fits_path)
    pc_per_pixel = (pixel_scale_arcsec / 206265.0) * distance_pc
    print(f"\n[物理幅換算] distance_pc={distance_pc}, "
          f"pixel_scale_arcsec={pixel_scale_arcsec:.3f}\", "
          f"pc_per_pixel={pc_per_pixel:.5f}, WIDTH_PC={WIDTH_PC}pc")

    # --- footprint境界の侵食(erosion)マスクを計算 ---
    footprint_eroded = compute_eroded_footprint(img_raw, erode_px=FOOTPRINT_ERODE_PX)
    print(f"[footprint] 境界{FOOTPRINT_ERODE_PX}px分を侵食した安全な内側領域を計算しました "
          f"(有効画素数={footprint_eroded.sum()} / 全体footprint画素数={np.isfinite(img_raw).sum()})")

    dir_fil_strict = os.path.join(REGION_DIR, "cache")         # FilFinder strict
    dir_fil_loose  = os.path.join(REGION_DIR, "cache_loose")  # FilFinder loose
    dir_dis_strict = os.path.join(REGION_DIR, "full")         # DisPerSE strict
    dir_dis_loose  = os.path.join(REGION_DIR, "full_loose")   # DisPerSE loose

    print("\n[マスク読み込み状況]")
    print("1. FilFinder Strict  (cache/):")
    m_fil_strict = load_mask_from_dir(dir_fil_strict, shape)

    print("2. FilFinder Loose   (cache_loose/):")
    m_fil_loose  = load_mask_from_dir(dir_fil_loose, shape)

    print("3. DisPerSE Strict   (full/):")
    m_dis_strict = load_mask_from_dir(dir_dis_strict, shape)

    print("4. DisPerSE Loose    (full_loose/):")
    m_dis_loose  = load_mask_from_dir(dir_dis_loose, shape)

    # --- footprint境界のノイズ除去(全マスクに適用) ---
    m_fil_strict = m_fil_strict & footprint_eroded
    m_fil_loose  = m_fil_loose & footprint_eroded
    m_dis_strict = m_dis_strict & footprint_eroded
    m_dis_loose  = m_dis_loose & footprint_eroded

    # --- DisPerSEの短いノイズ的セグメントを除外 ---
    print(f"\n[DisPerSE 長さフィルタ] {MIN_FILAMENT_LENGTH_PC}pc未満のセグメントを除外中...")
    m_dis_strict = filter_short_segments(m_dis_strict, pc_per_pixel, label_name="dis_strict")
    m_dis_loose  = filter_short_segments(m_dis_loose, pc_per_pixel, label_name="dis_loose")

    # Step 1: 各1px骨格線を 0.1pc に膨張
    fil_strict_thick = dilate_to_width(m_fil_strict, pc_per_pixel)
    dis_strict_thick = dilate_to_width(m_dis_strict, pc_per_pixel)
    fil_loose_thick  = dilate_to_width(m_fil_loose, pc_per_pixel)
    dis_loose_thick  = dilate_to_width(m_dis_loose, pc_per_pixel)

    # Step 2: Class 1 / Class 2 判定
    class1_thick = fil_strict_thick | dis_strict_thick
    loose_thick_or = fil_loose_thick | dis_loose_thick
    class2_thick   = loose_thick_or & ~class1_thick

    # Step 3: ノイズ除去
    if MIN_BLOB_AREA_PX > 0:
        class1_thick = remove_small_blobs(class1_thick, min_size_px=MIN_BLOB_AREA_PX)
        class2_thick = remove_small_blobs(class2_thick, min_size_px=MIN_BLOB_AREA_PX)

    # Step 4: 3クラスラベル割り当て
    label_map = np.zeros(shape, dtype=np.uint8)
    label_map[class1_thick] = 1  # Class 1: 確定コア
    label_map[class2_thick] = 2  # Class 2: 候補

    # FITS保存
    os.makedirs(OUT_DIR, exist_ok=True)
    fits_out_path = os.path.join(OUT_DIR, f"{FITS_NAME}_label3class.fits")
    hdu = fits.PrimaryHDU(data=label_map, header=header)
    hdu.writeto(fits_out_path, overwrite=True)

    print(f"-> ラベルFITS出力完了: {fits_out_path}")
    return label_map


# ============================================================
# 3. 可視化処理 (全体合成マップ描画)
# ============================================================
def visualize_full_map(img_raw, label_map, fits_path):
    """元画像の上に Class 1 (確定コア: 赤色) と Class 2 (候補: 水色) を描画して保存"""
    print("\n可視化マップ (PNG) の描画・保存処理を開始...")
    ny, nx = label_map.shape

    # --- 【修正1】 aquila 以外の領域で画像端（外枠）を NaN 扱いにする ---
    data_display = img_raw.copy()
    if "aquila" not in FITS_NAME.lower():
        print(f"  [画像端処理] 'aquila' 以外の領域のため、画像枠線 {BORDER_NAN_PX}px を NaN 化します。")
        data_display[:BORDER_NAN_PX, :] = np.nan
        data_display[-BORDER_NAN_PX:, :] = np.nan
        data_display[:, :BORDER_NAN_PX] = np.nan
        data_display[:, -BORDER_NAN_PX:] = np.nan

    fig, ax = plt.subplots(figsize=(12, 12))

    # --- 背景数値計算（実数表示用） ---
    finite_mask = np.isfinite(data_display) & (data_display > 0)
    data_masked = np.ma.masked_invalid(data_display)

    if np.any(finite_mask):
        vmin, vmax = np.percentile(data_display[finite_mask], [0.1, 99.5])
    else:
        vmin, vmax = 1e21, 1e23  # フォールバック値

    cmap_obj = plt.get_cmap("inferno").copy()
    cmap_obj.set_bad("black")

    im = ax.imshow(data_masked, cmap=cmap_obj, origin="lower", vmin=vmin, vmax=vmax)

    overlay = np.zeros((ny, nx, 4), dtype=float)
    overlay[label_map == 2] = COLOR_CAND_CYAN  # 水色 (Class 2: 候補)
    overlay[label_map == 1] = COLOR_CORE_RED   # 赤色 (Class 1: 確定コア)

    ax.imshow(overlay, origin="lower")

    # --- 【修正3】 カラーバーの実数表示 (指数表記 1e21 / 1e22等) ---
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.formatter = ticker.ScalarFormatter(useMathText=True)
    cbar.formatter.set_powerlimits((0, 0))  # 常に10の累乗表示 (1eXX) に強制
    cbar.update_ticks()
    cbar.set_label(r"Column Density $N(H_2)$ [cm$^{-2}$]")

    # --- 【修正2】 スケールバー描画 (5pc) ---
    distance_pc = guess_distance_pc(FITS_NAME)
    _, pixel_scale_arcsec = load_wcs_and_pixelscale(fits_path)
    draw_scalebar(ax, nx, ny, pixel_scale_arcsec, distance_pc, scalebar_pc=5.0, color="white")

    legend_elements = [
        Patch(facecolor=COLOR_CORE_RED[:3], edgecolor='none', alpha=0.9, label='Class 1: Confirmed Core (Strict OR / Red)'),
        Patch(facecolor=COLOR_CAND_CYAN[:3], edgecolor='none', alpha=0.7, label='Class 2: Candidate (Loose OR - Class 1 / Cyan)')
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=11, facecolor='#1e293b', edgecolor='none', labelcolor='white')

    ax.set_title(f"{FITS_NAME} — 3-Class Label Map Overlay", fontsize=14)
    ax.set_xlabel("x [pixels]")
    ax.set_ylabel("y [pixels]")

    out_map_path = os.path.join(OUT_DIR, f"{FITS_NAME}_3class_full_map.png")
    fig.tight_layout()
    fig.savefig(out_map_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"-> 可視化マップ (PNG) 保存完了: {out_map_path}\n")


# ============================================================
# メイン実行部
# ============================================================
def main():
    fits_path = os.path.join(ROOT_DIR, FITS_NAME + ".fits")
    if not os.path.exists(fits_path):
        fits_path = os.path.join(REGION_DIR, FITS_NAME + ".fits")

    if not os.path.exists(fits_path):
        print(f"[エラー] 元FITSファイルが見つかりません: {fits_path}")
        return

    print(f"元FITS読み込み: {fits_path}")
    with fits.open(fits_path) as hdul:
        header  = hdul[0].header.copy()
        img_raw = np.squeeze(np.asarray(hdul[0].data, dtype=np.float64))

    label_map = build_3class_label(img_raw, header, fits_path)
    visualize_full_map(img_raw, label_map, fits_path)

    print("==================================================")
    print("可視化処理がすべて完了しました！")
    print(f"出力先フォルダ: {OUT_DIR}")
    print("==================================================")


if __name__ == "__main__":
    main()
