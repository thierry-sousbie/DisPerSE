import os
import time
import pickle
import numpy as np
import cv2
from scipy.ndimage import binary_dilation, label, find_objects
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

# ============================================================
# DisPerSE と FilFinder の結果を比較するスクリプト(修正版)
#   - skeleton(1px)同士、mask(太い領域)同士をそれぞれ可視化比較
#   - 色分けオーバーレイ(黄=両方一致、赤=DisPerSEのみ、青=FilFinderのみ)
#     + 背景は柱密度の物理値(img_raw, infernoカラーマップ)
#     + 右下にpcスケールバーを追加(RGB画像として出力)
#   - 0.3pc未満の短いセグメントは比較・可視化前に除外(フィラメント長カットオフ)
#   - ファイル名にタイムスタンプを付与し、ビューア側のキャッシュ表示を回避
#   - 定量比較(IoU/Precision/Recall/F1)・CSV保存は行わない(可視化のみ)
#
# 【今回の修正点】
#   1. save_comparison_overlay の背景を、正規化済みPNG(input_u8)から
#      柱密度の物理値(img_raw)に変更。カラーバーも cm^-2 単位表記に統一。
#   2. 0.3pc未満の連結成分(短いセグメント)を、比較・可視化の前に除外する
#      filter_short_segments() を追加。disperse/filfinderのskeleton・mask
#      それぞれに適用。
#   3. 上記に伴い、img_raw / WCS / 距離の取得を先頭付近に移動し、
#      これらが取得できない場合は明示的にエラーで止まるようにした
#      (背景描画とカットオフ計算の両方に必須のため)。
# ============================================================

root_dir = "/mnt/sp4t/fujiwara.r"
fits_name = "HGBS_polaris_hires_column_density_map"  # 比較したい領域に合わせて変更

# --- 領域名(fits_nameに含まれる文字列)と距離[pc]の対応 -----------------
# 距離は現状の暫定値。Aquila/IC5146は会議で議論継続中のため要更新。
DISTANCE_PC = {
    "aquila": 440,
    "ic5146": 600,
    "orionb": 400,
    "polaris": 150,
    "taun3": 140,   # Taurusの実際のfits_nameは "HGBS_tauN3_..." のため別名として追加
}

# フィラメントとみなす最小長さ[pc]。これ未満の連結成分は比較・可視化前に除外する。
MIN_FILAMENT_LENGTH_PC = 0.3

# フィラメント幅[pc]。DisPerSE側マスク生成スクリプトと同じ値に揃える。
# skeleton比較図の可視化用の膨張(dilation)幅の計算にも使用する。
WIDTH_PC = 0.1


def guess_distance_pc(fits_name, distance_table=DISTANCE_PC):
    name_lower = fits_name.lower()
    for key, dist in distance_table.items():
        if key in name_lower:
            return dist
    print(f"[警告] fits_name='{fits_name}' から距離を推定できませんでした。")
    return None


# --- WCS取得用のFITSパス ------------------------------------------------
# root_dir直下に {fits_name}.fits として置かれている想定。
# 例: /mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map.fits
original_fits_path = os.path.join(root_dir, fits_name + ".fits")


def load_column_density(fits_path):
    """柱密度FITSのデータ本体(2次元配列, 物理値cm^-2)を読み込む。"""
    if not os.path.exists(fits_path):
        raise FileNotFoundError(
            f"柱密度FITSが見つかりません: {fits_path}\n"
            f"original_fits_path を実際のファイル配置に合わせて修正してください。"
        )
    return fits.getdata(fits_path)


def load_wcs_and_pixelscale(fits_path):
    """FITSヘッダからWCSとpixelスケール[arcsec/pix]を取得する。"""
    if not os.path.exists(fits_path):
        raise FileNotFoundError(
            f"FITSファイルが見つかりません: {fits_path}\n"
            f"original_fits_path を確認してください。"
        )
    header = fits.getheader(fits_path)
    w = WCS(header)
    pixel_scale_deg = proj_plane_pixel_scales(w)[:2].mean()
    pixel_scale_arcsec = pixel_scale_deg * 3600.0
    return w, pixel_scale_arcsec


disperse_dir = os.path.join(root_dir, fits_name, "disperse")
cache_dir = os.path.join(root_dir, fits_name, "cache")
full_dir = os.path.join(root_dir, fits_name, "full")
compare_dir = os.path.join(disperse_dir, "compare")
vis_dir = os.path.join(disperse_dir, "vis")
os.makedirs(compare_dir, exist_ok=True)
os.makedirs(vis_dir, exist_ok=True)

# --- 科学図(save_science_figure)用の設定 --------------------------------
# NSIG: DisPerSE実行時に使ったnumber-of-sigma(パラメータ探索でnsig=1に決定済み)
NSIG = 1
# BEAM_ARCSEC: HGBS柱密度マップのビームサイズ。領域ごとに異なる場合は要確認・修正。
BEAM_ARCSEC = 18.2
# display_name: 図タイトルに使う表示名。fits_nameから領域名を抜き出して整形。
display_name = fits_name.replace("HGBS_", "").replace("_hires_column_density_map", "")

print("fits_name:", fits_name)
print("disperse_dir:", disperse_dir)


# ------------------------------------------------------------
# DisPerSE側の結果を読み込む(0/1のラベル画像)
# ファイル名は実際にfull/に保存されているものに合わせてください
# ------------------------------------------------------------
disperse_skeleton_path = os.path.join(full_dir, fits_name + "_skeleton_01.png")
disperse_mask_path = os.path.join(full_dir, fits_name + "_mask_01.png")

print("full/ の中身:", os.listdir(os.path.join(full_dir)))

disperse_skeleton_raw = cv2.imread(disperse_skeleton_path, cv2.IMREAD_GRAYSCALE)
disperse_mask_raw = cv2.imread(disperse_mask_path, cv2.IMREAD_GRAYSCALE)

if disperse_skeleton_raw is None:
    raise FileNotFoundError(f"見つかりません: {disperse_skeleton_path}")
if disperse_mask_raw is None:
    raise FileNotFoundError(f"見つかりません: {disperse_mask_path}")

disperse_skeleton = disperse_skeleton_raw.astype(bool)
disperse_mask = disperse_mask_raw.astype(bool)


# ------------------------------------------------------------
# FilFinder側の結果を読み込む(pickleキャッシュから)
# ------------------------------------------------------------
skel_cache_path = os.path.join(cache_dir, fits_name + "_FilFinder2D-medskel.pkl")
mask_cache_path = os.path.join(cache_dir, fits_name + "_FilFinder2D-create_mask.pkl")

with open(skel_cache_path, "rb") as f:
    fil_ff_skel_obj = pickle.load(f)
filfinder_skeleton = fil_ff_skel_obj.skeleton.astype(bool)

with open(mask_cache_path, "rb") as f:
    fil_ff_mask_obj = pickle.load(f)
filfinder_mask = fil_ff_mask_obj.mask.astype(bool)


# ------------------------------------------------------------
# 診断: 画素数とshapeを確認
# ------------------------------------------------------------
print("disperse_skeleton の True画素数(フィルタ前):", disperse_skeleton.sum())
print("disperse_mask     の True画素数(フィルタ前):", disperse_mask.sum())
print("filfinder_skeleton の True画素数(フィルタ前):", filfinder_skeleton.sum())
print("filfinder_mask     の True画素数(フィルタ前):", filfinder_mask.sum())

print("disperse_skeleton shape:", disperse_skeleton.shape)
print("filfinder_skeleton shape:", filfinder_skeleton.shape)

assert disperse_skeleton.shape == filfinder_skeleton.shape, "shapeが一致しません。"
assert disperse_mask.shape == filfinder_mask.shape, "shapeが一致しません(mask)。"


# ------------------------------------------------------------
# 柱密度データ(背景・スケールバー・カットオフ計算に使用)
# ------------------------------------------------------------
distance_pc = guess_distance_pc(fits_name)
if distance_pc is None:
    raise RuntimeError(
        f"距離が推定できませんでした。DISTANCE_PC辞書にfits_nameに対応する"
        f"キーがあるか確認してください。"
    )

wcs_obj, pixel_scale_arcsec = load_wcs_and_pixelscale(original_fits_path)
img_raw = load_column_density(original_fits_path)

if img_raw.shape != disperse_skeleton.shape:
    raise RuntimeError(
        f"img_rawのshape{img_raw.shape}とdisperse_skeletonのshape"
        f"{disperse_skeleton.shape}が一致しません。ダウンサンプル等で"
        f"解像度が変わっている場合は、どちらかを揃えてから実行してください。"
    )

pc_per_pixel = (pixel_scale_arcsec / 206265.0) * distance_pc
print(f"pixel_scale_arcsec={pixel_scale_arcsec:.3f}\", "
      f"distance_pc={distance_pc}, pc_per_pixel={pc_per_pixel:.5f}")

# WIDTH_PC(フィラメント幅)を、skeleton比較図の可視化用膨張(dilation)の
# 半径ピクセル数に変換する。DisPerSE側マスク生成スクリプトと同じWIDTH_PCを
# 使うことで、可視化上の線の太さが実際のフィラメント幅定義と一致する。
vis_dilate_px = max(1, int(round((WIDTH_PC / pc_per_pixel) / 2)))
print(f"WIDTH_PC={WIDTH_PC}pc -> vis_dilate_px={vis_dilate_px}px "
      f"(skeleton比較の可視化用膨張半径)")


# ============================================================
# フィラメント長カットオフ(0.3pc未満の連結成分を除外)
# ============================================================
def filter_short_segments(mask, pc_per_pixel, min_length_pc=MIN_FILAMENT_LENGTH_PC,
                           label_name=""):
    """
    2値マスクを連結成分に分解し、各成分のバウンディングボックス対角長を
    物理長[pc]に換算して、min_length_pc未満の成分を除去する。

    NOTE: 対角長は「連結成分の広がり」の近似値であり、曲がりくねった
    フィラメントの実際の経路長(弧長)そのものではない点に注意。
    より厳密な弧長が必要な場合は、スケルトングラフ上の経路長を
    計算する実装に置き換えること。

    NOTE2: 連結成分判定は8方向連結(対角も接続とみなす)を使用する。
    1px幅のスケルトンは斜めにしかつながっていない箇所が多く、
    デフォルトの4方向連結(scipy.ndimage.labelの既定)だと、本来ひと続きの
    フィラメントが細かく分断され、それぞれが短すぎるとして誤って
    大量除去されてしまうため。
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


print("\n--- フィラメント長カットオフ適用中 ---")
disperse_skeleton = filter_short_segments(
    disperse_skeleton, pc_per_pixel, label_name="disperse_skeleton")
disperse_mask = filter_short_segments(
    disperse_mask, pc_per_pixel, label_name="disperse_mask")
filfinder_skeleton = filter_short_segments(
    filfinder_skeleton, pc_per_pixel, label_name="filfinder_skeleton")
filfinder_mask = filter_short_segments(
    filfinder_mask, pc_per_pixel, label_name="filfinder_mask")

print("disperse_skeleton の True画素数(フィルタ後):", disperse_skeleton.sum())
print("disperse_mask     の True画素数(フィルタ後):", disperse_mask.sum())
print("filfinder_skeleton の True画素数(フィルタ後):", filfinder_skeleton.sum())
print("filfinder_mask     の True画素数(フィルタ後):", filfinder_mask.sum())


# ============================================================
# 色分けオーバーレイ可視化(RGB出力 + 右下スケールバー)
#   黄 = 両方一致(TP), 赤 = DisPerSEのみ(FP), 青 = FilFinderのみ(FN)
# ============================================================
def draw_scalebar(ax, nx, ny, pixel_scale_arcsec, distance_pc,
                   scalebar_pc=5.0, color="white"):
    """
    右下にpcスケールバー(横棒＋テキスト)を描画する。
    pixel_scale_arcsec または distance_pc が None の場合は何もしない。
    """
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


def save_comparison_overlay(background_data, pred, gt, out_path, title,
                             vis_dilate_iter=0, bg_brightness=0.6, overlay_alpha=0.85,
                             pixel_scale_arcsec=None, distance_pc=None,
                             scalebar_pc=5.0, bg_cmap="inferno", bg_vmax=None):
    """
    background_data : 柱密度の物理値(cm^-2)の2次元配列(img_raw)。
    bg_brightness    : 背景の明るさ倍率(0〜1)。vmaxは固定したまま、表示データを
                        スケールダウンすることで暗く見せる(オーバーレイ色を目立たせる)。
    overlay_alpha    : オーバーレイ色の不透明度(0〜1)。
    pixel_scale_arcsec, distance_pc : スケールバー描画に使用。どちらかがNoneなら描画スキップ。
    bg_cmap          : 背景のカラーマップ。科学図(save_science_figure)と揃えて"inferno"を既定値に。
    bg_vmax          : Noneなら背景データの99.5パーセンタイルを自動使用。
    """
    ny, nx = pred.shape
    fig, ax = plt.subplots(figsize=(16, 13))

    finite_vals = background_data[np.isfinite(background_data)]
    if bg_vmax is None:
        bg_vmax = np.nanpercentile(finite_vals, 99.5) if finite_vals.size > 0 else 1.0

    dimmed_background = np.nan_to_num(background_data, nan=0.0) * bg_brightness
    im_bg = ax.imshow(dimmed_background, cmap=bg_cmap, origin="lower",
                       vmin=0, vmax=bg_vmax)

    if vis_dilate_iter > 0:
        pred_vis = binary_dilation(pred, iterations=vis_dilate_iter)
        gt_vis = binary_dilation(gt, iterations=vis_dilate_iter)
    else:
        pred_vis, gt_vis = pred, gt

    both = pred_vis & gt_vis
    disperse_only = pred_vis & ~gt_vis
    filfinder_only = ~pred_vis & gt_vis

    # RGBA float配列(0-1)として明示的に構築 -> matplotlibでRGB画像として描画される
    overlay = np.zeros((*pred_vis.shape, 4), dtype=float)

    yellow = mcolors.to_rgb("yellow")
    red = mcolors.to_rgb("red")
    blue = mcolors.to_rgb("dodgerblue")

    for mask_, color in [(both, yellow), (disperse_only, red), (filfinder_only, blue)]:
        overlay[mask_, 0] = color[0]
        overlay[mask_, 1] = color[1]
        overlay[mask_, 2] = color[2]
        overlay[mask_, 3] = overlay_alpha

    ax.imshow(overlay, origin="lower")
    ax.set_title(title, fontsize=16)
    ax.set_xlabel("x [pixels]")
    ax.set_ylabel("y [pixels]")

    cbar = fig.colorbar(im_bg, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"Column Density $N(H_2)$ [cm$^{-2}$]")

    draw_scalebar(ax, nx, ny, pixel_scale_arcsec, distance_pc, scalebar_pc=scalebar_pc)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  [{os.path.basename(out_path)}] 黄={both.sum()}px, "
          f"赤={disperse_only.sum()}px, 青={filfinder_only.sum()}px")


# ============================================================
# 科学図(論文体裁): inferno colormap + colorbar + タイトル + pixel軸 +
# フィラメント重ね描き + 右下スケールバー
# ============================================================
def save_science_figure(data2d, overlay_mask, out_path, wcs, distance_pc,
                         display_name, beam_arcsec, n_filaments,
                         scalebar_pc=5.0, vmax=None, cmap="inferno",
                         filament_color="red", filament_alpha=1.0):
    """
    data2d        : 物理値(線形, cm^-2)の柱密度2次元配列。NaNは黒表示になる。
    overlay_mask  : data2dと同形状の bool配列。Trueの画素に色を重ねて描く。
    filament_color: matplotlibが解釈できる色名 (例: "red", "cyan")
    filament_alpha: 重ね描きの不透明度 (skeletonは1.0推奨、maskは0.6程度推奨)
    n_filaments   : タイトルに表示するフィラメント数(0.3pcカットオフ後の連結成分数)。
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

    color_rgb = mcolors.to_rgb(filament_color)
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
    pixel_scale_arcsec_local = pixel_scale_deg * 3600.0
    draw_scalebar(ax, nx, ny, pixel_scale_arcsec_local, distance_pc, scalebar_pc=scalebar_pc)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"科学図を保存しました: {out_path}")


# ============================================================
# 比較オーバーレイの出力(カットオフ適用後のデータを使用)
# ============================================================
ts = int(time.time())

save_comparison_overlay(
    img_raw, disperse_skeleton, filfinder_skeleton,
    os.path.join(compare_dir, f"{fits_name}_compare_skeleton_{ts}.png"),
    "Yellow=Both / Red=DisPerSE only / Blue=FilFinder only (skeleton)",
    vis_dilate_iter=vis_dilate_px, bg_brightness=0.6, overlay_alpha=0.9,
    pixel_scale_arcsec=pixel_scale_arcsec, distance_pc=distance_pc, scalebar_pc=5.0,
)

save_comparison_overlay(
    img_raw, disperse_mask, filfinder_mask,
    os.path.join(compare_dir, f"{fits_name}_compare_mask_{ts}.png"),
    "Yellow=Both / Red=DisPerSE only / Blue=FilFinder only (mask)",
    vis_dilate_iter=0, bg_brightness=0.6, overlay_alpha=0.7,   # maskは面積が広いので少し透過を強めに
    pixel_scale_arcsec=pixel_scale_arcsec, distance_pc=distance_pc, scalebar_pc=5.0,
)

print("DONE (comparison overlays)")


# ============================================================
# 科学図の生成(DisPerSE結果を使用、カットオフ適用後)
# ============================================================
w = wcs_obj

# フィラメント本数(0.3pcカットオフ後の連結成分数)
# フィラメント本数(0.3pcカットオフ後の連結成分数)。filter_short_segmentsと
# 同じ8方向連結で揃える。
_, n_filaments = label(disperse_skeleton, structure=np.ones((3, 3), dtype=int))

# 科学図(1): 背景の上にDisPerSE skeleton(赤い細線)を重ねたもの
save_science_figure(
    data2d=img_raw,
    overlay_mask=disperse_skeleton,
    out_path=os.path.join(vis_dir, f"{fits_name}_science_figure_skeleton_nsig{NSIG}.png"),
    wcs=w,
    distance_pc=distance_pc,
    display_name=display_name,
    beam_arcsec=BEAM_ARCSEC,
    n_filaments=n_filaments,
    filament_color="red",
    filament_alpha=1.0,
)

# 科学図(2): 背景の上にDisPerSE mask(赤の太い半透明領域)を重ねたもの
save_science_figure(
    data2d=img_raw,
    overlay_mask=disperse_mask,
    out_path=os.path.join(vis_dir, f"{fits_name}_science_figure_mask_nsig{NSIG}.png"),
    wcs=w,
    distance_pc=distance_pc,
    display_name=display_name,
    beam_arcsec=BEAM_ARCSEC,
    n_filaments=n_filaments,
    filament_color="red",
    filament_alpha=0.6,
)

print("DONE (science figures)")
