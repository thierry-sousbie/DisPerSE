import os
import time
import pickle
import numpy as np
import cv2
from scipy.ndimage import binary_dilation, label
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
#     + 右下にpcスケールバーを追加(RGB画像として出力)
#   - ファイル名にタイムスタンプを付与し、ビューア側のキャッシュ表示を回避
#   - 定量比較(IoU/Precision/Recall/F1)・CSV保存は行わない(可視化のみ)
#
# 【前回からの修正点】
#   1. proj_plane_pixel_scales の import 漏れを追加
#   2. 元スクリプト後半にあった save_science_figure の呼び出しブロックは、
#      img_raw / skeleton / mask / w / distance_pc / display_name /
#      beam_arcsec / filaments / vis_dir / NSIG のいずれもこのスクリプト内で
#      定義されておらず、実行するとNameErrorになるコードだったため削除しました。
#      (別のDisPerSE実行スクリプト側のコードが混入したものと思われます。
#       もしそちらの科学図保存が必要であれば、そのスクリプトを見せてください)
#   3. save_comparison_overlay にスケールバー描画を統合し、RGBオーバーレイ画像
#      として出力するようにしました。
# ============================================================

root_dir = "/mnt/sp4t/fujiwara.r"
fits_name = "HGBS_ic5146_hires_column_density_map"  # 比較したい領域に合わせて変更

# --- 領域名(fits_nameに含まれる文字列)と距離[pc]の対応 -----------------
# 距離は現状の暫定値。Aquila/IC5146は会議で議論継続中のため要更新。
DISTANCE_PC = {
    "aquila": 440,
    "ic5146": 6000,
    "orionb": 400,
    "polaris": 150,
    "taurus": 140,
}

def guess_distance_pc(fits_name, distance_table=DISTANCE_PC):
    name_lower = fits_name.lower()
    for key, dist in distance_table.items():
        if key in name_lower:
            return dist
    print(f"[警告] fits_name='{fits_name}' から距離を推定できませんでした。"
          f"スケールバーは描画されません。")
    return None


# --- WCS取得用のFITSパス ------------------------------------------------
# root_dir直下に {fits_name}.fits として置かれている想定。
# 例: /mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map.fits
original_fits_path = os.path.join(root_dir, fits_name + ".fits")


def load_column_density(fits_path):
    """柱密度FITSのデータ本体(2次元配列)を読み込む。科学図の背景描画に使う。"""
    if not os.path.exists(fits_path):
        raise FileNotFoundError(
            f"柱密度FITSが見つかりません: {fits_path}\n"
            f"original_fits_path を実際のファイル配置に合わせて修正してください。"
        )
    return fits.getdata(fits_path)


def load_wcs_and_pixelscale(fits_path):
    """FITSヘッダからWCSとpixelスケール[arcsec/pix]を取得する。"""
    if not os.path.exists(fits_path):
        print(f"[警告] FITSファイルが見つかりません: {fits_path}\n"
              f"       スケールバーは描画されません。")
        return None, None
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
# 正規化済み入力画像(背景表示用)
# ------------------------------------------------------------
input_u8 = cv2.imread(
    os.path.join(full_dir, fits_name + "_input_full_ml.png"),
    cv2.IMREAD_GRAYSCALE
)
if input_u8 is None:
    raise FileNotFoundError("入力画像が読み込めません。パスを確認してください。")


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
print("disperse_skeleton の True画素数:", disperse_skeleton.sum())
print("disperse_mask     の True画素数:", disperse_mask.sum())
print("filfinder_skeleton の True画素数:", filfinder_skeleton.sum())
print("filfinder_mask     の True画素数:", filfinder_mask.sum())

print("input_u8 shape:", input_u8.shape)
print("disperse_skeleton shape:", disperse_skeleton.shape)
print("filfinder_skeleton shape:", filfinder_skeleton.shape)

assert disperse_skeleton.shape == filfinder_skeleton.shape, "shapeが一致しません。"
assert disperse_mask.shape == filfinder_mask.shape, "shapeが一致しません(mask)。"


# ============================================================
# 色分けオーバーレイ可視化(RGB出力 + 右下スケールバー)
#   黄 = 両方一致(TP), 赤 = DisPerSEのみ(FP), 青 = FilFinderのみ(FN)
# ============================================================
def draw_scalebar(ax, nx, ny, pixel_scale_arcsec, distance_pc,
                   scalebar_pc=5.0, color="white"):
    """
    右下にLPc棒(横棒＋テキスト)を描画する。
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


def save_comparison_overlay(background_u8, pred, gt, out_path, title,
                             vis_dilate_iter=0, bg_brightness=0.6, overlay_alpha=0.85,
                             pixel_scale_arcsec=None, distance_pc=None,
                             scalebar_pc=5.0, bg_cmap="inferno"):
    """
    bg_brightness : 背景の明るさ倍率(0〜1)。1.0で元の明るさのまま、0.3のようにすると暗くなる
    overlay_alpha : オーバーレイ色の不透明度(0〜1)。1.0だと背景が完全に隠れる、
                    0.85くらいだと背景がうっすら透けて見える
    pixel_scale_arcsec, distance_pc : スケールバー描画に使用。どちらかがNoneなら描画スキップ。
    bg_cmap       : 背景のカラーマップ。科学図(save_science_figure)と揃えて"inferno"を既定値に。
    """
    ny, nx = pred.shape
    fig, ax = plt.subplots(figsize=(16, 13))

    dimmed_background = (background_u8.astype(float) * bg_brightness).astype(np.uint8)
    im_bg = ax.imshow(dimmed_background, cmap=bg_cmap, origin="lower", vmin=0, vmax=255)

    if vis_dilate_iter > 0:
        pred_vis = binary_dilation(pred, iterations=vis_dilate_iter)
        gt_vis = binary_dilation(gt, iterations=vis_dilate_iter)
    else:
        pred_vis, gt_vis = pred, gt

    both = pred_vis | gt_vis
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
    cbar.set_label("Normalized intensity (background)")

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
                    skeleton(1px芯線)を渡せば細線オーバーレイ、mask(太らせ済み)を
                    渡せば太い領域オーバーレイになる。
    filament_color: matplotlibが解釈できる色名 (例: "red", "cyan")
    filament_alpha: 重ね描きの不透明度 (skeletonは1.0推奨、maskは0.6程度推奨)
    n_filaments   : タイトルに表示するフィラメント数。
                    NOTE: ここではdisperse_skeletonの連結成分数を簡易カウントした
                    ものを渡している。DisPerSEのarc-assembly(隣接アークの接続)を
                    経た「科学的に意味のあるフィラメント本数」とは異なり、
                    過大にカウントされやすい点に注意。正式な本数が別途あれば
                    そちらに置き換えること。
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

    # フィラメント(skeletonまたはmask)を指定色で重ね描き
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


# --- スケールバー用の情報を取得 -----------------------------------------
distance_pc = guess_distance_pc(fits_name)
wcs_obj, pixel_scale_arcsec = load_wcs_and_pixelscale(original_fits_path)

# input_u8とFITSのpixel数が一致しているか確認(ダウンサンプル等で異なる場合は
# スケールバーがずれるため、一致していない場合は警告のうえスケールバーなしで進める)
if wcs_obj is not None:
    fits_shape = fits.getdata(original_fits_path).shape
    if fits_shape != input_u8.shape:
        print(f"[警告] FITSのshape{fits_shape}とinput_u8のshape{input_u8.shape}が"
              f"一致していません。ダウンサンプル等で解像度が変わっている可能性があります。"
              f"スケールバーが不正確になる場合は pixel_scale_arcsec を"
              f"手動で補正してください。")

# タイムスタンプを付けて、ビューア側の古いキャッシュ表示を確実に回避する
ts = int(time.time())

save_comparison_overlay(
    input_u8, disperse_skeleton, filfinder_skeleton,
    os.path.join(compare_dir, f"{fits_name}_compare_skeleton_{ts}.png"),
    "Yellow=Both / Red=DisPerSE only / Blue=FilFinder only (skeleton)",
    vis_dilate_iter=8, bg_brightness=0.6, overlay_alpha=0.9,
    pixel_scale_arcsec=pixel_scale_arcsec, distance_pc=distance_pc, scalebar_pc=5.0,
)

save_comparison_overlay(
    input_u8, disperse_mask, filfinder_mask,
    os.path.join(compare_dir, f"{fits_name}_compare_mask_{ts}.png"),
    "Yellow=Both / Red=DisPerSE only / Blue=FilFinder only (mask)",
    vis_dilate_iter=0, bg_brightness=0.6, overlay_alpha=0.7,   # maskは面積が広いので少し透過を強めに
    pixel_scale_arcsec=pixel_scale_arcsec, distance_pc=distance_pc, scalebar_pc=5.0,
)

print("DONE")


# ============================================================
# 科学図の生成(DisPerSE結果を使用)
# ============================================================
if wcs_obj is None:
    raise RuntimeError(
        f"WCSが取得できなかったため科学図を生成できません。"
        f"original_fits_path='{original_fits_path}' を確認してください。"
    )
if distance_pc is None:
    raise RuntimeError(
        f"距離が推定できなかったため科学図を生成できません。"
        f"DISTANCE_PC辞書にfits_nameに対応するキーがあるか確認してください。"
    )

img_raw = load_column_density(original_fits_path)
w = wcs_obj

if img_raw.shape != disperse_skeleton.shape:
    raise RuntimeError(
        f"img_rawのshape{img_raw.shape}とdisperse_skeletonのshape"
        f"{disperse_skeleton.shape}が一致しません。ダウンサンプル等で"
        f"解像度が変わっている場合は、どちらかを揃えてから実行してください。"
    )

# フィラメント本数の簡易カウント(連結成分数)。上記docstring注記のとおり簡易値。
_, n_filaments = label(disperse_skeleton)

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

# 科学図(2): 背景の上にDisPerSE mask(赤い太い半透明領域)を重ねたもの
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
