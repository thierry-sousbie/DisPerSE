#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DisPerSE 実行 + FITS 処理・可視化スクリプト (Orion B対応・完全NaN排除＆NxN縮小・最終統合版)
- 画像が巨大かつ複雑な場合(Orion B等)にNxN縮小を実行(DOWNSAMPLE_FACTORで指定)
- バッファオーバーフローを確実に回避するため、縮小後の最終入力データから NaN/Inf を完全に排除(0.0置換)
- 境界をなだらかに平滑化してトポロジーの複雑さを極小化
- トポロジー解析後にスケールを4倍に正確に復元
- タイル切り出し(in, out)を完全廃止
- 端の赤枠ノイズを cv2.erode で完全除去
- 背景に全体のフル解像度グレースケール画像を採用
- 可視化の向きを FilFinder と同じ origin="lower" に完全統一

[修正 2026-07-17 その1]
- mse実行前のAv4未満→0.0マスキングを廃止（冗長のため）。
  → NaNのみ0.0に置換し、Av4カットはこれまで通りvalid_mask_erodedで事後適用する。

[修正 2026-07-17 その2]
- 上記修正後も同一箇所(Pass 1: building)でbuffer overflowが再現。ただし臨界点数は
  1616→35258と20倍以上変化したにもかかわらずクラッシュ箇所・症状は不変だった。
- 一方で、総セル数(31,757,184)・1-cells数(15,878,592)などグリッド解像度のみに
  依存する値は2回の実行で完全に一致しており、これはデータ内容とは無関係。
  → 原因はデータの中身ではなく、ダウンサンプル後もグリッドの総ピクセル数
  (2042x2592 ≈ 529万px)がこのmseビルドの内部バッファ容量を超えていることだと推定。
  → DOWNSAMPLE_FACTOR を4→8に引き上げ、総ピクセル数を約1/4(≈130万px)に削減して再試行する。
"""

import os
import glob
import subprocess
import warnings
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
import matplotlib.patheffects as pe
import matplotlib.colors

# ============================================================
# 設定 (CONFIG)
# ============================================================
root_dir = "/mnt/sp4t/fujiwara.r"

# --- 天体ごとのプリセット ---
# display_name : タイトルに表示する天体名
# fits_name    : FITSファイル名（拡張子なし、root_dir直下に "<fits_name>.fits" がある前提）
# distance_pc  : 距離[pc]（スケールバー・physical幅計算に使用）
# beam_arcsec  : ビームサイズ["]（タイトル表示用。hires列密度マップの典型値18.2"をデフォルトに、
#                実際の値と異なる場合は要調整）
REGIONS = {
    "OrionB":  dict(display_name="OrionB",  fits_name="HGBS_orionB_hires_column_density_map",  distance_pc=400.0, beam_arcsec=18.2),
    "Aquila":  dict(display_name="Aquila",  fits_name="HGBS_aquila_hires_column_density_map",  distance_pc=260.0, beam_arcsec=18.2),
    "IC5146":  dict(display_name="IC5146",  fits_name="HGBS_ic5146_hires_column_density_map",  distance_pc=460.0, beam_arcsec=18.2),
    "Taurus":  dict(display_name="Taurus",  fits_name="HGBS_taurus_hires_column_density_map",  distance_pc=140.0, beam_arcsec=18.2),
    "Polaris": dict(display_name="Polaris", fits_name="HGBS_polaris_hires_column_density_map", distance_pc=150.0, beam_arcsec=18.2),
}

# ここを変更するだけで対象天体を切り替えられる
TARGET_REGION = "OrionB"

_region_cfg = REGIONS[TARGET_REGION]
display_name = _region_cfg["display_name"]
fits_name = _region_cfg["fits_name"]
distance_pc = _region_cfg["distance_pc"]
beam_arcsec = _region_cfg["beam_arcsec"]

fits_path = os.path.join(root_dir, fits_name + ".fits")

# DisPerSE の mse コマンドがあるバイナリディレクトリ
DISPERSE_BIN_DIR = "/home1/fujiwara.r/PycharmProjects/DisPerSE/build/src/mse"

def bin_path(name):
    return os.path.join(DISPERSE_BIN_DIR, name) if DISPERSE_BIN_DIR else name

# DisPerSE 固有パラメータ
NSIG = 1

# マスク幅の方針
WIDTH_MODE = "skeleton_only"
WIDTH_PC = 0.1
# distance_pc は REGIONS[TARGET_REGION] から既に設定済み（上のCONFIG参照）

# Av = 4 相当の柱密度しきい値 (1.88e21 x 2 = 3.76e21 H2 cm-2)
# ※ mse実行前の事前マスキングには使わず、ステップ4のvalid_mask_erodedでのみ使用する
AV4_THRESH = 3.76e+21

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
# ステップ1: FITS読み込み・安全な正規化・完全NaN排除平滑化
# ============================================================
# DisPerSE用の一時的な「NaN排除済みFITSファイル」のパス
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

# RuntimeWarning: divide by zero / invalid value 対策
with np.errstate(divide='ignore', invalid='ignore'):
    img_log_all = np.where(finite_mask_all & (img_raw > 0), np.log10(img_raw), np.nan)

des_log_all = pd.DataFrame(pd.Series(img_log_all[np.isfinite(img_log_all)].ravel()).describe(
    percentiles=[0.005, 0.995]))

data_bak_ml = img_log_all.copy()
data_bak_ml -= des_log_all[0]["0.5%"]

# ゼロ除算を避ける安全ガード
scale_val = (des_log_all[0]["99.5%"] - des_log_all[0]["0.5%"]) / 255
if scale_val == 0:
    scale_val = 1.0
data_bak_ml /= scale_val

data_bak_ml = np.clip(data_bak_ml, 0, 255)
data_bak_ml[~finite_mask_all] = np.nan

# 背景となる「全体グレースケール画像（フルスケール）」
data_bak_ml_cleaned = np.nan_to_num(data_bak_ml, nan=0.0)
input_u8 = np.where(finite_mask_all, data_bak_ml_cleaned, 0).astype("uint8")
cv2.imwrite(os.path.join(full_dir, fits_name + "_input_full_ml.png"), input_u8)
print("全体グレースケール（背景用）入力を保存しました。")


# --- DisPerSE用：NaNのみ0.0に置換したデータの作成 ---
# [修正] Av4未満を0.0にする処理を廃止。NaNのみ0.0に置換する。
# Av4カットはステップ4のvalid_mask_erodedで事後適用されるため、
# ここで事前マスキングすると画像の大部分が厳密に同値(0.0)の巨大な
# 平坦プラトーになり、mseのpersistence pair計算でバッファオーバーフローを
# 引き起こす原因になり得るため。
img_raw_masked = img_raw.copy()
print("NaN領域のみ0.0に置換しています（Av4カットはmse実行後にvalid_mask_erodedで適用）...")
img_raw_masked[~np.isfinite(img_raw_masked)] = 0.0

# 【バッファオーバーフロー回避】巨大・複雑な画像は自動でNxN縮小(DOWNSAMPLE_FACTORで指定)
h_orig, w_orig = img_raw.shape
USE_DOWNSAMPLE = False
DOWNSAMPLE_FACTOR = 8  # 8x8縮小（4x4ではグリッド総px数≈529万でbuffer overflowしたため引き上げ）

if h_orig > 5000 or w_orig > 5000:
    print(f"画像サイズが巨大かつ複雑なため ({w_orig}x{h_orig})、DisPerSEのバッファオーバーフローを確実に回避するために{DOWNSAMPLE_FACTOR}x{DOWNSAMPLE_FACTOR}縮小および平滑化を行います。")
    USE_DOWNSAMPLE = True
    
    # 縮小単位の倍数にクロップ
    h_even = (h_orig // DOWNSAMPLE_FACTOR) * DOWNSAMPLE_FACTOR
    w_even = (w_orig // DOWNSAMPLE_FACTOR) * DOWNSAMPLE_FACTOR
    cropped = img_raw_masked[:h_even, :w_even]
    
    # NxN平均（リバニング、N=DOWNSAMPLE_FACTOR）を計算
    reshaped = cropped.reshape(h_even // DOWNSAMPLE_FACTOR, DOWNSAMPLE_FACTOR, w_even // DOWNSAMPLE_FACTOR, DOWNSAMPLE_FACTOR)
    img_disperse_input = np.mean(reshaped, axis=(1, 3)) # NaNを事前に排除したため通常のnp.meanで安全に計算可能
    
    # --------------------------------------------------------
    # 【追加・強化】縮小した後の画像にNaNが絶対に混入しないよう再度クリーンアップ
    # --------------------------------------------------------
    img_disperse_input[~np.isfinite(img_disperse_input)] = 0.0
    
    # --------------------------------------------------------
    # 【境界バリア除去】トポロジーを劇的に単純化するためのガウシアン平滑化
    # --------------------------------------------------------
    # カーネルサイズ 9x9, 標準偏差 2.0 でしきい値境界の「崖」をなだらかなスロープにします
    img_disperse_input = cv2.GaussianBlur(img_disperse_input, (9, 9), 2.0)
    
    # 平滑化処理中や境界部で奇跡的に発生したNaNも完全にシャットアウト
    img_disperse_input[~np.isfinite(img_disperse_input)] = 0.0
    
    # WCSヘッダーの座標情報もDOWNSAMPLE_FACTOR倍縮小に追従
    hdu_header_disperse = hdu_header.copy()
    if "CDELT1" in hdu_header_disperse: hdu_header_disperse["CDELT1"] *= DOWNSAMPLE_FACTOR
    if "CDELT2" in hdu_header_disperse: hdu_header_disperse["CDELT2"] *= DOWNSAMPLE_FACTOR
    if "CRPIX1" in hdu_header_disperse: hdu_header_disperse["CRPIX1"] = (hdu_header_disperse["CRPIX1"] + (DOWNSAMPLE_FACTOR - 1)/2.0) / float(DOWNSAMPLE_FACTOR)
    if "CRPIX2" in hdu_header_disperse: hdu_header_disperse["CRPIX2"] = (hdu_header_disperse["CRPIX2"] + (DOWNSAMPLE_FACTOR - 1)/2.0) / float(DOWNSAMPLE_FACTOR)
else:
    img_disperse_input = img_raw_masked
    hdu_header_disperse = hdu_header

masked_hdu = fits.PrimaryHDU(data=img_disperse_input, header=hdu_header_disperse)
masked_hdu.writeto(masked_fits_path, overwrite=True)
print(f"DisPerSE用のFITSを保存しました（縮小フラグ={USE_DOWNSAMPLE}、倍率={DOWNSAMPLE_FACTOR}、崖の平滑化＆完全NaN排除済み、Av4事前マスキングなし）: {masked_fits_path}")


# ============================================================
# ステップ2: mse 実行 (縮小＆平滑化・完全NaN排除された安全なFITSを入力にする)
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
# ステップ4: VTK読み込み・ラスタ化 (縮小座標をフルサイズに正確に引き戻す)
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

skeleton = np.zeros((h_orig, w_orig), dtype=np.uint8)

# 白いフィラメントの上に重ねる赤線を太く描くため、thicknessを 3 に設定
for pts in filaments:
    xs, ys = pts[:, 0], pts[:, 1]
    
    # NxN縮小状態で計算されていた場合、元の座標サイズに「DOWNSAMPLE_FACTOR倍に正確に復元」
    if USE_DOWNSAMPLE:
        xs = xs * float(DOWNSAMPLE_FACTOR)
        ys = ys * float(DOWNSAMPLE_FACTOR)
        
    for k in range(len(pts) - 1):
        # [修正] VTKのy座標は img_raw/img_disperse_input と同じ「row0=下端」慣習のまま。
        # 以前の (h_orig - 1 - y) 反転は skeleton 配列だけを img_raw と逆向きにしてしまい、
        # origin="lower" で重ね描きした際に skeleton が上下反転して表示される原因になっていた。
        p1 = (int(round(xs[k])), int(round(ys[k])))
        p2 = (int(round(xs[k + 1])), int(round(ys[k + 1])))
        cv2.line(skeleton, p1, p2, 1, thickness=3)

# 境界の偽フィラメント（赤線ノイズ）を消すためのマスク処理
# (本来のデータで物理的な Av >= 4 だった位置のみを真の領域とする)
# [このステップで初めてAv4カットを適用する]
# [修正] skeletonの反転を削除したため、valid_maskもimg_rawと同じ行慣習のまま使う
# (以前の [::-1, :] は skeleton側の反転に合わせるためのものだったが、その反転自体を廃止した)
original_finite_mask = np.isfinite(img_raw)
original_av4_mask = (img_raw >= AV4_THRESH) & original_finite_mask
valid_mask = original_av4_mask

# エッジ部分のノイズ赤線を強力に除去するため、cv2.erode で判定領域を2ピクセル縮小(Erode)
kernel = np.ones((5, 5), np.uint8)
valid_mask_eroded = cv2.erode(valid_mask.astype(np.uint8), kernel, iterations=2).astype(bool)

# マスク外（元のデータで Av < 4）の領域のスケルトンを完全にシャットアウト
skeleton[~valid_mask_eroded] = 0

cv2.imwrite(os.path.join(full_dir, fits_name + "_skeleton_01.png"), skeleton)
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


# ============================================================
# ステップ6: 可視化 (余分な余白を排除して FilFinder と完全同期)
# ============================================================
def add_scalebar(ax, w_px, h_px, wcs, distance_pc, scalebar_pc=5.0,
                  margin_frac=0.04, color="white"):
    """
    右下にL字型のスケールバー（横棒＋縦棒）を描画する。
    - wcs: 元のフル解像度画像に対応するWCS(astropy.wcs.WCS)
    - distance_pc: 対象天体までの距離[pc]（distance_pc変数を使用）
    - scalebar_pc: バーが表す物理長[pc]（デフォルト5pc）
    横・縦それぞれの軸のピクセルスケールを個別に使うため、
    投影の歪みでx/y方向のpc/pixelが多少異なっていても正しい長さになる。
    """
    pixel_scale_deg = proj_plane_pixel_scales(wcs)[:2]  # [x方向, y方向]
    pixel_scale_rad = (pixel_scale_deg * u.deg).to(u.rad).value
    pc_per_pixel_x, pc_per_pixel_y = distance_pc * pixel_scale_rad
    bar_length_px_x = scalebar_pc / pc_per_pixel_x
    bar_length_px_y = scalebar_pc / pc_per_pixel_y

    margin_x = w_px * margin_frac
    margin_y = h_px * margin_frac

    # L字の角(コーナー)を右下に固定し、そこから横棒は左へ、縦棒は上へ伸ばす
    x_corner = w_px - margin_x
    y_corner = margin_y

    x_start = x_corner - bar_length_px_x
    y_top = y_corner + bar_length_px_y

    # 視認性確保のため、黒縁取り(太め)＋白線(細め)の二重描画
    for lw, col, z in [(6, "black", 5), (3, color, 6)]:
        # 横棒
        ax.plot([x_start, x_corner], [y_corner, y_corner],
                color=col, linewidth=lw, solid_capstyle="butt", zorder=z)
        # 縦棒
        ax.plot([x_corner, x_corner], [y_corner, y_top],
                color=col, linewidth=lw, solid_capstyle="butt", zorder=z)

    fontsize = max(14, h_px * 0.018)
    text_kwargs = dict(color=color, fontsize=fontsize, zorder=6,
                        path_effects=[pe.withStroke(linewidth=3, foreground="black")])

    # 横棒のラベル（棒の下）
    ax.text((x_start + x_corner) / 2, y_corner - h_px * 0.012,
            f"{scalebar_pc:g} pc", ha="center", va="top", **text_kwargs)
    # 縦棒のラベル（棒の右、90度回転）
    ax.text(x_corner + w_px * 0.012, (y_corner + y_top) / 2,
            f"{scalebar_pc:g} pc", ha="left", va="center", rotation=90,
            **text_kwargs)


def save_overlay_png(background_u8, overlay_mask, out_path, wcs, distance_pc,
                      color=(1, 0, 0), alpha=1.0, scalebar_pc=5.0):
    h_px, w_px = background_u8.shape
    
    aspect_ratio = w_px / h_px
    fig, ax = plt.subplots(figsize=(12 * aspect_ratio, 12))
    
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
    # 1. origin="lower" で FITS 本来の向き（FilFinderのプロット）に合わせる
    # 2. aspect="equal" で縦横の歪み（引き伸ばし）を完全に防止
    ax.imshow(background_u8, cmap="gray", origin="lower", vmin=0, vmax=255, aspect='equal')
    
    # 重ね合わせレイヤーの作成
    overlay_rgba = np.zeros((*overlay_mask.shape, 4), dtype=float)
    overlay_rgba[..., 0] = color[0]  # R (赤)
    overlay_rgba[..., 1] = color[1]  # G
    overlay_rgba[..., 2] = color[2]  # B
    overlay_rgba[..., 3] = overlay_mask.astype(float) * alpha
    
    ax.imshow(overlay_rgba, origin="lower", aspect='equal')
    
    # 画像の物理的なピクセル範囲に座標軸を強制固定（余白の発生を防止）
    ax.set_xlim(0, w_px)
    ax.set_ylim(0, h_px)
    
    ax.axis("off")  # 軸線・枠線を消去

    # 右下にスケールバー（distance_pcを使ってpc→pixel換算）を追加
    add_scalebar(ax, w_px, h_px, wcs, distance_pc, scalebar_pc=scalebar_pc)
    
    fig.savefig(out_path, dpi=150, bbox_inches='tight', pad_inches=0)
    plt.close(fig)

# 全体グレースケール画像 (input_u8) の上に、スケルトン（赤線）を重ね合わせ
save_overlay_png(
    input_u8, skeleton.astype(bool),
    os.path.join(vis_dir, f"{fits_name}_disperse_skeleton_overlay_nsig{NSIG}.png"),
    wcs=w,
    distance_pc=distance_pc,
    color=(1, 0, 0),
    alpha=1.0
)

# 全体グレースケール画像 (input_u8) の上に、太いマスク（赤領域）を重ね合わせ
save_overlay_png(
    input_u8, mask.astype(bool),
    os.path.join(vis_dir, f"{fits_name}_disperse_mask_overlay_nsig{NSIG}.png"),
    wcs=w,
    distance_pc=distance_pc,
    color=(1, 0, 0),
    alpha=0.6
)
print("FilFinder版と完全に一致した重ね合わせ画像を保存しました。")


# ------------------------------------------------------------
# 論文・FilFinder出力と同じ体裁の科学図（inferno + colorbar + タイトル + pixel軸 +
# シアンのフィラメント線 + 右下の横棒スケールバー）を追加で保存する。
# 天体名/距離/beamは冒頭のREGIONS[TARGET_REGION]で切り替わる。
# ------------------------------------------------------------
def save_science_figure(data2d, skeleton_mask, out_path, wcs, distance_pc,
                         display_name, beam_arcsec, n_filaments,
                         scalebar_pc=5.0, vmax=None, cmap="inferno",
                         filament_color="red", filament_alpha=1.0):
    """
    data2d          : 物理値(線形, cm^-2)の柱密度2次元配列。NaNは黒表示になる。
    skeleton_mask   : data2dと同形状の bool配列。Trueの画素にフィラメント線を描く。
                      skeleton(1px芯線)を渡せば細線オーバーレイ、mask(太らせ済み)を
                      渡せば太い領域オーバーレイになる。
    filament_color  : matplotlibが解釈できる色名 (例: "red", "cyan")
    filament_alpha  : フィラメント重ね描きの不透明度 (skeletonは1.0推奨、maskは0.6程度推奨)
    n_filaments     : タイトルに表示するフィラメント数。
                      注意: DisPerSEのVTK出力をそのままread_vtk_skeleton()で数えた
                      セグメント数は、arc-assembly（隣接アークの接続）を経ていないため
                      「科学的に意味のあるフィラメント本数」より過大になりやすい。
                      正式な本数を使いたい場合は、assemble後の本数を別途計算して渡すこと。
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
    fil_rgba = np.zeros((*skeleton_mask.shape, 4), dtype=float)
    fil_rgba[..., 0] = color_rgb[0]
    fil_rgba[..., 1] = color_rgb[1]
    fil_rgba[..., 2] = color_rgb[2]
    fil_rgba[..., 3] = skeleton_mask.astype(float) * filament_alpha
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
    skeleton_mask=skeleton.astype(bool),
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
    skeleton_mask=mask.astype(bool),
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
