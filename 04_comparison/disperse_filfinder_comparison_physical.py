import os
import time
import pickle
import numpy as np
import pandas as pd
import cv2
from scipy.ndimage import binary_dilation
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# ============================================================
# DisPerSE(physical幅マスク) と FilFinder の結果を比較するスクリプト
#   - skeleton(1px)同士は従来通り
#   - mask同士は physical(0.1pc幅)版で比較
# ============================================================

root_dir = "/mnt/sp4t/fujiwara.r"
fits_name = "HGBS_aquilaM2_hires_column_density_map"

WIDTH_PC = 0.1  # DisPerSE側マスク生成スクリプトと同じ値に揃える

disperse_dir = os.path.join(root_dir, fits_name, "disperse")
filfinder_cache_dir = os.path.join(root_dir, fits_name, "cache")
compare_dir = os.path.join(disperse_dir, "compare")
os.makedirs(compare_dir, exist_ok=True)

print("fits_name:", fits_name)
print("disperse_dir:", disperse_dir)


# ------------------------------------------------------------
# 正規化済み入力画像(背景表示用)
# ------------------------------------------------------------
input_u8 = cv2.imread(
    os.path.join(disperse_dir, "full", fits_name + "_input_full_ml.png"),
    cv2.IMREAD_GRAYSCALE
)
if input_u8 is None:
    raise FileNotFoundError("入力画像が読み込めません。パスを確認してください。")


# ------------------------------------------------------------
# DisPerSE側の結果を読み込む(0/1のラベル画像)
#   skeletonは従来通り、maskはphysical版のファイル名に変更
# ------------------------------------------------------------
disperse_skeleton_path = os.path.join(disperse_dir, "full", fits_name + "_skeleton_01.png")
disperse_mask_path = os.path.join(disperse_dir, "full", fits_name + f"_mask_physical{WIDTH_PC}pc_01.png")

print("full/ の中身:", os.listdir(os.path.join(disperse_dir, "full")))

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
skel_cache_path = os.path.join(filfinder_cache_dir, fits_name + "_FilFinder2D-medskel.pkl")
mask_cache_path = os.path.join(filfinder_cache_dir, fits_name + "_FilFinder2D-create_mask.pkl")

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
print("disperse_mask(physical)の True画素数:", disperse_mask.sum())
print("filfinder_skeleton の True画素数:", filfinder_skeleton.sum())
print("filfinder_mask     の True画素数:", filfinder_mask.sum())

print("input_u8 shape:", input_u8.shape)
print("disperse_skeleton shape:", disperse_skeleton.shape)
print("filfinder_skeleton shape:", filfinder_skeleton.shape)

assert disperse_skeleton.shape == filfinder_skeleton.shape, "shapeが一致しません。"
assert disperse_mask.shape == filfinder_mask.shape, "shapeが一致しません(mask)。"


# ============================================================
# 定量比較: IoU / Precision / Recall / F1
# DisPerSE=pred, FilFinder=gt(暫定的な正解)として計算
# ============================================================
def compute_metrics(pred, gt):
    tp = np.logical_and(pred, gt).sum()
    fp = np.logical_and(pred, ~gt).sum()
    fn = np.logical_and(~pred, gt).sum()

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    iou = tp / (tp + fp + fn) if (tp + fp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {"TP": int(tp), "FP": int(fp), "FN": int(fn),
            "Precision": precision, "Recall": recall, "IoU": iou, "F1": f1}

metrics_skeleton = compute_metrics(disperse_skeleton, filfinder_skeleton)
metrics_mask = compute_metrics(disperse_mask, filfinder_mask)

df = pd.DataFrame([
    {"comparison": "skeleton(1px)", **metrics_skeleton},
    {"comparison": f"mask(physical {WIDTH_PC}pc)", **metrics_mask},
])
print(df)

csv_path = os.path.join(compare_dir, fits_name + f"_disperse_physical{WIDTH_PC}pc_vs_filfinder_metrics.csv")
df.to_csv(csv_path, index=False)
print("メトリクス保存:", csv_path)


# ============================================================
# 色分けオーバーレイ可視化
#   黄 = 両方一致(TP), 赤 = DisPerSEのみ(FP), 青 = FilFinderのみ(FN)
# ============================================================
def save_comparison_overlay(background_u8, pred, gt, out_path, title,
                              vis_dilate_iter=0, bg_brightness=0.6, overlay_alpha=0.85):
    fig, ax = plt.subplots(figsize=(16, 16))

    dimmed_background = (background_u8.astype(float) * bg_brightness).astype(np.uint8)
    ax.imshow(dimmed_background, cmap="gray", origin="upper", vmin=0, vmax=255)

    if vis_dilate_iter > 0:
        pred_vis = binary_dilation(pred, iterations=vis_dilate_iter)
        gt_vis = binary_dilation(gt, iterations=vis_dilate_iter)
    else:
        pred_vis, gt_vis = pred, gt

    both = pred_vis & gt_vis
    disperse_only = pred_vis & ~gt_vis
    filfinder_only = ~pred_vis & gt_vis

    overlay = np.zeros((*pred_vis.shape, 4), dtype=float)

    overlay[both, 0] = 1.0; overlay[both, 1] = 1.0; overlay[both, 2] = 0.0; overlay[both, 3] = overlay_alpha
    overlay[disperse_only, 0] = 1.0; overlay[disperse_only, 1] = 0.0; overlay[disperse_only, 2] = 0.0; overlay[disperse_only, 3] = overlay_alpha
    overlay[filfinder_only, 0] = 0.0; overlay[filfinder_only, 1] = 0.4; overlay[filfinder_only, 2] = 1.0; overlay[filfinder_only, 3] = overlay_alpha

    ax.imshow(overlay, origin="upper")
    ax.axis("off")
    ax.set_title(title, color="white", fontsize=16)
    fig.savefig(out_path, dpi=150, bbox_inches="tight", pad_inches=0, facecolor="black")
    plt.close(fig)

    print(f"  [{os.path.basename(out_path)}] 黄={both.sum()}px, 赤={disperse_only.sum()}px, 青={filfinder_only.sum()}px")

ts = int(time.time())

save_comparison_overlay(
    input_u8, disperse_skeleton, filfinder_skeleton,
    os.path.join(compare_dir, f"{fits_name}_compare_skeleton_{WIDTH_PC}pc_{ts}.png"),
    f"Yellow=Both / Red=DisPerSE only / Blue=FilFinder only (skeleton, {WIDTH_PC}pc)",
    vis_dilate_iter=8, bg_brightness=0.6, overlay_alpha=0.9
)

save_comparison_overlay(
    input_u8, disperse_mask, filfinder_mask,
    os.path.join(compare_dir, f"{fits_name}_compare_mask_{WIDTH_PC}pc_{ts}.png"),
    f"Yellow=Both / Red=DisPerSE only / Blue=FilFinder only (mask, {WIDTH_PC}pc)",
    vis_dilate_iter=0, bg_brightness=0.6, overlay_alpha=0.7
)

print("DONE")
