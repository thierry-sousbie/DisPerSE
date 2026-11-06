import os
import pickle
import numpy as np
import pandas as pd
import cv2
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# ============================================================
# DisPerSE と FilFinder の結果を比較するスクリプト
#   - skeleton(1px)同士、mask(太い領域)同士をそれぞれ比較
#   - IoU / Precision / Recall / F1 を算出
#   - 色分けオーバーレイ(白=両方一致、赤=DisPerSEのみ、緑=FilFinderのみ)
# ============================================================

root_dir = "/mnt/sp4t/fujiwara.r"
fits_name = "HGBS_aquilaM2_hires_column_density_map"

disperse_dir = os.path.join(root_dir, fits_name, "disperse")
filfinder_cache_dir = os.path.join(root_dir, fits_name, "cache")
compare_dir = os.path.join(disperse_dir, "compare")
os.makedirs(compare_dir, exist_ok=True)

# ------------------------------------------------------------
# 正規化済み入力画像(背景表示用、DisPerSEスクリプトで既に作成済み)
# ------------------------------------------------------------
input_u8 = cv2.imread(
    os.path.join(disperse_dir, "full", fits_name + "_input_full_ml.png"),
    cv2.IMREAD_GRAYSCALE
)

# ------------------------------------------------------------
# DisPerSE側の結果を読み込む(0/1のラベル画像)
# ------------------------------------------------------------
disperse_skeleton = cv2.imread(
    os.path.join(disperse_dir, "full", fits_name + "_skeleton_label01.png"),
    cv2.IMREAD_GRAYSCALE
).astype(bool)

disperse_mask = cv2.imread(
    os.path.join(disperse_dir, "full", fits_name + "_mask_label01.png"),
    cv2.IMREAD_GRAYSCALE
).astype(bool)

# ------------------------------------------------------------
# FilFinder側の結果を読み込む(pickleキャッシュから)
# ------------------------------------------------------------
with open(os.path.join(filfinder_cache_dir, fits_name + "_FilFinder2D-medskel.pkl"), "rb") as f:
    fil_ff_skel_obj = pickle.load(f)
filfinder_skeleton = fil_ff_skel_obj.skeleton.astype(bool)

with open(os.path.join(filfinder_cache_dir, fits_name + "_FilFinder2D-create_mask.pkl"), "rb") as f:
    fil_ff_mask_obj = pickle.load(f)
filfinder_mask = fil_ff_mask_obj.mask.astype(bool)

# shapeが一致しているか確認(異なるならここで気づける)
print("input_u8 shape:", input_u8.shape)
print("disperse_skeleton shape:", disperse_skeleton.shape)
print("filfinder_skeleton shape:", filfinder_skeleton.shape)
assert disperse_skeleton.shape == filfinder_skeleton.shape, "shapeが一致しません。再確認してください。"


# ============================================================
# 定量比較: IoU / Precision / Recall / F1
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

# DisPerSEをpredict、FilFinderをground truthとみなして計算(逆でも数値は対称ではないので注意)
metrics_skeleton = compute_metrics(disperse_skeleton, filfinder_skeleton)
metrics_mask = compute_metrics(disperse_mask, filfinder_mask)

df = pd.DataFrame([
    {"comparison": "skeleton(1px)", **metrics_skeleton},
    {"comparison": "mask(太い領域)", **metrics_mask},
])
print(df)

csv_path = os.path.join(compare_dir, fits_name + "_disperse_vs_filfinder_metrics.csv")
df.to_csv(csv_path, index=False)
print("メトリクス保存:", csv_path)


# ============================================================
# 色分けオーバーレイ可視化
#   白 = 両方一致(TP), 赤 = DisPerSEのみ(FP), 緑 = FilFinderのみ(FN)
# ============================================================
def save_comparison_overlay(background_u8, pred, gt, out_path, title):
    fig, ax = plt.subplots(figsize=(16, 16))
    ax.imshow(background_u8, cmap="gray", origin="upper", vmin=0, vmax=255)

    overlay = np.zeros((*pred.shape, 4), dtype=float)

    both = pred & gt
    disperse_only = pred & ~gt
    filfinder_only = ~pred & gt

   overlay[both, 0] = 0; overlay[both, 1] = 1; overlay[both, 2] = 1; overlay[both, 3] = 0.9   # シアン
   overlay[disperse_only, 0] = 1; overlay[disperse_only, 3] = 0.9                              # 赤
   overlay[filfinder_only, 0] = 1; overlay[filfinder_only, 1] = 0.5; overlay[filfinder_only, 2] = 0; overlay[filfinder_only, 3] = 0.9  # オレンジ
    ax.imshow(overlay, origin="upper")
    ax.axis("off")
    ax.set_title(title, color="white", fontsize=14)
    fig.savefig(out_path, dpi=150, bbox_inches="tight", pad_inches=0, facecolor="black")
    plt.close(fig)

save_comparison_overlay(
    input_u8, disperse_skeleton, filfinder_skeleton,
    os.path.join(compare_dir, fits_name + "_compare_skeleton.png"),
    "White=Both / Red=DisPerSE only / Green=FilFinder only (skeleton)"
)

save_comparison_overlay(
    input_u8, disperse_mask, filfinder_mask,
    os.path.join(compare_dir, fits_name + "_compare_mask.png"),
    "White=Both / Red=DisPerSE only / Green=FilFinder only (mask)"
)

print("比較オーバーレイを保存しました:")
print(" -", os.path.join(compare_dir, fits_name + "_compare_skeleton.png"))
print(" -", os.path.join(compare_dir, fits_name + "_compare_mask.png"))
print("DONE")
