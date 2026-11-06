import os
import numpy as np
import cv2
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

# ============================================================
# セマンティックセグメンテーション説明用: 3クラス可視化スクリプト(修正版)
#   label_3class_path の .npy が存在しないため、
#   装飾済みの3クラス重ね合わせPNG(Image1)から色を頼りにラベルを復元する
# ============================================================

root_dir = "/mnt/sp4t/fujiwara.r"
fits_name = "HGBS_aquilaM2_hires_column_density_map"

# 【変更点】.npyではなく、実際にアップロードされた2枚のPNGを直接指定
input_full_ml_path = os.path.join(root_dir, fits_name, "full", fits_name + "_vis_input.png")
label_overlay_png_path = os.path.join(root_dir,fits_name, "vis_3class_map",fits_name + "_3class_full_map.png")

out_path = os.path.join(root_dir, fits_name, fits_name + "_semseg_explainer.png")

matplotlib.rcParams["font.family"] = "IPAexGothic"

FLIP_LABEL_VERTICALLY = True  # Image1とImage2で上下が逆に見えたらこれを切り替える

CLASS_COLORS = {
    0: (0.05, 0.05, 0.08),
    1: (1.00, 0.85, 0.10),
    2: (0.20, 0.75, 0.95),
}
CLASS_LABELS = {
    0: "背景 (Background)",
    1: "Class 1: Confirmed Core",
    2: "Class 2: Candidate",
}

# ============================================================
# 入力画像を読み込む
# ============================================================
input_u8 = cv2.imread(input_full_ml_path, cv2.IMREAD_GRAYSCALE)
if input_u8 is None:
    raise FileNotFoundError(f"見つかりません: {input_full_ml_path}")


# ============================================================
# 装飾済みPNG(Image1)から、実データ部分だけを自動検出してクロップする
#   白い余白(タイトル・軸目盛りの外側)を除いた、最大の非白色領域を探す
# ============================================================
def crop_plot_area(bgr_img, white_thresh=250):
    gray = cv2.cvtColor(bgr_img, cv2.COLOR_BGR2GRAY)
    non_white_mask = (gray < white_thresh).astype(np.uint8)

    # 連結成分の中で最大のもの(=カラーバーや軸目盛りの数字ではなく、本体の図)を選ぶ
    n_labels, labels, stats, _ = cv2.connectedComponentsWithStats(non_white_mask, connectivity=8)
    if n_labels <= 1:
        return bgr_img  # 何も検出できなければそのまま返す

    largest_label = 1 + np.argmax(stats[1:, cv2.CC_STAT_AREA])
    x, y, w, h = stats[largest_label, cv2.CC_STAT_LEFT], stats[largest_label, cv2.CC_STAT_TOP], \
                 stats[largest_label, cv2.CC_STAT_WIDTH], stats[largest_label, cv2.CC_STAT_HEIGHT]
    return bgr_img[y:y + h, x:x + w]

label_overlay_bgr = cv2.imread(label_overlay_png_path)
if label_overlay_bgr is None:
    raise FileNotFoundError(f"見つかりません: {label_overlay_png_path}")

cropped = crop_plot_area(label_overlay_bgr)
print("クロップ後のサイズ:", cropped.shape)


# ============================================================
# 色分類でクラスラベルを復元する
#   黄色(Class1) / 水色(Class2) / それ以外(背景、magmaカラーマップの色調含む)
# ============================================================
def classify_by_color(bgr_img):
    b = bgr_img[:, :, 0].astype(int)
    g = bgr_img[:, :, 1].astype(int)
    r = bgr_img[:, :, 2].astype(int)

    label = np.zeros(bgr_img.shape[:2], dtype=np.uint8)

    # Class1(黄色): R,G高い かつ B低い(magma背景の白っぽい高輝度部分はBも高いので除外される)
    class1_mask = (r > 180) & (g > 150) & (b < 120)

    # Class2(水色): G,B高い かつ R低い(magmaには青系の色が無いので誤検出しにくい)
    class2_mask = (g > 120) & (b > 150) & (r < 150)

    label[class2_mask] = 2
    label[class1_mask] = 1  # Class1を後で上書き(優先度を高くする)
    return label

label_map = classify_by_color(cropped)

if FLIP_LABEL_VERTICALLY:
    label_map = np.flipud(label_map)

print("復元したラベルのクラス比率:")
for c in [0, 1, 2]:
    print(f"  class {c}: {100 * (label_map == c).mean():.2f}%")


# ============================================================
# input_u8 のサイズにラベルをリサイズして合わせる(最近傍補間でクラス値を維持)
# ============================================================
label_map_resized = cv2.resize(
    label_map, (input_u8.shape[1], input_u8.shape[0]),
    interpolation=cv2.INTER_NEAREST
)

H, W = input_u8.shape


# ============================================================
# 3クラスを単色で塗り分けたRGB画像を作る
# ============================================================
seg_rgb = np.zeros((H, W, 3), dtype=float)
for class_id, color in CLASS_COLORS.items():
    mask = (label_map_resized == class_id)
    seg_rgb[mask] = color


# ============================================================
# 描画
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(20, 10))

axes[0].imshow(input_u8, cmap="gray", origin="upper")
axes[0].set_title("① 入力画像(柱密度マップ)", fontsize=16)
axes[0].axis("off")

axes[1].imshow(seg_rgb, origin="upper")
axes[1].set_title("② セマンティックセグメンテーション結果(ピクセル単位のクラス分類)", fontsize=16)
axes[1].axis("off")

legend_elements = [Patch(facecolor=CLASS_COLORS[c], label=CLASS_LABELS[c]) for c in sorted(CLASS_COLORS)]
fig.legend(handles=legend_elements, loc="lower center", ncol=3, fontsize=13, bbox_to_anchor=(0.5, -0.02))

fig.suptitle(f"{fits_name} — Semantic Segmentation Example", fontsize=18)
fig.tight_layout(rect=[0, 0.05, 1, 0.95])
fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("保存:", out_path)
print("DONE")
