import os
import numpy as np
import cv2
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ============================================================
# train/val/test 分割の可視化スクリプト(可視化専用)
#   入力: 入力画像PNG、3クラスラベル画像、split_id_map PNG(スクリプトAの出力)
#   分割の計算は一切行わない
# ============================================================

root_dir = "/mnt/sp4t/fujiwara.r"
fits_name = "HGBS_aquilaM2_hires_column_density_map"

input_full_path = os.path.join(root_dir, fits_name, "disperse", "full", fits_name + "_input_full_ml.png")
label_3class_path = os.path.join(root_dir, fits_name, "compare", fits_name + "_label_3class.npy")
split_map_png_path = os.path.join(root_dir, fits_name, "split_3class", fits_name + "_split_id_map.png")

vis_out_path = os.path.join(root_dir, fits_name, fits_name + "_split_visualization.png")

PATCH_SIZE = 256

# ============================================================
# 読み込み(すべて計算済みファイルを読むだけ)
# ============================================================
def load_image_any(path):
    if path.endswith(".npy"):
        return np.load(path)
    img = cv2.imread(path, cv2.IMREAD_UNCHANGED)
    if img is None:
        raise FileNotFoundError(f"読み込めません: {path}")
    return img

input_u8 = load_image_any(input_full_path)
label_map_raw = load_image_any(label_3class_path)
split_id_png = cv2.imread(split_map_png_path, cv2.IMREAD_UNCHANGED)

if split_id_png is None:
    raise FileNotFoundError(f"split_id_map.pngが見つかりません: {split_map_png_path}")

split_id_map = split_id_png.astype(int) - 1  # PNG保存時の+1を戻す(0→-1, 1→0, 2→1, 3→2)

if label_map_raw.max() > 2:
    uniq = np.unique(label_map_raw)
    remap = {v: i for i, v in enumerate(sorted(uniq))}
    label_map = np.vectorize(remap.get)(label_map_raw).astype(int)
else:
    label_map = label_map_raw.astype(int)

assert input_u8.shape[:2] == split_id_map.shape[:2] == label_map.shape[:2], "shapeが一致しません。"
H, W = split_id_map.shape[:2]

SPLIT_TO_ID = {"train": 0, "val": 1, "test": 2}
SPLIT_COLOR = {
    "train": (0.20, 0.45, 0.95),
    "val":   (0.20, 0.80, 0.30),
    "test":  (0.95, 0.25, 0.25),
}
SPLIT_LABEL_JP = {"train": "学習", "val": "検証", "test": "テスト"}


# ============================================================
# オーバーレイ作成
# ============================================================
def build_split_overlay(shape):
    overlay = np.zeros((*shape, 4), dtype=float)
    for split_name, split_id in SPLIT_TO_ID.items():
        mask = (split_id_map == split_id)
        color = SPLIT_COLOR[split_name]
        overlay[mask, 0] = color[0]
        overlay[mask, 1] = color[1]
        overlay[mask, 2] = color[2]
        overlay[mask, 3] = 0.45
    return overlay

overlay_rgba = build_split_overlay((H, W))


# ============================================================
# 代表パッチの位置を選ぶ(各splitの中でフィラメント密度最大の場所)
#   ※ split_id_mapから直接、走査して求める(blocks_dfに依存しない)
# ============================================================
def find_best_patch(split_id_target):
    mask = (split_id_map == split_id_target)
    if not mask.any():
        return None

    ys, xs = np.where(mask)
    y_min, y_max, x_min, x_max = ys.min(), ys.max(), xs.min(), xs.max()

    best_score = -1
    best_loc = None
    step = 128  # 走査の粗さ(処理時間短縮のため。細かくしたければ小さくする)
    for y0 in range(y_min, max(y_min, y_max - PATCH_SIZE), step):
        for x0 in range(x_min, max(x_min, x_max - PATCH_SIZE), step):
            y1, x1 = y0 + PATCH_SIZE, x0 + PATCH_SIZE
            if y1 > H or x1 > W:
                continue
            if not np.all(split_id_map[y0:y1, x0:x1] == split_id_target):
                continue
            score = (label_map[y0:y1, x0:x1] > 0).mean()
            if score > best_score:
                best_score = score
                best_loc = (y0, y1, x0, x1)
    return best_loc

example_locations = {}
for split_name in ["test", "train", "val"]:
    loc = find_best_patch(SPLIT_TO_ID[split_name])
    if loc is not None:
        example_locations[split_name] = loc

print("代表パッチ位置:", example_locations)


# ============================================================
# 描画
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(22, 12))
panels = [
    {"ax": axes[0], "background": input_u8, "cmap": "gray", "title": "入力画像"},
    {"ax": axes[1], "background": (label_map > 0).astype(np.uint8) * 255, "cmap": "gray_r", "title": "暫定ラベル画像"},
]

inset_x_positions = [0.05, 0.38, 0.71]
inset_width, inset_height = 0.26, 0.30

for panel in panels:
    ax = panel["ax"]
    ax.imshow(panel["background"], cmap=panel["cmap"], origin="upper")
    ax.imshow(overlay_rgba, origin="upper")
    ax.set_title(panel["title"], fontsize=16)
    ax.axis("off")

    for i, split_name in enumerate(["test", "train", "val"]):
        if split_name not in example_locations:
            continue
        y0, y1, x0, x1 = example_locations[split_name]

        rect = Rectangle((x0, y0), x1 - x0, y1 - y0,
                          linewidth=2.5, edgecolor="navy", facecolor="none", zorder=5)
        ax.add_patch(rect)

        inset_ax = inset_axes(ax, width="100%", height="100%",
                               bbox_to_anchor=(inset_x_positions[i], 1.02, inset_width, inset_height),
                               bbox_transform=ax.transAxes, loc="lower left", borderpad=0)
        crop = panel["background"][y0:y1, x0:x1]
        inset_ax.imshow(crop, cmap=panel["cmap"], origin="upper")
        inset_ax.set_xticks([]); inset_ax.set_yticks([])
        for spine in inset_ax.spines.values():
            spine.set_edgecolor(SPLIT_COLOR[split_name])
            spine.set_linewidth(3)

        inset_ax.text(0.05, 0.92, SPLIT_LABEL_JP[split_name],
                       transform=inset_ax.transAxes, fontsize=14, fontweight="bold",
                       color="black", bbox=dict(facecolor=SPLIT_COLOR[split_name], alpha=0.9, edgecolor="none"),
                       va="top", ha="left")

        con = ConnectionPatch(
            xyA=((x0 + x1) / 2, y0), coordsA=ax.transData,
            xyB=(0.5, 0), coordsB=inset_ax.transAxes,
            linestyle="dotted", linewidth=2, color="navy", zorder=4
        )
        fig.add_artist(con)

legend_patches = [
    plt.Rectangle((0, 0), 1, 1, facecolor=SPLIT_COLOR[s], alpha=0.7, label=f"{SPLIT_LABEL_JP[s]} ({s})")
    for s in ["train", "val", "test"]
]
fig.legend(handles=legend_patches, loc="lower center", ncol=3, fontsize=13, bbox_to_anchor=(0.5, -0.02))

fig.suptitle(fits_name + " — train/val/test split visualization", fontsize=18)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
fig.savefig(vis_out_path, dpi=150, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("保存:", vis_out_path)
print("DONE")
