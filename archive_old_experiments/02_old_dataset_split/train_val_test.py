import os
import numpy as np
import pandas as pd
import cv2
from scipy.ndimage import binary_erosion

# ============================================================
# Aquila 3クラスラベルの train/val/test 分割スクリプト(決定的版)
#   - ランダムシャッフルなし(空間順のまま層内で先頭から割り当て)
#   - 背景/Class1(黄)/Class2(水色)の比率で層化
#   - ジグザグした観測範囲の端を侵食して除外
# ============================================================

root_dir = "/mnt/sp4t/fujiwara.r"
fits_name = "HGBS_aquilaM2_hires_column_density_map"

# 観測範囲の判定に使う(0=範囲外という前提のグレースケール画像)
footprint_ref_path = os.path.join(root_dir, fits_name, "disperse", "full", fits_name + "_input_full_ml.png")
# 3クラス重ね合わせの完成図(タイトル・凡例つき)。npyが無い場合は色から復元する
label_overlay_png_path = os.path.join(root_dir, fits_name, fits_name + "_3class_full_map.png")

out_root = os.path.join(root_dir, fits_name, "split_3class_balanced")
os.makedirs(out_root, exist_ok=True)
split_map_png_path = os.path.join(out_root, fits_name + "_split_id_map.png")
label_map_npy_path = os.path.join(out_root, fits_name + "_label_map_recovered.npy")

BLOCK_SIZE = 768
MARGIN_PX = 32          # ブロック境界の緩衝帯
EDGE_ERODE_PX = 32      # 観測範囲の外周(ジグザグ境界)から侵食して除外する幅
PATCH_SIZE = 256
STRIDE = 16
TILES_PER_DIR = 1000
N_STRAT_BINS = 5
TRAIN_RATIO, VAL_RATIO, TEST_RATIO = 0.7, 0.15, 0.15


# ============================================================
# ステップ1: 観測範囲マスクを作り、端を侵食して除外する
# ============================================================
footprint_gray = cv2.imread(footprint_ref_path, cv2.IMREAD_GRAYSCALE)
if footprint_gray is None:
    raise FileNotFoundError(f"見つかりません: {footprint_ref_path}")

footprint_mask = footprint_gray > 0  # 0=観測範囲外という前提
H, W = footprint_mask.shape

struct = np.ones((3, 3), dtype=bool)
footprint_mask_eroded = binary_erosion(footprint_mask, structure=struct, iterations=EDGE_ERODE_PX)
print(f"侵食前の有効画素割合: {footprint_mask.mean():.3f} → 侵食後: {footprint_mask_eroded.mean():.3f}")


# ============================================================
# ステップ2: 3クラスラベルを、装飾済みPNGの色から復元する
#   (label_3class.npy が既にあるなら、この関数呼び出しをnp.load()に差し替えてOK)
# ============================================================
def crop_plot_area(bgr_img, white_thresh=250):
    gray = cv2.cvtColor(bgr_img, cv2.COLOR_BGR2GRAY)
    non_white_mask = (gray < white_thresh).astype(np.uint8)
    n_labels, labels, stats, _ = cv2.connectedComponentsWithStats(non_white_mask, connectivity=8)
    if n_labels <= 1:
        return bgr_img, (0, 0)
    largest_label = 1 + np.argmax(stats[1:, cv2.CC_STAT_AREA])
    x, y, w, h = stats[largest_label, cv2.CC_STAT_LEFT], stats[largest_label, cv2.CC_STAT_TOP], \
                 stats[largest_label, cv2.CC_STAT_WIDTH], stats[largest_label, cv2.CC_STAT_HEIGHT]
    return bgr_img[y:y + h, x:x + w], (x, y)

def classify_by_color(bgr_img):
    b = bgr_img[:, :, 0].astype(int)
    g = bgr_img[:, :, 1].astype(int)
    r = bgr_img[:, :, 2].astype(int)
    label = np.zeros(bgr_img.shape[:2], dtype=np.uint8)
    class1_mask = (r > 180) & (g > 150) & (b < 120)   # 黄色
    class2_mask = (g > 120) & (b > 150) & (r < 150)   # 水色
    label[class2_mask] = 2
    label[class1_mask] = 1  # 黄色を優先(重複時は黄色を採用)
    return label

label_overlay_bgr = cv2.imread(label_overlay_png_path)
if label_overlay_bgr is None:
    raise FileNotFoundError(f"見つかりません: {label_overlay_png_path}")

cropped, _ = crop_plot_area(label_overlay_bgr)
label_map_small = classify_by_color(cropped)

# footprintと同じ解像度にリサイズ(最近傍補間でクラス値を保持)
label_map = cv2.resize(label_map_small, (W, H), interpolation=cv2.INTER_NEAREST)
np.save(label_map_npy_path, label_map)
print("復元ラベルを保存:", label_map_npy_path)

for c, name in [(0, "背景"), (1, "Class1(黄)"), (2, "Class2(水色)")]:
    frac = (label_map[footprint_mask_eroded] == c).mean()
    print(f"  全体でのクラス比率 [{name}]: {frac*100:.2f}%")


# ============================================================
# ステップ3: BLOCK_SIZE単位でブロック分割し、各クラスの比率を計算
# ============================================================
n_by = int(np.ceil(H / BLOCK_SIZE))
n_bx = int(np.ceil(W / BLOCK_SIZE))

blocks = []
for by in range(n_by):
    for bx in range(n_bx):
        y0, y1 = by * BLOCK_SIZE, min((by + 1) * BLOCK_SIZE, H)
        x0, x1 = bx * BLOCK_SIZE, min((bx + 1) * BLOCK_SIZE, W)

        block_footprint = footprint_mask_eroded[y0:y1, x0:x1]
        valid_frac = block_footprint.mean()
        if valid_frac < 0.3:
            continue  # ほぼ観測範囲外のブロックはスキップ

        block_label = label_map[y0:y1, x0:x1]
        frac_class1 = float((block_label == 1).mean())
        frac_class2 = float((block_label == 2).mean())

        blocks.append({
            "by": by, "bx": bx, "y0": y0, "y1": y1, "x0": x0, "x1": x1,
            "frac_class1": frac_class1, "frac_class2": frac_class2,
        })

blocks_df = pd.DataFrame(blocks)
print(f"有効ブロック数: {len(blocks_df)} / 全体: {n_by * n_bx}")


# ============================================================
# ステップ4: 3クラスの比率で層化し、乱数を使わず空間順のまま割り当てる
#   Class1(希少)を重視した複合スコアでビン分けする
# ============================================================
blocks_df["strat_score"] = blocks_df["frac_class1"] * 2.0 + blocks_df["frac_class2"] * 1.0
blocks_df["strat_bin"] = pd.qcut(
    blocks_df["strat_score"].rank(method="first"), q=N_STRAT_BINS, labels=False
)
blocks_df["split"] = "unassigned"

for b in range(N_STRAT_BINS):
    # 空間順(行優先)のまま並べる。shuffleは一切行わない。
    sub = blocks_df[blocks_df["strat_bin"] == b].sort_values(["by", "bx"])
    idx = sub.index.to_numpy()
    n = len(idx)
    n_train = int(round(n * TRAIN_RATIO))
    n_val = int(round(n * VAL_RATIO))
    blocks_df.loc[idx[:n_train], "split"] = "train"
    blocks_df.loc[idx[n_train:n_train + n_val], "split"] = "val"
    blocks_df.loc[idx[n_train + n_val:], "split"] = "test"

print(blocks_df["split"].value_counts())

# 各splitのクラス構成比を確認(ここが均等に近いほど良い)
for split_name in ["train", "val", "test"]:
    sub = blocks_df[blocks_df["split"] == split_name]
    pixels = [label_map[row.y0:row.y1, row.x0:row.x1][footprint_mask_eroded[row.y0:row.y1, row.x0:row.x1]]
              for _, row in sub.iterrows()]
    if pixels:
        concat = np.concatenate(pixels)
        vals, counts = np.unique(concat, return_counts=True)
        ratio = counts / counts.sum()
        print(f"[{split_name}] class比率:", dict(zip(vals.astype(int), np.round(ratio, 4))))


# ============================================================
# ステップ5: split_id_map を構築(端の侵食 + ブロック境界マージンの両方を反映)
# ============================================================
SPLIT_TO_ID = {"train": 0, "val": 1, "test": 2, "unassigned": -1}
split_id_map = np.full((H, W), -1, dtype=int)

for _, row in blocks_df.iterrows():
    split_id_map[row.y0:row.y1, row.x0:row.x1] = SPLIT_TO_ID[row.split]

# 観測範囲の端(侵食で除外された領域)を無効化
split_id_map[~footprint_mask_eroded] = -1

# ブロック境界のマージンも無効化
margin_excluded = np.zeros((H, W), dtype=bool)
for _, row in blocks_df.iterrows():
    y0, y1, x0, x1 = row.y0, row.y1, row.x0, row.x1
    margin_excluded[y0:min(y0 + MARGIN_PX, y1), x0:x1] = True
    margin_excluded[max(y1 - MARGIN_PX, y0):y1, x0:x1] = True
    margin_excluded[y0:y1, x0:min(x0 + MARGIN_PX, x1)] = True
    margin_excluded[y0:y1, max(x1 - MARGIN_PX, x0):x1] = True
split_id_map[margin_excluded] = -1

print("最終的な有効画素割合:", (split_id_map >= 0).mean())

split_id_png = (split_id_map + 1).astype(np.uint8)  # -1→0, 0→1, 1→2, 2→3
cv2.imwrite(split_map_png_path, split_id_png)
print("split_id_map(PNG)を保存:", split_map_png_path)


# ============================================================
# ステップ6: パッチ切り出し
# ============================================================
def check_make_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)

ID_TO_SPLIT = {v: k for k, v in SPLIT_TO_ID.items() if v != -1}
step = STRIDE
tile_counts = {"train": 0, "val": 0, "test": 0, "skipped": 0}

y = 0
while y + PATCH_SIZE <= H:
    x = 0
    while x + PATCH_SIZE <= W:
        patch_footprint = footprint_gray[y:y + PATCH_SIZE, x:x + PATCH_SIZE]
        patch_label = label_map[y:y + PATCH_SIZE, x:x + PATCH_SIZE]
        patch_split_id = split_id_map[y:y + PATCH_SIZE, x:x + PATCH_SIZE]

        x += step

        if (patch_footprint == 0).any():
            tile_counts["skipped"] += 1
            continue

        unique_ids = np.unique(patch_split_id)
        if len(unique_ids) != 1 or unique_ids[0] == -1:
            tile_counts["skipped"] += 1
            continue

        split_name = ID_TO_SPLIT[unique_ids[0]]
        o_name = fits_name + "_" + str(x - step).zfill(5) + "-" + str(y).zfill(5)

        in_dir = os.path.join(out_root, split_name, "in")
        out_dir = os.path.join(out_root, split_name, "out")
        tile_idx = tile_counts[split_name]
        sub_in = os.path.join(in_dir, "tiles_{:03d}".format(tile_idx // TILES_PER_DIR))
        sub_out = os.path.join(out_dir, "tiles_{:03d}".format(tile_idx // TILES_PER_DIR))
        check_make_dir(sub_in)
        check_make_dir(sub_out)

        cv2.imwrite(os.path.join(sub_in, o_name + "_in.png"), patch_footprint)
        cv2.imwrite(os.path.join(sub_out, o_name + "_out.png"), patch_label.astype("uint8"))

        tile_counts[split_name] += 1
    y += step

print("パッチ生成結果:", tile_counts)
print("DONE")


