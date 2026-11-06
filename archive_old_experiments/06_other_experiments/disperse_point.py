import os
import glob
import numpy as np
import cv2
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

root_dir = "/mnt/sp4t/fujiwara.r"
fits_name = "HGBS_aquilaM2_hires_column_density_map"
NSIG = 1  # これまで確定させた値

disperse_dir = os.path.join(root_dir, fits_name, "disperse")
raw_dir = os.path.join(disperse_dir, "raw")
vis_dir = os.path.join(disperse_dir, "vis")

input_full_ml_path = os.path.join(disperse_dir, "full", fits_name + "_input_full_ml.png")
input_u8 = cv2.imread(input_full_ml_path, cv2.IMREAD_GRAYSCALE)

# ============================================================
# ASCII出力を探す(nsig=1版のNDskl_asciiが必要。無ければskelconvで作る)
# ============================================================
ndskl_i = sorted(glob.glob(os.path.join(raw_dir, f"*_s{NSIG}.up.NDskl")))[-1]
ascii_candidates = sorted(glob.glob(ndskl_i + "*.a.NDskl"))

if not ascii_candidates:
    import subprocess
    DISPERSE_BIN_DIR = "/home1/fujiwara.r/PycharmProjects/DisPerSE/build/src/mse"
    cmd = [os.path.join(DISPERSE_BIN_DIR, "skelconv"), ndskl_i, "-to", "NDskl_ascii", "-outDir", raw_dir]
    subprocess.run(cmd, capture_output=True, text=True)
    ascii_candidates = sorted(glob.glob(ndskl_i + "*.a.NDskl"))

ascii_path = ascii_candidates[-1]
print("使用するASCII:", ascii_path)


# ============================================================
# [CRITICAL POINTS] セクションをパース(type + 座標のみ取得)
# ============================================================
def parse_critical_points(fpath):
    with open(fpath, "r") as f:
        lines = [l.rstrip("\n") for l in f]

    idx = lines.index("[CRITICAL POINTS]")
    nnodes = int(lines[idx + 1].strip())
    ptr = idx + 2

    nodes = []
    for _ in range(nnodes):
        header = lines[ptr].split()
        ptr += 1
        node_type = int(header[0])
        x, y = float(header[1]), float(header[2])
        nodes.append((node_type, x, y))

        nnext = int(lines[ptr].strip())
        ptr += 1
        ptr += nnext  # next/fil_indexの行はスキップ

    return nodes

nodes = parse_critical_points(ascii_path)
print("臨界点の総数:", len(nodes))

nodes_arr = np.array(nodes)
type_counts = {int(t): int((nodes_arr[:, 0] == t).sum()) for t in np.unique(nodes_arr[:, 0])}
print("種類ごとの個数(0=極小, 1=鞍点, 2=極大):", type_counts)


# ============================================================
# 可視化: 反転グレースケール背景に、種類ごとに色分けしてプロット
# ============================================================
inverted_u8 = 255 - input_u8

fig, ax = plt.subplots(figsize=(16, 16))
ax.imshow(inverted_u8, cmap="gray", origin="upper")

TYPE_STYLE = {
    0: {"color": "blue",   "marker": "v", "label": "極小 (Minimum)", "size": 15},
    1: {"color": "yellow", "marker": "^", "label": "鞍点 (Saddle)",  "size": 10},
    2: {"color": "red",    "marker": "*", "label": "極大 (Maximum)", "size": 20},
}

for t, style in TYPE_STYLE.items():
    pts = nodes_arr[nodes_arr[:, 0] == t]
    if len(pts) == 0:
        continue
    ax.scatter(pts[:, 1], pts[:, 2], c=style["color"], marker=style["marker"],
               s=style["size"], label=f'{style["label"]} (n={len(pts)})',
               edgecolors="black", linewidths=0.3, zorder=5)

ax.legend(loc="upper right", fontsize=13, framealpha=0.9)
ax.set_title(f"{fits_name} — DisPerSE Critical Points (nsig={NSIG})", fontsize=16)
ax.axis("off")

out_path = os.path.join(vis_dir, fits_name + f"_disperse_critical_points_nsig{NSIG}.png")
fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("保存:", out_path)
