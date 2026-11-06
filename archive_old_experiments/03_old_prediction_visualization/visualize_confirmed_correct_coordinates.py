#!/usr/bin/env python3
from __future__ import annotations

"""
元の柱密度FITS・元の3クラスラベル・AI予測を
すべて同じ「元画像座標系」で合成してからPNG保存する可視化スクリプト。

これにより、Matplotlibで保存済みのPNGに後から貼ったときに生じた
座標ずれを防ぐ。

暫定3クラスラベル:
    0   = background
    1   = confirmed_filament
    2   = candidate_filament
    255 = ignore

AI予測:
    background vs confirmed_filament
    pred == 1 を元座標 yXXXXX_xXXXXX に戻して重ねる。
"""

from pathlib import Path
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch
from astropy.io import fits

PROJECT_ROOT = Path(
    "/home1/fujiwara.r/PycharmProjects/mmsegmentation_research/mmsegmentation"
)

# =============================================================================
# USER SETTINGS
# 各領域ごとにここを変更
# =============================================================================

REGION_NAME = "ic516"

# 元の柱密度FITS
COLUMN_DENSITY_FITS_PATH = Path(
    "/mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map.fits"
)

# 元画像と同じ縦横サイズの3クラスラベル
# .fits または .npy を使用可能
LABEL_MAP_PATH = Path(
    "/mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map/vis_3class_map/HGBS_ic5146_hires_column_density_map_label3class.fits"
)

# その領域の test 入力パッチ
TEST_PATCH_DIR = Path(
    "/mnt/sp4t/fujiwara.r/split/test/in/ic5146"
)

# 出力画像
OUTPUT_PNG_PATH = Path(
    "/mnt/sp4t/fujiwara.r/split/visualization/"
    "ic5146_confirmed_full_overlay_corrected.png"
)

# AI予測だけの元サイズマスク
OUTPUT_PRED_MASK_PATH = Path(
    "/mnt/sp4t/fujiwara.r/split/visualization/"
    "ic5146_confirmed_full_pred_mask_corrected.png"
)

# =============================================================================
# Model
# =============================================================================

CONFIG_PATH = PROJECT_ROOT / (
    "configs/unet/unet_deeplabv3_hrf_2class_background_vs_confirmed.py"
)

CHECKPOINT_PATH = PROJECT_ROOT / (
    "work_dirs/unet_deeplabv3_hrf_2class_background_vs_confirmed/epoch_6.pth"
)

DEVICE = "cuda:0"

# =============================================================================
# Plot settings
# =============================================================================

# 柱密度表示の percentiles
VMIN_PERCENTILE = 1.0
VMAX_PERCENTILE = 99.5

# ラベルの透明度
LABEL_ALPHA = 0.90

# AI予測の透明度
PRED_ALPHA = 0.75

# line-like labelを少し見やすくするための表示色
CONFIRMED_COLOR = np.array([1.0, 1.0, 0.25])   # yellow
CANDIDATE_COLOR = np.array([0.15, 0.65, 0.90]) # cyan
PRED_COLOR = np.array([255, 60, 180], dtype=np.float32) / 255.0

PATCH_RE = re.compile(r".*y(\d+)_x(\d+)_in\.png$")


def load_fits_2d(path: Path) -> np.ndarray:
    with fits.open(path, memmap=True) as hdul:
        # 最初に2次元データが見つかったHDUを使う
        for hdu in hdul:
            if hdu.data is None:
                continue
            arr = np.asarray(hdu.data)
            arr = np.squeeze(arr)
            if arr.ndim == 2:
                return arr.astype(np.float32)

    raise ValueError(f"No 2D image found in FITS: {path}")


def load_label_map(path: Path) -> np.ndarray:
    suffix = path.suffix.lower()

    if suffix == ".npy":
        arr = np.load(path)
    elif suffix in (".fits", ".fit", ".fts"):
        with fits.open(path, memmap=True) as hdul:
            arr = None
            for hdu in hdul:
                if hdu.data is None:
                    continue
                candidate = np.squeeze(np.asarray(hdu.data))
                if candidate.ndim == 2:
                    arr = candidate
                    break
            if arr is None:
                raise ValueError(f"No 2D label map in {path}")
    else:
        raise ValueError("LABEL_MAP_PATH must be .fits or .npy")

    arr = np.squeeze(np.asarray(arr))
    if arr.ndim != 2:
        raise ValueError(f"Label map must be 2D: {arr.shape}")
    return arr.astype(np.uint8)


def parse_origin(path: Path):
    m = PATCH_RE.match(path.name)
    if m is None:
        raise ValueError(
            "Patch filename must contain yXXXXX_xXXXXX_in.png: "
            f"{path.name}"
        )
    return int(m.group(1)), int(m.group(2))


def extract_prediction(result):
    pred = result.pred_sem_seg.data
    if hasattr(pred, "detach"):
        pred = pred.detach().cpu().numpy()
    pred = np.squeeze(np.asarray(pred))
    if pred.ndim != 2:
        raise ValueError(f"Prediction must be 2D: {pred.shape}")
    return pred.astype(np.uint8)


def main():
    if str(PROJECT_ROOT) not in sys.path:
        sys.path.insert(0, str(PROJECT_ROOT))

    from mmseg.apis import init_model, inference_model

    for p in [
        COLUMN_DENSITY_FITS_PATH,
        LABEL_MAP_PATH,
        TEST_PATCH_DIR,
        CONFIG_PATH,
        CHECKPOINT_PATH,
    ]:
        if not p.exists():
            raise FileNotFoundError(p)

    column_density = load_fits_2d(COLUMN_DENSITY_FITS_PATH)
    label_map = load_label_map(LABEL_MAP_PATH)

    if column_density.shape != label_map.shape:
        raise ValueError(
            "Column-density and label-map shapes differ: "
            f"column={column_density.shape}, label={label_map.shape}"
        )

    full_h, full_w = label_map.shape

    patch_paths = sorted(TEST_PATCH_DIR.rglob("*_in.png"))
    if not patch_paths:
        raise RuntimeError(f"No test patches found in {TEST_PATCH_DIR}")

    print("=" * 92)
    print("Corrected full-region visualization: background vs confirmed_filament")
    print("=" * 92)
    print("Region      :", REGION_NAME)
    print("Map shape   :", label_map.shape)
    print("Num patches :", len(patch_paths))
    print("Checkpoint  :", CHECKPOINT_PATH)

    model = init_model(
        str(CONFIG_PATH),
        str(CHECKPOINT_PATH),
        device=DEVICE,
    )

    model.cfg.test_pipeline = [
        dict(type="LoadImageFromFile"),
        dict(type="PackSegInputs"),
    ]

    # ---------------------------------------------------------
    # 元画像と全く同じshapeでAI予測を復元
    # ---------------------------------------------------------
    full_pred = np.zeros((full_h, full_w), dtype=np.uint8)

    for i, patch_path in enumerate(patch_paths, start=1):
        y0, x0 = parse_origin(patch_path)

        result = inference_model(model, str(patch_path))
        pred = extract_prediction(result)

        ph, pw = pred.shape
        y1 = y0 + ph
        x1 = x0 + pw

        if y0 < 0 or x0 < 0 or y1 > full_h or x1 > full_w:
            raise ValueError(
                f"Patch outside original image: {patch_path}\n"
                f"origin=(y={y0}, x={x0}), patch={pred.shape}, "
                f"full={label_map.shape}"
            )

        region_view = full_pred[y0:y1, x0:x1]
        region_view[pred == 1] = 1

        if i % 50 == 0 or i == len(patch_paths):
            print(f"  processed {i} / {len(patch_paths)}")

    # ---------------------------------------------------------
    # AI予測だけの元サイズPNG
    # ---------------------------------------------------------
    OUTPUT_PRED_MASK_PATH.parent.mkdir(parents=True, exist_ok=True)

    from PIL import Image
    Image.fromarray((full_pred * 255).astype(np.uint8)).save(
        OUTPUT_PRED_MASK_PATH
    )

    # ---------------------------------------------------------
    # 可視化
    # すべて同じ array 座標上で imshow するので位置がずれない
    # ---------------------------------------------------------
    finite_positive = column_density[
        np.isfinite(column_density) & (column_density > 0)
    ]

    if finite_positive.size == 0:
        raise ValueError("Column density contains no positive finite pixels")

    vmin = np.percentile(finite_positive, VMIN_PERCENTILE)
    vmax = np.percentile(finite_positive, VMAX_PERCENTILE)

    if vmax <= vmin:
        vmax = float(finite_positive.max())

    fig, ax = plt.subplots(figsize=(12, 12))

    image = ax.imshow(
        column_density,
        origin="lower",
        cmap="inferno",
        norm=LogNorm(vmin=max(vmin, 1e-30), vmax=vmax),
        interpolation="nearest",
    )

    # confirmed GT
    confirmed_rgba = np.zeros((full_h, full_w, 4), dtype=np.float32)
    confirmed_mask = label_map == 1
    confirmed_rgba[confirmed_mask, :3] = CONFIRMED_COLOR
    confirmed_rgba[confirmed_mask, 3] = LABEL_ALPHA
    ax.imshow(
        confirmed_rgba,
        origin="lower",
        interpolation="nearest",
    )

    # candidate GT
    candidate_rgba = np.zeros((full_h, full_w, 4), dtype=np.float32)
    candidate_mask = label_map == 2
    candidate_rgba[candidate_mask, :3] = CANDIDATE_COLOR
    candidate_rgba[candidate_mask, 3] = LABEL_ALPHA
    ax.imshow(
        candidate_rgba,
        origin="lower",
        interpolation="nearest",
    )

    # AI prediction
    pred_rgba = np.zeros((full_h, full_w, 4), dtype=np.float32)
    pred_mask = full_pred == 1
    pred_rgba[pred_mask, :3] = PRED_COLOR
    pred_rgba[pred_mask, 3] = PRED_ALPHA

    ax.imshow(
        pred_rgba,
        origin="lower",
        interpolation="nearest",
    )

    ax.set_xlim(0, full_w)
    ax.set_ylim(0, full_h)

    ax.set_xlabel("x [pixels]")
    ax.set_ylabel("y [pixels]")
    ax.set_title(
        f"{REGION_NAME} — 3-Class provisional labels + AI confirmed_filament prediction"
    )

    legend_handles = [
        Patch(
            facecolor=CONFIRMED_COLOR,
            label="Provisional confirmed filament",
        ),
        Patch(
            facecolor=CANDIDATE_COLOR,
            label="Provisional candidate filament",
        ),
        Patch(
            facecolor=PRED_COLOR,
            label="AI predicted confirmed_filament",
        ),
    ]

    ax.legend(
        handles=legend_handles,
        loc="upper right",
    )

    cbar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.046,
        pad=0.04,
    )
    cbar.set_label(r"Column Density $N(H_2)$")

    fig.tight_layout()

    OUTPUT_PNG_PATH.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        OUTPUT_PNG_PATH,
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)

    print()
    print("Saved overlay :", OUTPUT_PNG_PATH)
    print("Saved AI mask :", OUTPUT_PRED_MASK_PATH)
    print("Predicted pixels:", int(pred_mask.sum()))
    print("Prediction color: magenta")


if __name__ == "__main__":
    main()
