#!/usr/bin/env python3
from __future__ import annotations

"""
test/in と test/out の256x256パッチを元座標へ戻して、
testされた場所・GT・AI予測を可視化する。

領域全体ラベルは使わない。
そのため、柱密度FITSと全体ラベルのshape不一致問題を回避できる。

出力:
  1) confirmed_test_overview.png
  2) confirmed_test_zoom.png
  3) test_patch_details/*.png
  4) confirmed_pred_full_mask.png
  5) confirmed_gt_full_mask.png
"""

from pathlib import Path
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch, Rectangle
from astropy.io import fits
from PIL import Image

PROJECT_ROOT = Path(
    "/home1/fujiwara.r/PycharmProjects/mmsegmentation_research/mmsegmentation"
)

# =============================================================================
# USER SETTINGS: 領域ごとにここを変更
# =============================================================================

REGION_NAME = "Aquila"

COLUMN_DENSITY_FITS_PATH = Path(
    "/mnt/sp4t/fujiwara.r/HGBS_aquilaM2_hires_column_density_map/"
    "HGBS_aquilaM2_hires_column_density_map.fits"
)

TEST_INPUT_DIR = Path(
    "/mnt/sp4t/fujiwara.r/split/test/in/aquila"
)

TEST_MASK_DIR = Path(
    "/mnt/sp4t/fujiwara.r/split/test/out/aquila"
)

OUTPUT_ROOT = Path(
    "/mnt/sp4t/fujiwara.r/split/visualization/"
    "test_detail_aquila_confirmed"
)

ZOOM_MARGIN = 256

# =============================================================================
# MODEL
# =============================================================================

CONFIG_PATH = PROJECT_ROOT / (
    "configs/unet/unet_deeplabv3_hrf_2class_background_vs_confirmed.py"
)

CHECKPOINT_PATH = PROJECT_ROOT / (
    "work_dirs/unet_deeplabv3_hrf_2class_background_vs_confirmed/epoch_6.pth"
)

DEVICE = "cuda:0"

# =============================================================================
# COLORS
# =============================================================================

BASE_CMAP = "gray"

GT_COLOR = np.array(
    [255, 215, 0], dtype=np.float32
) / 255.0

OTHER_CLASS_COLOR = np.array(
    [60, 220, 255], dtype=np.float32
) / 255.0

PRED_COLOR = np.array(
    [255, 0, 180], dtype=np.float32
) / 255.0

TEST_AREA_COLOR = np.array(
    [70, 130, 255], dtype=np.float32
) / 255.0

GT_ALPHA = 0.90
OTHER_ALPHA = 0.45
PRED_ALPHA = 0.80
TEST_AREA_ALPHA = 0.15

PATCH_RE = re.compile(
    r".*y(\d+)_x(\d+)_in\.png$"
)


def load_fits_2d(path):
    with fits.open(path, memmap=True) as hdul:
        for hdu in hdul:
            if hdu.data is None:
                continue

            arr = np.squeeze(
                np.asarray(hdu.data)
            )

            if arr.ndim == 2:
                return arr.astype(
                    np.float32
                )

    raise ValueError(
        f"No 2D data in FITS: {path}"
    )


def parse_origin(path):
    m = PATCH_RE.match(path.name)

    if m is None:
        raise ValueError(
            "Could not parse patch coordinates: "
            f"{path.name}"
        )

    return (
        int(m.group(1)),
        int(m.group(2)),
    )


def find_matching_mask(input_path):
    """
    test/in/aquila/0000/xxx_in.png
                    ↓
    test/out/aquila/0000/xxx_out.png
    """
    rel = input_path.relative_to(
        TEST_INPUT_DIR
    )

    mask_name = rel.name.replace(
        "_in.png",
        "_out.png",
    )

    mask_path = (
        TEST_MASK_DIR
        / rel.parent
        / mask_name
    )

    if not mask_path.is_file():
        raise FileNotFoundError(
            "Matching mask not found:\n"
            f"input={input_path}\n"
            f"mask={mask_path}"
        )

    return mask_path


def read_mask(path):
    arr = np.squeeze(
        np.asarray(
            Image.open(path)
        )
    )

    if arr.ndim != 2:
        raise ValueError(
            f"Mask must be 2D: {path}, "
            f"shape={arr.shape}"
        )

    return arr.astype(
        np.uint8
    )


def extract_prediction(result):
    pred = result.pred_sem_seg.data

    if hasattr(pred, "detach"):
        pred = (
            pred.detach()
            .cpu()
            .numpy()
        )

    pred = np.squeeze(
        np.asarray(pred)
    )

    if pred.ndim != 2:
        raise ValueError(
            f"Prediction must be 2D: "
            f"{pred.shape}"
        )

    return pred.astype(
        np.uint8
    )


def make_temp_config(path):
    text = path.read_text(
        encoding="utf-8"
    )

    text = re.sub(
        r"custom_imports\s*=\s*dict\([^\)]*\)",
        (
            "custom_imports = "
            "dict(imports=['mmseg.models'], "
            "allow_failed_imports=False)"
        ),
        text,
        count=1,
        flags=re.DOTALL,
    )

    temp = path.parent / (
        "_tmp_infer_" + path.name
    )

    temp.write_text(
        text,
        encoding="utf-8",
    )

    return temp


def normalized_gray(arr):
    finite = arr[
        np.isfinite(arr)
        & (arr > 0)
    ]

    if finite.size == 0:
        return np.zeros(
            arr.shape,
            dtype=np.float32,
        )

    lo = np.percentile(
        finite, 1.0
    )
    hi = np.percentile(
        finite, 99.5
    )

    if hi <= lo:
        hi = float(
            finite.max()
        )

    clipped = np.clip(
        arr,
        max(lo, 1e-30),
        hi,
    )

    disp = np.log10(
        clipped
    )

    mn = np.nanmin(disp)
    mx = np.nanmax(disp)

    return (
        (disp - mn)
        / max(mx - mn, 1e-8)
    )


def make_rgb(
    base,
    target_mask,
    other_mask,
    pred_mask,
):
    gray = normalized_gray(
        base
    )

    rgb = np.stack(
        [gray, gray, gray],
        axis=-1,
    ).astype(np.float32)

    rgb[other_mask] = (
        (1 - OTHER_ALPHA)
        * rgb[other_mask]
        + OTHER_ALPHA
        * OTHER_CLASS_COLOR
    )

    rgb[target_mask] = (
        (1 - GT_ALPHA)
        * rgb[target_mask]
        + GT_ALPHA
        * GT_COLOR
    )

    rgb[pred_mask] = (
        (1 - PRED_ALPHA)
        * rgb[pred_mask]
        + PRED_ALPHA
        * PRED_COLOR
    )

    return np.clip(
        rgb, 0, 1
    )


def save_patch_detail(
    index,
    patch_path,
    column_patch,
    mask_patch,
    pred_patch,
    output_dir,
):
    target_gt = (
        mask_patch == 1
    )

    other_gt = (
        (mask_patch > 0)
        & (mask_patch != 1)
        & (mask_patch != 255)
    )

    pred_fg = (
        pred_patch == 1
    )

    gt_rgb = make_rgb(
        column_patch,
        target_gt,
        other_gt,
        np.zeros_like(
            pred_fg,
            dtype=bool,
        ),
    )

    result_rgb = make_rgb(
        column_patch,
        target_gt,
        other_gt,
        pred_fg,
    )

    fig, axes = plt.subplots(
        1, 3,
        figsize=(12, 4),
    )

    finite = column_patch[
        np.isfinite(column_patch)
        & (column_patch > 0)
    ]

    if finite.size:
        lo = np.percentile(
            finite, 1
        )
        hi = np.percentile(
            finite, 99.5
        )
        if hi <= lo:
            hi = finite.max()
    else:
        lo, hi = 1, 2

    axes[0].imshow(
        column_patch,
        origin="lower",
        cmap=BASE_CMAP,
        norm=LogNorm(
            vmin=max(lo, 1e-30),
            vmax=hi,
        ),
    )

    axes[0].set_title(
        "Input"
    )

    axes[1].imshow(
        gt_rgb,
        origin="lower",
    )

    axes[1].set_title(
        "GT"
    )

    axes[2].imshow(
        result_rgb,
        origin="lower",
    )

    axes[2].set_title(
        "GT + Prediction"
    )

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(
        patch_path.name
    )

    fig.tight_layout()

    out = output_dir / (
        f"{index:04d}_"
        f"{patch_path.stem}.png"
    )

    fig.savefig(
        out,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)


def add_overlay_layer(
    ax,
    mask,
    color,
    alpha,
):
    rgba = np.zeros(
        (
            mask.shape[0],
            mask.shape[1],
            4,
        ),
        dtype=np.float32,
    )

    rgba[mask, :3] = color
    rgba[mask, 3] = alpha

    ax.imshow(
        rgba,
        origin="lower",
        interpolation="nearest",
    )


def main():
    if str(PROJECT_ROOT) not in sys.path:
        sys.path.insert(
            0,
            str(PROJECT_ROOT),
        )

    from mmseg.apis import (
        init_model,
        inference_model,
    )

    column = load_fits_2d(
        COLUMN_DENSITY_FITS_PATH
    )

    full_h, full_w = (
        column.shape
    )

    patch_paths = sorted(
        TEST_INPUT_DIR.rglob(
            "*_in.png"
        )
    )

    if not patch_paths:
        raise RuntimeError(
            f"No test patches: "
            f"{TEST_INPUT_DIR}"
        )

    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    detail_dir = (
        OUTPUT_ROOT
        / "test_patch_details"
    )

    detail_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # 255 = 未評価領域
    full_gt = np.full(
        (full_h, full_w),
        255,
        dtype=np.uint8,
    )

    full_pred = np.zeros(
        (full_h, full_w),
        dtype=np.uint8,
    )

    test_coverage = np.zeros(
        (full_h, full_w),
        dtype=np.uint8,
    )

    temp_config = (
        make_temp_config(
            CONFIG_PATH
        )
    )

    model = init_model(
        str(temp_config),
        str(CHECKPOINT_PATH),
        device=DEVICE,
    )

    model.cfg.test_pipeline = [
        dict(
            type="LoadImageFromFile"
        ),
        dict(
            type="PackSegInputs"
        ),
    ]

    ymins = []
    ymaxs = []
    xmins = []
    xmaxs = []

    print(
        "Column shape:",
        column.shape,
    )

    print(
        "Test patches:",
        len(patch_paths),
    )

    for i, patch_path in enumerate(
        patch_paths,
        start=1,
    ):
        y0, x0 = parse_origin(
            patch_path
        )

        mask_path = (
            find_matching_mask(
                patch_path
            )
        )

        gt_patch = read_mask(
            mask_path
        )

        result = inference_model(
            model,
            str(patch_path),
        )

        pred_patch = (
            extract_prediction(
                result
            )
        )

        ph, pw = (
            pred_patch.shape
        )

        y1 = y0 + ph
        x1 = x0 + pw

        if (
            y1 > full_h
            or x1 > full_w
        ):
            raise ValueError(
                "Patch outside FITS canvas:\n"
                f"{patch_path}\n"
                f"end=(y={y1},x={x1}), "
                f"full={column.shape}"
            )

        if (
            gt_patch.shape
            != pred_patch.shape
        ):
            raise ValueError(
                f"GT/pred shape mismatch: "
                f"gt={gt_patch.shape}, "
                f"pred={pred_patch.shape}"
            )

        full_gt[
            y0:y1,
            x0:x1
        ] = gt_patch

        full_pred[
            y0:y1,
            x0:x1
        ][
            pred_patch == 1
        ] = 1

        test_coverage[
            y0:y1,
            x0:x1
        ] = 1

        ymins.append(y0)
        ymaxs.append(y1)
        xmins.append(x0)
        xmaxs.append(x1)

        column_patch = column[
            y0:y1,
            x0:x1
        ]

        save_patch_detail(
            i,
            patch_path,
            column_patch,
            gt_patch,
            pred_patch,
            detail_dir,
        )

        if (
            i % 25 == 0
            or i == len(patch_paths)
        ):
            print(
                f"processed "
                f"{i}/"
                f"{len(patch_paths)}"
            )

    Image.fromarray(
        (
            full_pred * 255
        ).astype(np.uint8)
    ).save(
        OUTPUT_ROOT
        / "confirmed_pred_full_mask.png"
    )

    Image.fromarray(
        full_gt
    ).save(
        OUTPUT_ROOT
        / "confirmed_gt_full_mask.png"
    )

    finite = column[
        np.isfinite(column)
        & (column > 0)
    ]

    vmin = np.percentile(
        finite, 1
    )

    vmax = np.percentile(
        finite, 99.5
    )

    target_mask = (
        full_gt == 1
    )

    other_mask = (
        (full_gt > 0)
        & (full_gt != 1)
        & (full_gt != 255)
    )

    pred_mask = (
        full_pred == 1
    )

    test_mask = (
        test_coverage == 1
    )

    legend_handles = [
        Patch(
            facecolor=TEST_AREA_COLOR,
            alpha=TEST_AREA_ALPHA,
            label="test area",
        ),
        Patch(
            facecolor=GT_COLOR,
            alpha=GT_ALPHA,
            label="confirmed filament GT",
        ),
        Patch(
            facecolor=OTHER_CLASS_COLOR,
            alpha=OTHER_ALPHA,
            label="other class GT",
        ),
        Patch(
            facecolor=PRED_COLOR,
            alpha=PRED_ALPHA,
            label="AI prediction",
        ),
    ]

    # =========================================================
    # FULL OVERVIEW
    # =========================================================

    fig, ax = plt.subplots(
        figsize=(11, 11)
    )

    im = ax.imshow(
        column,
        origin="lower",
        cmap=BASE_CMAP,
        norm=LogNorm(
            vmin=max(vmin, 1e-30),
            vmax=vmax,
        ),
        interpolation="nearest",
    )

    add_overlay_layer(
        ax,
        test_mask,
        TEST_AREA_COLOR,
        TEST_AREA_ALPHA,
    )

    add_overlay_layer(
        ax,
        other_mask,
        OTHER_CLASS_COLOR,
        OTHER_ALPHA,
    )

    add_overlay_layer(
        ax,
        target_mask,
        GT_COLOR,
        GT_ALPHA,
    )

    add_overlay_layer(
        ax,
        pred_mask,
        PRED_COLOR,
        PRED_ALPHA,
    )

    for p in patch_paths:
        py, px = parse_origin(p)

        ax.add_patch(
            Rectangle(
                (px, py),
                256,
                256,
                linewidth=0.8,
                edgecolor="white",
                facecolor="none",
                alpha=0.85,
            )
        )

    ax.legend(
        handles=legend_handles,
        loc="upper right",
    )

    ax.set_title(
        f"{REGION_NAME} "
        f"— test overview "
        f"(confirmed filament)"
    )

    ax.set_xlabel(
        "x [pixels]"
    )

    ax.set_ylabel(
        "y [pixels]"
    )

    fig.colorbar(
        im,
        ax=ax,
        fraction=0.046,
        pad=0.04,
    )

    fig.tight_layout()

    fig.savefig(
        OUTPUT_ROOT
        / "confirmed_test_overview.png",
        dpi=220,
        bbox_inches="tight",
    )

    plt.close(fig)

    # =========================================================
    # ZOOM TEST AREA
    # =========================================================

    zy0 = max(
        min(ymins) - ZOOM_MARGIN,
        0,
    )

    zy1 = min(
        max(ymaxs) + ZOOM_MARGIN,
        full_h,
    )

    zx0 = max(
        min(xmins) - ZOOM_MARGIN,
        0,
    )

    zx1 = min(
        max(xmaxs) + ZOOM_MARGIN,
        full_w,
    )

    crop_column = column[
        zy0:zy1,
        zx0:zx1,
    ]

    crop_target = target_mask[
        zy0:zy1,
        zx0:zx1,
    ]

    crop_other = other_mask[
        zy0:zy1,
        zx0:zx1,
    ]

    crop_pred = pred_mask[
        zy0:zy1,
        zx0:zx1,
    ]

    crop_test = test_mask[
        zy0:zy1,
        zx0:zx1,
    ]

    fig, ax = plt.subplots(
        figsize=(12, 12)
    )

    ax.imshow(
        crop_column,
        origin="lower",
        cmap=BASE_CMAP,
        norm=LogNorm(
            vmin=max(vmin, 1e-30),
            vmax=vmax,
        ),
        interpolation="nearest",
    )

    add_overlay_layer(
        ax,
        crop_test,
        TEST_AREA_COLOR,
        TEST_AREA_ALPHA,
    )

    add_overlay_layer(
        ax,
        crop_other,
        OTHER_CLASS_COLOR,
        OTHER_ALPHA,
    )

    add_overlay_layer(
        ax,
        crop_target,
        GT_COLOR,
        GT_ALPHA,
    )

    add_overlay_layer(
        ax,
        crop_pred,
        PRED_COLOR,
        PRED_ALPHA,
    )

    for p in patch_paths:
        py, px = parse_origin(p)

        if (
            px + 256 <= zx0
            or px >= zx1
            or py + 256 <= zy0
            or py >= zy1
        ):
            continue

        ax.add_patch(
            Rectangle(
                (
                    px - zx0,
                    py - zy0,
                ),
                256,
                256,
                linewidth=1.2,
                edgecolor="white",
                facecolor="none",
                alpha=0.9,
            )
        )

    ax.legend(
        handles=legend_handles,
        loc="upper right",
    )

    ax.set_title(
        f"{REGION_NAME} "
        f"— enlarged test area "
        f"(confirmed filament)"
    )

    fig.tight_layout()

    fig.savefig(
        OUTPUT_ROOT
        / "confirmed_test_zoom.png",
        dpi=240,
        bbox_inches="tight",
    )

    plt.close(fig)

    try:
        temp_config.unlink()
    except Exception:
        pass

    print()
    print("Saved:")
    print(
        OUTPUT_ROOT
        / "confirmed_test_overview.png"
    )
    print(
        OUTPUT_ROOT
        / "confirmed_test_zoom.png"
    )
    print(
        detail_dir
    )


if __name__ == "__main__":
    main()
