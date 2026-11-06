#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
visualize_shared_aspp_dual_head_regions.py

HRF事前学習済み DeepLabV3（UNet Backbone）+ 共有ASPP + 最終1x1 Conv 2ヘッド
の学習済み checkpoint から Class 2 / Class 1 を推論し、線可視化を作る。

表示:
    Red    = provisional label
    Cyan   = AI prediction
    Yellow = overlap

重要:
- フィラメントを囲む contour は使わない。
- thick mask を scipy morphological skeleton に変換して、
  「フィラメントの線そのもの」を表示する。
- 線は不透明。背景だけ透明。
- 評価指標は eval threshold (default 0.5) で計算。
- Class 2 が eval threshold では0画素になる場合、
  可視化だけ probability の上位側から threshold を自動で下げる。
  GTは threshold 決定に使わない。
- 図には eval threshold と visualization threshold の両方を表示する。

前提:
- config:
  configs/unet/hrf_deeplabv3_unet_shared_aspp_dual_bce_dice_v2.py
- model:
  mmseg/models/segmentors/hrf_deeplabv3_unet_shared_aspp_dual_head_v2.py
- 出力は class1_prob / class2_prob を使用する。
"""

from pathlib import Path
import argparse
import csv
import json
import re
import sys

import numpy as np
from PIL import Image
from astropy.io import fits
from scipy import ndimage as ndi

try:
    from skimage.morphology import skeletonize
    HAS_SKIMAGE = True
except Exception:
    skeletonize = None
    HAS_SKIMAGE = False

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


# =============================================================================
# DEFAULT PATHS
# =============================================================================

PROJECT_ROOT = Path(
    "/home1/fujiwara.r/PycharmProjects/mmsegmentation_research/mmsegmentation"
)

DATASET_ROOT = Path(
    "/mnt/sp4t/fujiwara.r/split_v2"
)

CONFIG = (
    PROJECT_ROOT
    / "configs/unet/hrf_deeplabv3_unet_shared_aspp_dual_bce_dice_v2.py"
)

CHECKPOINT = (
    PROJECT_ROOT
    / "work_dirs/hrf_deeplabv3_unet_shared_aspp_dual_bce_dice_v2/epoch_17.pth"
)

OUTPUT_ROOT = (
    DATASET_ROOT
    / "shared_aspp_dual_head_visualization_class12"
)

REGIONS = (
    "aquila",
    "orionb",
    "ic5146",
    "taurus",
    "polaris",
)

REGION_FITS = {
    "aquila": Path(
        "/mnt/sp4t/fujiwara.r/HGBS_aquilaM2_hires_column_density_map.fits"
    ),
    "orionb": Path(
        "/mnt/sp4t/fujiwara.r/HGBS_orionB_hires_column_density_map.fits"
    ),
    "ic5146": Path(
        "/mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map.fits"
    ),
    "taurus": Path(
        "/mnt/sp4t/fujiwara.r/HGBS_tauN3_hires_column_density_map.fits"
    ),
    "polaris": Path(
        "/mnt/sp4t/fujiwara.r/HGBS_polaris_hires_column_density_map.fits"
    ),
}

REGION_LABEL_MAP = {
    "aquila": Path(
        "/mnt/sp4t/fujiwara.r/HGBS_aquilaM2_hires_column_density_map/"
        "vis_3class_map/"
        "HGBS_aquilaM2_hires_column_density_map_label3class.fits"
    ),
    "orionb": Path(
        "/mnt/sp4t/fujiwara.r/HGBS_orionB_hires_column_density_map/"
        "vis_3class_map/"
        "HGBS_orionB_hires_column_density_map_label3class.fits"
    ),
    "ic5146": Path(
        "/mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map/"
        "vis_3class_map/"
        "HGBS_ic5146_hires_column_density_map_label3class.fits"
    ),
    "taurus": Path(
        "/mnt/sp4t/fujiwara.r/HGBS_tauN3_hires_column_density_map/"
        "vis_3class_map/"
        "HGBS_tauN3_hires_column_density_map_label3class.fits"
    ),
    "polaris": Path(
        "/mnt/sp4t/fujiwara.r/HGBS_polaris_hires_column_density_map/"
        "vis_3class_map/"
        "HGBS_polaris_hires_column_density_map_label3class.fits"
    ),
}


# =============================================================================
# LABELS / DISPLAY SETTINGS
# =============================================================================

BACKGROUND_LABEL = 0
CONFIRMED_LABEL = 1
CANDIDATE_LABEL = 2
IGNORE_LABEL = 255

FIG_DPI = 240
ZOOM_MARGIN = 180

PATCH_EDGE_WIDTH = 0.30
PATCH_EDGE_ALPHA = 0.60

# 0 = skeletonそのまま（最も細い）
# 1 = 表示だけ1回dilation
LINE_THICKNESS = 0

# AIの小さい孤立成分を可視化だけで削除。
# 0 = 削除しない
AI_MIN_COMPONENT_PIXELS = 0

# ぶちぶち切れた線を減らすための前処理・後処理
# 予測マスクを少し close して細かい隙間を埋める
AI_MASK_CLOSE_ITERS = 1
# skeleton化したあとも 1px 程度の小さな切れ目をつなぐ
AI_LINE_BRIDGE_ITERS = 1
# provisional label は原形を崩しすぎないよう弱めにする
PROVISIONAL_MASK_CLOSE_ITERS = 0
# 線の極小成分を削除（0なら削除しない）
LINE_MIN_COMPONENT_PIXELS = 0

# Class 2が0.5で全く出ない場合の可視化fallback。
# probability上位0.5%程度を表示する。
DEFAULT_CANDIDATE_VIS_PERCENTILE = 99.5

# Class 1も万一0画素の場合。
DEFAULT_CONFIRMED_VIS_PERCENTILE = 99.7

DEFAULT_GALLERY_PATCH_COUNT = 6

# 完全不透明の線色
RED = [1.0, 0.0, 0.0, 1.0]
CYAN = [0.0, 1.0, 1.0, 1.0]
YELLOW = [1.0, 1.0, 0.0, 1.0]

PATCH_RE = re.compile(
    r".*y(\d+)_x(\d+)_in\.png$"
)


# =============================================================================
# ARGUMENTS
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
    )

    parser.add_argument(
        "--dataset-root",
        type=Path,
        default=DATASET_ROOT,
    )

    parser.add_argument(
        "--config",
        type=Path,
        default=CONFIG,
    )

    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=CHECKPOINT,
    )

    parser.add_argument(
        "--output-root",
        type=Path,
        default=OUTPUT_ROOT,
    )

    parser.add_argument(
        "--device",
        type=str,
        default="cuda:0",
    )

    parser.add_argument(
        "--split",
        choices=("train", "val", "test"),
        default="val",
        help="Visualize/evaluate this split. Default: val",
    )

    # 指標用
    parser.add_argument(
        "--candidate-threshold",
        type=float,
        default=0.5,
    )

    parser.add_argument(
        "--confirmed-threshold",
        type=float,
        default=0.5,
    )

    # 可視化専用。Noneなら自動。
    parser.add_argument(
        "--candidate-vis-threshold",
        type=float,
        default=None,
    )

    parser.add_argument(
        "--confirmed-vis-threshold",
        type=float,
        default=None,
    )

    parser.add_argument(
        "--candidate-vis-percentile",
        type=float,
        default=DEFAULT_CANDIDATE_VIS_PERCENTILE,
    )

    parser.add_argument(
        "--confirmed-vis-percentile",
        type=float,
        default=DEFAULT_CONFIRMED_VIS_PERCENTILE,
    )

    parser.add_argument(
        "--gallery-patch-count",
        type=int,
        default=DEFAULT_GALLERY_PATCH_COUNT,
    )

    parser.add_argument(
        "--regions",
        nargs="+",
        default=list(REGIONS),
    )

    parser.add_argument(
        "--no-save-each-test-patch",
        action="store_true",
    )

    parser.add_argument(
        "--no-save-patch-gallery",
        action="store_true",
    )

    args = parser.parse_args()

    for name in (
        "candidate_threshold",
        "confirmed_threshold",
    ):
        v = getattr(args, name)
        if not 0.0 <= v <= 1.0:
            raise ValueError(
                f"{name} must be within [0,1]"
            )

    for name in (
        "candidate_vis_threshold",
        "confirmed_vis_threshold",
    ):
        v = getattr(args, name)
        if (
            v is not None
            and not 0.0 <= v <= 1.0
        ):
            raise ValueError(
                f"{name} must be within [0,1]"
            )

    return args


# =============================================================================
# I/O
# =============================================================================

def load_fits(path):
    with fits.open(
        path,
        memmap=True,
    ) as hdul:
        arr = np.squeeze(
            np.asarray(
                hdul[0].data
            )
        )

    if arr.ndim != 2:
        raise ValueError(
            f"FITS must be 2D: {path}, shape={arr.shape}"
        )

    return arr.astype(
        np.float64
    )


def load_label_map(path):
    with fits.open(
        path,
        memmap=True,
    ) as hdul:
        arr = np.squeeze(
            np.asarray(
                hdul[0].data
            )
        )

    if arr.ndim != 2:
        raise ValueError(
            f"Label map must be 2D: {path}, shape={arr.shape}"
        )

    if np.issubdtype(
        arr.dtype,
        np.floating,
    ):
        arr = np.where(
            np.isfinite(arr),
            np.rint(arr),
            IGNORE_LABEL,
        )

    arr = arr.astype(
        np.uint8
    )

    allowed = {
        0, 1, 2, 255
    }

    bad = (
        set(
            int(v)
            for v in np.unique(arr)
        )
        - allowed
    )

    if bad:
        raise ValueError(
            f"Unexpected labels: {sorted(bad)}"
        )

    return arr


def read_mask(path):
    arr = np.squeeze(
        np.asarray(
            Image.open(path)
        )
    )

    if arr.ndim != 2:
        raise ValueError(
            f"Mask must be 2D: {path}"
        )

    return arr.astype(
        np.uint8
    )


def parse_origin(path):
    m = PATCH_RE.match(
        path.name
    )

    if not m:
        raise ValueError(
            f"Cannot parse patch coordinates: {path.name}"
        )

    return (
        int(m.group(1)),
        int(m.group(2)),
    )


def collect_pairs(
    dataset_root,
    region,
    split,
):
    inp = (
        dataset_root
        / split
        / "in"
        / region
    )

    out = (
        dataset_root
        / split
        / "out"
        / region
    )

    if not inp.is_dir():
        raise FileNotFoundError(
            inp
        )

    if not out.is_dir():
        raise FileNotFoundError(
            out
        )

    pairs = []

    for image_path in sorted(
        inp.rglob(
            "*_in.png"
        )
    ):
        rel = image_path.relative_to(
            inp
        )

        mask_path = (
            out
            / rel.parent
            / rel.name.replace(
                "_in.png",
                "_out.png",
            )
        )

        if not mask_path.is_file():
            raise FileNotFoundError(
                mask_path
            )

        pairs.append(
            (
                image_path,
                mask_path,
            )
        )

    if not pairs:
        raise RuntimeError(
            f"No patches: {region}"
        )

    return pairs


# =============================================================================
# MODEL OUTPUT
# =============================================================================

def _to_numpy(x):
    if (
        hasattr(x, "data")
        and not isinstance(
            x,
            np.ndarray,
        )
    ):
        try:
            x = x.data
        except Exception:
            pass

    if hasattr(
        x,
        "detach",
    ):
        x = (
            x.detach()
            .float()
            .cpu()
            .numpy()
        )

    return np.asarray(x)


def logits_to_prob(logits):
    logits = _to_numpy(
        logits
    )

    if (
        logits.ndim == 4
        and logits.shape[0] == 1
    ):
        logits = logits[0]

    if logits.ndim != 3:
        raise ValueError(
            f"Unexpected logits shape: {logits.shape}"
        )

    if logits.shape[0] == 2:
        z = (
            logits.astype(
                np.float64
            )
        )

        z -= np.max(
            z,
            axis=0,
            keepdims=True,
        )

        e = np.exp(z)

        return (
            e[1]
            / np.sum(
                e,
                axis=0,
            )
        ).astype(
            np.float32
        )

    if logits.shape[0] == 1:
        z = (
            logits[0]
            .astype(
                np.float64
            )
        )

        return (
            1.0
            / (
                1.0
                + np.exp(-z)
            )
        ).astype(
            np.float32
        )

    raise ValueError(
        f"Expected 1 or 2 channels, got {logits.shape}"
    )


def extract_probability(
    result,
    head,
):
    """
    Current model output:
        class1_prob = Class 1 probability
        class2_prob = Class 2 probability

    Backward-compatible aliases are also accepted.
    """
    if head == "confirmed":
        names = (
            "class1_prob",
            "confirmed_prob",
            "confirmed_positive_probability",
            "confirmed_positive_prob",
        )
    elif head == "candidate":
        names = (
            "class2_prob",
            "candidate_prob",
            "candidate_positive_probability",
            "candidate_positive_prob",
        )
    else:
        raise ValueError(head)

    for name in names:
        if hasattr(result, name):
            arr = np.squeeze(
                _to_numpy(
                    getattr(result, name)
                )
            )

            if arr.ndim == 2:
                return arr.astype(np.float32)

            if (
                arr.ndim == 3
                and arr.shape[0] in (1, 2)
            ):
                if arr.shape[0] == 1:
                    return arr[0].astype(np.float32)
                return logits_to_prob(arr)

    if isinstance(result, dict):
        for name in names:
            if name in result:
                arr = np.squeeze(
                    _to_numpy(result[name])
                )

                if arr.ndim == 2:
                    return arr.astype(np.float32)

                if (
                    arr.ndim == 3
                    and arr.shape[0] in (1, 2)
                ):
                    if arr.shape[0] == 1:
                        return arr[0].astype(np.float32)
                    return logits_to_prob(arr)

    attrs = [
        name
        for name in dir(result)
        if (
            "class1" in name
            or "class2" in name
            or "prob" in name
            or "logit" in name
        )
        and not name.startswith("_")
    ]

    raise AttributeError(
        f"{head} probability not found for shared-ASPP dual-head model. "
        f"Available={attrs}"
    )


# =============================================================================
# METRICS
# =============================================================================

def binary_gt(
    mask3,
    positive_label,
):
    out = np.full(
        mask3.shape,
        IGNORE_LABEL,
        dtype=np.uint8,
    )

    out[
        mask3
        == BACKGROUND_LABEL
    ] = 0

    out[
        mask3
        == positive_label
    ] = 1

    return out


def update_conf(
    conf,
    gt,
    pred,
):
    valid = (
        gt
        != IGNORE_LABEL
    )

    encoded = (
        2
        * gt[
            valid
        ].astype(
            np.int64
        )
        + pred[
            valid
        ].astype(
            np.int64
        )
    )

    conf += np.bincount(
        encoded,
        minlength=4,
    ).reshape(
        2,
        2,
    )


def metrics(conf):
    tn = int(
        conf[0, 0]
    )

    fp = int(
        conf[0, 1]
    )

    fn = int(
        conf[1, 0]
    )

    tp = int(
        conf[1, 1]
    )

    def div(a, b):
        if b == 0:
            return float(
                "nan"
            )

        return float(
            a / b
        )

    return {
        "precision": div(
            tp,
            tp + fp,
        ),
        "recall": div(
            tp,
            tp + fn,
        ),
        "iou": div(
            tp,
            tp + fp + fn,
        ),
        "f1_dice": div(
            2 * tp,
            2 * tp + fp + fn,
        ),
        "accuracy": div(
            tp + tn,
            conf.sum(),
        ),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
    }


def fmt_pct(v):
    if not np.isfinite(v):
        return "nan"

    return (
        f"{v * 100:.2f}"
    )


# =============================================================================
# FILAMENT LINE GENERATION
# =============================================================================

def remove_small_components(
    mask,
    min_pixels,
):
    mask = np.asarray(
        mask,
        dtype=bool,
    )

    if (
        min_pixels <= 1
        or not np.any(mask)
    ):
        return mask

    labels, n = ndi.label(
        mask,
        structure=np.ones(
            (3, 3),
            dtype=np.uint8,
        ),
    )

    if n == 0:
        return mask

    sizes = np.bincount(
        labels.ravel()
    )

    keep = np.where(
        sizes
        >= int(
            min_pixels
        )
    )[0]

    keep = keep[
        keep != 0
    ]

    if keep.size == 0:
        return np.zeros_like(
            mask,
            dtype=bool,
        )

    return np.isin(
        labels,
        keep,
    )


def _neighbors_8(img):
    P2 = img[:-2, 1:-1]
    P3 = img[:-2, 2:]
    P4 = img[1:-1, 2:]
    P5 = img[2:, 2:]
    P6 = img[2:, 1:-1]
    P7 = img[2:, :-2]
    P8 = img[1:-1, :-2]
    P9 = img[:-2, :-2]
    return P2, P3, P4, P5, P6, P7, P8, P9


def _zhang_suen_thinning(mask):
    img = np.asarray(mask, dtype=np.uint8).copy()

    if not np.any(img):
        return img.astype(bool)

    changed = True
    while changed:
        changed = False

        P2, P3, P4, P5, P6, P7, P8, P9 = _neighbors_8(img)
        neighbors = P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9
        transitions = ((P2 == 0) & (P3 == 1)).astype(np.uint8)
        transitions += ((P3 == 0) & (P4 == 1)).astype(np.uint8)
        transitions += ((P4 == 0) & (P5 == 1)).astype(np.uint8)
        transitions += ((P5 == 0) & (P6 == 1)).astype(np.uint8)
        transitions += ((P6 == 0) & (P7 == 1)).astype(np.uint8)
        transitions += ((P7 == 0) & (P8 == 1)).astype(np.uint8)
        transitions += ((P8 == 0) & (P9 == 1)).astype(np.uint8)
        transitions += ((P9 == 0) & (P2 == 1)).astype(np.uint8)

        core = img[1:-1, 1:-1] == 1
        m1 = (
            core
            & (neighbors >= 2)
            & (neighbors <= 6)
            & (transitions == 1)
            & ((P2 * P4 * P6) == 0)
            & ((P4 * P6 * P8) == 0)
        )
        if np.any(m1):
            img[1:-1, 1:-1][m1] = 0
            changed = True

        P2, P3, P4, P5, P6, P7, P8, P9 = _neighbors_8(img)
        neighbors = P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9
        transitions = ((P2 == 0) & (P3 == 1)).astype(np.uint8)
        transitions += ((P3 == 0) & (P4 == 1)).astype(np.uint8)
        transitions += ((P4 == 0) & (P5 == 1)).astype(np.uint8)
        transitions += ((P5 == 0) & (P6 == 1)).astype(np.uint8)
        transitions += ((P6 == 0) & (P7 == 1)).astype(np.uint8)
        transitions += ((P7 == 0) & (P8 == 1)).astype(np.uint8)
        transitions += ((P8 == 0) & (P9 == 1)).astype(np.uint8)
        transitions += ((P9 == 0) & (P2 == 1)).astype(np.uint8)

        core = img[1:-1, 1:-1] == 1
        m2 = (
            core
            & (neighbors >= 2)
            & (neighbors <= 6)
            & (transitions == 1)
            & ((P2 * P4 * P8) == 0)
            & ((P2 * P6 * P8) == 0)
        )
        if np.any(m2):
            img[1:-1, 1:-1][m2] = 0
            changed = True

    return img.astype(bool)


def _skeletonize_binary(mask):
    mask = np.asarray(mask, dtype=bool)
    if not np.any(mask):
        return mask
    if HAS_SKIMAGE:
        return np.asarray(skeletonize(mask), dtype=bool)
    return _zhang_suen_thinning(mask)


def _bridge_small_gaps(line_mask, iters=0):
    line_mask = np.asarray(line_mask, dtype=bool)
    if iters <= 0 or not np.any(line_mask):
        return line_mask
    # 小さな切れ目だけをつなぎ、最後に再度細線化する
    bridged = ndi.binary_closing(
        line_mask,
        structure=np.ones((3, 3), dtype=bool),
        iterations=int(iters),
    )
    bridged = ndi.binary_dilation(
        bridged,
        structure=np.ones((3, 3), dtype=bool),
        iterations=1,
    )
    return _skeletonize_binary(bridged)


def _mask_to_line(mask, close_iters=0, bridge_iters=0, min_component_pixels=0):
    mask = np.asarray(mask, dtype=bool)
    if not np.any(mask):
        return mask

    if min_component_pixels > 0:
        mask = remove_small_components(mask, min_component_pixels)

    if close_iters > 0:
        mask = ndi.binary_closing(
            mask,
            structure=np.ones((3, 3), dtype=bool),
            iterations=int(close_iters),
        )

    line = _skeletonize_binary(mask)
    line = _bridge_small_gaps(line, bridge_iters)

    if LINE_THICKNESS > 0:
        line = ndi.binary_dilation(line, iterations=LINE_THICKNESS)

    return np.asarray(line, dtype=bool)


def morphological_skeleton(mask):
    """Backward-compatible name.

    Convert a binary region to a centerline-like thin line.
    """
    return _mask_to_line(mask)


def provisional_line(mask):
    return _mask_to_line(
        mask,
        close_iters=PROVISIONAL_MASK_CLOSE_ITERS,
        bridge_iters=0,
        min_component_pixels=LINE_MIN_COMPONENT_PIXELS,
    )


def ai_line(mask):
    mask = remove_small_components(
        mask,
        AI_MIN_COMPONENT_PIXELS,
    )

    return _mask_to_line(
        mask,
        close_iters=AI_MASK_CLOSE_ITERS,
        bridge_iters=AI_LINE_BRIDGE_ITERS,
        min_component_pixels=LINE_MIN_COMPONENT_PIXELS,
    )


def provisional_overlay(

    candidate_mask,
    confirmed_mask,
):
    candidate = provisional_line(
        candidate_mask
    )

    confirmed = provisional_line(
        confirmed_mask
    )

    overlap = (
        candidate
        & confirmed
    )

    rgba = np.zeros(
        (
            candidate.shape[0],
            candidate.shape[1],
            4,
        ),
        dtype=np.float32,
    )

    rgba[
        candidate
        & ~overlap
    ] = RED

    rgba[
        confirmed
        & ~overlap
    ] = CYAN

    rgba[
        overlap
    ] = YELLOW

    return rgba


def comparison_overlay(
    provisional_mask,
    prediction_mask,
    roi=None,
):
    """Thin-line red/cyan/yellow overlay.

    No contour. No filled filament region.
    Only skeleton pixels are visible.
    """
    provisional_mask = np.asarray(
        provisional_mask,
        dtype=bool,
    ).copy()

    prediction_mask = np.asarray(
        prediction_mask,
        dtype=bool,
    ).copy()

    if roi is not None:
        roi = np.asarray(
            roi,
            dtype=bool,
        )

        provisional_mask &= roi
        prediction_mask &= roi

    red_line = provisional_line(
        provisional_mask
    )

    cyan_line = ai_line(
        prediction_mask
    )

    # AI centerline that falls inside the provisional region is overlap.
    yellow_line = (
        cyan_line
        & provisional_mask
    )

    cyan_only = (
        cyan_line
        & ~provisional_mask
    )

    red_only = (
        red_line
        & ~prediction_mask
    )

    rgba = np.zeros(
        (
            provisional_mask.shape[0],
            provisional_mask.shape[1],
            4,
        ),
        dtype=np.float32,
    )

    rgba[
        red_only
    ] = RED

    rgba[
        cyan_only
    ] = CYAN

    rgba[
        yellow_line
    ] = YELLOW

    return rgba


# =============================================================================
# PLOT HELPERS
# =============================================================================

def display_limits(
    column,
):
    values = column[
        np.isfinite(
            column
        )
        & (
            column > 0
        )
    ]

    lo, hi = np.percentile(
        values,
        [
            1.0,
            99.5,
        ],
    )

    return (
        max(
            float(lo),
            1e-30,
        ),
        float(hi),
    )


def show_column(
    ax,
    column,
    vmin,
    vmax,
):
    return ax.imshow(
        column,
        origin="lower",
        cmap="gray",
        norm=LogNorm(
            vmin=vmin,
            vmax=vmax,
        ),
        interpolation="nearest",
    )


def test_bbox(
    origins,
    shape,
):
    h, w = shape

    y0 = min(
        y
        for y, x, ph, pw
        in origins
    )

    x0 = min(
        x
        for y, x, ph, pw
        in origins
    )

    y1 = max(
        y + ph
        for y, x, ph, pw
        in origins
    )

    x1 = max(
        x + pw
        for y, x, ph, pw
        in origins
    )

    return (
        max(
            0,
            y0 - ZOOM_MARGIN,
        ),
        min(
            h,
            y1 + ZOOM_MARGIN,
        ),
        max(
            0,
            x0 - ZOOM_MARGIN,
        ),
        min(
            w,
            x1 + ZOOM_MARGIN,
        ),
    )


def apply_zoom(
    ax,
    origins,
    shape,
):
    y0, y1, x0, x1 = test_bbox(
        origins,
        shape,
    )

    ax.set_xlim(
        x0,
        x1,
    )

    ax.set_ylim(
        y0,
        y1,
    )


def patch_edges(
    ax,
    origins,
):
    for (
        y,
        x,
        h,
        w,
    ) in origins:
        ax.add_patch(
            Rectangle(
                (
                    x,
                    y,
                ),
                w,
                h,
                fill=False,
                edgecolor="white",
                linewidth=PATCH_EDGE_WIDTH,
                alpha=PATCH_EDGE_ALPHA,
            )
        )


# =============================================================================
# INFERENCE / RECONSTRUCTION
# =============================================================================

def reconstruct_region(
    region,
    model,
    inference_model,
    dataset_root,
    candidate_eval_threshold,
    confirmed_eval_threshold,
    split,
):
    column = load_fits(
        REGION_FITS[
            region
        ]
    )

    full_label = load_label_map(
        REGION_LABEL_MAP[
            region
        ]
    )

    if (
        full_label.shape
        != column.shape
    ):
        raise ValueError(
            f"{region}: column={column.shape}, label={full_label.shape}"
        )

    h, w = column.shape

    test_gt = np.full(
        (
            h,
            w,
        ),
        IGNORE_LABEL,
        dtype=np.uint8,
    )

    candidate_eval = np.zeros(
        (
            h,
            w,
        ),
        dtype=np.uint8,
    )

    confirmed_eval = np.zeros(
        (
            h,
            w,
        ),
        dtype=np.uint8,
    )

    candidate_cm = np.zeros(
        (
            2,
            2,
        ),
        dtype=np.int64,
    )

    confirmed_cm = np.zeros(
        (
            2,
            2,
        ),
        dtype=np.int64,
    )

    origins = []
    records = []

    pairs = collect_pairs(
        dataset_root,
        region,
        split,
    )

    for i, (
        image_path,
        mask_path,
    ) in enumerate(
        pairs,
        1,
    ):
        y, x = parse_origin(
            image_path
        )

        gt = read_mask(
            mask_path
        )

        result = inference_model(
            model,
            str(
                image_path
            ),
        )

        candidate_prob = extract_probability(
            result,
            "candidate",
        )

        confirmed_prob = extract_probability(
            result,
            "confirmed",
        )

        if (
            candidate_prob.shape
            != gt.shape
        ):
            raise ValueError(
                f"Candidate shape mismatch: {image_path}"
            )

        if (
            confirmed_prob.shape
            != gt.shape
        ):
            raise ValueError(
                f"Confirmed shape mismatch: {image_path}"
            )

        cand_eval = (
            candidate_prob
            >= candidate_eval_threshold
        ).astype(
            np.uint8
        )

        conf_eval = (
            confirmed_prob
            >= confirmed_eval_threshold
        ).astype(
            np.uint8
        )

        ph, pw = gt.shape

        test_gt[
            y:y + ph,
            x:x + pw,
        ] = gt

        candidate_eval[
            y:y + ph,
            x:x + pw,
        ] = cand_eval

        confirmed_eval[
            y:y + ph,
            x:x + pw,
        ] = conf_eval

        update_conf(
            candidate_cm,
            binary_gt(
                gt,
                CANDIDATE_LABEL,
            ),
            cand_eval,
        )

        update_conf(
            confirmed_cm,
            binary_gt(
                gt,
                CONFIRMED_LABEL,
            ),
            conf_eval,
        )

        origins.append(
            (
                y,
                x,
                ph,
                pw,
            )
        )

        records.append({
            "image_path": image_path,
            "mask_path": mask_path,
            "y": y,
            "x": x,
            "h": ph,
            "w": pw,
            "candidate_prob": candidate_prob.copy(),
            "confirmed_prob": confirmed_prob.copy(),
            "candidate_eval_pred": cand_eval.copy(),
            "confirmed_eval_pred": conf_eval.copy(),
        })

        if (
            i % 50 == 0
            or i == len(
                pairs
            )
        ):
            print(
                f"  {region}: {i}/{len(pairs)}"
            )

    return {
        "column": column,
        "full_label": full_label,
        "test_gt": test_gt,
        "candidate_eval": candidate_eval,
        "confirmed_eval": confirmed_eval,
        "candidate_cm": candidate_cm,
        "confirmed_cm": confirmed_cm,
        "origins": origins,
        "patch_records": records,
        "patch_count": len(
            records
        ),
    }


# =============================================================================
# VISUALIZATION THRESHOLD
# =============================================================================

def all_probabilities(
    d,
    head,
):
    key = (
        f"{head}_prob"
    )

    return np.concatenate(
        [
            record[
                key
            ].ravel()
            for record
            in d[
                "patch_records"
            ]
        ]
    ).astype(
        np.float32
    )


def choose_vis_threshold(
    d,
    head,
    requested,
    eval_threshold,
    fallback_percentile,
):
    """No GT is used here."""

    probs = all_probabilities(
        d,
        head,
    )

    max_prob = float(
        np.max(
            probs
        )
    )

    if requested is not None:
        return (
            float(
                requested
            ),
            "manual",
            max_prob,
        )

    eval_threshold = float(
        eval_threshold
    )

    if np.any(
        probs
        >= eval_threshold
    ):
        return (
            eval_threshold,
            "evaluation",
            max_prob,
        )

    percentile = float(
        np.clip(
            fallback_percentile,
            0.0,
            100.0,
        )
    )

    threshold = float(
        np.percentile(
            probs,
            percentile,
        )
    )

    if (
        threshold >= max_prob
        and max_prob > 0
    ):
        threshold = float(
            np.nextafter(
                max_prob,
                -np.inf,
            )
        )

    threshold = float(
        np.clip(
            threshold,
            0.0,
            1.0,
        )
    )

    return (
        threshold,
        f"auto_p{percentile:g}",
        max_prob,
    )


def apply_vis_predictions(
    d,
    candidate_threshold,
    confirmed_threshold,
):
    h, w = d[
        "column"
    ].shape

    candidate_full = np.zeros(
        (
            h,
            w,
        ),
        dtype=np.uint8,
    )

    confirmed_full = np.zeros(
        (
            h,
            w,
        ),
        dtype=np.uint8,
    )

    for record in d[
        "patch_records"
    ]:
        cand = (
            record[
                "candidate_prob"
            ]
            >= candidate_threshold
        ).astype(
            np.uint8
        )

        conf = (
            record[
                "confirmed_prob"
            ]
            >= confirmed_threshold
        ).astype(
            np.uint8
        )

        y = record["y"]
        x = record["x"]
        h0 = record["h"]
        w0 = record["w"]

        candidate_full[
            y:y + h0,
            x:x + w0,
        ] = cand

        confirmed_full[
            y:y + h0,
            x:x + w0,
        ] = conf

        record[
            "candidate_vis_pred"
        ] = cand.copy()

        record[
            "confirmed_vis_pred"
        ] = conf.copy()

    d[
        "candidate_vis"
    ] = candidate_full

    d[
        "confirmed_vis"
    ] = confirmed_full

    d[
        "candidate_vis_threshold"
    ] = float(
        candidate_threshold
    )

    d[
        "confirmed_vis_threshold"
    ] = float(
        confirmed_threshold
    )

    return d


# =============================================================================
# FULL REGION FIGURES
# =============================================================================

def save_region_figures(
    region,
    d,
    candidate_metric,
    confirmed_metric,
    output_root,
    candidate_eval_threshold,
    confirmed_eval_threshold,
):
    column = d[
        "column"
    ]

    label = d[
        "full_label"
    ]

    valid_roi = (
        d[
            "test_gt"
        ]
        != IGNORE_LABEL
    )

    cand_gt = (
        label
        == CANDIDATE_LABEL
    )

    conf_gt = (
        label
        == CONFIRMED_LABEL
    )

    cand_pred = (
        d[
            "candidate_vis"
        ]
        == 1
    )

    conf_pred = (
        d[
            "confirmed_vis"
        ]
        == 1
    )

    vmin, vmax = display_limits(
        column
    )

    region_dir = (
        output_root
        / region
    )

    region_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    saved = []

    # 01
    fig, ax = plt.subplots(
        figsize=(
            10,
            10,
        )
    )

    im = show_column(
        ax,
        column,
        vmin,
        vmax,
    )

    patch_edges(
        ax,
        d[
            "origins"
        ],
    )

    apply_zoom(
        ax,
        d[
            "origins"
        ],
        column.shape,
    )

    ax.set_title(
        f"{region}: Original H2 column density\n"
        f"white boxes = selected split patches ({d['patch_count']})"
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
    ).set_label(
        "Column Density"
    )

    fig.tight_layout()

    path = (
        region_dir
        / f"{region}_01_original_test_patches.png"
    )

    fig.savefig(
        path,
        dpi=FIG_DPI,
        bbox_inches="tight",
    )

    plt.close(fig)
    saved.append(path)

    # 02
    fig, ax = plt.subplots(
        figsize=(
            10,
            10,
        )
    )

    im = show_column(
        ax,
        column,
        vmin,
        vmax,
    )

    ax.imshow(
        provisional_overlay(
            cand_gt,
            conf_gt,
        ),
        origin="lower",
        interpolation="nearest",
    )

    patch_edges(
        ax,
        d[
            "origins"
        ],
    )

    apply_zoom(
        ax,
        d[
            "origins"
        ],
        column.shape,
    )

    ax.set_title(
        f"{region}: Provisional labels\n"
        "Red = Class 2, Cyan = Class 1"
    )

    ax.legend(
        handles=[
            Line2D(
                [0], [0],
                color="red",
                lw=1.5,
                label="Class 2 provisional label",
            ),
            Line2D(
                [0], [0],
                color="cyan",
                lw=1.5,
                label="Class 1 provisional label",
            ),
        ],
        loc="upper right",
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
    ).set_label(
        "Column Density"
    )

    fig.tight_layout()

    path = (
        region_dir
        / f"{region}_02_provisional_labels.png"
    )

    fig.savefig(
        path,
        dpi=FIG_DPI,
        bbox_inches="tight",
    )

    plt.close(fig)
    saved.append(path)

    task_rows = [
        (
            "03",
            "Class 2",
            cand_gt,
            cand_pred,
            candidate_metric,
            candidate_eval_threshold,
            d[
                "candidate_vis_threshold"
            ],
        ),
        (
            "04",
            "Class 1",
            conf_gt,
            conf_pred,
            confirmed_metric,
            confirmed_eval_threshold,
            d[
                "confirmed_vis_threshold"
            ],
        ),
    ]

    for (
        number,
        name,
        gt_mask,
        pred_mask,
        metric,
        eval_th,
        vis_th,
    ) in task_rows:

        fig, ax = plt.subplots(
            figsize=(
                10,
                10,
            )
        )

        im = show_column(
            ax,
            column,
            vmin,
            vmax,
        )

        ax.imshow(
            comparison_overlay(
                gt_mask,
                pred_mask,
                valid_roi,
            ),
            origin="lower",
            interpolation="nearest",
        )

        patch_edges(
            ax,
            d[
                "origins"
            ],
        )

        apply_zoom(
            ax,
            d[
                "origins"
            ],
            column.shape,
        )

        ax.set_title(
            f"{region}: {name} filament selected split\n"
            f"eval threshold={eval_th:.3f}, "
            f"visualization threshold={vis_th:.4f}  "
            f"P={fmt_pct(metric['precision'])}%  "
            f"R={fmt_pct(metric['recall'])}%  "
            f"IoU={fmt_pct(metric['iou'])}%  "
            f"F1={fmt_pct(metric['f1_dice'])}%\n"
            "Red = provisional, Cyan = AI prediction, Yellow = overlap"
        )

        ax.legend(
            handles=[
                Line2D(
                    [0], [0],
                    color="red",
                    lw=1.5,
                    label="provisional",
                ),
                Line2D(
                    [0], [0],
                    color="cyan",
                    lw=1.5,
                    label="AI prediction",
                ),
                Line2D(
                    [0], [0],
                    color="yellow",
                    lw=1.5,
                    label="overlap",
                ),
            ],
            loc="upper right",
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
        ).set_label(
            "Column Density"
        )

        fig.tight_layout()

        suffix = (
            "candidate"
            if number == "03"
            else "confirmed"
        )

        path = (
            region_dir
            / f"{region}_{number}_{suffix}_test_prediction.png"
        )

        fig.savefig(
            path,
            dpi=FIG_DPI,
            bbox_inches="tight",
        )

        plt.close(fig)
        saved.append(path)

    return saved


# =============================================================================
# PATCH FIGURES
# =============================================================================

def crop(
    arr,
    record,
):
    y = record["y"]
    x = record["x"]
    h = record["h"]
    w = record["w"]

    return arr[
        y:y + h,
        x:x + w,
    ]


def save_patch(
    region,
    d,
    record,
    index,
    output_root,
):
    column = crop(
        d[
            "column"
        ],
        record,
    )

    label = crop(
        d[
            "full_label"
        ],
        record,
    )

    cand_gt = (
        label
        == CANDIDATE_LABEL
    )

    conf_gt = (
        label
        == CONFIRMED_LABEL
    )

    cand_pred = (
        record[
            "candidate_vis_pred"
        ]
        == 1
    )

    conf_pred = (
        record[
            "confirmed_vis_pred"
        ]
        == 1
    )

    vmin, vmax = display_limits(
        d[
            "column"
        ]
    )

    fig, axes = plt.subplots(
        4,
        1,
        figsize=(
            4.3,
            14,
        ),
    )

    for ax in axes:
        ax.imshow(
            column,
            origin="lower",
            cmap="gray",
            norm=LogNorm(
                vmin=vmin,
                vmax=vmax,
            ),
            interpolation="nearest",
        )

        ax.set_xticks([])
        ax.set_yticks([])

    axes[
        0
    ].set_title(
        "(a) Input"
    )

    axes[
        1
    ].imshow(
        provisional_overlay(
            cand_gt,
            conf_gt,
        ),
        origin="lower",
        interpolation="nearest",
    )

    axes[
        1
    ].set_title(
        "(b) Provisional"
    )

    axes[
        2
    ].imshow(
        comparison_overlay(
            cand_gt,
            cand_pred,
        ),
        origin="lower",
        interpolation="nearest",
    )

    axes[
        2
    ].set_title(
        "(c) Class 2"
    )

    axes[
        3
    ].imshow(
        comparison_overlay(
            conf_gt,
            conf_pred,
        ),
        origin="lower",
        interpolation="nearest",
    )

    axes[
        3
    ].set_title(
        "(d) Class 1"
    )

    fig.suptitle(
        f"{region} TEST patch {index:04d}\n"
        "Red=provisional, Cyan=AI prediction, Yellow=overlap\n"
        f"vis thresholds: candidate={d['candidate_vis_threshold']:.4f}, "
        f"confirmed={d['confirmed_vis_threshold']:.4f}",
        fontsize=10,
    )

    fig.tight_layout(
        rect=(
            0,
            0,
            1,
            0.95,
        )
    )

    patch_dir = (
        output_root
        / region
        / "patch_visualization"
    )

    patch_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    path = (
        patch_dir
        / (
            f"patch_{index:04d}_"
            f"y{record['y']:05d}_"
            f"x{record['x']:05d}.png"
        )
    )

    fig.savefig(
        path,
        dpi=FIG_DPI,
        bbox_inches="tight",
    )

    plt.close(fig)

    return path


# =============================================================================
# GALLERY
# =============================================================================

def choose_gallery(
    d,
    count,
):
    scored = []

    for record in d[
        "patch_records"
    ]:
        label = crop(
            d[
                "full_label"
            ],
            record,
        )

        cand_gt = (
            label
            == CANDIDATE_LABEL
        )

        conf_gt = (
            label
            == CONFIRMED_LABEL
        )

        cand_pred = (
            record[
                "candidate_vis_pred"
            ]
            == 1
        )

        conf_pred = (
            record[
                "confirmed_vis_pred"
            ]
            == 1
        )

        cand_overlap = int(
            np.count_nonzero(
                cand_gt
                & cand_pred
            )
        )

        conf_overlap = int(
            np.count_nonzero(
                conf_gt
                & conf_pred
            )
        )

        cand_pred_n = int(
            np.count_nonzero(
                cand_pred
            )
        )

        conf_pred_n = int(
            np.count_nonzero(
                conf_pred
            )
        )

        cand_gt_n = int(
            np.count_nonzero(
                cand_gt
            )
        )

        conf_gt_n = int(
            np.count_nonzero(
                conf_gt
            )
        )

        # Class 2が見えるpatchを強く優先
        score = (
            40 * cand_overlap
            + 10 * cand_pred_n
            + 8 * conf_overlap
            + 2 * conf_pred_n
            + cand_gt_n
            + conf_gt_n
        )

        scored.append(
            (
                score,
                record,
            )
        )

    scored.sort(
        key=lambda x: x[0],
        reverse=True,
    )

    return [
        record
        for _, record
        in scored[
            :min(
                count,
                len(scored),
            )
        ]
    ]


def save_gallery(
    region,
    d,
    output_root,
    count,
):
    records = choose_gallery(
        d,
        count,
    )

    if not records:
        return None

    ncols = len(
        records
    )

    fig, axes = plt.subplots(
        4,
        ncols,
        figsize=(
            3.0 * ncols,
            11.0,
        ),
        squeeze=False,
    )

    vmin, vmax = display_limits(
        d[
            "column"
        ]
    )

    for col, record in enumerate(
        records
    ):
        column = crop(
            d[
                "column"
            ],
            record,
        )

        label = crop(
            d[
                "full_label"
            ],
            record,
        )

        cand_gt = (
            label
            == CANDIDATE_LABEL
        )

        conf_gt = (
            label
            == CONFIRMED_LABEL
        )

        cand_pred = (
            record[
                "candidate_vis_pred"
            ]
            == 1
        )

        conf_pred = (
            record[
                "confirmed_vis_pred"
            ]
            == 1
        )

        panels = [
            None,
            provisional_overlay(
                cand_gt,
                conf_gt,
            ),
            comparison_overlay(
                cand_gt,
                cand_pred,
            ),
            comparison_overlay(
                conf_gt,
                conf_pred,
            ),
        ]

        for row in range(
            4
        ):
            ax = axes[
                row,
                col,
            ]

            ax.imshow(
                column,
                origin="lower",
                cmap="gray",
                norm=LogNorm(
                    vmin=vmin,
                    vmax=vmax,
                ),
                interpolation="nearest",
            )

            if (
                panels[
                    row
                ]
                is not None
            ):
                ax.imshow(
                    panels[
                        row
                    ],
                    origin="lower",
                    interpolation="nearest",
                )

            ax.set_xticks([])
            ax.set_yticks([])

        axes[
            0,
            col,
        ].set_title(
            f"y={record['y']}, x={record['x']}",
            fontsize=9,
        )

    labels = [
        "(a) Input",
        "(b) Provisional",
        "(c) Class 2",
        "(d) Class 1",
    ]

    for row, label in enumerate(
        labels
    ):
        axes[
            row,
            0,
        ].set_ylabel(
            label,
            fontsize=10,
        )

    fig.suptitle(
        f"{region}: representative selected-split patches\n"
        "Red=provisional, Cyan=AI prediction, Yellow=overlap\n"
        f"visualization thresholds: "
        f"candidate={d['candidate_vis_threshold']:.4f}, "
        f"confirmed={d['confirmed_vis_threshold']:.4f}",
        fontsize=11,
    )

    fig.tight_layout(
        rect=(
            0,
            0,
            1,
            0.95,
        )
    )

    region_dir = (
        output_root
        / region
    )

    region_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    path = (
        region_dir
        / f"{region}_05_patch_gallery.png"
    )

    fig.savefig(
        path,
        dpi=FIG_DPI,
        bbox_inches="tight",
    )

    plt.close(fig)

    return path


# =============================================================================
# SAVE TABLES
# =============================================================================

def save_metrics(
    rows,
    output_root,
):
    fields = [
        "region",
        "task",
        "patch_count",
        "threshold",
        "precision",
        "recall",
        "iou",
        "f1_dice",
        "accuracy",
        "tp",
        "fp",
        "fn",
        "tn",
    ]

    csv_path = (
        output_root
        / "shared_aspp_dual_head_metrics_all_regions.csv"
    )

    json_path = (
        output_root
        / "shared_aspp_dual_head_metrics_all_regions.json"
    )

    with csv_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fields,
        )
        writer.writeheader()
        writer.writerows(
            rows
        )

    with json_path.open(
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(
            rows,
            f,
            ensure_ascii=False,
            indent=2,
            allow_nan=True,
        )

    return (
        csv_path,
        json_path,
    )


def save_thresholds(
    rows,
    output_root,
):
    path = (
        output_root
        / "visualization_thresholds.csv"
    )

    fields = [
        "region",
        "candidate_eval_threshold",
        "candidate_vis_threshold",
        "candidate_vis_source",
        "candidate_max_probability",
        "candidate_vis_pixels",
        "confirmed_eval_threshold",
        "confirmed_vis_threshold",
        "confirmed_vis_source",
        "confirmed_max_probability",
        "confirmed_vis_pixels",
    ]

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fields,
        )

        writer.writeheader()
        writer.writerows(
            rows
        )

    return path


# =============================================================================
# SANITY CHECK
# =============================================================================

def verify_model(
    model,
):
    expected = (
        "decode_head",
        "class1_classifier",
        "class2_classifier",
    )

    missing = [
        name
        for name in expected
        if not hasattr(
            model,
            name,
        )
    ]

    if missing:
        raise RuntimeError(
            "Wrong model implementation loaded. "
            f"Missing: {missing}. "
            "Use hrf_deeplabv3_unet_shared_aspp_dual_head_v2.py "
            "with its matching checkpoint."
        )


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()

    if (
        str(
            args.project_root
        )
        not in sys.path
    ):
        sys.path.insert(
            0,
            str(
                args.project_root
            ),
        )

    from mmengine.config import Config
    from mmseg.apis import (
        init_model,
        inference_model,
    )

    if not args.config.is_file():
        raise FileNotFoundError(
            args.config
        )

    if not args.checkpoint.is_file():
        raise FileNotFoundError(
            args.checkpoint
        )

    for region in args.regions:
        if region not in REGION_FITS:
            raise KeyError(
                region
            )

        if not REGION_FITS[
            region
        ].is_file():
            raise FileNotFoundError(
                REGION_FITS[
                    region
                ]
            )

        if not REGION_LABEL_MAP[
            region
        ].is_file():
            raise FileNotFoundError(
                REGION_LABEL_MAP[
                    region
                ]
            )

    args.output_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    cfg = Config.fromfile(
        str(
            args.config
        )
    )

    print(
        "Config:",
        args.config,
    )

    print(
        "Checkpoint:",
        args.checkpoint,
    )

    model = init_model(
        cfg,
        str(
            args.checkpoint
        ),
        device=args.device,
    )

    verify_model(
        model
    )

    model.cfg.test_pipeline = [
        dict(
            type="LoadImageFromFile"
        ),
        dict(
            type="PackSegInputs"
        ),
    ]

    metric_rows = []
    threshold_rows = []

    all_candidate = np.zeros(
        (
            2,
            2,
        ),
        dtype=np.int64,
    )

    all_confirmed = np.zeros(
        (
            2,
            2,
        ),
        dtype=np.int64,
    )

    for region in args.regions:
        print(
            "\n"
            + "=" * 90
        )

        print(
            "Region:",
            region,
        )

        d = reconstruct_region(
            region,
            model,
            inference_model,
            args.dataset_root,
            args.candidate_threshold,
            args.confirmed_threshold,
            args.split,
        )

        candidate_metric = metrics(
            d[
                "candidate_cm"
            ]
        )

        confirmed_metric = metrics(
            d[
                "confirmed_cm"
            ]
        )

        all_candidate += d[
            "candidate_cm"
        ]

        all_confirmed += d[
            "confirmed_cm"
        ]

        (
            candidate_vis_threshold,
            candidate_vis_source,
            candidate_max,
        ) = choose_vis_threshold(
            d,
            "candidate",
            args.candidate_vis_threshold,
            args.candidate_threshold,
            args.candidate_vis_percentile,
        )

        (
            confirmed_vis_threshold,
            confirmed_vis_source,
            confirmed_max,
        ) = choose_vis_threshold(
            d,
            "confirmed",
            args.confirmed_vis_threshold,
            args.confirmed_threshold,
            args.confirmed_vis_percentile,
        )

        d = apply_vis_predictions(
            d,
            candidate_vis_threshold,
            confirmed_vis_threshold,
        )

        candidate_vis_pixels = int(
            np.count_nonzero(
                d[
                    "candidate_vis"
                ]
            )
        )

        confirmed_vis_pixels = int(
            np.count_nonzero(
                d[
                    "confirmed_vis"
                ]
            )
        )

        print(
            f"Class 2 eval_th={args.candidate_threshold:.3f}, "
            f"vis_th={candidate_vis_threshold:.6f} "
            f"({candidate_vis_source}), "
            f"max_prob={candidate_max:.6f}, "
            f"vis_pixels={candidate_vis_pixels}"
        )

        print(
            f"Class 1 eval_th={args.confirmed_threshold:.3f}, "
            f"vis_th={confirmed_vis_threshold:.6f} "
            f"({confirmed_vis_source}), "
            f"max_prob={confirmed_max:.6f}, "
            f"vis_pixels={confirmed_vis_pixels}"
        )

        threshold_rows.append({
            "region": region,
            "candidate_eval_threshold": args.candidate_threshold,
            "candidate_vis_threshold": candidate_vis_threshold,
            "candidate_vis_source": candidate_vis_source,
            "candidate_max_probability": candidate_max,
            "candidate_vis_pixels": candidate_vis_pixels,
            "confirmed_eval_threshold": args.confirmed_threshold,
            "confirmed_vis_threshold": confirmed_vis_threshold,
            "confirmed_vis_source": confirmed_vis_source,
            "confirmed_max_probability": confirmed_max,
            "confirmed_vis_pixels": confirmed_vis_pixels,
        })

        metric_rows.extend([
            {
                "region": region,
                "task": "candidate",
                "patch_count": d["patch_count"],
                "threshold": args.candidate_threshold,
                **candidate_metric,
            },
            {
                "region": region,
                "task": "confirmed",
                "patch_count": d["patch_count"],
                "threshold": args.confirmed_threshold,
                **confirmed_metric,
            },
        ])

        for path in save_region_figures(
            region,
            d,
            candidate_metric,
            confirmed_metric,
            args.output_root,
            args.candidate_threshold,
            args.confirmed_threshold,
        ):
            print(
                "Saved:",
                path,
            )

        if not args.no_save_each_test_patch:
            for i, record in enumerate(
                d[
                    "patch_records"
                ],
                1,
            ):
                path = save_patch(
                    region,
                    d,
                    record,
                    i,
                    args.output_root,
                )

                if (
                    i == 1
                    or i % 50 == 0
                    or i == d[
                        "patch_count"
                    ]
                ):
                    print(
                        "Saved patch:",
                        path,
                    )

        if not args.no_save_patch_gallery:
            path = save_gallery(
                region,
                d,
                args.output_root,
                args.gallery_patch_count,
            )

            if path is not None:
                print(
                    "Saved gallery:",
                    path,
                )

    candidate_all = metrics(
        all_candidate
    )

    confirmed_all = metrics(
        all_confirmed
    )

    total_patches = sum(
        row[
            "patch_count"
        ]
        for row
        in metric_rows
        if (
            row[
                "task"
            ]
            == "candidate"
        )
    )

    metric_rows.extend([
        {
            "region": "ALL_REGIONS",
            "task": "candidate",
            "patch_count": total_patches,
            "threshold": args.candidate_threshold,
            **candidate_all,
        },
        {
            "region": "ALL_REGIONS",
            "task": "confirmed",
            "patch_count": total_patches,
            "threshold": args.confirmed_threshold,
            **confirmed_all,
        },
    ])

    metrics_csv, metrics_json = save_metrics(
        metric_rows,
        args.output_root,
    )

    thresholds_csv = save_thresholds(
        threshold_rows,
        args.output_root,
    )

    print(
        "\nALL candidate:",
        candidate_all,
    )

    print(
        "ALL confirmed:",
        confirmed_all,
    )

    print(
        "Saved:",
        metrics_csv,
    )

    print(
        "Saved:",
        metrics_json,
    )

    print(
        "Saved:",
        thresholds_csv,
    )

    print(
        "Output:",
        args.output_root,
    )


if __name__ == "__main__":
    main()
