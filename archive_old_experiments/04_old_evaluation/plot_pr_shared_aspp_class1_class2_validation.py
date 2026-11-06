#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_shared_aspp_dual_head_pr_curves_val.py

HRF事前学習済み DeepLabV3（UNet Backbone）+ 共有ASPP + 最終1x1 Conv 2ヘッド
の学習済み checkpoint を使って、split_v2 VALIDATION パッチに対する
Class 1 / Class 2 のPR曲線を領域ごとに作成するスクリプト。

前提:
- 学習時 config が
    configs/unet/hrf_deeplabv3_unet_shared_aspp_dual_bce_dice_v2.py
  であること
- モデル本体が
    mmseg/models/segmentors/hrf_deeplabv3_unet_shared_aspp_dual_head_v2.py
  であること
- validation PR曲線から threshold を選び、testでthreshold調整しないこと

出力:
  OUTPUT_DIR/
    class2_pr_curve_5regions_validation.png
    class1_pr_curve_5regions_validation.png
    {region}_class1_vs_class2_pr.png
    shared_aspp_dual_head_pr_summary_validation.csv

使い方例:
python tools/analysis_tools/plot_shared_aspp_dual_head_pr_curves_val.py
python tools/analysis_tools/plot_shared_aspp_dual_head_pr_curves_val.py --device cpu
python tools/analysis_tools/plot_shared_aspp_dual_head_pr_curves_val.py \
  --checkpoint /path/to/epoch_10.pth
"""

from pathlib import Path
import argparse
import csv
import sys

import numpy as np
from PIL import Image

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =============================================================================
# DEFAULT SETTINGS
# =============================================================================

DEFAULT_PROJECT_ROOT = Path(
    "/home1/fujiwara.r/PycharmProjects/mmsegmentation_research/mmsegmentation"
)
DEFAULT_DATASET_ROOT = Path("/mnt/sp4t/fujiwara.r/split_v2")
DEFAULT_CONFIG = (
    DEFAULT_PROJECT_ROOT
    / "configs/unet/hrf_deeplabv3_unet_shared_aspp_dual_bce_dice_v2.py"
)
DEFAULT_CHECKPOINT = (
    DEFAULT_PROJECT_ROOT
    / "work_dirs/hrf_deeplabv3_unet_shared_aspp_dual_bce_dice_v2/epoch_17.pth"
)
DEFAULT_OUTPUT_DIR = (
    DEFAULT_DATASET_ROOT / "shared_aspp_dual_head_pr_curves_validation"
)
DEFAULT_REGIONS = ["aquila", "orionb", "ic5146", "taurus", "polaris"]

IGNORE_LABEL = 255
CONFIRMED_LABEL = 1
CANDIDATE_LABEL = 2


# =============================================================================
# ARGUMENTS
# =============================================================================


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=DEFAULT_PROJECT_ROOT)
    parser.add_argument("--dataset-root", type=Path, default=DEFAULT_DATASET_ROOT)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--device", type=str, default="cuda:0")
    parser.add_argument(
        "--regions",
        nargs="+",
        default=DEFAULT_REGIONS,
        help="例: --regions aquila orionb"
    )
    return parser.parse_args()


# =============================================================================
# VALIDATION PATCH I/O
# =============================================================================


def collect_pairs(dataset_root: Path, region: str):
    inp = dataset_root / "val" / "in" / region
    out = dataset_root / "val" / "out" / region

    if not inp.is_dir():
        raise FileNotFoundError(f"VALIDATION input directory not found: {inp}")
    if not out.is_dir():
        raise FileNotFoundError(f"VALIDATION mask directory not found: {out}")

    pairs = []
    for image_path in sorted(inp.rglob("*_in.png")):
        rel = image_path.relative_to(inp)
        mask_path = out / rel.parent / rel.name.replace("_in.png", "_out.png")
        if not mask_path.is_file():
            raise FileNotFoundError(mask_path)
        pairs.append((image_path, mask_path))

    if not pairs:
        raise RuntimeError(f"No VALIDATION patches found for {region}: {inp}")

    return pairs



def read_mask(path: Path):
    arr = np.asarray(Image.open(path))
    arr = np.squeeze(arr)
    if arr.ndim != 2:
        raise ValueError(f"Mask must be 2D: {path}, shape={arr.shape}")
    return arr.astype(np.uint8)


# =============================================================================
# DUAL-BINARY OUTPUT EXTRACTION
# =============================================================================


def _to_numpy(x):
    if hasattr(x, "data") and not isinstance(x, np.ndarray):
        try:
            x = x.data
        except Exception:
            pass

    if hasattr(x, "detach"):
        x = x.detach().float().cpu().numpy()

    return np.asarray(x)



def logits_to_positive_prob(logits):
    logits = _to_numpy(logits)

    if logits.ndim == 4 and logits.shape[0] == 1:
        logits = logits[0]

    if logits.ndim != 3:
        raise ValueError(f"Unexpected logits shape: {logits.shape}")

    if logits.shape[0] == 2:
        z = logits.astype(np.float64)
        z = z - np.max(z, axis=0, keepdims=True)
        ez = np.exp(z)
        prob = ez[1] / np.sum(ez, axis=0)
        return prob.astype(np.float32)

    if logits.shape[0] == 1:
        z = logits[0].astype(np.float64)
        prob = 1.0 / (1.0 + np.exp(-z))
        return prob.astype(np.float32)

    raise ValueError(f"Binary head must have 1 or 2 channels, got {logits.shape}")



def extract_head_probability(result, head: str):
    """
    Current model output:
        result.class1_prob -> PixelData(data=[1,H,W])
        result.class2_prob -> PixelData(data=[1,H,W])

    head:
        "confirmed" -> Class 1 -> class1_prob
        "candidate" -> Class 2 -> class2_prob
    """
    if head == "confirmed":
        names = ("class1_prob", "confirmed_prob")
    elif head == "candidate":
        names = ("class2_prob", "candidate_prob")
    else:
        raise ValueError(head)

    # SegDataSample-like output
    for name in names:
        if hasattr(result, name):
            obj = getattr(result, name)
            arr = np.squeeze(_to_numpy(obj))
            if arr.ndim == 2:
                return arr.astype(np.float32)
            if arr.ndim == 3 and arr.shape[0] in (1, 2):
                if arr.shape[0] == 1:
                    return arr[0].astype(np.float32)
                return logits_to_positive_prob(arr)

    # dict output
    if isinstance(result, dict):
        for name in names:
            if name in result:
                arr = np.squeeze(_to_numpy(result[name]))
                if arr.ndim == 2:
                    return arr.astype(np.float32)
                if arr.ndim == 3 and arr.shape[0] in (1, 2):
                    if arr.shape[0] == 1:
                        return arr[0].astype(np.float32)
                    return logits_to_positive_prob(arr)

    attrs = [
        x for x in dir(result)
        if ("class1" in x or "class2" in x or "prob" in x or "logit" in x)
        and not x.startswith("_")
    ]
    raise AttributeError(
        f"{head} probability not found in current shared-ASPP dual-head output.\n"
        f"Available related attributes: {attrs}"
    )


# =============================================================================
# PR CALCULATION
# =============================================================================


def pr_curve(y_true, y_score):
    y_true = np.asarray(y_true, dtype=np.uint8)
    y_score = np.asarray(y_score, dtype=np.float64)

    n_pos = int(np.sum(y_true == 1))
    if n_pos == 0:
        raise ValueError("No positive pixels in this evaluation set.")

    order = np.argsort(-y_score, kind="mergesort")
    yt = y_true[order]
    ys = y_score[order]

    tp = np.cumsum(yt == 1, dtype=np.float64)
    fp = np.cumsum(yt == 0, dtype=np.float64)

    distinct = np.where(np.diff(ys) != 0)[0]
    idx = np.r_[distinct, yt.size - 1]

    tp = tp[idx]
    fp = fp[idx]
    thresholds = ys[idx]

    precision_raw = tp / np.maximum(tp + fp, 1.0)
    recall_raw = tp / float(n_pos)

    precision_plot = np.r_[1.0, precision_raw]
    recall_plot = np.r_[0.0, recall_raw]

    return precision_plot, recall_plot, thresholds, precision_raw, recall_raw



def average_precision(precision_raw, recall_raw):
    previous_recall = np.r_[0.0, recall_raw[:-1]]
    return float(np.sum((recall_raw - previous_recall) * precision_raw))



def best_f1_point(thresholds, precision_raw, recall_raw):
    f1 = 2.0 * precision_raw * recall_raw / np.maximum(precision_raw + recall_raw, 1e-12)
    i = int(np.nanargmax(f1))
    return {
        "best_threshold": float(thresholds[i]),
        "precision_at_best_f1": float(precision_raw[i]),
        "recall_at_best_f1": float(recall_raw[i]),
        "best_f1": float(f1[i]),
    }


# =============================================================================
# EVALUATION
# =============================================================================


def evaluate_region(model, inference_model, dataset_root: Path, region: str):
    pairs = collect_pairs(dataset_root, region)

    data = {
        "candidate_true": [],
        "candidate_score": [],
        "confirmed_true": [],
        "confirmed_score": [],
    }

    for i, (image_path, mask_path) in enumerate(pairs, 1):
        gt = read_mask(mask_path)
        result = inference_model(model, str(image_path))

        candidate_score = extract_head_probability(result, "candidate")
        confirmed_score = extract_head_probability(result, "confirmed")

        if candidate_score.shape != gt.shape:
            raise ValueError(
                f"Candidate shape mismatch: {image_path}\n"
                f"prob={candidate_score.shape}, gt={gt.shape}"
            )
        if confirmed_score.shape != gt.shape:
            raise ValueError(
                f"Confirmed shape mismatch: {image_path}\n"
                f"prob={confirmed_score.shape}, gt={gt.shape}"
            )

        valid_c = (gt == 0) | (gt == CANDIDATE_LABEL)
        data["candidate_true"].append((gt[valid_c] == CANDIDATE_LABEL).astype(np.uint8))
        data["candidate_score"].append(candidate_score[valid_c].astype(np.float32))

        valid_f = (gt == 0) | (gt == CONFIRMED_LABEL)
        data["confirmed_true"].append((gt[valid_f] == CONFIRMED_LABEL).astype(np.uint8))
        data["confirmed_score"].append(confirmed_score[valid_f].astype(np.float32))

        if i % 50 == 0 or i == len(pairs):
            print(f"  {region}: {i}/{len(pairs)}")

    result = {"patch_count": len(pairs)}

    for head in ("candidate", "confirmed"):
        y_true = np.concatenate(data[f"{head}_true"])
        y_score = np.concatenate(data[f"{head}_score"])

        precision_plot, recall_plot, thresholds, precision_raw, recall_raw = pr_curve(y_true, y_score)
        ap = average_precision(precision_raw, recall_raw)
        best = best_f1_point(thresholds, precision_raw, recall_raw)

        positive_pixels = int(np.sum(y_true == 1))
        negative_pixels = int(np.sum(y_true == 0))

        result[head] = {
            "precision": precision_plot,
            "recall": recall_plot,
            "auprc": ap,
            "positive_pixels": positive_pixels,
            "negative_pixels": negative_pixels,
            "positive_fraction": positive_pixels / max(positive_pixels + negative_pixels, 1),
            **best,
        }

    return result


# =============================================================================
# PLOTS / CSV
# =============================================================================

TITLE_FONTSIZE = 24
AXIS_LABEL_FONTSIZE = 18
TICK_FONTSIZE = 15
LEGEND_FONTSIZE = 18
LEGEND_TITLE_FONTSIZE = 18

# Display colors for Class 1 / Class 2 comparison figures
CLASS1_COLOR = "#FFD700"  # yellow
CLASS2_COLOR = "#0072BD"  # blue


def plot_regions(results, regions, head, output_path: Path):
    fig, ax = plt.subplots(figsize=(8.5, 7.0))

    for region in regions:
        d = results[region][head]
        ax.plot(
            d["recall"],
            d["precision"],
            linewidth=1.8,
            label=f"{region} (AUPRC={d['auprc']:.3f})",
        )

    display = "Class 2" if head == "candidate" else "Class 1"
    ax.set_xlabel("Recall", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("Precision", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title(f"{display} filament: Precision-Recall curves by region", fontsize=TITLE_FONTSIZE)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.25)
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE)
    ax.legend(loc="best", fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)

    fig.tight_layout()
    fig.savefig(output_path, dpi=250, bbox_inches="tight")
    plt.close(fig)



def plot_same_region(results, region, output_path: Path):
    fig, ax = plt.subplots(figsize=(8.0, 7.0))

    for head, display, color in (
        ("confirmed", "Class 1", CLASS1_COLOR),
        ("candidate", "Class 2", CLASS2_COLOR),
    ):
        d = results[region][head]
        ax.plot(
            d["recall"],
            d["precision"],
            linewidth=2.6,
            color=color,
            label=f"{display} (AUPRC={d['auprc']:.3f})",
        )

    ax.set_xlabel("Recall", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("Precision", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title(f"{region}: Class 1 vs Class 2 PR curve", fontsize=TITLE_FONTSIZE)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.25)
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE)
    ax.legend(loc="best", fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)

    fig.tight_layout()
    fig.savefig(output_path, dpi=250, bbox_inches="tight")
    plt.close(fig)



def save_summary(results, regions, output_dir: Path):
    path = output_dir / "shared_aspp_dual_head_pr_summary_validation.csv"
    fields = [
        "region", "head", "patch_count",
        "positive_pixels", "negative_pixels", "positive_fraction",
        "auprc", "best_threshold", "precision_at_best_f1",
        "recall_at_best_f1", "best_f1",
    ]

    rows = []
    for region in regions:
        for head in ("candidate", "confirmed"):
            d = results[region][head]
            rows.append({
                "region": region,
                "head": head,
                "patch_count": results[region]["patch_count"],
                "positive_pixels": d["positive_pixels"],
                "negative_pixels": d["negative_pixels"],
                "positive_fraction": d["positive_fraction"],
                "auprc": d["auprc"],
                "best_threshold": d["best_threshold"],
                "precision_at_best_f1": d["precision_at_best_f1"],
                "recall_at_best_f1": d["recall_at_best_f1"],
                "best_f1": d["best_f1"],
            })

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    return path


# =============================================================================
# MODEL SANITY CHECK
# =============================================================================

def verify_current_model(model):
    required = (
        "decode_head",
        "class1_classifier",
        "class2_classifier",
    )

    missing = [
        name
        for name in required
        if not hasattr(model, name)
    ]

    if missing:
        raise RuntimeError(
            "This PR script is for the shared-ASPP dual-head model. "
            f"Missing model attributes: {missing}. "
            "Use hrf_deeplabv3_unet_shared_aspp_dual_head_v2.py "
            "with its matching checkpoint."
        )


# =============================================================================
# MAIN
# =============================================================================


def main():
    args = parse_args()

    if str(args.project_root) not in sys.path:
        sys.path.insert(0, str(args.project_root))

    from mmengine.config import Config
    from mmseg.apis import init_model, inference_model

    if not args.config.is_file():
        raise FileNotFoundError(args.config)
    if not args.checkpoint.is_file():
        raise FileNotFoundError(args.checkpoint)

    cfg = Config.fromfile(str(args.config))

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Config     :", args.config)
    print("Checkpoint :", args.checkpoint)
    print("Output     :", args.output_dir)
    print("Regions    :", args.regions)

    model = init_model(cfg, str(args.checkpoint), device=args.device)
    verify_current_model(model)
    model.cfg.test_pipeline = [
        dict(type="LoadImageFromFile"),
        dict(type="PackSegInputs"),
    ]

    results = {}
    for region in args.regions:
        print("\n" + "=" * 80)
        print("Region:", region)
        results[region] = evaluate_region(model, inference_model, args.dataset_root, region)
        for head in ("candidate", "confirmed"):
            d = results[region][head]
            print(
                f"  {head:9s}: AUPRC={d['auprc']:.6f}, "
                f"best_th={d['best_threshold']:.6f}, "
                f"best_F1={d['best_f1']:.6f}"
            )

    candidate_png = args.output_dir / "class2_pr_curve_5regions_validation.png"
    confirmed_png = args.output_dir / "class1_pr_curve_5regions_validation.png"

    plot_regions(results, args.regions, "candidate", candidate_png)
    plot_regions(results, args.regions, "confirmed", confirmed_png)

    for region in args.regions:
        plot_same_region(
            results,
            region,
            args.output_dir / f"{region}_class1_vs_class2_pr.png",
        )

    csv_path = save_summary(results, args.regions, args.output_dir)

    print("\nSaved:", candidate_png)
    print("Saved:", confirmed_png)
    print("Saved:", csv_path)
    for region in args.regions:
        print("Saved:", args.output_dir / f"{region}_class1_vs_class2_pr.png")


if __name__ == "__main__":
    main()
