#!/usr/bin/env python3
from __future__ import annotations

"""
Read the saved outputs from make_all_regions_lshape_preferred_no_leakage_v4.py
and visualize the train / val / test split boundaries on the original FITS image.

What this script uses
---------------------
For each region, it reads:
  1) the original column-density FITS
  2) metadata/<region>/patch_manifest.csv

Then it saves a PNG like:
  split_boundary_visualization/<region>_split_boundaries.png

Display style
-------------
- train : blue
- val   : green
- test  : red

Each split is shown with:
- semi-transparent patch fill
- patch boundary lines
- a large "TRAIN / VAL / TEST" label near the split center

This makes the split much easier to understand than the old thin-line overview.
"""

import csv
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------

DATASET_ROOT = Path("/mnt/sp4t/fujiwara.r/split_v2")

RUN_ALL_CONFIGURED_REGIONS = False
SINGLE_REGION_NAME = "taurus"

REGION_DATASETS = {
    "aquila": {
        "enabled": True,
        "column_density_fits": (
            "/mnt/sp4t/fujiwara.r/HGBS_aquilaM2_hires_column_density_map.fits"
        ),
    },
    "orionb": {
        "enabled": True,
        "column_density_fits": (
            "/mnt/sp4t/fujiwara.r/HGBS_orionB_hires_column_density_map.fits"
        ),
    },
    "ic5146": {
        "enabled": True,
        "column_density_fits": (
            "/mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map.fits"
        ),
    },
    "taurus": {
        "enabled": True,
        "column_density_fits": (
            "/mnt/sp4t/fujiwara.r/HGBS_tauN3_hires_column_density_map.fits"
        ),
    },
    "polaris": {
        "enabled": True,
        "column_density_fits": (
            "/mnt/sp4t/fujiwara.r/HGBS_polaris_hires_column_density_map.fits"
        ),
    },
}

PATCH_SIZE = 256

OUTPUT_DIR = DATASET_ROOT / "split_boundary_visualization"

# split colors
SPLIT_COLORS = {
    # Same color convention as make_all_regions_lshape_preferred_no_leakage_v4.py
    # train = tab:blue, val = tab:green, test = tab:red
    "train": {
        "edge": "#1f77b4",   # matplotlib tab:blue
        "fill": "#1f77b4",
    },
    "val": {
        "edge": "#2ca02c",   # matplotlib tab:green
        "fill": "#2ca02c",
    },
    "test": {
        "edge": "#d62728",   # matplotlib tab:red
        "fill": "#d62728",
    },
}

FILL_ALPHA = 0.18
EDGE_ALPHA = 0.85
EDGE_LINEWIDTH = 1.0

TITLE_FONTSIZE = 18
LABEL_FONTSIZE = 16
LABEL_BBOX_ALPHA = 0.75

FIG_DPI = 200
FIGSIZE = (10, 10)
CMAP = "inferno"

# if True, also saves a patch-color-only image on white background
SAVE_SEPARATE_PLAIN_MAP = True

# -----------------------------------------------------------------------------
# IMPORTS THAT MAY NOT EXIST IN SOME ENVIRONMENTS
# -----------------------------------------------------------------------------

try:
    from astropy.io import fits
except ImportError as exc:
    raise ImportError(
        "astropy is required to read FITS files.\n"
        "Install it with: pip install astropy"
    ) from exc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle, Patch


# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------

def normalize_region_name(region_name: str) -> str:
    region = region_name.strip().lower().replace(" ", "_")
    if region not in REGION_DATASETS:
        raise ValueError(f"Unknown region: {region_name}")
    return region


def region_list() -> List[str]:
    if RUN_ALL_CONFIGURED_REGIONS:
        return [
            region
            for region, info in REGION_DATASETS.items()
            if bool(info.get("enabled", True))
        ]
    return [normalize_region_name(SINGLE_REGION_NAME)]


def load_fits_2d(path: Path) -> np.ndarray:
    with fits.open(path, memmap=True) as hdul:
        for hdu in hdul:
            if hdu.data is None:
                continue
            arr = np.squeeze(np.asarray(hdu.data))
            if arr.ndim == 2:
                return arr.astype(np.float32)
    raise ValueError(f"No 2D image found in FITS: {path}")


def load_patch_manifest(csv_path: Path) -> List[Dict[str, object]]:
    if not csv_path.is_file():
        raise FileNotFoundError(f"patch_manifest.csv not found: {csv_path}")

    records: List[Dict[str, object]] = []

    with csv_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"split", "y", "x"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(
                f"patch_manifest.csv is missing required columns. "
                f"Found={reader.fieldnames}"
            )

        for row in reader:
            records.append(
                {
                    "split": str(row["split"]).strip(),
                    "y": int(row["y"]),
                    "x": int(row["x"]),
                    "stem": str(row.get("stem", "")),
                }
            )

    if not records:
        raise ValueError(f"No patch records found in: {csv_path}")

    return records


def make_display_image(col_density: np.ndarray) -> Tuple[np.ndarray, float, float]:
    positive = np.isfinite(col_density) & (col_density > 0)

    vals = col_density[positive]
    if vals.size == 0:
        raise ValueError("FITS image has no positive finite pixels.")

    vmin = np.percentile(vals, 1.0)
    vmax = np.percentile(vals, 99.5)
    if vmax <= vmin:
        vmax = float(vals.max())

    return col_density, float(vmin), float(vmax)


def build_split_map(
    image_shape: Tuple[int, int],
    records: List[Dict[str, object]],
) -> Dict[str, np.ndarray]:
    h, w = image_shape

    out: Dict[str, np.ndarray] = {
        "train": np.zeros((h, w), dtype=bool),
        "val": np.zeros((h, w), dtype=bool),
        "test": np.zeros((h, w), dtype=bool),
    }

    for record in records:
        split = str(record["split"])
        y = int(record["y"])
        x = int(record["x"])
        y1 = y + PATCH_SIZE
        x1 = x + PATCH_SIZE

        if split not in out:
            raise ValueError(f"Unknown split name in manifest: {split}")

        if y < 0 or x < 0 or y1 > h or x1 > w:
            raise ValueError(
                "Patch from manifest is outside the FITS image.\n"
                f"split={split}, y={y}, x={x}, patch={PATCH_SIZE}, "
                f"image_shape={(h, w)}"
            )

        out[split][y:y1, x:x1] = True

    return out


def build_rgba(mask: np.ndarray, hex_color: str, alpha: float) -> np.ndarray:
    rgb = np.array(
        [int(hex_color[1:3], 16), int(hex_color[3:5], 16), int(hex_color[5:7], 16)],
        dtype=np.float32,
    ) / 255.0

    rgba = np.zeros((mask.shape[0], mask.shape[1], 4), dtype=np.float32)
    rgba[mask, :3] = rgb
    rgba[mask, 3] = alpha
    return rgba


def compute_centroid(mask: np.ndarray) -> Tuple[float, float]:
    ys, xs = np.where(mask)
    if len(xs) == 0:
        return (0.0, 0.0)
    return (float(xs.mean()), float(ys.mean()))


def draw_patch_edges(ax, records: List[Dict[str, object]]) -> None:
    for record in records:
        split = str(record["split"])
        y = int(record["y"])
        x = int(record["x"])

        ax.add_patch(
            Rectangle(
                (x, y),
                PATCH_SIZE,
                PATCH_SIZE,
                fill=False,
                linewidth=EDGE_LINEWIDTH,
                edgecolor=SPLIT_COLORS[split]["edge"],
                alpha=EDGE_ALPHA,
            )
        )


def draw_split_labels(ax, split_masks: Dict[str, np.ndarray]) -> None:
    names = {
        "train": "TRAIN",
        "val": "VAL",
        "test": "TEST",
    }

    for split_name in ("train", "val", "test"):
        mask = split_masks[split_name]
        if not np.any(mask):
            continue

        cx, cy = compute_centroid(mask)

        ax.text(
            cx,
            cy,
            names[split_name],
            color="white",
            fontsize=LABEL_FONTSIZE,
            fontweight="bold",
            ha="center",
            va="center",
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor=SPLIT_COLORS[split_name]["edge"],
                edgecolor="none",
                alpha=LABEL_BBOX_ALPHA,
            ),
        )


def save_visualization(
    region: str,
    col_density: np.ndarray,
    records: List[Dict[str, object]],
    output_png: Path,
    plain_png: Path | None = None,
) -> None:
    _, vmin, vmax = make_display_image(col_density)
    split_masks = build_split_map(col_density.shape, records)

    fig, ax = plt.subplots(figsize=FIGSIZE)

    im = ax.imshow(
        col_density,
        origin="lower",
        cmap=CMAP,
        norm=LogNorm(vmin=max(vmin, 1e-30), vmax=vmax),
        interpolation="nearest",
    )

    # semi-transparent split fill
    for split_name in ("train", "val", "test"):
        rgba = build_rgba(
            split_masks[split_name],
            SPLIT_COLORS[split_name]["fill"],
            FILL_ALPHA,
        )
        ax.imshow(rgba, origin="lower", interpolation="nearest")

    # patch edges
    draw_patch_edges(ax, records)

    # large labels
    draw_split_labels(ax, split_masks)

    ax.set_title(
        f"{region}  split visualization\n"
        "train=blue, val=green, test=red",
        fontsize=TITLE_FONTSIZE,
    )
    ax.set_xlabel("x [pixels]")
    ax.set_ylabel("y [pixels]")

    legend_handles = [
        Patch(facecolor=SPLIT_COLORS["train"]["fill"], edgecolor=SPLIT_COLORS["train"]["edge"], alpha=FILL_ALPHA, label="train"),
        Patch(facecolor=SPLIT_COLORS["val"]["fill"], edgecolor=SPLIT_COLORS["val"]["edge"], alpha=FILL_ALPHA, label="val"),
        Patch(facecolor=SPLIT_COLORS["test"]["fill"], edgecolor=SPLIT_COLORS["test"]["edge"], alpha=FILL_ALPHA, label="test"),
    ]
    ax.legend(handles=legend_handles, loc="upper right")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Column Density")

    fig.tight_layout()
    fig.savefig(output_png, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

    if plain_png is not None:
        fig, ax = plt.subplots(figsize=FIGSIZE)
        ax.set_facecolor("white")

        for split_name in ("train", "val", "test"):
            rgba = build_rgba(
                split_masks[split_name],
                SPLIT_COLORS[split_name]["fill"],
                0.35,
            )
            ax.imshow(rgba, origin="lower", interpolation="nearest")

        draw_patch_edges(ax, records)
        draw_split_labels(ax, split_masks)

        ax.set_title(
            f"{region}  split map only\n"
            "train=blue, val=green, test=red",
            fontsize=TITLE_FONTSIZE,
        )
        ax.set_xlabel("x [pixels]")
        ax.set_ylabel("y [pixels]")
        ax.legend(handles=legend_handles, loc="upper right")
        fig.tight_layout()
        fig.savefig(plain_png, dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)


def process_region(region: str) -> None:
    region = normalize_region_name(region)

    fits_path = Path(REGION_DATASETS[region]["column_density_fits"])
    manifest_path = DATASET_ROOT / "metadata" / region / "patch_manifest.csv"

    if not fits_path.is_file():
        raise FileNotFoundError(f"FITS not found: {fits_path}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 90)
    print(f"Processing region: {region}")
    print(f"FITS    : {fits_path}")
    print(f"Manifest: {manifest_path}")

    col_density = load_fits_2d(fits_path)
    records = load_patch_manifest(manifest_path)

    output_png = OUTPUT_DIR / f"{region}_split_boundaries.png"
    plain_png = OUTPUT_DIR / f"{region}_split_boundaries_plain.png" if SAVE_SEPARATE_PLAIN_MAP else None

    save_visualization(
        region=region,
        col_density=col_density,
        records=records,
        output_png=output_png,
        plain_png=plain_png,
    )

    split_counts = {"train": 0, "val": 0, "test": 0}
    for r in records:
        split_counts[str(r["split"])] += 1

    print("Saved:")
    print(f"  {output_png}")
    if plain_png is not None:
        print(f"  {plain_png}")
    print("Patch counts:")
    print(
        f"  train={split_counts['train']}, "
        f"val={split_counts['val']}, "
        f"test={split_counts['test']}"
    )


def main() -> None:
    regions = region_list()

    print("=" * 90)
    print("Saved split-boundary visualization")
    print(f"Dataset root: {DATASET_ROOT}")
    print(f"Output dir  : {OUTPUT_DIR}")
    print("Regions     : " + ", ".join(regions))
    print("=" * 90)

    for region in regions:
        process_region(region)

    print("=" * 90)
    print("Finished.")
    print(f"Output dir: {OUTPUT_DIR}")
    print("=" * 90)


if __name__ == "__main__":
    main()
