#!/usr/bin/env python3
from __future__ import annotations

"""
Visualize the saved train / val / test patch split as PNG files and highlight
disconnected islands / isolated patches.

This script reads:
    <DATASET_ROOT>/metadata/<region>/patch_manifest.csv

and saves:
    <DATASET_ROOT>/metadata/spatial_island_check_png/<region>_split_islands.png

What is shown in the PNG:
- train patches : blue
- val patches   : green
- test patches  : red
- each patch boundary is drawn
- disconnected components are identified by 4-neighbour connectivity
- isolated single-patch components are marked strongly
- small islands are annotated

This script does NOT rewrite the dataset.
It only reads the saved manifest and writes visualization PNGs + summary CSV/JSON.
"""

import csv
import json
from collections import deque
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch


# =============================================================================
# USER SETTINGS
# =============================================================================

DATASET_ROOT = Path("/mnt/sp4t/fujiwara.r/split_v2")

REGIONS = (
    "aquila",
    "orionb",
    "ic5146",
    "taurus",
    "polaris",
)

PATCH_SIZE = 256

# Components with <= this many patches (except the largest component) are
# marked as "small islands".
SMALL_ISLAND_MAX_PATCHES = 3

OUTPUT_DIR = DATASET_ROOT / "metadata" / "spatial_island_check_png"

# same convention as your current split-visualization
SPLIT_COLORS = {
    "train": "#1f77b4",  # blue
    "val":   "#2ca02c",  # green
    "test":  "#d62728",  # red
}

PATCH_FILL_ALPHA = 0.22
PATCH_EDGE_ALPHA = 0.95
PATCH_EDGE_LINEWIDTH = 1.0

# highlight for suspicious components
ISLAND_EDGE_COLOR = "black"
ISLAND_EDGE_LINEWIDTH = 2.0
ISOLATED_FILL_COLOR = "#FFD700"  # gold
ISOLATED_FILL_ALPHA = 0.55

FIG_DPI = 200
DEFAULT_FIGSIZE = (8, 8)

SAVE_SUMMARY_JSON = True
SAVE_SUMMARY_CSV = True


# =============================================================================
# IO
# =============================================================================

def read_manifest(path: Path) -> List[Dict[str, object]]:
    if not path.is_file():
        raise FileNotFoundError(f"Manifest not found: {path}")

    rows: List[Dict[str, object]] = []

    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"region", "split", "y", "x", "stem"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")

        for row in reader:
            rows.append(
                {
                    "region": str(row["region"]).strip(),
                    "split": str(row["split"]).strip(),
                    "stem": str(row["stem"]).strip(),
                    "y": int(row["y"]),
                    "x": int(row["x"]),
                }
            )

    if not rows:
        raise ValueError(f"No patch rows found in: {path}")

    return rows


# =============================================================================
# CONNECTIVITY
# =============================================================================

def neighbors4(coord: Tuple[int, int]) -> Tuple[Tuple[int, int], ...]:
    y, x = coord
    s = PATCH_SIZE
    return (
        (y - s, x),
        (y + s, x),
        (y, x - s),
        (y, x + s),
    )


def connected_components4(coords: Set[Tuple[int, int]]) -> List[Set[Tuple[int, int]]]:
    remaining = set(coords)
    components: List[Set[Tuple[int, int]]] = []

    while remaining:
        seed = next(iter(remaining))
        remaining.remove(seed)

        component = {seed}
        q = deque([seed])

        while q:
            current = q.popleft()
            for nb in neighbors4(current):
                if nb in remaining:
                    remaining.remove(nb)
                    component.add(nb)
                    q.append(nb)

        components.append(component)

    components.sort(
        key=lambda comp: (
            -len(comp),
            min(y for y, _ in comp),
            min(x for _, x in comp),
        )
    )
    return components


# =============================================================================
# ANALYSIS
# =============================================================================

def analyze_region(region: str, rows: Sequence[Dict[str, object]]) -> Dict[str, object]:
    splits_out: Dict[str, object] = {}
    region_review = False

    for split in ("train", "val", "test"):
        coords = {
            (int(row["y"]), int(row["x"]))
            for row in rows
            if str(row["split"]) == split
        }

        comps = connected_components4(coords)

        comp_infos = []
        for i, comp in enumerate(comps, start=1):
            sorted_coords = sorted(comp)
            isolated = (len(comp) == 1)
            small_island = (i > 1 and len(comp) <= SMALL_ISLAND_MAX_PATCHES)

            comp_infos.append(
                {
                    "component_id": i,
                    "patch_count": len(comp),
                    "isolated": isolated,
                    "small_island": small_island,
                    "coordinates": [{"y": y, "x": x} for y, x in sorted_coords],
                }
            )

        if len(comps) > 1:
            region_review = True

        splits_out[split] = {
            "patch_count": len(coords),
            "connected_component_count": len(comps),
            "isolated_patch_count": sum(len(c) == 1 for c in comps),
            "components": comp_infos,
        }

    return {
        "region": region,
        "status": "REVIEW" if region_review else "OK",
        "splits": splits_out,
    }


# =============================================================================
# VISUALIZATION
# =============================================================================

def calc_bounds(rows: Sequence[Dict[str, object]]) -> Tuple[int, int, int, int]:
    xs = [int(r["x"]) for r in rows]
    ys = [int(r["y"]) for r in rows]
    min_x = min(xs)
    min_y = min(ys)
    max_x = max(xs) + PATCH_SIZE
    max_y = max(ys) + PATCH_SIZE
    return min_x, min_y, max_x, max_y


def choose_figsize(width_px: int, height_px: int) -> Tuple[float, float]:
    if width_px <= 0 or height_px <= 0:
        return DEFAULT_FIGSIZE
    scale = max(width_px, height_px)
    base = 8.0
    if scale > 8000:
        base = 10.0
    if scale > 11000:
        base = 12.0
    aspect = width_px / max(height_px, 1)
    if aspect >= 1.0:
        return (base * aspect, base)
    return (base, base / max(aspect, 1e-6))


def component_centroid(coords: List[Dict[str, int]]) -> Tuple[float, float]:
    xs = [c["x"] + PATCH_SIZE / 2 for c in coords]
    ys = [c["y"] + PATCH_SIZE / 2 for c in coords]
    return sum(xs) / len(xs), sum(ys) / len(ys)


def draw_patch(ax, x: int, y: int, split: str) -> None:
    ax.add_patch(
        Rectangle(
            (x, y),
            PATCH_SIZE,
            PATCH_SIZE,
            facecolor=SPLIT_COLORS[split],
            edgecolor=SPLIT_COLORS[split],
            linewidth=PATCH_EDGE_LINEWIDTH,
            alpha=PATCH_FILL_ALPHA,
        )
    )
    ax.add_patch(
        Rectangle(
            (x, y),
            PATCH_SIZE,
            PATCH_SIZE,
            fill=False,
            edgecolor=SPLIT_COLORS[split],
            linewidth=PATCH_EDGE_LINEWIDTH,
            alpha=PATCH_EDGE_ALPHA,
        )
    )


def draw_island_highlight(
    ax,
    x: int,
    y: int,
    isolated: bool,
) -> None:
    if isolated:
        ax.add_patch(
            Rectangle(
                (x, y),
                PATCH_SIZE,
                PATCH_SIZE,
                facecolor=ISOLATED_FILL_COLOR,
                edgecolor=ISLAND_EDGE_COLOR,
                linewidth=ISLAND_EDGE_LINEWIDTH,
                alpha=ISOLATED_FILL_ALPHA,
            )
        )
    else:
        ax.add_patch(
            Rectangle(
                (x, y),
                PATCH_SIZE,
                PATCH_SIZE,
                fill=False,
                edgecolor=ISLAND_EDGE_COLOR,
                linewidth=ISLAND_EDGE_LINEWIDTH,
                alpha=1.0,
            )
        )


def save_png(
    region: str,
    rows: Sequence[Dict[str, object]],
    analysis: Dict[str, object],
    out_png: Path,
) -> None:
    min_x, min_y, max_x, max_y = calc_bounds(rows)
    width_px = max_x - min_x
    height_px = max_y - min_y

    fig, ax = plt.subplots(figsize=choose_figsize(width_px, height_px))
    ax.set_facecolor("white")

    # draw all patches
    for row in rows:
        split = str(row["split"])
        draw_patch(ax, int(row["x"]), int(row["y"]), split)

    # draw disconnected component highlights
    splits = analysis["splits"]
    for split in ("train", "val", "test"):
        components = splits[split]["components"]

        if len(components) <= 1:
            continue

        # largest component = component_id 1, normally not highlighted as suspicious
        for comp in components[1:]:
            isolated = bool(comp["isolated"])
            small_island = bool(comp["small_island"])

            for coord in comp["coordinates"]:
                draw_island_highlight(
                    ax,
                    int(coord["x"]),
                    int(coord["y"]),
                    isolated=isolated,
                )

            cx, cy = component_centroid(comp["coordinates"])
            label = f"{split}\nC{comp['component_id']}\n{comp['patch_count']} patch"
            if comp["patch_count"] > 1:
                label += "es"
            if isolated:
                label += "\nISOLATED"
            elif small_island:
                label += "\nSMALL ISLAND"
            else:
                label += "\nDISCONNECTED"

            ax.text(
                cx,
                cy,
                label,
                ha="center",
                va="center",
                fontsize=9,
                color="black",
                bbox=dict(
                    boxstyle="round,pad=0.25",
                    facecolor="white",
                    edgecolor="black",
                    alpha=0.85,
                ),
            )

    # main split labels on largest component
    for split in ("train", "val", "test"):
        components = splits[split]["components"]
        if not components:
            continue
        largest = components[0]
        cx, cy = component_centroid(largest["coordinates"])

        ax.text(
            cx,
            cy,
            split.upper(),
            ha="center",
            va="center",
            fontsize=16,
            fontweight="bold",
            color="white",
            bbox=dict(
                boxstyle="round,pad=0.35",
                facecolor=SPLIT_COLORS[split],
                edgecolor="none",
                alpha=0.75,
            ),
        )

    ax.set_title(
        f"{region} split map\n"
        f"train=blue, val=green, test=red"
        + ("  (REVIEW: disconnected split exists)" if analysis["status"] == "REVIEW" else "  (OK)"),
        fontsize=16,
    )
    ax.set_xlabel("x [pixels]")
    ax.set_ylabel("y [pixels]")

    # add margin
    margin = PATCH_SIZE
    ax.set_xlim(min_x - margin, max_x + margin)
    ax.set_ylim(min_y - margin, max_y + margin)
    ax.set_aspect("equal")

    legend_items = [
        Patch(facecolor=SPLIT_COLORS["train"], edgecolor=SPLIT_COLORS["train"], alpha=PATCH_FILL_ALPHA, label="train"),
        Patch(facecolor=SPLIT_COLORS["val"], edgecolor=SPLIT_COLORS["val"], alpha=PATCH_FILL_ALPHA, label="val"),
        Patch(facecolor=SPLIT_COLORS["test"], edgecolor=SPLIT_COLORS["test"], alpha=PATCH_FILL_ALPHA, label="test"),
        Patch(facecolor=ISOLATED_FILL_COLOR, edgecolor=ISLAND_EDGE_COLOR, alpha=ISOLATED_FILL_ALPHA, label="isolated patch"),
        Patch(facecolor="white", edgecolor=ISLAND_EDGE_COLOR, alpha=1.0, label="small / disconnected island"),
    ]
    ax.legend(handles=legend_items, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# SUMMARY FILES
# =============================================================================

def save_summary_json(results: Sequence[Dict[str, object]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)


def save_summary_csv(results: Sequence[Dict[str, object]], path: Path) -> None:
    fields = [
        "region",
        "status",
        "split",
        "patch_count",
        "connected_component_count",
        "isolated_patch_count",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()

        for res in results:
            for split in ("train", "val", "test"):
                info = res["splits"][split]
                writer.writerow(
                    {
                        "region": res["region"],
                        "status": res["status"],
                        "split": split,
                        "patch_count": info["patch_count"],
                        "connected_component_count": info["connected_component_count"],
                        "isolated_patch_count": info["isolated_patch_count"],
                    }
                )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []

    print("=" * 88)
    print("Visualize saved split map + disconnected islands")
    print(f"Dataset root: {DATASET_ROOT}")
    print(f"Output dir  : {OUTPUT_DIR}")
    print("=" * 88)

    for region in REGIONS:
        manifest = DATASET_ROOT / "metadata" / region / "patch_manifest.csv"
        rows = read_manifest(manifest)
        analysis = analyze_region(region, rows)
        all_results.append(analysis)

        out_png = OUTPUT_DIR / f"{region}_split_islands.png"
        save_png(region, rows, analysis, out_png)

        print(f"[{region}] status={analysis['status']}  saved -> {out_png}")

    if SAVE_SUMMARY_JSON:
        save_summary_json(all_results, OUTPUT_DIR / "spatial_island_summary.json")
    if SAVE_SUMMARY_CSV:
        save_summary_csv(all_results, OUTPUT_DIR / "spatial_island_summary.csv")

    print("=" * 88)
    print("Finished.")
    print(f"PNG files saved in: {OUTPUT_DIR}")
    print("=" * 88)


if __name__ == "__main__":
    main()
