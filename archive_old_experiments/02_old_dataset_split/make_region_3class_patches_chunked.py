#!/usr/bin/env python3
"""Build a chunked multi-region 3-class dataset for MMSegmentation.

Expected label definition in LABEL_MAP_PATH (.npy, .fits, .fit, or .fts):
    0: background
    1: FilFinder / DisPerSE normal-threshold union
    2: relaxed-threshold detection minus class 1
  255: ignore

Processing order:
1. Read one region's column-density FITS and 3-class label map.
2. Search a spatial train/val/test diagonal split (default 70/15/15).
3. Merge classes 1 and 2 for 8-neighbour connected-component analysis and
   ignore whole filament components that cross split boundaries.
4. Convert log10(column density) to an 8-bit, three-channel PNG image. The
   normalization limits are estimated from the current region's train split.
5. Save 256x256 image/mask patch pairs under a region directory, grouped into
   subdirectories containing at most FILES_PER_SUBDIR pairs each.

Output layout example (region=aquila, files_per_subdir=1000):
    DATASET_ROOT/
      train/in/aquila/0000/*_in.png
      train/out/aquila/0000/*_out.png
      train/in/aquila/0001/*_in.png
      train/out/aquila/0001/*_out.png
      val/in/aquila/0000/*_in.png
      val/out/aquila/0000/*_out.png
      test/in/aquila/0000/*_in.png
      test/out/aquila/0000/*_out.png
      metadata/aquila/*.npy, *.json, *.csv, *.png

Run the script once per region. Cleaning removes only that region's existing
output, so Aquila, Orion B, IC 5146, Taurus and Polaris can coexist.

Mask PNG files are single-channel uint8 images preserving 0, 1, 2 and 255.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from dataclasses import asdict, dataclass, replace
from pathlib import Path
import re
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from PIL import Image
from scipy.ndimage import binary_dilation, label as ndi_label

# =============================================================================
# Default example paths
# =============================================================================
REGION_NAME = "ic5146"
COLUMN_DENSITY_FITS = (
    "/mnt/sp4t/fujiwara.r/"
    "HGBS_ic5146_hires_column_density_map.fits"
)
LABEL_MAP_PATH = "/mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map/vis_3class_map/HGBS_ic5146_hires_column_density_map_label3class.fits"
DATASET_ROOT = "/mnt/sp4t/fujiwara.r/split"

# =============================================================================
# Labels and spatial split
# =============================================================================
BACKGROUND_LABEL = 0
CLASS1_LABEL = 1
CLASS2_LABEL = 2
IGNORE_LABEL = 255
ALLOWED_LABELS = {BACKGROUND_LABEL, CLASS1_LABEL, CLASS2_LABEL, IGNORE_LABEL}

# Invalid FITS borders/no-data are ignored. Dilation uses an 8-neighbourhood.
MARGIN_PX = 32

TARGET_SPLIT_RATIOS: Mapping[str, float] = {
    "train": 0.70,
    "val": 0.15,
    "test": 0.15,
}
MIN_SPLIT_RATIO = 0.10
RATIO_SEARCH_STEPS = 21
RATIO_SCORE_WEIGHT = 5.0
THETA_RANGE_DEG = (15.0, 165.0)
THETA_STEP_DEG = 2.0
BALANCE_TARGET_CLASSES = (CLASS1_LABEL, CLASS2_LABEL)

# Ignore a narrow zone around a connected filament that crosses split borders.
# Set to 0 to ignore only the filament pixels themselves.
CUT_COMPONENT_DILATION_PX = 2

# =============================================================================
# Patch generation
# =============================================================================
PATCH_SIZE = 256
TRAIN_STRIDE = 128
EVAL_STRIDE = 256

# A patch must lie completely inside one split geometrically. Within it, at
# least this fraction must remain valid after no-data/cut-component filtering.
MIN_PATCH_VALID_RATIO = 0.95

# Preserve all train patches by default. Set a float such as 1.0 to keep at
# most one background-only patch for every filament-containing train patch.
MAX_TRAIN_BACKGROUND_TO_POSITIVE: Optional[float] = None
MIN_FILAMENT_PIXELS_FOR_POSITIVE = 8

# log10 column-density display range, calculated on the train split only.
NORMALIZE_LOW_PERCENTILE = 1.0
NORMALIZE_HIGH_PERCENTILE = 99.0

RANDOM_SEED = 0
FILES_PER_SUBDIR = 1000
CLEAN_REGION_OUTPUT = True


@dataclass(frozen=True)
class PatchRecord:
    split: str
    region: str
    stem: str
    y: int
    x: int
    valid_ratio: float
    background_pixels: int
    class1_pixels: int
    class2_pixels: int
    ignore_pixels: int
    storage_subdir: str = ""

    @property
    def filament_pixels(self) -> int:
        return self.class1_pixels + self.class2_pixels


# =============================================================================
# I/O and validation
# =============================================================================
def _read_fits_array(path: str, hdu: int = 0) -> np.ndarray:
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise ImportError(
            "astropy is required to read FITS files. Install it in the "
            "MMSegmentation environment, for example: pip install astropy"
        ) from exc

    with fits.open(path, memmap=True) as hdul:
        if hdu >= len(hdul):
            raise IndexError(
                f"FITS HDU index {hdu} does not exist in {path}; "
                f"available HDUs: 0..{len(hdul) - 1}"
            )
        if hdul[hdu].data is None:
            raise ValueError(f"FITS HDU {hdu} contains no image data: {path}")
        return np.array(hdul[hdu].data, copy=True)


def _coerce_label_array(raw: np.ndarray, source: str) -> np.ndarray:
    """Convert a label array safely to uint8 without wrapping bad values.

    Non-finite floating-point pixels are treated as IGNORE_LABEL. All finite
    values must be integer-valued and one of 0, 1, 2 or 255.
    """
    array = np.squeeze(np.asarray(raw))
    if array.ndim != 2:
        raise ValueError(f"Label map must be 2-D after squeeze: {source}, shape={array.shape}")

    if np.issubdtype(array.dtype, np.floating):
        finite = np.isfinite(array)
        finite_values = array[finite]
        rounded = np.rint(finite_values)
        if finite_values.size and not np.allclose(finite_values, rounded, atol=1e-6):
            bad = finite_values[np.abs(finite_values - rounded) > 1e-6][:10]
            raise ValueError(
                f"Label FITS/array contains non-integer finite values in {source}: "
                f"examples={bad.tolist()}"
            )
        label_map = np.full(array.shape, IGNORE_LABEL, dtype=np.uint8)
        if rounded.size:
            numeric = rounded.astype(np.int64)
            unexpected = set(int(v) for v in np.unique(numeric)) - ALLOWED_LABELS
            if unexpected:
                raise ValueError(
                    f"Unexpected labels {sorted(unexpected)} in {source}. "
                    f"Allowed labels are {sorted(ALLOWED_LABELS)}."
                )
            label_map[finite] = numeric.astype(np.uint8)
        return label_map

    numeric = array.astype(np.int64, copy=False)
    unexpected = set(int(v) for v in np.unique(numeric)) - ALLOWED_LABELS
    if unexpected:
        raise ValueError(
            f"Unexpected labels {sorted(unexpected)} in {source}. "
            f"Allowed labels are {sorted(ALLOWED_LABELS)}."
        )
    return numeric.astype(np.uint8)


def load_data(
    column_density_fits: str,
    label_map_path: str,
    column_density_hdu: int = 0,
    label_hdu: int = 0,
) -> Tuple[np.ndarray, np.ndarray]:
    col_density = np.squeeze(
        np.asarray(
            _read_fits_array(column_density_fits, hdu=column_density_hdu),
            dtype=np.float64,
        )
    )

    label_path = Path(label_map_path)
    suffix = label_path.suffix.lower()
    if suffix == ".npy":
        raw_label = np.load(label_path)
    elif suffix in {".fits", ".fit", ".fts"}:
        raw_label = _read_fits_array(str(label_path), hdu=label_hdu)
    else:
        raise ValueError(
            f"Unsupported label file type: {label_path}. "
            "Use .npy, .fits, .fit or .fts."
        )
    label_map = _coerce_label_array(raw_label, str(label_path))

    if col_density.ndim != 2:
        raise ValueError(
            "Column-density input must be a 2-D array after squeeze: "
            f"shape={col_density.shape}"
        )
    if col_density.shape != label_map.shape:
        raise ValueError(
            f"Shape mismatch: column_density={col_density.shape}, "
            f"label_map={label_map.shape}"
        )

    return col_density, label_map



def build_three_class_label_map(
    filfinder_normal: np.ndarray,
    disperse_normal: np.ndarray,
    filfinder_relaxed: np.ndarray,
    disperse_relaxed: np.ndarray,
    invalid_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Build the requested mutually exclusive 3-class label map.

    class 1 = FilFinder(normal) OR DisPerSE(normal)
    class 2 = [FilFinder(relaxed) OR DisPerSE(relaxed)] minus class 1
    class 0 = everything else
    """
    arrays = [
        np.asarray(filfinder_normal, dtype=bool),
        np.asarray(disperse_normal, dtype=bool),
        np.asarray(filfinder_relaxed, dtype=bool),
        np.asarray(disperse_relaxed, dtype=bool),
    ]
    if len({array.shape for array in arrays}) != 1:
        raise ValueError("All FilFinder/DisPerSE masks must have the same shape.")

    class1 = arrays[0] | arrays[1]
    relaxed_union = arrays[2] | arrays[3]
    class2 = relaxed_union & ~class1

    label_map = np.zeros(class1.shape, dtype=np.uint8)
    label_map[class2] = CLASS2_LABEL
    label_map[class1] = CLASS1_LABEL
    if invalid_mask is not None:
        invalid = np.asarray(invalid_mask, dtype=bool)
        if invalid.shape != label_map.shape:
            raise ValueError("invalid_mask shape does not match the detection masks.")
        label_map[invalid] = IGNORE_LABEL
    return label_map

def normalize_region_name(region_name: str) -> str:
    region = region_name.strip().lower().replace(" ", "_")
    if not region or not re.fullmatch(r"[a-z0-9_-]+", region):
        raise ValueError(
            "region_name may contain only letters, digits, '_' and '-': "
            f"{region_name!r}"
        )
    return region


def prepare_output_dirs(root: Path, region: str, clean_region: bool) -> Path:
    """Prepare only one region; never remove other regions' patches."""
    root.mkdir(parents=True, exist_ok=True)
    for split in ("train", "val", "test"):
        for modality in ("in", "out"):
            region_dir = root / split / modality / region
            if clean_region and region_dir.exists():
                shutil.rmtree(region_dir)
            region_dir.mkdir(parents=True, exist_ok=True)

    metadata_dir = root / "metadata" / region
    if clean_region and metadata_dir.exists():
        shutil.rmtree(metadata_dir)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    return metadata_dir


def build_invalid_mask(col_density: np.ndarray, margin_px: int) -> np.ndarray:
    h, w = col_density.shape
    invalid = np.zeros((h, w), dtype=bool)

    if margin_px > 0:
        invalid[:margin_px, :] = True
        invalid[-margin_px:, :] = True
        invalid[:, :margin_px] = True
        invalid[:, -margin_px:] = True

    no_data = ~np.isfinite(col_density) | (col_density <= 0)
    if margin_px > 0 and no_data.any():
        no_data = binary_dilation(
            no_data,
            structure=np.ones((3, 3), dtype=bool),
            iterations=margin_px,
        )
    invalid |= no_data
    return invalid


# =============================================================================
# Diagonal split search
# =============================================================================
def _class_proportions(labels: np.ndarray, classes: Sequence[int]) -> np.ndarray:
    if labels.size == 0:
        return np.zeros(len(classes), dtype=np.float64)
    return np.asarray(
        [np.count_nonzero(labels == cls) / labels.size for cls in classes],
        dtype=np.float64,
    )


def _threshold_between(values: np.ndarray, index: int) -> float:
    """Return a threshold between sorted values[index-1] and values[index]."""
    n = values.size
    if n == 0:
        raise ValueError("Cannot create a threshold from an empty array.")
    if index <= 0:
        return float(np.nextafter(values[0], -np.inf))
    if index >= n:
        return float(np.nextafter(values[-1], np.inf))
    left = float(values[index - 1])
    right = float(values[index])
    if right > left:
        return (left + right) / 2.0
    return left


def _segment_class_proportions(
    cumulative_counts: np.ndarray,
    start: int,
    end: int,
) -> np.ndarray:
    size = end - start
    if size <= 0:
        return np.zeros(cumulative_counts.shape[0], dtype=np.float64)
    counts = cumulative_counts[:, end] - cumulative_counts[:, start]
    return counts.astype(np.float64) / size


def search_diagonal_split(
    label_map: np.ndarray,
    valid_mask: np.ndarray,
    target_ratios: Mapping[str, float] = TARGET_SPLIT_RATIOS,
    theta_range_deg: Tuple[float, float] = THETA_RANGE_DEG,
    theta_step_deg: float = THETA_STEP_DEG,
    balance_classes: Sequence[int] = BALANCE_TARGET_CLASSES,
    min_split_ratio: float = MIN_SPLIT_RATIO,
    steps: int = RATIO_SEARCH_STEPS,
    ratio_score_weight: float = RATIO_SCORE_WEIGHT,
) -> Dict[str, object]:
    if not np.isclose(sum(target_ratios.values()), 1.0):
        raise ValueError(f"Target split ratios must sum to 1: {target_ratios}")

    h, w = label_map.shape
    yy, xx = np.mgrid[0:h, 0:w]
    labels_valid = label_map[valid_mask]
    x_valid = xx[valid_mask].astype(np.float64)
    y_valid = yy[valid_mask].astype(np.float64)

    if labels_valid.size == 0:
        raise ValueError("No valid pixels remain after invalid-mask construction.")

    overall_props = _class_proportions(labels_valid, balance_classes)
    theta_values = np.arange(theta_range_deg[0], theta_range_deg[1], theta_step_deg)
    r1_values = np.linspace(min_split_ratio, 1.0 - 2.0 * min_split_ratio, steps)

    best: Optional[Dict[str, object]] = None
    eps = 1e-12

    for theta_deg in theta_values:
        theta = np.deg2rad(theta_deg)
        projection = x_valid * np.cos(theta) + y_valid * np.sin(theta)
        order = np.argsort(projection, kind="mergesort")
        projection_sorted = projection[order]
        labels_sorted = labels_valid[order]
        n = labels_sorted.size

        cumulative = np.zeros((len(balance_classes), n + 1), dtype=np.int64)
        for class_index, class_id in enumerate(balance_classes):
            cumulative[class_index, 1:] = np.cumsum(labels_sorted == class_id)

        for r1 in r1_values:
            r2_values = np.linspace(
                r1 + min_split_ratio,
                1.0 - min_split_ratio,
                steps,
            )
            for r2 in r2_values:
                idx1 = int(round(r1 * n))
                idx2 = int(round(r2 * n))
                idx1 = int(np.clip(idx1, 1, n - 2))
                idx2 = int(np.clip(idx2, idx1 + 1, n - 1))

                region_props = {
                    0: _segment_class_proportions(cumulative, 0, idx1),
                    1: _segment_class_proportions(cumulative, idx1, idx2),
                    2: _segment_class_proportions(cumulative, idx2, n),
                }
                region_ratios = {
                    0: idx1 / n,
                    1: (idx2 - idx1) / n,
                    2: (n - idx2) / n,
                }

                for orientation in ("low_to_high", "high_to_low"):
                    if orientation == "low_to_high":
                        split_regions = {"train": 0, "val": 1, "test": 2}
                    else:
                        split_regions = {"train": 2, "val": 1, "test": 0}

                    class_score = 0.0
                    ratio_score = 0.0
                    props: Dict[str, List[float]] = {}
                    actual_ratios: Dict[str, float] = {}

                    for split_name, region_index in split_regions.items():
                        split_props = region_props[region_index]
                        props[split_name] = split_props.tolist()
                        actual_ratio = region_ratios[region_index]
                        actual_ratios[split_name] = float(actual_ratio)

                        class_score += float(
                            np.sum(
                                np.abs(split_props - overall_props)
                                / (overall_props + eps)
                            )
                        )
                        ratio_score += abs(
                            actual_ratio - target_ratios[split_name]
                        ) / target_ratios[split_name]

                    score = class_score + ratio_score_weight * ratio_score
                    if best is None or score < float(best["score"]):
                        best = {
                            "score": float(score),
                            "class_score": float(class_score),
                            "ratio_score": float(ratio_score),
                            "theta_deg": float(theta_deg),
                            "c1": _threshold_between(projection_sorted, idx1),
                            "c2": _threshold_between(projection_sorted, idx2),
                            "orientation": orientation,
                            "props_before_cut_filter": props,
                            "actual_ratios_before_cut_filter": actual_ratios,
                            "overall_props": overall_props.tolist(),
                        }

    if best is None:
        raise RuntimeError("Failed to find a diagonal split.")

    theta = np.deg2rad(float(best["theta_deg"]))
    projection_full = xx * np.cos(theta) + yy * np.sin(theta)
    region_id = np.zeros((h, w), dtype=np.uint8)
    region_id[projection_full > float(best["c1"])] = 1
    region_id[projection_full > float(best["c2"])] = 2

    if best["orientation"] == "low_to_high":
        split_region_masks = {
            "train": region_id == 0,
            "val": region_id == 1,
            "test": region_id == 2,
        }
    else:
        split_region_masks = {
            "train": region_id == 2,
            "val": region_id == 1,
            "test": region_id == 0,
        }

    best["region_id"] = region_id
    best["split_region_masks"] = split_region_masks
    return best


def ignore_filaments_crossing_split_boundaries(
    label_map: np.ndarray,
    valid_mask: np.ndarray,
    region_id: np.ndarray,
    dilation_px: int = CUT_COMPONENT_DILATION_PX,
) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """Ignore whole physical filaments crossing train/val/test boundaries.

    Classes 1 and 2 are merged before connected-component labelling, because
    they can be different confidence regions of the same physical filament.
    Eight-neighbour connectivity preserves diagonal one-pixel structures.
    """
    filament = (
        ((label_map == CLASS1_LABEL) | (label_map == CLASS2_LABEL))
        & valid_mask
    )
    connectivity_8 = np.ones((3, 3), dtype=np.uint8)
    components, n_components = ndi_label(filament, structure=connectivity_8)

    crossing_mask = np.zeros(label_map.shape, dtype=bool)
    crossing_component_count = 0

    for component_id in range(1, n_components + 1):
        component = components == component_id
        if np.unique(region_id[component]).size > 1:
            crossing_mask |= component
            crossing_component_count += 1

    crossing_filament_pixels = int(np.count_nonzero(crossing_mask))
    ignore_mask = crossing_mask
    if dilation_px > 0 and crossing_mask.any():
        ignore_mask = binary_dilation(
            crossing_mask,
            structure=np.ones((3, 3), dtype=bool),
            iterations=dilation_px,
        )

    refined_valid_mask = valid_mask & ~ignore_mask
    return (
        refined_valid_mask,
        ignore_mask,
        crossing_component_count,
        crossing_filament_pixels,
    )


# =============================================================================
# Image conversion and patch extraction
# =============================================================================
def make_rgb_column_density_image(
    col_density: np.ndarray,
    train_valid_mask: np.ndarray,
    low_percentile: float = NORMALIZE_LOW_PERCENTILE,
    high_percentile: float = NORMALIZE_HIGH_PERCENTILE,
) -> Tuple[np.ndarray, Dict[str, float]]:
    positive = np.isfinite(col_density) & (col_density > 0)
    log_map = np.full(col_density.shape, np.nan, dtype=np.float64)
    log_map[positive] = np.log10(col_density[positive])

    train_values = log_map[train_valid_mask & positive]
    if train_values.size == 0:
        raise ValueError("No valid train pixels available for normalization.")

    value_low, value_high = np.percentile(
        train_values,
        [low_percentile, high_percentile],
    )
    if not np.isfinite(value_low) or not np.isfinite(value_high):
        raise ValueError("Non-finite normalization limits were calculated.")
    if value_high <= value_low:
        value_high = value_low + 1e-6

    scaled = (log_map - value_low) / (value_high - value_low)
    scaled = np.clip(scaled, 0.0, 1.0)
    scaled[~positive] = 0.0
    gray = np.rint(scaled * 255.0).astype(np.uint8)
    rgb = np.repeat(gray[..., None], 3, axis=2)

    normalization = {
        "transform": "log10",
        "source_split": "train",
        "low_percentile": float(low_percentile),
        "high_percentile": float(high_percentile),
        "log10_value_low": float(value_low),
        "log10_value_high": float(value_high),
        "output_min": 0,
        "output_max": 255,
        "channels": 3,
    }
    return rgb, normalization


def sliding_origins(length: int, patch_size: int, stride: int) -> List[int]:
    if patch_size > length:
        return []
    origins = list(range(0, length - patch_size + 1, stride))
    final_origin = length - patch_size
    if not origins or origins[-1] != final_origin:
        origins.append(final_origin)
    return origins


def enumerate_patch_records(
    split_name: str,
    region: str,
    split_region_mask: np.ndarray,
    label_map_masked: np.ndarray,
    refined_valid_mask: np.ndarray,
    patch_size: int,
    stride: int,
    min_valid_ratio: float,
) -> List[PatchRecord]:
    h, w = label_map_masked.shape
    y_origins = sliding_origins(h, patch_size, stride)
    x_origins = sliding_origins(w, patch_size, stride)
    records: List[PatchRecord] = []

    for y in y_origins:
        y2 = y + patch_size
        for x in x_origins:
            x2 = x + patch_size

            # Prevent input-context leakage across train/val/test.
            if not np.all(split_region_mask[y:y2, x:x2]):
                continue

            valid_ratio = float(np.mean(refined_valid_mask[y:y2, x:x2]))
            if valid_ratio < min_valid_ratio:
                continue

            mask_patch = label_map_masked[y:y2, x:x2]
            counts = {
                label: int(np.count_nonzero(mask_patch == label))
                for label in (0, 1, 2, 255)
            }
            stem = f"{region}_{split_name}_y{y:05d}_x{x:05d}"
            records.append(
                PatchRecord(
                    split=split_name,
                    region=region,
                    stem=stem,
                    y=y,
                    x=x,
                    valid_ratio=valid_ratio,
                    background_pixels=counts[0],
                    class1_pixels=counts[1],
                    class2_pixels=counts[2],
                    ignore_pixels=counts[255],
                )
            )

    return records


def optionally_limit_train_background_patches(
    records: List[PatchRecord],
    max_background_to_positive: Optional[float],
    seed: int,
) -> List[PatchRecord]:
    if max_background_to_positive is None:
        return records
    if max_background_to_positive < 0:
        raise ValueError("Background-to-positive ratio must be non-negative.")

    positives: List[PatchRecord] = []
    backgrounds: List[PatchRecord] = []
    for record in records:
        if record.filament_pixels >= MIN_FILAMENT_PIXELS_FOR_POSITIVE:
            positives.append(record)
        else:
            backgrounds.append(record)

    if not positives:
        return records

    max_background = int(round(len(positives) * max_background_to_positive))
    if len(backgrounds) <= max_background:
        return records

    rng = np.random.default_rng(seed)
    selected_indices = rng.choice(
        len(backgrounds),
        size=max_background,
        replace=False,
    )
    selected_backgrounds = [backgrounds[int(i)] for i in selected_indices]
    return sorted(positives + selected_backgrounds, key=lambda item: (item.y, item.x))


def assign_storage_subdirectories(
    records: Sequence[PatchRecord],
    files_per_subdir: int,
) -> List[PatchRecord]:
    if files_per_subdir <= 0:
        raise ValueError("files_per_subdir must be a positive integer.")

    ordered = sorted(records, key=lambda item: (item.y, item.x, item.stem))
    return [
        replace(record, storage_subdir=f"{index // files_per_subdir:04d}")
        for index, record in enumerate(ordered)
    ]


def save_patch_files(
    root: Path,
    rgb_image: np.ndarray,
    label_map_masked: np.ndarray,
    records: Iterable[PatchRecord],
    patch_size: int,
) -> None:
    for record in records:
        y1, x1 = record.y, record.x
        y2, x2 = y1 + patch_size, x1 + patch_size
        image_patch = rgb_image[y1:y2, x1:x2]
        mask_patch = label_map_masked[y1:y2, x1:x2]

        image_dir = (
            root / record.split / "in" / record.region / record.storage_subdir
        )
        mask_dir = (
            root / record.split / "out" / record.region / record.storage_subdir
        )
        image_dir.mkdir(parents=True, exist_ok=True)
        mask_dir.mkdir(parents=True, exist_ok=True)

        image_path = image_dir / f"{record.stem}_in.png"
        mask_path = mask_dir / f"{record.stem}_out.png"

        Image.fromarray(image_patch).save(image_path, compress_level=1)
        Image.fromarray(mask_patch.astype(np.uint8)).save(
            mask_path,
            compress_level=1,
        )


def records_to_stats(records: Sequence[PatchRecord]) -> Dict[str, object]:
    pixel_counts = {
        "background": int(sum(r.background_pixels for r in records)),
        "class1": int(sum(r.class1_pixels for r in records)),
        "class2": int(sum(r.class2_pixels for r in records)),
        "ignore": int(sum(r.ignore_pixels for r in records)),
    }
    valid_total = (
        pixel_counts["background"]
        + pixel_counts["class1"]
        + pixel_counts["class2"]
    )
    valid_fractions = {
        name: (count / valid_total if valid_total else 0.0)
        for name, count in pixel_counts.items()
        if name != "ignore"
    }
    return {
        "patch_count": len(records),
        "pixel_counts": pixel_counts,
        "valid_pixel_fractions": valid_fractions,
        "filament_patch_count": int(
            sum(
                r.filament_pixels >= MIN_FILAMENT_PIXELS_FOR_POSITIVE
                for r in records
            )
        ),
        "background_only_patch_count": int(
            sum(
                r.filament_pixels < MIN_FILAMENT_PIXELS_FOR_POSITIVE
                for r in records
            )
        ),
    }


def calculate_final_split_stats(
    label_map: np.ndarray,
    split_region_masks: Mapping[str, np.ndarray],
    refined_valid_mask: np.ndarray,
) -> Dict[str, object]:
    result: Dict[str, object] = {}
    total_valid = int(np.count_nonzero(refined_valid_mask))

    for split_name, region_mask in split_region_masks.items():
        mask = region_mask & refined_valid_mask
        labels = label_map[mask]
        counts = {str(c): int(np.count_nonzero(labels == c)) for c in (0, 1, 2)}
        n = labels.size
        result[split_name] = {
            "valid_pixels": int(n),
            "fraction_of_all_valid_pixels": n / total_valid if total_valid else 0.0,
            "class_counts": counts,
            "class_fractions": {
                str(c): counts[str(c)] / n if n else 0.0 for c in (0, 1, 2)
            },
        }
    return result


# =============================================================================
# Metadata and visualization
# =============================================================================
def save_manifest(path: Path, records: Sequence[PatchRecord]) -> None:
    fieldnames = list(asdict(records[0]).keys()) if records else [
        "split",
        "region",
        "stem",
        "y",
        "x",
        "valid_ratio",
        "background_pixels",
        "class1_pixels",
        "class2_pixels",
        "ignore_pixels",
        "storage_subdir",
    ]
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow(asdict(record))


def save_overview(
    col_density: np.ndarray,
    label_map_masked: np.ndarray,
    split_region_masks: Mapping[str, np.ndarray],
    out_path: Path,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; split overview was not generated.")
        return

    finite_positive = np.isfinite(col_density) & (col_density > 0)
    display = np.full(col_density.shape, np.nan, dtype=np.float64)
    display[finite_positive] = np.log10(col_density[finite_positive])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    axes[0].imshow(display, origin="lower", cmap="inferno")
    overlay = np.zeros(label_map_masked.shape + (4,), dtype=np.float32)
    overlay[split_region_masks["train"]] = (0.2, 0.6, 1.0, 0.30)
    overlay[split_region_masks["val"]] = (0.2, 1.0, 0.4, 0.30)
    overlay[split_region_masks["test"]] = (1.0, 0.3, 0.2, 0.30)
    axes[0].imshow(overlay, origin="lower")
    axes[0].set_title("Spatial split: blue=train, green=val, red=test")
    axes[0].axis("off")

    masked_for_display = label_map_masked.astype(np.float32)
    masked_for_display[masked_for_display == IGNORE_LABEL] = np.nan
    axes[1].imshow(masked_for_display, origin="lower", cmap="viridis", vmin=0, vmax=2)
    axes[1].set_title("3-class labels; ignored pixels are transparent")
    axes[1].axis("off")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================
def run(
    region_name: str = REGION_NAME,
    column_density_fits: str = COLUMN_DENSITY_FITS,
    label_map_path: str = LABEL_MAP_PATH,
    dataset_root: str = DATASET_ROOT,
    patch_size: int = PATCH_SIZE,
    train_stride: int = TRAIN_STRIDE,
    eval_stride: int = EVAL_STRIDE,
    files_per_subdir: int = FILES_PER_SUBDIR,
    clean_region: bool = CLEAN_REGION_OUTPUT,
    column_density_hdu: int = 0,
    label_hdu: int = 0,
) -> Dict[str, object]:
    if patch_size <= 0 or train_stride <= 0 or eval_stride <= 0:
        raise ValueError("Patch size and strides must be positive integers.")
    if files_per_subdir <= 0:
        raise ValueError("files_per_subdir must be a positive integer.")

    region = normalize_region_name(region_name)
    root = Path(dataset_root)
    metadata_dir = prepare_output_dirs(root, region, clean_region=clean_region)

    print(f"Loading region={region} FITS image and 3-class labels...")
    col_density, label_map = load_data(
        column_density_fits,
        label_map_path,
        column_density_hdu=column_density_hdu,
        label_hdu=label_hdu,
    )
    invalid_mask = build_invalid_mask(col_density, MARGIN_PX)
    valid_mask = ~invalid_mask

    print("Searching spatial diagonal split...")
    split_result = search_diagonal_split(label_map, valid_mask)
    region_id = np.asarray(split_result["region_id"])
    split_region_masks = split_result["split_region_masks"]

    print("Ignoring connected filaments crossing split boundaries...")
    (
        refined_valid_mask,
        crossing_ignore_mask,
        crossing_component_count,
        crossing_filament_pixels,
    ) = ignore_filaments_crossing_split_boundaries(
        label_map,
        valid_mask,
        region_id,
    )

    label_map_masked = label_map.copy()
    label_map_masked[~refined_valid_mask] = IGNORE_LABEL

    train_valid_mask = split_region_masks["train"] & refined_valid_mask
    rgb_image, normalization = make_rgb_column_density_image(
        col_density,
        train_valid_mask,
    )

    final_split_stats = calculate_final_split_stats(
        label_map,
        split_region_masks,
        refined_valid_mask,
    )

    all_records: List[PatchRecord] = []
    records_by_split: Dict[str, List[PatchRecord]] = {}
    for split_name in ("train", "val", "test"):
        stride = train_stride if split_name == "train" else eval_stride
        records = enumerate_patch_records(
            split_name=split_name,
            region=region,
            split_region_mask=split_region_masks[split_name],
            label_map_masked=label_map_masked,
            refined_valid_mask=refined_valid_mask,
            patch_size=patch_size,
            stride=stride,
            min_valid_ratio=MIN_PATCH_VALID_RATIO,
        )
        if split_name == "train":
            records = optionally_limit_train_background_patches(
                records,
                MAX_TRAIN_BACKGROUND_TO_POSITIVE,
                RANDOM_SEED,
            )
        records = assign_storage_subdirectories(records, files_per_subdir)
        if not records:
            raise RuntimeError(
                f"No {split_name} patches were generated. "
                "The split band may be narrower than PATCH_SIZE, or "
                "MIN_PATCH_VALID_RATIO may be too strict."
            )
        records_by_split[split_name] = records
        all_records.extend(records)

    print("Saving PNG image/mask patches...")
    save_patch_files(root, rgb_image, label_map_masked, all_records, patch_size)

    np.save(metadata_dir / "train_region_mask.npy", split_region_masks["train"])
    np.save(metadata_dir / "val_region_mask.npy", split_region_masks["val"])
    np.save(metadata_dir / "test_region_mask.npy", split_region_masks["test"])
    np.save(metadata_dir / "valid_mask_after_cut_filter.npy", refined_valid_mask)
    np.save(metadata_dir / "crossing_filament_ignore_mask.npy", crossing_ignore_mask)
    np.save(metadata_dir / "label_map_masked.npy", label_map_masked)

    save_manifest(metadata_dir / "patch_manifest.csv", all_records)
    save_overview(
        col_density,
        label_map_masked,
        split_region_masks,
        metadata_dir / "split_overview.png",
    )

    patch_stats = {
        split_name: records_to_stats(records)
        for split_name, records in records_by_split.items()
    }

    serializable_split = {
        key: value
        for key, value in split_result.items()
        if key not in {"region_id", "split_region_masks"}
    }
    summary: Dict[str, object] = {
        "input": {
            "region": region,
            "column_density_fits": column_density_fits,
            "column_density_hdu": column_density_hdu,
            "label_map_path": label_map_path,
            "label_hdu": label_hdu,
            "shape": list(label_map.shape),
            "label_definition": {
                "0": "background",
                "1": "normal-threshold FilFinder/DisPerSE union",
                "2": "relaxed-threshold result minus class 1",
                "255": "ignore",
            },
        },
        "output": {
            "dataset_root": str(root),
            "patch_size": patch_size,
            "train_stride": train_stride,
            "eval_stride": eval_stride,
            "files_per_subdir": files_per_subdir,
            "directory_pattern": "{split}/{in|out}/{region}/{chunk:04d}/",
            "image_suffix": "_in.png",
            "mask_suffix": "_out.png",
        },
        "split_search": serializable_split,
        "final_split_stats_after_cut_filter": final_split_stats,
        "crossing_filament_filter": {
            "crossing_component_count": crossing_component_count,
            "crossing_filament_pixels": crossing_filament_pixels,
            "ignored_pixels_after_dilation": int(np.count_nonzero(crossing_ignore_mask)),
            "dilation_px": CUT_COMPONENT_DILATION_PX,
        },
        "normalization": normalization,
        "patch_stats": patch_stats,
    }

    with (metadata_dir / "dataset_summary.json").open("w", encoding="utf-8") as file:
        json.dump(summary, file, ensure_ascii=False, indent=2)

    print("\nCompleted.")
    print(f"Dataset root: {root}")
    print(f"Region: {region}")
    print(f"Metadata: {metadata_dir}")
    print(
        "Best split: "
        f"theta={split_result['theta_deg']:.2f} deg, "
        f"orientation={split_result['orientation']}, "
        f"score={split_result['score']:.4f}"
    )
    print(
        "Crossing filaments ignored: "
        f"components={crossing_component_count}, "
        f"filament_pixels={crossing_filament_pixels}"
    )
    for split_name in ("train", "val", "test"):
        stats = patch_stats[split_name]
        fractions = stats["valid_pixel_fractions"]
        print(
            f"{split_name}: patches={stats['patch_count']}, "
            f"class fractions "
            f"0={fractions['background']:.6f}, "
            f"1={fractions['class1']:.6f}, "
            f"2={fractions['class2']:.6f}"
        )

    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", default=REGION_NAME)
    parser.add_argument("--column-density-fits", default=COLUMN_DENSITY_FITS)
    parser.add_argument("--label-map", default=LABEL_MAP_PATH)
    parser.add_argument("--dataset-root", default=DATASET_ROOT)
    parser.add_argument("--patch-size", type=int, default=PATCH_SIZE)
    parser.add_argument("--train-stride", type=int, default=TRAIN_STRIDE)
    parser.add_argument("--eval-stride", type=int, default=EVAL_STRIDE)
    parser.add_argument("--files-per-subdir", type=int, default=FILES_PER_SUBDIR)
    parser.add_argument("--column-density-hdu", type=int, default=0)
    parser.add_argument("--label-hdu", type=int, default=0)
    parser.add_argument(
        "--no-clean-region",
        action="store_true",
        help=(
            "Do not delete this region's existing output first. Other regions "
            "are never deleted, even without this option."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(
        region_name=args.region,
        column_density_fits=args.column_density_fits,
        label_map_path=args.label_map,
        dataset_root=args.dataset_root,
        patch_size=args.patch_size,
        train_stride=args.train_stride,
        eval_stride=args.eval_stride,
        files_per_subdir=args.files_per_subdir,
        clean_region=not args.no_clean_region,
        column_density_hdu=args.column_density_hdu,
        label_hdu=args.label_hdu,
    )
