#!/usr/bin/env python3
"""Split every cloud region into train/val/test without filament leakage.

The preferred split is one connected L/T-shaped spatial partition, for example:

    +-----------+-----+
    |           | val |
    |   train   +-----+
    |           |test |
    +-----------+-----+

The split unit is a non-overlapping 256x256 patch. Class-1 and class-2 pixels
are combined and labelled with 8-neighbour connected components. Every connected filament is assigned to exactly one split. A patch that
contains several independent filaments is kept only when all of those
filaments are assigned to the same split. If they are assigned to different
splits, only that conflicting patch is omitted; the independent filaments are
not merged permanently. Therefore the final layout remains as close as
possible to the L/T shape while avoiding transitive over-merging and leakage.

Each configured region must contribute class-1-positive patches to train,
validation and test.
"""

from __future__ import annotations

import csv
import json
import re
import shutil
import traceback
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np
from PIL import Image
from scipy.ndimage import binary_dilation, label as ndi_label


# =============================================================================
# USER SETTINGS
# =============================================================================

DATASET_ROOT = "/mnt/sp4t/fujiwara.r/split"
RUN_ALL_CONFIGURED_REGIONS = True
SINGLE_REGION_NAME = "aquila"
CONTINUE_ON_REGION_ERROR = True

REGION_DATASETS = {
    "aquila": {
        "enabled": True,
        "column_density_fits": (
            "/mnt/sp4t/fujiwara.r/HGBS_aquilaM2_hires_column_density_map.fits"
        ),
        "label_map_path": (
            "/mnt/sp4t/fujiwara.r/HGBS_aquilaM2_hires_column_density_map/vis_3class_map/HGBS_aquilaM2_hires_column_density_map_label3class.fits"
        ),
    },

    "orionb": {
        "enabled": True,
        "column_density_fits": (
            "/mnt/sp4t/fujiwara.r/HGBS_orionB_hires_column_density_map.fits"
        ),
        "label_map_path": (
            "/mnt/sp4t/fujiwara.r/HGBS_orionB_hires_column_density_map/vis_3class_map/HGBS_orionB_hires_column_density_map_label3class.fits"
        ),
    },

    "ic5146": {
        "enabled": True,
        "column_density_fits": (
            "/mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map.fits"
        ),
        "label_map_path": (
           "/mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map/vis_3class_map/HGBS_ic5146_hires_column_density_map_label3class.fits"
        ),
    },

    # 元画像名が HGBS_tauN3... でも、保存領域名は taurus にできます。
    "taurus": {
        "enabled": True,
        "column_density_fits": (
            "/mnt/sp4t/fujiwara.r/HGBS_tauN3_hires_column_density_map.fits"
        ),
        "label_map_path": (
            "/mnt/sp4t/fujiwara.r/HGBS_tauN3_hires_column_density_map/vis_3class_map/HGBS_tauN3_hires_column_density_map_label3class.fits"
        ),
    },

    "polaris": {
        "enabled": True,
        "column_density_fits": (
            "/mnt/sp4t/fujiwara.r/HGBS_polaris_hires_column_density_map.fits"
        ),
        "label_map_path": (
            "/mnt/sp4t/fujiwara.r/HGBS_polaris_hires_column_density_map/vis_3class_map/HGBS_polaris_hires_column_density_map_label3class.fits"
        ),
    },
}



# =============================================================================
# LABEL / PATCH SETTINGS
# =============================================================================

BACKGROUND_LABEL = 0
CLASS1_LABEL = 1
CLASS2_LABEL = 2
IGNORE_LABEL = 255
ALLOWED_LABELS = {0, 1, 2, 255}

PATCH_SIZE = 256
PATCH_STRIDE = 256  # Must equal PATCH_SIZE to prevent pixel overlap.

# Keep 256x256 non-overlapping patches, but search different grid origins when
# the ordinary x=0, y=0 grid cannot produce a valid three-way split.
PATCH_GRID_OFFSETS = (0, 64, 128, 192)

# Successful regions remain on the ordinary grid. Alternative offsets are
# searched only when offset (0, 0) fails.
SEARCH_GRID_OFFSETS_ON_FAILURE_ONLY = True
MIN_PATCH_VALID_RATIO = 0.80
MARGIN_PX = 32
FILES_PER_SUBDIR = 1000

TARGET_SPLIT_RATIOS: Mapping[str, float] = {
    "train": 0.60,
    "val": 0.20,
    "test": 0.20,
}
SPLIT_ORDER = ("train", "val", "test")
SPLIT_TO_ID = {"train": 0, "val": 1, "test": 2}
ID_TO_SPLIT = {value: key for key, value in SPLIT_TO_ID.items()}

# Every cloud region must satisfy these final conditions.
MIN_PATCHES_PER_SPLIT = {"train": 10, "val": 1, "test": 1}
MIN_CLASS1_PATCHES_PER_SPLIT = {"train": 1, "val": 1, "test": 1}
MIN_CLASS1_PIXELS_PER_SPLIT = {"train": 16, "val": 16, "test": 16}

REQUIRE_CLASS2_IN_EVERY_SPLIT = False
MIN_CLASS2_PIXELS_PER_SPLIT = 1

FILAMENT_COMPONENT_MODE = "union"  # class 1 + class 2
FILAMENT_COMPONENT_LINK_GAP_PX = 0
MIN_FILAMENT_COMPONENT_PIXELS = 1

ASSIGNMENT_TRIALS = 20_000
RANDOM_SEED = 0

# When a region has only a small number of independent components containing
# class 1, evaluate every train/val/test assignment exactly.
#
# Polaris has 5 class-1 components:
#     3 ** 5 = 243 assignments
# so exhaustive search is inexpensive and avoids greedy-search misses.
EXACT_CLASS1_ASSIGNMENT_MAX_GROUPS = 6

# =============================================================================
# L/T-SHAPE PREFERENCE
# =============================================================================
PREFERRED_SPATIAL_LAYOUT = "l_t_shape"

L_SHAPE_CUT_FRACTIONS = (
    0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75,
)
SECONDARY_CUT_FRACTIONS = (
    0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
)

SPATIAL_SHAPE_WEIGHT = 8.0
PATCH_RATIO_WEIGHT = 3.0
CLASS1_BALANCE_WEIGHT = 5.0
CLASS2_BALANCE_WEIGHT = 2.0

# None means there is no hard rejection threshold. The mismatch percentage is
# still minimized and written to dataset_summary.json.
MAX_POSITIVE_PATCH_MISMATCH_FRACTION = None


# A patch containing independent components assigned to different splits is
# omitted instead of merging those components into one permanent group.
CONFLICT_PATCH_DROP_WEIGHT = 12.0

# None disables a hard limit. The fraction is still minimized and recorded.
MAX_CONFLICT_PATCH_DROP_FRACTION = None

NORMALIZE_LOW_PERCENTILE = 1.0
NORMALIZE_HIGH_PERCENTILE = 99.0
CLEAN_REGION_OUTPUT = True


@dataclass(frozen=True)
class PatchRecord:
    region: str
    split: str
    stem: str
    y: int
    x: int
    valid_ratio: float
    background_pixels: int
    class1_pixels: int
    class2_pixels: int
    ignore_pixels: int
    component_ids: Tuple[int, ...]
    locked_group_id: int
    storage_subdir: str = ""

    @property
    def filament_pixels(self) -> int:
        return self.class1_pixels + self.class2_pixels

    @property
    def is_background_only(self) -> bool:
        return self.filament_pixels == 0


@dataclass
class LockedGroup:
    group_id: int
    component_ids: Tuple[int, ...]
    patch_indices: List[int]
    patch_count: int
    class1_patch_count: int
    class1_pixels: int
    class2_pixels: int
    valid_pixels: int


class DisjointSet:
    def __init__(self, size: int) -> None:
        self.parent = np.arange(size + 1, dtype=np.int32)
        self.rank = np.zeros(size + 1, dtype=np.int8)

    def find(self, value: int) -> int:
        parent = int(self.parent[value])
        if parent != value:
            self.parent[value] = self.find(parent)
        return int(self.parent[value])

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1


# =============================================================================
# INPUT
# =============================================================================

def normalize_region_name(region_name: str) -> str:
    region = region_name.strip().lower().replace(" ", "_")
    if not region or not re.fullmatch(r"[a-z0-9_-]+", region):
        raise ValueError(f"Invalid region name: {region_name!r}")
    return region


def read_fits_primary(path: str) -> np.ndarray:
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise ImportError("Install astropy with: pip install astropy") from exc

    with fits.open(path, memmap=True) as hdul:
        if not hdul or hdul[0].data is None:
            raise ValueError(f"Primary HDU contains no image data: {path}")
        return np.array(hdul[0].data, copy=True)


def coerce_label_array(raw: np.ndarray, source: str) -> np.ndarray:
    array = np.squeeze(np.asarray(raw))
    if array.ndim != 2:
        raise ValueError(
            f"Label map must be 2-D: {source}, shape={array.shape}"
        )

    if np.issubdtype(array.dtype, np.floating):
        finite = np.isfinite(array)
        values = array[finite]
        rounded = np.rint(values)
        if values.size and not np.allclose(values, rounded, atol=1e-6):
            raise ValueError(f"Non-integer label values in {source}")
        result = np.full(array.shape, IGNORE_LABEL, dtype=np.uint8)
        if rounded.size:
            numeric = rounded.astype(np.int64)
            unexpected = set(int(v) for v in np.unique(numeric)) - ALLOWED_LABELS
            if unexpected:
                raise ValueError(
                    f"Unexpected labels {sorted(unexpected)} in {source}"
                )
            result[finite] = numeric.astype(np.uint8)
        return result

    numeric = array.astype(np.int64, copy=False)
    unexpected = set(int(v) for v in np.unique(numeric)) - ALLOWED_LABELS
    if unexpected:
        raise ValueError(f"Unexpected labels {sorted(unexpected)} in {source}")
    return numeric.astype(np.uint8)


def load_data(
    column_density_fits: str,
    label_map_path: str,
) -> Tuple[np.ndarray, np.ndarray]:
    col_density = np.squeeze(
        np.asarray(read_fits_primary(column_density_fits), dtype=np.float64)
    )
    if col_density.ndim != 2:
        raise ValueError(f"Column-density image is not 2-D: {col_density.shape}")

    label_path = Path(label_map_path)
    suffix = label_path.suffix.lower()
    if suffix == ".npy":
        raw_label = np.load(label_path)
    elif suffix in {".fits", ".fit", ".fts"}:
        raw_label = read_fits_primary(str(label_path))
    else:
        raise ValueError(f"Unsupported label type: {label_path}")

    label_map = coerce_label_array(raw_label, str(label_path))
    if col_density.shape != label_map.shape:
        raise ValueError(
            f"Shape mismatch: image={col_density.shape}, label={label_map.shape}"
        )
    return col_density, label_map


def build_invalid_mask(col_density: np.ndarray) -> np.ndarray:
    h, w = col_density.shape
    invalid = np.zeros((h, w), dtype=bool)
    if MARGIN_PX > 0:
        invalid[:MARGIN_PX, :] = True
        invalid[-MARGIN_PX:, :] = True
        invalid[:, :MARGIN_PX] = True
        invalid[:, -MARGIN_PX:] = True

    no_data = ~np.isfinite(col_density) | (col_density <= 0)
    if no_data.any() and MARGIN_PX > 0:
        no_data = binary_dilation(
            no_data,
            structure=np.ones((3, 3), dtype=bool),
            iterations=MARGIN_PX,
        )
    invalid |= no_data
    return invalid


# =============================================================================
# COMPONENTS AND PATCHES
# =============================================================================

def build_filament_component_map(
    label_map: np.ndarray,
    valid_mask: np.ndarray,
) -> Tuple[np.ndarray, int, Dict[str, object]]:
    connectivity_8 = np.ones((3, 3), dtype=np.uint8)

    if FILAMENT_COMPONENT_MODE == "union":
        filament = (
            (label_map == CLASS1_LABEL) | (label_map == CLASS2_LABEL)
        ) & valid_mask
    elif FILAMENT_COMPONENT_MODE == "class1":
        filament = (label_map == CLASS1_LABEL) & valid_mask
    else:
        raise ValueError("FILAMENT_COMPONENT_MODE must be 'union' or 'class1'")

    if FILAMENT_COMPONENT_LINK_GAP_PX > 0:
        filament = binary_dilation(
            filament,
            structure=connectivity_8.astype(bool),
            iterations=FILAMENT_COMPONENT_LINK_GAP_PX,
        ) & valid_mask

    raw_components, raw_count = ndi_label(filament, structure=connectivity_8)
    if raw_count == 0:
        raise RuntimeError("No filament component was found.")

    sizes = np.bincount(raw_components.ravel(), minlength=raw_count + 1)
    kept = np.flatnonzero(sizes >= MIN_FILAMENT_COMPONENT_PIXELS)
    kept = kept[kept != 0]

    lookup = np.zeros(raw_count + 1, dtype=np.int32)
    lookup[kept] = np.arange(1, len(kept) + 1, dtype=np.int32)
    component_map = lookup[raw_components]
    component_count = int(len(kept))

    metadata = {
        "mode": FILAMENT_COMPONENT_MODE,
        "connectivity": 8,
        "link_gap_px": FILAMENT_COMPONENT_LINK_GAP_PX,
        "raw_component_count": int(raw_count),
        "kept_component_count": component_count,
        "component_pixels": int(np.count_nonzero(component_map)),
    }
    return component_map, component_count, metadata


def sliding_origins(
    length: int,
    offset: int,
) -> List[int]:
    if PATCH_SIZE > length:
        return []
    if offset < 0 or offset >= PATCH_STRIDE:
        raise ValueError(
            f"Grid offset must satisfy 0 <= offset < {PATCH_STRIDE}: "
            f"{offset}"
        )
    if offset + PATCH_SIZE > length:
        return []
    return list(
        range(
            offset,
            length - PATCH_SIZE + 1,
            PATCH_STRIDE,
        )
    )


def enumerate_candidate_patches(
    region: str,
    label_map: np.ndarray,
    valid_mask: np.ndarray,
    component_map: np.ndarray,
    offset_x: int = 0,
    offset_y: int = 0,
) -> List[PatchRecord]:
    if PATCH_STRIDE != PATCH_SIZE:
        raise ValueError(
            "PATCH_STRIDE must equal PATCH_SIZE to prevent input-pixel overlap."
        )

    h, w = label_map.shape
    records: List[PatchRecord] = []

    for y in sliding_origins(h, offset_y):
        y2 = y + PATCH_SIZE
        for x in sliding_origins(w, offset_x):
            x2 = x + PATCH_SIZE
            valid_patch = valid_mask[y:y2, x:x2]
            valid_ratio = float(np.mean(valid_patch))
            if valid_ratio < MIN_PATCH_VALID_RATIO:
                continue

            mask_patch = label_map[y:y2, x:x2].copy()
            mask_patch[~valid_patch] = IGNORE_LABEL
            counts = {
                value: int(np.count_nonzero(mask_patch == value))
                for value in (0, 1, 2, 255)
            }
            ids = np.unique(component_map[y:y2, x:x2])
            ids = tuple(int(v) for v in ids if int(v) > 0)

            records.append(
                PatchRecord(
                    region=region,
                    split="",
                    stem=f"y{y:05d}_x{x:05d}",
                    y=y,
                    x=x,
                    valid_ratio=valid_ratio,
                    background_pixels=counts[0],
                    class1_pixels=counts[1],
                    class2_pixels=counts[2],
                    ignore_pixels=counts[255],
                    component_ids=ids,
                    locked_group_id=0,
                )
            )

    if not records:
        raise RuntimeError("No valid 256x256 patches were generated.")
    return records


def create_component_groups(
    records: Sequence[PatchRecord],
    component_count: int,
    label_map: np.ndarray,
    component_map: np.ndarray,
) -> Tuple[List[PatchRecord], List[LockedGroup], Dict[str, object]]:
    """Keep every connected component independent.

    The previous implementation merged all components that happened to occur
    in the same 256x256 patch. That relation is transitive and can collapse
    dozens of independent filaments into only two groups. Here each connected
    component remains its own locked group. A multi-component patch is handled
    later: it is kept only if all of its components receive the same split.
    """
    component_patch_indices: Dict[int, List[int]] = {
        component_id: []
        for component_id in range(1, component_count + 1)
    }

    updated: List[PatchRecord] = []
    multi_component_patch_count = 0

    for patch_index, record in enumerate(records):
        for component_id in record.component_ids:
            component_patch_indices[component_id].append(patch_index)

        if not record.component_ids:
            group_id = 0
        elif len(record.component_ids) == 1:
            group_id = int(record.component_ids[0])
        else:
            # -1 explicitly marks a patch that contains multiple independent
            # components. It is not itself a permanent locked group.
            group_id = -1
            multi_component_patch_count += 1

        updated.append(replace(record, locked_group_id=group_id))

    groups: List[LockedGroup] = []

    for component_id in range(1, component_count + 1):
        patch_indices = sorted(set(component_patch_indices[component_id]))
        if not patch_indices:
            # The component is not represented in any retained patch for this
            # grid alignment, so it cannot contribute to this split.
            continue

        component_pixels = component_map == component_id

        class1_pixels = int(
            np.count_nonzero(
                component_pixels & (label_map == CLASS1_LABEL)
            )
        )
        class2_pixels = int(
            np.count_nonzero(
                component_pixels & (label_map == CLASS2_LABEL)
            )
        )

        class1_patch_count = 0
        valid_pixels = 0

        for patch_index in patch_indices:
            record = updated[patch_index]
            y1, x1 = record.y, record.x
            y2, x2 = y1 + PATCH_SIZE, x1 + PATCH_SIZE

            component_patch = (
                component_map[y1:y2, x1:x2] == component_id
            )
            label_patch = label_map[y1:y2, x1:x2]

            if np.any(component_patch & (label_patch == CLASS1_LABEL)):
                class1_patch_count += 1

            valid_pixels += int(np.count_nonzero(component_patch))

        groups.append(
            LockedGroup(
                group_id=component_id,
                component_ids=(component_id,),
                patch_indices=patch_indices,
                patch_count=len(patch_indices),
                class1_patch_count=class1_patch_count,
                class1_pixels=class1_pixels,
                class2_pixels=class2_pixels,
                valid_pixels=valid_pixels,
            )
        )

    metadata = {
        "connected_component_count": component_count,
        "locked_group_count": len(groups),
        "locking_unit": "one connected component",
        "background_only_patch_count": int(
            sum(not record.component_ids for record in updated)
        ),
        "single_component_patch_count": int(
            sum(len(record.component_ids) == 1 for record in updated)
        ),
        "multi_component_patch_count": int(multi_component_patch_count),
        "components_with_class1": int(
            sum(group.class1_patch_count > 0 for group in groups)
        ),
    }
    return updated, groups, metadata


@dataclass(frozen=True)
class LShapeCandidate:
    family: str
    primary_cut: float
    secondary_cut: float
    region_to_split: Tuple[int, int, int]

    @property
    def description(self) -> str:
        mapping = "/".join(
            ID_TO_SPLIT[split_id] for split_id in self.region_to_split
        )
        return (
            f"{self.family}(primary={self.primary_cut:.1f}, "
            f"secondary={self.secondary_cut:.1f}); "
            f"regions={mapping}"
        )


def l_shape_base_region(
    record: PatchRecord,
    candidate: LShapeCandidate,
) -> int:
    centre_x = record.x + PATCH_SIZE / 2.0
    centre_y = record.y + PATCH_SIZE / 2.0
    a = candidate.primary_cut
    b = candidate.secondary_cut

    if candidate.family == "left_large_right_stacked":
        if centre_x < a:
            return 0
        return 1 if centre_y < b else 2

    if candidate.family == "right_large_left_stacked":
        if centre_x >= a:
            return 0
        return 1 if centre_y < b else 2

    if candidate.family == "top_large_bottom_side_by_side":
        if centre_y < a:
            return 0
        return 1 if centre_x < b else 2

    if candidate.family == "bottom_large_top_side_by_side":
        if centre_y >= a:
            return 0
        return 1 if centre_x < b else 2

    raise ValueError(f"Unknown L-shape family: {candidate.family}")


def l_shape_preferred_split_id(
    record: PatchRecord,
    candidate: LShapeCandidate,
) -> int:
    base_region = l_shape_base_region(record, candidate)
    return int(candidate.region_to_split[base_region])


def generate_l_shape_candidates(
    records: Sequence[PatchRecord],
) -> List[LShapeCandidate]:
    from itertools import permutations

    if not records:
        return []

    centres_x = np.asarray(
        [record.x + PATCH_SIZE / 2.0 for record in records],
        dtype=np.float64,
    )
    centres_y = np.asarray(
        [record.y + PATCH_SIZE / 2.0 for record in records],
        dtype=np.float64,
    )

    x_primary = sorted(
        set(
            float(np.quantile(centres_x, fraction))
            for fraction in L_SHAPE_CUT_FRACTIONS
        )
    )
    y_primary = sorted(
        set(
            float(np.quantile(centres_y, fraction))
            for fraction in L_SHAPE_CUT_FRACTIONS
        )
    )
    x_secondary = sorted(
        set(
            float(np.quantile(centres_x, fraction))
            for fraction in SECONDARY_CUT_FRACTIONS
        )
    )
    y_secondary = sorted(
        set(
            float(np.quantile(centres_y, fraction))
            for fraction in SECONDARY_CUT_FRACTIONS
        )
    )

    mappings = list(permutations((0, 1, 2)))
    candidates: List[LShapeCandidate] = []

    for primary in x_primary:
        for secondary in y_secondary:
            for family in (
                "left_large_right_stacked",
                "right_large_left_stacked",
            ):
                for mapping in mappings:
                    candidates.append(
                        LShapeCandidate(
                            family=family,
                            primary_cut=primary,
                            secondary_cut=secondary,
                            region_to_split=tuple(int(v) for v in mapping),
                        )
                    )

    for primary in y_primary:
        for secondary in x_secondary:
            for family in (
                "top_large_bottom_side_by_side",
                "bottom_large_top_side_by_side",
            ):
                for mapping in mappings:
                    candidates.append(
                        LShapeCandidate(
                            family=family,
                            primary_cut=primary,
                            secondary_cut=secondary,
                            region_to_split=tuple(int(v) for v in mapping),
                        )
                    )

    return candidates


def preferred_split_counts_for_group(
    group: LockedGroup,
    records: Sequence[PatchRecord],
    candidate: LShapeCandidate,
) -> np.ndarray:
    counts = np.zeros(3, dtype=np.int64)
    for patch_index in group.patch_indices:
        split_id = l_shape_preferred_split_id(
            records[patch_index],
            candidate,
        )
        counts[split_id] += 1
    return counts


def choose_group_split_with_shape_preference(
    group: LockedGroup,
    preferred_counts: np.ndarray,
    patch_totals: np.ndarray,
    class1_totals: np.ndarray,
    class2_totals: np.ndarray,
    ideal_patches: np.ndarray,
    ideal_class1: np.ndarray,
    ideal_class2: np.ndarray,
) -> int:
    choices: List[Tuple[float, int]] = []

    for split_id in range(3):
        moved_patches = group.patch_count - int(preferred_counts[split_id])
        spatial_cost = SPATIAL_SHAPE_WEIGHT * (
            moved_patches / max(group.patch_count, 1)
        )

        p = patch_totals.copy()
        c1 = class1_totals.copy()
        c2 = class2_totals.copy()
        p[split_id] += group.patch_count
        c1[split_id] += group.class1_pixels
        c2[split_id] += group.class2_pixels

        balance_cost = float(
            PATCH_RATIO_WEIGHT
            * np.sum(((p - ideal_patches) / (ideal_patches + 1.0)) ** 2)
            + CLASS1_BALANCE_WEIGHT
            * np.sum(((c1 - ideal_class1) / (ideal_class1 + 1.0)) ** 2)
            + CLASS2_BALANCE_WEIGHT
            * np.sum(((c2 - ideal_class2) / (ideal_class2 + 1.0)) ** 2)
        )
        choices.append((spatial_cost + balance_cost, split_id))

    return min(choices)[1]


def assign_records_with_component_lock(
    records: Sequence[PatchRecord],
    component_assignment: Mapping[int, int],
    candidate: LShapeCandidate,
) -> Tuple[List[PatchRecord], List[PatchRecord], Dict[str, object]]:
    """Assign patches without merging independent components.

    Background-only patches follow the preferred L/T shape. A positive patch
    follows its component assignment when all contained components agree. If
    contained components belong to different splits, the patch is omitted.
    """
    kept: List[PatchRecord] = []
    dropped: List[PatchRecord] = []

    positive_patch_count = 0
    preferred_mismatch_count = 0

    for record in records:
        preferred_split_id = l_shape_preferred_split_id(
            record,
            candidate,
        )

        if not record.component_ids:
            split_id = preferred_split_id
        else:
            positive_patch_count += 1
            assigned_splits = {
                int(component_assignment[component_id])
                for component_id in record.component_ids
            }

            if len(assigned_splits) != 1:
                dropped.append(
                    replace(record, split="dropped_conflict")
                )
                continue

            split_id = next(iter(assigned_splits))
            if split_id != preferred_split_id:
                preferred_mismatch_count += 1

        kept.append(
            replace(record, split=ID_TO_SPLIT[split_id])
        )

    kept = sorted(
        kept,
        key=lambda record: (
            SPLIT_TO_ID[record.split],
            record.y,
            record.x,
        ),
    )
    dropped = sorted(dropped, key=lambda record: (record.y, record.x))

    metadata = {
        "candidate_patch_count": len(records),
        "kept_patch_count": len(kept),
        "dropped_conflicting_patch_count": len(dropped),
        "dropped_conflicting_patch_fraction": (
            len(dropped) / len(records) if records else 0.0
        ),
        "positive_patch_count": positive_patch_count,
        "positive_preferred_mismatch_count": preferred_mismatch_count,
        "positive_preferred_mismatch_fraction": (
            preferred_mismatch_count / positive_patch_count
            if positive_patch_count
            else 0.0
        ),
    }
    return kept, dropped, metadata


def final_record_assignment_score(
    records: Sequence[PatchRecord],
    assignment_metadata: Mapping[str, object],
) -> Optional[float]:
    """Evaluate the actual kept patches after conflict removal."""
    by_split = {
        split_name: [
            record for record in records
            if record.split == split_name
        ]
        for split_name in SPLIT_ORDER
    }

    patch_totals = np.asarray(
        [len(by_split[name]) for name in SPLIT_ORDER],
        dtype=np.float64,
    )
    class1_totals = np.asarray(
        [
            sum(record.class1_pixels for record in by_split[name])
            for name in SPLIT_ORDER
        ],
        dtype=np.float64,
    )
    class2_totals = np.asarray(
        [
            sum(record.class2_pixels for record in by_split[name])
            for name in SPLIT_ORDER
        ],
        dtype=np.float64,
    )

    for split_id, split_name in enumerate(SPLIT_ORDER):
        split_records = by_split[split_name]
        if len(split_records) < MIN_PATCHES_PER_SPLIT[split_name]:
            return None

        class1_patch_count = sum(
            record.class1_pixels > 0
            for record in split_records
        )
        if (
            class1_patch_count
            < MIN_CLASS1_PATCHES_PER_SPLIT[split_name]
        ):
            return None
        if (
            class1_totals[split_id]
            < MIN_CLASS1_PIXELS_PER_SPLIT[split_name]
        ):
            return None

        if (
            REQUIRE_CLASS2_IN_EVERY_SPLIT
            and class2_totals[split_id]
            < MIN_CLASS2_PIXELS_PER_SPLIT
        ):
            return None

    conflict_fraction = float(
        assignment_metadata["dropped_conflicting_patch_fraction"]
    )
    mismatch_fraction = float(
        assignment_metadata["positive_preferred_mismatch_fraction"]
    )

    if (
        MAX_CONFLICT_PATCH_DROP_FRACTION is not None
        and conflict_fraction > MAX_CONFLICT_PATCH_DROP_FRACTION
    ):
        return None

    if (
        MAX_POSITIVE_PATCH_MISMATCH_FRACTION is not None
        and mismatch_fraction
        > MAX_POSITIVE_PATCH_MISMATCH_FRACTION
    ):
        return None

    return float(
        PATCH_RATIO_WEIGHT * relative_target_error(patch_totals)
        + CLASS1_BALANCE_WEIGHT
        * relative_target_error(class1_totals)
        + CLASS2_BALANCE_WEIGHT
        * relative_target_error(class2_totals)
        + SPATIAL_SHAPE_WEIGHT * mismatch_fraction
        + CONFLICT_PATCH_DROP_WEIGHT * conflict_fraction
    )



def complete_component_assignment(
    initial_assignment: Mapping[int, int],
    candidate: LShapeCandidate,
    records: Sequence[PatchRecord],
    groups: Sequence[LockedGroup],
    preferred_counts: Mapping[int, np.ndarray],
    rng: np.random.Generator,
) -> Dict[int, int]:
    """Assign remaining non-seeded components greedily.

    Class-1 components may already be fixed by exhaustive enumeration. The
    remaining components, including class-2-only components, are assigned while
    balancing the requested split ratios and preserving the L/T shape.
    """
    assignment = {
        int(group_id): int(split_id)
        for group_id, split_id in initial_assignment.items()
    }

    total_group_patches = sum(group.patch_count for group in groups)
    total_class1 = sum(group.class1_pixels for group in groups)
    total_class2 = sum(group.class2_pixels for group in groups)
    target = np.asarray(
        [TARGET_SPLIT_RATIOS[name] for name in SPLIT_ORDER],
        dtype=np.float64,
    )
    ideal_patches = total_group_patches * target
    ideal_class1 = total_class1 * target
    ideal_class2 = total_class2 * target

    patch_totals = np.zeros(3, dtype=np.float64)
    class1_totals = np.zeros(3, dtype=np.float64)
    class2_totals = np.zeros(3, dtype=np.float64)
    group_by_id = {group.group_id: group for group in groups}

    for group_id, split_id in assignment.items():
        group = group_by_id[group_id]
        patch_totals[split_id] += group.patch_count
        class1_totals[split_id] += group.class1_pixels
        class2_totals[split_id] += group.class2_pixels

    remaining = [
        group
        for group in groups
        if group.group_id not in assignment
    ]

    # Deterministic ordering is preferable during exhaustive search. The small
    # random tie breaker remains seed-controlled and does not change validity.
    tie_breakers = {
        group.group_id: float(rng.random())
        for group in remaining
    }
    remaining.sort(
        key=lambda group: (
            -(
                group.class1_pixels * 5
                + group.class2_pixels * 2
                + group.patch_count * PATCH_SIZE
            ),
            tie_breakers[group.group_id],
            group.group_id,
        )
    )

    for group in remaining:
        split_id = choose_group_split_with_shape_preference(
            group=group,
            preferred_counts=preferred_counts[group.group_id],
            patch_totals=patch_totals,
            class1_totals=class1_totals,
            class2_totals=class2_totals,
            ideal_patches=ideal_patches,
            ideal_class1=ideal_class1,
            ideal_class2=ideal_class2,
        )
        assignment[group.group_id] = split_id
        patch_totals[split_id] += group.patch_count
        class1_totals[split_id] += group.class1_pixels
        class2_totals[split_id] += group.class2_pixels

    return assignment


def evaluate_complete_component_assignment(
    assignment: Mapping[int, int],
    candidate: LShapeCandidate,
    records: Sequence[PatchRecord],
) -> Optional[Tuple[float, Dict[str, object]]]:
    """Score the actual retained patches after conflict-patch removal."""
    preview_records, preview_dropped, patch_assignment_metadata = (
        assign_records_with_component_lock(
            records,
            assignment,
            candidate,
        )
    )

    score = final_record_assignment_score(
        preview_records,
        patch_assignment_metadata,
    )
    if score is None:
        return None

    metadata = {
        "candidate": candidate.description,
        **patch_assignment_metadata,
        "preview_dropped_stems": [
            record.stem for record in preview_dropped[:20]
        ],
        "final_score": float(score),
    }
    return float(score), metadata


def evaluate_l_shape_assignment_exhaustive(
    candidate: LShapeCandidate,
    records: Sequence[PatchRecord],
    groups: Sequence[LockedGroup],
    rng: np.random.Generator,
) -> Optional[Tuple[float, Dict[int, int], Dict[str, object]]]:
    """Exhaustively assign every class-1 component to the three splits."""
    from itertools import product

    class1_groups = sorted(
        (group for group in groups if group.class1_patch_count > 0),
        key=lambda group: group.group_id,
    )
    class1_count = len(class1_groups)

    if class1_count < 3:
        return None
    if class1_count > EXACT_CLASS1_ASSIGNMENT_MAX_GROUPS:
        raise ValueError(
            "Exhaustive evaluator called with too many class-1 groups: "
            f"{class1_count}"
        )

    preferred_counts = {
        group.group_id: preferred_split_counts_for_group(
            group,
            records,
            candidate,
        )
        for group in groups
    }

    best_score = float("inf")
    best_assignment: Optional[Dict[int, int]] = None
    best_metadata: Optional[Dict[str, object]] = None

    tested_assignment_count = 0
    surjective_assignment_count = 0
    valid_assignment_count = 0

    class1_group_ids = [group.group_id for group in class1_groups]

    for split_tuple in product(range(3), repeat=class1_count):
        tested_assignment_count += 1

        # Every split must receive at least one independent class-1 component.
        if set(split_tuple) != {0, 1, 2}:
            continue
        surjective_assignment_count += 1

        initial_assignment = {
            group_id: int(split_id)
            for group_id, split_id in zip(
                class1_group_ids,
                split_tuple,
            )
        }

        complete_assignment = complete_component_assignment(
            initial_assignment=initial_assignment,
            candidate=candidate,
            records=records,
            groups=groups,
            preferred_counts=preferred_counts,
            rng=rng,
        )

        evaluated = evaluate_complete_component_assignment(
            assignment=complete_assignment,
            candidate=candidate,
            records=records,
        )
        if evaluated is None:
            continue

        valid_assignment_count += 1
        score, metadata = evaluated

        if score < best_score:
            best_score = score
            best_assignment = complete_assignment
            best_metadata = metadata

    if best_assignment is None or best_metadata is None:
        return None

    best_metadata = {
        **best_metadata,
        "assignment_search_mode": "exhaustive_class1",
        "class1_group_count": class1_count,
        "theoretical_assignment_count": int(3 ** class1_count),
        "tested_assignment_count": tested_assignment_count,
        "surjective_assignment_count": surjective_assignment_count,
        "valid_assignment_count": valid_assignment_count,
    }
    return float(best_score), best_assignment, best_metadata


def evaluate_l_shape_assignment_greedy(
    candidate: LShapeCandidate,
    records: Sequence[PatchRecord],
    groups: Sequence[LockedGroup],
    rng: np.random.Generator,
) -> Optional[Tuple[float, Dict[int, int], Dict[str, object]]]:
    """Use the previous greedy seed method for larger class-1 sets."""
    class1_groups = [group for group in groups if group.class1_patch_count > 0]
    if len(class1_groups) < 3:
        return None

    preferred_counts = {
        group.group_id: preferred_split_counts_for_group(
            group,
            records,
            candidate,
        )
        for group in groups
    }

    initial_assignment: Dict[int, int] = {}
    unused_ids = {group.group_id for group in class1_groups}

    split_order = sorted(
        range(3),
        key=lambda split_id: sum(
            int(preferred_counts[group.group_id][split_id])
            for group in class1_groups
        ),
    )

    for split_id in split_order:
        choices = [
            group
            for group in class1_groups
            if group.group_id in unused_ids
        ]
        if not choices:
            return None

        selected = max(
            choices,
            key=lambda group: (
                int(preferred_counts[group.group_id][split_id]),
                -(
                    group.patch_count
                    - int(preferred_counts[group.group_id][split_id])
                ),
                group.class1_pixels,
            ),
        )
        initial_assignment[selected.group_id] = split_id
        unused_ids.remove(selected.group_id)

    assignment = complete_component_assignment(
        initial_assignment=initial_assignment,
        candidate=candidate,
        records=records,
        groups=groups,
        preferred_counts=preferred_counts,
        rng=rng,
    )

    evaluated = evaluate_complete_component_assignment(
        assignment=assignment,
        candidate=candidate,
        records=records,
    )
    if evaluated is None:
        return None

    score, metadata = evaluated
    metadata = {
        **metadata,
        "assignment_search_mode": "greedy_class1",
        "class1_group_count": len(class1_groups),
    }
    return float(score), assignment, metadata


def evaluate_l_shape_assignment(
    candidate: LShapeCandidate,
    records: Sequence[PatchRecord],
    groups: Sequence[LockedGroup],
    rng: np.random.Generator,
) -> Optional[Tuple[float, Dict[int, int], Dict[str, object]]]:
    """Choose exhaustive or greedy class-1 assignment search."""
    class1_group_count = sum(
        group.class1_patch_count > 0
        for group in groups
    )

    if class1_group_count <= EXACT_CLASS1_ASSIGNMENT_MAX_GROUPS:
        return evaluate_l_shape_assignment_exhaustive(
            candidate=candidate,
            records=records,
            groups=groups,
            rng=rng,
        )

    return evaluate_l_shape_assignment_greedy(
        candidate=candidate,
        records=records,
        groups=groups,
        rng=rng,
    )


def assign_locked_groups_l_shape_preferred(
    records: Sequence[PatchRecord],
    groups: Sequence[LockedGroup],
) -> Tuple[
    Dict[int, int],
    LShapeCandidate,
    Dict[str, object],
]:
    candidates = generate_l_shape_candidates(records)
    if not candidates:
        raise RuntimeError("No L/T-shaped split candidate was generated.")

    rng = np.random.default_rng(RANDOM_SEED)
    best_score = float("inf")
    best_assignment: Optional[Dict[int, int]] = None
    best_candidate: Optional[LShapeCandidate] = None
    best_metadata: Optional[Dict[str, object]] = None
    accepted = 0

    for candidate in candidates:
        result = evaluate_l_shape_assignment(
            candidate,
            records,
            groups,
            rng,
        )
        if result is None:
            continue
        accepted += 1
        score, assignment, metadata = result
        if score < best_score:
            best_score = score
            best_assignment = assignment
            best_candidate = candidate
            best_metadata = metadata

    if (
        best_assignment is None
        or best_candidate is None
        or best_metadata is None
    ):
        raise RuntimeError(
            "No L/T-shaped preferred assignment can keep whole connected "
            "filaments in one split while retaining class 1 in train, val "
            "and test."
        )

    best_metadata = {
        **best_metadata,
        "candidate_count": len(candidates),
        "accepted_candidate_count": accepted,
        "best_score": float(best_score),
        "group_assignment": {
            str(group_id): ID_TO_SPLIT[split_id]
            for group_id, split_id in sorted(best_assignment.items())
        },
    }
    return best_assignment, best_candidate, best_metadata


def assign_all_patches_from_l_shape(
    records: Sequence[PatchRecord],
    group_assignment: Mapping[int, int],
    candidate: LShapeCandidate,
) -> Tuple[List[PatchRecord], List[PatchRecord], Dict[str, object]]:
    return assign_records_with_component_lock(
        records,
        group_assignment,
        candidate,
    )



# =============================================================================
# GROUP ASSIGNMENT
# =============================================================================

def relative_target_error(values: np.ndarray) -> float:
    target = np.asarray(
        [TARGET_SPLIT_RATIOS[name] for name in SPLIT_ORDER],
        dtype=np.float64,
    )
    total = float(values.sum())
    if total <= 0:
        return 0.0
    expected = total * target
    return float(np.sum(((values - expected) / (expected + 1.0)) ** 2))


def assignment_score(
    assignment: Mapping[int, int],
    groups: Sequence[LockedGroup],
) -> float:
    patch_totals = np.zeros(3, dtype=np.float64)
    class1_totals = np.zeros(3, dtype=np.float64)
    class2_totals = np.zeros(3, dtype=np.float64)

    for group in groups:
        split_id = assignment[group.group_id]
        patch_totals[split_id] += group.patch_count
        class1_totals[split_id] += group.class1_pixels
        class2_totals[split_id] += group.class2_pixels

    score = (
        3.0 * relative_target_error(patch_totals)
        + 5.0 * relative_target_error(class1_totals)
        + 2.0 * relative_target_error(class2_totals)
    )

    for split_id, split_name in enumerate(SPLIT_ORDER):
        class1_patch_count = sum(
            group.class1_patch_count
            for group in groups
            if assignment[group.group_id] == split_id
        )
        if class1_patch_count < MIN_CLASS1_PATCHES_PER_SPLIT[split_name]:
            score += 1e12
        if class1_totals[split_id] < MIN_CLASS1_PIXELS_PER_SPLIT[split_name]:
            score += 1e12

    return score


def assign_locked_groups(
    groups: Sequence[LockedGroup],
) -> Tuple[Dict[int, int], Dict[str, object]]:
    class1_groups = [group for group in groups if group.class1_patch_count > 0]
    if len(class1_groups) < 3:
        raise RuntimeError(
            "Fewer than three independent class-1 locked groups exist. "
            "Class 1 cannot be placed in train, val and test without splitting "
            "a filament."
        )

    rng = np.random.default_rng(RANDOM_SEED)
    group_by_id = {group.group_id: group for group in groups}
    group_ids = [group.group_id for group in groups]
    class1_ids = [group.group_id for group in class1_groups]

    total_patches = sum(group.patch_count for group in groups)
    total_class1 = sum(group.class1_pixels for group in groups)
    total_class2 = sum(group.class2_pixels for group in groups)
    target = np.asarray(
        [TARGET_SPLIT_RATIOS[name] for name in SPLIT_ORDER],
        dtype=np.float64,
    )
    ideal_patches = total_patches * target
    ideal_class1 = total_class1 * target
    ideal_class2 = total_class2 * target

    strongest = sorted(
        class1_groups,
        key=lambda g: (g.class1_pixels, g.class1_patch_count, g.patch_count),
        reverse=True,
    )
    seed_triplets: List[Tuple[int, int, int]] = [
        (strongest[0].group_id, strongest[1].group_id, strongest[2].group_id)
    ]
    for _ in range(ASSIGNMENT_TRIALS):
        selected = rng.choice(class1_ids, size=3, replace=False)
        seed_triplets.append(tuple(int(v) for v in selected))

    best_assignment: Optional[Dict[int, int]] = None
    best_score = float("inf")

    for seed_triplet in seed_triplets:
        assignment: Dict[int, int] = {}
        split_permutation = rng.permutation(3)
        patch_totals = np.zeros(3, dtype=np.float64)
        class1_totals = np.zeros(3, dtype=np.float64)
        class2_totals = np.zeros(3, dtype=np.float64)

        for position, group_id in enumerate(seed_triplet):
            split_id = int(split_permutation[position])
            group = group_by_id[group_id]
            assignment[group_id] = split_id
            patch_totals[split_id] += group.patch_count
            class1_totals[split_id] += group.class1_pixels
            class2_totals[split_id] += group.class2_pixels

        remaining = [
            group_by_id[group_id]
            for group_id in group_ids
            if group_id not in assignment
        ]
        rng.shuffle(remaining)
        remaining.sort(
            key=lambda g: (
                g.class1_pixels * 5 + g.class2_pixels * 2 + g.patch_count * 256
            ),
            reverse=True,
        )

        for group in remaining:
            choices: List[Tuple[float, int]] = []
            for split_id in range(3):
                p = patch_totals.copy()
                c1 = class1_totals.copy()
                c2 = class2_totals.copy()
                p[split_id] += group.patch_count
                c1[split_id] += group.class1_pixels
                c2[split_id] += group.class2_pixels
                score = float(
                    np.sum(((p - ideal_patches) / (ideal_patches + 1.0)) ** 2)
                    + 2.0 * np.sum(
                        ((c1 - ideal_class1) / (ideal_class1 + 1.0)) ** 2
                    )
                    + np.sum(((c2 - ideal_class2) / (ideal_class2 + 1.0)) ** 2)
                )
                choices.append((score, split_id))

            _, selected_split = min(choices)
            assignment[group.group_id] = selected_split
            patch_totals[selected_split] += group.patch_count
            class1_totals[selected_split] += group.class1_pixels
            class2_totals[selected_split] += group.class2_pixels

        score = assignment_score(assignment, groups)
        if score < best_score:
            best_score = score
            best_assignment = dict(assignment)

    if best_assignment is None or best_score >= 1e12:
        raise RuntimeError(
            "Could not assign whole filament groups to all three splits while "
            "retaining class 1 in each split."
        )

    metadata = {
        "trials": len(seed_triplets),
        "best_score": float(best_score),
        "group_assignment": {
            str(group_id): ID_TO_SPLIT[split_id]
            for group_id, split_id in sorted(best_assignment.items())
        },
    }
    return best_assignment, metadata


def assign_background_patches(
    records: Sequence[PatchRecord],
    group_assignment: Mapping[int, int],
) -> List[PatchRecord]:
    assigned: List[PatchRecord] = []
    backgrounds: List[PatchRecord] = []
    counts = np.zeros(3, dtype=np.int64)

    for record in records:
        if record.locked_group_id == 0:
            backgrounds.append(record)
            continue
        split_id = int(group_assignment[record.locked_group_id])
        assigned.append(replace(record, split=ID_TO_SPLIT[split_id]))
        counts[split_id] += 1

    total = len(records)
    target_counts = np.asarray(
        [TARGET_SPLIT_RATIOS[name] * total for name in SPLIT_ORDER],
        dtype=np.float64,
    )

    for record in sorted(backgrounds, key=lambda r: (r.y, r.x)):
        deficits = target_counts - counts
        if np.max(deficits) > 0:
            split_id = int(np.argmax(deficits))
        else:
            ratios = counts / max(int(counts.sum()), 1)
            split_id = int(np.argmax(target_counts / total - ratios))
        assigned.append(replace(record, split=ID_TO_SPLIT[split_id]))
        counts[split_id] += 1

    return sorted(assigned, key=lambda r: (SPLIT_TO_ID[r.split], r.y, r.x))


# =============================================================================
# VALIDATION
# =============================================================================

def validate_final_assignment(
    records: Sequence[PatchRecord],
    component_count: int,
) -> Dict[str, object]:
    by_split = {
        name: [record for record in records if record.split == name]
        for name in SPLIT_ORDER
    }

    for split_name in SPLIT_ORDER:
        split_records = by_split[split_name]
        if len(split_records) < MIN_PATCHES_PER_SPLIT[split_name]:
            raise RuntimeError(
                f"{split_name} has {len(split_records)} patches; minimum is "
                f"{MIN_PATCHES_PER_SPLIT[split_name]}."
            )
        class1_patch_count = sum(r.class1_pixels > 0 for r in split_records)
        class1_pixels = sum(r.class1_pixels for r in split_records)
        class2_pixels = sum(r.class2_pixels for r in split_records)
        if class1_patch_count < MIN_CLASS1_PATCHES_PER_SPLIT[split_name]:
            raise RuntimeError(f"{split_name} lacks a class-1-positive patch.")
        if class1_pixels < MIN_CLASS1_PIXELS_PER_SPLIT[split_name]:
            raise RuntimeError(
                f"{split_name} has only {class1_pixels} class-1 pixels."
            )
        if (
            REQUIRE_CLASS2_IN_EVERY_SPLIT
            and class2_pixels < MIN_CLASS2_PIXELS_PER_SPLIT
        ):
            raise RuntimeError(f"{split_name} lacks class-2 pixels.")

    component_splits: Dict[int, Set[str]] = {
        component_id: set()
        for component_id in range(1, component_count + 1)
    }
    for record in records:
        for component_id in record.component_ids:
            component_splits[component_id].add(record.split)

    leaking = {
        component_id: sorted(splits)
        for component_id, splits in component_splits.items()
        if len(splits) > 1
    }
    if leaking:
        raise AssertionError(
            "Connected filament leakage detected: "
            f"{list(leaking.items())[:10]}"
        )

    coordinates = {
        name: {(record.y, record.x) for record in by_split[name]}
        for name in SPLIT_ORDER
    }
    for left_index, left in enumerate(SPLIT_ORDER):
        for right in SPLIT_ORDER[left_index + 1:]:
            if coordinates[left] & coordinates[right]:
                raise AssertionError(f"Patch overlap between {left} and {right}")

    return {
        "strict_component_leakage_count": 0,
        "split_patch_counts": {name: len(by_split[name]) for name in SPLIT_ORDER},
        "split_class1_patch_counts": {
            name: sum(r.class1_pixels > 0 for r in by_split[name])
            for name in SPLIT_ORDER
        },
        "split_class1_pixels": {
            name: sum(r.class1_pixels for r in by_split[name])
            for name in SPLIT_ORDER
        },
        "split_class2_pixels": {
            name: sum(r.class2_pixels for r in by_split[name])
            for name in SPLIT_ORDER
        },
    }


# =============================================================================
# OUTPUT
# =============================================================================

def prepare_output_dirs(root: Path, region: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    for split_name in SPLIT_ORDER:
        for modality in ("in", "out"):
            directory = root / split_name / modality / region
            if CLEAN_REGION_OUTPUT and directory.exists():
                shutil.rmtree(directory)
            directory.mkdir(parents=True, exist_ok=True)

    metadata_dir = root / "metadata" / region
    if CLEAN_REGION_OUTPUT and metadata_dir.exists():
        shutil.rmtree(metadata_dir)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    return metadata_dir


def make_rgb_image(
    col_density: np.ndarray,
    records: Sequence[PatchRecord],
) -> Tuple[np.ndarray, Dict[str, object]]:
    positive = np.isfinite(col_density) & (col_density > 0)
    log_map = np.full(col_density.shape, np.nan, dtype=np.float64)
    log_map[positive] = np.log10(col_density[positive])

    train_values: List[np.ndarray] = []
    for record in records:
        if record.split != "train":
            continue
        patch = log_map[
            record.y:record.y + PATCH_SIZE,
            record.x:record.x + PATCH_SIZE,
        ]
        values = patch[np.isfinite(patch)]
        if values.size:
            train_values.append(values)

    if not train_values:
        raise RuntimeError("No train pixels are available for normalization.")

    values = np.concatenate(train_values)
    low, high = np.percentile(
        values,
        [NORMALIZE_LOW_PERCENTILE, NORMALIZE_HIGH_PERCENTILE],
    )
    if high <= low:
        high = low + 1e-6

    scaled = np.clip((log_map - low) / (high - low), 0.0, 1.0)
    scaled[~positive] = 0.0
    gray = np.rint(scaled * 255.0).astype(np.uint8)
    rgb = np.repeat(gray[..., None], 3, axis=2)

    return rgb, {
        "transform": "log10",
        "source": "train patches only",
        "low_percentile": NORMALIZE_LOW_PERCENTILE,
        "high_percentile": NORMALIZE_HIGH_PERCENTILE,
        "log10_low": float(low),
        "log10_high": float(high),
    }


def assign_storage_subdirectories(
    records: Sequence[PatchRecord],
) -> List[PatchRecord]:
    output: List[PatchRecord] = []
    for split_name in SPLIT_ORDER:
        split_records = sorted(
            (record for record in records if record.split == split_name),
            key=lambda record: (record.y, record.x),
        )
        for index, record in enumerate(split_records):
            output.append(
                replace(
                    record,
                    storage_subdir=f"{index // FILES_PER_SUBDIR:04d}",
                )
            )
    return output


def save_patch_pairs(
    root: Path,
    rgb_image: np.ndarray,
    label_map: np.ndarray,
    valid_mask: np.ndarray,
    records: Iterable[PatchRecord],
) -> None:
    for record in records:
        y1, x1 = record.y, record.x
        y2, x2 = y1 + PATCH_SIZE, x1 + PATCH_SIZE
        image_patch = rgb_image[y1:y2, x1:x2]
        mask_patch = label_map[y1:y2, x1:x2].copy()
        mask_patch[~valid_mask[y1:y2, x1:x2]] = IGNORE_LABEL

        image_dir = (
            root / record.split / "in" / record.region / record.storage_subdir
        )
        mask_dir = (
            root / record.split / "out" / record.region / record.storage_subdir
        )
        image_dir.mkdir(parents=True, exist_ok=True)
        mask_dir.mkdir(parents=True, exist_ok=True)

        Image.fromarray(image_patch).save(
            image_dir / f"{record.stem}_in.png",
            compress_level=1,
        )
        Image.fromarray(mask_patch.astype(np.uint8)).save(
            mask_dir / f"{record.stem}_out.png",
            compress_level=1,
        )


def save_manifest(path: Path, records: Sequence[PatchRecord]) -> None:
    fieldnames = [
        "region", "split", "stem", "y", "x", "valid_ratio",
        "background_pixels", "class1_pixels", "class2_pixels",
        "ignore_pixels", "component_ids", "locked_group_id",
        "storage_subdir",
    ]
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            row = asdict(record)
            row["component_ids"] = ";".join(
                str(value) for value in record.component_ids
            )
            writer.writerow(row)


def records_stats(records: Sequence[PatchRecord]) -> Dict[str, object]:
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
    return {
        "patch_count": len(records),
        "class1_patch_count": int(sum(r.class1_pixels > 0 for r in records)),
        "class2_patch_count": int(sum(r.class2_pixels > 0 for r in records)),
        "background_only_patch_count": int(sum(r.is_background_only for r in records)),
        "pixel_counts": pixel_counts,
        "valid_pixel_fractions": {
            key: pixel_counts[key] / valid_total if valid_total else 0.0
            for key in ("background", "class1", "class2")
        },
    }


def save_split_map_and_overview(
    metadata_dir: Path,
    col_density: np.ndarray,
    label_map: np.ndarray,
    records: Sequence[PatchRecord],
) -> None:
    patch_split_map = np.full(label_map.shape, -1, dtype=np.int8)
    for record in records:
        patch_split_map[
            record.y:record.y + PATCH_SIZE,
            record.x:record.x + PATCH_SIZE,
        ] = SPLIT_TO_ID[record.split]
    np.save(metadata_dir / "patch_split_map.npy", patch_split_map)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
    except ImportError:
        print("matplotlib unavailable; split_overview.png was not created.")
        return

    positive = np.isfinite(col_density) & (col_density > 0)
    display = np.full(col_density.shape, np.nan, dtype=np.float64)
    display[positive] = np.log10(col_density[positive])

    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    axes[0].imshow(display, origin="lower", cmap="inferno")
    edge_colors = {"train": "tab:blue", "val": "tab:green", "test": "tab:red"}
    for record in records:
        axes[0].add_patch(
            Rectangle(
                (record.x, record.y),
                PATCH_SIZE,
                PATCH_SIZE,
                fill=False,
                linewidth=0.5,
                edgecolor=edge_colors[record.split],
            )
        )
    axes[0].set_title(
        "Patch split: blue=train, green=val, red=test\n"
        "L/T-shape preferred; whole filaments may locally bend the boundary"
    )
    axes[0].axis("off")

    labels_display = label_map.astype(np.float32)
    labels_display[labels_display == IGNORE_LABEL] = np.nan
    axes[1].imshow(labels_display, origin="lower", cmap="viridis", vmin=0, vmax=2)
    axes[1].set_title("Original 3-class label map")
    axes[1].axis("off")

    fig.tight_layout()
    fig.savefig(metadata_dir / "split_overview.png", dpi=160)
    plt.close(fig)


# =============================================================================
# REGION / BATCH PROCESSING
# =============================================================================


@dataclass
class GridSplitResult:
    offset_x: int
    offset_y: int
    candidate_records: List[PatchRecord]
    assigned_records: List[PatchRecord]
    dropped_records: List[PatchRecord]
    groups: List[LockedGroup]
    group_metadata: Dict[str, object]
    group_assignment: Dict[int, int]
    selected_l_shape: LShapeCandidate
    assignment_metadata: Dict[str, object]
    conflict_metadata: Dict[str, object]
    validation_metadata: Dict[str, object]

    @property
    def score(self) -> float:
        return float(self.assignment_metadata["best_score"])


def evaluate_patch_grid(
    region: str,
    label_map: np.ndarray,
    valid_mask: np.ndarray,
    component_map: np.ndarray,
    component_count: int,
    offset_x: int,
    offset_y: int,
) -> GridSplitResult:
    candidate_records = enumerate_candidate_patches(
        region=region,
        label_map=label_map,
        valid_mask=valid_mask,
        component_map=component_map,
        offset_x=offset_x,
        offset_y=offset_y,
    )

    candidate_records, groups, group_metadata = create_component_groups(
        candidate_records,
        component_count,
        label_map,
        component_map,
    )

    visible_class1_groups = int(
        group_metadata["components_with_class1"]
    )
    if visible_class1_groups < 3:
        raise RuntimeError(
            "This grid contains fewer than three usable independent "
            "class-1 components: "
            f"{visible_class1_groups}"
        )

    (
        group_assignment,
        selected_l_shape,
        assignment_metadata,
    ) = assign_locked_groups_l_shape_preferred(
        candidate_records,
        groups,
    )

    (
        assigned_records,
        dropped_records,
        conflict_metadata,
    ) = assign_all_patches_from_l_shape(
        candidate_records,
        group_assignment,
        selected_l_shape,
    )

    validation_metadata = validate_final_assignment(
        assigned_records,
        component_count,
    )

    return GridSplitResult(
        offset_x=offset_x,
        offset_y=offset_y,
        candidate_records=candidate_records,
        assigned_records=assigned_records,
        dropped_records=dropped_records,
        groups=groups,
        group_metadata=group_metadata,
        group_assignment=group_assignment,
        selected_l_shape=selected_l_shape,
        assignment_metadata=assignment_metadata,
        conflict_metadata=conflict_metadata,
        validation_metadata=validation_metadata,
    )


def candidate_grid_offsets() -> List[Tuple[int, int]]:
    return [(0, 0)] + [
        (offset_x, offset_y)
        for offset_y in PATCH_GRID_OFFSETS
        for offset_x in PATCH_GRID_OFFSETS
        if (offset_x, offset_y) != (0, 0)
    ]


def process_region(
    region_name: str,
    column_density_fits: str,
    label_map_path: str,
) -> Dict[str, object]:
    region = normalize_region_name(region_name)
    root = Path(DATASET_ROOT)
    metadata_dir = prepare_output_dirs(root, region)

    print(f"Loading region={region}...")
    col_density, label_map = load_data(column_density_fits, label_map_path)
    invalid_mask = build_invalid_mask(col_density)
    valid_mask = ~invalid_mask

    print("Building 8-neighbour filament connected components...")
    component_map, component_count, component_metadata = (
        build_filament_component_map(label_map, valid_mask)
    )
    print(f"Connected components: {component_count}")

    print(
        "Searching non-overlapping 256x256 patch-grid alignments..."
    )

    successful_grid_results: List[GridSplitResult] = []
    failed_grid_results: List[Dict[str, object]] = []
    offsets = candidate_grid_offsets()

    for grid_index, (offset_x, offset_y) in enumerate(offsets, start=1):
        if (
            SEARCH_GRID_OFFSETS_ON_FAILURE_ONLY
            and successful_grid_results
        ):
            break

        print(
            f"  Grid {grid_index}/{len(offsets)}: "
            f"offset_x={offset_x}, offset_y={offset_y}"
        )

        try:
            result = evaluate_patch_grid(
                region=region,
                label_map=label_map,
                valid_mask=valid_mask,
                component_map=component_map,
                component_count=component_count,
                offset_x=offset_x,
                offset_y=offset_y,
            )
            successful_grid_results.append(result)
            print(
                "    SUCCESS: "
                f"candidate_patches={len(result.candidate_records)}, "
                f"kept={len(result.assigned_records)}, "
                f"dropped={len(result.dropped_records)}, "
                f"class1_groups="
                f"{result.group_metadata['components_with_class1']}, "
                f"score={result.score:.4f}"
            )
        except Exception as exc:
            failed_grid_results.append(
                {
                    "offset_x": offset_x,
                    "offset_y": offset_y,
                    "exception_type": type(exc).__name__,
                    "message": str(exc),
                }
            )
            print(
                f"    rejected: {type(exc).__name__}: {exc}"
            )

    if not successful_grid_results:
        raise RuntimeError(
            "No tested non-overlapping patch-grid alignment produced a "
            "no-filament-leakage L/T-shaped split with class 1 in train, "
            "val and test."
        )

    selected_grid = min(
        successful_grid_results,
        key=lambda result: (
            result.score,
            len(result.dropped_records),
            -len(result.assigned_records),
        ),
    )

    records = selected_grid.assigned_records
    dropped_conflict_records = selected_grid.dropped_records
    groups = selected_grid.groups
    group_metadata = selected_grid.group_metadata
    group_assignment = selected_grid.group_assignment
    selected_l_shape = selected_grid.selected_l_shape
    assignment_metadata = selected_grid.assignment_metadata
    conflict_metadata = selected_grid.conflict_metadata
    validation_metadata = selected_grid.validation_metadata

    print(
        "Selected patch grid: "
        f"offset_x={selected_grid.offset_x}, "
        f"offset_y={selected_grid.offset_y}"
    )
    print(
        f"Independent component groups visible in grid: {len(groups)}, "
        f"class-1 components: "
        f"{group_metadata['components_with_class1']}, "
        f"multi-component patches: "
        f"{group_metadata['multi_component_patch_count']}"
    )
    print(f"Selected layout: {selected_l_shape.description}")
    print(
        "Class-1 assignment search: "
        f"{assignment_metadata.get('assignment_search_mode')}, "
        f"groups={assignment_metadata.get('class1_group_count')}, "
        f"valid assignments="
        f"{assignment_metadata.get('valid_assignment_count', 'n/a')}"
    )

    records = assign_storage_subdirectories(records)
    rgb_image, normalization_metadata = make_rgb_image(col_density, records)

    print("Saving PNG patch/mask pairs...")
    save_patch_pairs(root, rgb_image, label_map, valid_mask, records)

    np.save(metadata_dir / "filament_component_map.npy", component_map)
    save_manifest(metadata_dir / "patch_manifest.csv", records)
    save_manifest(
        metadata_dir / "dropped_conflicting_patches.csv",
        dropped_conflict_records,
    )
    save_split_map_and_overview(
        metadata_dir, col_density, label_map, records
    )

    patch_stats = {
        split_name: records_stats(
            [record for record in records if record.split == split_name]
        )
        for split_name in SPLIT_ORDER
    }

    summary = {
        "input": {
            "region": region,
            "column_density_fits": column_density_fits,
            "label_map_path": label_map_path,
            "shape": list(label_map.shape),
        },
        "method": {
            "name": (
                "L/T-shape preferred component-locked "
                "non-overlapping patch split"
            ),
            "patch_size": PATCH_SIZE,
            "patch_stride": PATCH_STRIDE,
            "target_split_ratios": dict(TARGET_SPLIT_RATIOS),
            "guarantees": {
                "each_region_contributes_to_all_three_splits": True,
                "same_connected_filament_in_multiple_splits": False,
                "independent_filaments_merged_by_patch_cooccurrence": False,
                "conflicting_multi_component_patches_are_saved": False,
                "input_pixel_overlap_between_splits": False,
                "class1_required_in_each_split": True,
            },
        },
        "patch_grid_search": {
            "selected_offset_x": selected_grid.offset_x,
            "selected_offset_y": selected_grid.offset_y,
            "patch_size": PATCH_SIZE,
            "patch_stride": PATCH_STRIDE,
            "configured_offsets": list(PATCH_GRID_OFFSETS),
            "failed_offsets": failed_grid_results,
            "success_count": len(successful_grid_results),
        },
        "component_metadata": component_metadata,
        "component_group_metadata": group_metadata,
        "conflicting_patch_handling": {
            **conflict_metadata,
            "policy": (
                "drop only a patch whose independent components "
                "are assigned to different splits"
            ),
        },
        "selected_spatial_layout": {
            "family": selected_l_shape.family,
            "primary_cut": selected_l_shape.primary_cut,
            "secondary_cut": selected_l_shape.secondary_cut,
            "region_to_split": [
                ID_TO_SPLIT[value]
                for value in selected_l_shape.region_to_split
            ],
            "description": selected_l_shape.description,
        },
        "assignment_metadata": assignment_metadata,
        "validation": validation_metadata,
        "normalization": normalization_metadata,
        "patch_stats": patch_stats,
    }
    with (metadata_dir / "dataset_summary.json").open(
        "w", encoding="utf-8"
    ) as file:
        json.dump(summary, file, ensure_ascii=False, indent=2)

    print("\nCompleted.")
    print(f"Region: {region}")
    print(f"Metadata: {metadata_dir}")
    print(
        "Patch-grid offset: "
        f"x={selected_grid.offset_x}, y={selected_grid.offset_y}"
    )
    print(
        "Dropped conflicting multi-component patches: "
        f"{len(dropped_conflict_records)}"
    )
    for split_name in SPLIT_ORDER:
        stats = patch_stats[split_name]
        fractions = stats["valid_pixel_fractions"]
        print(
            f"{split_name}: patches={stats['patch_count']}, "
            f"class1_patches={stats['class1_patch_count']}, "
            f"class fractions 0={fractions['background']:.6f}, "
            f"1={fractions['class1']:.6f}, "
            f"2={fractions['class2']:.6f}"
        )
    return summary


def validate_region_configuration() -> None:
    if not REGION_DATASETS:
        raise ValueError("REGION_DATASETS is empty.")
    for region, config in REGION_DATASETS.items():
        normalize_region_name(region)
        for key in ("column_density_fits", "label_map_path"):
            value = config.get(key)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(
                    f"REGION_DATASETS[{region!r}][{key!r}] is invalid."
                )


def run_all_regions() -> Dict[str, object]:
    validate_region_configuration()
    if RUN_ALL_CONFIGURED_REGIONS:
        regions = [
            region
            for region, config in REGION_DATASETS.items()
            if bool(config.get("enabled", True))
        ]
    else:
        regions = [normalize_region_name(SINGLE_REGION_NAME)]

    print("=" * 78)
    print("Regions to process:", ", ".join(regions))
    print("Method: L/T-shape preferred, independent-component locked split")
    print(f"Dataset root: {DATASET_ROOT}")
    print("=" * 78)

    completed: Dict[str, object] = {}
    failed: Dict[str, object] = {}

    for index, region in enumerate(regions, start=1):
        print()
        print("#" * 78)
        print(f"[{index}/{len(regions)}] Processing: {region}")
        print("#" * 78)
        config = REGION_DATASETS[region]
        try:
            completed[region] = process_region(
                region_name=region,
                column_density_fits=str(config["column_density_fits"]),
                label_map_path=str(config["label_map_path"]),
            )
        except Exception as exc:
            failed[region] = {
                "exception_type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            }
            print()
            print("!" * 78)
            print(f"FAILED: {region}")
            print(f"{type(exc).__name__}: {exc}")
            print("!" * 78)
            print(failed[region]["traceback"])
            if not CONTINUE_ON_REGION_ERROR:
                raise

    batch_summary = {
        "method": "L/T-shape preferred component-locked split with exact small-class1 assignment and patch-grid offset search",
        "dataset_root": DATASET_ROOT,
        "completed_regions": list(completed),
        "failed_regions": failed,
        "success_count": len(completed),
        "failure_count": len(failed),
    }
    metadata_root = Path(DATASET_ROOT) / "metadata"
    metadata_root.mkdir(parents=True, exist_ok=True)
    summary_path = metadata_root / "all_regions_batch_summary.json"
    with summary_path.open("w", encoding="utf-8") as file:
        json.dump(batch_summary, file, ensure_ascii=False, indent=2)

    print()
    print("=" * 78)
    print("ALL REGIONS FINISHED")
    print(f"Completed: {len(completed)} -> {', '.join(completed) or 'none'}")
    print(f"Failed:    {len(failed)} -> {', '.join(failed) or 'none'}")
    print(f"Summary:   {summary_path}")
    print("=" * 78)
    return batch_summary


if __name__ == "__main__":
    run_all_regions()


