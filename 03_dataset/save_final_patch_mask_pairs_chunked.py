#!/usr/bin/env python3
from __future__ import annotations

"""
Save image/mask patch pairs from the FINAL saved patch_manifest.csv files.

IMPORTANT
---------
This script DOES NOT split or reassign anything.

It reads the already-finalized:
    <MANIFEST_ROOT>/metadata/<region>/patch_manifest.csv

and writes PNG pairs exactly according to the saved split/y/x assignments.

Output structure:
    <OUTPUT_ROOT>/
      train/in/<region>/<chunk>/yXXXXX_xXXXXX_in.png
      train/out/<region>/<chunk>/yXXXXX_xXXXXX_out.png
      val/in/...
      val/out/...
      test/in/...
      test/out/...

The source 3-class mask values are preserved:
    0   = background
    1   = confirmed_filament
    2   = candidate_filament
    255 = ignore

Image preprocessing matches the dataset builder:
    column density
      -> log10
      -> normalize using TRAIN PATCHES ONLY
      -> uint8 grayscale
      -> duplicate to 3 channels

No train/val/test assignment is changed by this script.
"""

import csv
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
from PIL import Image

try:
    from astropy.io import fits
except ImportError as exc:
    raise ImportError(
        "astropy is required. Install with: pip install astropy"
    ) from exc


# =============================================================================
# USER SETTINGS
# =============================================================================

# Final manifests created by the connected/no-leakage split
MANIFEST_ROOT = Path("/mnt/sp4t/fujiwara.r/split_v2")

# Where the PNG image/mask pairs will be saved.
# If you want to rebuild directly inside split_v2, keep this equal to MANIFEST_ROOT.
OUTPUT_ROOT = Path("/mnt/sp4t/fujiwara.r/split_v2")

# Process all five finalized regions.
REGIONS = (
    "aquila",
    "orionb",
    "ic5146",
    "taurus",
    "polaris",
)

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
        "enabled": False,
        "column_density_fits": (
            "/mnt/sp4t/fujiwara.r/HGBS_orionB_hires_column_density_map.fits"
        ),
        "label_map_path": (
            "/mnt/sp4t/fujiwara.r/HGBS_orionB_hires_column_density_map/vis_3class_map/HGBS_orionB_hires_column_density_map_label3class.fits"
        ),
    },
    "ic5146": {
        "enabled": False,
        "column_density_fits": (
            "/mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map.fits"
        ),
        "label_map_path": (
            "/mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map/vis_3class_map/HGBS_ic5146_hires_column_density_map_label3class.fits"
        ),
    },
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
PATCH_SIZE = 256
FILES_PER_SUBDIR = 1000  # max image/mask pairs per chunk directory

BACKGROUND_LABEL = 0
CLASS1_LABEL = 1
CLASS2_LABEL = 2
IGNORE_LABEL = 255
ALLOWED_LABELS = {0, 1, 2, 255}

# Same settings as the split-builder
MARGIN_PX = 32
NORMALIZE_LOW_PERCENTILE = 1.0
NORMALIZE_HIGH_PERCENTILE = 99.0

# Delete ONLY image/mask folders for the region before rebuilding.
# metadata is never deleted by this script.
CLEAN_EXISTING_PATCH_OUTPUT = True

# Extra validation before saving.
VERIFY_NO_COORDINATE_DUPLICATION = True
VERIFY_MASK_VALUES = True


# =============================================================================
# INPUT
# =============================================================================

def read_fits_primary(path: str) -> np.ndarray:
    with fits.open(path, memmap=True) as hdul:
        if not hdul or hdul[0].data is None:
            raise ValueError(
                f"Primary HDU contains no image data: {path}"
            )
        return np.array(
            hdul[0].data,
            copy=True,
        )


def load_column_density(path: str) -> np.ndarray:
    arr = np.squeeze(
        np.asarray(
            read_fits_primary(path),
            dtype=np.float64,
        )
    )

    if arr.ndim != 2:
        raise ValueError(
            f"Column density must be 2D: {path}, shape={arr.shape}"
        )

    return arr


def load_label_map(path: str) -> np.ndarray:
    p = Path(path)
    suffix = p.suffix.lower()

    if suffix == ".npy":
        raw = np.load(p)
    elif suffix in {".fits", ".fit", ".fts"}:
        raw = read_fits_primary(str(p))
    else:
        raise ValueError(
            f"Unsupported label type: {p}"
        )

    arr = np.squeeze(
        np.asarray(raw)
    )

    if arr.ndim != 2:
        raise ValueError(
            f"Label map must be 2D: {p}, shape={arr.shape}"
        )

    if np.issubdtype(arr.dtype, np.floating):
        finite = np.isfinite(arr)
        result = np.full(
            arr.shape,
            IGNORE_LABEL,
            dtype=np.uint8,
        )

        values = arr[finite]
        rounded = np.rint(values)

        if values.size and not np.allclose(
            values,
            rounded,
            atol=1e-6,
        ):
            raise ValueError(
                f"Non-integer label values in {p}"
            )

        if rounded.size:
            numeric = rounded.astype(np.int64)
            unexpected = (
                set(int(v) for v in np.unique(numeric))
                - ALLOWED_LABELS
            )
            if unexpected:
                raise ValueError(
                    f"Unexpected labels {sorted(unexpected)} in {p}"
                )

            result[finite] = numeric.astype(
                np.uint8
            )

        return result

    numeric = arr.astype(
        np.int64,
        copy=False,
    )

    unexpected = (
        set(int(v) for v in np.unique(numeric))
        - ALLOWED_LABELS
    )

    if unexpected:
        raise ValueError(
            f"Unexpected labels {sorted(unexpected)} in {p}"
        )

    return numeric.astype(
        np.uint8
    )


def build_invalid_mask(
    col_density: np.ndarray,
) -> np.ndarray:
    h, w = col_density.shape

    invalid = np.zeros(
        (h, w),
        dtype=bool,
    )

    if MARGIN_PX > 0:
        invalid[:MARGIN_PX, :] = True
        invalid[-MARGIN_PX:, :] = True
        invalid[:, :MARGIN_PX] = True
        invalid[:, -MARGIN_PX:] = True

    no_data = (
        ~np.isfinite(col_density)
        | (col_density <= 0)
    )

    # To match the original builder as closely as possible,
    # use scipy dilation when available.
    if no_data.any() and MARGIN_PX > 0:
        try:
            from scipy.ndimage import binary_dilation

            no_data = binary_dilation(
                no_data,
                structure=np.ones(
                    (3, 3),
                    dtype=bool,
                ),
                iterations=MARGIN_PX,
            )
        except ImportError:
            raise ImportError(
                "scipy is required to reproduce the original invalid-mask logic."
            )

    invalid |= no_data
    return invalid


# =============================================================================
# MANIFEST
# =============================================================================

def read_manifest(
    region: str,
) -> List[Dict[str, object]]:
    path = (
        MANIFEST_ROOT
        / "metadata"
        / region
        / "patch_manifest.csv"
    )

    if not path.is_file():
        raise FileNotFoundError(
            f"Manifest not found: {path}"
        )

    rows: List[Dict[str, object]] = []

    with path.open(
        "r",
        newline="",
        encoding="utf-8",
    ) as f:
        reader = csv.DictReader(f)

        required = {
            "region",
            "split",
            "stem",
            "y",
            "x",
        }

        missing = required - set(
            reader.fieldnames or []
        )

        if missing:
            raise ValueError(
                f"{path} missing columns: {sorted(missing)}"
            )

        for row in reader:
            split = str(
                row["split"]
            ).strip()

            if split not in {
                "train",
                "val",
                "test",
            }:
                raise ValueError(
                    f"Unexpected split={split!r} in {path}"
                )

            rows.append(
                {
                    "region": region,
                    "split": split,
                    "stem": str(
                        row["stem"]
                    ).strip(),
                    "y": int(row["y"]),
                    "x": int(row["x"]),
                    "storage_subdir":
                        str(
                            row.get(
                                "storage_subdir",
                                "",
                            )
                        ).strip(),
                }
            )

    if not rows:
        raise ValueError(
            f"No rows in manifest: {path}"
        )

    return rows


def verify_manifest(
    rows: List[Dict[str, object]],
    image_shape: Tuple[int, int],
) -> None:
    h, w = image_shape

    coordinates = {}

    for row in rows:
        split = str(row["split"])
        y = int(row["y"])
        x = int(row["x"])

        y2 = y + PATCH_SIZE
        x2 = x + PATCH_SIZE

        if (
            y < 0
            or x < 0
            or y2 > h
            or x2 > w
        ):
            raise ValueError(
                "Manifest patch outside source image: "
                f"split={split}, y={y}, x={x}, "
                f"source_shape={image_shape}"
            )

        coord = (y, x)

        if (
            VERIFY_NO_COORDINATE_DUPLICATION
            and coord in coordinates
        ):
            raise ValueError(
                "Same patch coordinate appears multiple times: "
                f"{coord}, splits={coordinates[coord]} and {split}"
            )

        coordinates[coord] = split


# =============================================================================
# IMAGE NORMALIZATION
# =============================================================================

def make_rgb_image(
    col_density: np.ndarray,
    rows: List[Dict[str, object]],
) -> Tuple[np.ndarray, float, float]:

    positive = (
        np.isfinite(col_density)
        & (col_density > 0)
    )

    log_map = np.full(
        col_density.shape,
        np.nan,
        dtype=np.float64,
    )

    log_map[positive] = np.log10(
        col_density[positive]
    )

    train_values = []

    for row in rows:
        if row["split"] != "train":
            continue

        y = int(row["y"])
        x = int(row["x"])

        patch = log_map[
            y:y + PATCH_SIZE,
            x:x + PATCH_SIZE,
        ]

        values = patch[
            np.isfinite(patch)
        ]

        if values.size:
            train_values.append(
                values
            )

    if not train_values:
        raise RuntimeError(
            "No train pixels available for normalization."
        )

    values = np.concatenate(
        train_values
    )

    low, high = np.percentile(
        values,
        [
            NORMALIZE_LOW_PERCENTILE,
            NORMALIZE_HIGH_PERCENTILE,
        ],
    )

    if high <= low:
        high = low + 1e-6

    scaled = np.clip(
        (log_map - low)
        / (high - low),
        0.0,
        1.0,
    )

    scaled[~positive] = 0.0

    gray = np.rint(
        scaled * 255.0
    ).astype(np.uint8)

    rgb = np.repeat(
        gray[..., None],
        3,
        axis=2,
    )

    return (
        rgb,
        float(low),
        float(high),
    )


# =============================================================================
# OUTPUT
# =============================================================================

def clean_region_output(
    region: str,
) -> None:
    if not CLEAN_EXISTING_PATCH_OUTPUT:
        return

    for split in (
        "train",
        "val",
        "test",
    ):
        for modality in (
            "in",
            "out",
        ):
            path = (
                OUTPUT_ROOT
                / split
                / modality
                / region
            )

            if path.exists():
                shutil.rmtree(
                    path
                )


def storage_subdir_for_row(
    row: Dict[str, object],
    split_index: int,
) -> str:
    """
    Always assign a chunk directory from the current split-local index.

    FILES_PER_SUBDIR = 1000 means:
        index    0..999   -> 0000
        index 1000..1999  -> 0001
        index 2000..2999  -> 0002
        ...

    The storage_subdir value saved in an older manifest is intentionally
    ignored, so this export always guarantees at most FILES_PER_SUBDIR files
    per chunk directory.
    """
    return f"{split_index // FILES_PER_SUBDIR:04d}"


def save_region(
    region: str,
) -> None:

    config = REGION_DATASETS[
        region
    ]

    print()
    print("=" * 88)
    print(f"Region: {region}")
    print("=" * 88)

    col_density = load_column_density(
        config["column_density_fits"]
    )

    label_map = load_label_map(
        config["label_map_path"]
    )

    if col_density.shape != label_map.shape:
        raise ValueError(
            "Shape mismatch: "
            f"image={col_density.shape}, "
            f"label={label_map.shape}"
        )

    rows = read_manifest(
        region
    )

    verify_manifest(
        rows,
        col_density.shape,
    )

    invalid_mask = build_invalid_mask(
        col_density
    )

    valid_mask = ~invalid_mask

    rgb_image, log_low, log_high = (
        make_rgb_image(
            col_density,
            rows,
        )
    )

    clean_region_output(
        region
    )

    by_split = {
        split: [
            row
            for row in rows
            if row["split"] == split
        ]
        for split in (
            "train",
            "val",
            "test",
        )
    }

    total_saved = 0

    for split in (
        "train",
        "val",
        "test",
    ):
        split_rows = sorted(
            by_split[split],
            key=lambda r: (
                int(r["y"]),
                int(r["x"]),
            ),
        )

        for index, row in enumerate(
            split_rows
        ):
            y = int(row["y"])
            x = int(row["x"])
            y2 = y + PATCH_SIZE
            x2 = x + PATCH_SIZE

            image_patch = rgb_image[
                y:y2,
                x:x2,
            ]

            mask_patch = label_map[
                y:y2,
                x:x2,
            ].copy()

            mask_patch[
                ~valid_mask[
                    y:y2,
                    x:x2,
                ]
            ] = IGNORE_LABEL

            if VERIFY_MASK_VALUES:
                values = set(
                    int(v)
                    for v in np.unique(
                        mask_patch
                    )
                )

                unexpected = (
                    values
                    - ALLOWED_LABELS
                )

                if unexpected:
                    raise ValueError(
                        f"Unexpected mask labels "
                        f"{sorted(unexpected)} "
                        f"at region={region}, "
                        f"split={split}, "
                        f"y={y}, x={x}"
                    )

            subdir = (
                storage_subdir_for_row(
                    row,
                    index,
                )
            )

            image_dir = (
                OUTPUT_ROOT
                / split
                / "in"
                / region
                / subdir
            )

            mask_dir = (
                OUTPUT_ROOT
                / split
                / "out"
                / region
                / subdir
            )

            image_dir.mkdir(
                parents=True,
                exist_ok=True,
            )

            mask_dir.mkdir(
                parents=True,
                exist_ok=True,
            )

            stem = str(
                row["stem"]
            )

            Image.fromarray(
                image_patch
            ).save(
                image_dir
                / f"{stem}_in.png",
                compress_level=1,
            )

            Image.fromarray(
                mask_patch.astype(
                    np.uint8
                )
            ).save(
                mask_dir
                / f"{stem}_out.png",
                compress_level=1,
            )

            total_saved += 1

        print(
            f"{split}: "
            f"{len(split_rows)} pairs saved"
        )

    print(
        f"Total: {total_saved} pairs"
    )

    print(
        "Normalization: "
        f"log10_low={log_low:.6f}, "
        f"log10_high={log_high:.6f}"
    )


def main() -> None:
    print("=" * 88)
    print("Save finalized train/val/test patch + mask pairs")
    print(f"Manifest root: {MANIFEST_ROOT}")
    print(f"Output root  : {OUTPUT_ROOT}")
    print("=" * 88)

    for region in REGIONS:
        save_region(
            region
        )

    print()
    print("=" * 88)
    print("Finished.")
    print("=" * 88)


if __name__ == "__main__":
    main()
