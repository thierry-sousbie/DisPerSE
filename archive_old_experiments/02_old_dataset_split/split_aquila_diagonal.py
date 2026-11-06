"""
generate_all_regions_png_with_viz_fixed.py

修正版の目的
------------
1. train / val / test は乱数を使わず、決定論的な斜め帯分割にする。
2. 同一フィラメントが複数 split にまたがらないようにする。
   - まず、フィラメントを切らない境界候補を優先する。
   - そのような候補が存在しない場合だけ、切られるフィラメント全体を IGNORE にする。
3. 画像外周、NaN / inf、0以下の画素、その周辺を無効画素として扱う。
4. 帯幅を「面積 ÷ 投影長」ではなく、境界間の実際の垂直距離で評価する。
5. MAX_INVALID_RATIO と MIN_SPLIT_PURITY を矛盾させず、
   「無効画素は最大5%、有効画素は1つのsplitだけ」という条件でパッチを保存する。
6. 8近傍でフィラメント連結成分を求める。
7. 画素面積だけでなく、実際に取得可能なパッチ数も候補選択に利用する。
8. 最終分割後に、同一フィラメントの split 混在がないことを検証する。

注意
----
- 「フィラメントを途中で切らない」は train / val / test 間の分割について保証する。
- 256x256 パッチの端でフィラメントが切れることまで禁止したい場合は、
  REQUIRE_COMPLETE_FILAMENT_IN_PATCH = True にする。ただし取得可能パッチ数は大幅に減る。
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.ndimage import binary_dilation, label as ndi_label


# ============================================================
# CONFIG
# ============================================================

DATA_SOURCES = {
    #"aquila": {
    #    "fits": "/mnt/sp4t/fujiwara.r/HGBS_aquilaM2_hires_column_density_map.fits",
    #    "label": "/mnt/sp4t/fujiwara.r/HGBS_aquilaM2_hires_column_density_map/vis_3class_map/HGBS_aquilaM2_hires_column_density_map_label3class.fits",
    #},
    #"ic5146": {
    #    "fits": "/mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map.fits",
    #    "label": "/mnt/sp4t/fujiwara.r/HGBS_ic5146_hires_column_density_map/vis_3class_map/HGBS_ic5146_hires_column_density_map_label3class.fits",
    #},
    #"orionB": {
    #     "fits": "/mnt/sp4t/fujiwara.r/HGBS_orionB_hires_column_density_map.fits",
    #     "label": "/mnt/sp4t/fujiwara.r/HGBS_orionB_hires_column_density_map/vis_3class_map/HGBS_aorionB_hires_column_density_map_label3class.fits",
    # },
    # "taurus": {
    #     "fits": "/mnt/sp4t/fujiwara.r/HGBS_tauN3_hires_column_density_map.fits",
    #     "label": "/mnt/sp4t/fujiwara.r/HGBS_tauN3_hires_column_density_map/vis_3class_map/HGBS_tauN3_hires_column_density_map_label3class.fits"
    #},
     "polaris": {
         "fits": "/mnt/sp4t/fujiwara.r/HGBS_polaris_hires_column_density_map.fits",
         "label": "/mnt/sp4t/fujiwara.r/HGBS_polaris_hires_column_density_map/vis_3class_map/HGBS_polaris_hires_column_density_map_label3class.fits"
     },
}

OUTPUT_DIR = Path("/mnt/sp4t/fujiwara.r/split")
VIZ_DIR = OUTPUT_DIR / "split_visualizations"
PARAMS_DIR = OUTPUT_DIR / "split_parameters"

# 画像外周を無効にする幅
EDGE_MARGIN_PX = 32

# 元の NaN / inf / 0以下の周囲を追加で無効にする幅
NO_DATA_BUFFER_PX = 32

IGNORE_LABEL = 255
VALID_LABELS = (0, 1, 2)
BALANCE_TARGET_CLASSES = (1, 2)

# 目標分割比率。最低値ではなく、スコア上の目標。
TARGET_SPLIT_RATIOS = {
    "train": 0.60,
    "val": 0.20,
    "test": 0.20,
}

# 候補を棄却する最低面積比率
MIN_TRAIN_RATIO = 0.50
MIN_VAL_TEST_RATIO = 0.10

# 斜め分割探索
RATIO_SEARCH_STEPS = 15
THETA_RANGE_DEG = (15.0, 165.0)
THETA_STEP_DEG = 2.0

# パッチ設定
PATCH_SIZE = 256
STRIDE = 128

# パッチ内の無効画素を最大5%まで許可する。
# ただし残りの有効画素は、MIN_SPLIT_PURITY=1.0なら全て同一splitでなければならない。
MAX_INVALID_RATIO = 0.05
MIN_SPLIT_PURITY = 1.0

# 斜め帯の幅に対する安全係数。
# 軸平行の256x256パッチを斜め帯へ完全に入れるには、
# 投影方向に (PATCH_SIZE-1) * (|cos(theta)| + |sin(theta)|) が必要。
BAND_WIDTH_SAFETY_FACTOR = 1.05

# フィラメント境界処理
REQUIRE_NO_FILAMENT_CUT_CANDIDATE = True
ALLOW_CUT_FILAMENT_FALLBACK = True
FILAMENT_CONNECTIVITY = 8

# 候補探索後、上位候補だけ実際のパッチ数で再評価する。
PATCH_AWARE_TOP_K = 12

# パッチ数がこれ未満の候補には大きなペナルティを与える。
# データサイズに応じて調整可能。
MIN_PATCHES_PER_SPLIT = {
    "train": 1,
    "val": 1,
    "test": 1,
}

# パッチ端でフィラメントが切れるパッチまで拒否するか。
# True は非常に厳しく、特に長いフィラメントではパッチ数が激減する。
REQUIRE_COMPLETE_FILAMENT_IN_PATCH = False

# 1フォルダあたりの最大画像数
MAX_IMAGES_PER_DIR = 1000

# スコア重み
CLASS_BALANCE_WEIGHT = 1.0
AREA_RATIO_WEIGHT = 0.35
PATCH_RATIO_WEIGHT = 0.50
PATCH_YIELD_WEIGHT = 0.20
CUT_PIXEL_WEIGHT = 8.0
INSUFFICIENT_PATCH_PENALTY = 1_000_000.0

# 数値安定化
EPS = 1e-12


# ============================================================
# データ構造
# ============================================================


@dataclass(frozen=True)
class Candidate:
    theta_deg: float
    c1: float
    c2: float
    orientation: str
    preliminary_score: float
    class_score: float
    area_score: float
    cut_component_count: int
    cut_pixel_count: int
    width_low: float
    width_mid: float
    width_high: float
    required_band_width: float


@dataclass
class CandidateEvaluation:
    candidate: Candidate
    final_score: float
    split_masks: Dict[str, np.ndarray]
    refined_valid_mask: np.ndarray
    stats: Dict[str, object]


# ============================================================
# 基本ユーティリティ
# ============================================================


def validate_config() -> None:
    ratio_sum = sum(TARGET_SPLIT_RATIOS.values())
    if not np.isclose(ratio_sum, 1.0):
        raise ValueError(
            f"TARGET_SPLIT_RATIOS の合計は1.0である必要があります: {ratio_sum}"
        )

    if PATCH_SIZE <= 0 or STRIDE <= 0:
        raise ValueError("PATCH_SIZE と STRIDE は正の整数である必要があります。")

    if not (0.0 <= MAX_INVALID_RATIO < 1.0):
        raise ValueError("MAX_INVALID_RATIO は 0以上1未満にしてください。")

    if not (0.0 < MIN_SPLIT_PURITY <= 1.0):
        raise ValueError("MIN_SPLIT_PURITY は 0より大きく1以下にしてください。")

    if FILAMENT_CONNECTIVITY not in (4, 8):
        raise ValueError("FILAMENT_CONNECTIVITY は 4 または 8 にしてください。")


def connectivity_structure(connectivity: int) -> np.ndarray:
    if connectivity == 4:
        return np.array(
            [
                [0, 1, 0],
                [1, 1, 1],
                [0, 1, 0],
            ],
            dtype=bool,
        )
    if connectivity == 8:
        return np.ones((3, 3), dtype=bool)
    raise ValueError(f"Unsupported connectivity: {connectivity}")


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator > 0 else 0.0


def relative_distribution_error(
    observed: np.ndarray,
    target: np.ndarray,
) -> float:
    return float(np.sum(np.abs(observed - target) / (target + EPS)))


def area_ratio_error(actual: Mapping[str, float]) -> float:
    error = 0.0
    for split_name, target in TARGET_SPLIT_RATIOS.items():
        error += abs(actual[split_name] - target) / (target + EPS)
    return float(error)


def load_label_map(label_path: str) -> np.ndarray:
    from astropy.io import fits

    path = Path(label_path)
    if path.suffix.lower() == ".npy":
        label_map = np.load(path)
    else:
        with fits.open(path) as hdul:
            label_map = hdul[0].data

    label_map = np.asarray(label_map, dtype=np.int32)

    unexpected = ~np.isin(label_map, (*VALID_LABELS, IGNORE_LABEL))
    if unexpected.any():
        values = np.unique(label_map[unexpected])
        raise ValueError(
            "ラベルマップに想定外の値があります。"
            f" 許可値={(*VALID_LABELS, IGNORE_LABEL)}, 検出値={values[:20]}"
        )

    return label_map


# ============================================================
# 無効領域とフィラメント連結成分
# ============================================================


def build_invalid_mask(
    col_density: np.ndarray,
    edge_margin_px: int,
    no_data_buffer_px: int,
) -> np.ndarray:
    """画像外周と欠損領域、その周囲を無効化する。"""
    if col_density.ndim != 2:
        raise ValueError(f"col_density は2次元である必要があります: {col_density.shape}")

    h, w = col_density.shape
    invalid = np.zeros((h, w), dtype=bool)

    if edge_margin_px > 0:
        margin = min(edge_margin_px, h, w)
        invalid[:margin, :] = True
        invalid[-margin:, :] = True
        invalid[:, :margin] = True
        invalid[:, -margin:] = True

    no_data = ~np.isfinite(col_density) | (col_density <= 0)

    if no_data_buffer_px > 0 and no_data.any():
        # 3x3構造要素を指定回数膨張させる。
        # 巨大な (2r+1)x(2r+1) 構造要素を作るより省メモリ。
        no_data = binary_dilation(
            no_data,
            structure=np.ones((3, 3), dtype=bool),
            iterations=no_data_buffer_px,
        )

    invalid |= no_data
    return invalid


def build_filament_component_map(
    label_map: np.ndarray,
    target_classes: Sequence[int],
    connectivity: int,
) -> Tuple[np.ndarray, int]:
    """
    クラス1・2などのフィラメントを8近傍または4近傍で連結成分化する。

    component_id_map:
        0 はフィラメント以外。
        1以上は画像全体で一意なフィラメントID。
    """
    structure = connectivity_structure(connectivity)
    component_id_map = np.zeros(label_map.shape, dtype=np.int32)
    next_component_id = 1

    for cls in target_classes:
        binary_cls = label_map == cls
        labeled, n_components = ndi_label(binary_cls, structure=structure)

        if n_components == 0:
            continue

        positive = labeled > 0
        component_id_map[positive] = labeled[positive] + next_component_id - 1
        next_component_id += n_components

    return component_id_map, next_component_id - 1


def component_bounding_boxes(
    component_id_map: np.ndarray,
    n_components: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """各連結成分の y_min, y_max, x_min, x_max を返す。"""
    h, w = component_id_map.shape
    yy, xx = np.mgrid[0:h, 0:w]

    ids = component_id_map.ravel()
    positive = ids > 0
    ids0 = ids[positive] - 1

    y_min = np.full(n_components, h, dtype=np.int32)
    y_max = np.full(n_components, -1, dtype=np.int32)
    x_min = np.full(n_components, w, dtype=np.int32)
    x_max = np.full(n_components, -1, dtype=np.int32)

    if ids0.size > 0:
        np.minimum.at(y_min, ids0, yy.ravel()[positive])
        np.maximum.at(y_max, ids0, yy.ravel()[positive])
        np.minimum.at(x_min, ids0, xx.ravel()[positive])
        np.maximum.at(x_max, ids0, xx.ravel()[positive])

    return y_min, y_max, x_min, x_max


# ============================================================
# 積分画像によるパッチ数計算
# ============================================================


def window_sum_grid(
    mask: np.ndarray,
    patch_size: int,
    stride: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    すべてのスライド位置における mask のパッチ内合計を積分画像で返す。

    Returns
    -------
    sums: shape=(len(ys), len(xs))
    ys: パッチ左下ではなく、NumPy配列上の開始y
    xs: パッチ開始x
    """
    h, w = mask.shape
    if h < patch_size or w < patch_size:
        return (
            np.empty((0, 0), dtype=np.int64),
            np.empty((0,), dtype=np.int32),
            np.empty((0,), dtype=np.int32),
        )

    ys = np.arange(0, h - patch_size + 1, stride, dtype=np.int32)
    xs = np.arange(0, w - patch_size + 1, stride, dtype=np.int32)

    integral = np.pad(
        mask.astype(np.int64).cumsum(axis=0).cumsum(axis=1),
        ((1, 0), (1, 0)),
        mode="constant",
    )

    y0 = ys
    y1 = ys + patch_size
    x0 = xs
    x1 = xs + patch_size

    sums = (
        integral[np.ix_(y1, x1)]
        - integral[np.ix_(y0, x1)]
        - integral[np.ix_(y1, x0)]
        + integral[np.ix_(y0, x0)]
    )
    return sums, ys, xs


def count_extractable_patches(
    refined_valid_mask: np.ndarray,
    split_masks: Mapping[str, np.ndarray],
) -> Tuple[Dict[str, int], Dict[str, np.ndarray], np.ndarray, np.ndarray]:
    """
    実際の保存条件と同じ条件で、各splitから取れるパッチ数を数える。

    条件:
    - 無効画素率 <= MAX_INVALID_RATIO
    - 有効画素のうち、対象splitに属する割合 >= MIN_SPLIT_PURITY
    """
    total_pixels = PATCH_SIZE * PATCH_SIZE

    invalid_counts, ys, xs = window_sum_grid(
        ~refined_valid_mask,
        PATCH_SIZE,
        STRIDE,
    )

    if invalid_counts.size == 0:
        empty = {name: np.empty((0, 0), dtype=bool) for name in split_masks}
        return {name: 0 for name in split_masks}, empty, ys, xs

    valid_counts = total_pixels - invalid_counts
    invalid_ok = (invalid_counts / total_pixels) <= MAX_INVALID_RATIO
    valid_nonzero = valid_counts > 0

    counts: Dict[str, int] = {}
    eligible_grids: Dict[str, np.ndarray] = {}

    for split_name, split_mask in split_masks.items():
        split_counts, _, _ = window_sum_grid(split_mask, PATCH_SIZE, STRIDE)

        purity = np.divide(
            split_counts,
            valid_counts,
            out=np.zeros_like(split_counts, dtype=np.float64),
            where=valid_nonzero,
        )

        eligible = invalid_ok & valid_nonzero & (purity >= MIN_SPLIT_PURITY - EPS)
        eligible_grids[split_name] = eligible
        counts[split_name] = int(eligible.sum())

    return counts, eligible_grids, ys, xs


# ============================================================
# 斜め分割候補探索
# ============================================================


def projection_required_for_patch(theta: float) -> float:
    """軸平行PATCH_SIZE正方形の、theta方向への投影幅。"""
    return (
        (PATCH_SIZE - 1)
        * (abs(math.cos(theta)) + abs(math.sin(theta)))
        * BAND_WIDTH_SAFETY_FACTOR
    )


def build_component_projection_ranges(
    projection_valid: np.ndarray,
    component_ids_valid: np.ndarray,
    n_components: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """各フィラメント連結成分の投影最小値・最大値を求める。"""
    proj_min = np.full(n_components, np.inf, dtype=np.float64)
    proj_max = np.full(n_components, -np.inf, dtype=np.float64)

    positive = component_ids_valid > 0
    if positive.any():
        ids0 = component_ids_valid[positive] - 1
        values = projection_valid[positive]
        np.minimum.at(proj_min, ids0, values)
        np.maximum.at(proj_max, ids0, values)

    return proj_min, proj_max


def component_cut_info(
    c1: float,
    c2: float,
    component_proj_min: np.ndarray,
    component_proj_max: np.ndarray,
    component_valid_sizes: np.ndarray,
) -> Tuple[np.ndarray, int, int]:
    """境界c1/c2が内部を横切るフィラメントIDと画素数を返す。"""
    exists = np.isfinite(component_proj_min) & np.isfinite(component_proj_max)
    # region判定は projection > c で右側へ入るため、
    # 同一成分に projection <= c と projection > c の両方があれば切断。
    cut_by_c1 = exists & (component_proj_min <= c1) & (c1 < component_proj_max)
    cut_by_c2 = exists & (component_proj_min <= c2) & (c2 < component_proj_max)
    cut_flags = cut_by_c1 | cut_by_c2

    cut_component_indices = np.flatnonzero(cut_flags)
    cut_component_ids = cut_component_indices + 1
    cut_pixel_count = int(component_valid_sizes[cut_component_ids].sum())

    return cut_component_ids, int(cut_component_ids.size), cut_pixel_count


def region_class_counts_from_sorted(
    cumulative_class_counts: Mapping[int, np.ndarray],
    idx1: int,
    idx2: int,
    n_points: int,
) -> Dict[str, np.ndarray]:
    """low/mid/high各帯の対象クラス画素数を高速取得する。"""
    result: Dict[str, List[int]] = {"low": [], "mid": [], "high": []}

    for cls in BALANCE_TARGET_CLASSES:
        cumulative = cumulative_class_counts[cls]
        low_count = int(cumulative[idx1])
        mid_count = int(cumulative[idx2] - cumulative[idx1])
        high_count = int(cumulative[n_points] - cumulative[idx2])

        result["low"].append(low_count)
        result["mid"].append(mid_count)
        result["high"].append(high_count)

    return {
        name: np.asarray(values, dtype=np.float64)
        for name, values in result.items()
    }


def candidate_split_mapping(orientation: str) -> Dict[str, str]:
    if orientation == "low_is_train":
        return {"train": "low", "val": "mid", "test": "high"}
    if orientation == "high_is_train":
        return {"train": "high", "val": "mid", "test": "low"}
    raise ValueError(f"Unknown orientation: {orientation}")


def search_diagonal_candidates(
    label_map: np.ndarray,
    valid_mask: np.ndarray,
    component_id_map: np.ndarray,
    n_components: int,
) -> Tuple[List[Candidate], List[Candidate]]:
    """画素分布・帯幅・フィラメント切断条件で斜め分割候補を探索する。"""
    h, w = label_map.shape
    yy, xx = np.mgrid[0:h, 0:w]

    labels_valid = label_map[valid_mask]
    x_valid = xx[valid_mask].astype(np.float64)
    y_valid = yy[valid_mask].astype(np.float64)
    component_ids_valid = component_id_map[valid_mask]

    n_valid = labels_valid.size
    if n_valid == 0:
        raise ValueError("有効画素が0です。")

    overall_props = np.asarray(
        [(labels_valid == cls).sum() / n_valid for cls in BALANCE_TARGET_CLASSES],
        dtype=np.float64,
    )

    component_valid_sizes = np.bincount(
        component_ids_valid,
        minlength=n_components + 1,
    ).astype(np.int64)

    theta_values_deg = np.arange(
        THETA_RANGE_DEG[0],
        THETA_RANGE_DEG[1],
        THETA_STEP_DEG,
        dtype=np.float64,
    )

    r1_values = np.linspace(
        MIN_VAL_TEST_RATIO,
        1.0 - 2.0 * MIN_VAL_TEST_RATIO,
        RATIO_SEARCH_STEPS,
    )

    strict_candidates: List[Candidate] = []
    relaxed_candidates: List[Candidate] = []

    for theta_deg in theta_values_deg:
        theta = math.radians(float(theta_deg))
        cos_t = math.cos(theta)
        sin_t = math.sin(theta)

        projection_valid = x_valid * cos_t + y_valid * sin_t
        order = np.argsort(projection_valid, kind="mergesort")
        proj_sorted = projection_valid[order]
        labels_sorted = labels_valid[order]

        proj_min_all = float(proj_sorted[0])
        proj_max_all = float(proj_sorted[-1])
        required_width = projection_required_for_patch(theta)

        cumulative_class_counts: Dict[int, np.ndarray] = {}
        for cls in BALANCE_TARGET_CLASSES:
            flags = (labels_sorted == cls).astype(np.int64)
            cumulative_class_counts[cls] = np.concatenate(
                [np.zeros(1, dtype=np.int64), np.cumsum(flags)]
            )

        component_proj_min, component_proj_max = build_component_projection_ranges(
            projection_valid,
            component_ids_valid,
            n_components,
        )

        for r1 in r1_values:
            r2_values = np.linspace(
                r1 + MIN_VAL_TEST_RATIO,
                1.0 - MIN_VAL_TEST_RATIO,
                RATIO_SEARCH_STEPS,
            )

            for r2 in r2_values:
                c1 = float(np.quantile(proj_sorted, r1))
                c2 = float(np.quantile(proj_sorted, r2))

                if not c1 < c2:
                    continue

                idx1 = int(np.searchsorted(proj_sorted, c1, side="right"))
                idx2 = int(np.searchsorted(proj_sorted, c2, side="right"))

                n_low = idx1
                n_mid = idx2 - idx1
                n_high = n_valid - idx2

                if min(n_low, n_mid, n_high) <= 0:
                    continue

                area_by_band = {
                    "low": n_low / n_valid,
                    "mid": n_mid / n_valid,
                    "high": n_high / n_valid,
                }

                width_low = c1 - proj_min_all
                width_mid = c2 - c1
                width_high = proj_max_all - c2

                widths_ok = (
                    width_low >= required_width
                    and width_mid >= required_width
                    and width_high >= required_width
                )
                if not widths_ok:
                    continue

                cut_ids, cut_component_count, cut_pixel_count = component_cut_info(
                    c1,
                    c2,
                    component_proj_min,
                    component_proj_max,
                    component_valid_sizes,
                )
                del cut_ids

                band_class_counts = region_class_counts_from_sorted(
                    cumulative_class_counts,
                    idx1,
                    idx2,
                    n_valid,
                )

                for orientation in ("low_is_train", "high_is_train"):
                    mapping = candidate_split_mapping(orientation)
                    split_area = {
                        split_name: area_by_band[band_name]
                        for split_name, band_name in mapping.items()
                    }

                    if split_area["train"] < MIN_TRAIN_RATIO:
                        continue
                    if split_area["val"] < MIN_VAL_TEST_RATIO:
                        continue
                    if split_area["test"] < MIN_VAL_TEST_RATIO:
                        continue

                    class_score = 0.0
                    invalid_distribution = False
                    for split_name, band_name in mapping.items():
                        band_total = {
                            "low": n_low,
                            "mid": n_mid,
                            "high": n_high,
                        }[band_name]

                        if band_total <= 0:
                            invalid_distribution = True
                            break

                        observed_props = band_class_counts[band_name] / band_total
                        class_score += relative_distribution_error(
                            observed_props,
                            overall_props,
                        )

                    if invalid_distribution:
                        continue

                    current_area_score = area_ratio_error(split_area)
                    cut_ratio = cut_pixel_count / n_valid

                    preliminary_score = (
                        CLASS_BALANCE_WEIGHT * class_score
                        + AREA_RATIO_WEIGHT * current_area_score
                        + CUT_PIXEL_WEIGHT * cut_ratio
                    )

                    candidate = Candidate(
                        theta_deg=float(theta_deg),
                        c1=c1,
                        c2=c2,
                        orientation=orientation,
                        preliminary_score=float(preliminary_score),
                        class_score=float(class_score),
                        area_score=float(current_area_score),
                        cut_component_count=cut_component_count,
                        cut_pixel_count=cut_pixel_count,
                        width_low=float(width_low),
                        width_mid=float(width_mid),
                        width_high=float(width_high),
                        required_band_width=float(required_width),
                    )

                    relaxed_candidates.append(candidate)
                    if cut_component_count == 0:
                        strict_candidates.append(candidate)

    if not relaxed_candidates:
        raise RuntimeError(
            "帯幅・面積条件を満たす斜め分割候補がありません。"
            " THETA_RANGE_DEG、比率条件、PATCH_SIZE、"
            "BAND_WIDTH_SAFETY_FACTORを見直してください。"
        )

    strict_candidates.sort(key=lambda item: item.preliminary_score)
    relaxed_candidates.sort(key=lambda item: item.preliminary_score)

    top_k = max(PATCH_AWARE_TOP_K, 1)
    strict_shortlist = strict_candidates[:top_k]

    # フォールバック候補では、strict候補との重複を除き、
    # 実際にフィラメント全体のIGNOREが必要な候補だけを残す。
    cut_shortlist = [
        candidate
        for candidate in relaxed_candidates
        if candidate.cut_component_count > 0
    ][:top_k]

    return strict_shortlist, cut_shortlist


# ============================================================
# 候補の最終評価
# ============================================================


def build_split_masks_for_candidate(
    candidate: Candidate,
    valid_mask: np.ndarray,
    component_id_map: np.ndarray,
    n_components: int,
) -> Tuple[Dict[str, np.ndarray], np.ndarray, np.ndarray]:
    """
    候補境界を画像全体へ適用し、境界で切られるフィラメント全体を無効化する。

    Returns
    -------
    split_masks
    refined_valid_mask
    cut_component_ids
    """
    h, w = valid_mask.shape
    yy, xx = np.mgrid[0:h, 0:w]

    theta = math.radians(candidate.theta_deg)
    projection_full = xx.astype(np.float64) * math.cos(theta) + yy.astype(
        np.float64
    ) * math.sin(theta)

    region_id = np.zeros((h, w), dtype=np.int8)
    region_id[projection_full > candidate.c1] = 1
    region_id[projection_full > candidate.c2] = 2

    # 各フィラメントがどのregion_idに属するかを調べる。
    component_ids_valid = component_id_map[valid_mask]
    region_ids_valid = region_id[valid_mask]

    region_bits = np.zeros(n_components + 1, dtype=np.uint8)
    positive = component_ids_valid > 0
    if positive.any():
        # region 0,1,2 を bit 1,2,4 としてORする。
        bits = (1 << region_ids_valid[positive]).astype(np.uint8)
        np.bitwise_or.at(region_bits, component_ids_valid[positive], bits)

    # 1,2,4 のいずれか単独なら切られていない。
    # 3,5,6,7 なら複数領域にまたがっている。
    cut_flags = np.zeros(n_components + 1, dtype=bool)
    active_bits = region_bits[1:]
    cut_flags[1:] = (active_bits != 0) & ((active_bits & (active_bits - 1)) != 0)
    cut_component_ids = np.flatnonzero(cut_flags)

    refined_valid_mask = valid_mask.copy()
    if cut_component_ids.size > 0:
        cut_lookup = np.zeros(n_components + 1, dtype=bool)
        cut_lookup[cut_component_ids] = True
        cut_mask = cut_lookup[component_id_map] & valid_mask
        refined_valid_mask[cut_mask] = False

    if candidate.orientation == "low_is_train":
        split_masks = {
            "train": refined_valid_mask & (region_id == 0),
            "val": refined_valid_mask & (region_id == 1),
            "test": refined_valid_mask & (region_id == 2),
        }
    elif candidate.orientation == "high_is_train":
        split_masks = {
            "train": refined_valid_mask & (region_id == 2),
            "val": refined_valid_mask & (region_id == 1),
            "test": refined_valid_mask & (region_id == 0),
        }
    else:
        raise ValueError(f"Unknown orientation: {candidate.orientation}")

    return split_masks, refined_valid_mask, cut_component_ids


def calculate_post_split_statistics(
    label_map: np.ndarray,
    valid_mask_before: np.ndarray,
    refined_valid_mask: np.ndarray,
    split_masks: Mapping[str, np.ndarray],
    cut_component_ids: np.ndarray,
) -> Dict[str, object]:
    n_before = int(valid_mask_before.sum())
    n_after = int(refined_valid_mask.sum())

    overall_labels = label_map[refined_valid_mask]
    overall_props = np.asarray(
        [
            safe_ratio(int((overall_labels == cls).sum()), overall_labels.size)
            for cls in BALANCE_TARGET_CLASSES
        ],
        dtype=np.float64,
    )

    area_ratios: Dict[str, float] = {}
    class_props: Dict[str, Dict[str, float]] = {}
    class_score = 0.0

    for split_name, split_mask in split_masks.items():
        split_count = int(split_mask.sum())
        area_ratios[split_name] = safe_ratio(split_count, n_after)

        split_labels = label_map[split_mask]
        observed = np.asarray(
            [
                safe_ratio(int((split_labels == cls).sum()), split_labels.size)
                for cls in BALANCE_TARGET_CLASSES
            ],
            dtype=np.float64,
        )

        class_props[split_name] = {
            str(cls): float(prop)
            for cls, prop in zip(BALANCE_TARGET_CLASSES, observed)
        }
        class_score += relative_distribution_error(observed, overall_props)

    patch_counts, eligible_grids, ys, xs = count_extractable_patches(
        refined_valid_mask,
        split_masks,
    )

    total_patches = sum(patch_counts.values())
    patch_ratios = {
        split_name: safe_ratio(count, total_patches)
        for split_name, count in patch_counts.items()
    }

    patch_ratio_score = area_ratio_error(patch_ratios) if total_patches > 0 else math.inf
    post_area_score = area_ratio_error(area_ratios)

    total_possible_locations = int(len(ys) * len(xs))
    patch_yield = safe_ratio(total_patches, total_possible_locations)
    cut_pixel_count = n_before - n_after
    cut_pixel_ratio = safe_ratio(cut_pixel_count, n_before)

    return {
        "valid_pixels_before": n_before,
        "valid_pixels_after": n_after,
        "cut_component_count": int(cut_component_ids.size),
        "cut_pixel_count": int(cut_pixel_count),
        "cut_pixel_ratio": float(cut_pixel_ratio),
        "area_ratios": area_ratios,
        "class_props": class_props,
        "post_class_score": float(class_score),
        "post_area_score": float(post_area_score),
        "patch_counts": patch_counts,
        "patch_ratios": patch_ratios,
        "patch_ratio_score": float(patch_ratio_score),
        "patch_yield": float(patch_yield),
        "total_possible_patch_locations": total_possible_locations,
        "eligible_grids": eligible_grids,
        "grid_ys": ys,
        "grid_xs": xs,
    }


def score_evaluated_candidate(stats: Mapping[str, object]) -> float:
    patch_counts = stats["patch_counts"]
    assert isinstance(patch_counts, dict)

    score = (
        CLASS_BALANCE_WEIGHT * float(stats["post_class_score"])
        + AREA_RATIO_WEIGHT * float(stats["post_area_score"])
        + PATCH_RATIO_WEIGHT * float(stats["patch_ratio_score"])
        + PATCH_YIELD_WEIGHT * (1.0 - float(stats["patch_yield"]))
        + CUT_PIXEL_WEIGHT * float(stats["cut_pixel_ratio"])
    )

    for split_name, minimum in MIN_PATCHES_PER_SPLIT.items():
        actual = int(patch_counts.get(split_name, 0))
        if actual < minimum:
            score += INSUFFICIENT_PATCH_PENALTY * (minimum - actual)

    if not np.isfinite(score):
        return math.inf
    return float(score)


def candidate_meets_minimum_patch_counts(stats: Mapping[str, object]) -> bool:
    patch_counts = stats["patch_counts"]
    assert isinstance(patch_counts, dict)
    return all(
        int(patch_counts.get(split_name, 0)) >= minimum
        for split_name, minimum in MIN_PATCHES_PER_SPLIT.items()
    )


def select_best_candidate(
    label_map: np.ndarray,
    valid_mask: np.ndarray,
    component_id_map: np.ndarray,
    n_components: int,
    candidates: Sequence[Candidate],
) -> CandidateEvaluation:
    best_evaluation: Optional[CandidateEvaluation] = None

    for rank, candidate in enumerate(candidates, start=1):
        split_masks, refined_valid_mask, cut_component_ids = build_split_masks_for_candidate(
            candidate,
            valid_mask,
            component_id_map,
            n_components,
        )

        stats = calculate_post_split_statistics(
            label_map,
            valid_mask,
            refined_valid_mask,
            split_masks,
            cut_component_ids,
        )
        final_score = score_evaluated_candidate(stats)

        print(
            f"[candidate {rank:02d}] theta={candidate.theta_deg:.1f}, "
            f"orientation={candidate.orientation}, "
            f"final_score={final_score:.4f}, "
            f"patches={stats['patch_counts']}, "
            f"cut_components={stats['cut_component_count']}"
        )

        evaluation = CandidateEvaluation(
            candidate=candidate,
            final_score=final_score,
            split_masks=split_masks,
            refined_valid_mask=refined_valid_mask,
            stats=stats,
        )

        if best_evaluation is None or final_score < best_evaluation.final_score:
            best_evaluation = evaluation

    if best_evaluation is None or not np.isfinite(best_evaluation.final_score):
        raise RuntimeError(
            "パッチ数まで考慮した有効な分割候補が見つかりませんでした。"
            " MIN_PATCHES_PER_SPLIT、PATCH_SIZE、STRIDE、比率条件を見直してください。"
        )

    return best_evaluation


def search_diagonal_split(
    label_map: np.ndarray,
    valid_mask: np.ndarray,
    component_id_map: np.ndarray,
    n_components: int,
) -> CandidateEvaluation:
    strict_candidates, cut_candidates = search_diagonal_candidates(
        label_map,
        valid_mask,
        component_id_map,
        n_components,
    )

    if REQUIRE_NO_FILAMENT_CUT_CANDIDATE and strict_candidates:
        print(
            f"[info] フィラメントを1本も切らない候補をパッチ数で再評価: "
            f"{len(strict_candidates)} candidates"
        )
        strict_best = select_best_candidate(
            label_map,
            valid_mask,
            component_id_map,
            n_components,
            strict_candidates,
        )
        if candidate_meets_minimum_patch_counts(strict_best.stats):
            return strict_best

        if not ALLOW_CUT_FILAMENT_FALLBACK:
            raise RuntimeError(
                "フィラメント非切断候補は見つかりましたが、"
                "必要なパッチ数を確保できませんでした。"
            )

        print(
            "[warning] フィラメント非切断候補では必要パッチ数を確保できないため、"
            "切断対象フィラメント全体をIGNOREにする候補へフォールバックします。"
        )

    elif REQUIRE_NO_FILAMENT_CUT_CANDIDATE:
        if not ALLOW_CUT_FILAMENT_FALLBACK:
            raise RuntimeError(
                "フィラメントを切らない候補がありません。"
                " THETA_RANGE_DEG、比率条件、PATCH_SIZE、"
                "BAND_WIDTH_SAFETY_FACTORを見直してください。"
            )
        print(
            "[warning] フィラメントを切らない候補が見つからなかったため、"
            "切断対象フィラメント全体をIGNOREにする候補へフォールバックします。"
        )

    if REQUIRE_NO_FILAMENT_CUT_CANDIDATE:
        if not cut_candidates:
            raise RuntimeError(
                "フォールバック可能な分割候補がありません。"
                " 分割条件を見直してください。"
            )
        print(
            f"[info] IGNOREフォールバック候補をパッチ数で再評価: "
            f"{len(cut_candidates)} candidates"
        )
        fallback_best = select_best_candidate(
            label_map,
            valid_mask,
            component_id_map,
            n_components,
            cut_candidates,
        )
        if not candidate_meets_minimum_patch_counts(fallback_best.stats):
            raise RuntimeError(
                "フィラメント全体をIGNOREにする候補でも、"
                "必要なパッチ数を確保できませんでした。"
                " MIN_PATCHES_PER_SPLIT、PATCH_SIZE、STRIDE、比率条件を見直してください。"
            )
        return fallback_best

    # REQUIRE_NO_FILAMENT_CUT_CANDIDATE=False の場合は、
    # strict候補とcut候補をまとめてスコア順に評価する。
    candidates = sorted(
        [*strict_candidates, *cut_candidates],
        key=lambda item: item.preliminary_score,
    )[: max(PATCH_AWARE_TOP_K, 1)]
    if not candidates:
        raise RuntimeError("条件を満たす斜め分割候補がありません。")

    print(f"[info] パッチ数で再評価する候補数: {len(candidates)}")
    best = select_best_candidate(
        label_map,
        valid_mask,
        component_id_map,
        n_components,
        candidates,
    )
    if not candidate_meets_minimum_patch_counts(best.stats):
        raise RuntimeError(
            "必要なパッチ数を確保できる候補がありませんでした。"
        )
    return best


# ============================================================
# 分割検証
# ============================================================


def validate_final_split(
    split_masks: Mapping[str, np.ndarray],
    refined_valid_mask: np.ndarray,
    component_id_map: np.ndarray,
    n_components: int,
) -> None:
    """最終分割に重複・欠落・同一フィラメント混在がないか検証する。"""
    names = ("train", "val", "test")

    for name in names:
        if name not in split_masks:
            raise AssertionError(f"split_masks に {name} がありません。")

    overlap = (
        split_masks["train"].astype(np.uint8)
        + split_masks["val"].astype(np.uint8)
        + split_masks["test"].astype(np.uint8)
    )

    if np.any(overlap > 1):
        raise AssertionError("train / val / test のマスクが重複しています。")

    union = split_masks["train"] | split_masks["val"] | split_masks["test"]
    if not np.array_equal(union, refined_valid_mask):
        missing = int((refined_valid_mask & ~union).sum())
        extra = int((union & ~refined_valid_mask).sum())
        raise AssertionError(
            f"splitマスクの和と refined_valid_mask が一致しません: "
            f"missing={missing}, extra={extra}"
        )

    # 各フィラメントが複数splitに存在しないことを確認する。
    split_bits = np.zeros(n_components + 1, dtype=np.uint8)
    for bit, split_name in zip((1, 2, 4), names):
        ids = component_id_map[split_masks[split_name]]
        ids = ids[ids > 0]
        if ids.size > 0:
            np.bitwise_or.at(
                split_bits,
                ids,
                np.full(ids.shape, bit, dtype=np.uint8),
            )

    active = split_bits[1:]
    mixed = np.flatnonzero((active != 0) & ((active & (active - 1)) != 0)) + 1
    if mixed.size > 0:
        raise AssertionError(
            "同一フィラメントが複数splitに混在しています。"
            f" component_ids={mixed[:20].tolist()}"
        )

    print("[check] split重複なし・欠落なし・同一フィラメント混在なし")


# ============================================================
# 可視化とパラメータ保存
# ============================================================


def save_split_visualization(
    region_name: str,
    col_density: np.ndarray,
    split_masks: Mapping[str, np.ndarray],
    raw_valid_mask: np.ndarray,
    refined_valid_mask: np.ndarray,
    evaluation: CandidateEvaluation,
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    h, w = col_density.shape

    split_map = np.zeros((h, w), dtype=np.uint8)
    split_map[split_masks["train"]] = 1
    split_map[split_masks["val"]] = 2
    split_map[split_masks["test"]] = 3

    cut_filaments = raw_valid_mask & (~refined_valid_mask)
    split_map[cut_filaments] = 4

    finite_positive = np.isfinite(col_density) & (col_density > 0)
    log_density = np.full(col_density.shape, np.nan, dtype=np.float64)
    log_density[finite_positive] = np.log10(col_density[finite_positive])

    plt.figure(figsize=(12, 10))
    plt.imshow(log_density, cmap="gray", origin="lower", alpha=0.65)

    cmap = matplotlib.colormaps["Set1"].resampled(4)
    masked_split = np.ma.masked_where(split_map == 0, split_map)
    plt.imshow(
        masked_split,
        cmap=cmap,
        vmin=1,
        vmax=4,
        origin="lower",
        alpha=0.50,
    )

    import matplotlib.patches as mpatches

    legend_patches = [
        mpatches.Patch(color=cmap(0.0), label="Train Region"),
        mpatches.Patch(color=cmap(1.0 / 3.0), label="Val Region"),
        mpatches.Patch(color=cmap(2.0 / 3.0), label="Test Region"),
        mpatches.Patch(color=cmap(1.0), label="Cut Filaments / Ignored"),
    ]
    plt.legend(handles=legend_patches, loc="upper right", fontsize=11)

    candidate = evaluation.candidate
    stats = evaluation.stats
    plt.title(
        f"Diagonal Split Map: {region_name}\n"
        f"Theta={candidate.theta_deg:.1f} deg | "
        f"Score={evaluation.final_score:.4f} | "
        f"Patches={stats['patch_counts']}",
        fontsize=13,
    )
    plt.tight_layout()

    save_path = out_dir / f"{region_name}_split_map.png"
    plt.savefig(save_path, dpi=200)
    plt.close()
    print(f"-> Visualization saved: {save_path}")


def json_serializable_stats(evaluation: CandidateEvaluation) -> Dict[str, object]:
    candidate = evaluation.candidate
    stats = evaluation.stats

    serializable_stats = {
        key: value
        for key, value in stats.items()
        if key not in {"eligible_grids", "grid_ys", "grid_xs"}
    }

    return {
        "candidate": {
            "theta_deg": candidate.theta_deg,
            "c1": candidate.c1,
            "c2": candidate.c2,
            "orientation": candidate.orientation,
            "preliminary_score": candidate.preliminary_score,
            "class_score": candidate.class_score,
            "area_score": candidate.area_score,
            "cut_component_count_before_refinement": candidate.cut_component_count,
            "cut_pixel_count_before_refinement": candidate.cut_pixel_count,
            "width_low": candidate.width_low,
            "width_mid": candidate.width_mid,
            "width_high": candidate.width_high,
            "required_band_width": candidate.required_band_width,
        },
        "final_score": evaluation.final_score,
        "stats": serializable_stats,
        "config": {
            "target_split_ratios": TARGET_SPLIT_RATIOS,
            "patch_size": PATCH_SIZE,
            "stride": STRIDE,
            "max_invalid_ratio": MAX_INVALID_RATIO,
            "min_split_purity": MIN_SPLIT_PURITY,
            "edge_margin_px": EDGE_MARGIN_PX,
            "no_data_buffer_px": NO_DATA_BUFFER_PX,
            "filament_connectivity": FILAMENT_CONNECTIVITY,
            "require_complete_filament_in_patch": REQUIRE_COMPLETE_FILAMENT_IN_PATCH,
        },
    }


def save_split_parameters(
    region_name: str,
    evaluation: CandidateEvaluation,
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    save_path = out_dir / f"{region_name}_split_params.json"

    with save_path.open("w", encoding="utf-8") as file:
        json.dump(
            json_serializable_stats(evaluation),
            file,
            ensure_ascii=False,
            indent=2,
        )

    print(f"-> Split parameters saved: {save_path}")


# ============================================================
# パッチ保存
# ============================================================


def save_patch_for_mmseg(
    image_patch: np.ndarray,
    label_patch: np.ndarray,
    base_split_dir: Path,
    region_name: str,
    patch_idx: int,
    file_prefix: str,
) -> None:
    part_idx = patch_idx // MAX_IMAGES_PER_DIR
    part_dir_name = f"part_{part_idx:03d}"

    img_dir = base_split_dir / "in" / region_name / part_dir_name
    lbl_dir = base_split_dir / "out" / region_name / part_dir_name

    img_dir.mkdir(parents=True, exist_ok=True)
    lbl_dir.mkdir(parents=True, exist_ok=True)

    # NaNを直接 uint8 化すると警告や不定値が出るため、有効画素だけ正規化する。
    norm_img = np.zeros(image_patch.shape, dtype=np.uint8)
    valid_pixel_mask = np.isfinite(image_patch)

    if valid_pixel_mask.any():
        values = image_patch[valid_pixel_mask]
        min_val = float(values.min())
        max_val = float(values.max())

        if max_val > min_val:
            scaled = (values - min_val) / (max_val - min_val) * 255.0
            norm_img[valid_pixel_mask] = np.clip(scaled, 0.0, 255.0).astype(
                np.uint8
            )

    label_3class = label_patch.astype(np.int32, copy=True)
    unexpected = ~np.isin(label_3class, (*VALID_LABELS, IGNORE_LABEL))
    label_3class[unexpected] = IGNORE_LABEL
    label_3class = label_3class.astype(np.uint8)

    Image.fromarray(norm_img, mode="L").save(
        img_dir / f"{file_prefix}_in.png"
    )
    Image.fromarray(label_3class, mode="L").save(
        lbl_dir / f"{file_prefix}_out.png"
    )


def patch_contains_only_complete_filaments(
    component_patch: np.ndarray,
    y: int,
    x: int,
    component_boxes: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
) -> bool:
    """
    パッチ内に現れる全フィラメントが、そのパッチ内へ完全に収まるか確認する。

    REQUIRE_COMPLETE_FILAMENT_IN_PATCH=True の場合だけ使用する。
    """
    component_ids = np.unique(component_patch)
    component_ids = component_ids[component_ids > 0]
    if component_ids.size == 0:
        return True

    y_min, y_max, x_min, x_max = component_boxes
    ids0 = component_ids - 1

    y_end = y + PATCH_SIZE - 1
    x_end = x + PATCH_SIZE - 1

    fully_inside = (
        (y_min[ids0] >= y)
        & (y_max[ids0] <= y_end)
        & (x_min[ids0] >= x)
        & (x_max[ids0] <= x_end)
    )
    return bool(np.all(fully_inside))


def extract_and_save_patches(
    region_name: str,
    col_density: np.ndarray,
    label_map: np.ndarray,
    split_masks: Mapping[str, np.ndarray],
    refined_valid_mask: np.ndarray,
    component_id_map: np.ndarray,
    n_components: int,
    out_base_dir: Path,
) -> Dict[str, int]:
    """最終条件に従ってパッチを切り出して保存する。"""
    h, w = col_density.shape
    total_pixels = PATCH_SIZE * PATCH_SIZE

    invalid_map = ~refined_valid_mask
    counts = {"train": 0, "val": 0, "test": 0}

    component_boxes = None
    if REQUIRE_COMPLETE_FILAMENT_IN_PATCH:
        component_boxes = component_bounding_boxes(component_id_map, n_components)

    for y in range(0, h - PATCH_SIZE + 1, STRIDE):
        for x in range(0, w - PATCH_SIZE + 1, STRIDE):
            y_end = y + PATCH_SIZE
            x_end = x + PATCH_SIZE

            invalid_count = int(invalid_map[y:y_end, x:x_end].sum())
            invalid_ratio = invalid_count / total_pixels
            if invalid_ratio > MAX_INVALID_RATIO:
                continue

            valid_count = total_pixels - invalid_count
            if valid_count <= 0:
                continue

            if REQUIRE_COMPLETE_FILAMENT_IN_PATCH:
                assert component_boxes is not None
                component_patch = component_id_map[y:y_end, x:x_end]
                if not patch_contains_only_complete_filaments(
                    component_patch,
                    y,
                    x,
                    component_boxes,
                ):
                    continue

            selected_split: Optional[str] = None
            for split_name in ("train", "val", "test"):
                split_count = int(
                    split_masks[split_name][y:y_end, x:x_end].sum()
                )
                purity = split_count / valid_count

                if purity >= MIN_SPLIT_PURITY - EPS:
                    selected_split = split_name
                    break

            if selected_split is None:
                continue

            file_prefix = f"{region_name}_y{y:04d}_x{x:04d}"
            save_patch_for_mmseg(
                image_patch=col_density[y:y_end, x:x_end],
                label_patch=label_map[y:y_end, x:x_end],
                base_split_dir=out_base_dir / selected_split,
                region_name=region_name,
                patch_idx=counts[selected_split],
                file_prefix=file_prefix,
            )
            counts[selected_split] += 1

    return counts


# ============================================================
# 領域処理
# ============================================================


def process_region(
    region_name: str,
    fits_path: str,
    label_path: str,
    out_base_dir: Path,
) -> Dict[str, int]:
    from astropy.io import fits

    print(f"\n--- Processing Region: {region_name} ---")

    with fits.open(fits_path) as hdul:
        col_density = np.asarray(hdul[0].data, dtype=np.float64)

    label_map = load_label_map(label_path)

    if col_density.shape != label_map.shape:
        raise ValueError(
            f"画像とラベルのshapeが一致しません: "
            f"image={col_density.shape}, label={label_map.shape}"
        )

    invalid_mask = build_invalid_mask(
        col_density,
        EDGE_MARGIN_PX,
        NO_DATA_BUFFER_PX,
    )

    # ラベル側のIGNOREも無効領域へ含める。
    invalid_mask |= label_map == IGNORE_LABEL
    raw_valid_mask = ~invalid_mask

    component_id_map, n_components = build_filament_component_map(
        label_map,
        BALANCE_TARGET_CLASSES,
        FILAMENT_CONNECTIVITY,
    )
    print(f"[info] filament components: {n_components}")

    evaluation = search_diagonal_split(
        label_map,
        raw_valid_mask,
        component_id_map,
        n_components,
    )

    split_masks = evaluation.split_masks
    refined_valid_mask = evaluation.refined_valid_mask

    validate_final_split(
        split_masks,
        refined_valid_mask,
        component_id_map,
        n_components,
    )

    save_split_visualization(
        region_name,
        col_density,
        split_masks,
        raw_valid_mask,
        refined_valid_mask,
        evaluation,
        VIZ_DIR,
    )
    save_split_parameters(region_name, evaluation, PARAMS_DIR)

    # 無効領域を最終的に NaN / IGNORE にする。
    col_density_masked = col_density.copy()
    label_map_masked = label_map.copy()
    col_density_masked[~refined_valid_mask] = np.nan
    label_map_masked[~refined_valid_mask] = IGNORE_LABEL

    counts = extract_and_save_patches(
        region_name,
        col_density_masked,
        label_map_masked,
        split_masks,
        refined_valid_mask,
        component_id_map,
        n_components,
        out_base_dir,
    )

    predicted_counts = evaluation.stats["patch_counts"]
    if counts != predicted_counts:
        print(
            "[warning] 候補評価時の予測パッチ数と実保存数が異なります。"
            f" predicted={predicted_counts}, actual={counts}"
        )

    candidate = evaluation.candidate
    print(
        f"[selected] theta={candidate.theta_deg:.1f} deg, "
        f"c1={candidate.c1:.3f}, c2={candidate.c2:.3f}, "
        f"orientation={candidate.orientation}"
    )
    print(
        f"[selected] band widths="
        f"({candidate.width_low:.1f}, {candidate.width_mid:.1f}, "
        f"{candidate.width_high:.1f}) px, "
        f"required>={candidate.required_band_width:.1f}px"
    )
    print(
        f"Patches generated [{region_name}]: "
        f"Train={counts['train']}, Val={counts['val']}, Test={counts['test']}"
    )

    return counts


# ============================================================
# メイン処理
# ============================================================


def main() -> None:
    validate_config()

    total_summary = {"train": 0, "val": 0, "test": 0}

    for region_name, paths in DATA_SOURCES.items():
        counts = process_region(
            region_name,
            paths["fits"],
            paths["label"],
            OUTPUT_DIR,
        )
        for split_name in total_summary:
            total_summary[split_name] += counts[split_name]

    print("\n============================================")
    print("      ALL REGIONS EXTRACTION COMPLETE       ")
    print("============================================")
    print(
        "Total Combined Patches -> "
        f"Train: {total_summary['train']}, "
        f"Val: {total_summary['val']}, "
        f"Test: {total_summary['test']}"
    )
    print(f"Visualizations saved in: {VIZ_DIR}")
    print(f"Split parameters saved in: {PARAMS_DIR}")


if __name__ == "__main__":
    main()


