
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DisPerSE Loose版

- 5領域対応
- NSIG=0.5
- Loose Av / LEN_CUT_PC
- VTK中心線を1pixel幅でFITS保存
- 元柱密度画像に黄色線を重ねてPNG保存
- 物理幅の膨張は別スクリプトで行う
"""

import os
import csv
import glob
import subprocess
import tempfile

import numpy as np
import cv2
import vtk

from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import binary_dilation
from PIL import Image


# ============================================================
# 設定
# ============================================================

ROOT = "/mnt/sp4t/fujiwara.r"

REGION = os.environ.get(
    "REGION", "aquila"
).lower()

MODE = "loose"
NSIG = 0.5

AV_UNIT = 0.94e21

DISPERSE_BIN_DIR = (
    "/home1/fujiwara.r/"
    "PycharmProjects/DisPerSE/build/src/mse"
)

TARGETS = {
    "aquila": {
        "fits_name": "HGBS_aquilaM2_hires_column_density_map",
        "distance_pc": 440.0,
        "av": 6,
        "len_cut_pc": 0.15,
    },
    "ic5146": {
        "fits_name": "HGBS_ic5146_hires_column_density_map",
        "distance_pc": 600.0,
        "av": 2,
        "len_cut_pc": 0.15,
    },
    "orionb": {
        "fits_name": "HGBS_orionB_hires_column_density_map",
        "distance_pc": 400.0,
        "av": 3,
        "len_cut_pc": 0.30,
    },
    "taurus": {
        "fits_name": "HGBS_tauN3_hires_column_density_map",
        "distance_pc": 140.0,
        "av": 3,
        "len_cut_pc": 0.30,
    },
    "polaris": {
        "fits_name": "HGBS_polaris_hires_column_density_map",
        "distance_pc": 150.0,
        "av": 1,
        "len_cut_pc": 0.10,
    },
}


# ============================================================
# 共通関数
# ============================================================

def find_image_path(fits_name):
    candidates = [
        os.path.join(ROOT, fits_name + ".fits"),
        os.path.join(ROOT, fits_name, fits_name + ".fits"),
    ]

    for path in candidates:
        if os.path.isfile(path):
            return path

    raise FileNotFoundError(
        f"元FITSが見つかりません: {candidates}"
    )


def bin_path(name):
    if DISPERSE_BIN_DIR:
        return os.path.join(DISPERSE_BIN_DIR, name)
    return name


def run_command(cmd):
    print("RUN:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def find_single_file(pattern, label):
    files = sorted(glob.glob(pattern))

    if len(files) != 1:
        raise RuntimeError(
            f"{label}が一意に見つかりません: {files}"
        )

    return files[0]


def read_vtk_polylines(vtk_path):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(vtk_path)
    reader.Update()

    poly = reader.GetOutput()

    lines = poly.GetLines()
    lines.InitTraversal()

    ids = vtk.vtkIdList()
    filaments = []

    while lines.GetNextCell(ids):
        pts = [
            poly.GetPoint(ids.GetId(i))[:2]
            for i in range(ids.GetNumberOfIds())
        ]

        if len(pts) >= 2:
            filaments.append(
                np.asarray(pts, dtype=float)
            )

    return filaments


def make_gray_background(image):
    valid = np.isfinite(image) & (image > 0)

    if not np.any(valid):
        raise ValueError("有効な画像画素がありません")

    log_values = np.log10(image[valid])

    lo, hi = np.percentile(
        log_values, [1.0, 99.5]
    )

    if hi <= lo:
        hi = lo + 1e-6

    gray = np.zeros(image.shape, dtype=np.uint8)

    gray[valid] = np.clip(
        (log_values - lo) / (hi - lo) * 255,
        0,
        255
    ).astype(np.uint8)

    return gray


# ============================================================
# メイン
# ============================================================

def main():
    if REGION not in TARGETS:
        raise ValueError(f"不明な領域: {REGION}")

    cfg = TARGETS[REGION]

    fits_name = cfg["fits_name"]
    distance_pc = cfg["distance_pc"]
    av = cfg["av"]
    len_cut_pc = cfg["len_cut_pc"]

    av_thresh = av * AV_UNIT

    image_path = find_image_path(fits_name)

    with fits.open(image_path, memmap=True) as hdul:
        header = hdul[0].header.copy()
        image = np.squeeze(
            np.asarray(hdul[0].data, dtype=np.float64)
        ).copy()

    if image.ndim != 2:
        raise ValueError(
            f"元画像が2次元ではありません: {image.shape}"
        )

    # --------------------------------------------------------
    # pc/pixel
    # --------------------------------------------------------

    wcs = WCS(header).celestial

    scales_deg = np.asarray(
        proj_plane_pixel_scales(wcs)[:2],
        dtype=float
    )

    if (
        np.any(~np.isfinite(scales_deg))
        or np.any(scales_deg <= 0)
    ):
        raise ValueError("WCS pixel scaleが不正です")

    pc_per_pixel = (
        np.deg2rad(float(scales_deg.mean()))
        * distance_pc
    )

    print("================================")
    print("処理開始:", REGION, "/ Loose")
    print("距離:", distance_pc, "pc")
    print("Av:", av)
    print("Av threshold:", av_thresh)
    print("NSIG:", NSIG)
    print("LEN_CUT_PC:", len_cut_pc)
    print("pc_per_pixel:", pc_per_pixel)
    print("================================")

    # --------------------------------------------------------
    # Avマスク
    # --------------------------------------------------------

    valid = np.isfinite(image) & (image > 0)

    masked = image.copy()

    masked[
        ~valid | (image < av_thresh)
    ] = np.nan

    # --------------------------------------------------------
    # 保存先
    # --------------------------------------------------------

    base_dir = os.path.join(ROOT, fits_name)

    run_tag = (
        f"av{av}_nsig{NSIG:g}_len{len_cut_pc:g}"
    )

    raw_root = os.path.join(
        base_dir,
        "disperse",
        f"raw_loose_fixed_{run_tag}"
    )

    full_dir = os.path.join(
        base_dir,
        "full_disperse_loose_fixed"
    )

    vis_dir = os.path.join(
        base_dir,
        "vis_disperse_loose_fixed"
    )

    for directory in (raw_root, full_dir, vis_dir):
        os.makedirs(directory, exist_ok=True)

    # 中間結果の取り違えを防ぐため毎回新しいフォルダ
    run_dir = tempfile.mkdtemp(
        prefix="run_",
        dir=raw_root
    )

    masked_path = os.path.join(
        run_dir,
        f"{fits_name}_av{av}_masked.fits"
    )

    fits.writeto(
        masked_path,
        masked,
        header,
        overwrite=False
    )

    print("Masked FITS:", masked_path)

    # --------------------------------------------------------
    # Step 1: mse
    # --------------------------------------------------------

    run_command([
        bin_path("mse"),
        masked_path,
        "-upSkl",
        "-nsig", str(NSIG),
        "-outDir", run_dir,
    ])

    ndskl = find_single_file(
        os.path.join(run_dir, "*.NDskl"),
        "NDskl"
    )

    print("NDskl:", ndskl)

    # --------------------------------------------------------
    # Step 2: skelconv -> VTK
    # --------------------------------------------------------

    run_command([
        bin_path("skelconv"),
        ndskl,
        "-to", "vtk",
        "-outDir", run_dir,
    ])

    vtk_path = find_single_file(
        os.path.join(run_dir, "*.vtk"),
        "VTK"
    )

    print("VTK:", vtk_path)

    # --------------------------------------------------------
    # Step 3: VTKの長さ条件
    # --------------------------------------------------------

    filaments = read_vtk_polylines(vtk_path)

    skeleton = np.zeros(
        image.shape,
        dtype=np.uint8
    )

    length_records = []

    for fil_id, pts in enumerate(filaments):
        segment_vectors = np.diff(
            pts[:, :2],
            axis=0
        )

        length_px = float(
            np.linalg.norm(
                segment_vectors, axis=1
            ).sum()
        )

        length_pc = length_px * pc_per_pixel

        keep = length_pc >= len_cut_pc

        length_records.append((
            fil_id,
            length_pc,
            int(keep)
        ))

        if not keep:
            continue

        for p1, p2 in zip(pts[:-1], pts[1:]):
            xy1 = tuple(
                np.rint(p1[:2]).astype(int)
            )
            xy2 = tuple(
                np.rint(p2[:2]).astype(int)
            )

            cv2.line(
                skeleton,
                xy1,
                xy2,
                1,
                thickness=1
            )

    # --------------------------------------------------------
    # Step 4: Avマスク境界のノイズ除去
    # --------------------------------------------------------

    valid_mask = np.isfinite(masked)

    valid_eroded = cv2.erode(
        valid_mask.astype(np.uint8),
        np.ones((5, 5), dtype=np.uint8),
        iterations=2
    ).astype(bool)

    skeleton[~valid_eroded] = 0

    # --------------------------------------------------------
    # Step 5: Skeleton FITS保存
    # --------------------------------------------------------

    skeleton_path = os.path.join(
        full_dir,
        f"skeleton_{REGION}_filtered.fits"
    )

    fits.writeto(
        skeleton_path,
        skeleton,
        header,
        overwrite=True
    )

    # --------------------------------------------------------
    # Step 6: 長さCSV保存
    # --------------------------------------------------------

    csv_path = os.path.join(
        full_dir,
        f"filaments_{REGION}_length.csv"
    )

    with open(
        csv_path,
        "w",
        newline="",
        encoding="utf-8"
    ) as f:
        writer = csv.writer(f)

        writer.writerow([
            "vtk_polyline_id",
            "length_pc",
            "keep"
        ])

        writer.writerows(length_records)

    # --------------------------------------------------------
    # Step 7: 通常PNG（表示専用3x3）
    # --------------------------------------------------------

    gray = make_gray_background(image)

    skeleton_vis = binary_dilation(
        skeleton.astype(bool),
        structure=np.ones((3, 3), dtype=bool)
    )

    rgb = np.repeat(
        gray[:, :, None],
        3,
        axis=2
    )

    rgb[skeleton_vis] = (255, 255, 0)

    overlay_path = os.path.join(
        vis_dir,
        f"overlay_{REGION}.png"
    )

    Image.fromarray(rgb).save(overlay_path)

    n_kept = sum(r[2] for r in length_records)

    print("================================")
    print("検出VTKポリライン:", len(filaments))
    print("採用:", n_kept)
    print("Skeleton画素数:", np.count_nonzero(skeleton))
    print("Skeleton FITS:", skeleton_path)
    print("黄色PNG:", overlay_path)
    print("中間ファイル:", run_dir)
    print("DONE (loose)")


if __name__ == "__main__":
    main()

