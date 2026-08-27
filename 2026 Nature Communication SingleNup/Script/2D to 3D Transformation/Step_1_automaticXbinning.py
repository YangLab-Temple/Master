# -*- coding: utf-8 -*-
"""
Step 1 - Automatic X binning
============================

Reads every data file (Excel / CSV, with or without a header) from a ``data``
subfolder, automatically bins each file along the X dimension so that every bin
holds at least ``MIN_COUNT`` localizations, and writes the resulting cut points
to an editable settings CSV in a ``settings`` subfolder.

The settings file produced here is the input to Step 2 (the actual 2D -> 3D
transformation). It is plain CSV so it can be opened in Excel and adjusted by
hand (e.g. tweak an X center or a range) before running Step 2.

Folder layout (the main folder stays clean - only this script lives there)::

    <main folder>/
        Step_1_automaticXbinning.py
        data/        <- put your .csv / .xlsx / .xls files here
        settings/    <- this script writes (and re-reads) the settings CSV here

Original author: Wenlan Yu. Modernized June 2026.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# All tunable settings live in the shared config module so every pipeline step
# reads the same values. Edit config.py to change behaviour.
from config import (
    DATA_DIR,
    SETTINGS_DIR,
    SETTINGS_FILE,
    MIN_COUNT,
    MAX_BIN_WIDTH,
    MAX_BINS_PER_SIDE,
    EXCEL_EXTS,
    CSV_EXTS,
)

# The data loader and logging helpers are shared with Step 2 via data_io.
from data_io import warn, info, load_xy, check_integrity


# ---------------------------------------------------------------------------
# Automatic X binning
# ---------------------------------------------------------------------------
def auto_x_bin(
    x_values: np.ndarray,
    x_center: float,
    min_count: int = MIN_COUNT,
    max_bin_width: int = MAX_BIN_WIDTH,
    max_bins_per_side: int = MAX_BINS_PER_SIDE,
) -> tuple[list[float], list[int]] | None:
    """Build bin edges along X, growing outward from ``x_center``.

    Starting from a symmetric central bin that contains at least ``min_count``
    points, bins are appended to the left and right. Each new bin grows in width
    until it reaches ``min_count`` (or is abandoned at ``max_bin_width``).

    Returns ``(edges, counts)`` where ``edges`` has one more element than
    ``counts``, or ``None`` if even the widest central bin cannot reach
    ``min_count``.
    """
    x = np.asarray(x_values, dtype=float)
    x = x[~np.isnan(x)]

    def count(lo: float, hi: float) -> int:
        return int(np.count_nonzero((x >= lo) & (x < hi)))

    # 1. Central bin: smallest symmetric bin around the center with >= min_count.
    edges: list[float] | None = None
    counts: list[int] | None = None
    for half in range(1, max_bin_width // 2 + 1):
        lo, hi = x_center - half, x_center + half
        c = count(lo, hi)
        if c >= min_count:
            edges, counts = [lo, hi], [c]
            break
    if edges is None or counts is None:
        return None

    # 2. Extend left.
    for _ in range(max_bins_per_side):
        left = edges[0]
        for width in range(1, max_bin_width + 1):
            lo = left - width
            c = count(lo, left)
            if c >= min_count:
                edges.insert(0, lo)
                counts.insert(0, c)
                break
        else:
            break  # could not fill another bin on the left -> stop

    # 3. Extend right (mirror of the left search).
    for _ in range(max_bins_per_side):
        right = edges[-1]
        for width in range(1, max_bin_width + 1):
            hi = right + width
            c = count(right, hi)
            if c >= min_count:
                edges.append(hi)
                counts.append(c)
                break
        else:
            break

    return edges, counts


# ---------------------------------------------------------------------------
# Settings reuse
# ---------------------------------------------------------------------------
def load_existing_centers(settings_path: Path) -> dict[str, float]:
    """Read previously-saved X centers so manual edits survive a re-run.

    Returns a mapping of file name -> X center for every row whose ``X Center``
    is a usable number. Missing file / non-numeric values are skipped.
    """
    if not settings_path.exists():
        return {}
    try:
        prev = pd.read_csv(settings_path)
    except Exception as exc:  # noqa: BLE001
        warn(f"Could not read existing settings '{settings_path.name}': {exc}")
        return {}
    if "File" not in prev.columns or "X Center" not in prev.columns:
        return {}
    centers: dict[str, float] = {}
    for _, row in prev.iterrows():
        value = pd.to_numeric(row["X Center"], errors="coerce")
        if pd.notna(value):
            centers[str(row["File"])] = float(value)
    return centers


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir / DATA_DIR
    settings_dir = script_dir / SETTINGS_DIR
    settings_path = settings_dir / SETTINGS_FILE

    # --- verify the pipeline files/folders are intact --------------------------
    if not check_integrity(script_dir, step=1):
        return 1

    # --- locate the data folder ------------------------------------------------
    if not data_dir.is_dir():
        warn(
            f"No '{DATA_DIR}' subfolder found next to this script "
            f"(expected at: {data_dir}). Create it and put your CSV/Excel files "
            f"inside, then run again."
        )
        return 1

    data_files = sorted(
        p for p in data_dir.iterdir()
        if p.is_file() and p.suffix.lower() in (EXCEL_EXTS | CSV_EXTS)
    )
    if not data_files:
        if not any(data_dir.iterdir()):
            warn(f"The '{DATA_DIR}' folder is empty. Add CSV or Excel files and re-run.")
        else:
            warn(
                f"The '{DATA_DIR}' folder contains no CSV or Excel files "
                f"(looking for: {sorted(EXCEL_EXTS | CSV_EXTS)})."
            )
        return 1

    info(f"Found {len(data_files)} data file(s) in '{DATA_DIR}'.")

    # --- reuse manually-edited X centers if a settings file already exists ----
    existing_centers = load_existing_centers(settings_path)

    records = []
    for path in data_files:
        df = load_xy(path)
        if df is None:
            continue  # error already reported; skip this file

        x = df["X"].to_numpy(dtype=float)

        # X center: use a saved value if present, otherwise the median of X.
        # Rounded to an integer so that all bin edges stay integers (the bin
        # widths are integers), matching the original design.
        if path.name in existing_centers:
            x_center = int(round(existing_centers[path.name]))
            center_source = "settings"
        else:
            x_center = int(round(np.median(x)))
            center_source = "median"

        result = auto_x_bin(x, x_center)
        if result is None:
            warn(
                f"'{path.name}': no bin around X center {x_center:g} reached "
                f"{MIN_COUNT} points (checked widths up to {MAX_BIN_WIDTH} nm). "
                f"Skipped."
            )
            continue

        edges, counts = result
        records.append(
            {
                "File": path.name,
                "X Center": x_center,
                "Range_min": edges[0],
                "Range_max": edges[-1],
                "Localization Count": counts,
                "Cutpoint": edges,
            }
        )
        info(
            f"'{path.name}': center={x_center:g} ({center_source}), "
            f"{len(counts)} bins, range [{edges[0]:g}, {edges[-1]:g}]."
        )

    if not records:
        warn("No files could be binned successfully. Nothing written.")
        return 1

    # --- write the settings CSV ----------------------------------------------
    settings_dir.mkdir(exist_ok=True)
    out_df = pd.DataFrame(
        records,
        columns=["File", "X Center", "Range_min", "Range_max",
                 "Localization Count", "Cutpoint"],
    )
    out_df.to_csv(settings_path, index=False)
    info(f"Settings written to '{SETTINGS_DIR}/{SETTINGS_FILE}' ({len(out_df)} file(s)).")
    info("You can edit that CSV (e.g. X Center or ranges) before running Step 2.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
