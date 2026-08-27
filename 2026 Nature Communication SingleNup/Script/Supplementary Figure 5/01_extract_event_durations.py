# -*- coding: utf-8 -*-
"""
Step 1 of 3 -- Supplementary Figure S5B/S5C (Nup153C, 2 ms/frame).

Extract unfiltered event durations from the raw single-molecule localization
tables and write them to a single CSV.

Input
-----
A folder of .xlsx localization tables (one file per acquisition, one or more
sheets per file), exported headerless from the localization software. Column
layout (0-indexed):

    col 0   source image identifier (string)
    col 1   frame number (integer)
    col 2   x, nearest pixel (integer)
    col 3   y, nearest pixel (integer)
    col 4+  intensity, sigma, sub-pixel x/y, timestamp, ... (not used here)

Method
------
For each sheet:
  1. Keep only localizations inside a fixed square ROI centred on the pore
     (see ROI_X / ROI_Y below).
  2. Sort by frame number.
  3. An "event" is a run of consecutive frames (frame-to-frame difference
     of exactly 1). Its duration is the number of frames in the run.
  4. Optionally discard runs shorter than MIN_EVENT_FRAMES.

Note that events are defined by temporal contiguity within the ROI only;
there is no spatial linking between frames and no gap-closing. See README.

Output
------
data/combined_group_lengths.csv -- one row per event, duration in FRAMES.
Multiply by the frame time (2.0 ms) to convert to milliseconds.

Usage
-----
    python 01_extract_event_durations.py
"""

import os
import sys

import numpy as np
import pandas as pd

# ==================== SETTINGS ====================
data_dir = "Data_Nup153"                  # folder of .xlsx localization tables
output_csv = os.path.join("data", "combined_group_lengths.csv")

# Square ROI centred on the pore, in pixels, inclusive on both ends.
# The acquisitions are 128 x 128 crops already centred on the target pore,
# so the same window is applied to every file.
ROI_X = (54, 74)
ROI_Y = (54, 74)

FRAME_COL, X_COL, Y_COL = 1, 2, 3

# Minimum run length to keep, in frames. 1 = no filter (the setting used for
# the published S5B/S5C panels). Set to 3 to match the minimum three-
# localization tracking requirement used for the primary 50 ms datasets.
# See "Known limitations" in README.md before changing this.
MIN_EVENT_FRAMES = 1

# The original implementation discarded the first localization of every sheet
# as a side effect of how it handled the leading NaN produced by .diff(). This
# cost 18 single-frame events across the 35 sheets. True reproduces the
# published n = 29720 exactly; False gives the corrected n = 29738, which
# leaves every reported statistic unchanged at the precision quoted in the
# manuscript (delta-AIC 377.8 -> 378.8; component weights and means identical).
LEGACY_DROP_FIRST_LOCALIZATION = True

verbose = False        # per-sheet progress output
# ==================================================


def event_durations(frames):
    """Run lengths of consecutive frame numbers, in frames.

    `frames` is an unsorted 1-D array of frame indices for one sheet. Repeated
    frame numbers (two localizations detected in the same frame inside the ROI)
    break a run, because their frame-to-frame difference is 0, not 1.
    """
    f = np.sort(np.asarray(frames, dtype=float))
    if f.size == 0:
        return np.empty(0, dtype=int)
    # A new run starts wherever the step to the previous frame is not exactly 1.
    starts = np.concatenate(([True], np.diff(f) != 1))
    run_id = np.cumsum(starts)
    _, lengths = np.unique(run_id, return_counts=True)
    lengths = lengths.astype(int)

    # The original code dropped the leading NaN from .diff() positionally,
    # which silently discarded the sheet's first localization whenever it was
    # isolated -- i.e. whenever the first run had length 1. Runs of length > 1
    # were unaffected. Reproduced here so the repository regenerates the
    # published event list exactly.
    if LEGACY_DROP_FIRST_LOCALIZATION and lengths.size and lengths[0] == 1:
        lengths = lengths[1:]

    return lengths


def main():
    folder = os.path.join(os.getcwd(), data_dir)
    if not os.path.isdir(folder):
        sys.exit(f"Data folder not found: {folder}")

    excel_files = sorted(f for f in os.listdir(folder)
                         if f.lower().endswith((".xlsx", ".xls"))
                         and not f.startswith("~$"))
    if not excel_files:
        sys.exit(f"No Excel files found in {folder}")

    combined = []
    n_sheets = 0

    for name in excel_files:
        path = os.path.join(folder, name)
        try:
            book = pd.ExcelFile(path)
        except Exception as err:                       # noqa: BLE001
            print(f"  !! skipping {name}: {err}")
            continue

        for sheet in book.sheet_names:
            df = book.parse(sheet, header=None)
            in_roi = (
                df[X_COL].between(*ROI_X) & df[Y_COL].between(*ROI_Y)
            )
            lengths = event_durations(df.loc[in_roi, FRAME_COL].to_numpy())
            lengths = lengths[lengths >= MIN_EVENT_FRAMES]
            combined.extend(lengths.tolist())
            n_sheets += 1
            if verbose:
                print(f"  {name} :: {sheet} -- {in_roi.sum()} localizations "
                      f"in ROI, {len(lengths)} events")

    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    pd.DataFrame({"Group_Length": combined}).to_csv(output_csv, index=False)

    arr = np.asarray(combined)
    print(f"Files: {len(excel_files)}  |  sheets: {n_sheets}  |  "
          f"events: {arr.size}")
    print(f"Duration in frames -- min {arr.min()}, median {np.median(arr):.0f}, "
          f"max {arr.max()}")
    print(f"Saved: {output_csv}")


if __name__ == "__main__":
    main()
