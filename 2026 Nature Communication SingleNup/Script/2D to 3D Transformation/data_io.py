# -*- coding: utf-8 -*-
"""
Shared data input/output helpers for the 2D -> 3D transformation pipeline.
=========================================================================

Both Step 1 (automatic X binning) and Step 2 (transformation) read the same raw
data files, so the loader lives here in one place. ``load_xy`` accepts Excel and
CSV, with or without a header, and returns a clean two-column ``X`` / ``Y``
DataFrame (or ``None`` after printing an error if the file cannot be used).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from config import (
    X_HEADERS,
    Y_HEADERS,
    EXCEL_EXTS,
    DATA_DIR,
    SETTINGS_DIR,
    SETTINGS_FILE,
    RESULTS_DIR,
)


# ---------------------------------------------------------------------------
# Simple logging helpers (no external dependency)
# ---------------------------------------------------------------------------
def warn(msg: str) -> None:
    print(f"[WARNING] {msg}")


def error(msg: str) -> None:
    print(f"[ERROR]   {msg}")


def info(msg: str) -> None:
    print(f"[INFO]    {msg}")


# ---------------------------------------------------------------------------
# Loading and validating one data file
# ---------------------------------------------------------------------------
def _normalize(name: object) -> str:
    """Lower-case and remove all whitespace from a header cell for matching."""
    return "".join(str(name).split()).lower()


def _looks_like_header(first_row: pd.Series) -> bool:
    """A row is a header if at least one cell cannot be read as a number."""
    for cell in first_row:
        try:
            float(cell)
        except (TypeError, ValueError):
            return True
    return False


def load_xy(path: Path) -> "pd.DataFrame | None":
    """Load a data file into a two-column DataFrame with columns ``X`` and ``Y``.

    Handles Excel and CSV, with or without a header row. Returns ``None`` (after
    printing an error) if the file cannot be interpreted as clean X/Y data.
    """
    path = Path(path)
    ext = path.suffix.lower()
    try:
        if ext in EXCEL_EXTS:
            raw = pd.read_excel(path, header=None)
        else:
            raw = pd.read_csv(path, header=None)
    except Exception as exc:  # noqa: BLE001 - report any read failure to the user
        error(f"Could not read '{path.name}': {exc}")
        return None

    if raw.empty:
        error(f"'{path.name}' is empty.")
        return None

    has_header = _looks_like_header(raw.iloc[0])

    if has_header:
        header = [_normalize(c) for c in raw.iloc[0]]
        x_cols = [i for i, h in enumerate(header) if h in X_HEADERS]
        y_cols = [i for i, h in enumerate(header) if h in Y_HEADERS]
        if len(x_cols) != 1 or len(y_cols) != 1:
            error(
                f"'{path.name}': need exactly one X column and one Y column, but "
                f"found {len(x_cols)} X and {len(y_cols)} Y columns. "
                f"Accepted X headers: X / X(nm) / X (nm) (any case); "
                f"likewise for Y."
            )
            return None
        body = raw.iloc[1:, [x_cols[0], y_cols[0]]].copy()
    else:
        # No header -> the data itself must be exactly two columns (X, Y).
        if raw.shape[1] != 2:
            error(
                f"'{path.name}' has no header and {raw.shape[1]} columns. "
                f"Headerless files must have exactly two columns (X, Y)."
            )
            return None
        body = raw.copy()

    body.columns = ["X", "Y"]
    # Coerce to numeric; anything non-numeric becomes NaN so the row check below
    # can flag it.
    body["X"] = pd.to_numeric(body["X"], errors="coerce")
    body["Y"] = pd.to_numeric(body["Y"], errors="coerce")
    body = body.reset_index(drop=True)

    # Row check: each row must have BOTH numbers or BOTH NaN (a blank separator
    # line). Exactly one NaN means a corrupt row -> report it and reject the file.
    x_na = body["X"].isna()
    y_na = body["Y"].isna()
    one_na = x_na ^ y_na
    if one_na.any():
        # +2: one for the 1-based count, one for the header line (if present).
        offset = 2 if has_header else 1
        bad_rows = [int(i) + offset for i in body.index[one_na]]
        error(
            f"'{path.name}': {len(bad_rows)} row(s) have a value in only one of "
            f"X/Y (one number, one blank). Offending line number(s) in the file: "
            f"{bad_rows}. Please fix the data and re-run."
        )
        return None

    # Drop the all-NaN separator rows; keep only real data.
    clean = body[~(x_na & y_na)].reset_index(drop=True)
    if clean.empty:
        error(f"'{path.name}' contains no numeric data rows.")
        return None

    return clean


# ---------------------------------------------------------------------------
# Pipeline integrity check
# ---------------------------------------------------------------------------
# All four scripts must live together in the same folder for the pipeline to
# work. The check below guards against missing/deleted files and against running
# the steps out of order.
REQUIRED_SCRIPTS = [
    "config.py",
    "data_io.py",
    "Step_1_automaticXbinning.py",
    "Step_2_2D_to_3D_transformation.py",
]


def check_integrity(script_dir, step: int) -> bool:
    """Verify the pipeline's files and folders before a step runs.

    Confirms that all four pipeline scripts sit together in ``script_dir`` and
    that the input ``data`` folder exists. ``step`` is 1 or 2; Step 2 also
    requires the settings file produced by Step 1, which guards against running
    the steps in the wrong order.

    Returns True if it is safe to proceed, or False (after printing exactly what
    is missing) so the caller can abort.
    """
    script_dir = Path(script_dir)
    ok = True
    info("Checking pipeline integrity ...")

    # 1. all four scripts must be present together
    for fname in REQUIRED_SCRIPTS:
        if (script_dir / fname).is_file():
            info(f"  [ok]      {fname}")
        else:
            ok = False
            error(f"  [MISSING] {fname}")
    if not ok:
        error("One or more pipeline scripts are missing. All four must stay "
              "together in the same folder:")
        error("  " + ", ".join(REQUIRED_SCRIPTS))

    # 2. the input data folder must exist
    if (script_dir / DATA_DIR).is_dir():
        info(f"  [ok]      {DATA_DIR}/ (input data folder)")
    else:
        ok = False
        error(f"  [MISSING] {DATA_DIR}/ - create it and put your data files inside.")

    # 3. these are created automatically; just report their status
    for folder, note in ((SETTINGS_DIR, "created by Step 1"),
                         (RESULTS_DIR, "created by Step 2")):
        if (script_dir / folder).is_dir():
            info(f"  [ok]      {folder}/")
        else:
            info(f"  [--]      {folder}/ not present yet ({note}).")

    # 4. Step 2 cannot run until Step 1 has produced the settings file
    if step == 2:
        settings_file = script_dir / SETTINGS_DIR / SETTINGS_FILE
        if settings_file.is_file():
            info(f"  [ok]      {SETTINGS_DIR}/{SETTINGS_FILE}")
        else:
            ok = False
            error(f"  [MISSING] {SETTINGS_DIR}/{SETTINGS_FILE}")
            error("Step 2 needs the settings file produced by Step 1. "
                  "Please run Step_1_automaticXbinning.py first.")

    if ok:
        info("Integrity check passed.")
    else:
        error("Integrity check failed; aborting. See the messages above.")
    return ok
