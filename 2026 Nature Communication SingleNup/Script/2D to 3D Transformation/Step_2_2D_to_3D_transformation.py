# -*- coding: utf-8 -*-
"""
Step 2 - 2D -> 3D transformation
================================

Reads the binning settings produced by Step 1 (``settings/Xbinning_settings.csv``),
re-reads each data file from ``data/``, and reconstructs the radial (3D) density
for every X subregion and bin size via the area-matrix inversion.

Output layout::

    transformation_results/
        result_<timestamp>/
            raw/            <- density divided by X-range width (negatives kept)
            zero_negative/  <- negatives clipped to 0, then divided by max
            min_max/        <- negatives clipped to 0, then min-max scaled to [0,1]
            summary.xlsx    <- chi-square + diagnostics, one sheet per data file
            run_meta.json   <- archived settings + data fingerprint (for dedup)
            warnings.txt

Each run is fingerprinted from the settings and the data-file contents. If a
previous result folder has the same fingerprint, this script reports where that
result is and does nothing, so every result folder is guaranteed to differ.

All tunable values live in config.py. Original author: Wenlan Yu. Modernized 2026.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import math
import sys
from collections import namedtuple
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist
from scipy.stats import mannwhitneyu

from config import (
    DATA_DIR,
    SETTINGS_DIR,
    SETTINGS_FILE,
    RESULTS_DIR,
    R_upperlimit,
    bin_min,
    bin_max,
    bin_delta,
    BIN_DELTA_FLOOR,
    MWU_threshold,
    background_filter,
    background_filter_range,
    save_format,
)
from data_io import warn, info, load_xy, check_integrity

# The three normalization versions written to separate subfolders.
NORM_MODES = ["raw", "zero_negative", "min_max"]

# Decimal places implied by the bin-size floor (0.001 -> 3 decimals). Used both
# for the integer-arithmetic grid and for naming the output columns cleanly.
_BIN_DECIMALS = max(0, -int(round(math.log10(BIN_DELTA_FLOOR))))
_BIN_SCALE = round(1.0 / BIN_DELTA_FLOOR)


def make_bin_sizes(bmin: float, bmax: float, delta: float) -> list[float]:
    """Build the list of R bin sizes, supporting fractional steps.

    Works for integer or float bin_min/bin_max/bin_delta. A bin_delta below
    BIN_DELTA_FLOOR is raised to the floor (with a warning). Everything is snapped
    to the floor's grid using integer arithmetic, so there is no floating-point
    drift; bin_max is included when it lands on the grid.
    """
    if delta < BIN_DELTA_FLOOR:
        warn(f"bin_delta = {delta} is below the floor {BIN_DELTA_FLOOR}; "
             f"using {BIN_DELTA_FLOOR} instead.")
        delta = BIN_DELTA_FLOOR
    lo = round(bmin * _BIN_SCALE)
    hi = round(bmax * _BIN_SCALE)
    step = max(round(delta * _BIN_SCALE), 1)
    return [m / _BIN_SCALE for m in range(lo, hi + 1, step)]


def format_bin(b: float) -> str:
    """Clean string label for a (possibly fractional) bin size: 4.0 -> '4'."""
    return f"{round(b, _BIN_DECIMALS):.{_BIN_DECIMALS}f}".rstrip("0").rstrip(".")


# ---------------------------------------------------------------------------
# X-range separation
# ---------------------------------------------------------------------------
def parse_cutpoint(value: object) -> list[float]:
    """Parse a Cutpoint cell like '[-100, -67, ...]' into a list of floats."""
    text = str(value).replace("[", "").replace("]", "")
    return [float(p) for p in text.split(",") if p.strip() != ""]


def x_separation(range_min: float, range_max: float,
                 separate_points: list[float]) -> list[tuple[float, float]]:
    """Split [range_min, range_max] at every interior cut point.

    Returns consecutive (low, high) subregions. Cut points equal to or outside
    the endpoints are ignored. (Cleaned-up, order-stable replacement for the
    original X_sepration, which mutated the list while iterating.)
    """
    interior = sorted(p for p in separate_points if range_min < p < range_max)
    bounds = [range_min] + interior + [range_max]
    return [(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)]


# ---------------------------------------------------------------------------
# Symmetry check + radial data extraction
# ---------------------------------------------------------------------------
def data_filter(df: pd.DataFrame, x_lo: float, x_hi: float,
                file_label: str, range_label: str, warnings: list[str]) -> pd.Series:
    """Return |Y| for rows whose X falls in [x_lo, x_hi], with a symmetry check.

    A Mann-Whitney U test compares the negative and positive Y populations; if
    they differ more than MWU_threshold allows, a warning is recorded.
    """
    sub = df.loc[(df["X"] >= x_lo) & (df["X"] <= x_hi)]
    y = sub["Y"]
    y_pos = y[y >= 0]
    y_neg = y[y <= 0]
    if len(y_pos) > 0 and len(y_neg) > 0:
        try:
            result = mannwhitneyu(abs(y_neg), y_pos, method="auto")
            if result.pvalue < MWU_threshold:
                warnings.append(
                    f"Caution: Y data for {file_label} range {range_label} is not "
                    f"symmetric over 0 (Mann-Whitney U p = {result.pvalue:.3g})."
                )
        except ValueError:
            pass  # not enough points on one side -> skip the test
    return abs(y)


# ---------------------------------------------------------------------------
# Core area-matrix inversion (numerics preserved from the original)
# ---------------------------------------------------------------------------
# Per-inversion result: the physical 3D density plus reconstruction diagnostics.
InversionResult = namedtuple(
    "InversionResult",
    ["density", "n_compared", "chi_square", "chi_square_dof",
     "chi_square_pvalue", "poisson_deviance"],
)


def build_area_matrix(binsize: int, upperlimit: float) -> "tuple[np.ndarray, int]":
    """Build the (binnumber x binnumber) area matrix for one R bin size."""
    binnumber = math.ceil(upperlimit / binsize) + 1

    def calculate_S(i: int, j: int, binsize: int) -> float:
        return (
            (np.arccos((j - 1) / i) * (i ** 2) * (binsize ** 2)) / 2
            - ((math.sqrt((i ** 2) - ((j - 1) ** 2)) / i) * ((j - 1) / i)
               * (i ** 2) * (binsize ** 2)) / 2
        )

    area_matrix = np.zeros([binnumber, binnumber], dtype="f8")
    for i in range(binnumber):
        for j in range(binnumber):
            if j > i or i == binnumber - 1:
                break
            elif i == j:
                area_matrix[i, j] = calculate_S(i + 1, j + 1, binsize)
            else:
                area_matrix[i, j] = (
                    calculate_S(i + 1, j + 1, binsize)
                    - calculate_S(i, j + 1, binsize)
                    - calculate_S(i + 1, j + 2, binsize)
                    + calculate_S(i, j + 2, binsize)
                )
    area_matrix_sum = np.sum(area_matrix, axis=0)
    area_matrix_bg = np.empty(binnumber)
    area_matrix_bg.fill(binsize * background_filter_range)
    area_matrix_bg[binnumber - 1] = (
        (background_filter_range - (math.ceil(upperlimit / binsize) * binsize))
        * background_filter_range
    )
    area_matrix[binnumber - 1, ] = area_matrix_bg - area_matrix_sum
    return area_matrix, binnumber


def poisson_deviance(observed: np.ndarray, expected: np.ndarray) -> float:
    """EXPERIMENTAL - full Poisson deviance  G2 = 2 * sum( o*ln(o/e) - (o - e) ).

    Provided next to the chi-square for future use. The '-(o - e)' term keeps it
    well defined when the observed and expected totals differ (which they do
    here: clipping negative densities injects mass, so the reconstructed total
    exceeds the original). IMPORTANT: like the chi-square below, this is computed
    in-sample on a deterministic reconstruction, so it is a *relative* diagnostic
    only -- its value does not follow a chi-square distribution and the usual
    p-value interpretation does not apply. Not yet validated; here for comparison.
    """
    o = np.asarray(observed, dtype=float)
    e = np.asarray(expected, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_term = np.where(o > 0, o * np.log(o / e), 0.0)
    return float(2.0 * np.sum(log_term - (o - e)))


def invert_2d_to_3d(binsize: int, upperlimit: float, r_values) -> InversionResult:
    """Invert the 2D radial histogram into a 3D density and score the fit.

    Returns the physical 3D density profile plus Andrew's back-calculation check
    (done before any normalization): clip negative densities to zero, forward-
    project through the area matrix to reconstruct the 2D counts, then compare
    bin-by-bin to the original counts over bins that actually contain data.
    """
    area_matrix, binnumber = build_area_matrix(binsize, upperlimit)
    hist, _ = np.histogram(
        r_values,
        bins=math.ceil(upperlimit / binsize),
        range=(0, binsize * math.ceil(upperlimit / binsize)),
        density=False,
    )
    # zeros (not np.empty) so the background slot is well-defined even when the
    # background filter is off -- fixes a latent garbage-value bug.
    frequency_count = np.zeros(binnumber)
    frequency_count[0:binnumber - 1] = hist
    if background_filter and background_filter_range > R_upperlimit:
        hist_bg, _ = np.histogram(
            r_values,
            bins=1,
            range=(binsize * math.ceil(upperlimit / binsize), background_filter_range),
            density=False,
        )
        frequency_count[binnumber - 1] = hist_bg[0]

    density_full = np.linalg.inv(area_matrix.T) @ frequency_count

    # --- reconstruction check (before normalization) -----------------------
    density_clipped = np.clip(density_full, 0.0, None)       # negatives -> 0
    reconstructed = area_matrix.T @ density_clipped           # back-calc counts

    expected = frequency_count[0:binnumber - 1]               # original 2D counts
    observed = reconstructed[0:binnumber - 1]                 # reconstructed counts
    mask = expected > 0                                        # only bins with data
    e_m = expected[mask]
    o_m = observed[mask]

    if e_m.size > 0:
        chi_square = float(np.sum((o_m - e_m) ** 2 / e_m))
        deviance = poisson_deviance(o_m, e_m)
        dof = max(int(e_m.size) - 1, 1)
        # p-value with df = (bins with data) - 1. This df is a convention; the
        # in-sample reconstruction means it is not a rigorous test (see notes).
        pvalue = float(chi2_dist.sf(chi_square, dof))
    else:
        chi_square = deviance = pvalue = float("nan")
        dof = 0

    return InversionResult(
        density=density_full[0:binnumber - 1],
        n_compared=int(e_m.size),
        chi_square=chi_square,
        chi_square_dof=dof,
        chi_square_pvalue=pvalue,
        poisson_deviance=deviance,
    )


# ---------------------------------------------------------------------------
# Build the per-file result table (subregion-ordered, width-normalized)
# ---------------------------------------------------------------------------
def transform_file(df: pd.DataFrame, x_ranges: list[tuple[float, float]],
                   bin_sizes: list[int], file_label: str,
                   warnings: list[str], summary_rows: list[dict]) -> pd.DataFrame:
    """Reconstruct the radial density for every (subregion, bin size) pair.

    For each pair, the reconstruction chi-square (and experimental Poisson
    deviance) is recorded into ``summary_rows`` before any normalization.
    """
    result_df = pd.DataFrame()
    for (lo, hi) in x_ranges:
        r_values = data_filter(df, lo, hi, file_label, f"[{lo}, {hi}]", warnings)
        width = hi - lo
        for b in bin_sizes:
            inv = invert_2d_to_3d(b, R_upperlimit, r_values)
            prob = inv.density / width  # density per unit X-range (multi_region)
            bin_center = np.arange(b, b * (1 + math.ceil(R_upperlimit / b)), b) - b / 2
            label = format_bin(b)
            prob_df = pd.DataFrame([bin_center, prob]).T
            prob_df.columns = [
                f"[{lo}, {hi}]_Bin{label}_BinCenter",
                f"[{lo}, {hi}]_Bin{label}_Probability",
            ]
            result_df = pd.concat([result_df, prob_df], axis=1)
            summary_rows.append({
                "File": file_label,
                "X_region": f"[{lo}, {hi}]",
                "Region_width": width,
                "Bin_size": round(b, _BIN_DECIMALS),
                "N_localizations": int(len(r_values)),
                "N_bins_compared": inv.n_compared,
                "Chi_square": inv.chi_square,
                "Chi_square_dof": inv.chi_square_dof,
                "Chi_square_pvalue": inv.chi_square_pvalue,
                "Poisson_deviance_experimental": inv.poisson_deviance,
            })
    return result_df


# ---------------------------------------------------------------------------
# The three normalization versions
# ---------------------------------------------------------------------------
def normalize_result(result_df: pd.DataFrame, num_bin: int, mode: str,
                     file_label: str, warnings: list[str]) -> pd.DataFrame:
    """Return a normalized copy of result_df for the requested mode.

    Probability columns sit at odd positions; for a given bin size they recur
    every ``2*num_bin`` columns (one block per subregion). Each bin-size group
    is scaled together so subregions stay comparable (the multi_region scheme).
    """
    df = result_df.copy()
    if mode == "raw":
        return df

    ncols = df.shape[1]
    prob_pos = list(range(1, ncols, 2))
    df.iloc[:, prob_pos] = df.iloc[:, prob_pos].clip(lower=0)  # negatives -> 0

    for i in range(num_bin):
        group_pos = list(range(2 * i + 1, ncols, 2 * num_bin))
        if not group_pos:
            continue
        sub = df.iloc[:, group_pos]
        values = sub.to_numpy(dtype=float)
        if not np.isfinite(values).any():
            continue
        max_v = np.nanmax(values)
        min_v = np.nanmin(values)
        if mode == "zero_negative":
            if np.isfinite(max_v) and max_v != 0:
                df.iloc[:, group_pos] = sub / max_v
            else:
                warnings.append(
                    f"{file_label}: bin-size group {i} is all zero; "
                    f"max-normalization skipped."
                )
        elif mode == "min_max":
            if np.isfinite(max_v) and np.isfinite(min_v) and max_v != min_v:
                df.iloc[:, group_pos] = (sub - min_v) / (max_v - min_v)
            else:
                warnings.append(
                    f"{file_label}: bin-size group {i} is constant; "
                    f"min-max normalization skipped."
                )
    return df


def reorder_by_bin(df: pd.DataFrame, num_bin: int, num_regions: int) -> pd.DataFrame:
    """Regroup columns so the same bin size across subregions sits together."""
    out = pd.DataFrame()
    for x in range(num_bin):
        for y in range(num_regions):
            idx = num_bin * 2 * y + 2 * x
            out = pd.concat([out, df.iloc[:, idx:idx + 2]], axis=1)
    return out


# ---------------------------------------------------------------------------
# Run fingerprint (settings + data contents) for duplicate detection
# ---------------------------------------------------------------------------
def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def build_signature(settings_df: pd.DataFrame, data_dir: Path) -> tuple[dict, str]:
    """Build a JSON-serializable signature and its hash for this run."""
    params = {
        "R_upperlimit": R_upperlimit,
        "bin_min": bin_min,
        "bin_max": bin_max,
        "bin_delta": bin_delta,
        "MWU_threshold": MWU_threshold,
        "background_filter": background_filter,
        "background_filter_range": background_filter_range,
        "save_format": save_format,
    }
    files = []
    for _, row in settings_df.iterrows():
        name = str(row["File"])
        dpath = data_dir / name
        files.append({
            "name": name,
            "x_center": row.get("X Center"),
            "range_min": row.get("Range_min"),
            "range_max": row.get("Range_max"),
            "cutpoint": str(row.get("Cutpoint")),
            "data_sha256": file_sha256(dpath) if dpath.exists() else None,
        })
    signature = {"parameters": params, "files": files}
    canonical = json.dumps(signature, sort_keys=True, default=str)
    sig_hash = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    return signature, sig_hash


def find_duplicate(results_root: Path, sig_hash: str) -> "Path | None":
    """Return the folder of a previous identical run, or None."""
    if not results_root.is_dir():
        return None
    for meta_path in sorted(results_root.glob("*/run_meta.json")):
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001 - ignore unreadable/old meta files
            continue
        if meta.get("signature_hash") == sig_hash:
            return meta_path.parent
    return None


# ---------------------------------------------------------------------------
# Summary output (one Excel sheet per data file)
# ---------------------------------------------------------------------------
def _sheet_name(name: str, used: set) -> str:
    """Make an Excel-safe, unique sheet name (<=31 chars, no forbidden chars)."""
    forbidden = set(r"[]:*?/\\")
    safe = "".join("_" if c in forbidden else c for c in str(name)).strip()
    safe = (safe or "Sheet")[:31]
    base, n = safe, 1
    while safe.lower() in used:
        suffix = f"_{n}"
        safe = base[:31 - len(suffix)] + suffix
        n += 1
    used.add(safe.lower())
    return safe


def write_summary(summary_df: pd.DataFrame, result_dir: Path) -> None:
    """Write the per-(region, bin size) diagnostics, one sheet per data file.

    Falls back to a single ``summary.csv`` if openpyxl is not installed, so the
    pipeline never hard-fails just because the Excel writer is unavailable.
    """
    out = result_dir / "summary.xlsx"
    try:
        with pd.ExcelWriter(out, engine="openpyxl") as writer:
            if summary_df.empty:
                summary_df.to_excel(writer, sheet_name="summary", index=False)
            else:
                used: set = set()
                for fname, group in summary_df.groupby("File", sort=False):
                    sheet = _sheet_name(Path(str(fname)).stem, used)
                    # Keep the File column so the full (untruncated) name is
                    # preserved inside the sheet as well as in its title.
                    group.to_excel(writer, sheet_name=sheet, index=False)
        info(f"Summary written to '{out.name}' (one sheet per data file).")
    except Exception as exc:  # noqa: BLE001 - openpyxl missing or write failure
        try:
            out.unlink(missing_ok=True)
        except OSError:
            pass
        warn(f"Could not write Excel summary ({exc}); writing summary.csv instead.")
        summary_df.to_csv(result_dir / "summary.csv", index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir / DATA_DIR
    settings_path = script_dir / SETTINGS_DIR / SETTINGS_FILE
    results_root = script_dir / RESULTS_DIR

    # --- verify pipeline files/folders + correct run order ---------------------
    if not check_integrity(script_dir, step=2):
        return 1

    if not settings_path.exists():
        warn(f"Settings file not found at '{SETTINGS_DIR}/{SETTINGS_FILE}'. "
             f"Run Step 1 first.")
        return 1
    settings_df = pd.read_csv(settings_path)
    if settings_df.empty:
        warn("Settings file is empty. Nothing to transform.")
        return 1

    bin_sizes = make_bin_sizes(bin_min, bin_max, bin_delta)
    num_bin = len(bin_sizes)
    if num_bin == 0:
        warn("No bin sizes to analyze; check bin_min / bin_max / bin_delta.")
        return 1

    # --- duplicate detection ------------------------------------------------
    signature, sig_hash = build_signature(settings_df, data_dir)
    duplicate = find_duplicate(results_root, sig_hash)
    if duplicate is not None:
        info("The same settings and the same data files have already been "
             "transformed in a previous run.")
        info(f"You can find that result here: {duplicate}")
        return 0

    # --- run the transformation --------------------------------------------
    now = datetime.datetime.now()
    result_dir = results_root / f"result_{now.strftime('%Y_%m_%d_%H_%M_%S')}"

    warnings: list[str] = []
    summary_rows: list[dict] = []
    written = 0
    for _, row in settings_df.iterrows():
        name = str(row["File"])
        dpath = data_dir / name
        if not dpath.exists():
            warnings.append(f"Data file '{name}' not found in '{DATA_DIR}'. Skipped.")
            continue
        df = load_xy(dpath)
        if df is None:
            warnings.append(f"Data file '{name}' could not be loaded. Skipped.")
            continue

        range_min = row["Range_min"]
        range_max = row["Range_max"]
        if range_max < range_min:
            warnings.append(
                f"'{name}': Range_max < Range_min in the settings file. Skipped."
            )
            continue

        cutpoints = parse_cutpoint(row["Cutpoint"])
        x_ranges = x_separation(range_min, range_max, cutpoints)
        result_df = transform_file(df, x_ranges, bin_sizes, name, warnings,
                                   summary_rows)
        num_regions = len(x_ranges)

        out_name = Path(name).stem + ".csv"
        for mode in NORM_MODES:
            version = normalize_result(result_df, num_bin, mode, name, warnings)
            if save_format == "bin":
                version = reorder_by_bin(version, num_bin, num_regions)
            mode_dir = result_dir / mode
            mode_dir.mkdir(parents=True, exist_ok=True)
            version.to_csv(mode_dir / out_name, index=False)
        written += 1
        info(f"'{name}': {num_regions} subregion(s) x {num_bin} bin size(s) "
             f"-> raw / zero_negative / min_max.")

    if written == 0:
        warn("No files were transformed; no result folder written.")
        for w in warnings:
            warn(w)
        if result_dir.exists() and not any(result_dir.iterdir()):
            result_dir.rmdir()
        return 1

    # --- archive settings fingerprint + warnings inside the result folder ---
    meta = {
        "signature_hash": sig_hash,
        "created": now.strftime("%Y-%m-%d %H:%M:%S"),
        "normalization_versions": NORM_MODES,
        "bin_sizes": bin_sizes,
        **signature,
    }
    summary_df = pd.DataFrame(summary_rows, columns=[
        "File", "X_region", "Region_width", "Bin_size", "N_localizations",
        "N_bins_compared", "Chi_square", "Chi_square_dof", "Chi_square_pvalue",
        "Poisson_deviance_experimental",
    ])
    write_summary(summary_df, result_dir)

    (result_dir / "run_meta.json").write_text(
        json.dumps(meta, indent=2, default=str), encoding="utf-8"
    )
    (result_dir / "warnings.txt").write_text(
        ("\n".join(warnings) + "\n") if warnings else "No warnings.\n",
        encoding="utf-8",
    )

    for w in warnings:
        warn(w)
    info(f"Wrote {written} file(s) x 3 versions to "
         f"'{result_dir.relative_to(script_dir)}'.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
