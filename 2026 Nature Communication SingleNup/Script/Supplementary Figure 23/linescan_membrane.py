"""
linescan_membrane.py
=====================
Line-scan / membrane extraction for NE (nuclear-envelope) text images, plus an
alignment-error tilt-angle distribution for the SingleNup tilt simulation.

HOW TO USE (Spyder)
-------------------
1. Put your text-image CSVs in a subfolder called  data/  next to this script.
2. Edit the SETTINGS block just below (paths, fit, camera, plot style).
3. Press Run. Results are written to the  output/  folder.

Everything tunable lives in SETTINGS. There is no command line -- editing those
variables and hitting Run is all you need.

WHAT IT DOES
------------
For each CSV (columns = image X, rows = image Y; NE membrane runs ~vertically):
  1. Find the membrane peak X per row with sub-pixel precision (local Gaussian
     fit) using a tracking window so the trace stays on the membrane.
  2. Reject weak/ambiguous rows (SNR + prominence).
  3. Robustly fit  x = a*y^2 + b*y + c  with iterative sigma-clipping. The
     coordinate swap you described is implicit: we fit image-X (the peak) vs
     image-Y (the row), i.e. polyfit(y, x) -- no manual flipping.
  4. Overlay the fitted membrane on the image; report tangent angle vs vertical
     (theta = arctan(2*a*y + b)).

ALIGNMENT-ERROR TILT DISTRIBUTION
---------------------------------
The most-vertical part of the membrane (parabola APEX, dx/dy = 0) is aligned to
the laser crosshair, but alignment is imperfect (~500 nm = +/-3 px on the 160 nm
Andor EMCCD, +/-2 px on the 240 nm Cascade). The NPC actually imaged sits a
distance delta off the apex, where its tilt is theta = arctan(2*a*delta). Each
cell has a different curvature a, so pooling theta across all cells gives the
realistic EXPERIMENTAL tilt distribution to feed the Monte Carlo simulation.

  - The tilt ANGLE is scale-invariant (X and Y share a pixel size, it's a ratio),
    so px->nm does not change it. Pixel size enters ONLY via the offset:
    delta_px = ALIGN_OFFSET_NM / PIXEL_SIZE_NM.
  - This is the IN-PLANE membrane tilt at the imaged spot. Confirm with your
    mentor that this is the tilt axis the simulation expects (vs out-of-plane
    tilt toward the optical/Z axis).

OUTPUTS (in output/)
  <stem>_membrane.csv          y, x_peak, prominence, used, x_fit, tangent (+nm)
  <stem>_overlay.png           image with fitted membrane curve
  <stem>_angle.png             tangent angle vs y
  cells_summary.csv            per-cell a, apex, curvature radius, max tilt
  pooled_tilt_distribution.csv all sampled tilt angles across cells
  pooled_tilt_hist.png         histogram of the pooled distribution
"""

from __future__ import annotations

import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")          # safe for batch saving; comment out to show in Spyder
import matplotlib.pyplot as plt


# ========================================================================== #
#                              SETTINGS                                       #
#         Edit these, then press Run. This is the only thing to change.       #
# ========================================================================== #

# --- folders (relative to this script) ---------------------------------------
DATA_DIR = "data"          # where your input CSVs live
OUT_DIR  = "output"        # where results are written (kept separate from data)

# --- membrane fit ------------------------------------------------------------
SIDE          = "left"     # "left", "right", or "auto" -- which arc the membrane is on
TRACK_WIN     = 4          # peak-search half-window per row (px)
FIT_PTS       = 2          # +/- points around local max used for sub-pixel fit
MIN_SNR       = 2.0        # reject rows with peak SNR below this
MIN_PROM_FRAC = 0.15       # reject rows below this fraction of the max prominence
POLY_DEGREE   = 2          # polynomial degree (2 = parabola)
CLIP_SIGMA    = 2.5        # sigma-clip threshold for outlier rejection
EDGE_TRIM     = 0          # drop this many rows at top AND bottom before fitting

# --- camera / alignment error ------------------------------------------------
PIXEL_SIZE_NM   = 160.0    # Andor EMCCD = 160, Cascade = 240
ALIGN_OFFSET_NM = 500.0    # alignment uncertainty in nm (~500). Converted to px.
ALIGN_OFFSET_PX = None     # set a number to override the offset directly in px
OFFSET_MODEL    = "uniform"  # "uniform" (conservative, all offset along arc) or
                             # "disk" (2-D offset projected onto the membrane)
N_SAMPLES       = 2000     # Monte Carlo samples per cell for the tilt distribution
SEED            = 0        # RNG seed (reproducible)

# --- plot style --------------------------------------------------------------
SHOW_TITLE      = False
SHOW_LEGEND     = True
TITLE_FONTSIZE  = 20
TITLE_BOLD      = False
LABEL_FONTSIZE  = 18
LABEL_BOLD      = False
TICK_FONTSIZE   = 18
LEGEND_FONTSIZE = 15
FIT_COLOR       = "red"    # membrane line colour on the overlay
FIT_LW          = 2      # membrane line width
HIST_COLOR      = "#4477aa"
PLOT_UNITS      = "px"     # "px" or "nm" for overlay / angle axes
DPI             = 600

# ========================================================================== #
#                          END OF SETTINGS                                    #
# ========================================================================== #


def _build_style():
    return dict(
        show_title=SHOW_TITLE, show_legend=SHOW_LEGEND,
        title_fontsize=TITLE_FONTSIZE, title_bold=TITLE_BOLD,
        label_fontsize=LABEL_FONTSIZE, label_bold=LABEL_BOLD,
        tick_fontsize=TICK_FONTSIZE, legend_fontsize=LEGEND_FONTSIZE,
        fit_color=FIT_COLOR, fit_lw=FIT_LW, hist_color=HIST_COLOR,
        units=PLOT_UNITS, dpi=DPI,
    )


def _weight(flag):
    return "bold" if flag else "normal"


def style_axis(ax, style, title=None, xlabel=None, ylabel=None):
    """Apply font sizes / weights / title toggle to one axis."""
    if title is not None and style["show_title"]:
        ax.set_title(title, fontsize=style["title_fontsize"],
                     fontweight=_weight(style["title_bold"]))
    else:
        ax.set_title("")
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=style["label_fontsize"],
                      fontweight=_weight(style["label_bold"]))
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=style["label_fontsize"],
                      fontweight=_weight(style["label_bold"]))
    ax.tick_params(axis="both", labelsize=style["tick_fontsize"])


# --------------------------------------------------------------------------- #
#  Sub-pixel peak in a single row
# --------------------------------------------------------------------------- #
def gaussian_subpixel(x_idx, y_val):
    """Sub-pixel peak via parabola fit to log(intensity) (= Gaussian centre)."""
    y_val = np.clip(y_val, 1e-6, None)
    ly = np.log(y_val)
    try:
        A, B, _ = np.polyfit(x_idx, ly, 2)
        if A >= 0:
            return float(x_idx[np.argmax(y_val)])
        vtx = -B / (2.0 * A)
        if vtx < x_idx.min() or vtx > x_idx.max():
            return float(x_idx[np.argmax(y_val)])
        return float(vtx)
    except (np.linalg.LinAlgError, ValueError):
        return float(x_idx[np.argmax(y_val)])


def find_row_peak(row, center, half_win, fit_pts=2):
    """Membrane peak in one row within a tracking window."""
    n = row.size
    lo = int(max(0, np.floor(center - half_win)))
    hi = int(min(n - 1, np.ceil(center + half_win)))
    if hi <= lo:
        return np.nan, 0.0, 0.0
    seg = row[lo:hi + 1].astype(float)
    seg_idx = np.arange(lo, hi + 1)
    baseline = np.median(seg)
    k = int(np.argmax(seg))
    prominence = seg[k] - baseline
    resid = seg - baseline
    sigma = 1.4826 * np.median(np.abs(resid - np.median(resid))) + 1e-9
    snr = prominence / sigma
    g_lo = max(0, k - fit_pts)
    g_hi = min(seg.size - 1, k + fit_pts)
    peak_x = gaussian_subpixel(seg_idx[g_lo:g_hi + 1], seg[g_lo:g_hi + 1])
    return peak_x, float(prominence), float(snr)


# --------------------------------------------------------------------------- #
#  Trace membrane across rows
# --------------------------------------------------------------------------- #
def trace_membrane(img, half_win=4, fit_pts=2, side="auto"):
    n_rows, n_cols = img.shape
    if side == "left":
        cols = slice(0, n_cols // 2 + 2)
    elif side == "right":
        cols = slice(n_cols // 2 - 2, n_cols)
    else:
        cols = slice(0, n_cols)
    sub = img[:, cols]
    anchor_row = int(np.argmax(sub.max(axis=1)))
    anchor_x = float(np.argmax(img[anchor_row]) if side == "auto"
                     else (np.argmax(sub[anchor_row]) + (cols.start or 0)))

    peak_x = np.full(n_rows, np.nan)
    prom = np.zeros(n_rows)
    snr = np.zeros(n_rows)

    center = anchor_x
    for r in range(anchor_row, n_rows):
        px, pr, sn = find_row_peak(img[r], center, half_win, fit_pts)
        peak_x[r], prom[r], snr[r] = px, pr, sn
        if np.isfinite(px):
            center = px
    center = anchor_x
    for r in range(anchor_row - 1, -1, -1):
        px, pr, sn = find_row_peak(img[r], center, half_win, fit_pts)
        peak_x[r], prom[r], snr[r] = px, pr, sn
        if np.isfinite(px):
            center = px
    return peak_x, prom, snr, anchor_row


# --------------------------------------------------------------------------- #
#  Robust polynomial fit + tilt angle
# --------------------------------------------------------------------------- #
def robust_polyfit(y, x, weights=None, degree=2, clip_sigma=2.5, n_iter=5):
    y = np.asarray(y, float); x = np.asarray(x, float)
    w = np.ones_like(x) if weights is None else np.asarray(weights, float)
    mask = np.isfinite(x) & np.isfinite(y) & (w > 0)
    coeffs = None
    for _ in range(n_iter):
        if mask.sum() <= degree + 1:
            break
        coeffs = np.polyfit(y[mask], x[mask], degree, w=np.sqrt(w[mask]))
        resid = x - np.polyval(coeffs, y)
        s = np.std(resid[mask])
        if s < 1e-9:
            break
        new_mask = mask & (np.abs(resid) <= clip_sigma * s)
        if new_mask.sum() == mask.sum():
            mask = new_mask; break
        mask = new_mask
    if coeffs is None:
        coeffs = np.polyfit(y[mask], x[mask], degree)
    return np.poly1d(coeffs), mask


def tangent_angle_deg(coeffs, y):
    return np.degrees(np.arctan(np.polyder(coeffs)(y)))


# --------------------------------------------------------------------------- #
#  IO
# --------------------------------------------------------------------------- #
def load_image(path):
    with open(path) as fh:
        first = fh.readline()
    delim = "," if first.count(",") >= 1 else None
    return np.loadtxt(path, delimiter=delim)


def find_csv_files(data_dir):
    return sorted(glob.glob(os.path.join(data_dir, "*.csv")) +
                  glob.glob(os.path.join(data_dir, "*.CSV")))


# --------------------------------------------------------------------------- #
#  Per-file driver
# --------------------------------------------------------------------------- #
def run_one(path, out_dir, style, pixel_nm):
    img = load_image(path)
    n_rows, n_cols = img.shape
    rows = np.arange(n_rows)

    peak_x, prom, snr, anchor = trace_membrane(img, TRACK_WIN, FIT_PTS, SIDE)

    prom_thresh = MIN_PROM_FRAC * np.nanmax(prom)
    used = (snr >= MIN_SNR) & (prom >= prom_thresh) & np.isfinite(peak_x)
    if EDGE_TRIM > 0:
        used[:EDGE_TRIM] = False
        used[n_rows - EDGE_TRIM:] = False

    coeffs, _ = robust_polyfit(rows[used], peak_x[used], weights=prom[used],
                               degree=POLY_DEGREE, clip_sigma=CLIP_SIGMA)
    x_fit = coeffs(rows)
    theta = tangent_angle_deg(coeffs, rows)

    a = float(coeffs.c[0]) if POLY_DEGREE >= 2 else 0.0    # quadratic coeff (1/px)
    b = float(coeffs.c[1]) if POLY_DEGREE >= 1 else 0.0
    apex_y = (-b / (2.0 * a)) if a != 0 else np.nan
    radius_px = (1.0 / abs(2.0 * a)) if a != 0 else np.inf
    radius_nm = radius_px * pixel_nm

    stem = os.path.splitext(os.path.basename(path))[0]
    rmse_px = float(np.sqrt(np.mean((peak_x[used] - coeffs(rows[used])) ** 2)))

    print(f"\n=== {os.path.basename(path)} ===")
    print(f"  {n_rows}x{n_cols} px, anchor row {anchor}, rows used "
          f"{int(used.sum())}/{n_rows}")
    print(f"  poly: a={a:+.5f}  b={b:+.5f}  c={coeffs.c[-1]:+.3f}  (x in px)")
    print(f"  apex (vertical) at y={apex_y:.2f} px;  apex radius of curvature="
          f"{radius_px:.1f} px = {radius_nm:.0f} nm")
    print(f"  fit RMSE {rmse_px:.3f} px = {rmse_px*pixel_nm:.1f} nm")
    print(f"  tangent over crop (deg): mean={np.mean(theta):+.2f} "
          f"std={np.std(theta):.2f} range=[{theta.min():+.2f},{theta.max():+.2f}]")

    table = np.column_stack([rows, peak_x, prom, used.astype(int), x_fit,
                             theta, rows * pixel_nm, x_fit * pixel_nm])
    header = ("y_px,x_peak_px,prominence,used,x_fit_px,tangent_deg,"
              "y_nm,x_fit_nm")
    np.savetxt(os.path.join(out_dir, f"{stem}_membrane.csv"), table,
               delimiter=",", header=header, comments="", fmt="%.4f")

    scale = pixel_nm if style["units"] == "nm" else 1.0
    ulab = "nm" if style["units"] == "nm" else "px"
    yy = np.linspace(0, n_rows - 1, 300)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(img, cmap="gray", origin="upper", interpolation="nearest",
              extent=[0, n_cols * scale, n_rows * scale, 0])
    ax.plot(coeffs(yy) * scale, yy * scale, "-",
            color=style["fit_color"], lw=style["fit_lw"], label="membrane fit")
    style_axis(ax, style, title=os.path.basename(path),
               xlabel=f"X ({ulab})", ylabel=f"Y ({ulab})")
    if style["show_legend"]:
        ax.legend(fontsize=style["legend_fontsize"], loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"{stem}_overlay.png"), dpi=style["dpi"])
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(rows * scale, theta, "-o", ms=3, color="darkorange")
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    style_axis(ax, style, title=f"{os.path.basename(path)} - tilt",
               xlabel=f"Y ({ulab})", ylabel="tangent vs vertical (deg)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"{stem}_angle.png"), dpi=style["dpi"])
    plt.close(fig)

    return dict(stem=stem, a=a, b=b, apex_y=apex_y, radius_nm=radius_nm,
                rmse_px=rmse_px, theta=theta, used=used)


# --------------------------------------------------------------------------- #
#  Alignment-error tilt distribution (pooled across cells)
# --------------------------------------------------------------------------- #
def sample_tilt_for_cell(a, offset_px, n_samples, model, rng):
    """Tilt angles (deg) when the apex is mis-aligned by up to offset_px.
    theta = arctan(2*a*delta), delta = along-membrane offset."""
    if model == "disk":
        r = offset_px * np.sqrt(rng.uniform(0, 1, n_samples))
        phi = rng.uniform(0, 2 * np.pi, n_samples)
        delta = r * np.sin(phi)            # component along ~vertical membrane
    else:                                  # "uniform" -- conservative
        delta = rng.uniform(-offset_px, offset_px, n_samples)
    return np.degrees(np.arctan(2.0 * a * delta))


def alignment_tilt_analysis(results, out_dir, style, pixel_nm, offset_px,
                            n_samples, model, seed=0):
    rng = np.random.default_rng(seed)
    pooled, rows_summary = [], []
    for r in results:
        ang = sample_tilt_for_cell(r["a"], offset_px, n_samples, model, rng)
        pooled.append(ang)
        max_tilt = np.degrees(np.arctan(abs(2.0 * r["a"]) * offset_px))
        rows_summary.append((r["stem"], r["a"], r["apex_y"], r["radius_nm"],
                             r["rmse_px"] * pixel_nm, max_tilt))
    pooled = np.concatenate(pooled)

    with open(os.path.join(out_dir, "cells_summary.csv"), "w") as fh:
        fh.write("cell,a_per_px,apex_y_px,apex_radius_nm,rmse_nm,"
                 "max_tilt_at_offset_deg\n")
        for s in rows_summary:
            fh.write(f"{s[0]},{s[1]:.6f},{s[2]:.3f},{s[3]:.1f},{s[4]:.1f},"
                     f"{s[5]:.3f}\n")

    np.savetxt(os.path.join(out_dir, "pooled_tilt_distribution.csv"), pooled,
               delimiter=",", header="tilt_deg", comments="", fmt="%.4f")

    print("\n=========== POOLED alignment-error tilt distribution ===========")
    print(f"  cells: {len(results)}   offset: +/-{offset_px:.2f} px "
          f"(~{offset_px*pixel_nm:.0f} nm, {pixel_nm:.0f} nm/px, model={model})")
    print(f"  samples: {pooled.size}")
    print(f"  mean |tilt|      = {np.mean(np.abs(pooled)):.2f} deg")
    print(f"  std              = {np.std(pooled):.2f} deg")
    print(f"  95th pct |tilt|  = {np.percentile(np.abs(pooled), 95):.2f} deg")
    print(f"  max |tilt|       = {np.max(np.abs(pooled)):.2f} deg")

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.hist(pooled, bins=40, color=style["hist_color"], edgecolor="white")
    ax.axvline(0, color="gray", lw=0.8, ls="--")
    style_axis(ax, style,
               title=f"Pooled tilt (+/-{offset_px:.1f}px alignment)",
               xlabel="tilt angle vs vertical (deg)", ylabel="count")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "pooled_tilt_hist.png"), dpi=style["dpi"])
    plt.close(fig)
    return pooled


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #
def main():
    here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() \
        else os.getcwd()
    data_dir = os.path.join(here, DATA_DIR)
    out_dir = os.path.join(here, OUT_DIR)

    if not os.path.isdir(data_dir):
        raise SystemExit(f"data folder not found: {data_dir}\n"
                         f"Create a '{DATA_DIR}' subfolder and put your CSVs there.")
    os.makedirs(out_dir, exist_ok=True)
    files = find_csv_files(data_dir)
    if not files:
        raise SystemExit(f"no CSV files found in {data_dir}")

    style = _build_style()
    offset_px = (ALIGN_OFFSET_PX if ALIGN_OFFSET_PX is not None
                 else ALIGN_OFFSET_NM / PIXEL_SIZE_NM)

    print(f"Found {len(files)} CSV(s) in {data_dir}")
    print(f"Writing outputs to {out_dir}")
    print(f"Pixel size {PIXEL_SIZE_NM:.0f} nm  ->  alignment offset "
          f"+/-{offset_px:.2f} px (~{offset_px*PIXEL_SIZE_NM:.0f} nm)")

    results = []
    for i, f in enumerate(files, 1):
        print(f"\n[{i}/{len(files)}] {os.path.basename(f)} ...")
        try:
            results.append(run_one(f, out_dir, style, PIXEL_SIZE_NM))
        except Exception as e:                       # keep the batch going
            print(f"  !! skipped ({e})")

    if results:
        alignment_tilt_analysis(results, out_dir, style, PIXEL_SIZE_NM,
                                offset_px, N_SAMPLES, OFFSET_MODEL, seed=SEED)
        print(f"\nDone. See {out_dir}")


if __name__ == "__main__":
    main()
