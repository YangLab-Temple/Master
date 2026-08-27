# -*- coding: utf-8 -*-
"""
Bimodal-recovery / resampling validation for the 2D -> 3D transformation
========================================================================

Why this exists
---------------
Reviewer concern (round 2): the published Monte Carlo validates only the
recovery of a *single* peak radial position, not the *full shape* of the 3D
probability distribution. The open/closed-state claims for Nup214C and Nup153C
rest on resolving TWO distinct radial peaks (one central, one peripheral) within
a single population, so what must be validated is that the analysis reliably
recovers a two-component structure -- and, just as important, does NOT invent a
second peak when the truth is one component.

What this script does
---------------------
It reproduces the exact estimator used in the paper on a KNOWN ground-truth
mixture:

  * Build a "region" pool of localizations: N1 peripheral points on a ring of
    radius R (radial coordinate ~ R) and N2 central points (radial coordinate
    ~ 0). Each is broadened by localization precision. The ratio N1:N2 sets the
    relative weight of the two states.
  * Randomly subsample `n_sample` points (without replacement) from the pool,
    run the real 2D->3D transform, and read the SINGLE peak position (argmax of
    the deprojected density, with parabolic sub-bin refinement).
  * Repeat the subsample many times -> a histogram of recovered peak positions.
    Because a random subsample is, by chance, dominated by either the central or
    the peripheral component, the single recovered peak lands near 0 or near R,
    and the histogram is bimodal -- exactly the behaviour reported in the paper.
  * Regenerate the pool `n_pool_reps` times to confirm the result is not an
    artefact of one particular random pool.

Negative controls (essential for the reviewer)
-----------------------------------------------
The same procedure is run on two UNIMODAL ground truths -- all-peripheral and
all-central (same total count). If those histograms stay unimodal while the
mixture histogram is bimodal, then bimodality is diagnostic of a genuine
two-component truth and is NOT manufactured by the resampling procedure. This is
the control that answers the reviewer's strongest possible objection.

Faithfulness
------------
The transform under test is NOT reimplemented -- `invert_2d_to_3d` is imported
from `Step_2_2D_to_3D_transformation.py`, so simulated data passes through the
exact code used on real data.

Outputs (in ./bimodal_resampling_output/)
------------------------------------------
  peaks_raw.csv     one row per subsample: condition, pool_rep, peak_position
  summary.csv       per condition: mode positions, mode weights, input ratio,
                    fraction of peaks in the central vs peripheral band
  histograms.png    overlaid peak-position histograms: mixture vs controls
  per_pool.png      mixture histogram split by pool_rep (RNG-stability check)

Transformation: Andrew Ruba / Wenlan Yu.  Resampling validation harness: 2026.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Progress bar: use tqdm if available, else a lightweight stderr fallback so a
# long batch run always shows it is moving and not stuck in an unknown place.
try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - fallback when tqdm is not installed
    import sys, time

    class tqdm:  # minimal drop-in: supports total, desc, update(), close()
        def __init__(self, total=None, desc="", **_):
            self.total, self.desc, self.n = total, desc, 0
            self._t0 = time.time()
            self._last = -1
        def update(self, k=1):
            self.n += k
            if not self.total:
                return
            pct = int(100 * self.n / self.total)
            if pct != self._last:
                self._last = pct
                el = time.time() - self._t0
                rate = self.n / el if el > 0 else 0
                eta = (self.total - self.n) / rate if rate > 0 else 0
                sys.stderr.write(
                    f"\r{self.desc} {pct:3d}%  ({self.n}/{self.total}, "
                    f"{el:5.1f}s elapsed, ~{eta:5.1f}s left)   ")
                sys.stderr.flush()
        def close(self):
            sys.stderr.write("\n"); sys.stderr.flush()

# The transformation under test -- imported, not reimplemented.
from Step_2_2D_to_3D_transformation import invert_2d_to_3d


# ===========================================================================
# CONFIGURATION  (edit these)
# ===========================================================================
CFG = dict(
    R=30,            # peripheral ring radius (nm).
    R_central=0.0,     # central-component radial position (nm). 0 = on axis.
    precision=7.0,     # localization precision, sigma (nm).
    N1=170,            # peripheral points in the pool.
    N2=30,            # central points in the pool.  (pool size = N1 + N2)
    n_sample=100,      # points drawn per subsample (without replacement).
    n_resample=10000,  # subsamples per pool -> one histogram. Use ~30 to test.
    n_pool_reps=10,    # independent pool regenerations. Use ~2 to test.
    bin_size=10.0,      # R bin size for the transform (nm). Finer than 10 so the
                       # argmax peak position is not coarsely quantized.
    R_upper=150.0,     # radial upper limit passed to the transform (nm).
    refine_peak=True,  # parabolic sub-bin interpolation of the argmax.
    run_controls=True, # also run all-peripheral and all-central unimodal controls.
    seed=20260624,
)


# ===========================================================================
# Data generation
# ===========================================================================
def make_pool(rng, R, R_central, precision, n_periph, n_central):
    """Return |Y| for a pooled 'region': n_periph points projected from a ring of
    radius R, plus n_central points at radial R_central. Each point gets a true
    radial coordinate, its Y-projection, then Gaussian detection noise. Z dropped.
    """
    parts = []
    if n_periph > 0:
        theta = rng.uniform(0.0, 2.0 * np.pi, n_periph)
        y_periph = R * np.cos(theta) + rng.normal(0.0, precision, n_periph)
        parts.append(y_periph)
    if n_central > 0:
        # central component as a (near-)point source at radius R_central; for
        # R_central = 0 this is just precision-broadened noise about the axis.
        theta = rng.uniform(0.0, 2.0 * np.pi, n_central)
        y_central = R_central * np.cos(theta) + rng.normal(0.0, precision, n_central)
        parts.append(y_central)
    return np.abs(np.concatenate(parts))


# ===========================================================================
# Transform + single-peak readout
# ===========================================================================
def transform_density(r_values, bin_size, R_upper):
    """Run the actual pipeline transform; return (bin_centers, raw_density)."""
    inv = invert_2d_to_3d(bin_size, R_upper, r_values)
    density = np.asarray(inv.density, dtype=float)
    centers = (np.arange(1, len(density) + 1) * bin_size) - bin_size / 2.0
    return centers, density


def peak_position(centers, density, bin_size, refine=True):
    """Single recovered peak = argmax of the negative-clipped density. Optional
    parabolic interpolation of the argmax gives a sub-bin position so the
    histogram is not quantized to the bin grid. Returns nan if the density is
    empty / all <= 0.
    """
    d = np.clip(np.asarray(density, dtype=float), 0.0, None)
    if not np.isfinite(d).any() or d.max() <= 0.0:
        return float("nan")
    i = int(np.argmax(d))
    pos = centers[i]
    if refine and 0 < i < len(d) - 1:
        a, b, c = d[i - 1], d[i], d[i + 1]
        denom = a - 2.0 * b + c
        if denom != 0.0:
            delta = 0.5 * (a - c) / denom          # in units of bins, |delta|<=0.5
            pos = centers[i] + np.clip(delta, -0.5, 0.5) * bin_size
    return float(pos)


# ===========================================================================
# Monte Carlo over conditions
# ===========================================================================
def run_pool(rng, cfg, n_periph, n_central, progress=None):
    """One pool -> array of recovered peak positions over n_resample subsamples."""
    pool = make_pool(rng, cfg["R"], cfg["R_central"], cfg["precision"],
                     n_periph, n_central)
    n = len(pool)
    k = min(cfg["n_sample"], n)
    peaks = np.empty(cfg["n_resample"], dtype=float)
    for j in range(cfg["n_resample"]):
        idx = rng.choice(n, size=k, replace=False)
        centers, density = transform_density(pool[idx], cfg["bin_size"], cfg["R_upper"])
        peaks[j] = peak_position(centers, density, cfg["bin_size"], cfg["refine_peak"])
        if progress is not None:
            progress.update(1)
    return peaks


def run_condition(rng, cfg, n_periph, n_central, progress=None):
    """Run all pool reps for one condition. Returns (peaks_all, peaks_by_rep)."""
    by_rep = []
    for _ in range(cfg["n_pool_reps"]):
        by_rep.append(run_pool(rng, cfg, n_periph, n_central, progress))
    return np.concatenate(by_rep), by_rep


def band_fractions(peaks, R, R_central):
    """Split peaks at the midpoint between the two true modes; report the
    fraction in the central band and the peripheral band (ignoring nans)."""
    p = peaks[np.isfinite(peaks)]
    if p.size == 0:
        return float("nan"), float("nan"), 0
    mid = 0.5 * (R_central + R)
    frac_central = float(np.mean(p <= mid))
    return frac_central, 1.0 - frac_central, int(p.size)


def mode_positions(peaks, R, R_central):
    """Median peak position within each band -> recovered mode locations."""
    p = peaks[np.isfinite(peaks)]
    mid = 0.5 * (R_central + R)
    lo = p[p <= mid]
    hi = p[p > mid]
    m_central = float(np.median(lo)) if lo.size else float("nan")
    m_periph = float(np.median(hi)) if hi.size else float("nan")
    return m_central, m_periph


def main(cfg=None):
    cfg = cfg or CFG
    out_dir = Path(__file__).resolve().parent / "bimodal_resampling_output"
    out_dir.mkdir(exist_ok=True)
    rng = np.random.default_rng(cfg["seed"])

    pool_n = cfg["N1"] + cfg["N2"]
    conditions = [("mixture", cfg["N1"], cfg["N2"])]
    if cfg.get("run_controls", True):
        conditions += [
            ("all_peripheral", pool_n, 0),   # same total count, one component
            ("all_central", 0, pool_n),
        ]

    rows, summary, hist_data, per_pool = [], [], {}, {}
    input_ratio = cfg["N2"] / cfg["N1"] if cfg["N1"] else float("inf")  # central:periph

    total_draws = len(conditions) * cfg["n_pool_reps"] * cfg["n_resample"]
    bar = tqdm(total=total_draws, desc="subsamples")

    for label, n_periph, n_central in conditions:
        bar.desc = f"subsamples [{label}]"
        peaks_all, by_rep = run_condition(rng, cfg, n_periph, n_central, bar)
        hist_data[label] = peaks_all
        per_pool[label] = by_rep

        frac_c, frac_p, n_valid = band_fractions(peaks_all, cfg["R"], cfg["R_central"])
        m_c, m_p = mode_positions(peaks_all, cfg["R"], cfg["R_central"])
        summary.append(dict(
            condition=label, n_periph=n_periph, n_central=n_central,
            n_peaks=n_valid, frac_nan=float(np.mean(~np.isfinite(peaks_all))),
            mode_central=m_c, mode_peripheral=m_p,
            frac_in_central_band=frac_c, frac_in_peripheral_band=frac_p,
            input_central_over_periph=(n_central / n_periph) if n_periph else float("inf"),
        ))
        for rep_i, pk in enumerate(by_rep):
            for v in pk:
                rows.append(dict(condition=label, pool_rep=rep_i, peak_position=v))

    bar.close()

    pd.DataFrame(rows).to_csv(out_dir / "peaks_raw.csv", index=False)
    sum_df = pd.DataFrame(summary)
    sum_df.to_csv(out_dir / "summary.csv", index=False)

    # ---- overlaid histograms: mixture vs controls ----
    allv = np.concatenate([v[np.isfinite(v)] for v in hist_data.values()])
    lo, hi = float(np.min(allv)), float(np.max(allv))
    bins = np.linspace(lo, hi, 60)
    colors = {"mixture": "#DD8452", "all_peripheral": "#4C72B0",
              "all_central": "#55A868"}
    fig, ax = plt.subplots(figsize=(8, 5))
    for label in hist_data:
        v = hist_data[label]
        v = v[np.isfinite(v)]
        ax.hist(v, bins=bins, alpha=0.5, density=True,
                label=f"{label} (n={v.size})", color=colors.get(label))
    ax.axvline(cfg["R_central"], color="k", ls=":", lw=1.0)
    ax.axvline(cfg["R"], color="k", ls="--", lw=1.0)
    ax.set_xlabel("recovered peak position (nm)")
    ax.set_ylabel("normalized frequency")
    ax.set_title(
        f"Peak-position histogram  (R={cfg['R']:g}, central={cfg['R_central']:g}, "
        f"prec={cfg['precision']:g} nm, pool {cfg['N1']}+{cfg['N2']}, "
        f"sample {cfg['n_sample']}, {cfg['n_resample']}x{cfg['n_pool_reps']} draws)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "histograms.png", dpi=150)
    plt.close(fig)

    # ---- mixture histogram split by pool rep (RNG-stability check) ----
    fig, ax = plt.subplots(figsize=(8, 5))
    for rep_i, pk in enumerate(per_pool["mixture"]):
        v = pk[np.isfinite(pk)]
        ax.hist(v, bins=bins, histtype="step", lw=1.2, density=True,
                label=f"pool {rep_i}" if cfg["n_pool_reps"] <= 12 else None)
    ax.axvline(cfg["R_central"], color="k", ls=":", lw=1.0)
    ax.axvline(cfg["R"], color="k", ls="--", lw=1.0)
    ax.set_xlabel("recovered peak position (nm)")
    ax.set_ylabel("normalized frequency")
    ax.set_title("Mixture peak-position histogram per pool regeneration")
    if cfg["n_pool_reps"] <= 12:
        ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out_dir / "per_pool.png", dpi=150)
    plt.close(fig)

    pd.set_option("display.width", 170, "display.max_columns", 20)
    print("\n=== Per-condition summary ===")
    print(sum_df.to_string(index=False))
    print(f"\nInput central:peripheral ratio (mixture) = {input_ratio:.3g}")
    print(f"Outputs written to: {out_dir}")


if __name__ == "__main__":
    main()
