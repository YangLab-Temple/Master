# -*- coding: utf-8 -*-
"""
Ground-truth reference density + resample-histogram comparison
==============================================================

Purpose
-------
Sensitivity test for the reviewer: does the random-subsampling / peak-position
method reproduce the TWO-component structure of a known mixed distribution?

This script builds a transform-INDEPENDENT ground truth and compares it to the
resample-peak histogram produced by `bimodal_resampling_simulation.py`.

  1. Reference (true spatial density).  Generate a huge population (each N1, N2
     scaled by REF_SCALE) of TRUE 2D radial positions -- peripheral points on a
     ring of radius R, central points at R_central, each blurred by localization
     precision -- then histogram along the radial coordinate r and divide each
     bin's count by its exact ring area  pi*(r_hi^2 - r_lo^2).  The result is the
     molecules-per-unit-area (areal) density of the input mixture.  Because the
     huge N suppresses sampling noise, the bin can be small.  Crucially this uses
     NO 2D->3D transform, so it cannot be circular: it is the known input density.

  2. Resample histogram.  Read the peak positions the main script already wrote
     (bimodal_resampling_output/peaks_raw.csv, `mixture` condition) and bin them.

  3. Compare.  Overlay the two, each normalized to unit maximum.

What the comparison does and does NOT show
------------------------------------------
The reference is an areal DENSITY; the resample histogram is the distribution of
the single argmax peak across subsamples. They are different objects:
  * POSITIONS should match -- the histogram modes should sit on the reference
    peaks (r ~ 0 and r ~ R). This validates two-peak recovery / resolvability.
  * HEIGHTS will NOT match -- argmax is winner-take-all, so the taller (central)
    peak wins disproportionately and the histogram under-represents the
    peripheral peak relative to its true areal weight. Compare shape/positions,
    not amplitudes.

IMPORTANT: the reference is built on the TRUE radial coordinate r = sqrt(Y^2+Z^2),
not on the projected |Y| values fed to the transform. Area-normalizing the
projection would give the projection, not the spatial density.

Outputs (in ./bimodal_resampling_output/)
------------------------------------------
  reference_density.csv    r_lo, r_hi, r_center, count, ring_area, areal_density,
                           areal_density_norm   (small bins; for OriginLab)
  resample_histogram.csv   r_lo, r_hi, r_center, count, frequency, frequency_norm
  reference_vs_resample.png  overlay, each normalized to unit max

Author: Wenlan Yu / harness 2026.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Pull the SAME parameters used by the main simulation so the reference matches
# the resampling run (R, R_central, precision, N1, N2, seed).
try:
    from bimodal_resampling_simulation import CFG
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "Could not import CFG from bimodal_resampling_simulation.py -- run this "
        "script from the same folder as that file.\n  (%s)" % exc)

# ---- knobs for this script -------------------------------------------------
REF_SCALE = 20000      # multiply N1, N2 by this for the huge ground-truth sample.
REF_BIN = 0.5         # radial bin (nm) for the reference; small is fine at huge N.
HIST_BIN = 10.0        # radial bin (nm) for the resample-peak histogram.
R_MAX = None          # plot/CSV radial limit (nm); None -> auto from R.
PEAKS_CSV = "bimodal_resampling_output/peaks_raw.csv"
CONDITION = "mixture"


def sample_true_radii(rng, R, R_central, precision, n_periph, n_central):
    """TRUE radial coordinate r = sqrt(Y^2 + Z^2) for each molecule: peripheral on
    a ring of radius R, central at R_central, each blurred isotropically by the
    localization precision. This is the molecules' real spatial distribution (the
    best-case recoverable density), independent of the 2D->3D transform.
    """
    parts = []
    if n_periph > 0:
        phi = rng.uniform(0.0, 2.0 * np.pi, n_periph)
        x = R * np.cos(phi) + rng.normal(0.0, precision, n_periph)
        y = R * np.sin(phi) + rng.normal(0.0, precision, n_periph)
        parts.append(np.hypot(x, y))
    if n_central > 0:
        phi = rng.uniform(0.0, 2.0 * np.pi, n_central)
        x = R_central * np.cos(phi) + rng.normal(0.0, precision, n_central)
        y = R_central * np.sin(phi) + rng.normal(0.0, precision, n_central)
        parts.append(np.hypot(x, y))
    return np.concatenate(parts)


def build_reference(cfg, out_dir, r_max):
    rng = np.random.default_rng(cfg["seed"] + 1)
    r = sample_true_radii(rng, cfg["R"], cfg["R_central"], cfg["precision"],
                          cfg["N1"] * REF_SCALE, cfg["N2"] * REF_SCALE)
    edges = np.arange(0.0, r_max + REF_BIN, REF_BIN)
    counts, _ = np.histogram(r, bins=edges)
    r_lo, r_hi = edges[:-1], edges[1:]
    centers = 0.5 * (r_lo + r_hi)
    ring_area = np.pi * (r_hi ** 2 - r_lo ** 2)          # exact annulus area
    areal = counts / ring_area
    df = pd.DataFrame(dict(
        r_lo=r_lo, r_hi=r_hi, r_center=centers, count=counts,
        ring_area=ring_area, areal_density=areal,
        areal_density_norm=areal / areal.max() if areal.max() > 0 else areal))
    df.to_csv(out_dir / "reference_density.csv", index=False)
    return df


def build_resample_hist(out_dir, r_max):
    p = Path(PEAKS_CSV)
    if not p.is_absolute():
        p = out_dir.parent / PEAKS_CSV
    if not p.exists():
        print(f"[skip] {p} not found -- run bimodal_resampling_simulation.py "
              f"first to produce the resample peaks.")
        return None
    raw = pd.read_csv(p)
    peaks = raw.loc[raw["condition"] == CONDITION, "peak_position"].dropna().values
    edges = np.arange(0.0, r_max + HIST_BIN, HIST_BIN)
    counts, _ = np.histogram(peaks, bins=edges)
    r_lo, r_hi = edges[:-1], edges[1:]
    centers = 0.5 * (r_lo + r_hi)
    total = counts.sum()
    df = pd.DataFrame(dict(
        r_lo=r_lo, r_hi=r_hi, r_center=centers, count=counts,
        frequency=counts / total if total else counts,
        frequency_norm=counts / counts.max() if counts.max() > 0 else counts))
    df.to_csv(out_dir / "resample_histogram.csv", index=False)
    return df


def main(cfg=None):
    cfg = cfg or CFG
    out_dir = Path(__file__).resolve().parent / "bimodal_resampling_output"
    out_dir.mkdir(exist_ok=True)
    r_max = R_MAX if R_MAX is not None else float(max(60.0, 2.2 * cfg["R"]))

    ref = build_reference(cfg, out_dir, r_max)
    hist = build_resample_hist(out_dir, r_max)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ref["r_center"], ref["areal_density_norm"], "-", color="#333333",
            lw=2, label="true areal density (reference)")
    if hist is not None:
        ax.bar(hist["r_center"], hist["frequency_norm"], width=HIST_BIN,
               alpha=0.5, color="#DD8452", label="resample-peak histogram")
    ax.axvline(cfg["R_central"], color="k", ls=":", lw=1.0)
    ax.axvline(cfg["R"], color="k", ls="--", lw=1.0)
    ax.set_xlim(0, r_max)
    ax.set_xlabel("radial position r (nm)")
    ax.set_ylabel("normalized to unit max")
    ax.set_title(
        f"Reference density vs resample histogram  "
        f"(R={cfg['R']:g}, central={cfg['R_central']:g}, prec={cfg['precision']:g}, "
        f"N1:N2={cfg['N1']}:{cfg['N2']}, ref x{REF_SCALE})")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "reference_vs_resample.png", dpi=150)
    plt.close(fig)

    # quick text readout of where the reference peaks sit
    rr = ref["r_center"].values
    dd = ref["areal_density"].values
    central_pk = rr[np.argmax(np.where(rr <= 0.5 * cfg["R"], dd, -np.inf))]
    periph_pk = rr[np.argmax(np.where(rr > 0.5 * cfg["R"], dd, -np.inf))]
    print(f"reference peaks: central ~ {central_pk:.1f} nm, "
          f"peripheral ~ {periph_pk:.1f} nm (true R = {cfg['R']:g})")
    print(f"outputs written to: {out_dir}")


if __name__ == "__main__":
    main()
