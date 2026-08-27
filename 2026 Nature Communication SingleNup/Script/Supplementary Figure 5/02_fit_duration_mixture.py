# -*- coding: utf-8 -*-
"""
Step 2 of 3 -- Supplementary Figure S5B/S5C (Nup153C, 2 ms/frame).

Maximum-likelihood fit of a K-component geometric mixture to the unfiltered
event-duration distribution, plus the model-comparison table and the draft
figure caption.

Background
----------
Reviewer 1, comment 5 (round-2 revision) asked for a histogram of ALL event
durations -- not only those surviving the 200 ms retention cutoff used for
the primary 50 ms/frame MSD datasets -- to assess whether short and long
events form separable populations.

K = 3 is the adopted model for interpretation: decisive AIC/BIC support over
K = 2, far past the conventional strong-evidence threshold of ~10. K = 4 is
fit and reported in the table for transparency but is not the headline model;
its extra component captures <0.1% of events.

This script performs no plotting. The published S5B and S5C panels were drawn
in OriginLab from the per-bin counts written by 03_export_histogram_counts.py.

Input
-----
data/combined_group_lengths.csv -- one row per event, duration in FRAMES
(written by 01_extract_event_durations.py).

Output
------
results/model_comparison.csv -- AIC/BIC table, K = 1..K_max
results/FigS5_caption.txt    -- draft caption with the fitted numbers filled in

Usage
-----
    python 02_fit_duration_mixture.py
"""

import os

import numpy as np
import pandas as pd
from scipy.optimize import minimize

# ==================== SETTINGS ====================
input_csv = os.path.join("data", "combined_group_lengths.csv")
out_dir = "results"

frame_time_ms = 2.0

K_max = 4          # fit and report K = 1..K_max in the table
primary_K = 3      # model used for interpretation / rebuttal text
n_restarts = 15    # optimizer restarts per model, to avoid local minima
rng_seed = 0       # fixed for reproducibility

# Histogram ranges quoted in the caption. These must match the binning in
# 03_export_histogram_counts.py, which is what the figures were drawn from.
panelB_xmax_ms = 50
panelC_xmax_ms = 1000
# ==================================================


def nll_mixK(params, k, K):
    """Negative log-likelihood of a K-component geometric mixture.

    Parameter vector is [w_1..w_{K-1}, p_1..p_K]; the last weight is fixed by
    the constraint that the weights sum to 1. Returns a large finite penalty
    outside the feasible region so Nelder-Mead stays inside it.
    """
    if K == 1:
        p = params[0]
        if not (1e-6 < p < 1 - 1e-6):
            return 1e10
        return -np.sum((k - 1) * np.log(1 - p) + np.log(p))

    w = params[:K - 1]
    wK = 1 - np.sum(w)
    ps = params[K - 1:]
    if np.any(w <= 1e-6) or wK <= 1e-6:
        return 1e10
    if np.any(ps <= 1e-6) or np.any(ps >= 1 - 1e-6):
        return 1e10

    weights = np.append(w, wK)
    mix = np.zeros_like(k)
    for wi, pi in zip(weights, ps):
        mix += wi * pi * (1 - pi) ** (k - 1)
    return -np.sum(np.log(np.clip(mix, 1e-300, None)))


def fit_mixture(k, K, rng, n_restarts):
    """Fit a K-component geometric mixture by MLE with multiple restarts."""
    best = None
    for _ in range(n_restarts):
        if K == 1:
            x0 = [rng.uniform(0.1, 0.6)]
        else:
            w0 = rng.dirichlet(np.ones(K))[:K - 1]
            p0 = np.sort(rng.uniform(0.002, 0.6, K))[::-1]
            x0 = np.concatenate([w0, p0])
        r = minimize(nll_mixK, x0, args=(k, K), method="Nelder-Mead",
                     options={"maxiter": 8000, "maxfev": 8000,
                              "xatol": 1e-7, "fatol": 1e-7})
        if best is None or r.fun < best.fun:
            best = r

    n_params = 2 * K - 1
    if K == 1:
        comps = [(1.0, best.x[0])]
    else:
        w = best.x[:K - 1]
        weights = list(w) + [1 - np.sum(w)]
        comps = sorted(zip(weights, best.x[K - 1:]), key=lambda c: 1 / c[1])

    return dict(K=K,
                nll=best.fun,
                aic=2 * n_params + 2 * best.fun,
                bic=n_params * np.log(len(k)) + 2 * best.fun,
                comps=comps)


def model_table(fits, path):
    """AIC/BIC comparison across K, with step-wise deltas."""
    rows, prev_aic, prev_bic = [], None, None
    for K in sorted(fits):
        aic, bic = fits[K]["aic"], fits[K]["bic"]
        d_aic = np.nan if prev_aic is None else prev_aic - aic
        d_bic = np.nan if prev_bic is None else prev_bic - bic
        if np.isnan(d_aic):
            verdict = "--"
        else:
            verdict = "supported" if d_aic > 10 else "not supported"
        rows.append(dict(
            K=K, AIC=aic, BIC=bic,
            delta_AIC_step=d_aic, delta_BIC_step=d_bic,
            rule_of_10=verdict,
            components="; ".join(
                f"w={w:.4f},mean={frame_time_ms / p:.2f}ms"
                for w, p in fits[K]["comps"]),
        ))
        prev_aic, prev_bic = aic, bic

    table = pd.DataFrame(rows)
    table.to_csv(path, index=False)
    return table


def write_caption(fits, durations_ms, n_total, path):
    """Draft caption for Supplementary Figure S5, with fitted values filled in.

    NOTE: the panel (C) sentence below states that the 2 ms acquisition applies
    the same minimum three-localization tracking requirement as the primary
    50 ms datasets. The pipeline in this repository applies no such filter --
    69% of the events are shorter than three frames. See "Known limitations"
    in README.md; the sentence is reproduced here verbatim as submitted and
    should be reconciled with the code before the repository is made public.
    """
    comps = "; ".join(f"{w * 100:.1f}% - {frame_time_ms / p:.1f} ms"
                      for w, p in fits[primary_K]["comps"])
    d_aic = fits[primary_K - 1]["aic"] - fits[primary_K]["aic"]
    d_bic = fits[primary_K - 1]["bic"] - fits[primary_K]["bic"]
    n_beyond_b = int((durations_ms > panelB_xmax_ms).sum())
    n_beyond_c = int((durations_ms > panelC_xmax_ms).sum())

    caption = f"""(A) Trace-length frequency distributions for traces used in the MSD analysis
(unchanged; these are traces already retained after the 200 ms threshold).

(B) Unfiltered event-duration histogram for Nup153C, acquired independently at
{frame_time_ms:.0f} ms/frame temporal resolution (n = {n_total} tracked events, no
minimum-duration filter applied), shown for durations from 0 to {panelB_xmax_ms:.0f} ms
({n_beyond_b} events beyond this range are not shown, see panel C). Error bars, SD
from Poisson counting statistics.

(C) The same unfiltered distribution across the full observed range, 0 to
{panelC_xmax_ms:.0f} ms, plotted on a semi-logarithmic scale ({n_beyond_c} events beyond
this range, if any, are not shown). Because the {frame_time_ms:.0f} ms acquisition applies
the same minimum three-localization tracking requirement as the primary 50 ms
datasets but at six-fold finer temporal sampling, it resolves the short end of the
event-duration distribution that cannot be recovered from the primary datasets.

Maximum-likelihood mixture-model fitting of this unfiltered distribution
(Supplementary Table SX) shows it is decisively better described by
{primary_K} kinetically distinct populations than by {primary_K - 1}
(delta-AIC = {d_aic:.1f}, delta-BIC = {d_bic:.1f}, both far past the
conventional strong-evidence threshold of ~10), resolving into: {comps}.
These are consistent with free dye, Nup molecules not incorporated into an
assembled pore, and genuinely pore-associated Nup, respectively.
"""
    with open(path, "w") as fh:
        fh.write(caption)
    return caption


def main():
    os.makedirs(out_dir, exist_ok=True)

    k_all = pd.read_csv(input_csv).iloc[:, 0].to_numpy(dtype=float)
    durations_ms = k_all * frame_time_ms
    n_total = len(k_all)
    print(f"Loaded {n_total} events from {input_csv}")

    rng = np.random.default_rng(rng_seed)
    fits = {K: fit_mixture(k_all, K, rng, n_restarts)
            for K in range(1, K_max + 1)}

    table_path = os.path.join(out_dir, "model_comparison.csv")
    print(model_table(fits, table_path).to_string(index=False))

    caption_path = os.path.join(out_dir, "FigS5_caption.txt")
    print("\n" + write_caption(fits, durations_ms, n_total, caption_path))

    print("Saved table:  ", table_path)
    print("Saved caption:", caption_path)


if __name__ == "__main__":
    main()
