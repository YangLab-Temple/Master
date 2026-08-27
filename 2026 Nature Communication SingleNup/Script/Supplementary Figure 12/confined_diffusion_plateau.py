# -*- coding: utf-8 -*-
"""
Confined-diffusion plateau analysis for FG-Nup single-molecule trajectories.

Purpose (Reviewer 1, comment 11):
    Replace the sample-size-dependent "max extension" metric with a
    well-defined, sample-size-INDEPENDENT confinement measure by fitting a
    confined-diffusion model to the mean-square displacement (MSD) and
    reporting the PLATEAU (equivalently the RMS extension / position SD).

What it does, per protein and per spatial dimension (x, y):
    1. Re-extracts the time-averaged MSD separately for x and y from the raw
       trajectories (msdanalyzer exports a *combined* 2D <dr^2>; the Kusumi
       equation is 1D, so we need per-dimension curves).
    2. Fits TWO confinement models over the WELL-SAMPLED lag range only:
         - OU / harmonic:  MSD = plateau*(1-exp(-t/tau)) + 2*sig_loc^2
                           plateau = 2*sigma^2  (sigma = position SD)
         - Kusumi box (Eq. 11 of Kusumi 1993, Biophys J 65:2021):
                           MSD = L^2/6 - (16 L^2/pi^4) * sum_{n odd} (1/n^4)
                                 * exp(-n^2 pi^2 D t / L^2) + 2*sig_loc^2
                           plateau = L^2/6
       The 2*sig_loc^2 localization-error offset is FIXED from Table S3, not fitted.
    3. Bootstraps over whole traces to get confidence intervals on every estimate.
    4. Combines x,y into an RMS radial extension and tabulates everything
       against the OLD Gaussian-fit "max extension" mean for comparison.

Author: Wenlan Yu  (analysis assembled 2026-07)
Run in Spyder: edit the SETTINGS block below, then run the whole file.
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore", category=RuntimeWarning)

# progress bar (falls back to a tiny shim if tqdm is absent)
try:
    from tqdm import tqdm
except Exception:                                            # pragma: no cover
    def tqdm(it, **kw):
        it = list(it)
        tot = len(it)
        for i, x in enumerate(it):
            if i % max(1, tot // 20) == 0:
                print(f"  ...{i}/{tot}", flush=True)
            yield x

# =====================================================================
# SETTINGS  (edit these)
# =====================================================================
# Base folder = the folder that CONTAINS the two data subfolders below.
BASE_DIR = os.path.dirname(os.path.abspath(__file__))     # default: this script's folder

# The datasets to process in one run: {label: subfolder name under BASE_DIR}.
DATASETS = {
    "FG":    "Supplemental MSD FG",
    "NonFG": "Supplemental MSD NonFG",
}

# Everything (combined table + all figures) is written here.
OUTPUT_DIR = os.path.join(BASE_DIR, "Confinement_Analysis_Output")

DT         = 0.05        # s, time between frames
N_DIM_LABELS = ("x", "y")

# --- fit range (this is the crux: the far plateau is thinly sampled) ---
MIN_TRACES_PER_LAG = 5   # drop lags supported by fewer than this many traces
MAX_LAG_S          = 0.30 # s, hard cap on lag used in the fit (short-trace reliability)
MIN_LAGS_TO_FIT    = 4   # need at least this many lag points to attempt a fit

# --- bootstrap ---
N_BOOT   = 1000
RNG_SEED = 0

# --- Kusumi series ---
KUSUMI_N_ODD = 25        # number of odd terms in the series (1/n^4 -> converges fast)

# --- output ---
SAVE_FIGURES   = True

# Single-molecule localization precision sigma_loc (nm), per protein,
# from Supplementary Table S3 (and S2). Keyed by trajectory-file basename.
SIGMA_LOC = {
    "hCG1C": 7.7, "hCG1N": 6.5, "Nup153C": 8.7, "Nup214C": 7.3, "Nup358C": 8.8,
    "Nup50N": 8.9, "Nup54N": 6.3, "Nup58C": 6.2, "Nup58N": 8.0, "Nup62N": 8.2,
    "Nup98N": 6.5, "POM121C": 8.1, "TPRC": 9.4,
    "Nup153N": 6.1, "Nup214N": 8.4, "Nup358N": 8.4, "Nup50C": 7.9, "Nup54C": 6.3,
    "Nup62C": 7.2, "Nup98C": 6.8, "POM121N": 7.2, "TPRN": 9.2,
}
DEFAULT_SIGMA_LOC = 8.0   # nm, fallback if a protein is missing above

# =====================================================================
# plotting config (per-figure)
# =====================================================================
mpl.rcParams["figure.dpi"] = 130
mpl.rcParams["font.size"] = 12
mpl.rcParams["font.family"] = ["DejaVu Sans"]   # Arial-like; change to Arial locally
PLOTCFG = {
    "data_color":   "black",
    "band_alpha":   0.25,
    "band_color":   "grey",
    "ou_color":     "#1f77b4",
    "kusumi_color": "#d62728",
    "fitspan_color":"#fff2cc",
}

_ODD = np.arange(1, 2 * KUSUMI_N_ODD, 2)      # 1,3,5,...


# =====================================================================
# data loading
# =====================================================================
def load_traces(path):
    """Return list of (x, y) arrays. Traces are separated by all-NaN rows."""
    df = pd.read_csv(path, header=None)
    df = df.iloc[:, :2].apply(pd.to_numeric, errors="coerce")
    grp = df.isnull().all(axis=1).cumsum()
    traces = []
    for _, block in df.groupby(grp):
        b = block.dropna().values
        if b.shape[0] >= 3:                    # need >=3 points to be useful
            traces.append((b[:, 0].astype(float), b[:, 1].astype(float)))
    return traces


def msd_one_dim(traces, dim, max_lag):
    """
    Time-averaged, ensemble-pooled MSD for one dimension.
    Returns dict of arrays keyed by lag: lag_s, msd, sem, n_pairs, n_traces.
    """
    per_lag_sq = {k: [] for k in range(1, max_lag + 1)}
    per_lag_tr = {k: 0 for k in range(1, max_lag + 1)}
    for (x, y) in traces:
        pos = x if dim == 0 else y
        L = len(pos)
        for k in range(1, min(max_lag, L - 1) + 1):
            d = pos[k:] - pos[:-k]
            per_lag_sq[k].append(d * d)
            per_lag_tr[k] += 1
    lags, msd, sem, npair, ntr = [], [], [], [], []
    for k in range(1, max_lag + 1):
        if per_lag_sq[k]:
            allsq = np.concatenate(per_lag_sq[k])
            lags.append(k * DT)
            msd.append(allsq.mean())
            sem.append(allsq.std(ddof=1) / np.sqrt(len(allsq)) if len(allsq) > 1 else np.nan)
            npair.append(len(allsq))
            ntr.append(per_lag_tr[k])
    return dict(lag_s=np.array(lags), msd=np.array(msd), sem=np.array(sem),
                n_pairs=np.array(npair), n_traces=np.array(ntr))


def fit_mask(m):
    """Boolean mask of lags to include in the fit."""
    return (m["n_traces"] >= MIN_TRACES_PER_LAG) & (m["lag_s"] <= MAX_LAG_S + 1e-9)


# =====================================================================
# models  (offset = 2*sigma_loc^2 is fixed, passed via closure)
# =====================================================================
def make_ou(offset):
    def ou(t, plateau, tau):
        return plateau * (1.0 - np.exp(-t / tau)) + offset
    return ou


def make_kusumi(offset):
    def kusumi(t, L, D):
        t = np.asarray(t, float)
        plateau = L * L / 6.0
        n = _ODD[:, None]
        expo = np.exp(-(n ** 2) * (np.pi ** 2) * D * t[None, :] / (L * L))
        series = (16.0 * L * L / np.pi ** 4) * np.sum(expo / n ** 4, axis=0)
        return plateau - series + offset
    return kusumi


def _r2(y, yhat, w):
    ybar = np.average(y, weights=w)
    ss_res = np.sum(w * (y - yhat) ** 2)
    ss_tot = np.sum(w * (y - ybar) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan


def fit_dim(m, sigma_loc):
    """Fit OU and Kusumi to one dimension's MSD. Returns dict of results."""
    off = 2.0 * sigma_loc ** 2
    msk = fit_mask(m)
    out = {"n_fit_lags": int(msk.sum()), "max_lag_fit_s": float(m["lag_s"][msk].max()) if msk.any() else np.nan}
    if msk.sum() < MIN_LAGS_TO_FIT:
        out.update(dict(ou_plateau=np.nan, ou_sigma=np.nan, ou_tau=np.nan, ou_r2=np.nan,
                        ku_plateau=np.nan, ku_L=np.nan, ku_D=np.nan, ku_r2=np.nan))
        return out
    t, y = m["lag_s"][msk], m["msd"][msk]
    sem = m["sem"][msk]
    sem = np.where(np.isfinite(sem) & (sem > 0), sem, np.nanmedian(sem[sem > 0]) if np.any(sem > 0) else 1.0)
    w = 1.0 / sem ** 2
    plateau0 = max(np.nanmax(y) - off, 1.0)

    # OU
    try:
        ou = make_ou(off)
        p, _ = curve_fit(ou, t, y, p0=[plateau0, 0.1], sigma=sem, absolute_sigma=True,
                         bounds=([0, 1e-3], [np.inf, 10]), maxfev=20000)
        out["ou_plateau"], out["ou_tau"] = p[0], p[1]
        out["ou_sigma"] = np.sqrt(max(p[0], 0) / 2.0)          # position SD per dim
        out["ou_r2"] = _r2(y, ou(t, *p), w)
    except Exception:
        out.update(ou_plateau=np.nan, ou_sigma=np.nan, ou_tau=np.nan, ou_r2=np.nan)

    # Kusumi box
    try:
        ku = make_kusumi(off)
        L0 = np.sqrt(6 * plateau0)
        D0 = plateau0 / (2 * 0.1)
        p, _ = curve_fit(ku, t, y, p0=[L0, D0], sigma=sem, absolute_sigma=True,
                         bounds=([1, 1e-2], [1e4, 1e7]), maxfev=20000)
        out["ku_L"], out["ku_D"] = p[0], p[1]
        out["ku_plateau"] = p[0] ** 2 / 6.0
        out["ku_r2"] = _r2(y, ku(t, *p), w)
    except Exception:
        out.update(ku_plateau=np.nan, ku_L=np.nan, ku_D=np.nan, ku_r2=np.nan)
    return out


# =====================================================================
# per-protein driver
# =====================================================================
def process_protein(name, traces, sigma_loc, rng):
    max_lag = max(len(x) - 1 for x, _ in traces)
    max_lag = min(max_lag, int(np.ceil((MAX_LAG_S * 2) / DT)) + 5)
    _npts = np.array([len(x) for x, _ in traces])
    res = {"protein": name, "n_traces": len(traces),
           "mean_pts_per_track": float(_npts.mean()), "sd_pts_per_track": float(_npts.std()),
           "sigma_loc_nm": sigma_loc}
    curves = {}
    # point estimate
    for dim, lab in enumerate(N_DIM_LABELS):
        m = msd_one_dim(traces, dim, max_lag)
        curves[lab] = m
        f = fit_dim(m, sigma_loc)
        for k, v in f.items():
            res[f"{lab}_{k}"] = v

    # RMS radial extension from point estimate (Var_dim = plateau_dim / 2)
    def rms_radial(pref):
        px, py = res.get(f"x_{pref}_plateau", np.nan), res.get(f"y_{pref}_plateau", np.nan)
        return np.sqrt((px + py) / 2.0)
    res["rms_ext_ou_nm"] = rms_radial("ou")
    res["rms_ext_ku_nm"] = rms_radial("ku")

    # bootstrap over traces
    boot = {"ou": [], "ku": []}
    tr = np.array(traces, dtype=object)
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(tr), len(tr))
        bs = list(tr[idx])
        vals = {}
        ok = True
        for dim, lab in enumerate(N_DIM_LABELS):
            mb = msd_one_dim(bs, dim, max_lag)
            fb = fit_dim(mb, sigma_loc)
            vals[lab] = fb
            if not np.isfinite(fb.get("ou_plateau", np.nan)):
                ok = False
        if ok:
            for pref, key in (("ou", "ou_plateau"), ("ku", "ku_plateau")):
                px, py = vals["x"].get(key, np.nan), vals["y"].get(key, np.nan)
                boot[pref].append(np.sqrt((px + py) / 2.0))
    for pref in ("ou", "ku"):
        arr = np.array([b for b in boot[pref] if np.isfinite(b)])
        if arr.size:
            res[f"rms_ext_{pref}_lo"] = np.percentile(arr, 2.5)
            res[f"rms_ext_{pref}_hi"] = np.percentile(arr, 97.5)
            res[f"rms_ext_{pref}_boot_med"] = np.median(arr)
        else:
            res[f"rms_ext_{pref}_lo"] = res[f"rms_ext_{pref}_hi"] = res[f"rms_ext_{pref}_boot_med"] = np.nan
    return res, curves


def plot_protein(name, curves, res, fig_dir):
    off_txt = f"$\\sigma_{{loc}}$={res['sigma_loc_nm']:.1f} nm"
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for dim, lab in enumerate(N_DIM_LABELS):
        ax = axes[dim]
        m = curves[lab]
        msk = fit_mask(m)
        # full data + SEM band
        ax.fill_between(m["lag_s"], m["msd"] - m["sem"], m["msd"] + m["sem"],
                        color=PLOTCFG["band_color"], alpha=PLOTCFG["band_alpha"], label="SEM")
        ax.plot(m["lag_s"], m["msd"], "o-", color=PLOTCFG["data_color"], ms=3, lw=1, label="MSD data")
        # shade fit range
        if msk.any():
            ax.axvspan(0, m["lag_s"][msk].max(), color=PLOTCFG["fitspan_color"], zorder=0,
                       label="fit range")
        tt = np.linspace(0, m["lag_s"].max(), 200)
        off = 2 * res["sigma_loc_nm"] ** 2
        if np.isfinite(res[f"{lab}_ou_plateau"]):
            ou = make_ou(off)
            ax.plot(tt, ou(tt, res[f"{lab}_ou_plateau"], res[f"{lab}_ou_tau"]),
                    color=PLOTCFG["ou_color"], lw=2,
                    label=f"OU  plat={res[f'{lab}_ou_plateau']:.0f}")
            ax.axhline(res[f"{lab}_ou_plateau"] + off, color=PLOTCFG["ou_color"], ls=":", lw=1)
        if np.isfinite(res[f"{lab}_ku_plateau"]):
            ku = make_kusumi(off)
            ax.plot(tt, ku(tt, res[f"{lab}_ku_L"], res[f"{lab}_ku_D"]),
                    color=PLOTCFG["kusumi_color"], lw=2, ls="--",
                    label=f"Kusumi plat={res[f'{lab}_ku_plateau']:.0f}")
        ax.set_xlabel("Delay (s)")
        ax.set_ylabel(f"MSD$_{lab}$ (nm$^2$)")
        ax.set_title(f"{name} — {lab}")
        ax.spines[["right", "top"]].set_visible(False)
        ax.legend(fontsize=8, loc="lower right")
    fig.suptitle(f"{name}   ({res['n_traces']} traces, {off_txt})", fontweight="bold")
    fig.tight_layout()
    if SAVE_FIGURES:
        os.makedirs(fig_dir, exist_ok=True)
        fig.savefig(os.path.join(fig_dir, f"{name}_plateau_fit.png"), bbox_inches="tight")
    plt.close(fig)


# =====================================================================
# max-extension (msdanalyzer output) -- kept as a third reported approach
# =====================================================================
def load_max_extension(data_dir):
    """Return {protein: (mean, mean_SE, SD)} from the Setting CSV, if present."""
    out = {}
    for sf in ("Setting1.csv", "Setting2.csv"):
        sp = os.path.join(data_dir, sf)
        if os.path.exists(sp):
            s = pd.read_csv(sp)
            s.columns = [c.strip().lstrip("﻿") for c in s.columns]
            if "Nups" in s.columns and "Fitted Frequency mean" in s.columns:
                for _, r in s.iterrows():
                    out[str(r["Nups"]).strip()] = (
                        r.get("Fitted Frequency mean", np.nan),
                        r.get("Fitted Frequency mean SE", np.nan),
                        r.get("Fitted Frequency SD", np.nan),
                    )
            break
    return out


def process_dataset(label, data_dir, rng):
    """Process every protein in one folder; return list of result dicts."""
    files = sorted(glob.glob(os.path.join(data_dir, "*_Trajectory.csv")))
    files = [f for f in files if "_Joe" not in f]
    maxext = load_max_extension(data_dir)
    fig_dir = os.path.join(OUTPUT_DIR, f"{label}_plateau_fits")
    print(f"\n[{label}] {len(files)} trajectory files in {data_dir}")

    rows = []
    for f in tqdm(files, desc=f"{label} proteins"):
        name = os.path.basename(f).replace("_Trajectory.csv", "")
        sig = SIGMA_LOC.get(name, DEFAULT_SIGMA_LOC)
        traces = load_traces(f)
        if len(traces) < 2:
            print(f"  [skip] {name}: <2 usable traces")
            continue
        res, curves = process_protein(name, traces, sig, rng)
        res["dataset"] = label
        m, mse, sd = maxext.get(name, (np.nan, np.nan, np.nan))
        res["maxext_mean_nm"] = m
        res["maxext_mean_se_nm"] = mse
        res["maxext_sd_nm"] = sd
        rows.append(res)
        if SAVE_FIGURES:
            plot_protein(name, curves, res, fig_dir)
    return rows


def summary_figure(df, path):
    """Bar chart: three extension approaches side by side, all proteins."""
    d = df.sort_values(["dataset", "rms_ext_ou_nm"]).reset_index(drop=True)
    x = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(max(11, 0.55 * len(d)), 5.5))
    ou_err = np.vstack([d.rms_ext_ou_nm - d.rms_ext_ou_lo, d.rms_ext_ou_hi - d.rms_ext_ou_nm])
    ax.bar(x - 0.27, d.rms_ext_ou_nm, 0.26, yerr=ou_err, capsize=2,
           color="#1f77b4", label="RMS extension — OU plateau (95% CI)")
    ax.bar(x, d.rms_ext_ku_nm, 0.26, color="#2ca02c",
           label="RMS extension — Kusumi plateau")
    ax.bar(x + 0.27, d.maxext_mean_nm, 0.26, yerr=d.maxext_sd_nm, capsize=2,
           color="#bbbbbb", label="Max extension mean ± SD (msdanalyzer)")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{p}\n({g})" for p, g in zip(d.protein, d.dataset)],
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Extension (nm)")
    ax.set_title("FG-Nup confinement: three complementary extension metrics",
                 fontweight="bold")
    ax.spines[["right", "top"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# =====================================================================
# main
# =====================================================================
def main():
    rng = np.random.default_rng(RNG_SEED)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_rows = []
    for label, sub in DATASETS.items():
        data_dir = os.path.join(BASE_DIR, sub)
        if not os.path.isdir(data_dir):
            print(f"[warn] dataset folder not found, skipping: {data_dir}")
            continue
        all_rows += process_dataset(label, data_dir, rng)

    df = pd.DataFrame(all_rows)
    # order columns: identity first, then the three approaches, then details
    lead = ["dataset", "protein", "n_traces", "mean_pts_per_track", "sd_pts_per_track",
            "sigma_loc_nm",
            "rms_ext_ou_nm", "rms_ext_ou_lo", "rms_ext_ou_hi",
            "rms_ext_ku_nm", "rms_ext_ku_lo", "rms_ext_ku_hi",
            "maxext_mean_nm", "maxext_mean_se_nm", "maxext_sd_nm"]
    lead = [c for c in lead if c in df.columns]
    df = df[lead + [c for c in df.columns if c not in lead]]

    combined_csv = os.path.join(OUTPUT_DIR, "confinement_results_ALL.csv")
    df.to_csv(combined_csv, index=False)

    # compact per-approach summary table
    summ = df[[c for c in lead if c in df.columns]].copy()
    summ.to_csv(os.path.join(OUTPUT_DIR, "confinement_summary_table.csv"), index=False)

    if SAVE_FIGURES:
        summary_figure(df, os.path.join(OUTPUT_DIR, "SUMMARY_three_metrics.png"))

    pd.set_option("display.width", 200, "display.max_columns", 40)
    print("\n===== Extension metrics (nm), localization-error corrected =====")
    show = ["dataset", "protein", "n_traces", "rms_ext_ou_nm",
            "rms_ext_ou_lo", "rms_ext_ou_hi", "rms_ext_ku_nm", "maxext_mean_nm"]
    show = [c for c in show if c in df.columns]
    print(df[show].round(1).to_string(index=False))
    print(f"\nAll outputs -> {OUTPUT_DIR}")
    print("  - confinement_results_ALL.csv    (every fitted parameter)")
    print("  - confinement_summary_table.csv  (the three approaches side by side)")
    print("  - SUMMARY_three_metrics.png      (comparison bar chart)")
    print("  - <FG|NonFG>_plateau_fits/*.png  (per-protein MSD fits)")
    return df


if __name__ == "__main__":
    main()
