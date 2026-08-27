# -*- coding: utf-8 -*-
"""
MINFLUX 3D MSD pipeline: segment-level QC -> log-binned MSD -> model fits.

Reads the raw per-construct .xlsx folders and does NOT modify them. Replaces
MSD_MINFLUX.py + MSD_MINFLUX_robust.py for this purpose and is much faster:
it evaluates ~N_BINS target lag times instead of every integer lag, so it never
writes the multi-GB per-trace CSVs, and it processes .xlsx files in parallel.

LAYOUT
    Put every construct folder you want analysed inside a folder called `data`:

        Calculation/
            data/
                GFP-Nup54-HaloTag/    *.xlsx
                GFP-Nup58-HaloTag/    *.xlsx
                ...
            MINFLUX_MSD_pipeline.py
            MSD_v2/                   <- created by this script

    Anything not inside `data` is ignored, so retired constructs can simply be
    moved out. If no `data` folder is found the script falls back to scanning
    for <construct>/*.xlsx next to itself.

STAGE 1  Segment-level QC
    Terminal (and mid-trace) tracking failures are found with a centred rolling
    median / rolling SD on each axis. A localization is flagged when its local
    median departs from the trajectory median by more than QC_DEV_NM, or when
    its local SD exceeds QC_SD_FOLD x the trajectory's robust baseline SD. The
    longest contiguous run of unflagged localizations is kept. Because the
    rolling window is centred, the cut lands roughly QC_WINDOW/2 localizations
    before the true onset, i.e. the screen is slightly conservative. Every
    trimmed trajectory is recorded in QC_report.csv and plotted to
    QC_flagged/, so the screen is auditable and the raw files stay untouched.

STAGE 2  Log-binned MSD
    Per-axis (x, y, z) and 3D MSD on log-spaced lag times, using the measured
    dt of each displacement pair. Construct-level curves are reported three
    ways: pair-count weighted mean (as in @msdanalyzer, Tinevez; Tarantino et
    al. J Cell Biol 2014), median across trajectories, and 10% trimmed mean.
    The weighted SD and the effective sample size
    N_freedom = (sum w)^2 / sum(w^2) follow the @msdanalyzer definitions.

STAGE 3  Model fits
    (a) Anomalous / short-lag:  MSD = 2*s_loc^2 + 2*D_a*t^alpha
    (b) Confined, Kusumi et al., Biophys J 65:2021-2040 (1993), 1D per axis:
            MSD(t) = 2*s_loc^2 + (L^2/6)*[1 - (96/pi^4) * S(t)]
            S(t)   = sum_{n odd} (1/n^4) * exp(-n^2 * pi^2 * D * t / L^2)
        so MSD(0) = 2*s_loc^2 and MSD(inf) = 2*s_loc^2 + L^2/6.
        NB the exponent is negative. Written as exp(+(1/2)(n*pi*sigma/L)^2 * t)
        with sigma^2 = 2D the series diverges; with the sign corrected it is
        identical to the expression above.
    (c) Ornstein-Uhlenbeck: MSD = 2*s_loc^2 + 2*sig_c^2*(1 - exp(-t/tau)),
        D = sig_c^2 / tau, equivalent corral size L = sqrt(12)*sig_c so that
        both models report the plateau on the same footing.

OUTPUTS (in MSD_v2/)
    QC_report.csv                  every trimmed / dropped / skipped trajectory
    QC_flagged/<trace>.png         diagnostic plot per trimmed trajectory
    <construct>_per_trace_MSD.csv  per-trajectory binned curves
    <construct>_MSD_curves.csv     construct-level curves
    MSD_fit_summary.csv            s_loc, alpha, D, L, plateau per construct
"""

import os
import glob
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore")

try:
    from tqdm import tqdm
except ImportError:                                   # graceful fallback
    def tqdm(it=None, **kw):
        return it if it is not None else None

# ============================== CONFIGURATION ==============================
DATA_DIRNAME = "data"          # folder holding the construct subfolders
RAW_ROOT = None                # set explicitly to override auto-location
OUT_DIRNAME = "MSD_v2"
N_WORKERS = 0                  # 0 = auto (cpu_count - 1, capped at 16)

# ---- Stage 1: segment-level QC -------------------------------------------
QC_ENABLED = True
QC_WINDOW = 101          # localizations in the centred rolling window (odd)
QC_DEV_NM = 150.0        # local median may deviate this far from trace median
QC_SD_FOLD = 3.0         # local SD may exceed baseline SD by this factor
QC_MIN_KEEP = 200        # drop the trajectory if fewer localizations survive
QC_MIN_KEEP_FRAC = 0.30  # ...or if less than this fraction survives
QC_AXES = ("z", "x", "y")   # axes screened; a flag on any one flags the point
SKIP_REF_SHEETS = True   # skip sheets whose name contains "(ref)"
QC_PLOTS = True          # write a diagnostic PNG for each trimmed trajectory

# ---- Stage 2: MSD binning -------------------------------------------------
N_BINS = 40
DT_MIN, DT_MAX = 3e-4, 3.0     # s
MIN_PAIRS_PER_TRACE_BIN = 50
MIN_TRACES_PER_BIN = 20

# ---- Stage 3: fitting -----------------------------------------------------
FIT_CURVE = "weighted_mean"    # "weighted_mean" | "median" | "trimmed_mean"
# Short-lag window for the anomalous exponent. 0.002 s reproduces the published
# Table S3 values; alpha_z falls monotonically as this is widened (0.82 at
# 2 ms -> 0.47 at 50 ms for HaloTag-Nup62-GFP), so state the window in Methods.
ALPHA_FIT_MAX_S = 0.002
KUSUMI_TERMS = 25              # odd terms in the series

# Short-time model used for the alpha / D columns of Table S3:
#   False -> MSD = 2*D*t^alpha              (no localization-error term; this is
#            the convention used in the published Table S3. The static offset is
#            absorbed into the power law, which pushes alpha toward 0.)
#   True  -> MSD = 2*s_loc^2 + 2*D*t^alpha  (alpha describes motion only)
# Both are always computed and written to MSD_fit_summary.csv; this switch only
# selects which one populates the formatted table.
ALPHA_INCLUDE_LOC_ERROR = False
TABLE_AXES = ("x", "y", "z")   # axis columns of the formatted table

# ---- plateau window -------------------------------------------------------
# The plateau is averaged from PLATEAU_MIN_S up to the lag where the number of
# contributing trajectories falls below PLATEAU_MIN_TRACE_FRAC of its maximum.
# Without the second cut, the sparse tail inflates L: for HaloTag-Nup62-GFP,
# including lags past 0.5 s (where only 5-23 trajectories remain) moves L_z
# from 87 nm -- the published value -- to 93 nm.
PLATEAU_MIN_S = 0.03
PLATEAU_MIN_TRACE_FRAC = 0.40

# ---- trajectory-level axial-spread screen ---------------------------------
# The segment screen above finds tracking failures that occupy part of a
# trajectory. It is blind to trajectories that wander axially throughout, which
# have no "good run" to keep -- e.g. HaloTag-Nup58-GFP tid 2029 / 3903 / 1740,
# with SD_z = 104 / 85 / 74 nm against a population median of 24 nm, and axial
# excursions 1.8-3.2x larger than lateral (so not stage drift). A trajectory is
# dropped if its SD on any screened axis exceeds EITHER threshold.
QC_TRACE_SD_ENABLED = True
QC_TRACE_SD_FOLD = 3.0          # x the construct median SD on that axis
QC_TRACE_SD_ABS_NM = 50.0       # or this absolute value
QC_TRACE_SD_AXES = ("z",)
# ===========================================================================

EDGES = np.logspace(np.log10(DT_MIN), np.log10(DT_MAX), N_BINS)
CENTERS = np.sqrt(EDGES[:-1] * EDGES[1:])
AXIS_COL = {"x": "Points_0", "y": "Points_1", "z": "Points_2"}


# ------------------------------------------------------------------ stage 1
def rolling_stats(v, window):
    s = pd.Series(v)
    minp = max(5, window // 4)
    med = s.rolling(window, center=True, min_periods=minp).median()
    sd = s.rolling(window, center=True, min_periods=minp).std()
    return med.to_numpy(), sd.to_numpy()


def longest_good_run(good):
    if not good.any():
        return 0, 0
    idx = np.flatnonzero(np.diff(np.concatenate(([0], good.view(np.int8), [0]))))
    starts, stops = idx[0::2], idx[1::2]
    k = int(np.argmax(stops - starts))
    return int(starts[k]), int(stops[k])


def qc_trajectory(df):
    """Return (start, stop, {axis: max local deviation in nm})."""
    n = len(df)
    window = min(QC_WINDOW, max(11, (n // 5) | 1))
    good = np.ones(n, dtype=bool)
    max_dev = {}
    for ax in QC_AXES:
        v = df[AXIS_COL[ax]].to_numpy() * 1e6                 # um
        trace_med = np.median(v)
        mad = np.median(np.abs(v - trace_med))
        baseline_sd = 1.4826 * mad if mad > 0 else np.std(v)
        rmed, rsd = rolling_stats(v, window)
        dev_nm = np.abs(rmed - trace_med) * 1000.0
        max_dev[ax] = float(np.nanmax(dev_nm)) if n else 0.0
        bad = (dev_nm > QC_DEV_NM) | (rsd > QC_SD_FOLD * baseline_sd)
        good &= ~np.nan_to_num(bad, nan=False).astype(bool)
    start, stop = longest_good_run(good)
    return start, stop, max_dev


def plot_qc(df, start, stop, construct, trace, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = os.path.join(out_dir, "QC_flagged")
    os.makedirs(d, exist_ok=True)
    t = (df["tim"] - df["tim"].iloc[0]).to_numpy()
    fig, axes = plt.subplots(3, 1, figsize=(7, 6), sharex=True)
    for ax, name in zip(axes, ("x", "y", "z")):
        v = df[AXIS_COL[name]].to_numpy() * 1e9
        ax.plot(t, v - np.median(v), lw=0.5, color="0.75")
        ax.plot(t[start:stop], v[start:stop] - np.median(v), lw=0.5, color="C0")
        ax.set_ylabel(f"{name} (nm)")
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("time (s)")
    axes[0].set_title(f"{construct}  {trace}   kept {start}:{stop} of {len(df)}"
                      "   (grey = removed)", fontsize=9)
    fig.tight_layout()
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in trace)
    fig.savefig(os.path.join(d, f"{construct}_{safe}.png"), dpi=120)
    plt.close(fig)


# ------------------------------------------------------------------ stage 2
def msd_one_trace(df):
    """Log-binned MSD for one trajectory; only ~N_BINS lags are evaluated."""
    x = df[AXIS_COL["x"]].to_numpy() * 1e6
    y = df[AXIS_COL["y"]].to_numpy() * 1e6
    z = df[AXIS_COL["z"]].to_numpy() * 1e6
    t = df["tim"].to_numpy()
    n = len(t)
    if n < 10:
        return pd.DataFrame()
    dt_typ = np.median(np.diff(t))
    if not np.isfinite(dt_typ) or dt_typ <= 0:
        return pd.DataFrame()

    offsets = np.unique(np.clip(np.round(CENTERS / dt_typ).astype(int), 1, n - 1))
    rows = []
    for off in offsets:
        i0 = np.arange(n - off)
        i1 = i0 + off
        frame = pd.DataFrame(dict(
            b=np.digitize(t[i1] - t[i0], EDGES),
            dx2=(x[i1] - x[i0]) ** 2,
            dy2=(y[i1] - y[i0]) ** 2,
            dz2=(z[i1] - z[i0]) ** 2))
        frame = frame[(frame.b >= 1) & (frame.b <= len(CENTERS))]
        if frame.empty:
            continue
        rows.append(frame.groupby("b").agg(
            x=("dx2", "mean"), y=("dy2", "mean"), z=("dz2", "mean"),
            w=("dx2", "size")).reset_index())
    if not rows:
        return pd.DataFrame()
    allrows = pd.concat(rows, ignore_index=True)
    for c in "xyz":
        allrows[c + "w"] = allrows[c] * allrows["w"]
    g = allrows.groupby("b")[["xw", "yw", "zw", "w"]].sum()
    for c in "xyz":
        g[c] = g[c + "w"] / g["w"]
    g["msd3d"] = g["x"] + g["y"] + g["z"]
    g = g[g["w"] >= MIN_PAIRS_PER_TRACE_BIN]
    return g.reset_index()[["b", "x", "y", "z", "msd3d", "w"]]


def process_file(args):
    """Worker: QC + MSD for every trajectory in one .xlsx. Returns small data."""
    construct, path, out_dir = args
    fname = os.path.basename(path)
    parts = fname.replace(".xlsx", "").split()
    if len(parts) >= 2 and parts[-2].isdigit():
        meas = f"M{parts[-2].zfill(2)}_{parts[-1].zfill(2)}"
    else:
        meas = f"M{parts[-1].zfill(2)}"

    records, msd_rows, stats = [], [], []
    try:
        sheets = pd.read_excel(path, sheet_name=None)
    except Exception as exc:
        return construct, pd.DataFrame(), [dict(
            construct=construct, file=fname, sheet="", trace="", n_total=0,
            n_kept=0, n_removed=0, action=f"read_error: {exc}")], []

    for sheet, df in sheets.items():
        base = dict(construct=construct, file=fname, sheet=sheet)
        if SKIP_REF_SHEETS and "(ref)" in sheet:
            records.append(dict(base, trace=f"{meas}_{sheet}", n_total=len(df),
                                n_kept=0, n_removed=len(df),
                                action="skipped_reference_sheet"))
            continue
        if not all(c in df.columns for c in AXIS_COL.values()) or "tim" not in df:
            continue
        tid = sheet.split(" ")[-1] if " " in sheet else sheet
        trace = f"{meas}_tid_{tid}"
        n = len(df)
        if n < QC_MIN_KEEP:
            records.append(dict(base, trace=trace, n_total=n, n_kept=0,
                                n_removed=n, action="dropped_too_short"))
            continue

        if QC_ENABLED:
            start, stop, max_dev = qc_trajectory(df)
        else:
            start, stop, max_dev = 0, n, {a: 0.0 for a in QC_AXES}
        devs = {f"max_dev_{a}_nm": max_dev[a] for a in QC_AXES}
        kept = stop - start

        if kept < QC_MIN_KEEP or kept < QC_MIN_KEEP_FRAC * n:
            records.append(dict(base, trace=trace, n_total=n, n_kept=0,
                                n_removed=n, action="dropped_qc", **devs))
            continue
        if kept < n:
            records.append(dict(base, trace=trace, n_total=n, n_kept=kept,
                                n_removed=n - kept, action="trimmed",
                                keep_start=start, keep_stop=stop, **devs))
            if QC_PLOTS:
                try:
                    plot_qc(df, start, stop, construct, trace, out_dir)
                except Exception:
                    pass

        kept_df = df.iloc[start:stop]
        stats.append(dict(
            construct=construct, file=fname, trace=trace, n_kept=kept,
            **{f"sd_{a}_nm": float(kept_df[AXIS_COL[a]].std() * 1e9)
               for a in ("x", "y", "z")}))

        m = msd_one_trace(kept_df)
        if not m.empty:
            m["trace"] = trace
            msd_rows.append(m)

    per_trace = pd.concat(msd_rows, ignore_index=True) if msd_rows \
        else pd.DataFrame()
    return construct, per_trace, records, stats


def weighted_stats(values, weights):
    """@msdanalyzer definitions: reliability-weighted SD and N_freedom."""
    v = np.asarray(values, float)
    w = np.asarray(weights, float)
    sw, ssw = w.sum(), (w ** 2).sum()
    mean = np.sum(w * v) / sw
    denom = sw ** 2 - ssw
    sd = np.sqrt(sw / denom * np.sum(w * (v - mean) ** 2)) if denom > 0 else np.nan
    nfree = sw ** 2 / ssw
    return mean, sd, (sd / np.sqrt(nfree) if nfree > 0 else np.nan), nfree


def trimmed_mean(v, frac=0.10):
    v = np.sort(np.asarray(v, float))
    k = int(len(v) * frac)
    return v[k:len(v) - k].mean() if len(v) - 2 * k > 0 else v.mean()


def summarise_construct(per_trace):
    out = []
    for b, grp in per_trace.groupby("b"):
        if grp["trace"].nunique() < MIN_TRACES_PER_BIN:
            continue
        row = {"dt": CENTERS[b - 1], "n_traces": grp["trace"].nunique(),
               "n_pairs": grp["w"].sum()}
        for ax in ("x", "y", "z", "msd3d"):
            m, sd, sem, nf = weighted_stats(grp[ax], grp["w"])
            row[f"{ax}_weighted_mean"] = m
            row[f"{ax}_weighted_sd"] = sd
            row[f"{ax}_sem"] = sem
            row[f"{ax}_median"] = grp[ax].median()
            row[f"{ax}_trimmed_mean"] = trimmed_mean(grp[ax])
            if ax == "z":
                row["N_freedom"] = nf
        out.append(row)
    return pd.DataFrame(out).sort_values("dt").reset_index(drop=True)


# ------------------------------------------------------------------ stage 3
def m_anomalous(t, s_loc, D_a, alpha):
    return 2 * s_loc ** 2 + 2 * D_a * t ** alpha


def m_anomalous_nooffset(t, D_a, alpha):
    return 2 * D_a * t ** alpha


def m_kusumi(t, s_loc, L, D):
    n = np.arange(1, 2 * KUSUMI_TERMS, 2)[:, None]
    S = np.sum(n ** -4.0 * np.exp(-(n ** 2) * np.pi ** 2 * D * t[None, :] / L ** 2),
               axis=0)
    return 2 * s_loc ** 2 + (L ** 2 / 6.0) * (1 - (96 / np.pi ** 4) * S)


def m_ou(t, s_loc, sig_c, tau):
    return 2 * s_loc ** 2 + 2 * sig_c ** 2 * (1 - np.exp(-t / tau))


def r_squared(y, yhat):
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else np.nan


def _se(pcov):
    """Standard error of each fitted parameter, from the covariance matrix."""
    with np.errstate(invalid="ignore"):
        return np.sqrt(np.diag(pcov))


def fit_loglog(t, y, ysem=None):
    """log10(MSD) = b + a*log10(t).  alpha = a,  D = 10^b / 2.

    This is the convention used in the original Origin analysis
    (Ver 3/Summary_MSD.xlsx stores Intercept / Slope and their errors).
    """
    x, ly = np.log10(t), np.log10(y)
    X = np.vstack([np.ones(len(x)), x]).T
    if ysem is not None and np.all(np.isfinite(ysem)) and np.all(ysem > 0):
        w = 1.0 / (ysem / y / np.log(10)) ** 2          # error on log10(y)
    else:
        w = np.ones(len(x))
    W = np.diag(w)
    cov = np.linalg.inv(X.T @ W @ X)
    beta = cov @ (X.T @ W @ ly)
    resid = ly - X @ beta
    dof = max(1, len(x) - 2)
    scale = float(resid.T @ W @ resid) / dof            # scale by reduced chi2
    cov = cov * scale
    se = np.sqrt(np.diag(cov))
    D = 10 ** beta[0] / 2.0
    return dict(intercept=beta[0], intercept_se=se[0],
                alpha_loglog=beta[1], alpha_loglog_se=se[1],
                D_loglog_um2_s=D, D_loglog_se=np.log(10) * D * se[0],
                R2_loglog=r_squared(ly, X @ beta))


def fit_plateau(dt, msd, sigma=None, ntraces=None):
    """Mean MSD over the flat region, and L = sqrt(6 * plateau)."""
    ok = dt >= PLATEAU_MIN_S
    if ntraces is not None and len(ntraces) == len(dt):
        n = np.asarray(ntraces, float)
        ok = ok & (n >= PLATEAU_MIN_TRACE_FRAC * np.nanmax(n))
    if ok.sum() < 3:
        ok = dt >= PLATEAU_MIN_S
    if ok.sum() < 2:
        return {}
    P = float(np.mean(msd[ok]))
    dP = (float(np.median(sigma[ok])) if sigma is not None
          else float(np.std(msd[ok], ddof=1) / np.sqrt(ok.sum())))
    L = np.sqrt(6 * P)
    return dict(plateau_fit_um2=P, plateau_fit_se=dP,
                L_plateau_nm=L * 1000,
                L_plateau_se_nm=0.5 * np.sqrt(6 / P) * dP * 1000,
                plateau_dt_lo=float(dt[ok].min()),
                plateau_dt_hi=float(dt[ok].max()),
                plateau_n_bins=int(ok.sum()))


def fit_axis(dt, msd, label, sigma=None, ntraces=None):
    res = {"axis": label, "MSD_initial_um2": msd[0],
           "plateau_empirical_um2": np.median(msd[dt >= 0.1])
           if (dt >= 0.1).any() else np.nan,
           "n_points_fit": len(dt)}
    res.update(fit_plateau(dt, msd, sigma, ntraces))
    kw = dict(maxfev=40000)
    if sigma is not None:
        kw.update(sigma=sigma, absolute_sigma=False)

    short = dt <= ALPHA_FIT_MAX_S
    res["n_points_alpha_fit"] = int(short.sum())
    skw = dict(kw)
    if sigma is not None:
        skw["sigma"] = sigma[short]

    # (a) short-time power law, no localization-error term (published convention)
    if short.sum() >= 3:
        try:
            p, c = curve_fit(m_anomalous_nooffset, dt[short], msd[short],
                             p0=[0.01, 0.5], bounds=([1e-10, 0.01], [1e3, 2.0]),
                             **skw)
            e = _se(c)
            res.update(alpha_nooffset=p[1], alpha_nooffset_se=e[1],
                       D_nooffset_um2_s=p[0], D_nooffset_se=e[0],
                       R2_alpha_nooffset=r_squared(
                           msd[short], m_anomalous_nooffset(dt[short], *p)))
        except Exception:
            pass
    # (b) short-time power law with localization-error offset
    if short.sum() >= 4:
        try:
            p, c = curve_fit(m_anomalous, dt[short], msd[short],
                             p0=[0.01, 0.01, 0.5],
                             bounds=([0, 1e-10, 0.01], [0.2, 1e3, 2.0]), **skw)
            e = _se(c)
            res.update(alpha_offset=p[2], alpha_offset_se=e[2],
                       D_offset_um2_s=p[1], D_offset_se=e[1],
                       s_loc_short_nm=p[0] * 1000, s_loc_short_se_nm=e[0] * 1000,
                       R2_alpha_offset=r_squared(
                           msd[short], m_anomalous(dt[short], *p)))
        except Exception:
            pass

    # (c) log-log linear regression -- the original Origin convention
    if short.sum() >= 3:
        try:
            res.update(fit_loglog(dt[short], msd[short],
                                  sigma[short] if sigma is not None else None))
        except Exception:
            pass

    tag = "offset" if ALPHA_INCLUDE_LOC_ERROR else "nooffset"
    res["alpha"] = res.get(f"alpha_{tag}", np.nan)
    res["alpha_se"] = res.get(f"alpha_{tag}_se", np.nan)
    res["D_alpha_um2_s"] = res.get(f"D_{tag}_um2_s", np.nan)
    res["D_alpha_se"] = res.get(f"D_{tag}_se", np.nan)

    plateau0 = res["plateau_empirical_um2"]
    L0 = np.sqrt(6 * plateau0) if np.isfinite(plateau0) and plateau0 > 0 else 0.05
    try:
        p, c = curve_fit(m_kusumi, dt, msd, p0=[0.015, L0, 0.05],
                         bounds=([0, 1e-3, 1e-6], [0.2, 5.0, 1e3]), **kw)
        e = _se(c)
        res.update(kusumi_s_loc_nm=p[0] * 1000, kusumi_s_loc_se_nm=e[0] * 1000,
                   kusumi_L_nm=p[1] * 1000, kusumi_L_se_nm=e[1] * 1000,
                   kusumi_D_um2_s=p[2], kusumi_D_se=e[2],
                   kusumi_plateau_um2=p[1] ** 2 / 6,
                   R2_kusumi=r_squared(msd, m_kusumi(dt, *p)))
    except Exception:
        res.update(kusumi_L_nm=np.nan, kusumi_D_um2_s=np.nan)
    try:
        p, c = curve_fit(m_ou, dt, msd, p0=[0.015, L0 / 3.46, 0.02],
                         bounds=([0, 1e-4, 1e-5], [0.2, 2.0, 1e3]), **kw)
        e = _se(c)
        res.update(ou_s_loc_nm=p[0] * 1000, ou_sigma_c_nm=p[1] * 1000,
                   ou_tau_s=p[2], ou_tau_se=e[2],
                   ou_D_um2_s=p[1] ** 2 / p[2],
                   ou_L_equiv_nm=np.sqrt(12) * p[1] * 1000,
                   ou_L_equiv_se_nm=np.sqrt(12) * e[1] * 1000,
                   ou_plateau_um2=2 * p[1] ** 2,
                   R2_ou=r_squared(msd, m_ou(dt, *p)))
    except Exception:
        res.update(ou_tau_s=np.nan, ou_D_um2_s=np.nan)
    return res


# ------------------------------------------------- formatted Table S3 output
def _fmt_plain(v, se, nd=2):
    if not np.isfinite(v):
        return "n.d."
    return f"{v:.{nd}f} ± {se:.{nd}f}" if np.isfinite(se) else f"{v:.{nd}f}"


def _fmt_sci(v, se):
    """1.81 ± 0.06 x10-4, matching the published table."""
    if not np.isfinite(v) or v <= 0:
        return "n.d."
    ex = int(np.floor(np.log10(abs(v))))
    m, s = v / 10 ** ex, (se / 10 ** ex if np.isfinite(se) else np.nan)
    return (f"{m:.2f} ± {s:.2f} x10{ex}" if np.isfinite(s)
            else f"{m:.2f} x10{ex}")


def to_markdown_simple(df):
    """Pipe-table markdown without the optional `tabulate` dependency."""
    cols = [str(c) for c in df.columns]
    rows = [[("" if pd.isna(v) else str(v)) for v in r]
            for r in df.itertuples(index=False)]
    width = [max(len(cols[i]), *(len(r[i]) for r in rows)) if rows
             else len(cols[i]) for i in range(len(cols))]
    out = ["| " + " | ".join(c.ljust(width[i]) for i, c in enumerate(cols)) + " |",
           "|" + "|".join("-" * (w + 2) for w in width) + "|"]
    out += ["| " + " | ".join(v.ljust(width[i]) for i, v in enumerate(r)) + " |"
            for r in rows]
    return "\n".join(out)


def build_table(fit_df, method="loglog", l_source="plateau"):
    """Wide Table S3: rows = constructs, columns = alpha/D/L per axis.

    method    "loglog" -> alpha/D from the log-log linear regression
                          (original Origin convention)
              "direct" -> alpha/D from the nonlinear 2*D*t^alpha fit
    l_source  "plateau" -> L = sqrt(6 * plateau)   (original convention)
              "kusumi"  -> L from the Kusumi confined-diffusion fit
    """
    a_k, ase_k, d_k, dse_k = (
        ("alpha_loglog", "alpha_loglog_se", "D_loglog_um2_s", "D_loglog_se")
        if method == "loglog" else
        ("alpha", "alpha_se", "D_alpha_um2_s", "D_alpha_se"))
    l_k, lse_k = (("L_plateau_nm", "L_plateau_se_nm") if l_source == "plateau"
                  else ("kusumi_L_nm", "kusumi_L_se_nm"))
    rows_fmt, rows_raw = [], []
    for construct, grp in fit_df.groupby("construct", sort=False):
        fmt = {"construct": construct}
        raw = {"construct": construct}
        for ax in TABLE_AXES:
            r = grp[grp.axis == ax]
            if r.empty:
                continue
            r = r.iloc[0]
            fmt[f"alpha_{ax}"] = _fmt_plain(r.get(a_k), r.get(ase_k))
            fmt[f"D_{ax}_um2_s"] = _fmt_sci(r.get(d_k), r.get(dse_k))
            fmt[f"L_{ax}_nm"] = _fmt_plain(r.get(l_k), r.get(lse_k), nd=0)
            for out, src in (("alpha", a_k), ("alpha_se", ase_k),
                             ("D", d_k), ("D_se", dse_k),
                             ("L_nm", l_k), ("L_se_nm", lse_k)):
                raw[f"{out}_{ax}"] = r.get(src)
        rows_fmt.append(fmt)
        rows_raw.append(raw)
    fmt_df = pd.DataFrame(rows_fmt)
    fmt_df = fmt_df.reindex(columns=[c for c in
                                     ["construct"] +
                                     [x for ax in TABLE_AXES for x in
                                      (f"alpha_{ax}", f"D_{ax}_um2_s",
                                       f"L_{ax}_nm")]
                                     if c in fmt_df.columns])
    return fmt_df, pd.DataFrame(rows_raw)


def fit_construct(curves, construct):
    rows = []
    for ax in ("x", "y", "z", "msd3d"):
        col, semcol = f"{ax}_{FIT_CURVE}", f"{ax}_sem"
        cols = ["dt", col] + ([semcol] if semcol in curves.columns else [])
        d = curves[cols].dropna(subset=["dt", col])
        if len(d) < 6:
            continue
        sig = None
        if semcol in d.columns:
            s = d[semcol].to_numpy(float)
            if np.isfinite(s).all() and (s > 0).all():
                sig = s
        ntr = (curves.loc[d.index, "n_traces"].to_numpy(float)
               if "n_traces" in curves.columns else None)
        r = fit_axis(d["dt"].to_numpy(), d[col].to_numpy(), ax,
                     sigma=sig, ntraces=ntr)
        r["construct"] = construct
        r["curve"] = FIT_CURVE
        rows.append(r)
    return rows


# ---------------------------------------------------------------------- main
def locate_raw_root():
    if RAW_ROOT:
        return RAW_ROOT
    here = os.path.dirname(os.path.abspath(__file__))
    for base in (os.getcwd(), here):
        cand = os.path.join(base, DATA_DIRNAME)
        if glob.glob(os.path.join(cand, "*", "*.xlsx")):
            return cand
    for base in (os.getcwd(), here):
        if glob.glob(os.path.join(base, "*", "*.xlsx")):
            print(f"No '{DATA_DIRNAME}' folder found; scanning {base} directly.")
            return base
    raise SystemExit(
        f"Could not find any construct folders. Expected "
        f"{DATA_DIRNAME}/<construct>/*.xlsx under {os.getcwd()} or {here}.")


def main():
    root = locate_raw_root()
    parent = os.path.dirname(root) if os.path.basename(root) == DATA_DIRNAME \
        else root
    out_dir = os.path.join(parent, OUT_DIRNAME)
    os.makedirs(out_dir, exist_ok=True)

    workers = N_WORKERS or min(16, max(1, (os.cpu_count() or 4) - 1))
    constructs = sorted(
        d for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d))
        and glob.glob(os.path.join(root, d, "*.xlsx")))
    jobs = [(c, p, out_dir)
            for c in constructs
            for p in sorted(glob.glob(os.path.join(root, c, "*.xlsx")))]

    print(f"Raw data   : {root}")
    print(f"Output     : {out_dir}")
    print(f"Constructs : {len(constructs)}  ({', '.join(constructs)})")
    print(f"Files      : {len(jobs)}   workers: {workers}\n")

    collected = {c: [] for c in constructs}
    records, all_stats = [], []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(process_file, j): j for j in jobs}
        bar = tqdm(total=len(futures), desc="QC + MSD", unit="file")
        for fut in as_completed(futures):
            construct, per_trace, recs, stats = fut.result()
            if per_trace is not None and len(per_trace):
                collected[construct].append(per_trace)
            records.extend(recs)
            all_stats.extend(stats)
            if bar is not None:
                bar.set_postfix_str(construct[:24])
                bar.update(1)
        if bar is not None:
            bar.close()

    # ---- trajectory-level axial-spread screen (needs all traces of a construct)
    stats_df = pd.DataFrame(all_stats)
    dropped = set()
    if QC_TRACE_SD_ENABLED and len(stats_df):
        stats_df["excluded"] = False
        stats_df["exclude_reason"] = ""
        for construct, grp in stats_df.groupby("construct"):
            for ax in QC_TRACE_SD_AXES:
                col = f"sd_{ax}_nm"
                med = grp[col].median()
                thr = min(QC_TRACE_SD_FOLD * med, QC_TRACE_SD_ABS_NM)
                bad = grp.index[grp[col] > thr]
                print(f"  {construct}: median SD_{ax} = {med:.0f} nm, "
                      f"threshold = {thr:.0f} nm, "
                      f"{len(bad)} of {len(grp)} trajectories excluded")
                for i in bad:
                    r = stats_df.loc[i]
                    stats_df.loc[i, "excluded"] = True
                    stats_df.loc[i, "exclude_reason"] = (
                        f"sd_{ax}={r[col]:.0f}nm > {thr:.0f}nm")
                    dropped.add((r["construct"], r["trace"]))
                    records.append(dict(
                        construct=r["construct"], file=r["file"], sheet="",
                        trace=r["trace"], n_total=r["n_kept"], n_kept=0,
                        n_removed=r["n_kept"], action="dropped_trace_sd",
                        **{f"max_dev_{a}_nm": np.nan for a in QC_AXES}))
        stats_df.to_csv(os.path.join(out_dir, "QC_trace_stats.csv"), index=False)
        for c in collected:
            collected[c] = [
                d[~d["trace"].isin({t for cc, t in dropped if cc == c})]
                for d in collected[c]]

    qc_cols = ["construct", "file", "sheet", "trace", "n_total", "n_kept",
               "n_removed", "action", "keep_start", "keep_stop"] + \
              [f"max_dev_{a}_nm" for a in QC_AXES]
    qc_df = pd.DataFrame(records) if records else pd.DataFrame(columns=qc_cols)
    qc_df = qc_df.reindex(columns=[c for c in qc_cols if c in qc_df.columns]
                          + [c for c in qc_df.columns if c not in qc_cols])
    qc_df.to_csv(os.path.join(out_dir, "QC_report.csv"), index=False)

    all_fits = []
    for construct in tqdm(constructs, desc="Summarise + fit", unit="construct"):
        if not collected[construct]:
            print(f"  {construct}: no usable trajectories")
            continue
        per_trace = pd.concat(collected[construct], ignore_index=True)
        per_trace["dt"] = CENTERS[per_trace["b"] - 1]
        per_trace.to_csv(
            os.path.join(out_dir, f"{construct}_per_trace_MSD.csv"), index=False)
        curves = summarise_construct(per_trace)
        curves.to_csv(
            os.path.join(out_dir, f"{construct}_MSD_curves.csv"), index=False)
        all_fits.extend(fit_construct(curves, construct))

    if all_fits:
        lead = ["construct", "axis", "curve", "n_points_fit",
                "n_points_alpha_fit", "MSD_initial_um2",
                "alpha_loglog", "alpha_loglog_se", "D_loglog_um2_s",
                "D_loglog_se", "R2_loglog", "intercept", "intercept_se",
                "alpha", "alpha_se", "D_alpha_um2_s", "D_alpha_se",
                "plateau_fit_um2", "plateau_fit_se", "L_plateau_nm",
                "L_plateau_se_nm", "plateau_dt_lo", "plateau_dt_hi",
                "plateau_n_bins", "plateau_empirical_um2",
                "kusumi_L_nm", "kusumi_L_se_nm", "kusumi_D_um2_s",
                "kusumi_D_se", "kusumi_s_loc_nm", "kusumi_plateau_um2",
                "R2_kusumi"]
        df = pd.DataFrame(all_fits)
        df = df[[c for c in lead if c in df.columns]
                + [c for c in df.columns if c not in lead]]
        df.to_csv(os.path.join(out_dir, "MSD_fit_summary.csv"), index=False)

        desc = {
            "loglog": "alpha = slope of log10(MSD) vs log10(t), D = 10^intercept/2"
                      "  [original Origin convention]",
            "direct": "alpha, D from nonlinear fit of 2*D*t^alpha",
        }
        for method in ("loglog", "direct"):
            fmt_df, raw_df = build_table(df, method=method, l_source="plateau")
            stem = f"MSD_Table_S3_{method}"
            fmt_df.to_csv(os.path.join(out_dir, stem + ".csv"), index=False)
            raw_df.to_csv(os.path.join(out_dir, stem + "_raw.csv"), index=False)
            try:
                with open(os.path.join(out_dir, stem + ".md"), "w",
                          encoding="utf-8") as fh:
                    fh.write(to_markdown_simple(fmt_df))
            except Exception as exc:
                print(f"  (markdown copy skipped: {exc})")
            print(f"\nTable S3 [{method}]  {desc[method]}")
            print(f"  alpha/D fit to dt <= {ALPHA_FIT_MAX_S} s;  "
                  f"L = sqrt(6*plateau) over dt >= {PLATEAU_MIN_S} s "
                  f"with n_traces >= {PLATEAU_MIN_TRACE_FRAC:.0%} of max")
            print(fmt_df.to_string(index=False))

        print(f"\nFits   -> {os.path.join(out_dir, 'MSD_fit_summary.csv')}")
        print(f"Tables -> {out_dir}\\MSD_Table_S3_{{loglog,direct}}.csv "
              f"(+ .md, _raw.csv)")

    n_trim = sum(1 for r in records if r.get("action") == "trimmed")
    n_sd = sum(1 for r in records if r.get("action") == "dropped_trace_sd")
    n_drop = sum(1 for r in records
                 if str(r.get("action", "")).startswith("dropped")) - n_sd
    n_skip = sum(1 for r in records if r.get("action") == "skipped_reference_sheet")
    print(f"QC: {n_trim} segment-trimmed, {n_sd} dropped on trajectory SD, "
          f"{n_drop} dropped otherwise, {n_skip} reference sheets skipped")
    print(f"    -> {os.path.join(out_dir, 'QC_report.csv')}, QC_trace_stats.csv")


if __name__ == "__main__":
    main()
