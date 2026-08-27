# Supplementary Figure 12 — Confined-diffusion plateau analysis

Fits a confined-diffusion model to the MSD of FG-Nup single-molecule
trajectories and reports the **plateau** — a sample-size-independent measure of
confinement — in place of the "max extension" metric, whose value grows with the
number of trajectories analysed.

Prepared in response to Reviewer 1, comment 11.

---

## 1. Folder layout

```
Supplementary Figure 12/
    confined_diffusion_plateau.py   <- the analysis (edit the SETTINGS block)
    README.md                       <- this file

    Supplemental MSD FG/            <- 13 FG-Nup trajectory sets
    Supplemental MSD NonFG/         <-  9 non-FG trajectory sets
    Confinement_Analysis_Output/    <- created by the script
```

The two data folders are named in the `DATASETS` dict at the top of the script.
Everything is resolved relative to the script's own location, so the folder can
be moved as a unit.

---

## 2. Requirements

- Python 3.9+
- `numpy`, `pandas`, `scipy`, `matplotlib`
- `tqdm` — optional; a built-in fallback prints coarse progress if it is absent

```
pip install numpy pandas scipy matplotlib tqdm
```

---

## 3. Input format

The script reads every `*_Trajectory.csv` in each dataset folder. Files whose
name contains `_Joe` are skipped (an alternative tracking of Nup214N kept for
comparison).

Each file is headerless, with the first two columns being **x and y in nm**.
Individual traces are separated by rows that are blank in both columns; a trace
needs at least 3 points to be used. The protein name is the file stem with
`_Trajectory.csv` removed, and that name is the key into the `SIGMA_LOC` table.

The other CSVs sitting alongside (`*_alpha.csv`, `*_gamma.csv`,
`*_inditrack.csv`, `*_sumtrack.csv`, `*_Extension.csv/.xlsx`) are msdanalyzer
exports. They are **not** read by this script and are kept as the upstream
record.

---

## 4. Method

For each protein, and separately for x and y:

1. **Re-extract a 1D MSD.** msdanalyzer exports a combined 2D ⟨Δr²⟩, but both
   confinement models below are 1D, so the time-averaged, ensemble-pooled MSD is
   recomputed per dimension from the raw traces.

2. **Fit two confined-diffusion models** over the well-sampled lag range only:

   | Model | Equation | Plateau |
   |---|---|---|
   | OU / harmonic | MSD = P·(1 − e^(−t/τ)) + 2σ_loc² | P = 2σ² |
   | Kusumi box | MSD = L²/6 − (16L²/π⁴)·Σ_{n odd} n⁻⁴·e^(−n²π²Dt/L²) + 2σ_loc² | L²/6 |

   Kusumi box is Eq. 11 of Kusumi et al. 1993, *Biophys J* **65**:2021. The
   localization-error offset 2σ_loc² is **fixed** from Supplementary Table S3,
   not fitted — see `SIGMA_LOC` in the script.

3. **Bootstrap** over whole traces (1000 resamples) for 95% confidence
   intervals.

4. **Combine x and y** into an RMS radial extension, √((P_x + P_y)/2), and
   tabulate it against the older Gaussian-fit "max extension" mean.

---

## 5. Settings

All in the `SETTINGS` block near the top of the script.

| Setting | Default | Meaning |
|---|---|---|
| `DATASETS` | FG, NonFG | label → data subfolder |
| `DT` | 0.05 s | frame interval |
| `MIN_TRACES_PER_LAG` | 5 | drop lags supported by fewer traces than this |
| `MAX_LAG_S` | 0.30 s | hard cap on the lag used in the fit |
| `MIN_LAGS_TO_FIT` | 4 | minimum lag points needed to attempt a fit |
| `N_BOOT` | 1000 | bootstrap resamples |
| `RNG_SEED` | 0 | fixed, for reproducibility |
| `KUSUMI_N_ODD` | 25 | odd terms in the Kusumi series |
| `SIGMA_LOC` | per protein | localization precision (nm), Table S3 |
| `DEFAULT_SIGMA_LOC` | 8.0 nm | fallback when a protein is not in the table |
| `SAVE_FIGURES` | True | write the diagnostic PNGs |

The fit range is the crux of the analysis: the far plateau is thinly sampled, so
`MIN_TRACES_PER_LAG` and `MAX_LAG_S` together decide which lags are trustworthy.
At `DT` = 0.05 s, `MAX_LAG_S` = 0.30 s admits at most 6 lag points.

---

## 6. Running it

```
python confined_diffusion_plateau.py
```

Both datasets are processed in one run. Outputs land in
`Confinement_Analysis_Output/`:

| Output | Contents |
|---|---|
| `confinement_results_ALL.csv` | every fitted parameter, per protein and dimension |
| `confinement_summary_table.csv` | the three extension metrics side by side |
| `SUMMARY_three_metrics.png` | comparison bar chart across all proteins |
| `FG_plateau_fits/*.png` | per-protein MSD fits, x and y panels |
| `NonFG_plateau_fits/*.png` | same, non-FG proteins |

---

## 7. Result as committed

22 proteins (13 FG, 9 non-FG), 44–51 traces each. RMS extension from the OU
plateau spans 12.2 nm (POM121N) to 62.6 nm (TPRC). Across every protein the
plateau estimate is roughly half the old max-extension mean, which is the point
of the comment: max extension keeps growing with sample size, the plateau does
not.