# Supplementary Figure 30 — Bimodal-recovery validation of the 2D → 3D transformation

Tests whether the 2D → 3D transformation reliably recovers a **two-component**
radial distribution, and — just as important — whether it invents a second peak
when the truth has only one.

Prepared in response to the round-2 reviewer concern that the published Monte
Carlo validated only the recovery of a single peak position, not the full shape
of the 3D probability distribution. The open/closed-state claims for Nup214C and
Nup153C rest on resolving two distinct radial peaks within one population, so
two-component recovery is what has to be validated.

---

## 1. Folder layout

```
Supplementary Figure 30/
    bimodal_resampling_simulation.py   <- main validation (run first)
    reference_density.py               <- transform-independent ground truth
    sweep_simulation.py                <- ratio sweep, simulation half (heavy)
    sweep_plotting.py                  <- ratio sweep, plotting half (cheap)
    README.md                          <- this file

    Step_2_2D_to_3D_transformation.py  }  the transform under test —
    config.py                          }  copied from ../2D to 3D Transformation/
    data_io.py                         }  see section 6

    bimodal_resampling_output/         <- created by scripts 1 and 2
    sweep_output/                      <- created by the sweep (not committed)
```

---

## 2. Requirements

- Python 3.9+
- `numpy`, `pandas`, `scipy`, `matplotlib`
- `tqdm` — optional; every script has a built-in stderr fallback

```
pip install numpy pandas scipy matplotlib tqdm
```

All scripts use the `Agg` backend and save figures without opening windows.

---

## 3. What the validation does

There is no experimental input. Everything is simulated against a known
ground truth, and the estimator under test is the real one.

**Pool.** Build a "region" of localizations: `N1` peripheral points on a ring of
radius `R`, plus `N2` central points at `R_central`. Each point is projected to
its Y coordinate and broadened by the localization precision. The ratio N1:N2
sets the relative weight of the two states.

**Resample.** Draw `n_sample` points without replacement, run the **real** 2D→3D
transform, and read the single peak position (argmax of the deprojected density,
with parabolic sub-bin refinement). Repeat `n_resample` times.

Because a random subsample is by chance dominated by one component or the other,
the single recovered peak lands near 0 or near R, and the histogram of peak
positions is bimodal — the behaviour reported in the paper. Regenerating the
pool `n_pool_reps` times confirms this is not an artefact of one random pool.

**Negative controls.** The same procedure runs on two unimodal ground truths at
the same total count: all-peripheral and all-central. If those stay unimodal
while the mixture is bimodal, bimodality is diagnostic of a genuine two-component
truth rather than a product of the resampling. This is the control that answers
the reviewer's strongest objection.

**Faithfulness.** The transform is not reimplemented. `invert_2d_to_3d` is
imported from `Step_2_2D_to_3D_transformation.py`, so simulated data passes
through exactly the code used on real data.

**Reference density** (`reference_density.py`) builds a transform-independent
ground truth: a huge population of *true* radial coordinates r = √(Y²+Z²),
histogrammed and divided by the exact ring area π(r_hi² − r_lo²) to give an areal
density. Because it uses no transform, the comparison cannot be circular.

> Compare **positions**, not amplitudes. The reference is an areal density; the
> resample histogram is the distribution of a winner-take-all argmax. The taller
> central peak wins disproportionately, so the histogram under-represents the
> peripheral peak relative to its true areal weight. Matching modes validate
> two-peak resolvability; matching heights are not expected.

**Ratio sweep** (`sweep_simulation.py` → `sweep_plotting.py`) repeats all of the
above across the full peripheral:central range at a fixed pool size of 200,
sampled finely near the crossover (N1 = 200 … 160 in steps of 5) and coarsely
elsewhere (150 … 0 in steps of 10). Simulation and plotting are split so figures
can be restyled without re-running the Monte Carlo.

---

## 4. Configuration

Physical parameters live in one place — the `CFG` dict at the top of
`bimodal_resampling_simulation.py` — and the other scripts import it.

| Key | Default | Meaning |
|---|---|---|
| `R` | 30 nm | peripheral ring radius |
| `R_central` | 0.0 nm | central-component radial position |
| `precision` | 7.0 nm | localization precision, σ |
| `N1` | 170 | peripheral points in the pool |
| `N2` | 30 | central points in the pool |
| `n_sample` | 100 | points drawn per subsample, without replacement |
| `n_resample` | 10000 | subsamples per pool |
| `n_pool_reps` | 10 | independent pool regenerations |
| `bin_size` | 10.0 nm | R bin size passed to the transform |
| `R_upper` | 150.0 nm | radial upper limit passed to the transform |
| `refine_peak` | True | parabolic sub-bin interpolation of the argmax |
| `run_controls` | True | also run the two unimodal controls |
| `seed` | 20260624 | fixed, for reproducibility |

Script-local knobs: `REF_SCALE` (20000), `REF_BIN` (0.5 nm), `HIST_BIN`
(10.0 nm) in `reference_density.py`; `POOL`, `FINE_N1`, `COARSE_N1`,
`SWEEP_N_RESAMPLE`, `SWEEP_N_POOL_REPS` in `sweep_simulation.py`; `STYLE`,
`SQUARE`, `PANEL` in `sweep_plotting.py`.

To test quickly, drop `n_resample` to ~30 and `n_pool_reps` to ~2.

---

## 5. Running it

In order:

```
python bimodal_resampling_simulation.py    # writes peaks_raw.csv — needed by the next step
python reference_density.py                # reads those peaks, adds the ground truth
python sweep_simulation.py                 # optional, heavy: the full ratio sweep
python sweep_plotting.py                   # optional: figures from the sweep CSVs
```

`bimodal_resampling_output/`:

| Output | Contents |
|---|---|
| `peaks_raw.csv` | one row per subsample: condition, pool_rep, peak position |
| `summary.csv` | per condition: mode positions, band fractions, input ratio |
| `histograms.png` | overlaid peak-position histograms, mixture vs both controls |
| `per_pool.png` | mixture histogram split by pool regeneration (RNG stability) |
| `reference_density.csv` | true areal density, fine bins, for OriginLab |
| `resample_histogram.csv` | binned resample peaks |
| `reference_vs_resample.png` | overlay, each normalized to unit max |

`sweep_output/` (regenerated, not committed): one folder per ratio
(`N1-200_N2-000/`, …) each holding `peaks_raw.csv`, `resample_histogram.csv` and
`reference_density.csv`, plus `sweep_index.csv` and a `plots/` folder with
per-ratio overlays, `montage.png` and `transition.png`.

---

## 6. The vendored transform

`Step_2_2D_to_3D_transformation.py`, `config.py` and `data_io.py` are byte-for-byte
copies of the files in `../2D to 3D Transformation/`. They are here so the
validation imports the same transform the paper used.

Two consequences:

- **They must be kept in sync.** Nothing enforces it. If the transform is ever
  changed, re-copy all three or this validation silently tests a stale version.
- **`Step_2_2D_to_3D_transformation.py` will not run standalone from this
  folder.** Its `main()` calls `data_io.check_integrity()`, which requires
  `Step_1_automaticXbinning.py` and a `settings/` file that live in the
  transformation folder, not here. That is expected and harmless: `check_integrity`
  is only called from `main()`, and these scripts import only the
  `invert_2d_to_3d` function, so importing works. Run the transform itself from
  `../2D to 3D Transformation/`.

---

## 7. Result as committed

100,000 recovered peaks per condition (10 pools × 10,000 subsamples), no failed
transforms:

| Condition | pool | in central band | in peripheral band | peripheral mode |
|---|---|---|---|---|
| mixture | 170 : 30 | 53.4% | 46.6% | 25.1 nm |
| all-peripheral (control) | 200 : 0 | 12.4% | 87.6% | 27.5 nm |
| all-central (control) | 0 : 200 | 100% | 0% | — |

The controls behave as required: an all-central truth never produces a
peripheral peak, and an all-peripheral truth puts the large majority of its
peaks in the peripheral band while the mixture splits roughly evenly.
