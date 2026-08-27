# 2D → 3D Transformation Pipeline

Reconstructs the 3D radial density of rotationally symmetric biological
structures (e.g. the nuclear pore complex) from 2D single-molecule localization
data, using an area-matrix inversion. The workflow is two steps: automatic X
binning, then the 2D→3D transformation.

---

## 1. Folder layout

Keep all four scripts together in one folder. The subfolders are created and
used automatically.

```
<main folder>/
    config.py                          <- all settings (edit this)
    data_io.py                         <- shared data loader + integrity check
    Step_1_automaticXbinning.py        <- Step 1
    Step_2_2D_to_3D_transformation.py  <- Step 2
    README.md                          <- this file

    data/                  <- YOU create this; put your data files here
    settings/              <- Step 1 writes the binning settings here
    transformation_results/<- Step 2 writes results here
```

> **Do not delete or rename any of the four `.py` files, and do not run Step 2
> before Step 1.** Each script checks for the others at startup and will stop
> with a clear message if something is missing or out of order.

---

## 2. Requirements

- Python 3.9+
- `numpy`, `pandas`, `scipy`
- `openpyxl` — for the Excel summary output (and for reading `.xlsx` data files).
  If it is missing, Step 2 falls back to writing `summary.csv` instead.

Install everything with:

```
pip install numpy pandas scipy openpyxl
```

---

## 3. How to run

### Step 0 — Prepare your data

Create a `data/` folder next to the scripts and put your localization files in
it. Accepted formats:

- **Excel** (`.xlsx`, `.xls`) or **CSV** (`.csv`).
- **With a header**: the file must have exactly one X column and one Y column.
  Accepted header spellings (any case, spaces optional): `X`, `X(nm)`, `X (nm)`
  and `Y`, `Y(nm)`, `Y (nm)`. Extra columns are ignored.
- **Without a header**: the file must have exactly two columns, read as X and Y.
- Blank rows (both X and Y empty) are allowed as separators and are dropped.
- A row with a value in only one of X/Y is treated as an error: the file is
  reported (with the offending line numbers) and skipped.

### Step 1 — Automatic X binning

```
python Step_1_automaticXbinning.py
```

For each data file this finds, along the X axis, a set of bins each containing
at least `MIN_COUNT` localizations, starting from an X center and growing
outward. It writes `settings/Xbinning_settings.csv`.

The **X center** is chosen as follows: if `settings/Xbinning_settings.csv`
already exists and lists a numeric `X Center` for that file, that value is
reused (so your manual edits survive a re-run); otherwise it defaults to the
median of X. Centers are rounded to integers so all bin edges are integers.

**You can edit `settings/Xbinning_settings.csv` by hand** (e.g. adjust an
`X Center`, a range, or the cut points) before running Step 2.

Settings file columns:

| Column | Meaning |
|---|---|
| `File` | data file name (in `data/`) |
| `X Center` | binning anchor along X |
| `Range_min`, `Range_max` | overall X range that was binned |
| `Localization Count` | counts per X bin |
| `Cutpoint` | the X bin edges (these become the subregion boundaries in Step 2) |

### Step 2 — 2D → 3D transformation

```
python Step_2_2D_to_3D_transformation.py
```

For each file and each X subregion (defined by the `Cutpoint` edges), this
reconstructs the radial density for every bin size in the configured sweep and
writes the results.

Output for each run:

```
transformation_results/
    result_<timestamp>/
        raw/            <- density divided by the X-range width (negatives kept)
        zero_negative/  <- negatives clipped to 0, then divided by max
        min_max/        <- negatives clipped to 0, then min-max scaled to [0, 1]
        summary.xlsx    <- chi-square + diagnostics, one sheet per data file
        run_meta.json   <- archived settings + data fingerprint
        warnings.txt    <- any warnings from the run
```

Each `raw/`, `zero_negative/`, `min_max/` folder contains one CSV per data file.

**Duplicate detection.** Every run is fingerprinted from the settings and the
contents of the data files. If you re-run with identical settings and identical
data, Step 2 detects the matching previous result, prints where it is, and does
nothing — so every result folder is guaranteed to differ in some way. Change a
setting or the data to produce a new result.

---

## 4. The `summary.xlsx` (bin-size diagnostics)

For every X region and bin size, Step 2 records a chi-square from a
back-calculation check (clip negative densities to zero, forward-project through
the area matrix to reconstruct the 2D counts, compare bin-by-bin to the original
counts). The workbook has **one sheet per data file** for easy navigation; if
`openpyxl` is unavailable it falls back to a single `summary.csv`. Columns:

| Column | Meaning |
|---|---|
| `File`, `X_region`, `Region_width`, `Bin_size` | which transformation |
| `N_localizations` | points in that region (the region "weight") |
| `N_bins_compared` | bins with data used in the comparison |
| `Chi_square` | reconstruction chi-square statistic |
| `Chi_square_dof`, `Chi_square_pvalue` | df = (bins with data) − 1, and its p-value |
| `Poisson_deviance_experimental` | experimental alternative statistic |

Intended use: for each region, choose the smallest bin size whose chi-square
falls into an acceptable range — small enough to keep resolution, large enough
to avoid the negative-density artifacts of under-sampling.

> **Caveats.** The reconstruction is computed in-sample on a deterministic
> transform, so this chi-square is a *relative* diagnostic, not a rigorous
> hypothesis test; read it comparatively across bin sizes, not as an absolute
> p-value. The Poisson deviance column is experimental and not yet validated.

---

## 5. Settings (`config.py`)

All tunable values live in `config.py`, shared by both steps so they never
disagree.

**Folders:** `DATA_DIR`, `SETTINGS_DIR`, `SETTINGS_FILE`, `RESULTS_DIR`.

**Step 1 (X binning):**
- `MIN_COUNT` — minimum localizations per X bin (default 200).
- `MAX_BIN_WIDTH` — largest width an X bin may grow to (nm).
- `MAX_BINS_PER_SIDE` — cap on bins added to each side of the center.

**Input parsing:** `X_HEADERS`, `Y_HEADERS` (accepted header spellings),
`EXCEL_EXTS`, `CSV_EXTS`.

**Step 2 (transformation):**
- `R_upperlimit` — soft upper limit of the radial dimension (nm).
- `bin_min`, `bin_max` — range of R bin sizes to analyze.
- `bin_delta` — step between bin sizes; may be fractional (e.g. `0.5` → 4, 4.5,
  5, …). Values below `BIN_DELTA_FLOOR` (0.001) are clamped to it, with a warning.
- `BIN_DELTA_FLOOR` — smallest allowed `bin_delta` and the rounding resolution
  for bin sizes in output column names (default 0.001).
- `MWU_threshold` — Mann-Whitney U threshold for the Y-symmetry warning.
- `background_filter`, `background_filter_range` — background-noise handling.
- `save_format` — `"subregion"` (group bins of a region together; good for
  fitting) or `"bin"` (group the same bin size across regions; good for
  comparison).

---

## 6. Troubleshooting

- **"Integrity check failed"** — a `.py` file or the `data/` folder is missing.
  Restore the missing file(s); all four scripts must sit in the same folder.
- **"Step 2 needs the settings file ... run Step 1 first"** — you ran Step 2
  before Step 1. Run `Step_1_automaticXbinning.py` first.
- **"row(s) have a value in only one of X/Y"** — a data file has a half-empty
  row at the listed line number(s). Fix the data and re-run.
- **"need exactly one X column and one Y column"** — the header doesn't match an
  accepted spelling, or there are duplicate/missing X or Y columns.
- **"The same settings and the same data files have already been transformed"** —
  not an error: an identical run already exists at the path shown. Change a
  setting or the data to produce a new result.
