# Supplementary Figure 23 — Nuclear-envelope line scan and alignment-error tilt distribution

Traces the nuclear envelope in line-scan images, fits its local curvature, and
converts the imperfect alignment of the imaged pore to the laser crosshair into
a realistic distribution of membrane tilt angles. That distribution is the
experimental input to the tilt Monte Carlo simulation.

---

## 1. Folder layout

```
Supplementary Figure 23/
    linescan_membrane.py   <- the analysis (edit the SETTINGS block)
    README.md              <- this file

    data/                  <- input text images, one CSV per cell
    output/                <- created by the script
```

Both folders are resolved next to the script, so the folder can be moved as a
unit.

---

## 2. Requirements

- Python 3.9+
- `numpy`, `matplotlib`

```
pip install numpy matplotlib
```

The script sets the `Agg` backend so it saves figures without opening windows.
Comment out that line to see plots interactively in Spyder.

---

## 3. Input format

One CSV per cell in `data/`, a text image exported from the line scan: columns
are image X, rows are image Y, values are intensity. The NE membrane must run
roughly **vertically** through the crop. Comma-delimited or whitespace-delimited
both work. The committed set is 8 cells, `AVG_Cell1-1.csv` … `AVG_Cell8-1.csv`.

---

## 4. Method

**Per cell:**

1. **Sub-pixel membrane peak per row.** Within a tracking window that follows the
   membrane from row to row, the local maximum is refined by fitting a parabola
   to log(intensity) — the vertex of that parabola is the Gaussian centre.
2. **Reject weak rows** below `MIN_SNR` or below `MIN_PROM_FRAC` of the maximum
   prominence. SNR uses a robust (MAD-based) noise estimate.
3. **Robust parabola fit** of x = a·y² + b·y + c with iterative σ-clipping,
   weighted by peak prominence. Note the fit is image-X against image-Y, i.e.
   `polyfit(y, x)`, so no manual coordinate flip is needed.
4. **Tangent angle vs vertical**, θ(y) = arctan(2a·y + b). The apex — the
   most-vertical point, dx/dy = 0 — is at y = −b/2a, and the radius of curvature
   there is 1/|2a|.

**Alignment-error tilt distribution (pooled across cells):**

The apex is what gets aligned to the laser crosshair, but alignment is imperfect
(~500 nm). The pore actually imaged therefore sits some distance δ off the apex,
where the membrane tilt is θ = arctan(2a·δ). Each cell has its own curvature a,
so sampling δ per cell and pooling gives the experimental tilt distribution.

- The tilt **angle** is scale-invariant (X and Y share a pixel size, and θ is a
  ratio), so px → nm conversion does not change it. Pixel size enters **only**
  through the offset: `delta_px = ALIGN_OFFSET_NM / PIXEL_SIZE_NM`.
- `OFFSET_MODEL = "uniform"` draws δ uniformly on ±offset — conservative, it puts
  the whole alignment error along the membrane arc. `"disk"` draws a 2D offset
  and projects it onto the membrane, giving a narrower distribution.

---

## 5. Settings

All in the `SETTINGS` block; there is no command line.

| Setting | Default | Meaning |
|---|---|---|
| `DATA_DIR`, `OUT_DIR` | `data`, `output` | input and output folders |
| `SIDE` | `"left"` | which arc the membrane is on (`left`/`right`/`auto`) |
| `TRACK_WIN` | 4 px | peak-search half-window per row |
| `FIT_PTS` | 2 | points either side of the local max used for the sub-pixel fit |
| `MIN_SNR` | 2.0 | reject rows below this peak SNR |
| `MIN_PROM_FRAC` | 0.15 | reject rows below this fraction of max prominence |
| `POLY_DEGREE` | 2 | 2 = parabola |
| `CLIP_SIGMA` | 2.5 | σ-clip threshold for outlier rejection |
| `EDGE_TRIM` | 0 | rows dropped at top and bottom before fitting |
| `PIXEL_SIZE_NM` | 160.0 | Andor EMCCD = 160, Cascade = 240 |
| `ALIGN_OFFSET_NM` | 500.0 | alignment uncertainty; ±3.1 px at 160 nm/px |
| `ALIGN_OFFSET_PX` | None | set a number to override the offset directly in px |
| `OFFSET_MODEL` | `"uniform"` | `uniform` (conservative) or `disk` |
| `N_SAMPLES` | 2000 | Monte Carlo samples per cell |
| `SEED` | 0 | fixed, for reproducibility |

Plot style (fonts, colours, `PLOT_UNITS`, `DPI`) is in the same block.

---

## 6. Running it

```
python linescan_membrane.py
```

Outputs land in `output/`:

| Output | Contents |
|---|---|
| `<stem>_membrane.csv` | per row: y, x_peak, prominence, used, x_fit, tangent, plus nm columns |
| `<stem>_overlay.png` | image with the fitted membrane curve |
| `<stem>_angle.png` | tangent angle vs y |
| `cells_summary.csv` | per cell: a, apex y, apex radius, fit RMSE, max tilt at the offset |
| `pooled_tilt_distribution.csv` | every sampled tilt angle across all cells |
| `pooled_tilt_hist.png` | histogram of the pooled distribution |

A cell that fails is reported and skipped; the batch continues.

---

## 7. Result as committed

8 cells. Apex radius of curvature 4.0–13.5 µm, fit RMSE 18–80 nm, and maximum
tilt at the ±500 nm offset 2.1–7.0° per cell.

---
