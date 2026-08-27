# Supplementary Figure S5B/S5C — Nup153C event-duration distribution (2 ms/frame)

Analysis code for the unfiltered event-duration distribution of Nup153C acquired at
2 ms/frame, and the geometric mixture-model fit used to argue for kinetically
distinct populations. 

## Pipeline

```
Data_Nup153/*.xlsx          raw localization tables (not in this repository)
        │
        │  01_extract_event_durations.py
        ▼
data/combined_group_lengths.csv  one row per event, duration in FRAMES
        │
        ├─ 02_fit_duration_mixture.py   → results/model_comparison.csv
```
