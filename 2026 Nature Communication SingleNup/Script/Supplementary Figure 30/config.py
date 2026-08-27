# -*- coding: utf-8 -*-
"""
Shared configuration for the 2D -> 3D transformation pipeline.
==============================================================

Every step (Step 1 automatic X binning, Step 2 transformation, ...) imports
its settings from here, so the folder names and analysis parameters are defined
in exactly one place and cannot drift between scripts.

Edit the values below to change behaviour. Folder names are resolved relative to
each script's own location, so you only set the *name* here, not a full path.
"""

# ---------------------------------------------------------------------------
# Folder layout (names only; resolved next to each script at run time)
# ---------------------------------------------------------------------------
DATA_DIR = "data"            # subfolder holding the input data files
SETTINGS_DIR = "settings"    # subfolder for the settings CSV (input + output)
SETTINGS_FILE = "Xbinning_settings.csv"

# ---------------------------------------------------------------------------
# Automatic X-binning parameters (Step 1)
# ---------------------------------------------------------------------------
# A bin is accepted once it contains at least this many localizations.
MIN_COUNT = 200
# Largest width (in nm) a single bin may grow to while searching for MIN_COUNT.
# If a bin cannot reach MIN_COUNT within this width, that side stops.
MAX_BIN_WIDTH = 100
# Safety cap on how many bins are added to each side of the center.
MAX_BINS_PER_SIDE = 20

# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------
# Accepted header spellings for the X and Y columns. Comparison is done after
# lower-casing and removing all spaces, so "X", "X (nm)" and "x(nm)" all match.
# Add a normalized spelling here if your data uses another convention.
X_HEADERS = {"x", "x(nm)"}
Y_HEADERS = {"y", "y(nm)"}

# File extensions treated as data files.
EXCEL_EXTS = {".xlsx", ".xls"}
CSV_EXTS = {".csv"}

# ---------------------------------------------------------------------------
# 2D -> 3D transformation parameters (Step 2)
# ---------------------------------------------------------------------------
# Parent folder holding all transformation results. Each run gets its own
# timestamped subfolder underneath it; identical reruns are detected and skipped.
RESULTS_DIR = "transformation_results"

# Soft upper limit (nm) of the radial (R) dimension. Rounded up to a whole
# number of bins, so the effective limit depends on the bin size.
R_upperlimit = 150

# Range of bin sizes (nm) to analyze. Set bin_min == bin_max for a single size.
bin_min = 4
bin_max = 12
# Step between successive bin sizes in the bin_min..bin_max sweep.
# e.g. bin_min=4, bin_max=12, bin_delta=2  ->  bin sizes 4, 6, 8, 10, 12.
# May be fractional, e.g. bin_delta = 0.5  ->  4, 4.5, 5, ... bin_min, bin_max,
# and bin_delta may all be floats.
bin_delta = 1
# Safety floor: any bin_delta smaller than this is treated as this value (with a
# warning). It also sets the resolution to which bin sizes are rounded for the
# output column names, so there is no floating-point drift in the headers.
BIN_DELTA_FLOOR = 0.001

# Mann-Whitney U p-value threshold for the Y-symmetry check. A region whose Y
# data is less symmetric than this (smaller p) raises a warning.
MWU_threshold = 0.001

# Background-noise handling for localizations beyond R_upperlimit.
background_filter = True
background_filter_range = 180

# Output column ordering:
#   "subregion" - bins of one subregion grouped together (ideal for fitting)
#   "bin"       - same bin size across subregions grouped together (ideal for comparison)
save_format = "subregion"
