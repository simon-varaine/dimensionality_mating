################################################################################
## build_dataset.R
## AUTHORS ONLY. Rebuilds the shared analysis-ready dataset from the raw sources
## in data/ and writes it to data-derived/analysis_data.csv.
##
## Anyone reproducing the analyses does NOT need this script or the raw data:
## they use the committed data-derived/analysis_data.csv via run_all.R.
##
## Requirements: all raw files present in data/ (see README data manifest).
## Usage: run from the project root -> source("build_dataset.R")
################################################################################

source("R/00_setup.R")       # libraries, paths, helpers, constants
source("R/01_load_raw.R")    # load + clean the six raw sources
source("R/02_build_dataset.R")  # assemble + write data-derived/analysis_data.csv

cat("\nDone. Shared dataset written to:\n  ", path_analysis_data, "\n")
