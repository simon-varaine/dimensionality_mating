################################################################################
## run_all.R
## Reproduces every analysis, figure and table in the manuscript
## "Three-dimensional lifestyles are associated with reduced contest and
##  increased display in birds".
##
## This entry point needs only:
##   - data-derived/analysis_data.csv  (shared in this repository)
##   - data/AllBirdsHackett1.tre       (downloaded from birdtree.org; see README)
## It does NOT require the raw third-party datasets. (Authors rebuild the shared
## dataset from raw sources with build_dataset.R.)
##
## Usage: run from the project root (the folder containing data-derived/ and R/):
##   source("run_all.R")          # in an R session
##   Rscript run_all.R            # from the command line
##
## Outputs (figures + tables) are written to output/.
################################################################################

t0 <- Sys.time()

source("R/00_setup.R")         # libraries, paths, helpers, constants
source("R/03_load_dataset.R")  # read analysis_data.csv -> merged, avo_base
source("R/04_trees.R")         # master tree -> 100 sampled trees (object `trees`)
source("R/05_models.R")        # all PGLS models (res_* objects)
source("R/06_figures.R")       # Figure 1, Figure 2, Figure 3
source("R/07_tables.R")        # Table 1, A1, A2, A3, A4

cat(sprintf("\n=== Finished in %.1f min. Outputs in: %s ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), dir_out))
