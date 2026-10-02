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
## Outputs:
##   output/figures/  Figures 1, S1, S2  (.pdf for the paper, .png preview,
##                    .tex = ready-to-paste figure environment with caption)
##   output/tables/   Tables 1, S1, S2, S3  (.tex = ready-to-paste table
##                    environment with caption; .csv = same numbers, full precision)
##   output/model_results.rds  cache of every fitted model summary
##
## REFIT = TRUE (default) refits every model across the 100 trees (~30-40 min).
## REFIT = FALSE reuses the cache and only redraws figures and tables (seconds):
##   REFIT <- FALSE; source("run_all.R")
################################################################################

if (!exists("REFIT")) REFIT <- TRUE

t0 <- Sys.time()

source("R/00_setup.R")           # libraries, paths, helpers, constants
source("R/03_load_dataset.R")    # read analysis_data.csv -> merged, avo_base

if (REFIT || !file.exists(path_model_cache)) {
  source("R/04_trees.R")         # master tree -> 100 sampled trees (object `trees`)
  source("R/05_models.R")        # all models -> output/model_results.rds
} else {
  model_results <- readRDS(path_model_cache)
  cat("REFIT = FALSE: using cached models from", format(model_results$date), "\n")
}

source("R/06_figures.R")         # Figures 1, S1, S2
source("R/07_tables.R")          # Tables 1, S1, S2, S3

cat(sprintf("\n=== Finished in %.1f min. Outputs in: %s ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), dir_out))
