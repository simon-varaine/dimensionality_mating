################################################################################
## 00_setup.R
## Arena dimensionality and the form of sexual selection in birds.
## Paths, libraries, helper functions and constants SHARED by every module.
## This file performs NO computation: it is sourced first (by build_dataset.R
## and by run_all.R); the numbered modules then rely on the objects defined here.
################################################################################

## ---------------------------------------------------------------- libraries ---
library(readr)
library(readxl)
library(data.table)
library(dplyr)
library(stringr)
library(janitor)
library(tidyr)
library(ape)
library(phylolm)
library(ggplot2)

## -------------------------------------------------------------------- paths ---
## Repository layout:
##   <root>/data          raw inputs from third parties  (git-ignored, NOT shared)
##   <root>/data-derived  analysis-ready dataset          (shared on GitHub)
##   <root>/output        generated figures + tables      (git-ignored)
## The root is detected via {here} if available, otherwise the working directory.
if (requireNamespace("here", quietly = TRUE)) {
  dir_root <- here::here()
} else {
  dir_root <- getwd()
  message("Package {here} not found: root = getwd() = ", dir_root,
          "\n  -> start R from the project root, or install {here}.")
}
dir_data    <- file.path(dir_root, "data")           # raw sources (not shared)
dir_derived <- file.path(dir_root, "data-derived")   # shared analysis dataset
dir_out     <- file.path(dir_root, "output")         # generated outputs
for (d in c(dir_derived, dir_out))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

## --- raw input files (original file names, as downloaded from each source) ---
## Only needed to REBUILD the analysis dataset (build_dataset.R). Not required
## to reproduce the analyses from the shared dataset (run_all.R). Drop each file
## into data/ under the name it ships with (see README data manifest).
path_lislevand <- file.path(dir_data, "avian_ssd_jan07.txt")                        # Lislevand et al. 2007
path_avonet    <- file.path(dir_data, "AVONET3_BirdTree.xlsx")                      # Tobias et al. 2022
avonet_sheet   <- "AVONET3_BirdTree"
path_blio      <- file.path(dir_data, "BLIOCPhyloMasterTax.csv")                    # BirdTree taxonomy
path_dale      <- file.path(dir_data, "plumage_scores.csv")                         # Dale et al. 2015
path_marcondes <- file.path(dir_data, "Mating_systems_master_datasheet_10nov2023.xlsx")  # Marcondes & Douvas 2024
marcondes_sheet<- "Species_data"
path_spur      <- file.path(dir_data, "species_spur_data.csv")                      # Menezes & Palaoro 2022

## --- shared analysis-ready dataset (committed to GitHub) ---------------------
## Written by build_dataset.R, read by run_all.R. One row per species, keyed by
## tip_label, holding every variable used in the analyses (see README).
path_analysis_data <- file.path(dir_derived, "analysis_data.csv")

## --- phylogeny (VertLife; large file, NOT committed -- see README) -----------
## 1000 Hackett all-species trees (9993 tips); we sample n_trees_use of them.
path_master    <- file.path(dir_data, "AllBirdsHackett1.tre")

## --- global parameters --------------------------------------------------------
n_trees_use   <- 100   # number of trees sampled from the 1000
SEP_THRESHOLD <- 8     # (quasi-)separation guard for rare binary responses

## --- dimensionality contrast convention --------------------------------------
## dim_bin is a factor with levels c("2D","3D"), reference = "2D".
## Hence the model term for the 3D effect is "dim_bin3D" everywhere.
DIM_LEVELS <- c("2D", "3D")
COEF_3D    <- "dim_bin3D"

################################################################################
## DATA HELPERS ----------------------------------------------------------------
################################################################################

## normalise a binomial (join key): lower-case, single spaces, no underscores
clean_binom <- function(x) {
  x %>% str_replace_all("_", " ") %>% str_squish() %>% str_to_lower()
}

## first column whose name matches a pattern (robust reading of AVONET)
get_col <- function(df, pattern) {
  hit <- names(df)[str_detect(names(df), regex(pattern, ignore_case = TRUE))]
  if (length(hit)) hit[1] else NA_character_
}

################################################################################
## PGLS ENGINE -----------------------------------------------------------------
## All fitters assume:
##   - a data.frame indexed by tip_label (rownames = tree tips);
##   - the global object `trees` (list of 100 trees) created by 04_trees.R.
## Two families:
##   .fit_binary   -> phyloglm logistic  (binary responses: harem, lek, ...)
##   .fit_gaussian -> phylolm lambda     (continuous responses: SSD, dichrom, ...)
################################################################################

## build a data.frame indexed by tip_label, filtered to non-NA response.
## `data` defaults to `merged` (created by 03_load_dataset.R).
.prep_dat <- function(data = merged, response, extra_filter = NULL) {
  d <- data %>%
    filter(!is.na(dim_bin), !is.na(.data[[response]]), !is.na(tip_label))
  if (!is.null(extra_filter)) d <- d %>% filter(!!extra_filter)
  d <- d %>% distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
  rownames(d) <- d$tip_label
  d
}

## --- binary PGLS (phyloglm) on one tree --------------------------------------
## Returns NULL if phyloglm errors OR fails to converge (convergence warning),
## so that coefficients from unstable fits are never reported.
.fit_binary <- function(tree, response, dat, formula = NULL) {
  keep <- intersect(tree$tip.label, rownames(dat))
  if (length(keep) < 10) return(NULL)
  f <- if (is.null(formula)) as.formula(paste(response, "~ dim_bin")) else formula
  fit <- tryCatch({
    tr <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
    d  <- dat[tr$tip.label, ]
    if (sum(d[[response]], na.rm = TRUE) < 5) return(NULL)
    withCallingHandlers(
      phylolm::phyloglm(f, data = d, phy = tr,
                        method = "logistic_MPLE", btol = 30),
      warning = function(w) {
        if (grepl("converge", conditionMessage(w), ignore.case = TRUE))
          invokeRestart("muffleWarning")   # swallow warning, handled via $convergence
      })
  }, error = function(e) NULL)
  if (!is.null(fit) && !is.null(fit$convergence) && fit$convergence != 0)
    return(NULL)
  fit
}

## --- continuous PGLS (phylolm lambda) on one tree ----------------------------
## `drop_na_response`: if TRUE, first drop species with NA on the response
## (useful when `dat` carries several z-scored responses at once).
.fit_gaussian <- function(tree, response, dat, formula = NULL,
                          drop_na_response = FALSE) {
  if (drop_na_response) dat <- dat[!is.na(dat[[response]]), , drop = FALSE]
  keep <- intersect(tree$tip.label, rownames(dat))
  if (length(keep) < 20) return(NULL)
  f <- if (is.null(formula)) as.formula(paste(response, "~ dim_bin")) else formula
  tryCatch({
    tr <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
    d  <- dat[tr$tip.label, ]
    phylolm::phylolm(f, data = d, phy = tr, model = "lambda")
  }, error = function(e) NULL)
}

## --- loop over the 100 trees + summarise one coefficient ---------------------
## Output: across-tree median + range, % predicted sign, % p<0.05 / p<0.10
## (per-tree Wald), and a Rubin-pooled p (within + between-tree uncertainty).
## Feeds Table A2 (robustness) of the manuscript.
run_pgls <- function(response, type = "gaussian",
                     coef_name = COEF_3D,
                     direction = c("negative", "positive"),
                     formula = NULL, dat = NULL, label = response) {
  direction <- match.arg(direction)
  if (is.null(dat)) dat <- .prep_dat(merged, response)
  fitter <- if (type == "binary") .fit_binary else .fit_gaussian

  fits <- lapply(trees, fitter, response = response, dat = dat, formula = formula)
  ok   <- Filter(Negate(is.null), fits)

  cat(sprintf("\n== %-22s | N=%d | trees OK=%d/%d | predicted dir.=%s ==\n",
              label, nrow(dat), length(ok), length(trees), direction))
  if (length(ok) == 0) { cat("  No tree converged.\n"); return(invisible(NULL)) }
  if (length(ok) < 80)
    cat(sprintf("  /!\\ only %d/%d trees converged: result NOT RELIABLE",
                length(ok), length(trees)),
        "(instability / quasi-separation); do not report.\n")

  coefs <- sapply(ok, function(f) coef(f)[coef_name])
  ses   <- sapply(ok, function(f) sqrt(diag(vcov(f)))[coef_name])
  sgn   <- if (direction == "negative") -1 else 1

  ## per-tree two-sided Wald p
  pvals     <- 2 * pnorm(-abs(coefs / ses))
  pct_sig   <- round(100 * mean(pvals < 0.05, na.rm = TRUE))
  pct_sig10 <- round(100 * mean(pvals < 0.10, na.rm = TRUE))

  ## Rubin (1987) pooling: W within-tree, B between-tree,
  ## T = W + (1 + 1/M) B -> single pooled z-test / p-value.
  M        <- length(ok)
  W        <- mean(ses^2, na.rm = TRUE)
  B        <- stats::var(coefs)
  T_var    <- W + (1 + 1/M) * B
  se_rubin <- sqrt(T_var)
  z_rubin  <- mean(coefs) / se_rubin
  p_rubin  <- 2 * pnorm(-abs(z_rubin))

  cat(sprintf("  median = %+.4f | 95%% range = [%+.4f, %+.4f] | %d%% in direction | %d%% p<0.05 | %d%% p<0.10 (Wald)\n",
              median(coefs), quantile(coefs, 0.025), quantile(coefs, 0.975),
              round(100 * mean(coefs * sgn > 0)), pct_sig, pct_sig10))
  cat(sprintf("  Rubin pooled: SE=%.4f | z=%.3f | p=%.4g\n",
              se_rubin, z_rubin, p_rubin))

  invisible(list(coefs = coefs, ses = ses, pvals = pvals,
                 pct_sig = pct_sig, pct_sig10 = pct_sig10,
                 se_rubin = se_rubin, z_rubin = z_rubin, p_rubin = p_rubin,
                 n = nrow(dat), direction = direction, label = label))
}

## --- generic control for ONE continuous z-scored covariate -------------------
## Generalises the body-size (Rensch) control to any covariate (size, HWI, ...).
## Returns the raw and covariate-adjusted 3D effects, and prints the covariate's
## own coefficient. Feeds Table A3 (Rensch, covar = body_size_log) and Table A4
## (display, covar = hwi).
run_pgls_ctrl <- function(response, covar, direction, label,
                          type = "gaussian", data = merged) {
  base <- data %>%
    filter(!is.na(dim_bin), !is.na(.data[[response]]),
           !is.na(tip_label), !is.na(.data[[covar]])) %>%
    distinct(tip_label, .keep_all = TRUE)
  covz <- paste0(covar, "_z")
  base[[covz]] <- as.numeric(scale(base[[covar]]))
  dat <- as.data.frame(base); rownames(dat) <- dat$tip_label

  res_raw <- run_pgls(response, type = type, direction = direction, dat = dat,
                      label = paste0(label, " [raw]"))
  res_ctl <- run_pgls(
    response, type = type, direction = direction, coef_name = COEF_3D,
    formula = as.formula(paste(response, "~ dim_bin +", covz)),
    dat = dat, label = paste0(label, " [+ ", covar, "]"))

  ## the covariate's own coefficient (expected sign: Rensch > 0; HWI < 0)
  if (!is.null(res_ctl)) {
    fitter <- if (type == "binary") .fit_binary else .fit_gaussian
    fits <- lapply(trees, fitter, response = response, dat = dat,
                   formula = as.formula(paste(response, "~ dim_bin +", covz)))
    okf  <- Filter(Negate(is.null), fits)
    if (length(okf)) {
      cc <- sapply(okf, function(f) coef(f)[covz])
      cat(sprintf("   -> coef %s: median = %+.4f [%+.4f, %+.4f]\n",
                  covar, median(cc), quantile(cc, .025), quantile(cc, .975)))
    }
  }
  invisible(list(raw = res_raw, ctl = res_ctl))
}

################################################################################
## RARE-BINARY HELPERS (spurs, Marcondes systems) ------------------------------
## 2D/3D contingency table, separation guard, exact Fisher test.
################################################################################

## 2x2 table (dim_bin x {0,1}) for a binary variable `var` in dataset `data`.
xtab_dim <- function(data, var) {
  d <- data %>%
    filter(dim_bin %in% DIM_LEVELS, !is.na(.data[[var]])) %>%
    distinct(key, .keep_all = TRUE) %>%
    as.data.frame()
  d$dim_bin <- droplevels(d$dim_bin)
  table(dim_bin = d$dim_bin, val = d[[var]])
}

## count of the rarest positive cell (2D vs 3D) in a 2x2 table.
min_pos_cell <- function(tab) {
  if (!("1" %in% colnames(tab))) return(0)
  min(tab[, "1"])
}

## exact Fisher test on the 2x2 table, with formatted printing.
fisher_dim <- function(tab, label = "") {
  ft <- fisher.test(tab)
  cat(sprintf("  [%s] Fisher exact: OR=%.3f, p=%.3g\n",
              label, unname(ft$estimate), ft$p.value))
  ft
}

cat("00_setup.R loaded: paths, helpers and constants ready.\n")
