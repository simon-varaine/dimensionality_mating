################################################################################
## 05_models.R
## All PGLS models. Every result object (res_*) is on the RAW response scale and
## feeds the tables (A2/A3/A4) and Figure 1A. The z-scored versions used for the
## standardised panels of Figures 1B/3B are refitted inside 06_figures.R.
##
## Requires 00_setup.R (engine + helpers), 03_load_dataset.R (merged, avo_base),
## 04_trees.R (trees). Runtime: a few to a few tens of minutes (100 trees each).
##
## Sections:
##   1. Mating system      (harem, lek | Marcondes RDP, lekking)
##   2. Size dimorphism    (5 SSD traits) + Rensch body-size control
##   3. Display agility    (main + excl. aerial + HWI control + interaction)
##   4. Dichromatism       (big-N main effect + interaction)
##   5. Bony spurs         (descriptive + Fisher; no PGLS -- separation)
##   6. Marcondes systems  (contingencies + aerial-lek family diagnostic)
##   7. Passeriformes-only robustness (for Figure 3A)
################################################################################

## helper: interaction dataset (response x resource) on the Lislevand backbone.
## resource is z-scored; note the interaction coefficient is invariant to this
## centring, so its value equals the raw-resource interaction.
.interaction_dat <- function(response) {
  d <- merged %>%
    filter(!is.na(dim_bin), !is.na(.data[[response]]),
           !is.na(resource), !is.na(tip_label)) %>%
    mutate(resource_c = as.numeric(scale(resource))) %>%
    distinct(tip_label, .keep_all = TRUE) %>%
    as.data.frame()
  rownames(d) <- d$tip_label
  d
}
## Interaction term names. NB: to reproduce the manuscript EXACTLY, the two
## interactions currently use DIFFERENT resource scalings -- a pre-existing
## inconsistency in the original code, flagged here for harmonisation in the
## revision phase (do NOT change it during the refactor):
##   - display      uses RAW resource     -> Table A4 / Fig 2 caption value -0.12
##   - dichromatism uses SCALED resource_c-> Fig 2 caption value            -1.12
## (The interaction coefficient is invariant to centring but not to scaling, so
##  the two are not on a comparable scale.)
COEF_INT_C   <- paste0(COEF_3D, ":resource_c")   # scaled  (dichromatism)
COEF_INT_RAW <- paste0(COEF_3D, ":resource")     # raw     (display)

################################################################################
## 1. MATING SYSTEM -----------------------------------------------------------
## Predictions: harem / resource-defense polygamy DOWN in 3D (contest);
##              lekking not predicted to fall (choice-based system).
################################################################################
cat("\n########## 1. MATING SYSTEM ##########\n")

res_harem <- run_pgls("harem", type = "binary", direction = "negative",
                      label = "Harem polygyny (Lislevand)")
res_lek   <- run_pgls("lek",   type = "binary", direction = "positive",
                      label = "Lek/promiscuity (Lislevand)")

## Independent replication (Marcondes & Douvas), on avo_base (decoupled from
## Lislevand), so N is not capped by Lislevand coverage.
res_marc_poly <- run_pgls("marc_poly", type = "binary", direction = "negative",
                          dat = .prep_dat(avo_base, "marc_poly"),
                          label = "Resource-defense polygamy (Marcondes)")
res_marc_lek  <- run_pgls("marc_lek",  type = "binary", direction = "positive",
                          dat = .prep_dat(avo_base, "marc_lek"),
                          label = "Lekking (Marcondes)")

## --- 1b. Monogamy-referenced contrasts (revision, §2.3) ----------------------
## The one-vs-rest contrasts above conflate the reference category (e.g. a
## negative marc_poly could mean less RDP, more monogamy, OR more lekking). Refit
## each polygamous system against MONOGAMY as a clean reference, dropping the
## third category, so RDP-vs-monogamy and lek-vs-monogamy are separable -- the
## dissociation that carries the paper.
cat("\n--- 1b. Marcondes contrasts vs monogamy (clean reference) ---\n")

## RDP vs monogamy (exclude leks)
dat_rdp_mono <- avo_base %>%
  filter(marc_system %in% c("M", "P"), !is.na(dim_bin), !is.na(tip_label)) %>%
  mutate(rdp = as.integer(marc_system == "P")) %>%
  distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
rownames(dat_rdp_mono) <- dat_rdp_mono$tip_label
res_marc_rdp_vs_mono <- run_pgls("rdp", type = "binary", direction = "negative",
                                 dat = dat_rdp_mono,
                                 label = "RDP vs monogamy (Marcondes)")

## lek vs monogamy (exclude RDP)
dat_lek_mono <- avo_base %>%
  filter(marc_system %in% c("M", "L"), !is.na(dim_bin), !is.na(tip_label)) %>%
  mutate(lek_m = as.integer(marc_system == "L")) %>%
  distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
rownames(dat_lek_mono) <- dat_lek_mono$tip_label
res_marc_lek_vs_mono <- run_pgls("lek_m", type = "binary", direction = "positive",
                                 dat = dat_lek_mono,
                                 label = "Lek vs monogamy (Marcondes)")

## lek vs RDP directly (drop monogamy): the two polygamous systems head-to-head.
## Most direct statement of the dissociation, but NOT independent of the two
## contrasts above (log-odds L-vs-P = L-vs-M minus P-vs-M). Expect a large
## positive median but wide across-tree range: the lek signal is concentrated in
## hummingbirds, which the phylogenetic model correctly discounts.
dat_lek_rdp <- avo_base %>%
  filter(marc_system %in% c("P", "L"), !is.na(dim_bin), !is.na(tip_label)) %>%
  mutate(lek_vs_rdp = as.integer(marc_system == "L")) %>%
  distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
rownames(dat_lek_rdp) <- dat_lek_rdp$tip_label
res_marc_lek_vs_rdp <- run_pgls("lek_vs_rdp", type = "binary", direction = "positive",
                                dat = dat_lek_rdp,
                                label = "Lek vs RDP (Marcondes)")

################################################################################
## 2. SIZE DIMORPHISM (SSD) + RENSCH CONTROL ----------------------------------
## All size traits predicted DOWN in 3D (contest/armament axis). Then the
## body-size control (Rensch's rule): refit each SSD with +body_size_log.
################################################################################
cat("\n########## 2. SIZE DIMORPHISM + RENSCH ##########\n")

res_ssd_mass   <- run_pgls("ssd_mass",   direction = "negative", label = "Mass dimorphism")
res_ssd_tarsus <- run_pgls("ssd_tarsus", direction = "negative", label = "Tarsus dimorphism")
res_ssd_wing   <- run_pgls("ssd_wing",   direction = "negative", label = "Wing dimorphism")
res_ssd_bill   <- run_pgls("ssd_bill",   direction = "negative", label = "Bill dimorphism")
res_ssd_tail   <- run_pgls("ssd_tail",   direction = "negative", label = "Tail dimorphism")

## Rensch: raw vs body-size-controlled 3D effect for each SSD trait (Table A3).
ssd_traits <- c("ssd_mass", "ssd_tarsus", "ssd_wing", "ssd_bill", "ssd_tail")
cat("\n--- Rensch body-size control (raw vs + body_size_log) ---\n")
res_rensch <- lapply(ssd_traits, function(v)
  run_pgls_ctrl(v, covar = "body_size_log", direction = "negative",
                label = v, type = "gaussian", data = merged))
names(res_rensch) <- ssd_traits

################################################################################
## 3. DISPLAY AGILITY ---------------------------------------------------------
## Main effect (predicted UP in 3D) + two robustness checks and the
## resource-sharing interaction (Table A4).
################################################################################
cat("\n########## 3. DISPLAY AGILITY ##########\n")

res_display <- run_pgls("display_num", direction = "positive",
                        label = "Display agility (main)")

## robustness 1: drop Aerial species (partial definitional overlap with 3D)
res_display_excl_aerial <- run_pgls(
  "display_num", direction = "positive",
  dat = .prep_dat(merged, "display_num",
                  extra_filter = rlang::quo(
                    primary_lifestyle %in% c("Insessorial", "Terrestrial"))),
  label = "Display agility (excl. aerial)")

## robustness 2: control for flight efficiency (hand-wing index)
res_display_hwi <- run_pgls_ctrl("display_num", covar = "hwi",
                                 direction = "positive",
                                 label = "Display agility ~ HWI",
                                 type = "gaussian", data = merged)

## interaction: display x between-sex resource sharing (RAW resource; see note)
res_display_interact <- run_pgls(
  "display_num", type = "gaussian", direction = "negative",
  coef_name = COEF_INT_RAW, formula = display_num ~ dim_bin * resource,
  dat = .interaction_dat("display_num"),
  label = "Display x resource (interaction)")

################################################################################
## 4. PLUMAGE DICHROMATISM ----------------------------------------------------
## Main effect at large N on avo_base (Dale, decoupled from Lislevand), plus the
## resource-sharing interaction (on the Lislevand backbone, which carries the
## resource variable). Fig 2 marginal plots are drawn in 06_figures.R.
################################################################################
cat("\n########## 4. DICHROMATISM ##########\n")

res_dichrom_bigN <- run_pgls(
  "dichromatism", type = "gaussian", direction = "positive",
  dat = .prep_dat(avo_base, "dichromatism"),
  label = "Plumage dichromatism (main, big-N)")

res_dichrom_interact <- run_pgls(
  "dichromatism", type = "gaussian", direction = "negative",
  coef_name = COEF_INT_C, formula = dichromatism ~ dim_bin * resource_c,
  dat = .interaction_dat("dichromatism"),
  label = "Dichromatism x resource (interaction)")

################################################################################
## 4b. PLUMAGE DECOMPOSITION (male / female / difference)   [revision, §2.5] ---
## Dichromatism = male - female nets out drivers shared by both sexes, so it is
## insensitive to MUTUAL ornamentation (both sexes elaborate). To see whether the
## weak dichromatism signal hides mutual ornamentation, fit male and female
## elaboration separately (big-N, avo_base), alongside the difference.
##   Internal check: beta(male) - beta(female) should approximate beta(dichrom).
## CAVEAT: absolute male/female levels REINTRODUCE shared ecological drivers
## (predation, light environment; cf. Dunn et al. 2015) that the difference
## removes. Treat these as exploratory; lean on the difference and, later, on the
## resource interaction for causal inference.
################################################################################
cat("\n########## 4b. PLUMAGE DECOMPOSITION (male / female) ##########\n")

res_male_plumage   <- run_pgls(
  "male_plumage", type = "gaussian", direction = "positive",
  dat = .prep_dat(avo_base, "male_plumage"),
  label = "Male plumage elaboration (big-N)")
res_female_plumage <- run_pgls(
  "female_plumage", type = "gaussian", direction = "positive",
  dat = .prep_dat(avo_base, "female_plumage"),
  label = "Female plumage elaboration (big-N)")

cat(sprintf("\nDecomposition check: beta(male) - beta(female) = %+.4f  vs  beta(dichrom) = %+.4f\n",
            median(res_male_plumage$coefs) - median(res_female_plumage$coefs),
            median(res_dichrom_bigN$coefs)))

################################################################################
## 5. BONY SPURS --------------------------------------------------------------
## Spurs are almost absent from 3D species -> (quasi-)complete separation, so a
## phylogenetic logistic regression is not identifiable. Report the descriptive
## contingency + exact Fisher test (as in the manuscript).
################################################################################
cat("\n########## 5. BONY SPURS ##########\n")

spur_tab <- xtab_dim(avo_base, "spur_hi")
cat("Spur presence by dimensionality:\n"); print(spur_tab)

if (min_pos_cell(spur_tab) < SEP_THRESHOLD) {
  cat(sprintf(">>> (quasi-)separation (rarest positive cell < %d): descriptive + Fisher only.\n",
              SEP_THRESHOLD))
  spur_fisher <- fisher_dim(spur_tab, "spur ~ 3D")
  res_spur    <- NULL
} else {
  res_spur    <- run_pgls("spur_hi", type = "binary", direction = "negative",
                          dat = .prep_dat(avo_base, "spur_hi"),
                          label = "Spur presence")
  spur_fisher <- NULL
}

################################################################################
## 6. MARCONDES SYSTEMS: contingencies + aerial-lek diagnostic ----------------
## The PGLS are in section 1; here, the raw contingencies and a check that the
## aerial lek signal is concentrated in hummingbirds (supports a Discussion
## statement), rather than spread across families.
################################################################################
cat("\n########## 6. MARCONDES: contingencies + diagnostic ##########\n")

for (v in c("marc_poly", "marc_lek")) {
  tb <- xtab_dim(avo_base, v)
  cat("\nContingency", v, "x dimensionality:\n"); print(tb)
  cat(sprintf("  rarest positive cell: %d (separation threshold %d)\n",
              min_pos_cell(tb), SEP_THRESHOLD))
}

aerial_lek_fam <- avo_base %>%
  filter(primary_lifestyle == "Aerial", marc_lek == 1) %>%
  distinct(key, .keep_all = TRUE) %>%
  count(family, order, sort = TRUE)
cat("\nAerial lekking species by family:\n")
print(as.data.frame(aerial_lek_fam), row.names = FALSE)
n_aerial_lek <- sum(aerial_lek_fam$n)
n_troch      <- sum(aerial_lek_fam$n[aerial_lek_fam$family == "Trochilidae"])
cat(sprintf("Total aerial lek: %d | Trochilidae: %d (%.0f%%) | distinct families: %d\n",
            n_aerial_lek, n_troch, 100 * n_troch / max(n_aerial_lek, 1),
            nrow(aerial_lek_fam)))

################################################################################
## 7. PASSERIFORMES-ONLY ROBUSTNESS (for Figure 3A) ---------------------------
## Panel A of Fig 3 needs harem (Lislevand) and RDP (Marcondes) refitted within
## Passeriformes. Lekking is not estimable there (near-complete separation), so
## it is omitted, as in the manuscript. The z-scored morphology/behaviour panel
## (Fig 3B) is refitted inside 06_figures.R.
################################################################################
cat("\n########## 7. PASSERIFORMES ROBUSTNESS ##########\n")

res_harem_p <- run_pgls(
  "harem", type = "binary", direction = "negative",
  dat = .prep_dat(merged, "harem", extra_filter = rlang::quo(order == "Passeriformes")),
  label = "Harem polygyny [Passeriformes]")

## RDP within Passeriformes, with the separation guard
d_rdp_pass <- avo_base %>%
  filter(order == "Passeriformes", !is.na(dim_bin), !is.na(marc_poly),
         !is.na(tip_label)) %>%
  distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
rownames(d_rdp_pass) <- d_rdp_pass$tip_label
tb_rdp_pass <- xtab_dim(d_rdp_pass, "marc_poly")
cat(sprintf("[RDP | Passeriformes] N=%d | rarest positive cell=%d\n",
            nrow(d_rdp_pass), min_pos_cell(tb_rdp_pass)))
if (min_pos_cell(tb_rdp_pass) < SEP_THRESHOLD) {
  cat("  -> (quasi-)separation: not reported.\n")
  res_marc_poly_p <- NULL
} else {
  res_marc_poly_p <- run_pgls("marc_poly", type = "binary", direction = "negative",
                              dat = d_rdp_pass,
                              label = "Resource-defense polygamy [Passeriformes]")
}

cat("\n05_models.R done. Result objects (res_*) are in memory.\n")
