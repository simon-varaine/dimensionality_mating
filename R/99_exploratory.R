################################################################################
## 99_exploratory.R
## Exploratory analyses that are NOT part of the manuscript, kept for the record
## and NOT sourced by run_all.R. Moved here from 05_models.R / 06_figures.R when
## the pipeline was restricted to the paper's analyses (Block A main effects,
## controls, continuous-care interaction). Results are summarised in
## "Echanges avec Jan/results_synthesis_for_Jan.md".
##
##   1. Resource-sharing interactions (display, dichromatism)   [retired: construct]
##   2. Care moderators: care_female (F vs P), dev_pc1            [drawer]
##   3. Intensity moderation, dim_bin x Barber intensity (B1)     [dropped]
##   4. Flightlessness (descriptive)                              [drawer]
##   5. HWI as an alternative dimensionality proxy                [negative result]
##   6. Old Figure 2 (resource-sharing) and Figure 4 (dev. mode)
##
## Usage (from the project root, after the main pipeline has created `trees`):
##   source("R/00_setup.R"); source("R/03_load_dataset.R"); source("R/04_trees.R")
##   source("R/99_exploratory.R")
################################################################################

################################################################################
## 1. RESOURCE-SHARING INTERACTIONS ------------------------------------------
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

## interaction: display x between-sex resource sharing (RAW resource; see note)
res_display_interact <- run_pgls(
  "display_num", type = "gaussian", direction = "negative",
  coef_name = COEF_INT_RAW, formula = display_num ~ dim_bin * resource,
  dat = .interaction_dat("display_num"),
  label = "Display x resource (interaction)")

res_dichrom_interact <- run_pgls(
  "dichromatism", type = "gaussian", direction = "negative",
  coef_name = COEF_INT_C, formula = dichromatism ~ dim_bin * resource_c,
  dat = .interaction_dat("dichromatism"),
  label = "Dichromatism x resource (interaction)")

################################################################################
## 2-3. CARE MODERATORS (binary / dev. mode) + INTENSITY MODERATION (B1) -----
################################################################################
## 8. CARE MODERATORS: interaction dim_bin x care on the choice traits  [rev.] -
## Prediction: the 3D -> male-ornament effect concentrates where care is LOW /
## female-borne (female choice) and attenuates where care is high / shared
## (mutual choice). Two moderators, from near-circular to distal:
##   care_female : F vs P            (PROXIMAL, near-circular; upper bound)     +
##   dev_pc1     : hatchling PCA     (DISTAL; HIGH = precocial = low care need) +
## (last column = predicted sign of the dim_bin3D:moderator interaction.)
## (Care duration from BirdBase set aside for now; dev_pc1 largely captures it.)
## Responses: male_plumage, dichromatism (avo_base); display_num (merged).
################################################################################
cat("\n########## 8. CARE MODERATORS (interactions) ##########\n")

## build the interaction dataset for one response x one moderator.
## The moderator is always renamed `mv`, so the interaction term is dim_bin3D:mv.
.interaction_care <- function(response, base, mod) {
  d <- base %>% filter(!is.na(dim_bin), !is.na(.data[[response]]), !is.na(tip_label))
  if (mod == "care_female") {
    d   <- d %>% filter(care_mode %in% c("F", "P")) %>%
      mutate(mv = as.integer(care_mode == "F"))
    dir <- "positive"
  } else if (mod == "dev_pc1") {
    d   <- d %>% filter(!is.na(dev_pc1)) %>% mutate(mv = as.numeric(scale(dev_pc1)))
    dir <- "positive"          # high dev_pc1 = precocial = low care
  }
  d <- d %>% distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
  rownames(d) <- d$tip_label
  list(dat = d, dir = dir)
}

care_grid <- expand.grid(
  response = c("male_plumage", "dichromatism", "display_num"),
  mod      = c("care_female", "dev_pc1"),
  stringsAsFactors = FALSE)

res_care_int <- list()
for (i in seq_len(nrow(care_grid))) {
  rsp  <- care_grid$response[i]; md <- care_grid$mod[i]
  base <- if (rsp == "display_num") merged else avo_base
  s    <- .interaction_care(rsp, base, md)
  lab  <- paste0(rsp, " x ", md)
  res_care_int[[lab]] <- run_pgls(
    rsp, type = "gaussian", direction = s$dir,
    coef_name = paste0(COEF_3D, ":mv"),
    formula   = as.formula(paste(rsp, "~ dim_bin * mv")),
    dat = s$dat, label = lab)
}

################################################################################
## 9. INTENSITY MODERATION (non-circular channel-switch test)          [rev.] --
## Does a given AMOUNT of sexual selection buy CONTEST traits in 2D but CHOICE
## traits in 3D? Fit  trait ~ dim_bin * intensity(Barber).
##   - intensity is measured independently of the morphology/plumage/display
##     traits, and is ~orthogonal to dimensionality (main effect -0.085), so the
##     interaction is NOT circular (unlike the main effect of mating system);
##   - applied to FORM traits only -- NOT to mating-system outcomes (circular);
##   - sex-role-reversed species EXCLUDED so intensity reflects male-biased
##     polygyny, not female-female competition (which would flip the prediction).
## Predicted dim_bin3D:intensity sign:
##   contest traits (SSD)               -> negative (intensity buys dimorphism in 2D)
##   choice traits (dichrom/male/display)-> positive (intensity buys ornament in 3D)
################################################################################
cat("\n########## 9. INTENSITY MODERATION (channel switch) ##########\n")

int_grid <- data.frame(
  response = c("ssd_mass", "ssd_tarsus", "ssd_wing", "ssd_bill", "ssd_tail",
               "dichromatism", "male_plumage", "female_plumage", "display_num"),
  base_nm  = c(rep("merged", 5), "avo_base", "avo_base", "avo_base", "merged"),
  dir      = c(rep("negative", 5), "positive", "positive", "positive", "positive"),
  stringsAsFactors = FALSE)

res_int_mod <- list()
for (i in seq_len(nrow(int_grid))) {
  rsp  <- int_grid$response[i]
  base <- get(int_grid$base_nm[i])
  d <- base %>%
    filter(!is.na(dim_bin), !is.na(.data[[rsp]]), !is.na(barber_ss), !is.na(tip_label),
           is.na(barber_srr) | barber_srr == 0) %>%   # exclude sex-role-reversed:
    mutate(ss_z = as.numeric(scale(barber_ss))) %>%    # intensity = male-biased polygyny
    distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
  rownames(d) <- d$tip_label
  res_int_mod[[rsp]] <- run_pgls(
    rsp, type = "gaussian", direction = int_grid$dir[i],
    coef_name = paste0(COEF_3D, ":ss_z"),
    formula   = as.formula(paste(rsp, "~ dim_bin * ss_z")),
    dat = d, label = paste0(rsp, " x intensity"))
}

################################################################################
## 4. FLIGHTLESSNESS (descriptive only) --------------------------------------
################################################################################
## (c) flightlessness -- descriptive only (rare, clade-concentrated)
fless_tab <- xtab_dim(avo_base, "flightless")
cat("Flightlessness by dimensionality:\n"); print(fless_tab)
if (min_pos_cell(fless_tab) < SEP_THRESHOLD) {
  cat(sprintf(">>> (quasi-)separation (rarest positive cell < %d): descriptive + Fisher only.\n",
              SEP_THRESHOLD))
  fless_fisher <- fisher_dim(fless_tab, "flightless ~ 3D")
} else {
  fless_fisher <- fisher_dim(fless_tab, "flightless ~ 3D")
}

################################################################################
## 5. HWI AS AN ALTERNATIVE DIMENSIONALITY PROXY (set aside) -----------------
################################################################################
## 11. HWI ROBUSTNESS: dimensionality proxied by hand-wing index       [rev.] --
## HWI (continuous flight efficiency ~ escape capacity) as an ALTERNATIVE,
## independent operationalization of dimensionality (HIGH HWI = more 3D-like).
## Triangulation with the foraging 2D/3D contrast. Displays EXCLUDED (HWI and
## display agility are near-tautological). Caveat: HWI also indexes dispersal/
## migration, so read as robustness, not a clean causal test. Covers all species
## with HWI (incl. aquatic/generalist), not just the 2D/3D categories.
## Predicted sign of the hwi coefficient mirrors the 3D effect:
##   contest (SSD, spurs, strong-poly/RDP, SRR, polyandry) -> negative
##   choice  (dichrom, male plumage, uv_cd, lek)           -> positive
################################################################################
## SET ASIDE (not run): HWI fails as a dimensionality proxy -- opposite sign on
## intensity (+0.20, p=1e-24, vs -0.085 for foraging), confounded by dispersal/
## migration. Kept for reference; wrapped in if(FALSE). Re-enable to reproduce
## the negative result that justifies choosing the foraging-lifestyle proxy.
if (FALSE) {
cat("\n########## 11. HWI ROBUSTNESS (alternative dimensionality proxy) ##########\n")

.hwi_dat <- function(response) {
  d <- avo_base %>%
    filter(!is.na(.data[[response]]), !is.na(hwi), !is.na(tip_label)) %>%
    mutate(hwi_z = as.numeric(scale(hwi))) %>%
    distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
  rownames(d) <- d$tip_label
  d
}

hwi_grid <- data.frame(
  response = c("ssd_mass", "ssd_tarsus", "ssd_wing", "ssd_bill", "ssd_tail",
               "spur_hi",
               "dichromatism", "male_plumage", "female_plumage", "uv_cd",
               "barber_poly", "marc_poly", "barber_lek", "marc_lek",
               "barber_ss", "allopreen", "barber_srr", "polyandry"),
  type = c("gaussian", "gaussian", "gaussian", "gaussian", "gaussian",
           "binary",
           "gaussian", "gaussian", "gaussian", "gaussian",
           "binary", "binary", "binary", "binary",
           "gaussian", "binary", "binary", "binary"),
  dir  = c("negative", "negative", "negative", "negative", "negative",
           "negative",
           "positive", "positive", "positive", "positive",
           "negative", "negative", "positive", "positive",
           "negative", "positive", "negative", "negative"),
  stringsAsFactors = FALSE)

res_hwi <- list()
for (i in seq_len(nrow(hwi_grid))) {
  rsp <- hwi_grid$response[i]
  res_hwi[[rsp]] <- run_pgls(
    rsp, type = hwi_grid$type[i], direction = hwi_grid$dir[i],
    coef_name = "hwi_z", formula = as.formula(paste(rsp, "~ hwi_z")),
    dat = .hwi_dat(rsp), label = paste0(rsp, " ~ HWI"))
}
}  # end if(FALSE): HWI set aside (fails as a dimensionality proxy; see synthesis)

################################################################################
## 6. OLD FIGURES 2 AND 4 -----------------------------------------------------
################################################################################
dir_out_explo <- file.path(dir_out, "exploratory")
if (!dir.exists(dir_out_explo)) dir.create(dir_out_explo, recursive = TRUE)

################################################################################
## FIGURE 2: resource-sharing interaction (marginal predictions, tree 1) ------
################################################################################
cat("\n=== Figure 2 ===\n")

## NB: this plot uses SCALED resource_c for both traits (so the two panels are
## drawn on a common resource axis). The manuscript's baked-in beta values were
## on inconsistent scales (see 05_models.R note); the caption here is qualitative
## to avoid re-embedding those numbers -- to be finalised in the revision.
marginal_grid <- function(response, trait_label) {
  d <- merged %>%
    filter(!is.na(dim_bin), !is.na(.data[[response]]),
           !is.na(resource), !is.na(tip_label)) %>%
    mutate(resource_c = as.numeric(scale(resource))) %>%
    distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
  rownames(d) <- d$tip_label
  tr <- ape::drop.tip(trees[[1]], setdiff(trees[[1]]$tip.label, rownames(d)))
  d1 <- d[tr$tip.label, ]
  fit <- phylolm::phylolm(as.formula(paste(response, "~ dim_bin * resource_c")),
                          data = d1, phy = tr, model = "lambda")
  co  <- coef(fit)
  res_levels <- quantile(d1$resource_c, c(0.15, 0.5, 0.85))
  g <- expand.grid(dim_bin = factor(DIM_LEVELS, levels = DIM_LEVELS),
                   resource_c = res_levels)
  g$dim3d <- as.integer(g$dim_bin == "3D")
  g$pred  <- co["(Intercept)"] + co[COEF_3D] * g$dim3d +
             co["resource_c"] * g$resource_c +
             co[COEF_INT_C] * g$dim3d * g$resource_c
  g$res_lab <- factor(rep(c("Low resource-sharing\n(no shared territory)",
                            "Medium",
                            "High resource-sharing\n(year-round territory)"),
                          each = 2),
                      levels = c("Low resource-sharing\n(no shared territory)",
                                 "Medium",
                                 "High resource-sharing\n(year-round territory)"))
  g$trait <- trait_label
  g
}

fig2_df <- bind_rows(
  marginal_grid("dichromatism", "Plumage dichromatism (male - female)"),
  marginal_grid("display_num",  "Display agility (1-5)")
)
fig2_df$trait <- factor(fig2_df$trait,
                        levels = c("Plumage dichromatism (male - female)",
                                   "Display agility (1-5)"))

fig2 <- ggplot(fig2_df, aes(dim_bin, pred, group = res_lab, colour = res_lab)) +
  geom_line(linewidth = 1.1) + geom_point(size = 3) +
  facet_wrap(~ trait, scales = "free_y") +
  scale_colour_manual(values = c("#D55E00", "grey55", "#0072B2"), name = NULL) +
  labs(
    title = "Both choice-related traits depend on between-sex resource sharing",
    subtitle = "Predicted trait value from the PGLS interaction model (tree 1)",
    x = "Mating-arena dimensionality (lifestyle proxy)",
    y = "Predicted trait value",
    caption = paste0(
      "In both panels, a 3D lifestyle raises the trait where the sexes share ",
      "little on a common territory, but not where they share resources ",
      "year-round.\nY-axes are trait-specific and not comparable in magnitude.")) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey50", size = 8),
        strip.text = element_text(face = "bold"),
        legend.position = "right", panel.grid.minor = element_blank())

ggsave(file.path(dir_out_explo, "fig2_interaction.pdf"), fig2, width = 9, height = 4.4)
cat("  saved fig2_interaction.pdf\n")


################################################################################
## FIGURE 4: care moderation — dim_bin x developmental mode (dev_pc1) ---------
## Exploratory. Predicted trait value in 2D vs 3D across the developmental-mode
## axis (hatchling PC1: low = altricial/high care, high = precocial/low care).
## Prediction: the 2D->3D gap in male ornament / display WIDENS toward the
## precocial (low-care) end. Fit on tree 1, raw dev_pc1 (prediction is invariant
## to moderator scaling). NB dev_pc1 is compressed for passerine-only plumage.
################################################################################
cat("\n=== Figure 4 (developmental-mode interaction) ===\n")

marginal_dev <- function(response, base, trait_label) {
  d <- base %>%
    filter(!is.na(dim_bin), !is.na(.data[[response]]),
           !is.na(dev_pc1), !is.na(tip_label)) %>%
    distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
  rownames(d) <- d$tip_label
  tr  <- ape::drop.tip(trees[[1]], setdiff(trees[[1]]$tip.label, rownames(d)))
  d1  <- d[tr$tip.label, ]
  fit <- phylolm::phylolm(as.formula(paste(response, "~ dim_bin * dev_pc1")),
                          data = d1, phy = tr, model = "lambda")
  co  <- coef(fit)
  gx  <- seq(quantile(d1$dev_pc1, 0.05), quantile(d1$dev_pc1, 0.95),
             length.out = 60)
  g   <- expand.grid(dim_bin = factor(DIM_LEVELS, levels = DIM_LEVELS), dev_pc1 = gx)
  g$d3   <- as.integer(g$dim_bin == "3D")
  g$pred <- co["(Intercept)"] + co[COEF_3D] * g$d3 +
            co["dev_pc1"] * g$dev_pc1 +
            co[paste0(COEF_3D, ":dev_pc1")] * g$d3 * g$dev_pc1
  g$trait <- trait_label
  g
}

fig4_df <- bind_rows(
  marginal_dev("male_plumage", avo_base, "Male plumage elaboration"),
  marginal_dev("dichromatism", avo_base, "Plumage dichromatism (M-F)"),
  marginal_dev("display_num",  merged,   "Display agility (1-5)")
)
fig4_df$trait <- factor(fig4_df$trait,
                        levels = c("Male plumage elaboration",
                                   "Plumage dichromatism (M-F)",
                                   "Display agility (1-5)"))

fig4 <- ggplot(fig4_df, aes(dev_pc1, pred, colour = dim_bin)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ trait, scales = "free") +
  scale_colour_manual(values = c("2D" = "#D55E00", "3D" = "#009E73"), name = NULL) +
  labs(
    title = "Does the 3D effect on choice traits depend on developmental mode?",
    subtitle = "Predicted trait value from the PGLS interaction model (tree 1); exploratory",
    x = "Developmental mode (hatchling PC1):  altricial / long care  \u2192  precocial / short care",
    y = "Predicted trait value",
    caption = paste0(
      "Prediction: the 2D\u21923D gap widens toward the precocial (low-care) end. ",
      "Direction is consistent for male plumage and display (not for the M-F difference),\n",
      "but not statistically robust. dev_pc1 is compressed for passerine-only plumage; ",
      "y-axes are trait-specific and not comparable in magnitude.")) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey50", size = 8),
        strip.text = element_text(face = "bold"),
        legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(dir_out_explo, "fig4_devmode_interaction.pdf"), fig4,
       width = 11, height = 4.2)
cat("  saved fig4_devmode_interaction.pdf\n")

cat("
99_exploratory.R done.
")
