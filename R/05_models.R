################################################################################
## 05_models.R
## Every model reported in the manuscript, fitted across the 100 trees, and saved
## to output/model_results.rds for 06_figures.R and 07_tables.R.
##
## Requires 00_setup.R (engine + helpers), 03_load_dataset.R (merged, avo_base),
## 04_trees.R (trees). Runtime: ~30-40 min on a laptop (100 trees x ~70 models).
##
## Sections:
##   A. Block A -- main (bivariate) effects of 3D, by block:     [main text]
##        sexual selection & mating system | contest | choice
##   C. Block A with controls / on subsets                       [appendix]
##        C2 Rensch body-size control,
##        C4 display excl. aerial / + HWI,
##        C5 Passeriformes only
##   B. Block B2 -- dim_bin x continuous parental care           [appendix]
##   D. Descriptive diagnostics quoted in the text (aerial leks)
##
## Analyses not in the paper (resource-sharing and intensity interactions,
## binary care / developmental mode, flightlessness, HWI proxy) are kept in
## R/99_exploratory.R, which run_all.R does not source.
################################################################################

## Short source keys, resolved to citations in 07_tables.R (`sources`).
## One spec row per Block A response. Columns:
##   block     : ss (sexual selection & mating system) | contest | choice
##   type      : gaussian (PGLS lambda) | binary (penalised phylogenetic logistic)
##   direction : a priori prediction for the 3D effect (none = no prediction)
##   base      : backbone the response comes from (avo_base = all species;
##               merged = Lislevand subset)
##   separated : rare binary response with (quasi-)complete separation between
##               2D and 3D -> fitted with the penalised model anyway, but also
##               summarised by an exact Fisher test (descriptive)
##   in_fig    : shown in Figure 1 (all rows appear in Table S1)
specA <- tibble::tribble(
  ~id,          ~block,    ~response,        ~type,      ~direction, ~base,      ~source,     ~separated, ~in_fig, ~label,
  "intensity",  "ss",      "barber_ss",      "gaussian", "none",     "avo_base", "barber",    FALSE,      TRUE,    "Sexual-selection intensity",
  "rdp_marc",   "ss",      "rdp",            "binary",   "negative", "avo_base", "marcondes", FALSE,      TRUE,    "Resource-defence polygamy vs monogamy",
  "poly_barb",  "ss",      "barber_poly",    "binary",   "negative", "avo_base", "barber",    FALSE,      TRUE,    "Strong polygamy vs monogamy",
  "lek_marc",   "ss",      "lek_m",          "binary",   "positive", "avo_base", "marcondes", FALSE,      TRUE,    "Lek vs monogamy",
  "lek_barb",   "ss",      "barber_lek",     "binary",   "positive", "avo_base", "barber",    FALSE,      TRUE,    "Lek vs monogamy",
  "lek_vs_rdp", "ss",      "lek_vs_rdp",     "binary",   "positive", "avo_base", "marcondes", FALSE,      FALSE,   "Lek vs resource-defence polygamy",
  "ssd_mass",   "contest", "ssd_mass",       "gaussian", "negative", "merged",   "lislevand", FALSE,      TRUE,    "Mass dimorphism",
  "ssd_tarsus", "contest", "ssd_tarsus",     "gaussian", "negative", "merged",   "lislevand", FALSE,      TRUE,    "Tarsus dimorphism",
  "ssd_wing",   "contest", "ssd_wing",       "gaussian", "negative", "merged",   "lislevand", FALSE,      TRUE,    "Wing dimorphism",
  "ssd_bill",   "contest", "ssd_bill",       "gaussian", "negative", "merged",   "lislevand", FALSE,      TRUE,    "Bill dimorphism",
  "ssd_tail",   "contest", "ssd_tail",       "gaussian", "negative", "merged",   "lislevand", FALSE,      TRUE,    "Tail dimorphism",
  "spur",       "contest", "spur_hi",        "binary",   "negative", "avo_base", "menezes",   TRUE,       FALSE,   "Bony spurs",
  "dichrom",    "choice",  "dichromatism",   "gaussian", "positive", "avo_base", "dale",      FALSE,      TRUE,    "Plumage dichromatism (M $-$ F)",
  "male_pl",    "choice",  "male_plumage",   "gaussian", "positive", "avo_base", "dale",      FALSE,      TRUE,    "Male plumage elaboration",
  "female_pl",  "choice",  "female_plumage", "gaussian", "positive", "avo_base", "dale",      FALSE,      TRUE,    "Female plumage elaboration",
  "dichro_sr",  "choice",  "dichro_sr",      "gaussian", "positive", "avo_base", "sexrole",   FALSE,      TRUE,    "Plumage dimorphism (all orders)",
  "uv",         "choice",  "uv_cd",          "gaussian", "positive", "avo_base", "uv",        FALSE,      TRUE,    "UV-inclusive colour dimorphism",
  "display",    "choice",  "display_num",    "gaussian", "positive", "merged",   "lislevand", FALSE,      TRUE,    "Display agility",
  "allopreen",  "choice",  "allopreen",      "binary",   "positive", "avo_base", "kenny",     FALSE,      TRUE,    "Allopreening"
)

## fit one spec row, optionally on a subset (quosure) of its backbone.
## Rare binary responses are skipped when separated (rarest positive cell below
## SEP_THRESHOLD) unless force = TRUE (the penalised model is then fitted anyway).
fit_spec <- function(s, subset = NULL, force = FALSE, tag = "") {
  base <- get(s$base)
  if (!is.null(subset)) base <- base %>% filter(!!subset)
  if (s$type == "binary" && !force) {
    tb <- xtab_dim(base, s$response)
    if (min_pos_cell(tb) < SEP_THRESHOLD) {
      cat(sprintf("\n== %s%s: (quasi-)separation (rarest positive cell = %d) -> not fitted ==\n",
                  s$label, tag, min_pos_cell(tb)))
      return(NULL)
    }
  }
  run_pgls(s$response, type = s$type, direction = s$direction,
           dat = .prep_dat(base, s$response), label = paste0(s$label, tag))
}

################################################################################
## A. BLOCK A: MAIN (BIVARIATE) EFFECTS OF 3D --------------------------------
## response ~ dim_bin, one model per response. Predictions:
##   contest traits (SSD, spurs, RDP) DOWN in 3D;
##   choice traits (plumage, display, lekking) UP in 3D;
##   intensity: no directional prediction (form, not intensity).
################################################################################
cat("\n########## A. BLOCK A: MAIN EFFECTS ##########\n")

resA <- lapply(seq_len(nrow(specA)), function(i)
  fit_spec(specA[i, ], force = specA$separated[i]))
names(resA) <- specA$id

## exact Fisher tests for the separated binaries (descriptive, non-phylogenetic)
fisherA <- lapply(specA$id[specA$separated], function(id) {
  s  <- specA[specA$id == id, ]
  tb <- xtab_dim(get(s$base), s$response)
  cat(sprintf("\n%s by dimensionality:\n", s$label)); print(tb)
  ft <- fisher_dim(tb, s$label)
  list(tab = tb, or = unname(ft$estimate), p = ft$p.value)
})
names(fisherA) <- specA$id[specA$separated]

cat(sprintf("\nDecomposition check: beta(male) - beta(female) = %+.4f  vs  beta(dichrom) = %+.4f\n",
            resA$male_pl$est - resA$female_pl$est, resA$dichrom$est))

################################################################################
## C. BLOCK A WITH CONTROLS / ON SUBSETS (appendix) --------------------------
## Each entry: list(ref = bivariate reference, ctl = controlled/subset model,
## control = description). For covariate controls, `ref` is the bivariate
## model refitted on the same sample (species with the covariate); for subset
## analyses, `ref` is the main Block A model on the full sample.
################################################################################
cat("\n########## C. CONTROLS ##########\n")
spec_row <- function(id) specA[specA$id == id, ]
resC <- list()

## (Sex-role reversal and polyandry, and the intensity model excluding
##  sex-role-reversed species, were dropped from the paper: expectations and
##  mechanism too uncertain. No C1 section.)

## C2. Rensch: body-size control (+ log mean body mass) for each SSD trait
for (v in c("ssd_mass", "ssd_tarsus", "ssd_wing", "ssd_bill", "ssd_tail")) {
  r <- run_pgls_ctrl(v, covar = "body_size_log", direction = "negative",
                     label = spec_row(v)$label, type = "gaussian", data = merged)
  resC[[paste0(v, "_rensch")]] <- list(id = v, control = "+ log body mass",
                                       ref = r$raw, ctl = r$ctl)
}

## (C3, tail dimorphism controlled for mass dimorphism, was dropped from the
##  paper.)

## C4. display agility: drop Aerial species (partial definitional overlap with
##     3D), and control for flight efficiency (hand-wing index)
resC$display_noaerial <- list(
  id = "display", control = "Excluding aerial species",
  ref = resA$display,
  ctl = fit_spec(spec_row("display"),
                 subset = rlang::quo(primary_lifestyle %in% c("Insessorial", "Terrestrial")),
                 tag = " [excl. aerial]"))
r <- run_pgls_ctrl("display_num", covar = "hwi", direction = "positive",
                   label = "Display agility", type = "gaussian", data = merged)
resC$display_hwi <- list(id = "display", control = "+ hand-wing index",
                         ref = r$raw, ctl = r$ctl)

## C5. Passeriformes only: every Figure-1 response refitted within the one
##     order large enough for a within-clade test (separated binaries skipped)
cat("\n--- C5. Passeriformes only ---\n")
resP <- lapply(which(specA$in_fig), function(i)
  fit_spec(specA[i, ], subset = rlang::quo(order == "Passeriformes"),
           tag = " [Passeriformes]"))
names(resP) <- specA$id[specA$in_fig]
for (id in names(resP))
  resC[[paste0(id, "_passer")]] <- list(id = id, control = "Passeriformes only",
                                        ref = resA[[id]], ctl = resP[[id]])

################################################################################
## B. BLOCK B2: dim_bin x CONTINUOUS PARENTAL CARE (appendix) -----------------
## trait ~ dim_bin * care_z, care_z = z-scored relative care investment of the
## sexes (care_cont > 0 = female-biased care; see sign check in 02_build).
## Channel-switch prediction for the interaction dim_bin3D:care_z:
##   contest traits (SSD)                  -> negative (3D drop larger under female care)
##   choice traits (plumage, display, ...) -> positive (3D ornament larger under female care)
################################################################################
cat("\n########## B. CONTINUOUS-CARE MODERATION (B2) ##########\n")

specB2 <- tibble::tribble(
  ~id,          ~block,    ~response,      ~type,      ~direction, ~base,
  "ssd_mass",   "contest", "ssd_mass",     "gaussian", "negative", "merged",
  "ssd_tarsus", "contest", "ssd_tarsus",   "gaussian", "negative", "merged",
  "ssd_wing",   "contest", "ssd_wing",     "gaussian", "negative", "merged",
  "ssd_bill",   "contest", "ssd_bill",     "gaussian", "negative", "merged",
  "ssd_tail",   "contest", "ssd_tail",     "gaussian", "negative", "merged",
  "dichrom",    "choice",  "dichromatism", "gaussian", "positive", "avo_base",
  "male_pl",    "choice",  "male_plumage", "gaussian", "positive", "avo_base",
  "dichro_sr",  "choice",  "dichro_sr",    "gaussian", "positive", "avo_base",
  "display",    "choice",  "display_num",  "gaussian", "positive", "merged",
  "allopreen",  "choice",  "allopreen",    "binary",   "positive", "avo_base"
) %>%
  left_join(specA %>% select(id, source, label), by = "id")

resB2 <- lapply(seq_len(nrow(specB2)), function(i) {
  s <- specB2[i, ]
  d <- get(s$base) %>%
    filter(!is.na(dim_bin), !is.na(.data[[s$response]]), !is.na(care_cont),
           !is.na(tip_label)) %>%
    mutate(care_z = as.numeric(scale(care_cont))) %>%
    distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
  rownames(d) <- d$tip_label
  run_pgls(s$response, type = s$type, direction = s$direction,
           coef_name = paste0(COEF_3D, ":care_z"),
           formula   = as.formula(paste(s$response, "~ dim_bin * care_z")),
           dat = d, label = paste0(s$label, " x care"))
})
names(resB2) <- specB2$id

################################################################################
## D. DESCRIPTIVE DIAGNOSTICS quoted in the text ------------------------------
## The aerial lek signal is concentrated in hummingbirds (supports a Discussion
## statement), rather than spread across families.
################################################################################
cat("\n########## D. DIAGNOSTICS ##########\n")

aerial_lek_fam <- avo_base %>%
  filter(primary_lifestyle == "Aerial", marc_lek == 1) %>%
  distinct(key, .keep_all = TRUE) %>%
  count(family, order, sort = TRUE)
cat("\nAerial lekking species by family (Marcondes):\n")
print(as.data.frame(aerial_lek_fam), row.names = FALSE)
cat(sprintf("Total aerial lek: %d | Trochilidae: %d | distinct families: %d\n",
            sum(aerial_lek_fam$n),
            sum(aerial_lek_fam$n[aerial_lek_fam$family == "Trochilidae"]),
            nrow(aerial_lek_fam)))

################################################################################
## SAVE -----------------------------------------------------------------------
################################################################################
model_results <- list(specA = specA, resA = resA, fisherA = fisherA,
                      resC = resC, specB2 = specB2, resB2 = resB2,
                      aerial_lek_fam = aerial_lek_fam,
                      n_trees = length(trees), date = Sys.time())
saveRDS(model_results, path_model_cache)
cat("\n05_models.R done. Results saved to", path_model_cache, "\n")
