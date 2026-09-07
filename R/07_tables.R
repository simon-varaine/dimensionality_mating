################################################################################
## 07_tables.R
## Table 1 (raw values by lifestyle) and Tables A1-A4 (sample sizes, robustness,
## Rensch control, display robustness). Each is printed and written to output/.
## Requires 00_setup.R, 03_load_dataset.R, 04_trees.R, 05_models.R.
################################################################################

lifestyle_order <- c("Terrestrial", "Insessorial", "Aerial",
                     "Generalist", "Aquatic")

################################################################################
## TABLE 1: raw values by foraging lifestyle ----------------------------------
## N and avo_base-derived columns (RDP, lek-Marcondes, dichromatism, spurs) from
## avo_base; Lislevand columns (% poly, % lek, SSD mass, display) from merged.
################################################################################
by_life_avo <- avo_base %>%
  filter(!is.na(primary_lifestyle)) %>% distinct(key, .keep_all = TRUE) %>%
  group_by(primary_lifestyle) %>%
  summarise(N            = n(),
            pct_RDP      = round(100 * mean(marc_poly == 1, na.rm = TRUE), 1),
            pct_lek_marc = round(100 * mean(marc_lek  == 1, na.rm = TRUE), 1),
            dichrom_med  = round(median(dichromatism, na.rm = TRUE), 2),
            pct_spur     = round(100 * mean(spur_hi == 1, na.rm = TRUE), 1),
            .groups = "drop")

by_life_lis <- merged %>%
  filter(!is.na(primary_lifestyle)) %>% distinct(key, .keep_all = TRUE) %>%
  group_by(primary_lifestyle) %>%
  summarise(pct_poly     = round(100 * mean(harem == 1, na.rm = TRUE), 1),
            pct_lek_lis  = round(100 * mean(lek   == 1, na.rm = TRUE), 1),
            ssd_mass_med = round(median(ssd_mass,    na.rm = TRUE), 3),
            display_med  = round(median(display_num, na.rm = TRUE), 1),
            .groups = "drop")

table1 <- by_life_avo %>%
  left_join(by_life_lis, by = "primary_lifestyle") %>%
  mutate(primary_lifestyle = factor(primary_lifestyle, levels = lifestyle_order)) %>%
  arrange(primary_lifestyle) %>%
  ## column order as in the manuscript Table 1
  select(primary_lifestyle, N, pct_poly, pct_lek_lis, pct_RDP, pct_lek_marc,
         ssd_mass_med, display_med, dichrom_med, pct_spur)

cat("\n=== TABLE 1: raw values by lifestyle ===\n")
print(as.data.frame(table1), row.names = FALSE)
write_csv(table1, file.path(dir_out, "table1_by_lifestyle.csv"))

################################################################################
## TABLE A1: sample sizes per model (N, N_2D, N_3D) ---------------------------
################################################################################
cnt <- function(data, resp, clade = NULL, extra = NULL) {
  d <- data
  if (!is.null(clade)) d <- dplyr::filter(d, order == clade)
  d <- d %>% dplyr::filter(!is.na(dim_bin), !is.na(.data[[resp]]), !is.na(tip_label))
  if (!is.null(extra)) d <- dplyr::filter(d, !is.na(.data[[extra]]))
  d <- dplyr::distinct(d, tip_label, .keep_all = TRUE)
  c(N = nrow(d),
    N2D = sum(d$dim_bin == "2D"),
    N3D = sum(d$dim_bin == "3D"))
}

tableA1 <- as.data.frame(do.call(rbind, list(
  ## mating system
  "Harem / lek (Lislevand)"             = cnt(merged,   "harem"),
  "Resource-defense polygamy / lekking" = cnt(avo_base, "marc_poly"),
  ## contest-related
  "Mass dimorphism"                     = cnt(merged,   "ssd_mass"),
  "Tarsus dimorphism"                   = cnt(merged,   "ssd_tarsus"),
  "Wing dimorphism"                     = cnt(merged,   "ssd_wing"),
  "Bill dimorphism"                     = cnt(merged,   "ssd_bill"),
  "Tail dimorphism"                     = cnt(merged,   "ssd_tail"),
  "Spur presence"                       = cnt(avo_base, "spur_hi"),
  ## display / choice
  "Display agility (main)"              = cnt(merged,   "display_num"),
  "Display agility x resource"          = cnt(merged,   "display_num", extra = "resource"),
  "Dichromatism (main)"                 = cnt(avo_base, "dichromatism"),
  "Dichromatism x resource"             = cnt(merged,   "dichromatism", extra = "resource"),
  ## within-Passeriformes robustness
  "Harem [Passeriformes]"               = cnt(merged,   "harem",        "Passeriformes"),
  "RDP [Passeriformes]"                 = cnt(avo_base, "marc_poly",    "Passeriformes"),
  "Mass dim. [Passeriformes]"           = cnt(merged,   "ssd_mass",     "Passeriformes"),
  "Tarsus dim. [Passeriformes]"         = cnt(merged,   "ssd_tarsus",   "Passeriformes"),
  "Wing dim. [Passeriformes]"           = cnt(merged,   "ssd_wing",     "Passeriformes"),
  "Bill dim. [Passeriformes]"           = cnt(merged,   "ssd_bill",     "Passeriformes"),
  "Tail dim. [Passeriformes]"           = cnt(merged,   "ssd_tail",     "Passeriformes"),
  "Display [Passeriformes]"             = cnt(merged,   "display_num",  "Passeriformes"),
  "Dichromatism [Passeriformes]"        = cnt(avo_base, "dichromatism", "Passeriformes")
)))
tableA1$response <- rownames(tableA1); rownames(tableA1) <- NULL
tableA1 <- tableA1[, c("response", "N", "N2D", "N3D")]

cat("\n=== TABLE A1: sample sizes per model ===\n")
print(tableA1, row.names = FALSE)
write_csv(tableA1, file.path(dir_out, "tableA1_sample_sizes.csv"))

################################################################################
## TABLE A2: statistical robustness of the Fig. 1 effects ---------------------
################################################################################
row_from_res <- function(res, label) {
  if (is.null(res)) return(NULL)
  sgn <- if (res$direction == "negative") -1 else 1
  data.frame(
    response = label, N = res$n,
    median   = round(median(res$coefs), 4),
    lo       = round(unname(quantile(res$coefs, 0.025)), 4),
    hi       = round(unname(quantile(res$coefs, 0.975)), 4),
    pct_sign = round(100 * mean(res$coefs * sgn > 0)),
    pct_p05  = res$pct_sig,
    pct_p10  = res$pct_sig10,
    rubin_p  = signif(res$p_rubin, 3)
  )
}
tableA2 <- bind_rows(
  row_from_res(res_harem,        "Harem polygyny (Lislevand)"),
  row_from_res(res_marc_poly,    "Resource-defense polygamy (Marcondes)"),
  row_from_res(res_lek,          "Lek / promiscuity (Lislevand)"),
  row_from_res(res_marc_lek,     "Lekking (Marcondes)"),
  row_from_res(res_ssd_mass,     "Mass dimorphism"),
  row_from_res(res_ssd_tarsus,   "Tarsus dimorphism"),
  row_from_res(res_ssd_wing,     "Wing dimorphism"),
  row_from_res(res_ssd_bill,     "Bill dimorphism"),
  row_from_res(res_ssd_tail,     "Tail dimorphism"),
  row_from_res(res_display,      "Display agility"),
  row_from_res(res_dichrom_bigN, "Plumage dichromatism (main effect)")
)
cat("\n=== TABLE A2: statistical robustness ===\n")
print(as.data.frame(tableA2), row.names = FALSE)
write_csv(tableA2, file.path(dir_out, "tableA2_robustness.csv"))

################################################################################
## TABLE A3: body-size (Rensch) control ---------------------------------------
################################################################################
rensch_row <- function(v, label) {
  r <- res_rensch[[v]]
  if (is.null(r) || is.null(r$raw) || is.null(r$ctl)) return(NULL)
  data.frame(
    trait   = label, N = r$raw$n,
    beta_raw = round(median(r$raw$coefs), 4),
    raw_lo   = round(unname(quantile(r$raw$coefs, 0.025)), 4),
    raw_hi   = round(unname(quantile(r$raw$coefs, 0.975)), 4),
    beta_size_ctrl = round(median(r$ctl$coefs), 4),
    ctl_lo   = round(unname(quantile(r$ctl$coefs, 0.025)), 4),
    ctl_hi   = round(unname(quantile(r$ctl$coefs, 0.975)), 4)
  )
}
tableA3 <- bind_rows(
  rensch_row("ssd_mass",   "Mass dimorphism"),
  rensch_row("ssd_tarsus", "Tarsus dimorphism"),
  rensch_row("ssd_tail",   "Tail dimorphism"),
  rensch_row("ssd_wing",   "Wing dimorphism"),
  rensch_row("ssd_bill",   "Bill dimorphism")
)
cat("\n=== TABLE A3: Rensch body-size control ===\n")
print(as.data.frame(tableA3), row.names = FALSE)
write_csv(tableA3, file.path(dir_out, "tableA3_rensch.csv"))

################################################################################
## TABLE A4: display agility robustness and modulation ------------------------
################################################################################
a4_row <- function(res, label) {
  if (is.null(res)) return(NULL)
  data.frame(model = label, N = res$n,
             beta = round(median(res$coefs), 3),
             lo   = round(unname(quantile(res$coefs, 0.025)), 3),
             hi   = round(unname(quantile(res$coefs, 0.975)), 3))
}
tableA4 <- bind_rows(
  a4_row(res_display,             "Display ~ 3D (main)"),
  a4_row(res_display_excl_aerial, "excluding aerial species"),
  a4_row(res_display_hwi$ctl,     "controlling for hand-wing index"),
  a4_row(res_display_interact,    "Display ~ 3D x resource (interaction)")
)
cat("\n=== TABLE A4: display agility robustness ===\n")
print(as.data.frame(tableA4), row.names = FALSE)
write_csv(tableA4, file.path(dir_out, "tableA4_display.csv"))

## spur descriptive result (reported in text / Table 1)
if (exists("spur_fisher") && !is.null(spur_fisher))
  cat(sprintf("\nSpurs 2D vs 3D: Fisher exact OR=%.3f, p=%.3g\n",
              unname(spur_fisher$estimate), spur_fisher$p.value))

cat("\n07_tables.R done. Tables written to", dir_out, "\n")
