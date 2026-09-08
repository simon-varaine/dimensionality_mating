################################################################################
## 02_build_dataset.R  [AUTHORS ONLY]
## Assemble the six cleaned sources into ONE analysis-ready table (one row per
## AVONET species) and write it to data-derived/analysis_data.csv.
##
## Design: a single wide table replaces the former two backbones. AVONET is the
## base (all species); the Lislevand block is left-joined and flagged with
## `lislevand`. Downstream (03_load_dataset.R) reconstructs:
##   avo_base = all species                (spurs, Marcondes, big-N dichromatism)
##   merged   = subset where lislevand==TRUE (harem/lek, SSD, display, resource)
## This keeps the original two-backbone behaviour while sharing a single file.
## Requires 00_setup.R and 01_load_raw.R.
################################################################################

stopifnot(exists("lis"), exists("avo"), exists("blio"),
          exists("dale"), exists("spur_sp"), exists("marc_sp"), exists("barber"))

## --- 2a. base = every AVONET species + taxonomy/tree, colour, spurs, systems -
analysis <- avo %>%
  transmute(
    key, scientific = sci_avo,
    order = order_avo, family = family_avo,
    primary_lifestyle, mass_avonet = avo_mass, hwi
  ) %>%
  left_join(blio %>% select(key_blio, tip_label), by = c("key" = "key_blio")) %>%
  left_join(dale %>% select(tip_label, male_plumage, female_plumage, dichromatism),
            by = "tip_label") %>%
  left_join(spur_sp %>% select(key, spur_hi), by = "key") %>%
  left_join(marc_sp %>% select(key, marc_poly, marc_lek, marc_system), by = "key") %>%
  left_join(barber %>% select(key, barber_ss, barber_srr), by = "key")

## --- 2b. left-join the Lislevand block (subset of species) -------------------
lis_block <- lis %>%
  transmute(
    key, english_name, mating_system,
    display_num = as.numeric(display), resource,
    m_mass, f_mass, m_wing, f_wing, m_tarsus, f_tarsus,
    m_tail, f_tail, m_bill, f_bill,
    lislevand = TRUE
  )
analysis <- analysis %>% left_join(lis_block, by = "key")

## --- 2c. derive all analysis variables --------------------------------------
analysis <- analysis %>%
  mutate(
    lislevand = coalesce(lislevand, FALSE),

    ## dimensionality proxy from AVONET foraging lifestyle.
    ## /!\ foraging lifestyle, not the mating arena itself (documented caveat).
    ##     Aquatic (waterfowl) = reduced-dimensionality arena -> kept apart.
    dim_proxy = case_when(
      primary_lifestyle %in% c("Aerial", "Insessorial") ~ "more_3D",
      primary_lifestyle == "Terrestrial"                ~ "more_2D",
      primary_lifestyle == "Aquatic"                    ~ "aquatic_ambiguous",
      primary_lifestyle == "Generalist"                 ~ "intermediate",
      TRUE ~ NA_character_),
    ## focal contrast (character here; factored in 03_load_dataset.R)
    dim_bin = case_when(dim_proxy == "more_3D" ~ "3D",
                        dim_proxy == "more_2D" ~ "2D",
                        TRUE ~ NA_character_),

    ## mating-system binary contrasts (NA where mating_system is NA)
    ##   harem (4) = polygyny >15%   ; lek (5) = lek/promiscuous
    ##   polyandry (1) = reversed roles
    harem     = if_else(mating_system == 4L, 1L, 0L),
    lek       = if_else(mating_system == 5L, 1L, 0L),
    polyandry = if_else(mating_system == 1L, 1L, 0L),

    ## sexual size dimorphism = log(male/female), trait by trait
    ssd_mass   = ifelse(m_mass   > 0 & f_mass   > 0, log(m_mass   / f_mass),   NA_real_),
    ssd_tarsus = ifelse(m_tarsus > 0 & f_tarsus > 0, log(m_tarsus / f_tarsus), NA_real_),
    ssd_wing   = ifelse(m_wing   > 0 & f_wing   > 0, log(m_wing   / f_wing),   NA_real_),
    ssd_tail   = ifelse(m_tail   > 0 & f_tail   > 0, log(m_tail   / f_tail),   NA_real_),
    ssd_bill   = ifelse(m_bill   > 0 & f_bill   > 0, log(m_bill   / f_bill),   NA_real_),

    ## mean body mass (Lislevand M/F; fallback to AVONET) for the Rensch control
    body_mass_mean = case_when(
      m_mass > 0 & f_mass > 0 ~ (m_mass + f_mass) / 2,
      m_mass > 0             ~ m_mass,
      f_mass > 0             ~ f_mass,
      !is.na(mass_avonet)    ~ mass_avonet,
      TRUE                   ~ NA_real_),
    body_size_log = ifelse(body_mass_mean > 0, log(body_mass_mean), NA_real_),

    ## Barber mating-system contrasts, monogamy-referenced (scores 0-1 = monogamy)
    ##   strong polygamy (2-3) vs monogamy ; lek (4) vs monogamy (mutually excl.).
    ## NB 'strong polygamy' is general polygamy, NOT purely resource-defense.
    barber_poly = case_when(barber_ss %in% c(2L, 3L) ~ 1L,
                            barber_ss %in% c(0L, 1L) ~ 0L,
                            TRUE ~ NA_integer_),
    barber_lek  = case_when(barber_ss == 4L          ~ 1L,
                            barber_ss %in% c(0L, 1L) ~ 0L,
                            TRUE ~ NA_integer_)
  )

## --- 2d. keep the columns used by the analyses (drop raw M/F, mass_avonet) ---
analysis <- analysis %>%
  select(
    ## identifiers
    key, tip_label, scientific, order, family,
    ## predictor
    primary_lifestyle, dim_proxy, dim_bin,
    ## whole-clade variables (avo_base analyses)
    hwi, dichromatism, male_plumage, female_plumage,
    spur_hi, marc_poly, marc_lek, marc_system,
    barber_ss, barber_srr, barber_poly, barber_lek,
    ## Lislevand block (merged analyses)
    lislevand, english_name, mating_system, display_num, resource,
    harem, lek, polyandry,
    ssd_mass, ssd_tarsus, ssd_wing, ssd_tail, ssd_bill,
    body_mass_mean, body_size_log
  )

## --- 2e. write the shared dataset -------------------------------------------
write_csv(analysis, path_analysis_data)
cat("\nShared dataset written:", path_analysis_data,
    "\n  ", nrow(analysis), "species x", ncol(analysis), "columns\n")

## --- 2f. validation anchors (compare to the previous pipeline) --------------
## Model N = species with dim_bin, response and tip_label, deduplicated by tip.
n_anchor <- function(resp, lislevand_only = FALSE) {
  d <- analysis
  if (lislevand_only) d <- dplyr::filter(d, lislevand)
  d %>% dplyr::filter(!is.na(dim_bin), !is.na(tip_label), !is.na(.data[[resp]])) %>%
    dplyr::distinct(tip_label) %>% nrow()
}
cat("\n--- validation anchors (expected values in parentheses) ---\n")
cat(sprintf("  total species            : %d\n", nrow(analysis)))
cat(sprintf("  with tip_label           : %d\n", sum(!is.na(analysis$tip_label))))
cat(sprintf("  lislevand == TRUE        : %d\n", sum(analysis$lislevand)))
cat(sprintf("  harem  (merged)          : %d   (~927)\n",  n_anchor("harem",  TRUE)))
cat(sprintf("  ssd_mass (merged)        : %d   (~2026)\n", n_anchor("ssd_mass", TRUE)))
cat(sprintf("  display_num (merged)     : %d   (~922)\n",  n_anchor("display_num", TRUE)))
cat(sprintf("  dichromatism (avo_base)  : %d   (~5220)\n", n_anchor("dichromatism")))
cat(sprintf("  marc_poly (avo_base)     : %d   (~8278)\n", n_anchor("marc_poly")))
cat(sprintf("  spur_hi (avo_base)       : %d   (~626)\n",  n_anchor("spur_hi")))

cat("\n02_build_dataset.R done.\n")
