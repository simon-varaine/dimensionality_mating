################################################################################
## 01_load_raw.R  [AUTHORS ONLY]
## Load and clean the six raw sources into tidy in-memory tables:
##   lis      Lislevand et al. 2007  : SSD, mating system, display, resource
##   avo      AVONET / BirdTree      : primary lifestyle (2D/3D proxy), HWI, mass
##   blio     BLIOCPhyloMasterTax    : taxonomy <-> tree tip labels
##   dale     Dale et al. 2015       : male/female plumage scores (dichromatism)
##   spur_sp  Menezes & Palaoro 2022 : bony-spur presence
##   marc_sp  Marcondes & Douvas 2024: mating system (M / P / L)
## No joins and no derived variables here -- that is 02_build_dataset.R.
## Requires 00_setup.R (paths, clean_binom, get_col).
################################################################################

## --- 1a. Lislevand 2007 (tab-delimited, NA coded as -999) --------------------
lis <- read_tsv(path_lislevand,
                na = c("", "NA", "-999", "-999.0"),
                show_col_types = FALSE) %>%
  clean_names() %>%
  transmute(
    sci_lis       = species_name,                 # binomial (Monroe & Sibley 1997)
    english_name  = english_name,
    mating_system = as.integer(mating_system),    # 1..5 social mating system
    display       = as.integer(display),          # 1..5 display agility
    resource      = as.integer(resource),         # 0..2 between-sex resource sharing
    ## sex-specific measures for the SSD ratios
    m_mass   = suppressWarnings(as.numeric(m_mass)),
    f_mass   = suppressWarnings(as.numeric(f_mass)),
    m_wing   = suppressWarnings(as.numeric(m_wing)),
    f_wing   = suppressWarnings(as.numeric(f_wing)),
    m_tarsus = suppressWarnings(as.numeric(m_tarsus)),
    f_tarsus = suppressWarnings(as.numeric(f_tarsus)),
    m_tail   = suppressWarnings(as.numeric(m_tail)),
    f_tail   = suppressWarnings(as.numeric(f_tail)),
    m_bill   = suppressWarnings(as.numeric(m_bill)),
    f_bill   = suppressWarnings(as.numeric(f_bill))
  ) %>%
  ## any residual negatives in the morphological measures -> NA
  mutate(across(c(starts_with("m_"), starts_with("f_")),
                ~ ifelse(. < 0, NA_real_, .))) %>%
  mutate(key = clean_binom(sci_lis)) %>%
  filter(!is.na(key), key != "") %>%
  distinct(key, .keep_all = TRUE)

cat("Lislevand:", nrow(lis), "binomials;",
    sum(!is.na(lis$mating_system)), "with mating_system\n")

## --- 1b. AVONET / BirdTree (robust column detection) -------------------------
avo_raw <- read_excel(path_avonet, sheet = avonet_sheet) %>% clean_names()

c_sp   <- get_col(avo_raw, "^species")    # species3
c_fam  <- get_col(avo_raw, "^family")
c_ord  <- get_col(avo_raw, "^order")
c_life <- get_col(avo_raw, "lifestyle")   # primary_lifestyle
c_mass <- get_col(avo_raw, "^mass$")
c_hwi  <- get_col(avo_raw, "hand.?wing")  # hand_wing_index
stopifnot(!is.na(c_sp), !is.na(c_life))

avo <- avo_raw %>%
  transmute(
    sci_avo           = .data[[c_sp]],
    family_avo        = if (!is.na(c_fam))  .data[[c_fam]]  else NA,
    order_avo         = if (!is.na(c_ord))  .data[[c_ord]]  else NA,
    primary_lifestyle = .data[[c_life]],
    avo_mass          = if (!is.na(c_mass)) as.numeric(.data[[c_mass]]) else NA_real_,
    hwi               = if (!is.na(c_hwi))  as.numeric(.data[[c_hwi]])  else NA_real_
  ) %>%
  mutate(key = clean_binom(sci_avo)) %>%
  filter(!is.na(key), key != "") %>%
  distinct(key, .keep_all = TRUE)

cat("AVONET:", nrow(avo), "binomials;",
    sum(!is.na(avo$hwi)), "with hand-wing index\n")

## --- 1c. BLIOCPhyloMasterTax (taxonomy <-> tree tip labels) ------------------
blio <- read_csv(path_blio, show_col_types = FALSE) %>%
  clean_names() %>%
  transmute(
    key_blio  = clean_binom(scientific),
    tip_label = tip_label     # "Genus_species" (matches the tree tips)
  )
cat("BLIOCPhyloMasterTax:", nrow(blio), "species\n")

## --- 1d. Dale 2015 (plumage scores; male/female kept for the M/F split) ------
dale <- read_csv(path_dale, show_col_types = FALSE) %>%
  clean_names() %>%
  transmute(
    tip_label      = tip_label,   # already "Genus_species"
    male_plumage   = male_plumage_score,
    female_plumage = female_plumage_score,
    dichromatism   = male_plumage_score - female_plumage_score
  ) %>%
  filter(!is.na(dichromatism), !is.na(tip_label)) %>%
  distinct(tip_label, .keep_all = TRUE)
cat("Dale:", nrow(dale), "species\n")

## --- 1e. Menezes & Palaoro 2022 (bony-spur presence) ------------------------
## The manuscript uses only the high-confidence presence (descriptive + Fisher).
if (!file.exists(path_spur))
  stop("Spur file not found: ", path_spur)
spur_sp <- read_csv(path_spur, show_col_types = FALSE) %>%
  clean_names() %>%
  transmute(
    key     = clean_binom(scientific),
    spur_hi = as.integer(any_spur_presence_high_conf)
  ) %>%
  filter(!is.na(key)) %>%
  distinct(key, .keep_all = TRUE)
cat("Spurs:", nrow(spur_sp), "species (",
    sum(spur_sp$spur_hi == 1, na.rm = TRUE), "with a spur )\n")

## --- 1f. Marcondes & Douvas 2024 (independent mating-system source) ----------
## Coding on the pair bond: M = monogamy, P = resource-defense polygamy,
## L = lekking; U / PU excluded (-> NA). Two binary contrasts:
##   marc_poly = P vs (M, L) ; marc_lek = L vs (M, P).
if (!file.exists(path_marcondes))
  stop("Marcondes file not found: ", path_marcondes)
marc_sp <- read_excel(path_marcondes, sheet = marcondes_sheet) %>%
  clean_names() %>%
  transmute(
    key         = clean_binom(species),
    ms_marc_raw = str_squish(as.character(mating_system))  # "M","P","L","U","PU"
  ) %>%
  filter(!is.na(key), key != "") %>%
  mutate(
    marc_poly = case_when(ms_marc_raw == "P"          ~ 1L,
                          ms_marc_raw %in% c("M", "L") ~ 0L,
                          TRUE ~ NA_integer_),
    marc_lek  = case_when(ms_marc_raw == "L"          ~ 1L,
                          ms_marc_raw %in% c("M", "P") ~ 0L,
                          TRUE ~ NA_integer_)
  ) %>%
  distinct(key, .keep_all = TRUE)
cat("Marcondes:", nrow(marc_sp), "species | P =",
    sum(marc_sp$marc_poly == 1, na.rm = TRUE), "| L =",
    sum(marc_sp$marc_lek == 1, na.rm = TRUE), "\n")

cat("01_load_raw.R done.\n")
