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
                          TRUE ~ NA_integer_),
    ## raw 3-state category kept for the monogamy-referenced contrasts (§2.3)
    marc_system = case_when(ms_marc_raw %in% c("M", "P", "L") ~ ms_marc_raw,
                            TRUE ~ NA_character_)
  ) %>%
  distinct(key, .keep_all = TRUE)
cat("Marcondes:", nrow(marc_sp), "species | P =",
    sum(marc_sp$marc_poly == 1, na.rm = TRUE), "| L =",
    sum(marc_sp$marc_lek == 1, na.rm = TRUE), "\n")

## --- 1g. Barber et al. 2024 (sexual-selection intensity, BirdTree) -----------
## Sexual selection: ordinal 0 (strict monogamy) .. 4 (lekking / display posts),
##   inferred from mating system + EPP + OSR. Broad coverage (~all birds).
## Sex-role reversal: 0/1, used for the "polyandry more frequent in 2D" test.
if (!file.exists(path_barber))
  stop("Barber file not found: ", path_barber)
barber <- read_excel(path_barber, sheet = barber_sheet) %>%
  clean_names() %>%
  transmute(
    key        = clean_binom(scientific_name_bird_tree),
    barber_ss  = suppressWarnings(as.integer(sexual_selection)),   # 0..4 intensity
    barber_srr = suppressWarnings(as.integer(sex_role_reversal))   # 0/1
  ) %>%
  filter(!is.na(key), key != "") %>%
  distinct(key, .keep_all = TRUE)
cat("Barber:", nrow(barber), "species | with SS score:",
    sum(!is.na(barber$barber_ss)), "| sex-role reversed:",
    sum(barber$barber_srr == 1, na.rm = TRUE), "\n")

## --- 1h. Who cares? parental-care type (Cockburn/HBW; Jetz taxonomy) ---------
## pc.mode: P pair/biparental, M male-only, F female-only, C cooperative, N none.
## This is the PROXIMAL (near-circular) care moderator -> upper-bound test.
if (!file.exists(path_carewho))
  stop("Who cares? file not found: ", path_carewho)
carew_raw <- read_csv(path_carewho, show_col_types = FALSE) %>% clean_names()
cw_sp   <- get_col(carew_raw, "species")
cw_mode <- get_col(carew_raw, "pc.?mode")
stopifnot(!is.na(cw_sp), !is.na(cw_mode))
carew <- carew_raw %>%
  transmute(key       = clean_binom(.data[[cw_sp]]),
            care_mode = str_squish(as.character(.data[[cw_mode]]))) %>%
  filter(!is.na(key), key != "", care_mode %in% c("P", "M", "F", "C", "N")) %>%
  distinct(key, .keep_all = TRUE)
cat("Who cares?:", nrow(carew), "species | ",
    paste(names(table(carew$care_mode)), table(carew$care_mode), sep = "=", collapse = " "), "\n")

## --- 1i. Developmental mode (Cooney et al. 2021; hatchling/chick PCA) --------
## Continuous altricial<->precocial axis. DISTAL, clean care moderator.
## /!\ sign (which end = precocial) to be verified after the build.
if (!file.exists(path_devmode))
  stop("Dev-mode file not found: ", path_devmode)
dev_raw <- read_excel(path_devmode, sheet = devmode_sheet) %>% clean_names()
dv_sp <- get_col(dev_raw, "species")
dv_h1 <- get_col(dev_raw, "hatchling.*pc.?1")
dv_c1 <- get_col(dev_raw, "chick.*pc.?1")
stopifnot(!is.na(dv_sp), !is.na(dv_h1))
devmode <- dev_raw %>%
  transmute(key          = clean_binom(.data[[dv_sp]]),
            dev_pc1      = suppressWarnings(as.numeric(.data[[dv_h1]])),   # hatchling PC1
            dev_chickpc1 = if (!is.na(dv_c1)) suppressWarnings(as.numeric(.data[[dv_c1]])) else NA_real_) %>%
  filter(!is.na(key), key != "") %>%
  distinct(key, .keep_all = TRUE)
cat("Dev mode:", nrow(devmode), "species | with hatchling PC1:",
    sum(!is.na(devmode$dev_pc1)), "\n")

## --- 1k. Volancy (Sayol et al. 2020): flightlessness (hard, rare 2D proxy) ---
if (!file.exists(path_volancy)) stop("Volancy file not found: ", path_volancy)
vol_raw <- read_excel(path_volancy) %>% clean_names()
volancy <- vol_raw %>%
  transmute(key        = clean_binom(.data[[get_col(vol_raw, "^species$")]]),
            flightless = as.integer(str_detect(str_to_lower(.data[[get_col(vol_raw, "volancy")]]),
                                                "flightless"))) %>%
  filter(!is.na(key), key != "") %>% distinct(key, .keep_all = TRUE)
cat("Volancy:", nrow(volancy), "species | flightless:",
    sum(volancy$flightless == 1, na.rm = TRUE), "\n")

## --- 1l. Allopreening + care traits (Kenny et al. 2017) ----------------------
## allopreen (0/1): a pair-bond / "patience" signal (mutual-choice channel).
## par_coop (continuous): parental cooperation. age_indep (d): care duration.
if (!file.exists(path_allopreen)) stop("Allopreening file not found: ", path_allopreen)
allo_raw <- read_excel(path_allopreen, sheet = allopreen_sheet) %>% clean_names()
allo <- allo_raw %>%
  transmute(key       = clean_binom(.data[[get_col(allo_raw, "^species$")]]),
            allopreen = suppressWarnings(as.integer(.data[[get_col(allo_raw, "allopreen_pairs")]])),
            par_coop  = suppressWarnings(as.numeric(.data[[get_col(allo_raw, "parental_cooperation_score")]])),
            age_indep = suppressWarnings(as.numeric(.data[[get_col(allo_raw, "age_of_independence")]]))) %>%
  filter(!is.na(key), key != "") %>% distinct(key, .keep_all = TRUE)
cat("Allopreening:", nrow(allo), "species | allopreen=1:",
    sum(allo$allopreen == 1, na.rm = TRUE), "| par_coop:",
    sum(!is.na(allo$par_coop)), "| age_indep:", sum(!is.na(allo$age_indep)), "\n")

## --- 1m. UV-inclusive dichromatism (colour discriminability; JZO 2025) -------
if (!file.exists(path_uvdichrom)) stop("UV dichromatism file not found: ", path_uvdichrom)
uv_raw <- read_excel(path_uvdichrom) %>% clean_names()
uvdi <- uv_raw %>%
  transmute(key   = clean_binom(.data[[get_col(uv_raw, "^species$")]]),
            uv_cd = suppressWarnings(as.numeric(.data[[get_col(uv_raw, "colour_discriminability_absolute")]]))) %>%
  filter(!is.na(key), key != "") %>% distinct(key, .keep_all = TRUE)
cat("UV dichromatism:", nrow(uvdi), "species | with CD:", sum(!is.na(uvdi$uv_cd)), "\n")

## --- 1n. Sex-role ecology (Tobias/Szekely): continuous care + broad dichro ---
## care_cont : relative investment of the sexes in parental care (continuous) --
##   a better proximal care moderator than the F-vs-P binary (more coverage/power).
## dichro_sr : plumage dimorphism score, near-complete coverage (~9960, all clades).
if (!file.exists(path_sexrole)) stop("Sex-role file not found: ", path_sexrole)
sr_raw <- read_excel(path_sexrole, sheet = sexrole_sheet) %>% clean_names()
sexrole <- sr_raw %>%
  transmute(key       = clean_binom(.data[[get_col(sr_raw, "^species$")]]),
            care_cont = suppressWarnings(as.numeric(.data[[get_col(sr_raw, "^care$")]])),
            dichro_sr = suppressWarnings(as.numeric(.data[[get_col(sr_raw, "^dichro$")]]))) %>%
  filter(!is.na(key), key != "") %>% distinct(key, .keep_all = TRUE)
cat("Sex-role ecology:", nrow(sexrole), "species | care:",
    sum(!is.na(sexrole$care_cont)), "| dichro:", sum(!is.na(sexrole$dichro_sr)), "\n")

cat("01_load_raw.R done.\n")

## --- 1j. BirdBase (care durations) -- set aside for now (large, slow to read;
## dev_pc1 largely captures the same axis). To re-enable, restore the loader and
## the joins/schema in 02/03, and pre-convert the .xlsx to a small CSV first.
