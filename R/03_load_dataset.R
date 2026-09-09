################################################################################
## 03_load_dataset.R
## Read the shared analysis dataset and reconstruct the two analysis backbones
## used throughout:
##   avo_base = all species                 (spurs, Marcondes, big-N dichromatism)
##   merged   = subset where lislevand==TRUE (harem/lek, SSD, display, resource)
## This is the entry point for reproducing the analyses: it needs only
## data-derived/analysis_data.csv (no raw data). Requires 00_setup.R.
################################################################################

if (!file.exists(path_analysis_data))
  stop("Shared dataset not found: ", path_analysis_data,
       "\n  Authors: run build_dataset.R first. Others: see README.")

## explicit column schema (also guards against a malformed/renamed file) -------
analysis_data <- read_csv(
  path_analysis_data,
  col_types = cols(
    key               = col_character(),
    tip_label         = col_character(),
    scientific        = col_character(),
    order             = col_character(),
    family            = col_character(),
    primary_lifestyle = col_character(),
    dim_proxy         = col_character(),
    dim_bin           = col_character(),   # factored just below
    hwi               = col_double(),
    dichromatism      = col_double(),
    male_plumage      = col_double(),
    female_plumage    = col_double(),
    spur_hi           = col_integer(),
    marc_poly         = col_integer(),
    marc_lek          = col_integer(),
    marc_system       = col_character(),
    barber_ss         = col_integer(),
    barber_srr        = col_integer(),
    barber_poly       = col_integer(),
    barber_lek        = col_integer(),
    care_mode         = col_character(),
    dev_pc1           = col_double(),
    dev_chickpc1      = col_double(),
    flightless        = col_integer(),
    allopreen         = col_integer(),
    par_coop          = col_double(),
    age_indep         = col_double(),
    uv_cd             = col_double(),
    care_cont         = col_double(),
    dichro_sr         = col_double(),
    lislevand         = col_logical(),
    english_name      = col_character(),
    mating_system     = col_integer(),
    display_num       = col_double(),
    resource          = col_integer(),
    harem             = col_integer(),
    lek               = col_integer(),
    polyandry         = col_integer(),
    ssd_mass          = col_double(),
    ssd_tarsus        = col_double(),
    ssd_wing          = col_double(),
    ssd_tail          = col_double(),
    ssd_bill          = col_double(),
    body_mass_mean    = col_double(),
    body_size_log     = col_double()
  )
) %>%
  ## restore the focal contrast as a factor with 2D as the reference level
  mutate(dim_bin = factor(dim_bin, levels = DIM_LEVELS))

## --- reconstruct the two backbones ------------------------------------------
avo_base <- analysis_data                          # every species
merged   <- analysis_data %>% filter(lislevand)    # Lislevand subset

## --- validation anchors ------------------------------------------------------
cat("Loaded analysis dataset:\n")
cat(sprintf("  avo_base : %d species (expected 9993)\n", nrow(avo_base)))
cat(sprintf("  merged   : %d species (expected ~3509)\n", nrow(merged)))
cat("  dim_bin distribution (avo_base):\n")
print(table(avo_base$dim_bin, useNA = "ifany"))
cat("  dim_bin distribution (merged):\n")
print(table(merged$dim_bin, useNA = "ifany"))

cat("03_load_dataset.R done.\n")
