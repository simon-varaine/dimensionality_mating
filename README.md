# Arena dimensionality and the form of sexual selection in birds

Reproducible code for the phylogenetic comparative analyses in the manuscript
*"Three-dimensional lifestyles are associated with reduced contest and increased
display in birds"*.

## Two ways to use this repository

**Reproduce the analyses (anyone).** You need only the shared analysis dataset
(`data-derived/analysis_data.csv`, included here) and the phylogeny file
(downloaded separately, see below). Then, from the project root:

```r
source("run_all.R")      # in an R session
# or:  Rscript run_all.R
```

**Rebuild the analysis dataset from raw sources (authors).** This requires the
original third-party datasets in `data/` (not shared here — see the manifest).

```r
source("build_dataset.R")   # writes data-derived/analysis_data.csv
```

A full `run_all.R` fits every model across 100 trees (~30-40 min on a laptop) and
caches the results in `output/model_results.rds`. To redraw figures and tables
from the cache without refitting (a few seconds):

```r
REFIT <- FALSE; source("run_all.R")
```

## Outputs (paper-ready)

| File | Content | In the paper |
|------|---------|--------------|
| `output/tables/table1_descriptives.tex` | every response by foraging lifestyle (% or mean (SD), with *n*) | Table 1 |
| `output/figures/fig1_main_effects.pdf` | main (bivariate) effects of 3D, 90% + 95% CIs | Figure 1 |
| `output/tables/tableS1_main_effects.tex` | main effects, full statistics (N 2D/3D, estimate, CI, % trees, one- and two-sided p) | Table S1 |
| `output/tables/tableS2_controls.tex` | main effects with controls / on subsets (Rensch, display excl. aerial / + HWI, Passeriformes) | Table S2 |
| `output/figures/figS1_passeriformes.pdf` | main effects within Passeriformes | Figure S1 |
| `output/tables/tableS3_care_interaction.tex` + `output/figures/figS2_care_interaction.pdf` | 3D × continuous parental-care interaction | Table S3, Figure S2 |

Each `.tex` is a complete `table`/`figure` environment (caption and label
included; tables need only `booktabs`) to paste or `\input{}` into the
manuscript; each table also has a `.csv` twin with full-precision numbers, and
each figure a `.png` preview.

**Statistics.** Each model is fitted on each of 100 trees and pooled with
Rubin's rules (within- + between-tree variance). Binary responses: penalised
phylogenetic logistic regression (`phyloglm`, `logistic_MPLE`); continuous
responses: PGLS with Pagel's λ (`phylolm`). Predictions are directional, so
one-sided p-values are reported alongside two-sided ones, and figures show both
90% (≈ one-sided 5%) and 95% CIs.

## Repository layout

```
project/
├── data/                 # raw third-party inputs (git-ignored; NOT shared)
├── data-derived/
│   └── analysis_data.csv # shared analysis-ready dataset (one row per species)
├── output/               # generated figures + tables (git-ignored)
│   ├── figures/          #   Figures 1, S1, S2 (.pdf, .png, .tex)
│   ├── tables/           #   Tables 1, S1, S2, S3 (.tex, .csv)
│   └── model_results.rds #   cache of all fitted models
├── R/
│   ├── 00_setup.R        # libraries, paths, model engine, helpers, constants
│   ├── 01_load_raw.R     # load + clean the raw sources             [authors]
│   ├── 02_build_dataset.R# assemble + write analysis_data.csv       [authors]
│   ├── 03_load_dataset.R # read analysis_data.csv -> merged, avo_base
│   ├── 04_trees.R        # master tree -> 100 sampled trees
│   ├── 05_models.R       # all models (main effects, controls, care interaction)
│   ├── 06_figures.R      # Figures 1, S1, S2
│   ├── 07_tables.R       # Tables 1, S1, S2, S3
│   └── 99_exploratory.R  # analyses NOT in the paper (not run by run_all.R)
├── build_dataset.R       # entry point for authors (raw -> derived)
├── run_all.R             # entry point for anyone (derived -> results)
└── README.md
```

## Phylogeny (download required)

The tree file is large and is **not** stored in this repository. Download the
Hackett backbone "All species" set of 1000 trees from
<https://birdtree.org>, and place it in `data/` as:

```
data/AllBirdsHackett1.tre
```

`04_trees.R` reads it and samples 100 trees.

## The shared analysis dataset

`data-derived/analysis_data.csv` has one row per species and contains every
variable used in the analyses, so the results can be reproduced without the raw
sources. Key columns:

| Column | Description |
|--------|-------------|
| `key`, `tip_label`, `scientific` | join key, BirdTree tip label, binomial |
| `order`, `family`, `primary_lifestyle` | taxonomy and AVONET foraging lifestyle |
| `dim_bin` | dimensionality contrast, `"2D"` (Terrestrial) vs `"3D"` (Insessorial+Aerial) |
| `lislevand` | logical flag: species present in the Lislevand join |
| `mating_system` | Lislevand 1–5 social mating-system score |
| `harem`, `lek`, `polyandry` | binary contrasts derived from `mating_system` |
| `marc_poly`, `marc_lek` | Marcondes & Douvas resource-defense-polygamy / lekking contrasts |
| `ssd_mass`, `ssd_tarsus`, `ssd_wing`, `ssd_tail`, `ssd_bill` | log(male/female) size dimorphism |
| `display_num` | Lislevand display-agility score (1–5) |
| `resource` | Lislevand between-sex resource-sharing index (0–2) |
| `dichromatism`, `male_plumage`, `female_plumage` | Dale plumage scores (male − female, and each sex) |
| `spur_hi` | bony-spur presence (high confidence) |
| `hwi` | hand-wing index |
| `body_mass_mean`, `body_size_log` | mean body mass and its log (for the Rensch control) |
| `marc_system` | Marcondes & Douvas category (`M`/`P`/`L`); `03_load_dataset.R` derives the monogamy-referenced contrasts `rdp`, `lek_m`, `lek_vs_rdp` |
| `barber_ss`, `barber_srr` | Barber sexual-selection intensity (0–4) and sex-role reversal |
| `barber_poly`, `barber_lek` | strong polygamy (2–3) / lek (4) vs monogamy (0–1), from `barber_ss` |
| `dichro_sr`, `care_cont` | broad plumage dimorphism (all orders) and relative care investment of the sexes (> 0 = female-biased), sex-role ecology dataset |
| `uv_cd` | UV-inclusive colour discriminability between the sexes |
| `allopreen` | allopreening between pair members (0/1), Kenny et al. 2017 |
| `care_mode`, `dev_pc1`, `dev_chickpc1`, `flightless`, `par_coop`, `age_indep` | used only in `R/99_exploratory.R` |

In `03_load_dataset.R` the dataset is split into the two analysis backbones used
throughout: `avo_base` (all species) and `merged` (the Lislevand subset,
`lislevand == TRUE`).

## Data manifest (raw sources, to rebuild the dataset)

Download each source and drop the file into `data/` **under the name it ships
with** — no renaming needed. These files are not redistributed here. Where a
source is delivered as an archive (e.g. AVONET's `ELEData/`, a Dryad `.zip`),
extract the single file named below into `data/`.

| File name in `data/` | Source / notes | Reference |
|---|---|---|
| `avian_ssd_jan07.txt` | Ecological Archives E088-096 | Lislevand et al. 2007, *Ecology* |
| `AVONET3_BirdTree.xlsx` (sheet `AVONET3_BirdTree`) | from the AVONET supplement (`ELEData/TraitData/`) | Tobias et al. 2022, *Ecol. Lett.* |
| `BLIOCPhyloMasterTax.csv` | BirdTree taxonomy ↔ tip labels | Jetz et al. 2012, *Nature* |
| `plumage_scores.csv` | male/female plumage scores | Dale et al. 2015, *Nature* |
| `Mating_systems_master_datasheet_10nov2023.xlsx` (sheet `Species_data`) | remove any ` (n)` download suffix so the name matches | Marcondes & Douvas 2024, *Evolution* |
| `species_spur_data.csv` | Dryad | Menezes & Palaoro 2022, *Ecol. Lett.* |
| `S1_Data.xlsx` (sheet `Data1 (BirdTree)`) | sexual-selection intensity (0–4) + sex-role reversal | Barber et al. 2024, *PLoS Biology* |
| `SexroleEcologyFinal.xlsx` (sheet `Data`) | continuous care (`Care`) + broad plumage dimorphism (`Dichro`) | sex-role ecology dataset (reference to complete) |
| `MergedCDLatSet.xlsx` | UV-inclusive colour discriminability | Villar et al. 2025, *J. Zool.* |
| `arx078_suppl_kenny_esm_tables1.xlsx` (sheet `Kenny_TableS1_ESM`) | allopreening, parental cooperation | Kenny et al. 2017, *Behav. Ecol.* |
| `GeneralDataFrame.csv` | care mode (exploratory only) | "Who cares?" dataset, Dryad |
| `evo14365-sup-0009-datasets2.xlsx` | developmental mode (exploratory only) | Cooney et al. 2021, *Evolution* |
| `DataFileS2.xlsx` | volancy (exploratory only) | Sayol et al. 2020 |
| `AllBirdsHackett1.tre` | 1000 Hackett "All species" trees (large; download from birdtree.org) | Jetz et al. 2012 |

If a future release of any source changes its file name, either keep the old
name or update the corresponding `path_*` line in `R/00_setup.R`.

## Data sources & licences

`data-derived/analysis_data.csv` contains variables **derived from** the sources
above. Before redistributing it, confirm that each source's licence permits
reuse of derived values, and cite all sources. AVONET is released under CC-BY;
Dryad deposits are typically CC0; journal/archive supplements generally permit
reuse with citation — but verify per source. If any source disallows
redistribution of derived data, share the code plus a species list instead and
have users obtain that source themselves.

## Column expectations (raw sources)

`01_load_raw.R` reads specific columns. If a source is a different release with
renamed columns, adjust the `transmute()` blocks there. Expected fields:

- **Lislevand**: `species_name`, `English_name`, `Mating_system`, `Display`,
  `Resource`, and `M_*/F_*` mass/wing/tarsus/tail/bill (NA coded as `-999`).
- **AVONET**: species, family, order, `Primary.Lifestyle`, habitat, mass,
  hand-wing index (detected by pattern via `get_col()`).
- **BLIOCPhyloMasterTax**: `Scientific`, `TipLabel`.
- **Dale**: `TipLabel`, `male_plumage_score`, `female_plumage_score`.
- **Marcondes**: `Species`, `Mating_system` coded `M`/`P`/`L` (`U`/`PU` = NA).
- **Menezes & Palaoro**: `Scientific`, `AnySpurPresence.HighConf`,
  `HandWingIndex`, `logBodyMass`.

## Scope

`run_all.R` reproduces exactly the analyses reported in the manuscript: the
main (bivariate) effects of 3D (main text), the same effects with controls or
on subsets, and the 3D × continuous parental-care interaction (appendix).
Exploratory analyses that are not in the paper are kept, but not run, in
`R/99_exploratory.R`: resource-sharing interactions, the 3D × sexual-selection
intensity interaction, binary care mode and developmental-mode interactions,
flightlessness, and the hand-wing index as an alternative dimensionality proxy.
Earlier branches (non-phylogenetic sanity checks, fuzzy-matching diagnostics,
an extra-pair-paternity module, an aquatic-vs-terrestrial regression) were
removed.
