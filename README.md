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

Outputs (figures + tables) are written to `output/`. A full `run_all.R` fits
every model across 100 trees and takes on the order of a few to a few tens of
minutes depending on the machine.

## Repository layout

```
project/
├── data/                 # raw third-party inputs (git-ignored; NOT shared)
├── data-derived/
│   └── analysis_data.csv # shared analysis-ready dataset (one row per species)
├── output/               # generated figures + tables (git-ignored)
├── R/
│   ├── 00_setup.R        # libraries, paths, helpers, constants
│   ├── 01_load_raw.R     # load + clean the six raw sources        [authors]
│   ├── 02_build_dataset.R# assemble + write analysis_data.csv       [authors]
│   ├── 03_load_dataset.R # read analysis_data.csv -> merged, avo_base
│   ├── 04_trees.R        # master tree -> 100 sampled trees
│   ├── 05_models.R       # all PGLS models
│   ├── 06_figures.R      # Figure 1, Figure 2, Figure 3
│   └── 07_tables.R       # Table 1, A1, A2, A3, A4
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

In `03_load_dataset.R` the dataset is split into the two analysis backbones used
throughout: `avo_base` (all species) and `merged` (the Lislevand subset,
`lislevand == TRUE`).

## Data manifest (raw sources, to rebuild the dataset)

Place each file directly in `data/` under the cleaned name below. These files
are **not** redistributed here; download them from the original publications.

| Cleaned name in `data/` | Original file | Reference |
|---|---|---|
| `lislevand_2007_avian_ssd.txt` | `avian_ssd_jan07.txt` (Ecological Archives E088-096) | Lislevand et al. 2007, *Ecology* |
| `avonet3_birdtree.xlsx` (sheet `AVONET3_BirdTree`) | AVONET supplementary dataset | Tobias et al. 2022, *Ecol. Lett.* |
| `birdtree_taxonomy.csv` | `BLIOCPhyloMasterTax.csv` | Jetz et al. 2012, *Nature* |
| `dale_2015_plumage_scores.csv` | plumage scores | Dale et al. 2015, *Nature* |
| `marcondes_douvas_2024_mating_systems.xlsx` (sheet `Species_data`) | mating-systems datasheet | Marcondes & Douvas 2024, *Evolution* |
| `menezes_palaoro_2022_spurs.csv` | `species_spur_data.csv` (Dryad) | Menezes & Palaoro 2022, *Ecol. Lett.* |
| `AllBirdsHackett1.tre` | 1000 Hackett all-species trees | Jetz et al. 2012 / birdtree.org |

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

This code reproduces exactly the analyses reported in the manuscript. Earlier
exploratory branches that were not part of the paper have been removed for
clarity (non-phylogenetic sanity checks, fuzzy-matching diagnostics, an
extra-pair-paternity module, an aquatic-vs-terrestrial regression, spur
demarcation regressions prevented by complete separation, and a superseded
summary figure).
