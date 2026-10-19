# Desert-siderophore-weathering

Analysis code for the manuscript **"Microbial siderophore activity can initiate weathering of barren desert sand"** (first author: Qiong Chen).

The study asks whether microbes can directly initiate weathering of barren desert sand, using a 17-year Tengger Desert restoration chronosequence, screening of 805 cultivated isolates, and the *bibC*-dependent desferrioxamine E (DFOE) siderophore pathway of strain F8 (*Pantoea alhagi* / *Mixta alhagi*) with knockout, complementation and overexpression genetics plus 84-day field plots.

## Repository contents

- `Global_weathering_potential_model.R` — global weathering potential model: function-only microbial suitability with V18 fixed-domain substrate scaling and AI+PET climate. Reads curated input tables (sand-matrix, climate, SoilGrids covariates, microbe tolerance table) and writes model outputs, tables and figures.
- `code/` — main analyses, split from the original `Desert.Rmd` into one file per part (parts numbered continuously 01–07):

| File | Part | Content |
|---|---|---|
| `00_setup.R` | — | Shared package loading (source this first) |
| `01_chronosequence_sand_properties.Rmd` | 1 | Field chronosequence and sand-matrix properties: physicochemical properties, D50 and hydraulic conductivity, mineral composition and PCA |
| `02a_community_diversity.Rmd` | 2a | Microbial diversity: beta-diversity, alpha diversity, phylogenetic tree, fuzzy C-means clusters |
| `02b_community_glm_lasso.Rmd` | 2b | Community modeling I: GLM regression, graph-based permutation test, AEM, LASSO prediction |
| `02c_community_dbrda_gbm.Rmd` | 2c | Community modeling II: dbRDA and variation partitioning, gradient boosting |
| `03_functional_reprogramming.Rmd` | 3 | Functional and biogeochemical reprogramming: KEGG pathways, functional maptree, mediation analysis |
| `04_isolate_screening.Rmd` | 4 | Cultivated-isolate screening and field representation: isolate taxonomy, isolate–MAG overlap, sterile-sand screening |
| `05_lab_weathering_validation.Rmd` | 5 | Laboratory validation of sand weathering: F8 incubation experiments |
| `06_F8_metabolism_siderophore.Rmd` | 6 | F8 metabolism and siderophore mechanism: metabolomics, DFOE identification and quantification, *bibC* mutant and overexpression strains |
| `07_siderophore_field_patterns.Rmd` | 7 | Siderophore-associated patterns across the field chronosequence |

Each part file starts with a setup chunk that sources `00_setup.R`. Run the parts in numbered order in the same R session (later parts may use objects created earlier).

## Requirements

R (≥ 4.0 recommended). All required packages are listed in `code/00_setup.R` (deduplicated from the original draft):

tidyverse, vegan, ggplot2, ggtext, patchwork, ggrepel, broom, phyloseq, ggtree, lme4, glmnet, caret, gbm, pheatmap, readxl, janitor, mgcv, and others — see `00_setup.R` for the full list.

## Data

Input data files are **not** included in this repository. The scripts reference local file paths (e.g. `~/Desktop/replot/...`, `~/Documents/...`); update the path variables at the top of each file to point to your own data locations before running.

Sequencing data are deposited at ENA. See the manuscript for accession details and data availability.

## Status

Manuscript under appeal. This repository is private until publication; code will be made public on acceptance.
