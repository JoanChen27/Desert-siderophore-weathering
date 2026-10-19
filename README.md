# Desert-siderophore-weathering

Analysis code for the manuscript **"Microbial siderophore activity can initiate weathering of barren desert sand"** (first author: Qiong Chen).

The study asks whether microbes can directly initiate weathering of barren desert sand, using a 17-year Tengger Desert restoration chronosequence, screening of 805 cultivated isolates, and the *bibC*-dependent desferrioxamine E (DFOE) siderophore pathway of strain F8 (*Pantoea alhagi* / *Mixta alhagi*) with knockout, complementation and overexpression genetics plus 84-day field plots.

## Repository contents

- `Desert.Rmd` — main analysis notebook (R Markdown), 8 parts:
  1. Field chronosequence and sand-matrix properties (physicochemical properties, D50 and hydraulic conductivity, mineral composition and PCA)
  2. Microbial diversity and community succession (beta-diversity, alpha diversity, phylogenetic tree, fuzzy C-means clusters, GLM, permutation tests, LASSO, dbRDA, gradient boosting)
  3. Functional and biogeochemical reprogramming (KEGG pathways, functional maptree, mediation analysis)
  4. Cultivated-isolate screening and field representation (isolate taxonomy, isolate–MAG overlap, sterile-sand screening)
  5. Laboratory validation of sand weathering (F8 incubation experiments)
  6. F8 metabolism and siderophore mechanism (metabolomics, DFOE identification and quantification, *bibC* mutant and overexpression strains)
  7. Siderophore-associated patterns across the field chronosequence
- `Global_weathering_potential_model.R` — global weathering potential model: function-only microbial suitability with V18 fixed-domain substrate scaling and AI+PET climate. Reads curated input tables (sand-matrix, climate, SoilGrids covariates, microbe tolerance table) and writes model outputs, tables and figures.

## Requirements

R (≥ 4.0 recommended) with the following packages (see `library()` calls at the top of each file):

- `Desert.Rmd`: tidyverse, vegan, ggplot2, ggtext, patchwork, ggrepel, broom, phyloseq, ggtree, lme4, glmnet, caret, gbm, pheatmap, and others
- Global model script: tidyverse, readxl, janitor, mgcv

## Data

Input data files are **not** included in this repository. Both scripts reference local file paths (e.g. `~/Desktop/replot/...`, `~/Documents/...`); update the path variables at the top of each file to point to your own data locations before running.

Sequencing data are deposited at ENA. See the manuscript for accession details and data availability.

## Running the code

- `Desert.Rmd`: open in RStudio and knit, or run chunks interactively in order. Sections are numbered and can be run independently once upstream objects exist.
- Global model script: set `input_dir` / `output_dir` at the top of the file, then `source()` the whole script or run sections in order.

## Status

Manuscript under appeal. This repository is private until publication; code will be made public on acceptance.
