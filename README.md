# Desert-siderophore-weathering

Analysis code for the manuscript **"Microbial siderophore activity can initiate weathering of barren desert sand"** (first author: Qiong Chen).

The study asks whether microbes can directly initiate weathering of barren desert sand, using a 17-year Tengger Desert restoration chronosequence, screening of 805 cultivated isolates, and the *bibC*-dependent desferrioxamine E (DFOE) siderophore pathway of strain F8 (*Pantoea alhagi* / *Mixta alhagi*) with knockout, complementation and overexpression genetics plus 84-day field plots.

## Repository contents

All analysis code is in `code/`, numbered in run order (later files may use objects created earlier, so run them in sequence in the same R session):

| File | Content |
|---|---|
| `00_setup.R` | Shared package loading — source this first |
| `01_chronosequence_sand_properties.Rmd` | Field chronosequence and sand-matrix properties: physicochemical properties, D50 and hydraulic conductivity, mineral composition and PCA |
| `02a_community_diversity.Rmd` | Microbial diversity: beta-diversity, alpha diversity, phylogenetic tree, fuzzy C-means clusters |
| `02b_community_glm_lasso.Rmd` | Community modeling I: GLM regression, graph-based permutation test, AEM, LASSO prediction |
| `02c_community_dbrda_gbm.Rmd` | Community modeling II: dbRDA and variation partitioning, gradient boosting |
| `03_functional_reprogramming.Rmd` | Functional and biogeochemical reprogramming: KEGG pathways, functional maptree, mediation analysis |
| `04_isolate_screening.Rmd` | Cultivated-isolate screening and field representation: isolate taxonomy, isolate–MAG overlap, sterile-sand screening |
| `05_lab_weathering_validation.Rmd` | Laboratory validation of sand weathering: F8 incubation experiments |
| `06_F8_metabolism_siderophore.Rmd` | F8 metabolism and siderophore mechanism: metabolomics, DFOE identification and quantification, *bibC* mutant and overexpression strains |
| `07_siderophore_field_patterns.Rmd` | Siderophore-associated patterns across the field chronosequence |
| `08_global_weathering_potential_model.R` | Global comparison assessing environmental suitability for F8-associated weathering across sandy drylands: a relative weathering-potential model parameterized by sand-matrix reactivity, climate suitability, and microbial suitability (the experimentally measured functional response of F8) |

## Requirements

R (≥ 4.0 recommended). All required packages are listed in `code/00_setup.R`:

tidyverse, vegan, ggplot2, ggtext, patchwork, ggrepel, broom, phyloseq, ggtree, lme4, glmnet, caret, gbm, pheatmap, readxl, janitor, mgcv, and others — see `00_setup.R` for the full list.

## Data

The underlying datasets are publicly available as described in the manuscript:

- Raw sequencing reads, metagenome-assembled genomes (MAGs), isolate 16S sequences and genomes: ENA under [PRJEB60065](https://www.ebi.ac.uk/ena/browser/view/PRJEB60065), [PRJEB128436](https://www.ebi.ac.uk/ena/browser/view/PRJEB128436) and [PRJEB128464](https://www.ebi.ac.uk/ena/browser/view/PRJEB128464)
- LC–MS metabolomics data: MetaboLights under [MTBLS15530](https://www.ebi.ac.uk/metabolights/MTBLS15530)

## Status

Manuscript under review. This repository is private until publication; code will be made public on acceptance.
