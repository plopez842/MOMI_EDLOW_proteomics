# MOMI-EDLOW proteomics

Reproducible R code for the main figures in:

> Lopez Zapana PA, Shook LL, Joughin BA, et al. Maternal proteome profiling reveals dynamic gestational age-specific responses to de novo vaccination. *The Journal of Immunology*. 2025. doi: [10.1093/jimmun/vkaf298](https://doi.org/10.1093/jimmun/vkaf298)

The study analyzed 1,497 SOMAmer analytes targeting 1,451 unique proteins in 466 plasma samples from 278 pregnant individuals. The scripts here reconstruct the five main manuscript figures from the analysis-ready data object.

## Repository layout

| Manuscript panel | Code | Analysis |
|---|---|---|
| Figure 1 | `figures/figure_01/figure_01.R` | Baseline trimester PLSDA, longitudinal GAMs, self-organizing map, Reactome enrichment |
| Figure 2 | `figures/figure_02/figure_02.R` | Trimester differential abundance and KEGG GSEA |
| Figure 3 | `figures/figure_03/figure_03.R` | Dose 1 versus baseline GAM-permutation analysis and KEGG GSEA |
| Figure 4 | `figures/figure_04/figure_04.R` | Dose 2 versus baseline GAM-permutation analysis and KEGG GSEA |
| Figure 5 | `figures/figure_05/figure_05.R` | Acute (days 1-7) post-vaccination MoM, LRT-permutation analysis, and KEGG GSEA |

Shared functions live in `R/`. Curated pathway panels used in the publication figures live in `config/`. Each figure script saves analysis tables to `results/tables/figure_XX/` and plots to `results/figures/figure_XX/`.

## Data

The participant-level proteomics data are not stored in Git. The manuscript data are available through ImmPort study **SDY2913**. See `data/README.md` for the required local files and schema.

Place these two files under `data/processed/`:

```text
data/processed/comb_v0_vax.rds
data/processed/aptamer_annotations.rds
```

Alternatively, point `MOMI_DATA_DIR` to a directory containing those files:

```bash
export MOMI_DATA_DIR=/path/to/processed/data
```

## Setup

R 4.3 or later is recommended. Install the declared dependencies with:

```r
install.packages(c(
  "BiocManager", "circlize", "dplyr", "fs", "ggplot2", "ggrepel",
  "ggsci", "glmnet", "here", "kohonen", "patchwork", "readr",
  "scales", "tidyr"
))

BiocManager::install(c(
  "clusterProfiler", "ComplexHeatmap", "limma", "org.Hs.eg.db",
  "ReactomePA", "ropls"
))
```

Create the output directories and validate the inputs:

```bash
Rscript scripts/create_folders.R
Rscript scripts/check_inputs.R
```

## Run

Run a single figure:

```bash
Rscript figures/figure_01/figure_01.R
```

Run all main figures:

```bash
Rscript scripts/run_all.R
```

The publication analyses use 100 repeated LASSO selections and 1,000 permutations. These are the defaults and can be changed for a quick test:

```bash
MOMI_LASSO_TRIALS=2 MOMI_N_PERM=2 MOMI_CORES=1 Rscript tests/smoke_test.R
```

The full GAM-permutation workflows are computationally intensive. Intermediate results are checkpointed in `results/cache/`, allowing interrupted runs to resume without repeating completed proteins.

## Reproducibility notes

- A fixed random seed is used throughout.
- Gestational trimesters are defined as first: `<14` weeks, second: `14-<28` weeks, and third: `>=28` weeks.
- Main-figure vaccine analyses compare V0 with samples collected more than 7 days after Dose 1 (V1) or Dose 2 (V2). Acute analyses use days 1-7 after vaccination.
- Differential abundance uses `limma`, with FDR `<0.05` and absolute log2 fold change `>log2(1.5)`.
- GAMs use REML-fitted cubic regression splines with `k=8`. Vaccine significance is evaluated against label-permuted empirical null distributions with Benjamini-Hochberg correction.
- Pathways are significant at FDR `<0.25`, consistent with the manuscript.

## Citation

Please cite the article above and this repository when reusing the workflow.
