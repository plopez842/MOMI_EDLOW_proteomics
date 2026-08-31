[![Article DOI](https://img.shields.io/badge/DOI-10.1093%2Fjimmun%2Fvkaf298-2F6F9F)](https://doi.org/10.1093/jimmun/vkaf298)

# MOMI-EDLOW Proteomics

R code for the five main figures in:

> Lopez Zapana PA, Shook LL, Joughin BA, et al. Maternal proteome profiling reveals dynamic gestational age-specific responses to de novo vaccination. *The Journal of Immunology*. 2025. [doi:10.1093/jimmun/vkaf298](https://doi.org/10.1093/jimmun/vkaf298)

The study includes 1,497 SOMAmer analytes measured in 466 plasma samples from 278 pregnant individuals.

## Where is the code?

Each manuscript figure has one top-level folder. Open the folder, then run the R script inside it. If a figure uses a selected list of pathways, that list is in the same folder as `pathways_shown.csv`.

| Figure | Script | Analysis |
|---|---|---|
| 1 | `Figure_1/figure_1.R` | Baseline trimester PLSDA, gestational GAMs, SOM clustering, Reactome enrichment |
| 2 | `Figure_2/figure_2.R` | Trimester differential abundance and KEGG GSEA |
| 3 | `Figure_3/figure_3.R` | Dose 1 response and KEGG GSEA |
| 4 | `Figure_4/figure_4.R` | Dose 2 response and KEGG GSEA |
| 5 | `Figure_5/figure_5.R` | Acute days 1-7 response and KEGG GSEA |

```text
MOMI_EDLOW_proteomics/
├── Figure_1/              Figure 1 script and pathways shown
├── Figure_2/              Figure 2 script and pathways shown
├── Figure_3/              Figure 3 script and pathways shown
├── Figure_4/              Figure 4 script and pathways shown
├── Figure_5/              Figure 5 script and pathways shown
├── helpful_functions/     Clearly named shared analysis and plotting functions
├── data/                  Data instructions; participant data are not in Git
├── results/               Generated plots, tables, and saved intermediate results
├── scripts/               Setup, data check, and run-all scripts
└── tests/                 Small end-to-end test
```

## Data

Study data are available from ImmPort under accession **SDY2913**. The participant-level data are not included in Git.

The analysis needs these two R objects:

```text
data/processed/comb_v0_vax.rds
data/processed/aptamer_annotations.rds
```

See [`data/README.md`](data/README.md) for the required columns. The files may also live elsewhere:

```bash
export MOMI_DATA_DIR=/path/to/the/folder/containing/the/rds/files
```

## Run the analysis

Use R 4.3 or later. From the repository folder:

```bash
# Install packages once
Rscript scripts/install_dependencies.R

# Check that the data are readable
Rscript scripts/check_inputs.R

# Run one figure
Rscript Figure_1/figure_1.R

# Or run all five figures
Rscript scripts/run_all.R
```

The full analysis uses 100 LASSO repetitions and 1,000 permutations and may take a long time. Completed calculations are saved in `results/cache/`, so the scripts can resume after an interruption.

For a quick test of the pipeline:

```bash
MOMI_LASSO_TRIALS=1 MOMI_N_PERM=1 MOMI_CORES=1 Rscript tests/smoke_test.R
```

## Outputs

- Plots: `results/figures/figure_01/` through `figure_05/`
- Tables: `results/tables/figure_01/` through `figure_05/`
- Saved intermediate calculations: `results/cache/`

The output folders are created automatically. `Rscript scripts/create_folders.R` can also create them before an analysis starts.

## How the scripts are organized

Every figure script begins with its purpose, manuscript panels, required inputs, outputs, and exact run command. Numbered section headings then follow the analysis in order. Shared functions are grouped by scientific purpose:

- `data_and_setup.R`: read data, check columns, define settings, and create folders
- `statistical_models.R`: LASSO, PLSDA, GAMs, permutations, SOM, and acute models
- `pathway_analysis.R`: KEGG GSEA and Reactome enrichment
- `plotting_functions.R`: manuscript plots and colors
- `vaccine_analysis.R`: trimester, vaccine, and acute comparison helpers

## Citation

Please cite the article above when using this workflow. GitHub can export the complete citation from [`CITATION.cff`](CITATION.cff).

## License and contact

No open-source license has been selected yet; see [`LICENSE`](LICENSE). For questions, please open a GitHub Issue.
