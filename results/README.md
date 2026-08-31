# Generated results

The analysis scripts create this directory tree automatically:

```text
results/
├── figures/
│   ├── figure_01/
│   ├── figure_02/
│   ├── figure_03/
│   ├── figure_04/
│   └── figure_05/
├── tables/
│   ├── figure_01/
│   ├── figure_02/
│   ├── figure_03/
│   ├── figure_04/
│   └── figure_05/
└── cache/
    ├── figure_01/
    ├── figure_02/
    ├── figure_03/
    ├── figure_04/
    └── figure_05/
```

- `figures/` contains manuscript-panel PDFs.
- `tables/` contains model statistics, protein lists, and enrichment results as TSV files.
- `cache/` contains restartable intermediate calculations for computationally intensive permutation analyses.

Generated outputs are excluded from Git. Run `Rscript scripts/create_folders.R` from the repository root to create the empty structure, or run any figure script and its required output directories will be created automatically.
