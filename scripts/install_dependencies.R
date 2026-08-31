#!/usr/bin/env Rscript

# INSTALL R PACKAGE DEPENDENCIES ----
#
# Purpose:
#   Install the CRAN and Bioconductor packages used by the analysis.
#   Packages that are already installed are left unchanged.
#
# Run from the repository folder with:
#   Rscript scripts/install_dependencies.R


# CRAN packages ----

cran_packages <- c(
  "BiocManager",
  "dplyr",
  "ggplot2",
  "ggrepel",
  "glmnet",
  "kohonen",
  "mgcv",
  "patchwork",
  "readr",
  "tidyr"
)

bioconductor_packages <- c(
  "clusterProfiler",
  "limma",
  "org.Hs.eg.db",
  "ReactomePA",
  "ropls"
)


# Install missing CRAN packages ----

installed <- rownames(installed.packages())
missing_cran <- setdiff(cran_packages, installed)
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}


# Install missing Bioconductor packages ----

installed <- rownames(installed.packages())
missing_bioconductor <- setdiff(bioconductor_packages, installed)
if (length(missing_bioconductor)) {
  BiocManager::install(missing_bioconductor, ask = FALSE, update = FALSE)
}


# Confirm that every dependency is available ----

still_missing <- setdiff(
  c(cran_packages[-1], bioconductor_packages),
  rownames(installed.packages())
)
if (length(still_missing)) {
  stop(
    "The following packages could not be installed: ",
    paste(still_missing, collapse = ", "),
    call. = FALSE
  )
}

message("All MOMI-EDLOW proteomics dependencies are installed.")
