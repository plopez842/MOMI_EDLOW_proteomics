#!/usr/bin/env Rscript

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

installed <- rownames(installed.packages())
missing_cran <- setdiff(cran_packages, installed)
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

installed <- rownames(installed.packages())
missing_bioconductor <- setdiff(bioconductor_packages, installed)
if (length(missing_bioconductor)) {
  BiocManager::install(missing_bioconductor, ask = FALSE, update = FALSE)
}

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
