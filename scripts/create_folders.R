#!/usr/bin/env Rscript

# CREATE OUTPUT FOLDERS ----
#
# Purpose:
#   Create the results/figures, results/tables, and results/cache folders used
#   by Figures 1-5. Figure scripts also create these folders automatically.
#
# Run from the repository folder with:
#   Rscript scripts/create_folders.R

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the MOMI_EDLOW_proteomics repository folder.", call. = FALSE)
}

source("helpful_functions/project_setup.R")
create_output_directories()

message("Output folders are ready under: results/")
