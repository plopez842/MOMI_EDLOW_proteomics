#!/usr/bin/env Rscript

# CREATE OUTPUT FOLDERS ----
#
# Purpose:
#   Create the results/figures, results/tables, and results/cache folders used
#   by Figures 1-5. Figure scripts also create these folders automatically.
#
# Run from the repository folder with:
#   Rscript scripts/create_folders.R

project_directory <- getwd()
if (!file.exists(file.path(project_directory, "DESCRIPTION"))) {
  stop("Run this script from the MOMI_EDLOW_proteomics repository folder.", call. = FALSE)
}

source(file.path(project_directory, "helpful_functions", "data_and_setup.R"))
create_output_directories(project_directory)

message("Output folders are ready under: ", file.path(project_directory, "results"))
