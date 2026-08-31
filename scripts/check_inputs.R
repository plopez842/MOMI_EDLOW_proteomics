#!/usr/bin/env Rscript

# CHECK THE ANALYSIS INPUT FILES ----
#
# Purpose:
#   Confirm that both RDS files are readable and contain the samples, protein
#   columns, timepoints, and collection trimesters expected by the analysis.
#
# Run from the repository folder with:
#   Rscript scripts/check_inputs.R

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the MOMI_EDLOW_proteomics repository folder.", call. = FALSE)
}

source("helpful_functions/project_setup.R")
study_data <- load_analysis_inputs()

cat("Analysis data:", study_data$data_path, "\n")
cat("Annotation data:", study_data$annotation_path, "\n")
cat("Samples:", nrow(study_data$samples), "\n")
cat("Protein analytes:", length(study_data$proteins), "\n")

cat("Timepoints:\n")
print(table(study_data$samples$Timepoint_v1, useNA = "ifany"))

cat("Baseline collection trimesters:\n")
baseline_rows <- study_data$samples$Timepoint_v1 == "V0"
print(table(study_data$samples$Trim_Collec[baseline_rows], useNA = "ifany"))
