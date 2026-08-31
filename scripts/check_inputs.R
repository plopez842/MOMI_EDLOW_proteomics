#!/usr/bin/env Rscript

# CHECK THE ANALYSIS INPUT FILES ----
#
# Purpose:
#   Confirm that both RDS files are readable and contain the samples, protein
#   columns, timepoints, and collection trimesters expected by the analysis.
#
# Run from the repository folder with:
#   Rscript scripts/check_inputs.R

project_directory <- getwd()
if (!file.exists(file.path(project_directory, "DESCRIPTION"))) {
  stop("Run this script from the MOMI_EDLOW_proteomics repository folder.", call. = FALSE)
}

source(file.path(project_directory, "helpful_functions", "data_and_setup.R"))
study_data <- load_analysis_inputs(project_directory)

cat("Analysis data:", study_data$data_path, "\n")
cat("Annotation data:", study_data$annotation_path, "\n")
cat("Samples:", nrow(study_data$samples), "\n")
cat("Protein analytes:", length(study_data$proteins), "\n")

cat("Timepoints:\n")
print(table(study_data$samples$Timepoint_v1, useNA = "ifany"))

cat("Baseline collection trimesters:\n")
baseline_rows <- study_data$samples$Timepoint_v1 == "V0"
print(table(study_data$samples$Trim_Collec[baseline_rows], useNA = "ifany"))
