#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
source(file.path(root, "R", "io.R"))
inputs <- load_analysis_inputs(root)

cat("Analysis data:", inputs$data_path, "\n")
cat("Annotation data:", inputs$annotation_path, "\n")
cat("Samples:", nrow(inputs$samples), "\n")
cat("Protein analytes:", length(inputs$proteins), "\n")
cat("Timepoints:\n")
print(table(inputs$samples$Timepoint_v1, useNA = "ifany"))
cat("Baseline collection trimesters:\n")
print(table(inputs$samples$Trim_Collec[inputs$samples$Timepoint_v1 == "V0"], useNA = "ifany"))
