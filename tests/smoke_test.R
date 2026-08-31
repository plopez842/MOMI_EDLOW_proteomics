#!/usr/bin/env Rscript

# REDUCED END-TO-END PIPELINE TEST ----
#
# Purpose:
#   Check the major analysis functions on 12 proteins and one permutation.
#   This confirms that the data format and installed packages work together;
#   it does not reproduce the publication results.
#
# Run from the repository folder with:
#   Rscript tests/smoke_test.R

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the MOMI_EDLOW_proteomics repository folder.", call. = FALSE)
}

source("helpful_functions/project_setup.R")
source_analysis_functions()

study_data <- load_analysis_inputs()
stopifnot(nrow(study_data$samples) > 0L, length(study_data$proteins) >= 12L)

# Use only 12 proteins so this check finishes quickly.
test_proteins <- study_data$proteins[seq_len(12L)]
baseline_samples <- study_data$samples[
  study_data$samples$Timepoint_v1 == "V0",
  ,
  drop = FALSE
]


# LASSO and PLSDA ----

collection_trimester <- droplevels(baseline_samples$Trim_Collec)
scaled_test_proteins <- scale(
  as.matrix(baseline_samples[, test_proteins, drop = FALSE])
)

lasso_test <- repeated_lasso_select(
  scaled_test_proteins,
  collection_trimester,
  trials = 1L,
  threshold = 0,
  seed = 1010L
)
stopifnot(length(lasso_test$selected) > 0L)

plsda_test <- fit_plsda(
  scaled_test_proteins[, lasso_test$selected, drop = FALSE],
  collection_trimester,
  n_components = 2L
)
stopifnot(nrow(plsda_test$scores) == nrow(baseline_samples))


# Trimester differential abundance ----

limma_test <- baseline_trimester_limma(study_data$samples, test_proteins)
stopifnot(length(limma_test$tables) == 3L)


# Dose 1 GAM permutation analysis ----

dose_1_test <- run_vaccine_gam_analysis(
  study_data$samples,
  test_proteins[1:2],
  dose = "V1",
  n_perm = 1L,
  k = 8L,
  cores = 1L,
  seed = 1010L
)
stopifnot(all(dim(dose_1_test$log2fc) == c(2L, 33L)))


# Baseline GAM and self-organizing map ----

baseline_gam_test <- run_baseline_gam_analysis(
  study_data$samples,
  test_proteins,
  n_perm = 1L,
  k = 8L,
  cores = 1L,
  seed = 1010L
)
som_test <- train_self_organizing_map(
  baseline_gam_test$fitted,
  clusters = 3L,
  seed = 1010L
)
stopifnot(length(som_test$protein_cluster) == length(test_proteins))


# Acute Dose 2 analysis ----

acute_dose_2_test <- run_acute_analysis(
  study_data$samples,
  test_proteins[1:2],
  dose = "V2",
  n_perm = 1L,
  k = 8L,
  cores = 1L,
  seed = 1010L
)
stopifnot(
  nrow(acute_dose_2_test$acute) == 45L,
  nrow(acute_dose_2_test$stats) == 2L
)


# Plotting ----

heatmap_test <- plot_log2fc_heatmap(
  dose_1_test$log2fc,
  study_data$annotations,
  rownames(dose_1_test$log2fc),
  "Smoke test"
)
stopifnot(inherits(heatmap_test, "ggplot"))

cat("Smoke test passed.\n")
