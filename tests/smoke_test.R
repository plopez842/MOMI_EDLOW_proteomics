#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
source(file.path(root, "R", "io.R"))
source_project_functions(root)

inputs <- load_analysis_inputs(root)
stopifnot(nrow(inputs$samples) > 0L, length(inputs$proteins) >= 12L)

baseline <- inputs$samples[inputs$samples$Timepoint_v1 == "V0", , drop = FALSE]
trimester <- droplevels(baseline$Trim_Collec)
test_proteins <- inputs$proteins[seq_len(12L)]
x <- scale(as.matrix(baseline[, test_proteins, drop = FALSE]))

lasso <- repeated_lasso_select(x, trimester, trials = 1L, threshold = 0, seed = 1010L)
stopifnot(length(lasso$selected) > 0L)
plsda <- fit_plsda(x[, lasso$selected, drop = FALSE], trimester, n_components = 2L)
stopifnot(nrow(plsda$scores) == nrow(baseline))

limma_results <- baseline_trimester_limma(inputs$samples, test_proteins)
stopifnot(length(limma_results$tables) == 3L)

vaccine <- run_vaccine_gam_analysis(
  inputs$samples, test_proteins[1:2], dose = "V1",
  n_perm = 1L, k = 8L, cores = 1L, seed = 1010L
)
stopifnot(all(dim(vaccine$log2fc) == c(2L, 33L)))

baseline_gam <- run_baseline_gam_analysis(
  inputs$samples, test_proteins,
  n_perm = 1L, k = 8L, cores = 1L, seed = 1010L
)
som <- train_self_organizing_map(baseline_gam$fitted, clusters = 3L, seed = 1010L)
stopifnot(length(som$protein_cluster) == length(test_proteins))

acute <- run_acute_analysis(
  inputs$samples, test_proteins[1:2], dose = "V2",
  n_perm = 1L, k = 8L, cores = 1L, seed = 1010L
)
stopifnot(nrow(acute$acute) == 45L, nrow(acute$stats) == 2L)

plot <- plot_log2fc_heatmap(
  vaccine$log2fc,
  inputs$annotations,
  rownames(vaccine$log2fc),
  "Smoke test"
)
stopifnot(inherits(plot, "ggplot"))

cat("Smoke test passed.\n")
