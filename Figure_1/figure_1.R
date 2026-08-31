#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
source(file.path(root, "R", "io.R"))
source_project_functions(root)
create_output_directories(root)

parameters <- analysis_parameters()
set.seed(parameters$seed)
inputs <- load_analysis_inputs(root)
paths <- figure_paths(1, root)
dir.create(paths$cache, recursive = TRUE, showWarnings = FALSE)

baseline <- inputs$samples[inputs$samples$Timepoint_v1 == "V0", , drop = FALSE]
trimester <- factor(
  baseline$Trim_Collec,
  levels = c("Trim_1st", "Trim_2nd", "Trim_3rd"),
  labels = c("1st", "2nd", "3rd")
)
x <- scale(as.matrix(baseline[, inputs$proteins, drop = FALSE]))
x <- x[, colSums(!is.finite(x)) == 0, drop = FALSE]

lasso_cache <- file.path(paths$cache, paste0("lasso_", parameters$lasso_trials, "_trials.rds"))
if (file.exists(lasso_cache)) {
  lasso <- readRDS(lasso_cache)
} else {
  lasso <- repeated_lasso_select(
    x, trimester,
    trials = parameters$lasso_trials,
    threshold = 0.8,
    seed = parameters$seed
  )
  saveRDS(lasso, lasso_cache)
}

x_selected <- x[, lasso$selected, drop = FALSE]
plsda <- fit_plsda(x_selected, trimester, n_components = 2L)
vip <- sort(plsda$vip, decreasing = TRUE)

save_table(
  data.frame(protein = names(lasso$frequency), selection_frequency = as.numeric(lasso$frequency)),
  file.path(paths$tables, "lasso_selection_frequency.csv")
)
save_table(
  data.frame(protein = names(vip), gene = gene_symbols(names(vip), inputs$annotations), VIP = as.numeric(vip)),
  file.path(paths$tables, "plsda_vip_scores.csv")
)
loading_table <- data.frame(protein = rownames(plsda$loadings), plsda$loadings, check.names = FALSE)
loading_table$gene <- gene_symbols(loading_table$protein, inputs$annotations)
save_table(loading_table, file.path(paths$tables, "plsda_loadings.csv"))

panel_a <- patchwork::wrap_plots(
  list(
    plot_plsda_scores(plsda$scores, trimester),
    plot_vip_scores(vip, x_selected, trimester, inputs$annotations)
  ),
  ncol = 1
)
panel_b <- plot_plsda_loadings(plsda$loadings, x_selected, trimester, inputs$annotations, 1L)
panel_c <- plot_top_protein_trajectories(baseline, names(head(vip, 6)), inputs$annotations)
save_plot(panel_a, file.path(paths$plots, "figure_01A_plsda_and_vip.pdf"), 7, 8)
save_plot(panel_b, file.path(paths$plots, "figure_01B_lv1_loadings.pdf"), 4.5, 7)
save_plot(panel_c, file.path(paths$plots, "figure_01C_top_protein_trajectories.pdf"), 7, 7)

gam_cache <- file.path(paths$cache, paste0("baseline_gam_nperm_", parameters$n_perm))
baseline_gam <- run_baseline_gam_analysis(
  inputs$samples,
  inputs$proteins,
  n_perm = parameters$n_perm,
  k = parameters$gam_k,
  cores = parameters$cores,
  seed = parameters$seed,
  cache_dir = gam_cache
)
save_table(baseline_gam$stats, file.path(paths$tables, "baseline_gam_permutation_results.csv"))

significant <- baseline_gam$stats$protein[baseline_gam$stats$fdr < 0.05]
if (length(significant) < 11L) {
  warning("Fewer than 11 proteins passed FDR <0.05; using the top 25 for a diagnostic SOM.")
  significant <- head(baseline_gam$stats$protein[order(baseline_gam$stats$pvalue)], 25L)
}
som_result <- train_self_organizing_map(
  baseline_gam$fitted[significant, , drop = FALSE],
  clusters = 11L,
  seed = parameters$seed
)
cluster_table <- data.frame(
  protein = names(som_result$protein_cluster),
  gene = gene_symbols(names(som_result$protein_cluster), inputs$annotations),
  cluster = as.integer(som_result$protein_cluster)
)
save_table(cluster_table, file.path(paths$tables, "som_protein_clusters.csv"))

save_som_map(som_result, file.path(paths$plots, "figure_01D_self_organizing_map.pdf"))
panel_e <- plot_som_trajectories(som_result)
save_plot(panel_e, file.path(paths$plots, "figure_01E_som_trajectories.pdf"), 8, 6)

representatives <- choose_representative_aptamers(inputs$samples, inputs$proteins, inputs$annotations)
reactome <- run_reactome_ora(som_result$protein_cluster, representatives)
save_table(reactome, file.path(paths$tables, "som_reactome_enrichment.csv"))
reactome_panel <- readr::read_csv(file.path(root, "Figure_1", "pathways_shown.csv"), show_col_types = FALSE)
panel_f <- plot_reactome_panel(reactome, reactome_panel)
save_plot(panel_f, file.path(paths$plots, "figure_01F_reactome_pathways.pdf"), 12, 5)

message("Figure 1 panels written to: ", paths$plots)
