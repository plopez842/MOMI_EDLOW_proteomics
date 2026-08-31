#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
source(file.path(root, "R", "io.R"))
source_project_functions(root)
create_output_directories(root)

parameters <- analysis_parameters()
set.seed(parameters$seed)
inputs <- load_analysis_inputs(root)
paths <- figure_paths(2, root)

trimester_analysis <- baseline_trimester_limma(inputs$samples, inputs$proteins)
comparison_titles <- c(
  "1st_vs_2nd" = "1st trimester vs 2nd trimester",
  "2nd_vs_3rd" = "2nd trimester vs 3rd trimester",
  "1st_vs_3rd" = "1st trimester vs 3rd trimester"
)

volcanoes <- lapply(names(trimester_analysis$tables), function(comparison) {
  table <- trimester_analysis$tables[[comparison]]
  output <- data.frame(protein = rownames(table), table, check.names = FALSE)
  output$gene <- gene_symbols(output$protein, inputs$annotations)
  save_table(output, file.path(paths$tables, paste0("limma_", comparison, ".csv")))
  volcano_plot(table, comparison_titles[[comparison]], inputs$annotations)
})
names(volcanoes) <- names(trimester_analysis$tables)
panel_abc <- patchwork::wrap_plots(volcanoes, nrow = 1)
save_plot(panel_abc, file.path(paths$plots, "figure_02A-C_trimester_volcanoes.pdf"), 12, 4)

representatives <- choose_representative_aptamers(inputs$samples, inputs$proteins, inputs$annotations)
trimester_matrix <- trimester_contrast_matrix(trimester_analysis$tables)
trimester_gsea <- run_gsea_matrix(trimester_matrix, representatives, parameters$seed)

weekly_fitted <- baseline_weekly_predictions(
  inputs$samples, inputs$proteins,
  k = parameters$gam_k,
  cores = parameters$cores
)
weekly_log2fc <- within_trimester_log2fc(weekly_fitted)
weekly_gsea <- run_gsea_matrix(weekly_log2fc, representatives, parameters$seed + 1000L)

all_gsea <- dplyr::bind_rows(trimester_gsea, weekly_gsea)
save_table(all_gsea, file.path(paths$tables, "kegg_gsea_all_comparisons.csv"))
panel <- readr::read_csv(file.path(root, "config", "figure2_pathway_panel.csv"), show_col_types = FALSE)
panel_results <- filter_pathway_panel(all_gsea, panel, fdr = 0.25)
save_table(panel_results, file.path(paths$tables, "kegg_gsea_main_figure_panel.csv"))

panel_d <- plot_pathway_heatmap(panel_results, "Pathway activity across and within pregnancy trimesters")
save_plot(panel_d, file.path(paths$plots, "figure_02D_kegg_gsea_heatmap.pdf"), 13, 10)

message("Figure 2 panels written to: ", paths$plots)
