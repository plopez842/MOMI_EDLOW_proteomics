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
paths <- figure_paths(4, root)
cache <- file.path(paths$cache, paste0("nperm_", parameters$n_perm, "_k", parameters$gam_k))

analysis <- run_vaccine_gam_analysis(
  inputs$samples,
  inputs$proteins,
  dose = "V2",
  n_perm = parameters$n_perm,
  k = parameters$gam_k,
  cores = parameters$cores,
  seed = parameters$seed,
  cache_dir = cache
)
statistics <- analysis$stats
statistics$gene <- gene_symbols(statistics$protein, inputs$annotations)
save_table(statistics, file.path(paths$tables, "dose2_vs_baseline_gam_permutation.csv"))

selected <- vaccine_panel_proteins(analysis, inputs$annotations)
if (!length(selected)) selected <- head(statistics$protein[order(statistics$pvalue)], 30L)
masked_log2fc <- mask_unobserved_weeks(analysis$log2fc, analysis$baseline, analysis$vaccinated)
panel_a <- plot_log2fc_heatmap(
  masked_log2fc,
  inputs$annotations,
  selected,
  "Differentially expressed proteins after Dose 2"
)
save_plot(panel_a, file.path(paths$plots, "figure_04A_dose2_log2fc_heatmap.pdf"), 8, 6)

panel_b_proteins <- proteins_for_genes(
  c("CXCL3", "FABP5", "MMP7", "ANXA1", "OLR1", "S100A9"),
  inputs$annotations,
  rownames(analysis$log2fc)
)
panel_b <- plot_vaccine_proteins(analysis$comparison, panel_b_proteins, inputs$annotations, "V2")
save_plot(panel_b, file.path(paths$plots, "figure_04B_dose2_protein_trajectories.pdf"), 11, 5)

representatives <- choose_representative_aptamers(inputs$samples, inputs$proteins, inputs$annotations)
observed_columns <- colSums(!is.na(masked_log2fc)) > 0
gsea <- run_gsea_matrix(analysis$log2fc[, observed_columns, drop = FALSE], representatives, parameters$seed)
save_table(gsea, file.path(paths$tables, "dose2_kegg_gsea_all_weeks.csv"))
pathway_panel <- readr::read_csv(file.path(root, "config", "vaccine_pathway_panel.csv"), show_col_types = FALSE)
panel_results <- filter_pathway_panel(gsea, pathway_panel, "figure4", fdr = 0.25)
save_table(panel_results, file.path(paths$tables, "dose2_kegg_gsea_main_figure_panel.csv"))
panel_c <- plot_gsea_bubbles(panel_results, "Dose 2 versus baseline KEGG GSEA")
save_plot(panel_c, file.path(paths$plots, "figure_04C_dose2_kegg_gsea.pdf"), 11, 10)

message("Figure 4 panels written to: ", paths$plots)
