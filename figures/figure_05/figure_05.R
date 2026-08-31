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
paths <- figure_paths(5, root)

v1 <- run_acute_analysis(
  inputs$samples,
  inputs$proteins,
  dose = "V1",
  n_perm = parameters$n_perm,
  k = parameters$gam_k,
  cores = parameters$cores,
  seed = parameters$seed,
  cache_dir = file.path(paths$cache, paste0("V1_nperm_", parameters$n_perm))
)
v2 <- run_acute_analysis(
  inputs$samples,
  inputs$proteins,
  dose = "V2",
  n_perm = parameters$n_perm,
  k = parameters$gam_k,
  cores = parameters$cores,
  seed = parameters$seed + 500000L,
  cache_dir = file.path(paths$cache, paste0("V2_nperm_", parameters$n_perm))
)

v1_stats <- transform(v1$stats, gene = gene_symbols(v1$stats$protein, inputs$annotations))
v2_stats <- transform(v2$stats, gene = gene_symbols(v2$stats$protein, inputs$annotations))
save_table(v1_stats, file.path(paths$tables, "acute_dose1_lrt_permutation.csv"))
save_table(v2_stats, file.path(paths$tables, "acute_dose2_lrt_permutation.csv"))

v1_panel_proteins <- proteins_for_genes(
  c("VWF", "NOG", "NAAA", "TIMP3", "MRC1", "CXCL12"),
  inputs$annotations,
  names(v1$protein_results)
)
v2_panel_proteins <- proteins_for_genes(
  c("TFF3", "C9", "TNFRSF1B", "IL1RN", "TLR5", "AGT", "CXCL10", "POR"),
  inputs$annotations,
  names(v2$protein_results)
)

panel_a <- plot_acute_proteins(v1, v1_panel_proteins, inputs$annotations)
panel_b <- plot_acute_proteins(v2, v2_panel_proteins, inputs$annotations)
save_plot(panel_a, file.path(paths$plots, "figure_05A_acute_dose1_proteins.pdf"), 7, 7)
save_plot(panel_b, file.path(paths$plots, "figure_05B_acute_dose2_proteins.pdf"), 7, 9)

representatives <- choose_representative_aptamers(inputs$samples, inputs$proteins, inputs$annotations)
v1_daily <- daily_acute_mom_matrix(v1)
v2_daily <- daily_acute_mom_matrix(v2)
v1_gsea <- run_gsea_matrix(v1_daily, representatives, parameters$seed)
v2_gsea <- run_gsea_matrix(v2_daily, representatives, parameters$seed + 1000L)
v1_gsea$dose <- "V1"
v2_gsea$dose <- "V2"
save_table(v1_gsea, file.path(paths$tables, "acute_dose1_kegg_gsea.csv"))
save_table(v2_gsea, file.path(paths$tables, "acute_dose2_kegg_gsea.csv"))

v1_counts <- summarize_gsea_directions(v1_gsea, 0.25)
v2_counts <- summarize_gsea_directions(v2_gsea, 0.25)
panel_c <- patchwork::wrap_plots(
  list(plot_gsea_counts(v1_counts, "V1"), plot_gsea_counts(v2_counts, "V2")),
  nrow = 1,
  guides = "collect"
)
save_plot(panel_c, file.path(paths$plots, "figure_05C_acute_gsea_counts.pdf"), 8, 3.5)

pathway_panel <- readr::read_csv(file.path(root, "config", "vaccine_pathway_panel.csv"), show_col_types = FALSE)
v1_panel <- filter_pathway_panel(v1_gsea, pathway_panel, "figure5", fdr = 0.25)
v2_panel <- filter_pathway_panel(v2_gsea, pathway_panel, "figure5", fdr = 0.25)
save_table(dplyr::bind_rows(V1 = v1_panel, V2 = v2_panel, .id = "dose"), file.path(paths$tables, "acute_kegg_gsea_main_figure_panel.csv"))
panel_d <- patchwork::wrap_plots(
  list(
    plot_gsea_bubbles(v1_panel, "V1 acute KEGG GSEA"),
    plot_gsea_bubbles(v2_panel, "V2 acute KEGG GSEA")
  ),
  nrow = 1,
  guides = "collect"
)
save_plot(panel_d, file.path(paths$plots, "figure_05D_acute_kegg_gsea.pdf"), 12, 9)

message("Figure 5 panels written to: ", paths$plots)
