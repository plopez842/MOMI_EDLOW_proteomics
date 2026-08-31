#!/usr/bin/env Rscript

# FIGURE 2: PROTEOMIC CHANGES ACROSS PREGNANCY ----
#
# Purpose:
#   Compare baseline protein abundance between trimesters and identify KEGG
#   pathways that change across and within pregnancy trimesters.
#
# Manuscript panels:
#   A-C  Trimester differential-abundance volcano plots
#   D    KEGG GSEA heatmap
#
# Input:
#   data/processed/comb_v0_vax.rds
#   data/processed/aptamer_annotations.rds
#   Figure_2/pathways_shown.csv
#
# Output:
#   results/figures/figure_02/
#   results/tables/figure_02/
#
# Run from the repository folder with:
#   Rscript Figure_2/figure_2.R


# 1. Load functions, settings, and data ----

project_directory <- getwd()
if (!file.exists(file.path(project_directory, "DESCRIPTION"))) {
  stop("Run this script from the MOMI_EDLOW_proteomics repository folder.", call. = FALSE)
}

source(file.path(project_directory, "helpful_functions", "data_and_setup.R"))
source_analysis_functions(project_directory)
create_output_directories(project_directory)

settings <- analysis_parameters()
set.seed(settings$seed)

study_data <- load_analysis_inputs(project_directory)
output <- figure_output_folders(2, project_directory)


# 2. Differential abundance between trimesters (Figure 2A-C) ----

trimester_results <- baseline_trimester_limma(
  study_data$samples,
  study_data$proteins
)

comparison_titles <- c(
  "1st_vs_2nd" = "1st trimester vs 2nd trimester",
  "2nd_vs_3rd" = "2nd trimester vs 3rd trimester",
  "1st_vs_3rd" = "1st trimester vs 3rd trimester"
)

volcano_plots <- lapply(names(trimester_results$tables), function(comparison) {
  comparison_results <- trimester_results$tables[[comparison]]

  comparison_table <- data.frame(
    protein = rownames(comparison_results),
    comparison_results,
    check.names = FALSE
  )
  comparison_table$gene <- gene_symbols(
    comparison_table$protein,
    study_data$annotations
  )

  save_table(
    comparison_table,
    file.path(output$tables, paste0("limma_", comparison, ".csv"))
  )

  volcano_plot(
    comparison_results,
    comparison_titles[[comparison]],
    study_data$annotations
  )
})

names(volcano_plots) <- names(trimester_results$tables)
figure_2abc <- patchwork::wrap_plots(volcano_plots, nrow = 1)
save_plot(
  figure_2abc,
  file.path(output$plots, "figure_02A-C_trimester_volcanoes.pdf"),
  12,
  4
)


# 3. Build protein rankings for KEGG GSEA ----

representative_aptamers <- choose_representative_aptamers(
  study_data$samples,
  study_data$proteins,
  study_data$annotations
)

# Rankings for direct comparisons between trimesters.
between_trimester_log2fc <- trimester_contrast_matrix(trimester_results$tables)
between_trimester_gsea <- run_gsea_matrix(
  between_trimester_log2fc,
  representative_aptamers,
  settings$seed
)

# Rankings for weekly changes within each trimester.
weekly_fitted_abundance <- baseline_weekly_predictions(
  study_data$samples,
  study_data$proteins,
  k = settings$gam_k,
  cores = settings$cores
)
within_trimester_changes <- within_trimester_log2fc(weekly_fitted_abundance)
within_trimester_gsea <- run_gsea_matrix(
  within_trimester_changes,
  representative_aptamers,
  settings$seed + 1000L
)

all_gsea_results <- dplyr::bind_rows(
  between_trimester_gsea,
  within_trimester_gsea
)
save_table(
  all_gsea_results,
  file.path(output$tables, "kegg_gsea_all_comparisons.csv")
)


# 4. Select manuscript pathways and make Figure 2D ----

pathways_shown <- readr::read_csv(
  file.path(project_directory, "Figure_2", "pathways_shown.csv"),
  show_col_types = FALSE
)
figure_pathway_results <- filter_pathway_panel(
  all_gsea_results,
  pathways_shown,
  fdr = 0.25
)
save_table(
  figure_pathway_results,
  file.path(output$tables, "kegg_gsea_main_figure_panel.csv")
)

figure_2d <- plot_pathway_heatmap(
  figure_pathway_results,
  "Pathway activity across and within pregnancy trimesters"
)
save_plot(
  figure_2d,
  file.path(output$plots, "figure_02D_kegg_gsea_heatmap.pdf"),
  13,
  10
)

message("Figure 2 is complete. Files were written to: ", output$plots)
