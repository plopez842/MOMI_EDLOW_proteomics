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
# Output:
#   results/figures/figure_02/
#   results/tables/figure_02/
#
# Run from the repository folder with:
#   Rscript Figure_2/figure_2.R


# 1. Project setup ----

source("helpful_functions/project_setup.R")
figure_setup <- prepare_figure(2)

settings <- figure_setup$settings
study_data <- figure_setup$study_data
output <- figure_setup$output


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

  save_results_table(
    comparison_table,
    output,
    paste0("limma_", comparison)
  )

  volcano_plot(
    comparison_results,
    comparison_titles[[comparison]],
    study_data$annotations
  )
})

names(volcano_plots) <- names(trimester_results$tables)
figure_2abc <- patchwork::wrap_plots(volcano_plots, nrow = 1)
save_results_plot(
  figure_2abc,
  output,
  "figure_02A-C_trimester_volcanoes",
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
save_results_table(
  all_gsea_results,
  output,
  "kegg_gsea_all_comparisons"
)


# 4. Select manuscript pathways and make Figure 2D ----

pathways_shown <- pathways_for_figure(2)
figure_pathway_results <- filter_pathway_panel(
  all_gsea_results,
  pathways_shown,
  fdr = 0.25
)
save_results_table(
  figure_pathway_results,
  output,
  "kegg_gsea_main_figure_panel"
)

figure_2d <- plot_pathway_heatmap(
  figure_pathway_results,
  "Pathway activity across and within pregnancy trimesters"
)
save_results_plot(
  figure_2d,
  output,
  "figure_02D_kegg_gsea_heatmap",
  13,
  10
)

message("Figure 2 is complete. Files were written to: ", output$plots)
