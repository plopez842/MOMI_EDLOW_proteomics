#!/usr/bin/env Rscript

# FIGURE 5: ACUTE RESPONSE DURING DAYS 1-7 AFTER VACCINATION ----
#
# Purpose:
#   Measure acute protein and pathway changes during the first 7 days after
#   Dose 1 and Dose 2. Each observed protein value is compared with the value
#   expected at the same gestational age from the baseline GAM.
#
# Manuscript panels:
#   A  Acute example proteins after Dose 1
#   B  Acute example proteins after Dose 2
#   C  Number of positively and negatively enriched pathways by day
#   D  KEGG pathway enrichment during days 1-7
#
# Output:
#   results/figures/figure_05/
#   results/tables/figure_05/
#
# Run from the repository folder with:
#   Rscript Figure_5/figure_5.R


# 1. Project setup ----

source("helpful_functions/project_setup.R")
figure_setup <- prepare_figure(5)

settings <- figure_setup$settings
study_data <- figure_setup$study_data
output <- figure_setup$output


# 2. Acute Dose 1 and Dose 2 permutation analyses ----

acute_dose_1 <- run_acute_analysis(
  study_data$samples,
  study_data$proteins,
  dose = "V1",
  n_perm = settings$n_perm,
  k = settings$gam_k,
  cores = settings$cores,
  seed = settings$seed,
  cache_dir = result_cache_directory(
    output,
    paste0("V1_nperm_", settings$n_perm)
  )
)

acute_dose_2 <- run_acute_analysis(
  study_data$samples,
  study_data$proteins,
  dose = "V2",
  n_perm = settings$n_perm,
  k = settings$gam_k,
  cores = settings$cores,
  seed = settings$seed + 500000L,
  cache_dir = result_cache_directory(
    output,
    paste0("V2_nperm_", settings$n_perm)
  )
)

dose_1_statistics <- transform(
  acute_dose_1$stats,
  gene = gene_symbols(acute_dose_1$stats$protein, study_data$annotations)
)
dose_2_statistics <- transform(
  acute_dose_2$stats,
  gene = gene_symbols(acute_dose_2$stats$protein, study_data$annotations)
)

save_results_table(
  dose_1_statistics,
  output,
  "acute_dose1_lrt_permutation"
)
save_results_table(
  dose_2_statistics,
  output,
  "acute_dose2_lrt_permutation"
)


# 3. Example acute protein trajectories (Figure 5A-B) ----

dose_1_example_proteins <- proteins_for_genes(
  c("VWF", "NOG", "NAAA", "TIMP3", "MRC1", "CXCL12"),
  study_data$annotations,
  names(acute_dose_1$protein_results)
)
dose_2_example_proteins <- proteins_for_genes(
  c("TFF3", "C9", "TNFRSF1B", "IL1RN", "TLR5", "AGT", "CXCL10", "POR"),
  study_data$annotations,
  names(acute_dose_2$protein_results)
)

figure_5a <- plot_acute_proteins(
  acute_dose_1,
  dose_1_example_proteins,
  study_data$annotations
)
figure_5b <- plot_acute_proteins(
  acute_dose_2,
  dose_2_example_proteins,
  study_data$annotations
)

save_results_plot(
  figure_5a,
  output,
  "figure_05A_acute_dose1_proteins",
  7,
  7
)
save_results_plot(
  figure_5b,
  output,
  "figure_05B_acute_dose2_proteins",
  7,
  9
)


# 4. Daily KEGG pathway enrichment ----

representative_aptamers <- choose_representative_aptamers(
  study_data$samples,
  study_data$proteins,
  study_data$annotations
)

dose_1_daily_mom <- daily_acute_mom_matrix(acute_dose_1)
dose_2_daily_mom <- daily_acute_mom_matrix(acute_dose_2)

dose_1_gsea <- run_gsea_matrix(
  dose_1_daily_mom,
  representative_aptamers,
  settings$seed
)
dose_2_gsea <- run_gsea_matrix(
  dose_2_daily_mom,
  representative_aptamers,
  settings$seed + 1000L
)
dose_1_gsea$dose <- "V1"
dose_2_gsea$dose <- "V2"

save_results_table(
  dose_1_gsea,
  output,
  "acute_dose1_kegg_gsea"
)
save_results_table(
  dose_2_gsea,
  output,
  "acute_dose2_kegg_gsea"
)


# 5. Number of enriched pathways by day (Figure 5C) ----

dose_1_pathway_counts <- summarize_gsea_directions(dose_1_gsea, fdr = 0.25)
dose_2_pathway_counts <- summarize_gsea_directions(dose_2_gsea, fdr = 0.25)

figure_5c <- patchwork::wrap_plots(
  list(
    plot_gsea_counts(dose_1_pathway_counts, "V1"),
    plot_gsea_counts(dose_2_pathway_counts, "V2")
  ),
  nrow = 1,
  guides = "collect"
)
save_results_plot(
  figure_5c,
  output,
  "figure_05C_acute_gsea_counts",
  8,
  3.5
)


# 6. Selected acute pathways (Figure 5D) ----

pathways_shown <- pathways_for_figure(5)
dose_1_pathways_shown <- filter_pathway_panel(
  dose_1_gsea,
  pathways_shown,
  fdr = 0.25
)
dose_2_pathways_shown <- filter_pathway_panel(
  dose_2_gsea,
  pathways_shown,
  fdr = 0.25
)

save_results_table(
  dplyr::bind_rows(
    V1 = dose_1_pathways_shown,
    V2 = dose_2_pathways_shown,
    .id = "dose"
  ),
  output,
  "acute_kegg_gsea_main_figure_panel"
)

figure_5d <- patchwork::wrap_plots(
  list(
    plot_gsea_bubbles(dose_1_pathways_shown, "V1 acute KEGG GSEA"),
    plot_gsea_bubbles(dose_2_pathways_shown, "V2 acute KEGG GSEA")
  ),
  nrow = 1,
  guides = "collect"
)
save_results_plot(
  figure_5d,
  output,
  "figure_05D_acute_kegg_gsea",
  12,
  9
)

message("Figure 5 is complete. Files were written to: ", output$plots)
