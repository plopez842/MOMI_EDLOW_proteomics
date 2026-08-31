#!/usr/bin/env Rscript

# FIGURE 4: PROTEOMIC RESPONSE AFTER VACCINE DOSE 2 ----
#
# Purpose:
#   Compare baseline samples (V0) with samples collected more than 7 days after
#   Dose 2 (V2), while modeling protein abundance across gestational age.
#
# Manuscript panels:
#   A  Heatmap of Dose 2-associated protein changes across gestation
#   B  Example protein trajectories
#   C  KEGG pathway enrichment across gestational weeks
#
# Input:
#   data/processed/comb_v0_vax.rds
#   data/processed/aptamer_annotations.rds
#   Figure_4/pathways_shown.csv
#
# Output:
#   results/figures/figure_04/
#   results/tables/figure_04/
#
# Run from the repository folder with:
#   Rscript Figure_4/figure_4.R


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
output <- figure_output_folders(4, project_directory)


# 2. Dose 2 versus baseline GAM permutation analysis ----

dose_2_cache <- file.path(
  output$cache,
  paste0("nperm_", settings$n_perm, "_k", settings$gam_k)
)

dose_2_results <- run_vaccine_gam_analysis(
  study_data$samples,
  study_data$proteins,
  dose = "V2",
  n_perm = settings$n_perm,
  k = settings$gam_k,
  cores = settings$cores,
  seed = settings$seed,
  cache_dir = dose_2_cache
)

dose_2_statistics <- dose_2_results$stats
dose_2_statistics$gene <- gene_symbols(
  dose_2_statistics$protein,
  study_data$annotations
)
save_table(
  dose_2_statistics,
  file.path(output$tables, "dose2_vs_baseline_gam_permutation.csv")
)


# 3. Protein heatmap across gestation (Figure 4A) ----

proteins_shown <- vaccine_panel_proteins(
  dose_2_results,
  study_data$annotations
)
if (!length(proteins_shown)) {
  proteins_shown <- head(
    dose_2_statistics$protein[order(dose_2_statistics$pvalue)],
    30L
  )
}

observed_log2fc <- mask_unobserved_weeks(
  dose_2_results$log2fc,
  dose_2_results$baseline,
  dose_2_results$vaccinated
)

figure_4a <- plot_log2fc_heatmap(
  observed_log2fc,
  study_data$annotations,
  proteins_shown,
  "Differentially expressed proteins after Dose 2"
)
save_plot(
  figure_4a,
  file.path(output$plots, "figure_04A_dose2_log2fc_heatmap.pdf"),
  8,
  6
)


# 4. Example protein trajectories (Figure 4B) ----

example_proteins <- proteins_for_genes(
  c("CXCL3", "FABP5", "MMP7", "ANXA1", "OLR1", "S100A9"),
  study_data$annotations,
  rownames(dose_2_results$log2fc)
)

figure_4b <- plot_vaccine_proteins(
  dose_2_results$comparison,
  example_proteins,
  study_data$annotations,
  dose = "V2"
)
save_plot(
  figure_4b,
  file.path(output$plots, "figure_04B_dose2_protein_trajectories.pdf"),
  11,
  5
)


# 5. KEGG pathway enrichment across gestation (Figure 4C) ----

representative_aptamers <- choose_representative_aptamers(
  study_data$samples,
  study_data$proteins,
  study_data$annotations
)

weeks_with_data <- colSums(!is.na(observed_log2fc)) > 0
dose_2_gsea <- run_gsea_matrix(
  dose_2_results$log2fc[, weeks_with_data, drop = FALSE],
  representative_aptamers,
  settings$seed
)
save_table(
  dose_2_gsea,
  file.path(output$tables, "dose2_kegg_gsea_all_weeks.csv")
)

pathways_shown <- readr::read_csv(
  file.path(project_directory, "Figure_4", "pathways_shown.csv"),
  show_col_types = FALSE
)
figure_pathway_results <- filter_pathway_panel(
  dose_2_gsea,
  pathways_shown,
  fdr = 0.25
)
save_table(
  figure_pathway_results,
  file.path(output$tables, "dose2_kegg_gsea_main_figure_panel.csv")
)

figure_4c <- plot_gsea_bubbles(
  figure_pathway_results,
  "Dose 2 versus baseline KEGG GSEA"
)
save_plot(
  figure_4c,
  file.path(output$plots, "figure_04C_dose2_kegg_gsea.pdf"),
  11,
  10
)

message("Figure 4 is complete. Files were written to: ", output$plots)
