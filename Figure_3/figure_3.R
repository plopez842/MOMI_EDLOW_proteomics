#!/usr/bin/env Rscript

# FIGURE 3: PROTEOMIC RESPONSE AFTER VACCINE DOSE 1 ----
#
# Purpose:
#   Compare baseline samples (V0) with samples collected more than 7 days after
#   Dose 1 (V1), while modeling protein abundance across gestational age.
#
# Manuscript panels:
#   A  Heatmap of Dose 1-associated protein changes across gestation
#   B  Example protein trajectories
#   C  KEGG pathway enrichment across gestational weeks
#
# Input:
#   data/processed/comb_v0_vax.rds
#   data/processed/aptamer_annotations.rds
#   Figure_3/pathways_shown.csv
#
# Output:
#   results/figures/figure_03/
#   results/tables/figure_03/
#
# Run from the repository folder with:
#   Rscript Figure_3/figure_3.R


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
output <- figure_output_folders(3, project_directory)


# 2. Dose 1 versus baseline GAM permutation analysis ----

dose_1_cache <- file.path(
  output$cache,
  paste0("nperm_", settings$n_perm, "_k", settings$gam_k)
)

dose_1_results <- run_vaccine_gam_analysis(
  study_data$samples,
  study_data$proteins,
  dose = "V1",
  n_perm = settings$n_perm,
  k = settings$gam_k,
  cores = settings$cores,
  seed = settings$seed,
  cache_dir = dose_1_cache
)

dose_1_statistics <- dose_1_results$stats
dose_1_statistics$gene <- gene_symbols(
  dose_1_statistics$protein,
  study_data$annotations
)
save_table(
  dose_1_statistics,
  file.path(output$tables, "dose1_vs_baseline_gam_permutation.csv")
)


# 3. Protein heatmap across gestation (Figure 3A) ----

proteins_shown <- vaccine_panel_proteins(
  dose_1_results,
  study_data$annotations
)
if (!length(proteins_shown)) {
  proteins_shown <- head(
    dose_1_statistics$protein[order(dose_1_statistics$pvalue)],
    13L
  )
}

observed_log2fc <- mask_unobserved_weeks(
  dose_1_results$log2fc,
  dose_1_results$baseline,
  dose_1_results$vaccinated
)

figure_3a <- plot_log2fc_heatmap(
  observed_log2fc,
  study_data$annotations,
  proteins_shown,
  "Differentially expressed proteins after Dose 1"
)
save_plot(
  figure_3a,
  file.path(output$plots, "figure_03A_dose1_log2fc_heatmap.pdf"),
  8,
  4
)


# 4. Example protein trajectories (Figure 3B) ----

example_proteins <- proteins_for_genes(
  c("IL1RL1", "CXCL3", "ZNF174"),
  study_data$annotations,
  rownames(dose_1_results$log2fc)
)

figure_3b <- plot_vaccine_proteins(
  dose_1_results$comparison,
  example_proteins,
  study_data$annotations,
  dose = "V1"
)
save_plot(
  figure_3b,
  file.path(output$plots, "figure_03B_dose1_protein_trajectories.pdf"),
  9,
  3
)


# 5. KEGG pathway enrichment across gestation (Figure 3C) ----

representative_aptamers <- choose_representative_aptamers(
  study_data$samples,
  study_data$proteins,
  study_data$annotations
)

weeks_with_data <- colSums(!is.na(observed_log2fc)) > 0
dose_1_gsea <- run_gsea_matrix(
  dose_1_results$log2fc[, weeks_with_data, drop = FALSE],
  representative_aptamers,
  settings$seed
)
save_table(
  dose_1_gsea,
  file.path(output$tables, "dose1_kegg_gsea_all_weeks.csv")
)

pathways_shown <- readr::read_csv(
  file.path(project_directory, "Figure_3", "pathways_shown.csv"),
  show_col_types = FALSE
)
figure_pathway_results <- filter_pathway_panel(
  dose_1_gsea,
  pathways_shown,
  fdr = 0.25
)
save_table(
  figure_pathway_results,
  file.path(output$tables, "dose1_kegg_gsea_main_figure_panel.csv")
)

figure_3c <- plot_gsea_bubbles(
  figure_pathway_results,
  "Dose 1 versus baseline KEGG GSEA"
)
save_plot(
  figure_3c,
  file.path(output$plots, "figure_03C_dose1_kegg_gsea.pdf"),
  10,
  6
)

message("Figure 3 is complete. Files were written to: ", output$plots)
