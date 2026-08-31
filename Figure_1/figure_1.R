#!/usr/bin/env Rscript

# FIGURE 1: BASELINE PROTEOME ACROSS PREGNANCY ----
#
# Purpose:
#   Reproduce the baseline gestational-age analyses shown in Figure 1.
#
# Manuscript panels:
#   A-C  Trimester PLSDA, VIP/loadings, and example protein trajectories
#   D-E  Self-organizing map of significant gestational protein trajectories
#   F    Reactome enrichment for the SOM clusters
#
# Input:
#   data/processed/comb_v0_vax.rds
#   data/processed/aptamer_annotations.rds
#   Figure_1/pathways_shown.csv
#
# Output:
#   results/figures/figure_01/
#   results/tables/figure_01/
#
# Run from the repository folder with:
#   Rscript Figure_1/figure_1.R


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
output <- figure_output_folders(1, project_directory)
dir.create(output$cache, recursive = TRUE, showWarnings = FALSE)


# 2. Keep baseline samples and define collection trimester ----

baseline_samples <- study_data$samples[
  study_data$samples$Timepoint_v1 == "V0",
  ,
  drop = FALSE
]

collection_trimester <- factor(
  baseline_samples$Trim_Collec,
  levels = c("Trim_1st", "Trim_2nd", "Trim_3rd"),
  labels = c("1st", "2nd", "3rd")
)

scaled_baseline_proteins <- scale(
  as.matrix(baseline_samples[, study_data$proteins, drop = FALSE])
)
finite_proteins <- colSums(!is.finite(scaled_baseline_proteins)) == 0
scaled_baseline_proteins <- scaled_baseline_proteins[, finite_proteins, drop = FALSE]


# 3. Repeated LASSO feature selection ----

lasso_cache_file <- file.path(
  output$cache,
  paste0("lasso_", settings$lasso_trials, "_trials.rds")
)

if (file.exists(lasso_cache_file)) {
  lasso_selection <- readRDS(lasso_cache_file)
} else {
  lasso_selection <- repeated_lasso_select(
    scaled_baseline_proteins,
    collection_trimester,
    trials = settings$lasso_trials,
    threshold = 0.80,
    seed = settings$seed
  )
  saveRDS(lasso_selection, lasso_cache_file)
}

selected_protein_matrix <- scaled_baseline_proteins[
  ,
  lasso_selection$selected,
  drop = FALSE
]

save_table(
  data.frame(
    protein = names(lasso_selection$frequency),
    selection_frequency = as.numeric(lasso_selection$frequency)
  ),
  file.path(output$tables, "lasso_selection_frequency.csv")
)


# 4. PLSDA and Figure 1A-C ----

plsda_model <- fit_plsda(
  selected_protein_matrix,
  collection_trimester,
  n_components = 2L
)
vip_scores <- sort(plsda_model$vip, decreasing = TRUE)

save_table(
  data.frame(
    protein = names(vip_scores),
    gene = gene_symbols(names(vip_scores), study_data$annotations),
    VIP = as.numeric(vip_scores)
  ),
  file.path(output$tables, "plsda_vip_scores.csv")
)

plsda_loadings_table <- data.frame(
  protein = rownames(plsda_model$loadings),
  plsda_model$loadings,
  check.names = FALSE
)
plsda_loadings_table$gene <- gene_symbols(
  plsda_loadings_table$protein,
  study_data$annotations
)
save_table(
  plsda_loadings_table,
  file.path(output$tables, "plsda_loadings.csv")
)

figure_1a <- patchwork::wrap_plots(
  list(
    plot_plsda_scores(plsda_model$scores, collection_trimester),
    plot_vip_scores(
      vip_scores,
      selected_protein_matrix,
      collection_trimester,
      study_data$annotations
    )
  ),
  ncol = 1
)

figure_1b <- plot_plsda_loadings(
  plsda_model$loadings,
  selected_protein_matrix,
  collection_trimester,
  study_data$annotations,
  component = 1L
)

figure_1c <- plot_top_protein_trajectories(
  baseline_samples,
  names(head(vip_scores, 6)),
  study_data$annotations
)

save_plot(figure_1a, file.path(output$plots, "figure_01A_plsda_and_vip.pdf"), 7, 8)
save_plot(figure_1b, file.path(output$plots, "figure_01B_lv1_loadings.pdf"), 4.5, 7)
save_plot(figure_1c, file.path(output$plots, "figure_01C_top_protein_trajectories.pdf"), 7, 7)


# 5. Gestational GAM permutation analysis ----

baseline_gam_cache <- file.path(
  output$cache,
  paste0("baseline_gam_nperm_", settings$n_perm)
)

baseline_gam_results <- run_baseline_gam_analysis(
  study_data$samples,
  study_data$proteins,
  n_perm = settings$n_perm,
  k = settings$gam_k,
  cores = settings$cores,
  seed = settings$seed,
  cache_dir = baseline_gam_cache
)

save_table(
  baseline_gam_results$stats,
  file.path(output$tables, "baseline_gam_permutation_results.csv")
)


# 6. Self-organizing map and Figure 1D-E ----

som_proteins <- baseline_gam_results$stats$protein[
  baseline_gam_results$stats$fdr < 0.05
]

if (length(som_proteins) < 11L) {
  warning("Fewer than 11 proteins passed FDR < 0.05; using the top 25 proteins.")
  som_proteins <- head(
    baseline_gam_results$stats$protein[
      order(baseline_gam_results$stats$pvalue)
    ],
    25L
  )
}

som_model <- train_self_organizing_map(
  baseline_gam_results$fitted[som_proteins, , drop = FALSE],
  clusters = 11L,
  seed = settings$seed
)

som_cluster_table <- data.frame(
  protein = names(som_model$protein_cluster),
  gene = gene_symbols(names(som_model$protein_cluster), study_data$annotations),
  cluster = as.integer(som_model$protein_cluster)
)
save_table(
  som_cluster_table,
  file.path(output$tables, "som_protein_clusters.csv")
)

save_som_map(
  som_model,
  file.path(output$plots, "figure_01D_self_organizing_map.pdf")
)
figure_1e <- plot_som_trajectories(som_model)
save_plot(figure_1e, file.path(output$plots, "figure_01E_som_trajectories.pdf"), 8, 6)


# 7. Reactome pathway enrichment and Figure 1F ----

representative_aptamers <- choose_representative_aptamers(
  study_data$samples,
  study_data$proteins,
  study_data$annotations
)

reactome_results <- run_reactome_ora(
  som_model$protein_cluster,
  representative_aptamers
)
save_table(
  reactome_results,
  file.path(output$tables, "som_reactome_enrichment.csv")
)

pathways_shown <- readr::read_csv(
  file.path(project_directory, "Figure_1", "pathways_shown.csv"),
  show_col_types = FALSE
)
figure_1f <- plot_reactome_panel(reactome_results, pathways_shown)
save_plot(
  figure_1f,
  file.path(output$plots, "figure_01F_reactome_pathways.pdf"),
  12,
  5
)

message("Figure 1 is complete. Files were written to: ", output$plots)
