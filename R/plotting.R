trimester_colors <- c("1st" = "#6A9BD4", "2nd" = "#F5C242", "3rd" = "#E76F51")
vaccine_colors <- c("V0" = "#A7A9AC", "V1" = "#1B9E77", "V2" = "#756BB1")
nes_colors <- c(low = "#3B2C85", mid = "#F7F7F7", high = "#D7301F")

theme_manuscript <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )
}

save_plot <- function(plot, path, width, height, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = dpi)
  invisible(path)
}

plot_plsda_scores <- function(scores, trimester) {
  columns <- names(scores)[seq_len(min(2L, ncol(scores)))]
  plot_data <- data.frame(
    LV1 = scores[[columns[1]]],
    LV2 = scores[[columns[2]]],
    trimester = factor(trimester, levels = c("1st", "2nd", "3rd"))
  )
  ggplot2::ggplot(plot_data, ggplot2::aes(.data$LV1, .data$LV2, colour = .data$trimester)) +
    ggplot2::stat_ellipse(level = 0.95, linewidth = 0.7) +
    ggplot2::geom_point(size = 2.3, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = trimester_colors, name = "Trimester") +
    ggplot2::labs(x = "Scores on LV1", y = "Scores on LV2") +
    theme_manuscript()
}

plot_vip_scores <- function(vip, x_scaled, trimester, annotations, minimum = 1) {
  vip_data <- data.frame(protein = names(vip), VIP = as.numeric(vip))
  vip_data <- vip_data[vip_data$VIP >= minimum, , drop = FALSE]
  if (!nrow(vip_data)) vip_data <- head(vip_data[order(-vip_data$VIP), , drop = FALSE], 20)
  vip_data$gene <- gene_symbols(vip_data$protein, annotations)
  vip_data$enriched <- vapply(vip_data$protein, function(protein) {
    medians <- tapply(x_scaled[, protein], trimester, stats::median, na.rm = TRUE)
    names(which.max(medians))
  }, character(1))
  vip_data <- vip_data[order(-vip_data$VIP), , drop = FALSE]
  vip_data$gene <- factor(vip_data$gene, levels = vip_data$gene)
  ggplot2::ggplot(vip_data, ggplot2::aes(.data$gene, .data$VIP, fill = .data$enriched)) +
    ggplot2::geom_col(colour = "black", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = trimester_colors, guide = "none") +
    ggplot2::labs(x = NULL, y = "VIP score") +
    theme_manuscript() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7))
}

plot_plsda_loadings <- function(loadings, x_scaled, trimester, annotations, component = 1L) {
  column <- names(loadings)[component]
  plot_data <- data.frame(
    protein = rownames(loadings),
    loading = loadings[[column]]
  )
  plot_data$gene <- gene_symbols(plot_data$protein, annotations)
  plot_data$enriched <- vapply(plot_data$protein, function(protein) {
    medians <- tapply(x_scaled[, protein], trimester, stats::median, na.rm = TRUE)
    names(which.max(medians))
  }, character(1))
  plot_data <- plot_data[order(plot_data$loading), , drop = FALSE]
  plot_data$gene <- factor(plot_data$gene, levels = plot_data$gene)
  ggplot2::ggplot(plot_data, ggplot2::aes(.data$loading, .data$gene, fill = .data$enriched)) +
    ggplot2::geom_col(colour = "black", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = trimester_colors, guide = "none") +
    ggplot2::labs(x = paste0("LV", component, " loadings"), y = NULL) +
    theme_manuscript(base_size = 8)
}

plot_top_protein_trajectories <- function(samples, proteins, annotations, n = 6L) {
  selected <- head(proteins, n)
  plots <- lapply(selected, function(protein) {
    data <- data.frame(
      ga = samples$GA_collection_days,
      value = samples[[protein]],
      trimester = trimester_from_weeks(samples$GA_collection_days / 7)
    )
    ggplot2::ggplot(data, ggplot2::aes(.data$ga, .data$value, colour = .data$trimester)) +
      ggplot2::geom_point(size = 1.4, alpha = 0.8) +
      ggplot2::geom_smooth(
        method = "gam", formula = y ~ s(x, bs = "cr", k = 8),
        se = FALSE, colour = "black", linewidth = 0.7
      ) +
      ggplot2::geom_vline(xintercept = c(14, 28) * 7, linetype = "dashed") +
      ggplot2::scale_colour_manual(values = trimester_colors, guide = "none") +
      ggplot2::labs(title = gene_symbols(protein, annotations), x = NULL, y = "log2 RFU") +
      theme_manuscript(base_size = 8)
  })
  patchwork::wrap_plots(plots, ncol = 2) + patchwork::plot_annotation()
}

save_som_map <- function(som_result, path, palette = NULL) {
  clusters <- sort(unique(som_result$unit_cluster))
  if (is.null(palette)) palette <- grDevices::hcl.colors(length(clusters), "Spectral")
  grDevices::pdf(path, width = 7, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(
    som_result$model,
    type = "codes",
    bgcol = palette[som_result$unit_cluster],
    main = "SOM clusters",
    shape = "straight"
  )
  kohonen::add.cluster.boundaries(som_result$model, som_result$unit_cluster)
  invisible(path)
}

plot_som_trajectories <- function(som_result, weeks = 8:40) {
  fitted <- som_result$scaled_fitted
  long <- as.data.frame(fitted)
  long$protein <- rownames(long)
  long$cluster <- factor(som_result$protein_cluster[long$protein])
  long <- tidyr::pivot_longer(long, -c(.data$protein, .data$cluster), names_to = "week", values_to = "scaled")
  long$week <- rep(weeks, times = nrow(fitted))
  ggplot2::ggplot(long, ggplot2::aes(.data$week, .data$scaled, group = .data$protein, colour = .data$protein)) +
    ggplot2::geom_smooth(se = FALSE, linewidth = 0.45) +
    ggplot2::facet_wrap(~cluster, ncol = 4) +
    ggplot2::geom_vline(xintercept = c(14, 28), linetype = "dashed") +
    ggplot2::guides(colour = "none") +
    ggplot2::labs(title = "SOM clustering", x = "Gestational age (weeks)", y = "Scaled fitted abundance") +
    theme_manuscript(base_size = 8)
}

plot_reactome_panel <- function(enrichment, panel) {
  data <- dplyr::inner_join(enrichment, panel, by = c("Description" = "pathway", "cluster" = "cluster_order"))
  data <- data[data$p.adjust < 0.25, , drop = FALSE]
  data$Description <- factor(data$Description, levels = rev(unique(data$Description)))
  ggplot2::ggplot(data, ggplot2::aes(.data$Description, -log10(.data$pvalue), fill = factor(.data$cluster))) +
    ggplot2::geom_col(colour = "black", linewidth = 0.25) +
    ggplot2::coord_flip() +
    ggplot2::facet_grid(. ~ cluster, scales = "free_x", space = "free_x") +
    ggplot2::labs(title = "SOM cluster Reactome pathway analysis", x = NULL, y = "-log10(p-value)", fill = "Cluster") +
    theme_manuscript(base_size = 8) +
    ggplot2::theme(legend.position = "none")
}

volcano_plot <- function(table, title, annotations, threshold = log2(1.5), fdr = 0.05, labels = 12L) {
  data <- as.data.frame(table)
  data$protein <- rownames(data)
  data$gene <- gene_symbols(data$protein, annotations)
  data$significant <- data$adj.P.Val < fdr & abs(data$logFC) > threshold
  rank_index <- order(data$adj.P.Val, -abs(data$logFC))
  data$label <- ""
  data$label[head(rank_index[data$significant[rank_index]], labels)] <- data$gene[head(rank_index[data$significant[rank_index]], labels)]

  ggplot2::ggplot(data, ggplot2::aes(.data$logFC, -log10(.data$adj.P.Val))) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$significant), size = 1.1) +
    ggrepel::geom_text_repel(ggplot2::aes(label = .data$label), size = 2.6, max.overlaps = Inf) +
    ggplot2::geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed", colour = "grey70") +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = "dashed", colour = "grey70") +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "#D7301F"), guide = "none") +
    ggplot2::labs(title = title, x = "log2 fold change", y = "-log10(FDR)") +
    theme_manuscript(base_size = 9)
}

plot_log2fc_heatmap <- function(log2fc, annotations, selected_proteins, title, limits = NULL) {
  selected_proteins <- intersect(selected_proteins, rownames(log2fc))
  matrix <- log2fc[selected_proteins, , drop = FALSE]
  long <- as.data.frame(matrix)
  long$protein <- rownames(long)
  long <- tidyr::pivot_longer(long, -"protein", names_to = "week", values_to = "log2fc")
  long$week_number <- as.integer(sub("W", "", long$week))
  long$gene <- gene_symbols(long$protein, annotations)
  long$gene <- factor(long$gene, levels = rev(gene_symbols(selected_proteins, annotations)))
  if (is.null(limits)) limits <- max(abs(long$log2fc), na.rm = TRUE) * c(-1, 1)

  ggplot2::ggplot(long, ggplot2::aes(.data$week_number, .data$gene, fill = .data$log2fc)) +
    ggplot2::geom_tile(colour = "grey55", linewidth = 0.15) +
    ggplot2::scale_fill_gradient2(
      low = nes_colors[["low"]], mid = nes_colors[["mid"]], high = nes_colors[["high"]],
      midpoint = 0, limits = limits, na.value = "grey80", name = "log2FC"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(8, 40, 4), expand = c(0, 0)) +
    ggplot2::labs(title = title, x = "Gestational age at collection (weeks)", y = NULL) +
    theme_manuscript(base_size = 8)
}

plot_pathway_heatmap <- function(panel_results, title = "KEGG pathway enrichment") {
  data <- panel_results
  data$comparison <- factor(data$comparison, levels = unique(data$comparison))
  data$Description <- factor(data$Description, levels = rev(unique(data$Description[order(data$display_order)])))
  ggplot2::ggplot(data, ggplot2::aes(.data$comparison, .data$Description, fill = .data$NES)) +
    ggplot2::geom_tile(colour = "grey60", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(.data$significant, "*", "")), size = 2.5) +
    ggplot2::facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient2(
      low = nes_colors[["low"]], mid = nes_colors[["mid"]], high = nes_colors[["high"]], midpoint = 0
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = "NES") +
    theme_manuscript(base_size = 7) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
}

plot_gsea_bubbles <- function(panel_results, title = "KEGG GSEA") {
  data <- panel_results[panel_results$significant, , drop = FALSE]
  data$comparison <- factor(data$comparison, levels = unique(data$comparison))
  data$Description <- factor(data$Description, levels = rev(unique(data$Description[order(data$display_order)])))
  ggplot2::ggplot(data, ggplot2::aes(.data$comparison, .data$Description)) +
    ggplot2::geom_point(ggplot2::aes(size = -log10(.data$pvalue), colour = .data$NES), stroke = 0.2) +
    ggplot2::facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_colour_gradient2(
      low = nes_colors[["low"]], mid = nes_colors[["mid"]], high = nes_colors[["high"]], midpoint = 0
    ) +
    ggplot2::scale_size_continuous(range = c(1.2, 5)) +
    ggplot2::labs(title = title, x = NULL, y = NULL, size = "-log10(p)", colour = "NES") +
    theme_manuscript(base_size = 7) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
}

plot_vaccine_proteins <- function(comparison, proteins, annotations, dose) {
  plots <- lapply(proteins, function(protein) {
    data <- data.frame(
      week = comparison$GA_collection_days / 7,
      value = comparison[[protein]],
      status = droplevels(factor(comparison$Timepoint_v1, levels = c("V0", dose)))
    )
    ggplot2::ggplot(data, ggplot2::aes(.data$week, .data$value, colour = .data$status)) +
      ggplot2::geom_point(size = 0.9, alpha = 0.55) +
      ggplot2::geom_smooth(
        method = "gam", formula = y ~ s(x, bs = "cr", k = 8),
        se = FALSE, linewidth = 0.9
      ) +
      ggplot2::scale_colour_manual(values = vaccine_colors, name = "Vaccine status") +
      ggplot2::scale_x_continuous(breaks = seq(8, 40, 5)) +
      ggplot2::labs(title = gene_symbols(protein, annotations), x = "GA at collection (weeks)", y = "log2 RFU") +
      theme_manuscript(base_size = 8)
  })
  patchwork::wrap_plots(plots, nrow = 1, guides = "collect")
}

plot_gsea_counts <- function(counts, dose) {
  ggplot2::ggplot(counts, ggplot2::aes(.data$comparison, .data$n_pathways, fill = .data$direction)) +
    ggplot2::geom_col(position = "identity", colour = "black") +
    ggplot2::scale_fill_manual(values = c("Negative" = "#24126A", "Positive" = "#C00000")) +
    ggplot2::labs(title = paste("GSEA analysis for", dose), x = paste("Days post", dose), y = "Number of pathways", fill = "NES") +
    theme_manuscript(base_size = 9)
}

plot_acute_proteins <- function(acute_result, proteins, annotations) {
  plots <- lapply(proteins, function(protein) {
    result <- acute_result$protein_results[[protein]]
    acute_rows <- result$keep
    data <- data.frame(
      day = result$day,
      mom = result$mom,
      trimester = factor(acute_result$acute$Vax1_Trim[acute_rows], levels = c("Trim_1st", "Trim_2nd", "Trim_3rd"))
    )
    levels(data$trimester) <- c("1st", "2nd", "3rd")
    ggplot2::ggplot(data, ggplot2::aes(.data$day, .data$mom, colour = .data$trimester)) +
      ggplot2::geom_point(size = 1.5, alpha = 0.9) +
      ggplot2::geom_smooth(
        method = "gam", formula = y ~ s(x, bs = "cr", k = min(6, length(unique(x)) - 1)),
        se = FALSE, colour = "black", linewidth = 0.7
      ) +
      ggplot2::scale_colour_manual(values = trimester_colors, name = "Trimester of Dose 1") +
      ggplot2::scale_x_continuous(breaks = 1:7) +
      ggplot2::labs(title = gene_symbols(protein, annotations), x = paste("Days post", acute_result$dose), y = "log2(MoM)") +
      theme_manuscript(base_size = 8)
  })
  patchwork::wrap_plots(plots, ncol = 2, guides = "collect")
}
