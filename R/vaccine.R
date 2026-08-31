baseline_trimester_limma <- function(samples, proteins) {
  baseline <- samples[samples$Timepoint_v1 == "V0", , drop = FALSE]
  expression <- t(as.matrix(baseline[, proteins, drop = FALSE]))
  group <- droplevels(baseline$Trim_Collec)
  design <- stats::model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- limma::lmFit(expression, design)
  contrasts <- limma::makeContrasts(
    `1st_vs_2nd` = Trim_1st - Trim_2nd,
    `2nd_vs_3rd` = Trim_2nd - Trim_3rd,
    `1st_vs_3rd` = Trim_1st - Trim_3rd,
    levels = design
  )
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrasts))
  tables <- lapply(seq_len(ncol(contrasts)), function(index) {
    limma::topTable(fit, coef = index, number = Inf, sort.by = "none")
  })
  names(tables) <- colnames(contrasts)
  list(baseline = baseline, tables = tables)
}

baseline_weekly_predictions <- function(samples, proteins, k = 8L, cores = 1L) {
  baseline <- samples[samples$Timepoint_v1 == "V0", , drop = FALSE]
  weeks <- 8:40
  worker <- function(protein) {
    model <- fit_gam_curve(baseline[[protein]], baseline$GA_collection_days, k)
    predict_gam_curve(model, weeks * 7)
  }
  fitted <- do.call(rbind, parallel_lapply(proteins, worker, cores))
  rownames(fitted) <- proteins
  colnames(fitted) <- paste0("W", weeks)
  fitted
}

within_trimester_log2fc <- function(fitted) {
  weeks <- as.integer(sub("W", "", colnames(fitted)))
  starts <- ifelse(weeks < 14, 8L, ifelse(weeks < 28, 14L, 28L))
  output <- vapply(seq_along(weeks), function(index) {
    fitted[, index] - fitted[, match(paste0("W", starts[index]), colnames(fitted))]
  }, numeric(nrow(fitted)))
  output <- t(output)
  rownames(output) <- rownames(fitted)
  colnames(output) <- paste0("W", weeks, "_vs_W", starts)
  output
}

trimester_contrast_matrix <- function(limma_tables) {
  matrix <- do.call(cbind, lapply(limma_tables, function(table) table$logFC))
  rownames(matrix) <- rownames(limma_tables[[1]])
  colnames(matrix) <- names(limma_tables)
  matrix
}

mask_unobserved_weeks <- function(log2fc, baseline, vaccinated) {
  weeks <- as.integer(sub("W", "", colnames(log2fc)))
  baseline_weeks <- floor(baseline$GA_collection_days / 7)
  vaccinated_weeks <- floor(vaccinated$GA_collection_days / 7)
  observed <- weeks %in% baseline_weeks & weeks %in% vaccinated_weeks
  output <- log2fc
  output[, !observed] <- NA_real_
  output
}

proteins_for_genes <- function(genes, annotations, available = NULL) {
  matches <- annotations$AptName[match(genes, annotations$EntrezGeneSymbol)]
  names(matches) <- genes
  matches <- matches[!is.na(matches)]
  if (!is.null(available)) matches <- matches[matches %in% available]
  unname(matches)
}

vaccine_panel_proteins <- function(analysis, annotations) {
  dose <- analysis$dose
  if (dose == "V1") {
    selected <- analysis$stats$protein[analysis$stats$fdr < 0.25]
  } else {
    selected <- analysis$stats$protein[analysis$stats$fdr < 0.10]
    carryover <- proteins_for_genes(c("CXCL3", "PDLIM3", "CASP10", "ZNF174"), annotations, rownames(analysis$log2fc))
    selected <- unique(c(selected, carryover))
  }
  selected
}

daily_acute_mom_matrix <- function(acute_result) {
  days <- sort(unique(acute_result$acute[[acute_result$days_column]]))
  proteins <- names(acute_result$protein_results)
  output <- matrix(NA_real_, nrow = length(proteins), ncol = length(days), dimnames = list(proteins, paste0("D", days)))
  for (protein in proteins) {
    result <- acute_result$protein_results[[protein]]
    for (day in days) {
      output[protein, paste0("D", day)] <- mean(result$mom[result$day == day], na.rm = TRUE)
    }
  }
  output
}

significant_protein_table <- function(stats, annotations, fdr = 0.25) {
  output <- stats[stats$fdr < fdr, , drop = FALSE]
  output$gene <- gene_symbols(output$protein, annotations)
  output[order(output$fdr, -output$statistic), , drop = FALSE]
}
