choose_representative_aptamers <- function(samples, proteins, annotations) {
  mapping <- annotations[match(proteins, annotations$AptName), , drop = FALSE]
  mapping <- mapping[!is.na(mapping$Single_UniPro) & nzchar(mapping$Single_UniPro), , drop = FALSE]
  if (!nrow(mapping)) stop("No proteins could be mapped to UniProt identifiers.", call. = FALSE)

  if ("F_stat" %in% names(mapping)) {
    mapping$selection_score <- as.numeric(mapping$F_stat)
  } else {
    mapping$selection_score <- vapply(mapping$AptName, function(protein) {
      stats::var(samples[[protein]], na.rm = TRUE)
    }, numeric(1))
  }
  mapping <- mapping[order(mapping$Single_UniPro, -mapping$selection_score), , drop = FALSE]
  mapping[!duplicated(mapping$Single_UniPro), , drop = FALSE]
}

prepare_uniprot_ranks <- function(aptamer_values, representative_mapping) {
  values <- as.numeric(aptamer_values)
  names(values) <- names(aptamer_values)
  mapping <- representative_mapping[
    match(names(values), representative_mapping$AptName),
    c("AptName", "Single_UniPro"),
    drop = FALSE
  ]
  keep <- !is.na(mapping$Single_UniPro) & is.finite(values)
  ranked <- values[keep]
  names(ranked) <- mapping$Single_UniPro[keep]
  ranked <- ranked[!duplicated(names(ranked))]
  sort(ranked, decreasing = TRUE)
}

run_kegg_gsea <- function(ranked_uniprot, seed = 1010L) {
  if (length(ranked_uniprot) < 50L) return(data.frame())
  set.seed(seed)
  result <- suppressMessages(clusterProfiler::gseKEGG(
    geneList = ranked_uniprot,
    organism = "hsa",
    keyType = "uniprot",
    minGSSize = 10,
    maxGSSize = 1400,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  ))
  as.data.frame(result)
}

run_gsea_matrix <- function(value_matrix, representative_mapping, seed = 1010L) {
  value_matrix <- as.matrix(value_matrix)
  results <- lapply(seq_len(ncol(value_matrix)), function(index) {
    values <- value_matrix[, index]
    names(values) <- rownames(value_matrix)
    result <- run_kegg_gsea(prepare_uniprot_ranks(values, representative_mapping), seed + index)
    if (!nrow(result)) return(result)
    result$comparison <- colnames(value_matrix)[index]
    result
  })
  dplyr::bind_rows(results)
}

run_reactome_ora <- function(protein_clusters, representative_mapping) {
  mapping <- representative_mapping[
    representative_mapping$AptName %in% names(protein_clusters),
    , drop = FALSE
  ]
  background <- unique(mapping$Single_UniPro)
  background_entrez <- suppressMessages(clusterProfiler::bitr(
    background,
    fromType = "UNIPROT",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  ))

  results <- lapply(sort(unique(protein_clusters)), function(cluster_number) {
    aptamers <- names(protein_clusters)[protein_clusters == cluster_number]
    uniprot <- unique(mapping$Single_UniPro[mapping$AptName %in% aptamers])
    entrez <- suppressMessages(clusterProfiler::bitr(
      uniprot,
      fromType = "UNIPROT",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    ))
    if (!nrow(entrez)) return(data.frame())
    enrichment <- suppressMessages(ReactomePA::enrichPathway(
      gene = unique(entrez$ENTREZID),
      universe = unique(background_entrez$ENTREZID),
      organism = "human",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = 10,
      maxGSSize = 1400,
      pAdjustMethod = "BH",
      readable = TRUE
    ))
    output <- as.data.frame(enrichment)
    if (nrow(output)) output$cluster <- cluster_number
    output
  })
  dplyr::bind_rows(results)
}

filter_pathway_panel <- function(gsea_results, panel, fdr = 0.25) {
  if (!nrow(gsea_results)) return(gsea_results)
  output <- dplyr::inner_join(
    gsea_results,
    panel,
    by = c("Description" = "pathway")
  )
  output$significant <- is.finite(output$p.adjust) & output$p.adjust < fdr
  output
}

summarize_gsea_directions <- function(gsea_results, fdr = 0.25) {
  gsea_results |>
    dplyr::filter(.data$p.adjust < fdr) |>
    dplyr::mutate(direction = ifelse(.data$NES >= 0, "Positive", "Negative")) |>
    dplyr::count(.data$comparison, .data$direction, name = "n_pathways")
}
