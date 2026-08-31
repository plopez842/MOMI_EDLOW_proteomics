# PATHWAYS SHOWN IN THE MANUSCRIPT ----
#
# The curated pathway selections used in Figures 1-5 are stored directly in R.
# This keeps the selections in one readable R file while preserving the
# exact pathway names, categories, and display order used by the plots.


# Figure 1: Reactome pathways by SOM cluster ----

figure_1_pathways <- data.frame(
  pathway = c(
    "Toll-like Receptor Cascades",
    "Extracellular matrix organization",
    "Collagen formation",
    "The citric acid (TCA) cycle and respiratory electron transport",
    "PI-3K cascade:FGFR2",
    "Signaling by Receptor Tyrosine Kinases",
    "Cellular responses to stress",
    "Autophagy",
    "Vesicle-mediated transport",
    "KEAP1-NFE2L2 pathway",
    "Signal Transduction",
    "Glycosaminoglycan metabolism",
    "EPH-ephrin mediated repulsion of cells",
    "Dissolution of Fibrin Clot",
    "Sphingolipid metabolism",
    "Toll-like Receptor Cascades",
    "Neutrophil degranulation",
    "Innate Immune System",
    "Metabolism of water-soluble vitamins and cofactors"
  ),
  cluster_order = c(
    rep(1L, 3),
    rep(2L, 8),
    rep(3L, 4),
    rep(9L, 4)
  )
)


# Figure 2: KEGG pathways across pregnancy ----

figure_2_pathways <- data.frame(
  pathway = c(
    "Cellular senescence",
    "Gonadal meiosis",
    "Adherens junction",
    "Focal adhesion",
    "Gap junction",
    "Regulation of actin cytoskeleton",
    "Signaling pathways regulating pluripotency of stem cells",
    "Tight junction",
    "Adrenergic signaling in cardiomyocytes",
    "Vascular smooth muscle contraction",
    "Adipocytokine signaling pathway",
    "Aldosterone-regulated sodium reabsorption",
    "Estrogen signaling pathway",
    "Glucagon signaling pathway",
    "Glycolysis / Gluconeogenesis",
    "GnRH signaling pathway",
    "Growth hormone synthesis secretion and action",
    "Insulin signaling pathway",
    "Oxytocin signaling pathway",
    "Prolactin signaling pathway",
    "Relaxin signaling pathway",
    "Thyroid hormone signaling pathway",
    "B cell receptor signaling pathway",
    "C-type lectin receptor signaling pathway",
    "Fc gamma R-mediated phagocytosis",
    "Natural killer cell mediated cytotoxicity",
    "NOD-like receptor signaling pathway",
    "Platelet activation",
    "T cell receptor signaling pathway",
    "Th17 cell differentiation",
    "Toll-like receptor signaling pathway",
    "Cholinergic synapse",
    "Dopaminergic synapse",
    "Long-term potentiation",
    "Neurotrophin signaling pathway",
    "Serotonergic synapse",
    "AMPK signaling pathway",
    "Apelin signaling pathway",
    "ErbB signaling pathway",
    "Hippo signaling pathway",
    "MAPK signaling pathway",
    "mTOR signaling pathway",
    "NF-kappa B signaling pathway",
    "Sphingolipid signaling pathway",
    "TNF signaling pathway",
    "VEGF signaling pathway",
    "Wnt signaling pathway"
  ),
  category = rep(
    c(
      "Cell growth and death",
      "Cellular motility and signaling",
      "Circulatory system",
      "Endocrine system",
      "Immune system",
      "Nervous system",
      "Signal transduction"
    ),
    times = c(2, 6, 2, 12, 9, 5, 11)
  ),
  display_order = 1:47
)


# Figures 3-5: shared vaccine-response KEGG pathways ----

vaccine_pathways <- data.frame(
  pathway = c(
    "Endocytosis",
    "Efferocytosis",
    "Cellular senescence",
    "Regulation of actin cytoskeleton",
    "Gap junction",
    "Focal adhesion",
    "Tight junction",
    "Thyroid hormone signaling pathway",
    "Relaxin signaling pathway",
    "Prolactin signaling pathway",
    "Oxytocin signaling pathway",
    "Insulin signaling pathway",
    "GnRH signaling pathway",
    "Glucagon signaling pathway",
    "Adipocytokine signaling pathway",
    "NOD-like receptor signaling pathway",
    "IL-17 signaling pathway",
    "Chemokine signaling pathway",
    "Toll-like receptor signaling pathway",
    "Th17 cell differentiation",
    "T cell receptor signaling pathway",
    "Platelet activation",
    "Neutrophil extracellular trap formation",
    "Natural killer cell mediated cytotoxicity",
    "Leukocyte transendothelial migration",
    "Fc gamma R-mediated phagocytosis",
    "Fc epsilon RI signaling pathway",
    "Complement and coagulation cascades",
    "C-type lectin receptor signaling pathway",
    "TNF signaling pathway",
    "Ras signaling pathway",
    "Phospholipase D signaling pathway",
    "NF-kappa B signaling pathway",
    "MAPK signaling pathway",
    "HIF-1 signaling pathway",
    "FoxO signaling pathway",
    "ErbB signaling pathway",
    "Wnt signaling pathway",
    "VEGF signaling pathway",
    "TGF-beta signaling pathway",
    "mTOR signaling pathway",
    "Apelin signaling pathway"
  ),
  category = rep(
    c(
      "Cell growth and death",
      "Cellular motility and signaling",
      "Endocrine system",
      "Immune system",
      "Signal transduction"
    ),
    times = c(3, 4, 8, 14, 13)
  ),
  display_order = 1:42
)

vaccine_pathway_rows <- list(
  figure_3 = c(4:6, 16:18, 30:37),
  figure_4 = c(4:6, 8:42),
  figure_5 = c(1:4, 6:15, 18:30, 32:33, 37:42)
)


# Return the pathway table used by one manuscript figure ----

pathways_for_figure <- function(figure_number) {
  if (!figure_number %in% 1:5) {
    stop("figure_number must be between 1 and 5.", call. = FALSE)
  }

  if (figure_number == 1L) return(figure_1_pathways)
  if (figure_number == 2L) return(figure_2_pathways)

  selected_rows <- vaccine_pathway_rows[[paste0("figure_", figure_number)]]
  selected_pathways <- vaccine_pathways[selected_rows, , drop = FALSE]
  rownames(selected_pathways) <- NULL
  selected_pathways
}
