# STATISTICAL MODELING FUNCTIONS ----
#
# This file contains the models shared by multiple manuscript figures:
#   - repeated LASSO feature selection and PLSDA
#   - gestational-age generalized additive models (GAMs)
#   - label-permutation significance tests
#   - self-organizing map (SOM) clustering
#   - acute post-vaccination median-of-medians (MoM) models
#
# Figure scripts provide the biological comparisons and call these functions.


# Parallel processing helper ----

# Use multiple macOS/Linux processes when MOMI_CORES is greater than one.
# Windows and single-core runs use ordinary lapply().
parallel_lapply <- function(x, fun, cores = 1L) {
  if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(x, fun, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(x, fun)
  }
}


# LASSO feature selection and PLSDA ----

# Repeat cross-validated LASSO and retain proteins selected in at least the
# requested fraction of trials.
repeated_lasso_select <- function(x, y, trials = 100L, threshold = 0.8, seed = 1010L) {
  x <- as.matrix(x)
  y <- droplevels(factor(y))
  if (nlevels(y) < 2L) stop("LASSO requires at least two outcome classes.", call. = FALSE)
  family <- if (nlevels(y) == 2L) "binomial" else "multinomial"
  counts <- stats::setNames(integer(ncol(x)), colnames(x))
  nfolds <- min(5L, nrow(x))

  for (trial in seq_len(trials)) {
    set.seed(seed + trial)
    foldid <- sample(rep(seq_len(nfolds), length.out = nrow(x)))
    fit <- glmnet::cv.glmnet(
      x, y,
      family = family,
      alpha = 1,
      nfolds = nfolds,
      foldid = foldid,
      type.measure = "mse",
      type.multinomial = "grouped",
      standardize = FALSE
    )
    coefficients <- stats::coef(fit, s = "lambda.min")
    if (!is.list(coefficients)) coefficients <- list(coefficients)
    selected <- unique(unlist(lapply(coefficients, function(coef_matrix) {
      nonzero <- which(as.numeric(coef_matrix[-1, 1]) != 0)
      rownames(coef_matrix)[-1][nonzero]
    })))
    counts[selected] <- counts[selected] + 1L
  }

  frequency <- counts / trials
  selected <- names(frequency)[frequency >= threshold]
  if (!length(selected)) selected <- names(frequency)[frequency == max(frequency)]
  list(selected = selected, frequency = sort(frequency, decreasing = TRUE))
}

# Fit a two-component PLSDA model and return the scores, loadings, and VIPs in
# plain R objects that are easy to save and plot.
fit_plsda <- function(x, y, n_components = 2L) {
  x <- as.matrix(x)
  y <- droplevels(factor(y))
  model <- ropls::opls(
    x, y,
    predI = n_components,
    orthoI = 0,
    permI = 0,
    crossvalI = 5,
    plotSubC = "none",
    fig.pdfC = "none",
    info.txtC = "none"
  )
  list(
    model = model,
    scores = as.data.frame(ropls::getScoreMN(model)),
    loadings = as.data.frame(ropls::getLoadingMN(model)),
    vip = ropls::getVipVn(model)
  )
}


# Generalized additive model helpers ----

# Limit spline complexity when only a few distinct time values are available.
effective_k <- function(x, requested = 8L) {
  max(3L, min(as.integer(requested), length(unique(stats::na.omit(x))) - 1L))
}

# Fit one REML cubic-regression spline across gestational age in days.
fit_gam_curve <- function(y, ga_days, k = 8L) {
  keep <- is.finite(y) & is.finite(ga_days)
  model_data <- data.frame(y = as.numeric(y[keep]), ga = as.numeric(ga_days[keep]))
  if (nrow(model_data) < 8L || length(unique(model_data$ga)) < 4L) {
    stop("Insufficient observations or gestational ages for GAM fitting.", call. = FALSE)
  }
  k_eff <- effective_k(model_data$ga, k)
  mgcv::gam(y ~ s(ga, bs = "cr", k = k_eff), method = "REML", data = model_data)
}

# Predict a fitted gestational GAM at specific gestational ages.
predict_gam_curve <- function(model, ga_days) {
  as.numeric(stats::predict(model, newdata = data.frame(ga = ga_days), type = "response"))
}

# Numerically integrate a curve using the trapezoid rule.
trapezoid_sum <- function(x, y) {
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}


# Vaccine response across gestation ----

# For one protein, fit separate baseline and vaccinated GAMs. The test statistic
# is the integrated absolute distance between curves. Label permutation creates
# the empirical null distribution.
fit_one_vaccine_protein <- function(y, ga_days, group, n_perm = 1000L, k = 8L, seed = 1010L) {
  keep <- is.finite(y) & is.finite(ga_days) & !is.na(group)
  y <- as.numeric(y[keep])
  ga_days <- as.numeric(ga_days[keep])
  group <- droplevels(factor(group[keep]))
  if (nlevels(group) != 2L) stop("Vaccine comparison requires exactly two groups.", call. = FALSE)

  levels_group <- levels(group)
  fit0 <- fit_gam_curve(y[group == levels_group[1]], ga_days[group == levels_group[1]], k)
  fit1 <- fit_gam_curve(y[group == levels_group[2]], ga_days[group == levels_group[2]], k)
  prediction_days <- seq(8 * 7, 40 * 7, by = 1)
  prediction_weeks <- 8:40
  weekly_days <- prediction_weeks * 7

  pred0_daily <- predict_gam_curve(fit0, prediction_days)
  pred1_daily <- predict_gam_curve(fit1, prediction_days)
  observed <- trapezoid_sum(prediction_days, abs(pred1_daily - pred0_daily))

  set.seed(seed)
  null <- numeric(n_perm)
  for (permutation in seq_len(n_perm)) {
    permuted_group <- sample(group)
    null[permutation] <- tryCatch({
      perm0 <- fit_gam_curve(y[permuted_group == levels_group[1]], ga_days[permuted_group == levels_group[1]], k)
      perm1 <- fit_gam_curve(y[permuted_group == levels_group[2]], ga_days[permuted_group == levels_group[2]], k)
      pred0 <- predict_gam_curve(perm0, prediction_days)
      pred1 <- predict_gam_curve(perm1, prediction_days)
      trapezoid_sum(prediction_days, abs(pred1 - pred0))
    }, error = function(error) NA_real_)
  }
  valid_null <- null[is.finite(null)]
  pvalue <- if (length(valid_null)) mean(valid_null >= observed) else NA_real_

  weekly0 <- predict_gam_curve(fit0, weekly_days)
  weekly1 <- predict_gam_curve(fit1, weekly_days)
  list(
    observed = observed,
    pvalue = pvalue,
    null = null,
    weeks = prediction_weeks,
    baseline = weekly0,
    vaccinated = weekly1,
    log2fc = weekly1 - weekly0
  )
}

# Apply the vaccine GAM permutation test to every protein for Dose 1 or Dose 2.
# Optional per-protein cache files allow long analyses to resume.
run_vaccine_gam_analysis <- function(samples, proteins, dose = c("V1", "V2"), n_perm = 1000L,
                                     k = 8L, cores = 1L, seed = 1010L, cache_dir = NULL) {
  dose <- match.arg(dose)
  days_column <- if (dose == "V1") "diff_collec_vax1" else "diff_collec_vax2"
  baseline <- samples[samples$Timepoint_v1 == "V0", , drop = FALSE]
  vaccinated <- samples[
    samples$Timepoint_v1 == dose & samples[[days_column]] > 7,
    , drop = FALSE
  ]
  comparison <- rbind(baseline, vaccinated)
  group <- factor(as.character(comparison$Timepoint_v1), levels = c("V0", dose))

  if (!is.null(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  worker <- function(index) {
    protein <- proteins[index]
    cache_file <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(protein, ".rds")) else NULL
    if (!is.null(cache_file) && file.exists(cache_file)) return(readRDS(cache_file))
    result <- fit_one_vaccine_protein(
      comparison[[protein]], comparison$GA_collection_days, group,
      n_perm = n_perm, k = k, seed = seed + index * 10000L
    )
    result$protein <- protein
    if (!is.null(cache_file)) saveRDS(result, cache_file)
    result
  }

  results <- parallel_lapply(seq_along(proteins), worker, cores)
  stats <- data.frame(
    protein = proteins,
    auc_difference = vapply(results, `[[`, numeric(1), "observed"),
    pvalue = vapply(results, `[[`, numeric(1), "pvalue")
  )
  stats$fdr <- stats::p.adjust(stats$pvalue, method = "BH")
  log2fc <- do.call(rbind, lapply(results, `[[`, "log2fc"))
  rownames(log2fc) <- proteins
  colnames(log2fc) <- paste0("W", 8:40)
  list(
    dose = dose,
    baseline = baseline,
    vaccinated = vaccinated,
    comparison = comparison,
    stats = stats,
    log2fc = log2fc,
    protein_results = stats::setNames(results, proteins)
  )
}


# Baseline gestational trajectories ----

# For one baseline protein, test whether its GAM explains more variation than
# expected after permuting protein abundance across gestational ages.
fit_one_baseline_gam <- function(y, ga_days, prediction_days, n_perm = 1000L, k = 8L, seed = 1010L) {
  model <- fit_gam_curve(y, ga_days, k)
  observed <- summary(model)$dev.expl
  set.seed(seed)
  null <- numeric(n_perm)
  for (permutation in seq_len(n_perm)) {
    null[permutation] <- tryCatch({
      permuted <- fit_gam_curve(sample(y), ga_days, k)
      summary(permuted)$dev.expl
    }, error = function(error) NA_real_)
  }
  valid_null <- null[is.finite(null)]
  list(
    deviance_explained = observed,
    pvalue = if (length(valid_null)) mean(valid_null >= observed) else NA_real_,
    fitted = predict_gam_curve(model, prediction_days)
  )
}

# Run the baseline gestational GAM permutation test for every protein.
run_baseline_gam_analysis <- function(samples, proteins, n_perm = 1000L, k = 8L,
                                      cores = 1L, seed = 1010L, cache_dir = NULL) {
  baseline <- samples[samples$Timepoint_v1 == "V0", , drop = FALSE]
  prediction_days <- seq(8 * 7, 40 * 7, by = 7)
  if (!is.null(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  worker <- function(index) {
    protein <- proteins[index]
    cache_file <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(protein, ".rds")) else NULL
    if (!is.null(cache_file) && file.exists(cache_file)) return(readRDS(cache_file))
    result <- fit_one_baseline_gam(
      baseline[[protein]], baseline$GA_collection_days, prediction_days,
      n_perm = n_perm, k = k, seed = seed + index * 10000L
    )
    result$protein <- protein
    if (!is.null(cache_file)) saveRDS(result, cache_file)
    result
  }

  results <- parallel_lapply(seq_along(proteins), worker, cores)
  stats <- data.frame(
    protein = proteins,
    deviance_explained = vapply(results, `[[`, numeric(1), "deviance_explained"),
    pvalue = vapply(results, `[[`, numeric(1), "pvalue")
  )
  stats$fdr <- stats::p.adjust(stats$pvalue, method = "BH")
  fitted <- do.call(rbind, lapply(results, `[[`, "fitted"))
  rownames(fitted) <- proteins
  colnames(fitted) <- paste0("W", 8:40)
  list(baseline = baseline, stats = stats, fitted = fitted)
}


# Self-organizing map ----

# Cluster standardized fitted gestational trajectories. The publication uses an
# 11 x 11 SOM; smaller maps are used automatically for reduced tests.
train_self_organizing_map <- function(fitted_matrix, clusters = 11L, seed = 1010L) {
  set.seed(seed)
  x <- t(scale(t(as.matrix(fitted_matrix))))
  grid_side <- min(11L, max(2L, floor(sqrt(nrow(x)))))
  grid <- kohonen::somgrid(grid_side, grid_side, topo = "hexagonal")
  radius_start <- 2 * grid_side / 3
  som1 <- kohonen::som(
    x, grid = grid, rlen = 1000, alpha = c(0.5, 0.5),
    radius = c(radius_start, 0), mode = "online", keep.data = TRUE
  )
  som2 <- kohonen::som(
    x, grid = grid, rlen = 2000, alpha = c(0.1, 0.1),
    radius = c(radius_start, 0), mode = "online", init = som1$codes, keep.data = TRUE
  )
  som3 <- kohonen::som(
    x, grid = grid, rlen = 5000, alpha = c(0.01, 0.01),
    radius = c(radius_start, 0), mode = "online", init = som2$codes, keep.data = TRUE
  )
  code_distance <- kohonen::object.distances(som3, type = "codes")
  tree <- stats::hclust(code_distance, method = "ward.D2")
  unit_cluster <- stats::cutree(tree, k = min(as.integer(clusters), grid_side^2))
  protein_cluster <- unit_cluster[som3$unit.classif]
  names(protein_cluster) <- rownames(x)
  list(
    model = som3,
    scaled_fitted = x,
    tree = tree,
    unit_cluster = unit_cluster,
    protein_cluster = protein_cluster
  )
}


# Acute response during days 1-7 ----

# For one protein, subtract the gestational-age-specific baseline expectation
# from each acute sample, then test whether MoM varies across days 1-7.
fit_acute_mom <- function(baseline_y, baseline_ga, acute_y, acute_ga, acute_day,
                          n_perm = 1000L, k = 8L, seed = 1010L) {
  baseline_model <- fit_gam_curve(baseline_y, baseline_ga, k)
  expected <- predict_gam_curve(baseline_model, acute_ga)
  mom <- acute_y - expected
  keep <- is.finite(mom) & is.finite(acute_day)
  mom <- mom[keep]
  day <- acute_day[keep]
  k_eff <- effective_k(day, min(k, 6L))
  model_data <- data.frame(mom = mom, day = day)
  full <- mgcv::gam(mom ~ s(day, bs = "cr", k = k_eff), method = "REML", data = model_data)
  null <- stats::lm(mom ~ 1, data = model_data)
  observed <- max(0, 2 * (as.numeric(stats::logLik(full)) - as.numeric(stats::logLik(null))))

  set.seed(seed)
  null_stats <- numeric(n_perm)
  for (permutation in seq_len(n_perm)) {
    permuted_data <- transform(model_data, mom = sample(mom))
    null_stats[permutation] <- tryCatch({
      perm_full <- mgcv::gam(mom ~ s(day, bs = "cr", k = k_eff), method = "REML", data = permuted_data)
      perm_null <- stats::lm(mom ~ 1, data = permuted_data)
      max(0, 2 * (as.numeric(stats::logLik(perm_full)) - as.numeric(stats::logLik(perm_null))))
    }, error = function(error) NA_real_)
  }
  valid_null <- null_stats[is.finite(null_stats)]
  list(
    statistic = observed,
    pvalue = if (length(valid_null)) mean(valid_null >= observed) else NA_real_,
    mom = mom,
    day = day,
    keep = which(keep)
  )
}

# Apply the acute MoM likelihood-ratio permutation test to every protein.
run_acute_analysis <- function(samples, proteins, dose = c("V1", "V2"), n_perm = 1000L,
                               k = 8L, cores = 1L, seed = 1010L, cache_dir = NULL) {
  dose <- match.arg(dose)
  days_column <- if (dose == "V1") "diff_collec_vax1" else "diff_collec_vax2"
  baseline <- samples[samples$Timepoint_v1 == "V0", , drop = FALSE]
  acute <- samples[
    samples$Timepoint_v1 == dose & samples[[days_column]] >= 1 & samples[[days_column]] <= 7,
    , drop = FALSE
  ]
  if (!is.null(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  worker <- function(index) {
    protein <- proteins[index]
    cache_file <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(protein, ".rds")) else NULL
    if (!is.null(cache_file) && file.exists(cache_file)) return(readRDS(cache_file))
    result <- fit_acute_mom(
      baseline[[protein]], baseline$GA_collection_days,
      acute[[protein]], acute$GA_collection_days, acute[[days_column]],
      n_perm = n_perm, k = k, seed = seed + index * 10000L
    )
    result$protein <- protein
    if (!is.null(cache_file)) saveRDS(result, cache_file)
    result
  }
  results <- parallel_lapply(seq_along(proteins), worker, cores)
  stats <- data.frame(
    protein = proteins,
    statistic = vapply(results, `[[`, numeric(1), "statistic"),
    pvalue = vapply(results, `[[`, numeric(1), "pvalue")
  )
  stats$fdr <- stats::p.adjust(stats$pvalue, method = "BH")
  list(
    dose = dose,
    baseline = baseline,
    acute = acute,
    days_column = days_column,
    stats = stats,
    protein_results = stats::setNames(results, proteins)
  )
}
