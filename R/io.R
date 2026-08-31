find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) && dir.exists(file.path(path, "figures"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the repository root from: ", start, call. = FALSE)
    }
    path <- parent
  }
}

source_project_functions <- function(root = find_project_root()) {
  files <- c("io.R", "modeling.R", "gsea.R", "plotting.R", "vaccine.R")
  for (file in files[files != "io.R"]) {
    source(file.path(root, "R", file), local = globalenv())
  }
  invisible(root)
}

env_integer <- function(name, default, minimum = 1L) {
  value <- Sys.getenv(name, unset = as.character(default))
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < minimum) {
    stop(name, " must be an integer >= ", minimum, call. = FALSE)
  }
  parsed
}

analysis_parameters <- function() {
  list(
    seed = env_integer("MOMI_SEED", 1010L, 1L),
    lasso_trials = env_integer("MOMI_LASSO_TRIALS", 100L, 1L),
    n_perm = env_integer("MOMI_N_PERM", 1000L, 1L),
    cores = env_integer("MOMI_CORES", 1L, 1L),
    gam_k = env_integer("MOMI_GAM_K", 8L, 4L)
  )
}

data_directory <- function(root = find_project_root()) {
  configured <- Sys.getenv("MOMI_DATA_DIR", unset = "")
  if (nzchar(configured)) normalizePath(configured, mustWork = TRUE) else file.path(root, "data", "processed")
}

locate_input <- function(directory, candidates) {
  paths <- file.path(directory, candidates)
  present <- paths[file.exists(paths)]
  if (!length(present)) {
    stop(
      "Missing input in ", directory, ". Expected one of: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  present[[1]]
}

load_analysis_inputs <- function(root = find_project_root()) {
  directory <- data_directory(root)
  data_path <- locate_input(directory, c("comb_v0_vax.rds", "analysis_data.rds"))
  annotation_path <- locate_input(
    directory,
    c("aptamer_annotations.rds", "MOMI_proteomics_Somalogic_AptInfo.rds", "somalogic_protein_info.rds")
  )

  samples <- readRDS(data_path)
  annotations <- readRDS(annotation_path)

  required_sample_columns <- c(
    "Sample_Label", "Timepoint_v1", "GA_collection_days", "Trim_Collec",
    "Vax1_Trim", "diff_collec_vax1", "diff_collec_vax2"
  )
  missing_sample_columns <- setdiff(required_sample_columns, names(samples))
  if (length(missing_sample_columns)) {
    stop("comb_v0_vax.rds is missing: ", paste(missing_sample_columns, collapse = ", "), call. = FALSE)
  }

  required_annotation_columns <- c("AptName", "EntrezGeneSymbol", "Single_UniPro")
  missing_annotation_columns <- setdiff(required_annotation_columns, names(annotations))
  if (length(missing_annotation_columns)) {
    stop("aptamer_annotations.rds is missing: ", paste(missing_annotation_columns, collapse = ", "), call. = FALSE)
  }

  proteins <- grep("^seq\\.", names(samples), value = TRUE)
  proteins <- intersect(proteins, annotations$AptName)
  if (!length(proteins)) stop("No shared seq.* protein columns were found.", call. = FALSE)

  samples$Timepoint_v1 <- factor(samples$Timepoint_v1, levels = c("V0", "V1", "V2"))
  samples$Trim_Collec <- factor(samples$Trim_Collec, levels = c("Trim_1st", "Trim_2nd", "Trim_3rd"))

  list(
    samples = samples,
    annotations = annotations,
    proteins = proteins,
    data_path = data_path,
    annotation_path = annotation_path
  )
}

create_output_directories <- function(root = find_project_root()) {
  directories <- c(
    file.path(root, "results", "cache"),
    unlist(lapply(sprintf("figure_%02d", 1:5), function(figure) {
      c(
        file.path(root, "results", "figures", figure),
        file.path(root, "results", "tables", figure)
      )
    }))
  )
  invisible(lapply(directories, dir.create, recursive = TRUE, showWarnings = FALSE))
}

figure_paths <- function(number, root = find_project_root()) {
  figure <- sprintf("figure_%02d", number)
  list(
    plots = file.path(root, "results", "figures", figure),
    tables = file.path(root, "results", "tables", figure),
    cache = file.path(root, "results", "cache", figure)
  )
}

gene_symbols <- function(aptamers, annotations) {
  symbols <- annotations$EntrezGeneSymbol[match(aptamers, annotations$AptName)]
  symbols[is.na(symbols) | !nzchar(symbols)] <- aptamers[is.na(symbols) | !nzchar(symbols)]
  make.unique(symbols)
}

trimester_from_weeks <- function(weeks) {
  factor(
    ifelse(weeks < 14, "1st", ifelse(weeks < 28, "2nd", "3rd")),
    levels = c("1st", "2nd", "3rd")
  )
}

save_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(as.data.frame(x), path, na = "")
  invisible(path)
}
