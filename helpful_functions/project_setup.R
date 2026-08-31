# PROJECT SETUP ----
#
# This is the one central file that knows where inputs and pathway tables live.
# Figure scripts do not contain file names or file-location plumbing.
#
# This file:
#   1. Loads the other helper files.
#   2. Reads and checks the two analysis input files.
#   3. Stores the analysis settings (seed, permutations, cores, and GAM k).
#   4. Creates consistent output folders for Figures 1-5.


# Load all analysis functions ----

source_analysis_functions <- function(project_directory = getwd()) {
  helper_files <- c(
    "manuscript_pathways.R",
    "statistical_models.R",
    "pathway_analysis.R",
    "plotting_functions.R",
    "vaccine_analysis.R"
  )

  for (helper_file in helper_files) {
    source(
      file.path(project_directory, "helpful_functions", helper_file),
      local = globalenv()
    )
  }

  invisible(project_directory)
}


# Analysis settings ----

# Read a positive integer from an environment variable. If the variable was
# not set, use the publication-scale default.
read_integer_setting <- function(variable_name, default, minimum = 1L) {
  value <- Sys.getenv(variable_name, unset = as.character(default))
  value <- suppressWarnings(as.integer(value))

  if (is.na(value) || value < minimum) {
    stop(variable_name, " must be an integer >= ", minimum, call. = FALSE)
  }

  value
}

# Central list of settings used throughout the analysis. These values can be
# changed temporarily from the command line without editing any R code.
analysis_parameters <- function() {
  list(
    seed = read_integer_setting("MOMI_SEED", 1010L),
    lasso_trials = read_integer_setting("MOMI_LASSO_TRIALS", 100L),
    n_perm = read_integer_setting("MOMI_N_PERM", 1000L),
    cores = read_integer_setting("MOMI_CORES", 1L),
    gam_k = read_integer_setting("MOMI_GAM_K", 8L, minimum = 4L)
  )
}


# Locate and read the private analysis inputs ----

# The analysis-ready objects were created locally from clinician-provided
# information and are not distributed with this repository. Their private local
# folder must be supplied through MOMI_INPUT_DIR.
private_input_directory <- function() {
  input_directory <- Sys.getenv("MOMI_INPUT_DIR", unset = "")

  if (!nzchar(input_directory)) {
    stop(
      "Set MOMI_INPUT_DIR to the private folder containing the local analysis objects.",
      call. = FALSE
    )
  }

  normalizePath(input_directory, mustWork = TRUE)
}

# Read the local sample-level proteomics object and SOMAmer annotation object,
# then check every column needed by the five figure workflows.
load_analysis_inputs <- function() {
  input_directory <- private_input_directory()
  sample_data_file <- file.path(input_directory, "comb_v0_vax.rds")
  annotation_file <- file.path(
    input_directory,
    "MOMI_proteomics_Somalogic_AptInfo.rds"
  )

  missing_files <- c(sample_data_file, annotation_file)[
    !file.exists(c(sample_data_file, annotation_file))
  ]
  if (length(missing_files)) {
    stop(
      "Missing private analysis object(s): ",
      paste(basename(missing_files), collapse = ", "),
      call. = FALSE
    )
  }

  samples <- readRDS(sample_data_file)
  annotations <- readRDS(annotation_file)

  required_sample_columns <- c(
    "Sample_Label",
    "Timepoint_v1",
    "GA_collection_days",
    "Trim_Collec",
    "Vax1_Trim",
    "diff_collec_vax1",
    "diff_collec_vax2"
  )
  missing_sample_columns <- setdiff(required_sample_columns, names(samples))
  if (length(missing_sample_columns)) {
    stop(
      "The sample data are missing: ",
      paste(missing_sample_columns, collapse = ", "),
      call. = FALSE
    )
  }

  required_annotation_columns <- c("AptName", "EntrezGeneSymbol", "Single_UniPro")
  missing_annotation_columns <- setdiff(required_annotation_columns, names(annotations))
  if (length(missing_annotation_columns)) {
    stop(
      "The annotation data are missing: ",
      paste(missing_annotation_columns, collapse = ", "),
      call. = FALSE
    )
  }

  protein_columns <- grep("^seq\\.", names(samples), value = TRUE)
  protein_columns <- intersect(protein_columns, annotations$AptName)
  if (!length(protein_columns)) {
    stop("No shared seq.* protein columns were found.", call. = FALSE)
  }

  samples$Timepoint_v1 <- factor(samples$Timepoint_v1, levels = c("V0", "V1", "V2"))
  samples$Trim_Collec <- factor(
    samples$Trim_Collec,
    levels = c("Trim_1st", "Trim_2nd", "Trim_3rd")
  )

  list(
    samples = samples,
    annotations = annotations,
    proteins = protein_columns,
    data_path = sample_data_file,
    annotation_path = annotation_file
  )
}


# Output folders and file helpers ----

# Create one plot folder and one table folder for every manuscript figure.
create_output_directories <- function(project_directory = getwd()) {
  output_directories <- c(
    file.path(project_directory, "results", "cache"),
    unlist(lapply(sprintf("figure_%02d", 1:5), function(figure_name) {
      c(
        file.path(project_directory, "results", "figures", figure_name),
        file.path(project_directory, "results", "tables", figure_name)
      )
    }))
  )

  invisible(lapply(
    output_directories,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
}

# Return the three output folders used by one manuscript figure.
figure_output_folders <- function(figure_number, project_directory = getwd()) {
  figure_name <- sprintf("figure_%02d", figure_number)

  list(
    plots = file.path(project_directory, "results", "figures", figure_name),
    tables = file.path(project_directory, "results", "tables", figure_name),
    cache = file.path(project_directory, "results", "cache", figure_name)
  )
}

# Perform all common setup for one figure. Figure scripts receive readable
# objects named settings, study_data, and output without knowing file locations.
prepare_figure <- function(figure_number, project_directory = getwd()) {
  if (!file.exists(file.path(project_directory, "DESCRIPTION"))) {
    stop(
      "Run this script from the MOMI_EDLOW_proteomics repository folder.",
      call. = FALSE
    )
  }

  source_analysis_functions(project_directory)
  create_output_directories(project_directory)

  settings <- analysis_parameters()
  set.seed(settings$seed)

  list(
    settings = settings,
    study_data = load_analysis_inputs(),
    output = figure_output_folders(figure_number, project_directory)
  )
}

# Build a plot path without exposing folders or file extensions in figure code.
result_plot_file <- function(output, file_name) {
  file.path(output$plots, paste0(file_name, ".pdf"))
}

# Build a cache path without exposing folders or file extensions in figure code.
result_cache_file <- function(output, file_name) {
  file.path(output$cache, paste0(file_name, ".rds"))
}

# Create and return a named subfolder for per-protein intermediate calculations.
result_cache_directory <- function(output, directory_name) {
  cache_directory <- file.path(output$cache, directory_name)
  dir.create(cache_directory, recursive = TRUE, showWarnings = FALSE)
  cache_directory
}

# Reuse a completed calculation when available; otherwise calculate and save it.
use_cached_result <- function(output, file_name, calculate) {
  cache_file <- result_cache_file(output, file_name)

  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  result <- calculate()
  saveRDS(result, cache_file)
  result
}

# Save one manuscript result table. Figure scripts provide only a readable name.
save_results_table <- function(table, output, file_name) {
  save_table(table, file.path(output$tables, paste0(file_name, ".tsv")))
}

# Save one manuscript plot. Figure scripts provide only a readable name and size.
save_results_plot <- function(plot, output, file_name, width, height) {
  save_plot(
    plot,
    result_plot_file(output, file_name),
    width = width,
    height = height
  )
}

# Match SOMAmer identifiers to readable gene symbols for tables and plots.
gene_symbols <- function(aptamers, annotations) {
  symbols <- annotations$EntrezGeneSymbol[match(aptamers, annotations$AptName)]
  missing_symbols <- is.na(symbols) | !nzchar(symbols)
  symbols[missing_symbols] <- aptamers[missing_symbols]
  make.unique(symbols)
}

# Assign collection week to the trimester definitions used in the manuscript.
trimester_from_weeks <- function(weeks) {
  factor(
    ifelse(weeks < 14, "1st", ifelse(weeks < 28, "2nd", "3rd")),
    levels = c("1st", "2nd", "3rd")
  )
}

# Save a data frame as a tab-separated table and create its parent folder.
save_table <- function(table, output_file) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(as.data.frame(table), output_file, na = "")
  invisible(output_file)
}
