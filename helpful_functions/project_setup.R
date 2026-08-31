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


# Locate and read the analysis data ----

# By default, data are read from data/processed/. MOMI_DATA_DIR can point to a
# different folder so participant-level data do not need to be copied into Git.
data_directory <- function(project_directory = getwd()) {
  external_data_directory <- Sys.getenv("MOMI_DATA_DIR", unset = "")

  if (nzchar(external_data_directory)) {
    normalizePath(external_data_directory, mustWork = TRUE)
  } else {
    file.path(project_directory, "data", "processed")
  }
}

# Accept the public file names documented in data/README.md and the original
# local file names used while developing the manuscript analysis.
locate_input_file <- function(directory, accepted_file_names) {
  candidate_paths <- file.path(directory, accepted_file_names)
  files_that_exist <- candidate_paths[file.exists(candidate_paths)]

  if (!length(files_that_exist)) {
    stop(
      "Missing input in ", directory, ". Expected one of: ",
      paste(accepted_file_names, collapse = ", "),
      call. = FALSE
    )
  }

  files_that_exist[[1]]
}

# Read the sample-level proteomics data and the SOMAmer annotation table, then
# check all columns needed by the five figure workflows.
load_analysis_inputs <- function(project_directory = getwd()) {
  input_directory <- data_directory(project_directory)

  sample_data_file <- locate_input_file(
    input_directory,
    c("comb_v0_vax.rds", "analysis_data.rds")
  )
  annotation_file <- locate_input_file(
    input_directory,
    c(
      "aptamer_annotations.rds",
      "MOMI_proteomics_Somalogic_AptInfo.rds",
      "somalogic_protein_info.rds"
    )
  )

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
    study_data = load_analysis_inputs(project_directory),
    output = figure_output_folders(figure_number, project_directory)
  )
}

# Read one curated pathway list. Figure scripts do not know its file location.
read_pathways_shown <- function(figure_number, project_directory = getwd()) {
  readr::read_csv(
    file.path(
      project_directory,
      "pathway_lists",
      sprintf("figure_%d.csv", figure_number)
    ),
    show_col_types = FALSE
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
  save_table(table, file.path(output$tables, paste0(file_name, ".csv")))
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

# Save a data frame as a CSV and create its parent folder if needed.
save_table <- function(table, output_file) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(as.data.frame(table), output_file, na = "")
  invisible(output_file)
}
