#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)

for (number in 1:5) {
  figure_script <- file.path(root, sprintf("Figure_%d", number), sprintf("figure_%d.R", number))
  message("Running ", basename(dirname(figure_script)), "...")
  status <- system2("Rscript", figure_script)
  if (status != 0L) stop("Figure ", number, " failed with status ", status, call. = FALSE)
}
message("All main-figure workflows completed.")
