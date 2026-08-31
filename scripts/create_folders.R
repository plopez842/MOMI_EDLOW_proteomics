#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]), ".."), mustWork = FALSE)
root <- sub("^--file=", "", root)
if (!file.exists(file.path(root, "DESCRIPTION"))) root <- getwd()
source(file.path(root, "R", "io.R"))
root <- find_project_root(root)
create_output_directories(root)
message("Created output folders under: ", file.path(root, "results"))
