#!/usr/bin/env Rscript

# RUN ALL FIVE MANUSCRIPT FIGURES ----
#
# Purpose:
#   Run Figure 1 through Figure 5 in manuscript order. Each figure saves its
#   plots and result tables before the next figure begins.
#
# Run from the repository folder with:
#   Rscript scripts/run_all.R

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the MOMI_EDLOW_proteomics repository folder.", call. = FALSE)
}

figure_scripts <- c(
  "Figure_1/figure_1.R",
  "Figure_2/figure_2.R",
  "Figure_3/figure_3.R",
  "Figure_4/figure_4.R",
  "Figure_5/figure_5.R"
)

for (figure_script in figure_scripts) {
  message("Running ", figure_script, " ...")
  exit_status <- system2("Rscript", figure_script)

  if (exit_status != 0L) {
    stop(figure_script, " failed with exit status ", exit_status, call. = FALSE)
  }
}

message("All five manuscript figure workflows are complete.")
