#!/usr/bin/env Rscript

# AOC Multiverse — Subject-Level: Run Analysis + Visualization
#
# Runs all subject-level multiverse scripts in order:
#   1. Sternberg analysis
#   2. N-back analysis
#   3. Sternberg visualization
#   4. N-back visualization
#
# Does NOT run: prep (MATLAB), model tests.
#
# Run from R: setwd to AOC root, then source("multiverse/subject-level/AOC_multiverse_runALL_subject-level.R")

subj_dir <- if (file.exists("multiverse/subject-level/AOC_multiverse_sternberg_analysis_subject.R")) {
  "multiverse/subject-level"
} else if (file.exists("AOC_multiverse_sternberg_analysis_subject.R")) {
  "."
} else {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
    setwd(script_dir)
    "."
  } else {
    stop("Cannot find scripts. Run from AOC project root or multiverse/subject-level/")
  }
}

message("=== AOC Subject-Level Multiverse: Analysis + Visualization ===\n")

message("--- Sternberg analysis ---")
source(file.path(subj_dir, "AOC_multiverse_sternberg_analysis_subject.R"), local = new.env())

message("\n--- N-back analysis ---")
source(file.path(subj_dir, "AOC_multiverse_nback_analysis_subject.R"), local = new.env())

message("\n--- Sternberg visualization ---")
source(file.path(subj_dir, "AOC_multiverse_sternberg_visualize_subject.R"), local = new.env())

message("\n--- N-back visualization ---")
source(file.path(subj_dir, "AOC_multiverse_nback_visualize_subject.R"), local = new.env())

message("\n=== Subject-level multiverse complete ===")
