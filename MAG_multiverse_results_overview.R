# MAG Multiverse — Results Overview Table
#
# Extracts key multiverse statistics from result CSVs into one table.
# Filters: excl. 0_500ms, excl. microsaccades, fooof_method == "welch"
#
# Run: Rscript MAG_multiverse_verify_numbers.R

library(tidyverse)

FOOOF_METHOD <- "welch"
GAZE_MEASURES <- c("gaze_deviation", "gaze_velocity", "scan_path_length", "BCEA")
csv_dir <- Sys.getenv("AOC_MULTIVERSE_DIR",
                      unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/MAG/data/multiverse")

filter_data <- function(df) {
  if ("fooof_method" %in% names(df)) df <- df %>% filter(fooof_method == FOOOF_METHOD)
  if ("latency_ms" %in% names(df)) df <- df %>% filter(latency_ms != "0_500ms")
  if ("gaze_measure" %in% names(df)) df <- df %>% filter(gaze_measure %in% GAZE_MEASURES)
  df
}

add_sig <- function(df) {
  df %>% mutate(
    sig = case_when(
      p.value < 0.05 & estimate > 0 ~ "pos",
      p.value < 0.05 & estimate < 0 ~ "neg",
      TRUE ~ "ns"
    )
  )
}

tbl <- tibble(section = character(), effect = character(), task = character(),
              subset = character(), metric = character(), n = integer())

# Gaze ~ condition
for (task in c("nback", "sternberg")) {
  p <- file.path(csv_dir, paste0("MAG_multiverse_", task, "_condition_gaze_results.csv"))
  if (file.exists(p)) {
    M <- read.csv(p, stringsAsFactors = FALSE) %>% filter_data() %>% add_sig()
    tbl <- bind_rows(tbl,
      tibble(section = "3.5.2 Gaze~cond", effect = "Gaze ~ Condition", task = task, subset = "", metric = "pos", n = sum(M$sig == "pos", na.rm = TRUE)),
      tibble(section = "3.5.2 Gaze~cond", effect = "Gaze ~ Condition", task = task, subset = "", metric = "neg", n = sum(M$sig == "neg", na.rm = TRUE)),
      tibble(section = "3.5.2 Gaze~cond", effect = "Gaze ~ Condition", task = task, subset = "", metric = "total", n = nrow(M))
    )
  }
}

# Alpha ~ condition
for (task in c("nback", "sternberg")) {
  p <- file.path(csv_dir, paste0("MAG_multiverse_", task, "_condition_results.csv"))
  if (file.exists(p)) {
    M <- read.csv(p, stringsAsFactors = FALSE) %>% filter_data() %>% add_sig()
    Mp <- M %>% filter(fooof == "FOOOFed")
    Mn <- M %>% filter(fooof == "nonFOOOFed")
    tbl <- bind_rows(tbl,
      tibble(section = "3.5.3 Alpha~cond", effect = "Alpha ~ Condition", task = task, subset = "param", metric = "neg", n = sum(Mp$sig == "neg", na.rm = TRUE)),
      tibble(section = "3.5.3 Alpha~cond", effect = "Alpha ~ Condition", task = task, subset = "param", metric = "total", n = nrow(Mp)),
      tibble(section = "3.5.3 Alpha~cond", effect = "Alpha ~ Condition", task = task, subset = "nonparam", metric = "neg", n = sum(Mn$sig == "neg", na.rm = TRUE)),
      tibble(section = "3.5.3 Alpha~cond", effect = "Alpha ~ Condition", task = task, subset = "nonparam", metric = "pos", n = sum(Mn$sig == "pos", na.rm = TRUE)),
      tibble(section = "3.5.3 Alpha~cond", effect = "Alpha ~ Condition", task = task, subset = "nonparam", metric = "total", n = nrow(Mn))
    )
  }
}

# Gaze -> Alpha (highest condition)
for (task in c("nback", "sternberg")) {
  p <- file.path(csv_dir, paste0("MAG_multiverse_", task, "_conditions_results.csv"))
  if (file.exists(p)) {
    M <- read.csv(p, stringsAsFactors = FALSE) %>% filter_data() %>% add_sig()
    hl <- unique(M$cond_label)[which.max(as.numeric(gsub("[^0-9]", "", unique(M$cond_label))))]
    Mh <- M %>% filter(cond_label == hl)
    Mp <- Mh %>% filter(fooof == "FOOOFed")
    Mn <- Mh %>% filter(fooof == "nonFOOOFed")
    tbl <- bind_rows(tbl,
      tibble(section = "3.5.4 Gaze->Alpha", effect = "Gaze -> Alpha", task = task, subset = "param", metric = "sig", n = sum(Mp$sig %in% c("pos","neg"), na.rm = TRUE)),
      tibble(section = "3.5.4 Gaze->Alpha", effect = "Gaze -> Alpha", task = task, subset = "param", metric = "total", n = nrow(Mp)),
      tibble(section = "3.5.4 Gaze->Alpha", effect = "Gaze -> Alpha", task = task, subset = "nonparam", metric = "sig", n = sum(Mn$sig %in% c("pos","neg"), na.rm = TRUE)),
      tibble(section = "3.5.4 Gaze->Alpha", effect = "Gaze -> Alpha", task = task, subset = "nonparam", metric = "total", n = nrow(Mn)),
      tibble(section = "3.5.4 Gaze->Alpha", effect = "Gaze -> Alpha", task = task, subset = "nonparam", metric = "pos", n = sum(Mn$sig == "pos", na.rm = TRUE)),
      tibble(section = "3.5.4 Gaze->Alpha", effect = "Gaze -> Alpha", task = task, subset = "nonparam", metric = "neg", n = sum(Mn$sig == "neg", na.rm = TRUE))
    )
  }
}

# Interaction
for (task in c("nback", "sternberg")) {
  p <- file.path(csv_dir, paste0("MAG_multiverse_", task, "_interaction_results.csv"))
  if (file.exists(p)) {
    M <- read.csv(p, stringsAsFactors = FALSE) %>% filter_data() %>% add_sig()
    tbl <- bind_rows(tbl,
      tibble(section = "3.5.4 Interaction", effect = "Gaze×Cond non-sig", task = task, subset = "", metric = "ns", n = sum(M$sig == "ns", na.rm = TRUE)),
      tibble(section = "3.5.4 Interaction", effect = "Gaze×Cond non-sig", task = task, subset = "", metric = "total", n = nrow(M))
    )
  }
}

# Aperiodic
ap_files <- c("MAG_multiverse_nback_aperiodic_gaze_results.csv", "MAG_multiverse_sternberg_aperiodic_gaze_results.csv")
if (all(file.exists(file.path(csv_dir, ap_files)))) {
  M_ap <- bind_rows(lapply(file.path(csv_dir, ap_files), read.csv, stringsAsFactors = FALSE)) %>% filter_data() %>% add_sig()
  tbl <- bind_rows(tbl,
    tibble(section = "3.5.5 Aperiodic", effect = "Gaze -> Offset", task = "both", subset = "", metric = "sig_pos", n = sum(M_ap$aperiodic_measure == "Offset" & M_ap$sig == "pos", na.rm = TRUE)),
    tibble(section = "3.5.5 Aperiodic", effect = "Gaze -> Offset", task = "both", subset = "", metric = "total", n = sum(M_ap$aperiodic_measure == "Offset", na.rm = TRUE)),
    tibble(section = "3.5.5 Aperiodic", effect = "Gaze -> Exponent", task = "both", subset = "", metric = "sig_pos", n = sum(M_ap$aperiodic_measure == "Exponent" & M_ap$sig == "pos", na.rm = TRUE)),
    tibble(section = "3.5.5 Aperiodic", effect = "Gaze -> Exponent", task = "both", subset = "", metric = "total", n = sum(M_ap$aperiodic_measure == "Exponent", na.rm = TRUE))
  )
}

# Output
message("=== MAG Multiverse Results Overview (excl. 0_500ms, microsaccades) ===\n")
print(tbl, n = Inf)
write.csv(tbl, file.path(csv_dir, "MAG_multiverse_results_overview.csv"), row.names = FALSE)
message(sprintf("\nSaved: %s", file.path(csv_dir, "MAG_multiverse_results_overview.csv")))
