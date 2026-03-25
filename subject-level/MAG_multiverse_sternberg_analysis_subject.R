# MAG Multiverse — Sternberg: Subject-Level Analysis (robustness)
# Loads pre-aggregated subject-level CSV (from MAG_multiverse_prep_subject.m),
# fits models per universe as an aggregation sensitivity analysis.
# Uses multiverse R package (Sarma et al., 2021).
#
# Models:
#   M (main, 7D):      alpha ~ gaze_value * Condition + (1|subjectID)
#                       alpha ~ gaze_value + (1|subjectID) [per condition]
#   M_eeg (EEG, 5D):   alpha ~ Condition + (1|subjectID)
#   M_gaze (gaze, 3D): gaze_value ~ Condition + (1|subjectID)
#
# Output CSVs (saved to csv_dir, with _subject_ prefix):
#   MAG_multiverse_sternberg_subject_results.csv
#   MAG_multiverse_sternberg_subject_conditions_results.csv
#   MAG_multiverse_sternberg_subject_condition_results.csv
#   MAG_multiverse_sternberg_subject_interaction_results.csv
#   MAG_multiverse_sternberg_subject_condition_gaze_results.csv
#
# Note:
#   This script is intended as a robustness analysis complementing the
#   primary trial-level models. Trial-level variance is already aggregated,
#   so adding "trial" random effects here is not appropriate.

library(lme4)
library(lmerTest)
library(broom.mixed)
library(tidyverse)
library(multiverse)

# ========== LOGGING MODE ==========
# Set AOC_VERBOSE_LOGGING=true to show full lmer warnings/messages.
VERBOSE_LOGGING <- tolower(Sys.getenv("AOC_VERBOSE_LOGGING", unset = "false")) %in% c("1", "true", "yes", "y")
MIN_N_PER_UNIVERSE <- 10

# ========== PATHS ==========
csv_dir  <- Sys.getenv("AOC_MULTIVERSE_DIR",
                        unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/MAG/data/multiverse")
csv_path <- file.path(csv_dir, "MAG_multiverse_sternberg_subject.csv")
if (!file.exists(csv_path)) stop("Subject-level CSV not found: ", csv_path,
                                  "\nRun MAG_multiverse_prep_subject.m first.")

# ========== LOAD & FILTER DATA ==========
dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
dat$subjectID <- as.factor(dat$subjectID)
dat$Condition <- as.factor(dat$Condition)
if (!("n_trials_subject_condition" %in% names(dat))) {
  dat$n_trials_subject_condition <- 1
}
dat$n_trials_subject_condition <- suppressWarnings(as.numeric(dat$n_trials_subject_condition))
dat$n_trials_subject_condition[!is.finite(dat$n_trials_subject_condition) | dat$n_trials_subject_condition <= 0] <- 1

message(sprintf("Loaded: %d subject-level rows, %d subjects.",
                nrow(dat), n_distinct(dat$subjectID)))

n_rows_0500 <- sum(dat$latency_ms == "0_500ms", na.rm = TRUE)
dat <- dat %>% filter(latency_ms != "0_500ms")
message(sprintf("Excluded 0_500ms rows: %d. Remaining: %d rows.",
                n_rows_0500, nrow(dat)))

if ("gaze_measure" %in% names(dat)) {
  dat <- dat %>% filter(gaze_measure != "gaze_density" | is.na(gaze_measure))
}
dat <- dat %>% filter(!electrodes %in% c("all", "parietal"))
dat <- dat %>% filter(baseline_eeg != "dB", baseline_gaze != "dB")
message(sprintf("After filtering: %d rows.", nrow(dat)))

# ========== HELPER FUNCTIONS ==========
robust_z <- function(x, winsor = 3) {
  med <- median(x, na.rm = TRUE)
  ma  <- mad(x, na.rm = TRUE)
  if (is.na(ma) || ma == 0) return(rep(NaN, length(x)))
  z <- (x - med) / ma
  z[z >  winsor] <-  winsor
  z[z < -winsor] <- -winsor
  z
}

add_sig <- function(df) {
  df %>% mutate(
    condition = factor(case_when(
      p.value < 0.05 & estimate > 0 ~ "Positive",
      p.value < 0.05 & estimate < 0 ~ "Negative",
      TRUE ~ "Non-significant"
    ), levels = c("Positive", "Negative", "Non-significant"))
  )
}

safe_extract <- function(results, var_name) {
  r <- results[[var_name]]
  if (is.null(r) || !is.data.frame(r) || nrow(r) == 0) return(NULL)
  r
}

run_lm <- function(formula_obj, data_obj) {
  use_weights <- "model_weight" %in% names(data_obj) && any(is.finite(data_obj$model_weight))
  lm_args <- list(formula = formula_obj, data = data_obj)
  if (use_weights) lm_args$weights <- data_obj$model_weight
  tryCatch(
    suppressWarnings(suppressMessages(do.call(lm, lm_args))),
    error = function(e) NULL
  )
}

run_lmer <- function(formula_obj, data_obj) {
  use_weights <- "model_weight" %in% names(data_obj) && any(is.finite(data_obj$model_weight))
  lmer_args <- list(
    formula = formula_obj,
    data = data_obj,
    control = lmerControl(optimizer = "bobyqa")
  )
  if (use_weights) lmer_args$weights <- data_obj$model_weight

  if (VERBOSE_LOGGING) {
    return(tryCatch(
      do.call(lmer, lmer_args),
      error = function(e) NULL
    ))
  }
  tryCatch(
    suppressWarnings(
      suppressMessages(
        do.call(lmer, lmer_args)
      )
    ),
    error = function(e) NULL
  )
}

prepare_model_data <- function(df, cols) {
  out <- df %>% filter(complete.cases(across(all_of(cols))))
  if (!("model_weight" %in% names(out))) out$model_weight <- 1
  out$model_weight <- suppressWarnings(as.numeric(out$model_weight))
  out$model_weight[!is.finite(out$model_weight) | out$model_weight <= 0] <- 1
  out
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

summarize_robustness <- function(df, model_family) {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  grouping_cols <- intersect(c("term", "cond_label", "aperiodic_measure"), names(df))
  if (!(".universe" %in% names(df)) || !("estimate" %in% names(df))) return(tibble())

  out <- df %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(
      n_rows = n(),
      n_universes = n_distinct(.universe),
      median_beta = median(estimate, na.rm = TRUE),
      pct_beta_positive = 100 * mean(estimate > 0, na.rm = TRUE),
      pct_beta_negative = 100 * mean(estimate < 0, na.rm = TRUE),
      dominant_sign = case_when(
        pct_beta_positive > pct_beta_negative ~ "Positive",
        pct_beta_negative > pct_beta_positive ~ "Negative",
        TRUE ~ "Mixed"
      ),
      common_overlap_low = safe_max(conf.low),
      common_overlap_high = safe_min(conf.high),
      has_common_overlap = common_overlap_high >= common_overlap_low,
      pct_cis_crossing_zero = 100 * mean(conf.low <= 0 & conf.high >= 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(model_family = model_family, .before = 1)

  out
}

# Weights: trial count per subject-condition (no trial-level variance merge)
branch_cols <- c("electrodes", "fooof", "latency_ms", "alpha_type",
                 "gaze_measure", "baseline_eeg", "baseline_gaze")

# Condition setup
cond_levels  <- sort(unique(dat$Condition))
cond_labels  <- setNames(paste0("Set size ", cond_levels), as.character(cond_levels))
highest_cond <- max(as.numeric(as.character(cond_levels)))
analysis_role <- "subject_level_robustness"

# ========== MAIN MULTIVERSE (7 dimensions) ==========
message("Setting up main multiverse (7 dimensions)...")

M <- multiverse()

inside(M, {
  .elec   <- branch(electrodes,    "posterior", "occipital")
  .fooof  <- branch(fooof,         "FOOOFed", "nonFOOOFed")
  .lat    <- branch(latency_ms,    "0_1000ms", "0_2000ms", "1000_2000ms")
  .alpha  <- branch(alpha_type,    "canonical", "IAF")
  .gaze   <- branch(gaze_measure,  "scan_path_length", "gaze_velocity", "microsaccades", "BCEA", "gaze_deviation")
  .bleeg  <- branch(baseline_eeg,  "raw", "pct_change")
  .blgaze <- branch(baseline_gaze, "raw", "pct_change")

  df <- dat %>%
    filter(electrodes == .elec, fooof == .fooof, latency_ms == .lat,
           alpha_type == .alpha, gaze_measure == .gaze,
           baseline_eeg == .bleeg, baseline_gaze == .blgaze) %>%
    filter(complete.cases(alpha, gaze_value, Condition, subjectID))

  df$model_weight <- df$n_trials_subject_condition
  df$alpha        <- robust_z(df$alpha)
  df$gaze_value   <- robust_z(df$gaze_value)
  df <- prepare_model_data(df, c("alpha", "gaze_value", "Condition", "subjectID", "model_weight"))
  valid <- nrow(df) >= MIN_N_PER_UNIVERSE && !any(is.nan(df$gaze_value)) && !any(is.nan(df$alpha))

  tid_int <- if (valid) {
    fit <- run_lmer(alpha ~ gaze_value * Condition + (1 | subjectID), df)
    if (!is.null(fit)) broom.mixed::tidy(fit, conf.int = TRUE) else tibble()
  } else tibble()

  tid_cond <- if (valid) {
    bind_rows(lapply(as.character(cond_levels), function(cl) {
      dc <- df[df$Condition == cl, ]
      dc <- prepare_model_data(dc, c("alpha", "gaze_value", "subjectID", "model_weight"))
      if (nrow(dc) < MIN_N_PER_UNIVERSE) return(tibble())
      fit_c <- run_lm(alpha ~ gaze_value, dc)
      if (is.null(fit_c)) return(tibble())
      tid_c <- broom::tidy(fit_c, conf.int = TRUE) %>% filter(term == "gaze_value")
      tid_c$cond_label <- cond_labels[cl]
      tid_c
    }))
  } else tibble()
})

message("Executing main multiverse...")
execute_multiverse(M)

# ========== EXTRACT RESULTS ==========
message("Extracting results from main multiverse...")
M_expanded <- expand(M)

M_int <- M_expanded %>%
  mutate(tid = map(.results, safe_extract, "tid_int")) %>%
  filter(!map_lgl(tid, is.null)) %>%
  select(.universe, all_of(branch_cols), tid) %>%
  unnest(tid) %>%
  filter(grepl("gaze_value", term)) %>%
  mutate(analysis_role = analysis_role)

M_cond <- M_expanded %>%
  mutate(tid = map(.results, safe_extract, "tid_cond")) %>%
  filter(!map_lgl(tid, is.null)) %>%
  select(.universe, all_of(branch_cols), tid) %>%
  unnest(tid) %>%
  mutate(analysis_role = analysis_role)

if (nrow(M_cond) == 0L) stop("No successful LMM fits.")

M_cond <- add_sig(M_cond)
M_int  <- add_sig(M_int)

write.csv(M_int, file.path(csv_dir, "MAG_multiverse_sternberg_subject_results.csv"), row.names = FALSE)
write.csv(M_cond, file.path(csv_dir, "MAG_multiverse_sternberg_subject_conditions_results.csv"), row.names = FALSE)
message("Saved: MAG_multiverse_sternberg_subject_results.csv, MAG_multiverse_sternberg_subject_conditions_results.csv")

# Extract interaction term
highest_int_term <- paste0("gaze_value:Condition", highest_cond)
M_interaction    <- M_int %>% filter(term == highest_int_term)
if (nrow(M_interaction) > 0) {
  write.csv(M_interaction, file.path(csv_dir, "MAG_multiverse_sternberg_subject_interaction_results.csv"), row.names = FALSE)
  message(sprintf("Saved: MAG_multiverse_sternberg_subject_interaction_results.csv (%d universes)", nrow(M_interaction)))
} else {
  message("WARNING: No interaction terms found for ", highest_int_term)
}

# ========== EEG-ONLY MULTIVERSE (5 dimensions) ==========
message("Setting up EEG-only multiverse (5 dimensions)...")

dat_eeg <- dat %>%
  filter(!electrodes %in% c("all", "parietal"), baseline_eeg != "dB") %>%
  select(subjectID, Condition, alpha, n_trials_subject_condition, electrodes, fooof, latency_ms, alpha_type, baseline_eeg) %>%
  distinct()

highest_alpha_term <- paste0("Condition", highest_cond)

M_eeg <- multiverse()

inside(M_eeg, {
  .elec  <- branch(electrodes,   "posterior", "occipital")
  .fooof <- branch(fooof,        "FOOOFed", "nonFOOOFed")
  .lat   <- branch(latency_ms,   "0_1000ms", "0_2000ms", "1000_2000ms")
  .alpha <- branch(alpha_type,   "canonical", "IAF")
  .bleeg <- branch(baseline_eeg, "raw", "pct_change")

  de <- dat_eeg %>%
    filter(electrodes == .elec, fooof == .fooof, latency_ms == .lat,
           alpha_type == .alpha, baseline_eeg == .bleeg) %>%
    filter(complete.cases(alpha, Condition, subjectID))

  de$model_weight <- de$n_trials_subject_condition
  de$alpha <- robust_z(de$alpha)
  de <- prepare_model_data(de, c("alpha", "Condition", "subjectID", "model_weight"))
  valid <- nrow(de) >= MIN_N_PER_UNIVERSE && !any(is.nan(de$alpha))

  tid_ca <- if (valid) {
    fit <- run_lmer(alpha ~ Condition + (1 | subjectID), de)
    if (!is.null(fit)) {
      broom.mixed::tidy(fit, conf.int = TRUE) %>% filter(term == highest_alpha_term)
    } else tibble()
  } else tibble()
})

message("Executing EEG-only multiverse...")
execute_multiverse(M_eeg)

eeg_branch_cols  <- c("electrodes", "fooof", "latency_ms", "alpha_type", "baseline_eeg")
M_eeg_expanded   <- expand(M_eeg)

M_ca <- M_eeg_expanded %>%
  mutate(tid = map(.results, safe_extract, "tid_ca")) %>%
  filter(!map_lgl(tid, is.null)) %>%
  select(.universe, all_of(eeg_branch_cols), tid) %>%
  unnest(tid) %>%
  mutate(analysis_role = analysis_role)

if (nrow(M_ca) > 0L) {
  M_ca <- add_sig(M_ca)
  write.csv(M_ca, file.path(csv_dir, "MAG_multiverse_sternberg_subject_condition_results.csv"), row.names = FALSE)
  message(sprintf("Saved: MAG_multiverse_sternberg_subject_condition_results.csv (%d universes)", nrow(M_ca)))
} else {
  message("WARNING: No successful condition \u2192 alpha fits.")
}

# ========== GAZE-ONLY MULTIVERSE (3 dimensions) ==========
message("Setting up gaze-only multiverse (3 dimensions)...")

dat_gaze <- dat %>%
  select(subjectID, Condition, gaze_value, n_trials_subject_condition, latency_ms, gaze_measure, baseline_gaze) %>%
  distinct()

highest_gaze_term <- paste0("Condition", highest_cond)

M_gaze <- multiverse()

inside(M_gaze, {
  .lat    <- branch(latency_ms,    "0_1000ms", "0_2000ms", "1000_2000ms")
  .gaze   <- branch(gaze_measure,  "scan_path_length", "gaze_velocity", "microsaccades", "BCEA", "gaze_deviation")
  .blgaze <- branch(baseline_gaze, "raw", "pct_change")

  dg <- dat_gaze %>%
    filter(latency_ms == .lat, gaze_measure == .gaze, baseline_gaze == .blgaze) %>%
    filter(complete.cases(gaze_value, Condition, subjectID))

  dg$model_weight <- dg$n_trials_subject_condition
  dg$gaze_value   <- robust_z(dg$gaze_value)
  dg <- prepare_model_data(dg, c("gaze_value", "Condition", "subjectID", "model_weight"))
  valid <- nrow(dg) >= MIN_N_PER_UNIVERSE && !any(is.nan(dg$gaze_value))

  tid_cg <- if (valid) {
    fit <- run_lmer(gaze_value ~ Condition + (1 | subjectID), dg)
    if (!is.null(fit)) {
      broom.mixed::tidy(fit, conf.int = TRUE) %>% filter(term == highest_gaze_term)
    } else tibble()
  } else tibble()
})

message("Executing gaze-only multiverse...")
execute_multiverse(M_gaze)

gaze_branch_cols <- c("latency_ms", "gaze_measure", "baseline_gaze")
M_gaze_expanded  <- expand(M_gaze)

M_cg <- M_gaze_expanded %>%
  mutate(tid = map(.results, safe_extract, "tid_cg")) %>%
  filter(!map_lgl(tid, is.null)) %>%
  select(.universe, all_of(gaze_branch_cols), tid) %>%
  unnest(tid) %>%
  mutate(analysis_role = analysis_role)

if (nrow(M_cg) > 0L) {
  M_cg <- add_sig(M_cg)
  write.csv(M_cg, file.path(csv_dir, "MAG_multiverse_sternberg_subject_condition_gaze_results.csv"), row.names = FALSE)
  message(sprintf("Saved: MAG_multiverse_sternberg_subject_condition_gaze_results.csv (%d universes)", nrow(M_cg)))
} else {
  message("WARNING: No successful condition \u2192 gaze fits.")
}

# ========== APERIODIC MULTIVERSE ==========
M_ap_gaze_results <- tibble()
M_ap_cond_results <- tibble()
if ("aperiodic_offset" %in% names(dat) && "aperiodic_exponent" %in% names(dat)) {
  message("Setting up aperiodic multiverse (subject-level)...")

  dat_ap <- dat %>%
    filter(fooof == "FOOOFed") %>%
    filter(complete.cases(aperiodic_offset, aperiodic_exponent))

  dat_ap_gaze <- dat_ap %>%
    select(subjectID, Condition, aperiodic_offset, aperiodic_exponent,
           gaze_value, n_trials_subject_condition,
           electrodes, latency_ms, gaze_measure, baseline_gaze) %>%
    distinct()

  M_ap_gaze <- multiverse()

  inside(M_ap_gaze, {
    .elec   <- branch(electrodes,    "posterior", "occipital")
    .lat    <- branch(latency_ms,    "0_1000ms", "0_2000ms", "1000_2000ms")
    .gaze   <- branch(gaze_measure,  "scan_path_length", "gaze_velocity", "microsaccades", "BCEA", "gaze_deviation")
    .blgaze <- branch(baseline_gaze, "raw", "pct_change")

    dap <- dat_ap_gaze %>%
      filter(electrodes == .elec, latency_ms == .lat,
             gaze_measure == .gaze, baseline_gaze == .blgaze) %>%
      filter(complete.cases(aperiodic_exponent, aperiodic_offset, gaze_value, Condition, subjectID))

    dap$aperiodic_exponent <- robust_z(dap$aperiodic_exponent)
    dap$aperiodic_offset <- robust_z(dap$aperiodic_offset)
    dap$gaze_value   <- robust_z(dap$gaze_value)
    dap$model_weight <- dap$n_trials_subject_condition
    dap <- prepare_model_data(dap, c("aperiodic_exponent", "aperiodic_offset", "gaze_value", "Condition", "subjectID", "model_weight"))
    valid <- nrow(dap) >= MIN_N_PER_UNIVERSE && !any(is.nan(dap$gaze_value)) &&
             !any(is.nan(dap$aperiodic_exponent)) && !any(is.nan(dap$aperiodic_offset))

    tid_exp_gaze <- if (valid) {
      fit <- run_lmer(aperiodic_exponent ~ gaze_value + (1 | subjectID), dap)
      if (!is.null(fit)) {
        broom.mixed::tidy(fit, conf.int = TRUE) %>% filter(term == "gaze_value") %>%
          mutate(aperiodic_measure = "Exponent")
      } else tibble()
    } else tibble()

    tid_off_gaze <- if (valid) {
      fit <- run_lmer(aperiodic_offset ~ gaze_value + (1 | subjectID), dap)
      if (!is.null(fit)) {
        broom.mixed::tidy(fit, conf.int = TRUE) %>% filter(term == "gaze_value") %>%
          mutate(aperiodic_measure = "Offset")
      } else tibble()
    } else tibble()
  })

  message("Executing aperiodic ~ gaze multiverse...")
  execute_multiverse(M_ap_gaze)

  ap_gaze_branch_cols <- c("electrodes", "latency_ms", "gaze_measure", "baseline_gaze")
  M_ap_gaze_exp <- expand(M_ap_gaze)

  M_ap_gaze_results <- bind_rows(
    M_ap_gaze_exp %>%
      mutate(tid = map(.results, safe_extract, "tid_exp_gaze")) %>%
      filter(!map_lgl(tid, is.null)) %>%
      select(.universe, all_of(ap_gaze_branch_cols), tid) %>%
      unnest(tid),
    M_ap_gaze_exp %>%
      mutate(tid = map(.results, safe_extract, "tid_off_gaze")) %>%
      filter(!map_lgl(tid, is.null)) %>%
      select(.universe, all_of(ap_gaze_branch_cols), tid) %>%
      unnest(tid)
  ) %>%
    mutate(analysis_role = analysis_role)

  if (nrow(M_ap_gaze_results) > 0) {
    M_ap_gaze_results <- add_sig(M_ap_gaze_results)
    write.csv(M_ap_gaze_results, file.path(csv_dir, "MAG_multiverse_sternberg_subject_aperiodic_gaze_results.csv"), row.names = FALSE)
    message(sprintf("Saved: aperiodic_gaze_results.csv (%d rows)", nrow(M_ap_gaze_results)))
  }

  dat_ap_eeg <- dat_ap %>%
    select(subjectID, Condition, aperiodic_offset, aperiodic_exponent, n_trials_subject_condition, electrodes, latency_ms) %>%
    distinct()

  M_ap_cond <- multiverse()

  inside(M_ap_cond, {
    .elec <- branch(electrodes, "posterior", "occipital")
    .lat  <- branch(latency_ms, "0_1000ms", "0_2000ms", "1000_2000ms")

    dap <- dat_ap_eeg %>%
      filter(electrodes == .elec, latency_ms == .lat) %>%
      filter(complete.cases(aperiodic_exponent, aperiodic_offset, Condition, subjectID))

    dap$aperiodic_exponent <- robust_z(dap$aperiodic_exponent)
    dap$aperiodic_offset <- robust_z(dap$aperiodic_offset)
    dap$model_weight <- dap$n_trials_subject_condition
    dap <- prepare_model_data(dap, c("aperiodic_exponent", "aperiodic_offset", "Condition", "subjectID", "model_weight"))
    valid <- nrow(dap) >= MIN_N_PER_UNIVERSE && !any(is.nan(dap$aperiodic_exponent)) && !any(is.nan(dap$aperiodic_offset))

    tid_exp_cond <- if (valid) {
      fit <- run_lmer(aperiodic_exponent ~ Condition + (1 | subjectID), dap)
      if (!is.null(fit)) {
        broom.mixed::tidy(fit, conf.int = TRUE) %>%
          filter(term == highest_alpha_term) %>%
          mutate(aperiodic_measure = "Exponent")
      } else tibble()
    } else tibble()

    tid_off_cond <- if (valid) {
      fit <- run_lmer(aperiodic_offset ~ Condition + (1 | subjectID), dap)
      if (!is.null(fit)) {
        broom.mixed::tidy(fit, conf.int = TRUE) %>%
          filter(term == highest_alpha_term) %>%
          mutate(aperiodic_measure = "Offset")
      } else tibble()
    } else tibble()
  })

  message("Executing aperiodic ~ condition multiverse...")
  execute_multiverse(M_ap_cond)

  ap_cond_branch_cols <- c("electrodes", "latency_ms")
  M_ap_cond_exp <- expand(M_ap_cond)

  M_ap_cond_results <- bind_rows(
    M_ap_cond_exp %>%
      mutate(tid = map(.results, safe_extract, "tid_exp_cond")) %>%
      filter(!map_lgl(tid, is.null)) %>%
      select(.universe, all_of(ap_cond_branch_cols), tid) %>%
      unnest(tid),
    M_ap_cond_exp %>%
      mutate(tid = map(.results, safe_extract, "tid_off_cond")) %>%
      filter(!map_lgl(tid, is.null)) %>%
      select(.universe, all_of(ap_cond_branch_cols), tid) %>%
      unnest(tid)
  ) %>%
    mutate(analysis_role = analysis_role)

  if (nrow(M_ap_cond_results) > 0) {
    M_ap_cond_results <- add_sig(M_ap_cond_results)
    write.csv(M_ap_cond_results, file.path(csv_dir, "MAG_multiverse_sternberg_subject_aperiodic_condition_results.csv"), row.names = FALSE)
    message(sprintf("Saved: aperiodic_condition_results.csv (%d rows)", nrow(M_ap_cond_results)))
  }

} else {
  message("Skipping aperiodic multiverse: columns not found in CSV.")
}

# ========== ROBUSTNESS SUMMARY OUTPUTS (CSV ONLY) ==========
robust_main <- summarize_robustness(M_int, "main_model_terms")
robust_cond <- summarize_robustness(M_cond, "per_condition_gaze_slopes")
robust_interaction <- summarize_robustness(M_interaction, "highest_condition_interaction")
robust_ca <- summarize_robustness(M_ca, "condition_to_alpha")
robust_cg <- summarize_robustness(M_cg, "condition_to_gaze")
robust_ap_gaze <- summarize_robustness(M_ap_gaze_results, "aperiodic_to_gaze")
robust_ap_cond <- summarize_robustness(M_ap_cond_results, "aperiodic_to_condition")

if (nrow(robust_main) > 0) write.csv(robust_main, file.path(csv_dir, "MAG_multiverse_sternberg_subject_results_robustness_summary.csv"), row.names = FALSE)
if (nrow(robust_cond) > 0) write.csv(robust_cond, file.path(csv_dir, "MAG_multiverse_sternberg_subject_conditions_results_robustness_summary.csv"), row.names = FALSE)
if (nrow(robust_interaction) > 0) write.csv(robust_interaction, file.path(csv_dir, "MAG_multiverse_sternberg_subject_interaction_results_robustness_summary.csv"), row.names = FALSE)
if (nrow(robust_ca) > 0) write.csv(robust_ca, file.path(csv_dir, "MAG_multiverse_sternberg_subject_condition_results_robustness_summary.csv"), row.names = FALSE)
if (nrow(robust_cg) > 0) write.csv(robust_cg, file.path(csv_dir, "MAG_multiverse_sternberg_subject_condition_gaze_results_robustness_summary.csv"), row.names = FALSE)
if (nrow(robust_ap_gaze) > 0) write.csv(robust_ap_gaze, file.path(csv_dir, "MAG_multiverse_sternberg_subject_aperiodic_gaze_results_robustness_summary.csv"), row.names = FALSE)
if (nrow(robust_ap_cond) > 0) write.csv(robust_ap_cond, file.path(csv_dir, "MAG_multiverse_sternberg_subject_aperiodic_condition_results_robustness_summary.csv"), row.names = FALSE)

robustness_summary <- bind_rows(
  robust_main, robust_cond, robust_interaction, robust_ca, robust_cg, robust_ap_gaze, robust_ap_cond
)
if (nrow(robustness_summary) > 0) {
  write.csv(robustness_summary, file.path(csv_dir, "MAG_multiverse_sternberg_subject_robustness_summary.csv"), row.names = FALSE)
  message(sprintf("Saved robustness summaries (combined rows: %d)", nrow(robustness_summary)))
}

message("=== Sternberg SUBJECT-LEVEL multiverse ANALYSIS complete ===")
message("Result CSVs saved to: ", csv_dir)
