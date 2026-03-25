# AOC Multiverse — Sternberg: Analysis
# Fits LMMs for all multiverse universes and saves result CSVs.
# Uses multiverse R package (Sarma et al., 2021).
#
# Models (primary, trial-level):
#   M (main, 7D):      alpha ~ gaze_value * Condition + (1 + Condition || subjectID) [fallback: (1|subjectID)]
#                       alpha ~ gaze_value + (1|subjectID) [per condition]
#   M_eeg (EEG, 5D):   alpha ~ Condition + (1 + Condition || subjectID) [fallback: (1|subjectID)]
#   M_gaze (gaze, 3D): gaze_value ~ Condition + (1 + Condition || subjectID) [fallback: (1|subjectID)]
#
# Output CSVs (saved to csv_dir):
#   MAG_multiverse_sternberg_results.csv               — full interaction model terms
#   MAG_multiverse_sternberg_conditions_results.csv     — per-condition gaze→alpha
#   MAG_multiverse_sternberg_condition_results.csv      — condition→alpha (EEG-only)
#   MAG_multiverse_sternberg_interaction_results.csv    — gaze×condition interaction term
#   MAG_multiverse_sternberg_condition_gaze_results.csv — condition→gaze (gaze-only)

library(lme4)
library(lmerTest)
library(broom.mixed)
library(tidyverse)
library(multiverse)

# ========== LOGGING MODE ==========
# Set AOC_VERBOSE_LOGGING=true to show full lmer warnings/messages.
VERBOSE_LOGGING <- tolower(Sys.getenv("AOC_VERBOSE_LOGGING", unset = "false")) %in% c("1", "true", "yes", "y")
# Allowed: "singleFFT", "both", "welch"
FOOOF_METHOD <- "welch"
R2_THRESHOLD <- 0.90
r2_dir <- Sys.getenv("AOC_R2_DIR",
                     unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/MAG/data/controls/multiverse")

# ========== PATHS ==========
csv_dir  <- Sys.getenv("AOC_MULTIVERSE_DIR",
                        unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/MAG/data/multiverse")
csv_path_method <- file.path(csv_dir, paste0("MAG_multiverse_sternberg_", FOOOF_METHOD, ".csv"))
csv_path <- if (file.exists(csv_path_method)) csv_path_method else file.path(csv_dir, "MAG_multiverse_sternberg.csv")
if (!file.exists(csv_path)) stop("CSV not found: ", csv_path)

# ========== LOAD & FILTER DATA ==========
dat <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
dat$subjectID <- as.factor(dat$subjectID)
dat$Condition <- as.factor(dat$Condition)

message(sprintf("Loaded: %d rows, %d universes.", nrow(dat), n_distinct(dat$universe_id)))
n_rows_0500 <- sum(dat$latency_ms == "0_500ms", na.rm = TRUE)
dat <- dat %>% filter(latency_ms != "0_500ms")
message(sprintf("Excluded 0_500ms rows: %d. Remaining: %d rows, %d universes.",
                n_rows_0500, nrow(dat), n_distinct(dat$universe_id)))

if (!("r_squared" %in% names(dat))) {
  r2_path_method <- file.path(r2_dir, paste0("MAG_fooof_r2_sternberg_", FOOOF_METHOD, ".csv"))
  r2_path <- if (file.exists(r2_path_method)) r2_path_method else file.path(r2_dir, "MAG_fooof_r2_sternberg.csv")

  if (file.exists(r2_path)) {
    r2 <- read.csv(r2_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
    join_cols <- c("subjectID", "Condition", "Trial", "latency_ms", "electrodes")

    if (all(c(join_cols, "r_squared") %in% names(r2)) && all(join_cols %in% names(dat))) {
      dat$subjectID_num <- suppressWarnings(as.numeric(as.character(dat$subjectID)))
      r2$subjectID_num <- suppressWarnings(as.numeric(as.character(r2$subjectID)))
      dat$Condition <- as.factor(dat$Condition)
      r2$Condition <- as.factor(r2$Condition)

      dat <- dat %>%
        left_join(
          r2 %>% select(subjectID_num, Condition, Trial, latency_ms, electrodes, r_squared),
          by = c("subjectID_num", "Condition", "Trial", "latency_ms", "electrodes")
        ) %>%
        select(-subjectID_num)
      message(sprintf("Joined trial-level FOOOF R2 from: %s", r2_path))
    } else {
      warning("Skipping FOOOF R2 join: required join columns missing in multiverse and/or R2 CSV.")
    }
  } else {
    warning("FOOOF R2 CSV not found; no trial-level R2 filter will be applied.")
  }
}

if (!("r_squared" %in% names(dat))) {
  dat$r_squared <- NA_real_
}

n_rows_before_r2 <- nrow(dat)
n_univ_before_r2 <- n_distinct(dat$universe_id)
n_fooof_rows <- sum(dat$fooof == "FOOOFed", na.rm = TRUE)
low_r2_mask <- dat$fooof == "FOOOFed" & is.finite(dat$r_squared) & dat$r_squared < R2_THRESHOLD
n_low_r2_removed <- sum(low_r2_mask, na.rm = TRUE)
n_fooof_no_r2 <- sum(dat$fooof == "FOOOFed" & !is.finite(dat$r_squared), na.rm = TRUE)

latency_exclusion_summary <- dat %>%
  mutate(low_r2_excluded = low_r2_mask) %>%
  group_by(latency_ms) %>%
  summarise(
    fooof_rows_total = sum(fooof == "FOOOFed", na.rm = TRUE),
    rows_removed_low_r2 = sum(low_r2_excluded, na.rm = TRUE),
    fooof_rows_missing_r2 = sum(fooof == "FOOOFed" & !is.finite(r_squared), na.rm = TRUE),
    pct_removed_low_r2 = 100 * rows_removed_low_r2 / pmax(fooof_rows_total, 1),
    .groups = "drop"
  ) %>%
  arrange(latency_ms)

dat <- dat %>% filter(!low_r2_mask)

message(sprintf(
  "FOOOF R2 filter (>= %.2f): removed %d/%d FOOOFed rows (%.2f%%). Rows missing R2 kept: %d.",
  R2_THRESHOLD, n_low_r2_removed, n_fooof_rows,
  100 * n_low_r2_removed / max(n_fooof_rows, 1), n_fooof_no_r2
))
message(sprintf(
  "After FOOOF R2 filtering: %d rows (from %d), %d universes (from %d).",
  nrow(dat), n_rows_before_r2, n_distinct(dat$universe_id), n_univ_before_r2
))
message("FOOOF R2 exclusions by latency_ms:")
for (i in seq_len(nrow(latency_exclusion_summary))) {
  row_i <- latency_exclusion_summary[i, ]
  message(sprintf(
    "  %s: removed %d/%d FOOOFed rows (%.2f%%), missing R2 kept: %d",
    row_i$latency_ms, row_i$rows_removed_low_r2, row_i$fooof_rows_total,
    row_i$pct_removed_low_r2, row_i$fooof_rows_missing_r2
  ))
}

write.csv(
  data.frame(
    task = "sternberg",
    fooof_method = FOOOF_METHOD,
    r2_threshold = R2_THRESHOLD,
    rows_before = n_rows_before_r2,
    rows_after = nrow(dat),
    rows_removed_low_r2 = n_low_r2_removed,
    fooof_rows_total = n_fooof_rows,
    fooof_rows_missing_r2 = n_fooof_no_r2,
    universes_before = n_univ_before_r2,
    universes_after = n_distinct(dat$universe_id)
  ),
  file.path(csv_dir, "MAG_multiverse_sternberg_fooof_trial_r2_filter_summary.csv"),
  row.names = FALSE
)

write.csv(
  latency_exclusion_summary %>%
    mutate(task = "sternberg", fooof_method = FOOOF_METHOD, r2_threshold = R2_THRESHOLD, .before = 1),
  file.path(csv_dir, "MAG_multiverse_sternberg_fooof_trial_r2_filter_summary_by_latency.csv"),
  row.names = FALSE
)

if ("gaze_measure" %in% names(dat)) {
  dat <- dat %>% filter(gaze_measure != "gaze_density" | is.na(gaze_measure))
  dat <- dat %>% filter(gaze_measure != "microsaccades")
}
dat <- dat %>% filter(!electrodes %in% c("all", "parietal"))
dat <- dat %>% filter(baseline_eeg != "dB", baseline_gaze != "dB")
message(sprintf("After filtering: %d rows, %d universes.", nrow(dat), n_distinct(dat$universe_id)))

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

save_with_method <- function(df, filename) {
  out_path <- file.path(csv_dir, filename)
  df_out <- df %>% mutate(fooof_method = FOOOF_METHOD, .before = 1)

  if (file.exists(out_path)) {
    old <- read.csv(out_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
    if (!("fooof_method" %in% names(old))) {
      old <- old %>% mutate(fooof_method = Sys.getenv("AOC_LEGACY_FOOOF_METHOD", unset = "singleFFT"), .before = 1)
    }
    old <- old %>% filter(fooof_method != FOOOF_METHOD)
    df_out <- bind_rows(old, df_out)
  }

  write.csv(df_out, out_path, row.names = FALSE)
}

run_lmer <- function(formula_obj, data_obj) {
  if (VERBOSE_LOGGING) {
    return(tryCatch(
      lmer(formula_obj, data = data_obj, control = lmerControl(optimizer = "bobyqa")),
      error = function(e) NULL
    ))
  }
  tryCatch(
    suppressWarnings(
      suppressMessages(
        lmer(formula_obj, data = data_obj, control = lmerControl(optimizer = "bobyqa"))
      )
    ),
    error = function(e) NULL
  )
}

fit_with_condition_slope <- function(formula_slope, formula_intercept, data) {
  fit_slope <- run_lmer(formula_slope, data)
  if (!is.null(fit_slope) && !isSingular(fit_slope, tol = 1e-4)) {
    return(list(fit = fit_slope, re_spec = "subject_intercept_condition_slope"))
  }

  fit_intercept <- run_lmer(formula_intercept, data)
  if (is.null(fit_intercept)) return(NULL)
  list(fit = fit_intercept, re_spec = "subject_intercept_only")
}

branch_cols <- c("electrodes", "fooof", "latency_ms", "alpha_type",
                 "gaze_measure", "baseline_eeg", "baseline_gaze")

# Condition setup
cond_levels  <- sort(unique(dat$Condition))
cond_labels  <- setNames(paste0("Set size ", cond_levels), as.character(cond_levels))
highest_cond <- max(as.numeric(as.character(cond_levels)))

# ========== MAIN MULTIVERSE (7 dimensions) ==========
message("Setting up main multiverse (7 dimensions)...")

M <- multiverse()

inside(M, {
  .elec   <- branch(electrodes,    "posterior", "occipital")
  .fooof  <- branch(fooof,         "FOOOFed", "nonFOOOFed")
  .lat    <- branch(latency_ms,    "0_1000ms", "0_2000ms", "1000_2000ms")
  .alpha  <- branch(alpha_type,    "canonical", "IAF")
  .gaze   <- branch(gaze_measure,  "scan_path_length", "gaze_velocity", "BCEA", "gaze_deviation")
  .bleeg  <- branch(baseline_eeg,  "raw", "pct_change")
  .blgaze <- branch(baseline_gaze, "raw", "pct_change")

  df <- dat %>%
    filter(electrodes == .elec, fooof == .fooof, latency_ms == .lat,
           alpha_type == .alpha, gaze_measure == .gaze,
           baseline_eeg == .bleeg, baseline_gaze == .blgaze) %>%
    filter(complete.cases(alpha, gaze_value, Condition, subjectID))

  df$gaze_value <- robust_z(df$gaze_value)
  df$alpha      <- robust_z(df$alpha)
  valid <- nrow(df) >= 10 && !any(is.nan(df$gaze_value)) && !any(is.nan(df$alpha))

  tid_int <- if (valid) {
    fit_obj <- fit_with_condition_slope(
      alpha ~ gaze_value * Condition + (1 + Condition || subjectID),
      alpha ~ gaze_value * Condition + (1 | subjectID),
      df
    )
    if (!is.null(fit_obj)) {
      broom.mixed::tidy(fit_obj$fit, conf.int = TRUE) %>%
        mutate(random_effects = fit_obj$re_spec)
    } else tibble()
  } else tibble()

  tid_cond <- if (valid) {
    bind_rows(lapply(as.character(cond_levels), function(cl) {
      dc <- df[df$Condition == cl, ]
      if (nrow(dc) < 5) return(tibble())
      fit_c <- run_lmer(alpha ~ gaze_value + (1 | subjectID), dc)
      if (is.null(fit_c)) return(tibble())
      tid_c <- broom.mixed::tidy(fit_c, conf.int = TRUE) %>% filter(term == "gaze_value")
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
  filter(grepl("gaze_value", term))

M_cond <- M_expanded %>%
  mutate(tid = map(.results, safe_extract, "tid_cond")) %>%
  filter(!map_lgl(tid, is.null)) %>%
  select(.universe, all_of(branch_cols), tid) %>%
  unnest(tid)

if (nrow(M_cond) == 0L) stop("No successful LMM fits.")

# Drop unstable universes (SE > 95th percentile)
se_thresh <- quantile(M_cond$std.error, 0.95, na.rm = TRUE)
bad_ids   <- M_cond %>% filter(std.error > se_thresh) %>% pull(.universe) %>% unique()
M_cond    <- M_cond %>% filter(!.universe %in% bad_ids)
M_int     <- M_int  %>% filter(!.universe %in% bad_ids)
message(sprintf("Dropped %d unstable universes (SE > %.4f). %d remain.",
                length(bad_ids), se_thresh, n_distinct(M_cond$.universe)))

M_cond <- add_sig(M_cond)
M_int  <- add_sig(M_int)

save_with_method(M_int, "MAG_multiverse_sternberg_results.csv")
save_with_method(M_cond, "MAG_multiverse_sternberg_conditions_results.csv")
message("Saved: MAG_multiverse_sternberg_results.csv, MAG_multiverse_sternberg_conditions_results.csv")

# Extract interaction term
highest_int_term <- paste0("gaze_value:Condition", highest_cond)
M_interaction    <- M_int %>% filter(term == highest_int_term)
if (nrow(M_interaction) > 0) {
  save_with_method(M_interaction, "MAG_multiverse_sternberg_interaction_results.csv")
  message(sprintf("Saved: MAG_multiverse_sternberg_interaction_results.csv (%d universes)", nrow(M_interaction)))
} else {
  message("WARNING: No interaction terms found for ", highest_int_term)
}

# ========== EEG-ONLY MULTIVERSE (5 dimensions) ==========
message("Setting up EEG-only multiverse (5 dimensions)...")

dat_eeg <- dat %>%
  filter(!electrodes %in% c("all", "parietal"), baseline_eeg != "dB") %>%
  select(subjectID, Condition, alpha, electrodes, fooof, latency_ms, alpha_type, baseline_eeg) %>%
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

  de$alpha <- robust_z(de$alpha)
  valid <- nrow(de) >= 10 && !any(is.nan(de$alpha))

  tid_ca <- if (valid) {
    fit_obj <- fit_with_condition_slope(
      alpha ~ Condition + (1 + Condition || subjectID),
      alpha ~ Condition + (1 | subjectID),
      de
    )
    if (!is.null(fit_obj)) {
      broom.mixed::tidy(fit_obj$fit, conf.int = TRUE) %>%
        filter(term == highest_alpha_term) %>%
        mutate(random_effects = fit_obj$re_spec)
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
  unnest(tid)

if (nrow(M_ca) > 0L) {
  M_ca <- add_sig(M_ca)
  save_with_method(M_ca, "MAG_multiverse_sternberg_condition_results.csv")
  message(sprintf("Saved: MAG_multiverse_sternberg_condition_results.csv (%d universes)", nrow(M_ca)))
} else {
  message("WARNING: No successful condition \u2192 alpha fits.")
}

# ========== GAZE-ONLY MULTIVERSE (3 dimensions) ==========
message("Setting up gaze-only multiverse (3 dimensions)...")

dat_gaze <- dat %>%
  select(subjectID, Condition, gaze_value, latency_ms, gaze_measure, baseline_gaze) %>%
  distinct()

highest_gaze_term <- paste0("Condition", highest_cond)

M_gaze <- multiverse()

inside(M_gaze, {
  .lat    <- branch(latency_ms,    "0_1000ms", "0_2000ms", "1000_2000ms")
  .gaze   <- branch(gaze_measure,  "scan_path_length", "gaze_velocity", "BCEA", "gaze_deviation")
  .blgaze <- branch(baseline_gaze, "raw", "pct_change")

  dg <- dat_gaze %>%
    filter(latency_ms == .lat, gaze_measure == .gaze, baseline_gaze == .blgaze) %>%
    filter(complete.cases(gaze_value, Condition, subjectID))

  dg$gaze_value <- robust_z(dg$gaze_value)
  valid <- nrow(dg) >= 10 && !any(is.nan(dg$gaze_value))

  tid_cg <- if (valid) {
    fit_obj <- fit_with_condition_slope(
      gaze_value ~ Condition + (1 + Condition || subjectID),
      gaze_value ~ Condition + (1 | subjectID),
      dg
    )
    if (!is.null(fit_obj)) {
      broom.mixed::tidy(fit_obj$fit, conf.int = TRUE) %>%
        filter(term == highest_gaze_term) %>%
        mutate(random_effects = fit_obj$re_spec)
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
  unnest(tid)

if (nrow(M_cg) > 0L) {
  M_cg <- add_sig(M_cg)
  save_with_method(M_cg, "MAG_multiverse_sternberg_condition_gaze_results.csv")
  message(sprintf("Saved: MAG_multiverse_sternberg_condition_gaze_results.csv (%d universes)", nrow(M_cg)))
} else {
  message("WARNING: No successful condition \u2192 gaze fits.")
}

# ========== APERIODIC MULTIVERSE ==========
# Aperiodic parameters (offset, exponent) from spectral parameterization.
# Only exists for FOOOFed universes. Reduced dimension space.
if ("aperiodic_offset" %in% names(dat) && "aperiodic_exponent" %in% names(dat)) {
  message("Setting up aperiodic multiverse...")

  dat_ap <- dat %>%
    filter(fooof == "FOOOFed") %>%
    filter(complete.cases(aperiodic_offset, aperiodic_exponent))

  # --- Aperiodic ~ gaze (4D: latency × electrodes × gaze_measure × gaze_baseline) ---
  dat_ap_gaze <- dat_ap %>%
    select(subjectID, Condition, aperiodic_offset, aperiodic_exponent,
           gaze_value, electrodes, latency_ms, gaze_measure, baseline_gaze) %>%
    distinct()

  M_ap_gaze <- multiverse()

  inside(M_ap_gaze, {
    .elec   <- branch(electrodes,    "posterior", "occipital")
    .lat    <- branch(latency_ms,    "0_1000ms", "0_2000ms", "1000_2000ms")
    .gaze   <- branch(gaze_measure,  "scan_path_length", "gaze_velocity", "BCEA", "gaze_deviation")
    .blgaze <- branch(baseline_gaze, "raw", "pct_change")

    dap <- dat_ap_gaze %>%
      filter(electrodes == .elec, latency_ms == .lat,
             gaze_measure == .gaze, baseline_gaze == .blgaze) %>%
      filter(complete.cases(aperiodic_exponent, aperiodic_offset, gaze_value, Condition, subjectID))

    dap$gaze_value <- robust_z(dap$gaze_value)
    dap$aperiodic_exponent <- robust_z(dap$aperiodic_exponent)
    dap$aperiodic_offset <- robust_z(dap$aperiodic_offset)
    valid <- nrow(dap) >= 10 && !any(is.nan(dap$gaze_value)) &&
             !any(is.nan(dap$aperiodic_exponent)) && !any(is.nan(dap$aperiodic_offset))

    tid_exp_gaze <- if (valid) {
      fit_obj <- fit_with_condition_slope(
        aperiodic_exponent ~ gaze_value + Condition + (1 + Condition || subjectID),
        aperiodic_exponent ~ gaze_value + Condition + (1 | subjectID),
        dap
      )
      if (!is.null(fit_obj)) {
        broom.mixed::tidy(fit_obj$fit, conf.int = TRUE) %>%
          filter(term == "gaze_value") %>%
          mutate(
            aperiodic_measure = "Exponent",
            random_effects = fit_obj$re_spec
          )
      } else tibble()
    } else tibble()

    tid_off_gaze <- if (valid) {
      fit_obj <- fit_with_condition_slope(
        aperiodic_offset ~ gaze_value + Condition + (1 + Condition || subjectID),
        aperiodic_offset ~ gaze_value + Condition + (1 | subjectID),
        dap
      )
      if (!is.null(fit_obj)) {
        broom.mixed::tidy(fit_obj$fit, conf.int = TRUE) %>%
          filter(term == "gaze_value") %>%
          mutate(
            aperiodic_measure = "Offset",
            random_effects = fit_obj$re_spec
          )
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
  )

  if (nrow(M_ap_gaze_results) > 0) {
    M_ap_gaze_results <- add_sig(M_ap_gaze_results)
    save_with_method(M_ap_gaze_results, "MAG_multiverse_sternberg_aperiodic_gaze_results.csv")
    message(sprintf("Saved: MAG_multiverse_sternberg_aperiodic_gaze_results.csv (%d rows)", nrow(M_ap_gaze_results)))
  }

  # --- Aperiodic ~ condition (2D: latency × electrodes) ---
  dat_ap_eeg <- dat_ap %>%
    select(subjectID, Condition, aperiodic_offset, aperiodic_exponent, electrodes, latency_ms) %>%
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
    valid <- nrow(dap) >= 10 && !any(is.nan(dap$aperiodic_exponent)) && !any(is.nan(dap$aperiodic_offset))

    tid_exp_cond <- if (valid) {
      fit_obj <- fit_with_condition_slope(
        aperiodic_exponent ~ Condition + (1 + Condition || subjectID),
        aperiodic_exponent ~ Condition + (1 | subjectID),
        dap
      )
      if (!is.null(fit_obj)) {
        broom.mixed::tidy(fit_obj$fit, conf.int = TRUE) %>%
          filter(term == highest_alpha_term) %>%
          mutate(
            aperiodic_measure = "Exponent",
            random_effects = fit_obj$re_spec
          )
      } else tibble()
    } else tibble()

    tid_off_cond <- if (valid) {
      fit_obj <- fit_with_condition_slope(
        aperiodic_offset ~ Condition + (1 + Condition || subjectID),
        aperiodic_offset ~ Condition + (1 | subjectID),
        dap
      )
      if (!is.null(fit_obj)) {
        broom.mixed::tidy(fit_obj$fit, conf.int = TRUE) %>%
          filter(term == highest_alpha_term) %>%
          mutate(
            aperiodic_measure = "Offset",
            random_effects = fit_obj$re_spec
          )
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
  )

  if (nrow(M_ap_cond_results) > 0) {
    M_ap_cond_results <- add_sig(M_ap_cond_results)
    save_with_method(M_ap_cond_results, "MAG_multiverse_sternberg_aperiodic_condition_results.csv")
    message(sprintf("Saved: MAG_multiverse_sternberg_aperiodic_condition_results.csv (%d rows)", nrow(M_ap_cond_results)))
  }

} else {
  message("Skipping aperiodic multiverse: columns not found in CSV.")
}

message("=== Sternberg multiverse ANALYSIS complete ===")
message("Result CSVs saved to: ", csv_dir)
