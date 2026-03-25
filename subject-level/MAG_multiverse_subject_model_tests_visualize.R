#!/usr/bin/env Rscript

# AOC Multiverse — Subject-Level Model Tests (Specification Curves)
#
# Goal:
#   Compare model-choice effects on subject-level multiverse figures.
#   - Gaze models: compare old gaze_value vs centered alternatives.
#   - Non-gaze models: test FOOOF subject-threshold scenarios only.
#
# Output directory:
#   /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/multiverse/subject-level/model_tests

library(tidyverse)
library(ggplot2)
library(patchwork)
library(lme4)
library(lmerTest)
library(broom.mixed)

csv_dir <- Sys.getenv(
  "AOC_MULTIVERSE_DIR",
  unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/MAG/data/multiverse"
)
r2_dir <- Sys.getenv(
  "AOC_R2_DIR",
  unset = "/Volumes/g_psyplafor_methlab$/Students/Arne/MAG/data/controls/multiverse"
)
FOOOF_METHOD <- Sys.getenv("AOC_FOOOF_METHOD", unset = "singleFFT")
MIN_SUBJECTS <- suppressWarnings(as.numeric(Sys.getenv("AOC_SUBJECT_MIN_N", unset = "20")))
if (!is.finite(MIN_SUBJECTS) || MIN_SUBJECTS < 5) MIN_SUBJECTS <- 20

storage_plot <- "/Volumes/g_psyplafor_methlab$/Students/Arne/MAG/figures/multiverse/subject-level/model_tests"
if (!dir.exists(storage_plot)) dir.create(storage_plot, recursive = TRUE)

sig_colors <- c("Positive" = "#33CC66", "Negative" = "#fe0000", "Non-significant" = "#d1d1d1")
sig_levels <- c("Positive", "Negative", "Non-significant")

elec_order <- c("posterior", "occipital")
lat_order  <- c("0_1000ms", "1000_2000ms", "0_2000ms")
gaze_order <- c("gaze_deviation", "gaze_velocity", "scan_path_length", "BCEA")

filter_labels <- c(
  "none" = "No FOOOF subject filter",
  "subject_min_r2_ge_0.90" = "Subject min FOOOF R2 >= .90",
  "subject_min_r2_ge_0.95" = "Subject min FOOOF R2 >= .95"
)

term_labels <- c(
  "gaze_value_z" = "gaze_value",
  "gaze_within_z" = "gaze_within",
  "gaze_between_z" = "gaze_between"
)

robust_z <- function(x, winsor = 3) {
  med <- median(x, na.rm = TRUE)
  ma <- mad(x, na.rm = TRUE)
  if (is.na(ma) || ma == 0) return(rep(NaN, length(x)))
  z <- (x - med) / ma
  z[z > winsor] <- winsor
  z[z < -winsor] <- -winsor
  z
}

add_sig <- function(df) {
  df %>%
    mutate(
      condition = factor(case_when(
        p.value < 0.05 & estimate > 0 ~ "Positive",
        p.value < 0.05 & estimate < 0 ~ "Negative",
        TRUE ~ "Non-significant"
      ), levels = sig_levels)
    )
}

run_lmer <- function(formula_obj, data_obj) {
  args <- list(
    formula = formula_obj,
    data = data_obj,
    control = lmerControl(optimizer = "bobyqa")
  )
  if ("model_weight" %in% names(data_obj) && any(is.finite(data_obj$model_weight))) {
    args$weights <- data_obj$model_weight
  }
  tryCatch(
    suppressWarnings(suppressMessages(do.call(lmer, args))),
    error = function(e) NULL
  )
}

valid_subjects <- function(df) {
  n_distinct(df$subjectID) >= MIN_SUBJECTS
}

plot_spec_curve <- function(df, panel_df, title_txt, subtitle_txt, out_path, h = 10) {
  if (nrow(df) == 0) return(invisible(NULL))
  ymax <- max(abs(c(df$conf.low, df$conf.high)), na.rm = TRUE) * 1.05
  if (!is.finite(ymax) || ymax == 0) ymax <- 0.1

  p_curve <- ggplot(df, aes(x = ordered_universe, y = estimate, color = condition)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha = 0.7) +
    geom_point(size = 1.1, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = sig_colors, name = "Significance") +
    labs(title = title_txt, subtitle = subtitle_txt, x = "Universe", y = expression(bold("Standardized " * beta))) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", plot.title = element_text(face = "bold")) +
    coord_cartesian(ylim = c(-ymax, ymax))

  p_panel <- ggplot(panel_df, aes(x = ordered_universe, y = value, fill = condition)) +
    geom_tile() +
    facet_grid(Variable ~ ., scales = "free_y", space = "free") +
    scale_fill_manual(values = sig_colors, name = "Significance") +
    labs(x = "Universe", y = "Analysis Decision") +
    theme_minimal(base_size = 11) +
    theme(strip.text.y = element_text(angle = 0, face = "bold"), legend.position = "bottom")

  g <- p_curve / p_panel + plot_layout(heights = c(0.8, 1.2))
  ggsave(out_path, g, width = 13, height = h, dpi = 300, bg = "white")
}

make_alpha_panel <- function(df) {
  df %>%
    transmute(ordered_universe, latency_ms, electrodes, fooof, baseline_eeg, alpha_type, gaze_measure, baseline_gaze, condition) %>%
    pivot_longer(cols = c(latency_ms, electrodes, fooof, baseline_eeg, alpha_type, gaze_measure, baseline_gaze),
                 names_to = "Variable", values_to = "value") %>%
    mutate(Variable = factor(Variable, levels = c("latency_ms", "electrodes", "fooof", "baseline_eeg", "alpha_type", "gaze_measure", "baseline_gaze")))
}

make_eeg_panel <- function(df) {
  df %>%
    transmute(ordered_universe, latency_ms, electrodes, fooof, baseline_eeg, alpha_type, condition) %>%
    pivot_longer(cols = c(latency_ms, electrodes, fooof, baseline_eeg, alpha_type),
                 names_to = "Variable", values_to = "value") %>%
    mutate(Variable = factor(Variable, levels = c("latency_ms", "electrodes", "fooof", "baseline_eeg", "alpha_type")))
}

make_gaze_panel <- function(df) {
  df %>%
    transmute(ordered_universe, latency_ms, gaze_measure, baseline_gaze, condition) %>%
    pivot_longer(cols = c(latency_ms, gaze_measure, baseline_gaze),
                 names_to = "Variable", values_to = "value") %>%
    mutate(Variable = factor(Variable, levels = c("latency_ms", "gaze_measure", "baseline_gaze")))
}

make_ap_panel <- function(df) {
  df %>%
    transmute(ordered_universe, latency_ms, electrodes, gaze_measure, baseline_gaze, condition) %>%
    pivot_longer(cols = c(latency_ms, electrodes, gaze_measure, baseline_gaze),
                 names_to = "Variable", values_to = "value") %>%
    mutate(Variable = factor(Variable, levels = c("latency_ms", "electrodes", "gaze_measure", "baseline_gaze")))
}

make_ap_cond_panel <- function(df) {
  df %>%
    transmute(ordered_universe, latency_ms, electrodes, fooof, baseline_eeg, condition) %>%
    pivot_longer(cols = c(latency_ms, electrodes, fooof, baseline_eeg), names_to = "Variable", values_to = "value") %>%
    mutate(Variable = factor(Variable, levels = c("latency_ms", "electrodes", "fooof", "baseline_eeg")))
}

get_subject_keep <- function(task, threshold) {
  if (!is.finite(threshold)) return(NULL)
  r2_path_method <- file.path(r2_dir, paste0("fooof_r2_", task, "_subject_", FOOOF_METHOD, ".csv"))
  r2_path <- if (file.exists(r2_path_method)) r2_path_method else file.path(r2_dir, paste0("fooof_r2_", task, "_subject.csv"))
  if (!file.exists(r2_path)) return(NULL)

  r2 <- read.csv(r2_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
  if (!all(c("subjectID", "r_squared") %in% names(r2))) return(NULL)
  r2 %>%
    mutate(subjectID = as.factor(subjectID)) %>%
    group_by(subjectID) %>%
    summarise(subject_min_r2 = min(r_squared, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(subject_min_r2), subject_min_r2 >= threshold) %>%
    pull(subjectID)
}

for (task in c("sternberg", "nback")) {
  # 1) Gaze model tests from ablation details (old vs centered gaze handling).
  alpha_detail_path <- file.path(csv_dir, paste0("multiverse_", task, "_subject_ablation_alpha_detail.csv"))
  ap_detail_path <- file.path(csv_dir, paste0("multiverse_", task, "_subject_ablation_aperiodic_detail.csv"))
  if (file.exists(alpha_detail_path) && file.exists(ap_detail_path)) {
    A <- read.csv(alpha_detail_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
    AP <- read.csv(ap_detail_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))

    if (nrow(A) > 0) {
      A <- add_sig(A) %>%
        filter((ablation_step == "S0_old_unweighted_gaze_value" & term == "gaze_value_z") |
                 (ablation_step == "S2_weighted_switch_to_within_between" & term %in% c("gaze_within_z", "gaze_between_z")))

      for (filt in unique(A$r2_subject_filter)) {
        for (term in unique(A$term)) {
          d <- A %>% filter(r2_subject_filter == .env$filt, term == .env$term)
          if (nrow(d) == 0) next
          conds <- suppressWarnings(as.numeric(unique(d$cond_label)))
          if (all(is.na(conds))) next
          high_cond <- as.character(max(conds, na.rm = TRUE))
          d <- d %>%
            filter(cond_label == high_cond) %>%
            mutate(.lat_ord = match(latency_ms, lat_order),
                   .elec_ord = match(electrodes, elec_order),
                   .gaze_ord = match(gaze_measure, gaze_order)) %>%
            arrange(ablation_step, fooof, .lat_ord, .elec_ord, baseline_eeg, alpha_type, .gaze_ord, baseline_gaze) %>%
            select(-.lat_ord, -.elec_ord, -.gaze_ord) %>%
            mutate(ordered_universe = row_number())
          panel_d <- make_alpha_panel(d)
          term_lab <- ifelse(term %in% names(term_labels), term_labels[[term]], term)
          filt_lab <- ifelse(filt %in% names(filter_labels), filter_labels[[filt]], filt)
          out <- file.path(storage_plot, paste0("AOC_", task, "_subject_alpha_gaze_test_", term_lab, "_", filt, ".png"))
          plot_spec_curve(d, panel_d,
                          paste0("Alpha ~ gaze (", toupper(task), ") — ", term_lab),
                          paste0(filt_lab, "; highest condition"),
                          out)
        }
      }
    }

    if (nrow(AP) > 0) {
      AP <- add_sig(AP) %>%
        filter((ablation_step == "A0_old_unweighted_gaze_value" & term == "gaze_value_z") |
                 (ablation_step %in% c("A2_weighted_switch_to_within_between", "A3_weighted_within_between_plus_condition") &
                    term %in% c("gaze_within_z", "gaze_between_z")))

      for (meas in unique(AP$aperiodic_measure)) {
        for (filt in unique(AP$r2_subject_filter)) {
          for (term in unique(AP$term)) {
            d <- AP %>% filter(aperiodic_measure == .env$meas, r2_subject_filter == .env$filt, term == .env$term)
            if (nrow(d) == 0) next
            d <- d %>%
              mutate(.lat_ord = match(latency_ms, lat_order),
                     .elec_ord = match(electrodes, elec_order),
                     .gaze_ord = match(gaze_measure, gaze_order)) %>%
              arrange(ablation_step, .lat_ord, .elec_ord, .gaze_ord, baseline_gaze) %>%
              select(-.lat_ord, -.elec_ord, -.gaze_ord) %>%
              mutate(ordered_universe = row_number())
            panel_d <- make_ap_panel(d)
            term_lab <- ifelse(term %in% names(term_labels), term_labels[[term]], term)
            filt_lab <- ifelse(filt %in% names(filter_labels), filter_labels[[filt]], filt)
            out <- file.path(storage_plot, paste0("AOC_", task, "_subject_aperiodic_", tolower(meas), "_gaze_test_", term_lab, "_", filt, ".png"))
            plot_spec_curve(d, panel_d,
                            paste0("Aperiodic ", tolower(meas), " ~ gaze (", toupper(task), ") — ", term_lab),
                            paste0(filt_lab, "; A0/A2/A3 variants"),
                            out)
          }
        }
      }
    }
  } else {
    message("Missing ablation detail CSVs for ", task, ".")
  }

  # 2) Non-gaze model tests: threshold-only comparisons.
  subject_path_method <- file.path(csv_dir, paste0("multiverse_", task, "_subject_", FOOOF_METHOD, ".csv"))
  subject_path <- if (file.exists(subject_path_method)) subject_path_method else file.path(csv_dir, paste0("multiverse_", task, "_subject.csv"))
  if (!file.exists(subject_path)) next

  dat <- read.csv(subject_path, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
  dat$subjectID <- as.factor(dat$subjectID)
  dat$Condition <- as.factor(dat$Condition)
  if (!("n_trials_subject_condition" %in% names(dat))) dat$n_trials_subject_condition <- 1
  dat$n_trials_subject_condition <- suppressWarnings(as.numeric(dat$n_trials_subject_condition))
  dat$n_trials_subject_condition[!is.finite(dat$n_trials_subject_condition) | dat$n_trials_subject_condition <= 0] <- 1
  if ("gaze_measure" %in% names(dat)) {
    dat <- dat %>% filter(gaze_measure != "gaze_density" | is.na(gaze_measure))
  }
  dat <- dat %>% filter(!electrodes %in% c("all", "parietal"), baseline_eeg != "dB", baseline_gaze != "dB")
  cond_levels <- sort(unique(dat$Condition))
  high_cond <- max(as.numeric(as.character(cond_levels)))
  high_term <- paste0("Condition", high_cond)

  scenarios <- tibble(
    filter_id = c("none", "subject_min_r2_ge_0.90", "subject_min_r2_ge_0.95"),
    threshold = c(NA_real_, 0.90, 0.95)
  )

  for (i in seq_len(nrow(scenarios))) {
    sid <- scenarios$filter_id[i]
    thr <- scenarios$threshold[i]
    keep_subj <- get_subject_keep(task, thr)
    ds <- dat
    if (is.finite(thr) && !is.null(keep_subj)) ds <- ds %>% filter(subjectID %in% keep_subj)
    filt_lab <- ifelse(sid %in% names(filter_labels), filter_labels[[sid]], sid)

    # condition -> alpha
    D_eeg <- ds %>% select(subjectID, Condition, alpha, n_trials_subject_condition, electrodes, fooof, latency_ms, alpha_type, baseline_eeg) %>% distinct()
    specs_eeg <- D_eeg %>% distinct(electrodes, fooof, latency_ms, alpha_type, baseline_eeg)
    out_eeg <- vector("list", nrow(specs_eeg))
    k <- 1
    for (r in seq_len(nrow(specs_eeg))) {
      key <- specs_eeg[r, ]
      d <- D_eeg %>% semi_join(key, by = names(key)) %>% filter(complete.cases(alpha, Condition, subjectID))
      if (!valid_subjects(d)) next
      d <- d %>% mutate(alpha = robust_z(alpha), model_weight = n_trials_subject_condition)
      fit <- run_lmer(alpha ~ Condition + (1 | subjectID), d)
      if (is.null(fit)) next
      td <- broom.mixed::tidy(fit, conf.int = TRUE) %>% filter(term == high_term)
      if (nrow(td) == 0) next
      out_eeg[[k]] <- bind_cols(td, key)
      k <- k + 1
    }
    M_ca <- bind_rows(out_eeg)
    if (nrow(M_ca) > 0) {
      M_ca <- add_sig(M_ca) %>%
        mutate(.lat_ord = match(latency_ms, lat_order), .elec_ord = match(electrodes, elec_order)) %>%
        arrange(.lat_ord, .elec_ord, fooof, baseline_eeg, alpha_type) %>%
        select(-.lat_ord, -.elec_ord) %>%
        mutate(ordered_universe = row_number())
      panel <- make_eeg_panel(M_ca)
      out <- file.path(storage_plot, paste0("AOC_", task, "_subject_condition_alpha_threshold_", sid, ".png"))
      plot_spec_curve(M_ca, panel, paste0("alpha ~ Condition (", toupper(task), ")"), filt_lab, out, h = 9)
    }

    # condition -> gaze
    D_gaze <- ds %>% select(subjectID, Condition, gaze_value, n_trials_subject_condition, latency_ms, gaze_measure, baseline_gaze) %>% distinct()
    specs_gaze <- D_gaze %>% distinct(latency_ms, gaze_measure, baseline_gaze)
    out_gaze <- vector("list", nrow(specs_gaze))
    k <- 1
    for (r in seq_len(nrow(specs_gaze))) {
      key <- specs_gaze[r, ]
      d <- D_gaze %>% semi_join(key, by = names(key)) %>% filter(complete.cases(gaze_value, Condition, subjectID))
      if (!valid_subjects(d)) next
      d <- d %>% mutate(gaze_value = robust_z(gaze_value), model_weight = n_trials_subject_condition)
      fit <- run_lmer(gaze_value ~ Condition + (1 | subjectID), d)
      if (is.null(fit)) next
      td <- broom.mixed::tidy(fit, conf.int = TRUE) %>% filter(term == high_term)
      if (nrow(td) == 0) next
      out_gaze[[k]] <- bind_cols(td, key)
      k <- k + 1
    }
    M_cg <- bind_rows(out_gaze)
    if (nrow(M_cg) > 0) {
      M_cg <- add_sig(M_cg) %>%
        mutate(.lat_ord = match(latency_ms, lat_order), .gaze_ord = match(gaze_measure, gaze_order)) %>%
        arrange(.lat_ord, .gaze_ord, baseline_gaze) %>%
        select(-.lat_ord, -.gaze_ord) %>%
        mutate(ordered_universe = row_number())
      panel <- make_gaze_panel(M_cg)
      out <- file.path(storage_plot, paste0("AOC_", task, "_subject_condition_gaze_threshold_", sid, ".png"))
      plot_spec_curve(M_cg, panel, paste0("gaze ~ Condition (", toupper(task), ")"), filt_lab, out, h = 8)
    }

    # aperiodic ~ condition
    D_ap <- ds %>%
      select(subjectID, Condition, aperiodic_offset, aperiodic_exponent, n_trials_subject_condition, electrodes, fooof, latency_ms, baseline_eeg) %>%
      distinct()
    specs_ap <- D_ap %>% distinct(electrodes, fooof, latency_ms, baseline_eeg)
    out_ap <- list()
    k <- 1
    for (r in seq_len(nrow(specs_ap))) {
      key <- specs_ap[r, ]
      d <- D_ap %>% semi_join(key, by = names(key))
      if (!valid_subjects(d)) next
      d <- d %>% mutate(
        aperiodic_exponent = robust_z(aperiodic_exponent),
        aperiodic_offset = robust_z(aperiodic_offset),
        model_weight = n_trials_subject_condition
      )
      fit_e <- run_lmer(aperiodic_exponent ~ Condition + (1 | subjectID), d)
      fit_o <- run_lmer(aperiodic_offset ~ Condition + (1 | subjectID), d)
      if (!is.null(fit_e)) {
        td <- broom.mixed::tidy(fit_e, conf.int = TRUE) %>% filter(term == high_term)
        if (nrow(td) > 0) { out_ap[[k]] <- bind_cols(td, key) %>% mutate(aperiodic_measure = "Exponent"); k <- k + 1 }
      }
      if (!is.null(fit_o)) {
        td <- broom.mixed::tidy(fit_o, conf.int = TRUE) %>% filter(term == high_term)
        if (nrow(td) > 0) { out_ap[[k]] <- bind_cols(td, key) %>% mutate(aperiodic_measure = "Offset"); k <- k + 1 }
      }
    }
    M_apc <- bind_rows(out_ap)
    if (nrow(M_apc) > 0) {
      for (meas in c("Exponent", "Offset")) {
        d <- M_apc %>% filter(aperiodic_measure == meas)
        if (nrow(d) == 0) next
        d <- add_sig(d) %>%
          mutate(.lat_ord = match(latency_ms, lat_order), .elec_ord = match(electrodes, elec_order)) %>%
          arrange(.lat_ord, .elec_ord, fooof, baseline_eeg) %>%
          select(-.lat_ord, -.elec_ord) %>%
          mutate(ordered_universe = row_number())
        panel <- make_ap_cond_panel(d)
        out <- file.path(storage_plot, paste0("AOC_", task, "_subject_aperiodic_", tolower(meas), "_condition_threshold_", sid, ".png"))
        plot_spec_curve(d, panel, paste0("aperiodic_", tolower(meas), " ~ Condition (", toupper(task), ")"), filt_lab, out, h = 8)
      }
    }
  }

  message("Saved model-test specification curves for ", task, " to: ", storage_plot)
}

message("=== Subject-level model tests specification-curve visualization complete ===")
