#!/usr/bin/env Rscript
# ==============================================================================
# 22_evaluate_forecasts.R
#
# Evaluate rolling-origin forecast results: CRPS, log-score, calibration,
# skill scores, and generate summary figures.
#
# Input:
#   results/forecasts/forecast_entropy{0,1}_origin{N}.rds
#
# Output:
#   results/figures/kernel_forecast_crps.png
#   results/figures/kernel_forecast_calibration.png
#   results/figures/kernel_forecast_trajectories.png
#   results/figures/kernel_forecast_entropy_effect.png
#   results/forecast_evaluation_summary.csv
# ==============================================================================

library(tidyverse)
library(patchwork)

# Set up paths
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- grep("--file=", args, value = TRUE)
  if (length(script_path) > 0) {
    setwd(dirname(normalizePath(sub("--file=", "", script_path))))
  }
}

dir.create("../results/figures", showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank()))

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("FORECAST EVALUATION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD ALL FORECAST RESULTS
# ==============================================================================

cat("Loading forecast results...\n")

forecast_files <- list.files("../results/forecasts",
                             pattern = "^forecast_entropy[01]_origin[0-9]+\\.rds$",
                             full.names = TRUE)

if (length(forecast_files) == 0) {
  stop("No forecast files found. Run 21_run_forecasts.R first.")
}

results <- lapply(forecast_files, readRDS)
names(results) <- basename(forecast_files)

cat(sprintf("  Loaded %d forecast files\n", length(results)))

# Extract configuration from first result
HORIZONS <- c(4, 8, 13)
S <- results[[1]]$S

# ==============================================================================
# 2. CRPS FUNCTION
# ==============================================================================

# CRPS for a single observation against a sample of predictions
# CRPS(F, y) = E|X - y| - 0.5 * E|X - X'|
crps_sample_fn <- function(y, samples) {
  # Use scoringRules if available, otherwise manual
  if (requireNamespace("scoringRules", quietly = TRUE)) {
    return(scoringRules::crps_sample(y = y, dat = samples))
  }

  # Manual implementation
  n <- length(samples)
  term1 <- mean(abs(samples - y))

  # E|X - X'| via sorted samples (efficient O(n log n) formula)
  sorted <- sort(samples)
  weights <- (2 * seq_len(n) - n - 1) / n
  term2 <- 2 * sum(weights * sorted) / n

  term1 - 0.5 * abs(term2)
}

# ==============================================================================
# 3. COMPUTE EVALUATION METRICS
# ==============================================================================

cat("\nComputing evaluation metrics...\n")

eval_rows <- list()

for (res_name in names(results)) {
  res <- results[[res_name]]
  n_train <- res$N_train
  n_forecast <- res$N_forecast
  n_model <- res$N_model
  use_ent <- res$use_entropy

  # cases_forward draws: matrix (n_draws x N_model)
  # Forecast period columns: (n_train+1):n_model
  cases_fwd <- res$cases_forward
  rt_draws <- res$Rt
  obs_cases <- res$observed_cases  # full cases vector (length N)

  for (h in HORIZONS) {
    # Model index for this horizon
    idx <- n_train + h
    if (idx > n_model) next

    # Absolute index into cases vector
    abs_idx <- S + idx
    if (abs_idx > length(obs_cases)) next

    y_obs <- obs_cases[abs_idx]

    # Column name for draws (1-indexed)
    col_name <- sprintf("cases_forward[%d]", idx)
    if (!(col_name %in% colnames(cases_fwd))) next

    samples <- cases_fwd[, col_name]

    # CRPS
    crps_val <- crps_sample_fn(y_obs, samples)

    # Log-score (predictive log density, approximated by KDE)
    log_score <- tryCatch({
      d <- density(samples, from = max(0, min(samples) - 50),
                   to = max(samples) + 50, n = 1024)
      log_pred <- approx(d$x, d$y, xout = y_obs)$y
      if (is.na(log_pred) || log_pred <= 0) -Inf else log(log_pred)
    }, error = function(e) NA_real_)

    # Interval coverage
    q025 <- quantile(samples, 0.025)
    q975 <- quantile(samples, 0.975)
    q10  <- quantile(samples, 0.10)
    q90  <- quantile(samples, 0.90)
    q25  <- quantile(samples, 0.25)
    q75  <- quantile(samples, 0.75)

    in_95 <- as.integer(y_obs >= q025 & y_obs <= q975)
    in_80 <- as.integer(y_obs >= q10 & y_obs <= q90)
    in_50 <- as.integer(y_obs >= q25 & y_obs <= q75)

    # Interval width (95%)
    width_95 <- unname(q975 - q025)

    # Bias (median prediction - observed)
    bias <- median(samples) - y_obs

    # --- Baseline 1: Rt = 1 -> cases = infectious_pressure ---
    ip_baseline <- 0
    gi <- res$gi
    for (s_lag in 1:length(gi)) {
      case_idx <- abs_idx - s_lag
      if (case_idx >= 1 && case_idx <= length(obs_cases)) {
        ip_baseline <- ip_baseline + obs_cases[case_idx] * gi[s_lag]
      }
    }
    crps_baseline_rt1 <- abs(y_obs - ip_baseline)

    # --- Baseline 2: Random walk (last observed Rt held constant) ---
    # Use last training-period cases to estimate Rt at the boundary
    # Last observed case ratio: cases[N_train] / infectious_pressure[N_train]
    last_train_idx <- n_train  # last training week (relative to model start)
    last_ip <- 0
    for (s_lag in 1:length(gi)) {
      case_idx_last <- last_train_idx + length(gi) - s_lag  # absolute index
      if (case_idx_last >= 1 && case_idx_last <= length(obs_cases)) {
        last_ip <- last_ip + obs_cases[case_idx_last] * gi[s_lag]
      }
    }
    # Last Rt ≈ cases[last] / infectious_pressure[last]
    last_cases <- obs_cases[last_train_idx + length(gi)]
    rt_last <- if (last_ip > 0) last_cases / last_ip else 1.0
    # Random walk forecast: Rt stays at rt_last
    rw_predicted <- rt_last * ip_baseline
    crps_baseline_rw <- abs(y_obs - rw_predicted)

    eval_rows[[length(eval_rows) + 1]] <- tibble(
      file = res_name,
      use_entropy = use_ent,
      N_train = n_train,
      horizon = h,
      date_forecast = res$dates_model[idx],
      y_obs = y_obs,
      median_pred = median(samples),
      crps = crps_val,
      log_score = log_score,
      in_50 = in_50,
      in_80 = in_80,
      in_95 = in_95,
      width_95 = width_95,
      bias = bias,
      crps_baseline_rt1 = crps_baseline_rt1,
      crps_baseline_rw = crps_baseline_rw,
      skill_vs_rt1 = 1 - crps_val / max(crps_baseline_rt1, 1e-6),
      skill_vs_rw = 1 - crps_val / max(crps_baseline_rw, 1e-6)
    )
  }
}

eval_df <- bind_rows(eval_rows)

cat(sprintf("  Total evaluation points: %d\n", nrow(eval_df)))
cat(sprintf("  Entropy options: %s\n",
            paste(sort(unique(eval_df$use_entropy)), collapse = ", ")))
cat(sprintf("  Horizons: %s\n",
            paste(sort(unique(eval_df$horizon)), collapse = ", ")))

# ==============================================================================
# 4. SUMMARY TABLE
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("FORECAST PERFORMANCE SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

summary_df <- eval_df |>
  group_by(use_entropy, horizon) |>
  summarize(
    n_origins = n(),
    mean_crps = mean(crps, na.rm = TRUE),
    mean_log_score = mean(log_score[is.finite(log_score)], na.rm = TRUE),
    coverage_50 = mean(in_50, na.rm = TRUE),
    coverage_80 = mean(in_80, na.rm = TRUE),
    coverage_95 = mean(in_95, na.rm = TRUE),
    mean_width_95 = mean(width_95, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    mean_crps_baseline_rt1 = mean(crps_baseline_rt1, na.rm = TRUE),
    mean_crps_baseline_rw = mean(crps_baseline_rw, na.rm = TRUE),
    mean_skill_vs_rt1 = mean(skill_vs_rt1, na.rm = TRUE),
    mean_skill_vs_rw = mean(skill_vs_rw, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    entropy_label = ifelse(use_entropy == 1, "With entropy", "Without entropy")
  )

print(summary_df |> select(entropy_label, horizon, n_origins, mean_crps,
                            mean_skill_vs_rt1, mean_skill_vs_rw,
                            coverage_50, coverage_80, coverage_95,
                            mean_bias))

# Save detailed evaluation
write_csv(eval_df, "../results/forecast_evaluation_detail.csv")
write_csv(summary_df, "../results/forecast_evaluation_summary.csv")
cat("\n  Saved: results/forecast_evaluation_summary.csv\n")
cat("  Saved: results/forecast_evaluation_detail.csv\n")

# ==============================================================================
# 5. FIGURE 1: CRPS BY HORIZON
# ==============================================================================

cat("\nGenerating figures...\n")

p_crps <- eval_df |>
  mutate(entropy_label = ifelse(use_entropy == 1, "With entropy", "Without entropy")) |>
  ggplot(aes(x = factor(horizon), y = crps, fill = entropy_label)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(values = c("Without entropy" = "#4575b4",
                                "With entropy" = "#d73027")) +
  labs(
    title = "Forecast CRPS by Horizon",
    subtitle = "Lower is better. Matern 3/2 kernel.",
    x = "Forecast horizon (weeks)",
    y = "CRPS (cases)",
    fill = NULL
  ) +
  theme(legend.position = "bottom")

ggsave("../results/figures/kernel_forecast_crps.png", p_crps,
       width = 7, height = 5, dpi = 300)
cat("  Saved: kernel_forecast_crps.png\n")

# ==============================================================================
# 6. FIGURE 2: CALIBRATION
# ==============================================================================

# Compute empirical coverage at multiple nominal levels
nominal_levels <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

calib_rows <- list()
for (res_name in names(results)) {
  res <- results[[res_name]]
  n_train <- res$N_train
  n_model <- res$N_model
  use_ent <- res$use_entropy
  cases_fwd <- res$cases_forward
  obs_cases <- res$observed_cases

  for (h in HORIZONS) {
    idx <- n_train + h
    if (idx > n_model) next
    abs_idx <- S + idx
    if (abs_idx > length(obs_cases)) next

    y_obs <- obs_cases[abs_idx]
    col_name <- sprintf("cases_forward[%d]", idx)
    if (!(col_name %in% colnames(cases_fwd))) next
    samples <- cases_fwd[, col_name]

    for (nom in nominal_levels) {
      lo <- quantile(samples, (1 - nom) / 2)
      hi <- quantile(samples, 1 - (1 - nom) / 2)
      covered <- as.integer(y_obs >= lo & y_obs <= hi)
      calib_rows[[length(calib_rows) + 1]] <- tibble(
        use_entropy = use_ent,
        horizon = h,
        nominal = nom,
        covered = covered
      )
    }
  }
}

calib_df <- bind_rows(calib_rows) |>
  group_by(use_entropy, horizon, nominal) |>
  summarize(empirical = mean(covered), .groups = "drop") |>
  mutate(entropy_label = ifelse(use_entropy == 1, "With entropy", "Without entropy"))

p_calib <- calib_df |>
  ggplot(aes(x = nominal, y = empirical, color = entropy_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~paste0(horizon, "-week horizon")) +
  scale_color_manual(values = c("Without entropy" = "#4575b4",
                                 "With entropy" = "#d73027")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Forecast Calibration",
    subtitle = "Nominal vs empirical coverage. Diagonal = perfect calibration.",
    x = "Nominal coverage",
    y = "Empirical coverage",
    color = NULL
  ) +
  theme(legend.position = "bottom")

ggsave("../results/figures/kernel_forecast_calibration.png", p_calib,
       width = 10, height = 4.5, dpi = 300)
cat("  Saved: kernel_forecast_calibration.png\n")

# ==============================================================================
# 7. FIGURE 3: EXAMPLE FORECAST TRAJECTORIES
# ==============================================================================

# Pick 2-3 origins to show trajectories
example_origins <- sort(unique(eval_df$N_train))
if (length(example_origins) > 3) {
  # Pick first, middle, last
  idx_sel <- c(1, ceiling(length(example_origins) / 2), length(example_origins))
  example_origins <- example_origins[idx_sel]
}

traj_plots <- list()

for (eo in example_origins) {
  # Use entropy=0 version for trajectories
  res_name <- sprintf("forecast_entropy0_origin%d.rds", eo)
  if (!(res_name %in% names(results))) next
  res <- results[[res_name]]

  cases_fwd <- res$cases_forward
  obs_cases <- res$observed_cases
  dates <- res$dates_model
  n_model <- res$N_model

  # Build trajectory summary for forecast period
  forecast_idx <- (eo + 1):min(eo + max(HORIZONS), n_model)
  if (length(forecast_idx) == 0) next

  traj_data <- tibble(
    date = dates[forecast_idx],
    observed = obs_cases[S + forecast_idx]
  )

  # Extract quantiles from posterior
  medians <- numeric(length(forecast_idx))
  q025 <- numeric(length(forecast_idx))
  q975 <- numeric(length(forecast_idx))
  q10 <- numeric(length(forecast_idx))
  q90 <- numeric(length(forecast_idx))

  for (j in seq_along(forecast_idx)) {
    col_name <- sprintf("cases_forward[%d]", forecast_idx[j])
    if (col_name %in% colnames(cases_fwd)) {
      samp <- cases_fwd[, col_name]
      medians[j] <- median(samp)
      q025[j] <- quantile(samp, 0.025)
      q975[j] <- quantile(samp, 0.975)
      q10[j] <- quantile(samp, 0.10)
      q90[j] <- quantile(samp, 0.90)
    }
  }

  traj_data <- traj_data |>
    mutate(median_pred = medians, q025 = q025, q975 = q975,
           q10 = q10, q90 = q90)

  # Also show some training context
  context_start <- max(1, eo - 12)
  context_idx <- context_start:eo
  context_data <- tibble(
    date = dates[context_idx],
    observed = obs_cases[S + context_idx]
  )

  p_traj <- ggplot() +
    # Training context
    geom_line(data = context_data, aes(x = date, y = observed),
              color = "grey40", linewidth = 0.5) +
    # Forecast intervals
    geom_ribbon(data = traj_data, aes(x = date, ymin = q025, ymax = q975),
                fill = "#4575b4", alpha = 0.2) +
    geom_ribbon(data = traj_data, aes(x = date, ymin = q10, ymax = q90),
                fill = "#4575b4", alpha = 0.3) +
    geom_line(data = traj_data, aes(x = date, y = median_pred),
              color = "#4575b4", linewidth = 0.8) +
    # Observed in forecast window
    geom_point(data = traj_data, aes(x = date, y = observed),
               color = "black", size = 1.5) +
    # Vertical line at origin
    geom_vline(xintercept = dates[eo], linetype = "dashed", color = "red") +
    labs(
      title = sprintf("Origin: %s (N_train=%d)", dates[eo], eo),
      x = NULL, y = "Cases"
    )

  traj_plots[[length(traj_plots) + 1]] <- p_traj
}

if (length(traj_plots) > 0) {
  p_trajectories <- wrap_plots(traj_plots, ncol = 1) +
    plot_annotation(
      title = "Example Forecast Trajectories",
      subtitle = "Red dashed line = forecast origin. Blue = predicted (median + 80/95% CI). Black dots = observed."
    )

  ggsave("../results/figures/kernel_forecast_trajectories.png", p_trajectories,
         width = 10, height = 3.5 * length(traj_plots), dpi = 300)
  cat("  Saved: kernel_forecast_trajectories.png\n")
}

# ==============================================================================
# 8. FIGURE 4: ENTROPY EFFECT (PAIRED COMPARISON)
# ==============================================================================

# Paired comparison: same origin, with vs without entropy
paired_df <- eval_df |>
  select(N_train, horizon, use_entropy, crps) |>
  pivot_wider(names_from = use_entropy, values_from = crps,
              names_prefix = "crps_ent") |>
  filter(!is.na(crps_ent0), !is.na(crps_ent1)) |>
  mutate(
    crps_diff = crps_ent1 - crps_ent0,  # negative = entropy helps
    pct_improvement = -100 * crps_diff / crps_ent0
  )

if (nrow(paired_df) > 0) {
  p_entropy <- paired_df |>
    ggplot(aes(x = factor(horizon), y = pct_improvement)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_boxplot(fill = "#fee090", alpha = 0.7, outlier.size = 1) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
    labs(
      title = "Effect of Shannon Entropy Covariate on Forecast Skill",
      subtitle = "Positive = entropy improves CRPS. Each point = one forecast origin.",
      x = "Forecast horizon (weeks)",
      y = "CRPS improvement (%)"
    )

  ggsave("../results/figures/kernel_forecast_entropy_effect.png", p_entropy,
         width = 7, height = 5, dpi = 300)
  cat("  Saved: kernel_forecast_entropy_effect.png\n")

  # Print summary
  cat("\n")
  cat("=" |> rep(70) |> paste(collapse = ""), "\n")
  cat("ENTROPY EFFECT SUMMARY\n")
  cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

  paired_summary <- paired_df |>
    group_by(horizon) |>
    summarize(
      n = n(),
      mean_improvement_pct = mean(pct_improvement),
      median_improvement_pct = median(pct_improvement),
      prop_improved = mean(crps_diff < 0),
      .groups = "drop"
    )
  print(paired_summary)
}

# ==============================================================================
# 9. FINAL SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("EVALUATION COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Figures saved:\n")
cat("  results/figures/kernel_forecast_crps.png\n")
cat("  results/figures/kernel_forecast_calibration.png\n")
cat("  results/figures/kernel_forecast_trajectories.png\n")
cat("  results/figures/kernel_forecast_entropy_effect.png\n")
cat("\nTables saved:\n")
cat("  results/forecast_evaluation_summary.csv\n")
cat("  results/forecast_evaluation_detail.csv\n")
