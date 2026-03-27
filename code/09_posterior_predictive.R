#!/usr/bin/env Rscript
# ==============================================================================
# 09_posterior_predictive.R
#
# Posterior predictive checks for dengue Rt estimation model
# Verifies the fitted model can reproduce key features of observed data
#
# Input: data/model_data.rds, results/fit_model3.rds (or fit_model2.rds)
# Output: results/figures/posterior_predictive_*.png
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
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

# Create directories
dir.create("../results", showWarnings = FALSE)
dir.create("../results/figures", showWarnings = FALSE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("POSTERIOR PREDICTIVE CHECKS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA AND FIT
# ==============================================================================

cat("Loading data and fitted model...\n")

model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data
metadata <- model_data$metadata

# Try to load the best available fit
fit_file <- NULL
for (f in c("../results/fit_model3.rds",
            "../results/fit_model2_short_gp.rds",
            "../results/fit_model2.rds")) {
  if (file.exists(f)) {
    fit_file <- f
    break
  }
}

if (is.null(fit_file)) {
  stop("No fitted model found. Run 13_fit_model3.R first.")
}

cat(sprintf("  Loading: %s\n", fit_file))
fit <- readRDS(fit_file)

# Extract observed cases
S <- stan_data$S
N <- stan_data$N
N_model <- stan_data$N_model
y <- stan_data$cases[(S + 1):N]

cat(sprintf("  N_model = %d weeks\n", N_model))
cat(sprintf("  Observed cases: min=%.0f, median=%.0f, max=%.0f\n",
            min(y), median(y), max(y)))

# ==============================================================================
# 2. EXTRACT POSTERIOR PREDICTIVE SAMPLES
# ==============================================================================

cat("\nExtracting posterior predictive samples...\n")

# Extract y_rep (cases_pred in our model)
y_rep <- fit$draws("cases_pred", format = "matrix")

cat(sprintf("  Posterior samples: %d draws x %d time points\n",
            nrow(y_rep), ncol(y_rep)))

# Ensure dimensions match (handle potential off-by-one issues)
if (ncol(y_rep) != length(y)) {
  cat(sprintf("  Note: Adjusting dimensions (y_rep has %d cols, y has %d elements)\n",
              ncol(y_rep), length(y)))
  # Take the minimum to ensure alignment
  n_common <- min(ncol(y_rep), length(y))
  y_rep <- y_rep[, 1:n_common]
  y <- y[1:n_common]
  N_model <- n_common
}

# Also extract Rt for context
Rt_draws <- fit$draws("Rt", format = "matrix")

# ==============================================================================
# 3. BAYESPLOT POSTERIOR PREDICTIVE CHECKS
# ==============================================================================

cat("\nGenerating posterior predictive visualizations...\n")

# Set theme
theme_set(theme_minimal(base_size = 12))
color_scheme_set("brightblue")

# --- PPC 1: Density overlay ---
# Compare distribution of replicated data to observed
p_dens <- ppc_dens_overlay(y, y_rep[1:100, ]) +
  labs(
    title = "Posterior Predictive: Density Overlay",
    subtitle = "100 replications (light blue) vs observed (dark blue)",
    x = "Weekly Cases"
  ) +
  coord_cartesian(xlim = c(0, quantile(y, 0.99) * 1.5))

# --- PPC 2: Histogram of replicated vs observed ---
p_hist <- ppc_hist(y, y_rep[1:8, ], binwidth = 50) +
  labs(
    title = "Posterior Predictive: Histograms",
    subtitle = "y (observed) vs 8 replicated datasets"
  )

# --- PPC 3: Summary statistics ---
# Mean
p_stat_mean <- ppc_stat(y, y_rep, stat = "mean") +
  labs(title = "Posterior Predictive: Mean")

# SD
p_stat_sd <- ppc_stat(y, y_rep, stat = "sd") +
  labs(title = "Posterior Predictive: SD")

# Max
p_stat_max <- ppc_stat(y, y_rep, stat = "max") +
  labs(title = "Posterior Predictive: Maximum")

# Min
p_stat_min <- ppc_stat(y, y_rep, stat = "min") +
  labs(title = "Posterior Predictive: Minimum")

# --- PPC 4: Time series with intervals ---
y_rep_summary <- tibble(
  week = 1:N_model,
  y_obs = y,
  median = apply(y_rep, 2, median),
  q05 = apply(y_rep, 2, quantile, 0.05),
  q25 = apply(y_rep, 2, quantile, 0.25),
  q75 = apply(y_rep, 2, quantile, 0.75),
  q95 = apply(y_rep, 2, quantile, 0.95)
)

# Add dates if available
if (!is.null(metadata$dates)) {
  dates <- metadata$dates[(S + 1):N]
  y_rep_summary$date <- dates
  x_var <- "date"
  x_lab <- "Date"
} else {
  x_var <- "week"
  x_lab <- "Week"
}

p_ts <- ggplot(y_rep_summary, aes_string(x = x_var)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "steelblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "steelblue", alpha = 0.4) +
  geom_line(aes(y = median), color = "steelblue", linewidth = 0.8) +
  geom_point(aes(y = y_obs), color = "darkred", size = 0.8, alpha = 0.7) +
  labs(
    title = "Posterior Predictive: Time Series",
    subtitle = "Observed (red points) vs predicted intervals (blue)",
    x = x_lab,
    y = "Weekly Cases"
  )

# --- PPC 5: Intervals (pointwise coverage) ---
p_intervals <- ppc_intervals(y, y_rep, x = 1:N_model, prob = 0.5, prob_outer = 0.9) +
  labs(
    title = "Posterior Predictive: Intervals",
    subtitle = "50% and 90% intervals; points are observed",
    x = "Week Index",
    y = "Weekly Cases"
  )

# --- PPC 6: Scatter plot (observed vs predicted median) ---
p_scatter <- ggplot(y_rep_summary, aes(x = y_obs, y = median)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_point(alpha = 0.5, color = "steelblue") +
  labs(
    title = "Posterior Predictive: Calibration",
    subtitle = "Observed vs predicted median",
    x = "Observed Cases",
    y = "Predicted Median"
  ) +
  coord_equal()

# ==============================================================================
# 4. QUANTITATIVE CHECKS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("POSTERIOR PREDICTIVE SUMMARIES\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Check 1: Coverage of intervals
in_50_interval <- y >= y_rep_summary$q25 & y <= y_rep_summary$q75
in_90_interval <- y >= y_rep_summary$q05 & y <= y_rep_summary$q95

coverage_50 <- mean(in_50_interval)
coverage_90 <- mean(in_90_interval)

cat("Interval Coverage:\n")
cat(sprintf("  50%% intervals contain observed: %.1f%% (nominal: 50%%)\n",
            100 * coverage_50))
cat(sprintf("  90%% intervals contain observed: %.1f%% (nominal: 90%%)\n",
            100 * coverage_90))

# Check 2: Summary statistics
obs_mean <- mean(y)
obs_sd <- sd(y)
obs_max <- max(y)
obs_min <- min(y)

rep_means <- apply(y_rep, 1, mean)
rep_sds <- apply(y_rep, 1, sd)
rep_maxs <- apply(y_rep, 1, max)
rep_mins <- apply(y_rep, 1, min)

cat("\nSummary Statistics (observed vs replicated):\n")
cat(sprintf("  Mean:   observed=%.0f, replicated median=%.0f (95%% CI: %.0f-%.0f)\n",
            obs_mean, median(rep_means),
            quantile(rep_means, 0.025), quantile(rep_means, 0.975)))
cat(sprintf("  SD:     observed=%.0f, replicated median=%.0f (95%% CI: %.0f-%.0f)\n",
            obs_sd, median(rep_sds),
            quantile(rep_sds, 0.025), quantile(rep_sds, 0.975)))
cat(sprintf("  Max:    observed=%.0f, replicated median=%.0f (95%% CI: %.0f-%.0f)\n",
            obs_max, median(rep_maxs),
            quantile(rep_maxs, 0.025), quantile(rep_maxs, 0.975)))
cat(sprintf("  Min:    observed=%.0f, replicated median=%.0f (95%% CI: %.0f-%.0f)\n",
            obs_min, median(rep_mins),
            quantile(rep_mins, 0.025), quantile(rep_mins, 0.975)))

# Check 3: Bayesian p-values (proportion of replications more extreme than observed)
p_mean <- mean(rep_means >= obs_mean)
p_sd <- mean(rep_sds >= obs_sd)
p_max <- mean(rep_maxs >= obs_max)
p_min <- mean(rep_mins <= obs_min)

cat("\nBayesian p-values (should be near 0.5 for well-calibrated model):\n")
cat(sprintf("  Mean: %.3f\n", p_mean))
cat(sprintf("  SD:   %.3f\n", p_sd))
cat(sprintf("  Max:  %.3f\n", p_max))
cat(sprintf("  Min:  %.3f\n", p_min))

# Check 4: Root mean squared error
rmse_per_draw <- apply(y_rep, 1, function(y_r) sqrt(mean((y - y_r)^2)))
cat(sprintf("\nRMSE: median=%.0f (95%% CI: %.0f-%.0f)\n",
            median(rmse_per_draw),
            quantile(rmse_per_draw, 0.025),
            quantile(rmse_per_draw, 0.975)))

# Check 5: Mean absolute error
mae_per_draw <- apply(y_rep, 1, function(y_r) mean(abs(y - y_r)))
cat(sprintf("MAE:  median=%.0f (95%% CI: %.0f-%.0f)\n",
            median(mae_per_draw),
            quantile(mae_per_draw, 0.025),
            quantile(mae_per_draw, 0.975)))

# ==============================================================================
# 5. SAVE PLOTS
# ==============================================================================

cat("\nSaving visualizations...\n")

# Combined summary plot
p_summary <- (p_dens | p_scatter) / (p_stat_mean | p_stat_sd)
ggsave("../results/figures/posterior_predictive_summary.png", p_summary,
       width = 12, height = 10, dpi = 150)
cat("  Saved: results/figures/posterior_predictive_summary.png\n")

# Time series plot
ggsave("../results/figures/posterior_predictive_timeseries.png", p_ts,
       width = 12, height = 6, dpi = 150)
cat("  Saved: results/figures/posterior_predictive_timeseries.png\n")

# Intervals plot
ggsave("../results/figures/posterior_predictive_intervals.png", p_intervals,
       width = 12, height = 6, dpi = 150)
cat("  Saved: results/figures/posterior_predictive_intervals.png\n")

# Density overlay
ggsave("../results/figures/posterior_predictive_density.png", p_dens,
       width = 8, height = 6, dpi = 150)
cat("  Saved: results/figures/posterior_predictive_density.png\n")

# Statistics plots
p_stats <- (p_stat_mean | p_stat_sd) / (p_stat_max | p_stat_min)
ggsave("../results/figures/posterior_predictive_stats.png", p_stats,
       width = 10, height = 8, dpi = 150)
cat("  Saved: results/figures/posterior_predictive_stats.png\n")

# Histograms
ggsave("../results/figures/posterior_predictive_histograms.png", p_hist,
       width = 10, height = 8, dpi = 150)
cat("  Saved: results/figures/posterior_predictive_histograms.png\n")

# ==============================================================================
# 6. ASSESSMENT
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("POSTERIOR PREDICTIVE ASSESSMENT\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

checks <- list(
  coverage_50 = abs(coverage_50 - 0.5) < 0.15,  # Within 15% of nominal
  coverage_90 = abs(coverage_90 - 0.9) < 0.1,   # Within 10% of nominal
  mean_captured = obs_mean >= quantile(rep_means, 0.025) &
                  obs_mean <= quantile(rep_means, 0.975),
  sd_captured = obs_sd >= quantile(rep_sds, 0.025) &
                obs_sd <= quantile(rep_sds, 0.975),
  p_values_ok = all(c(p_mean, p_sd) > 0.05 & c(p_mean, p_sd) < 0.95)
)

cat("Checks:\n")
cat(sprintf("  [%s] 50%% interval coverage reasonable (35-65%%)\n",
            ifelse(checks$coverage_50, "PASS", "WARN")))
cat(sprintf("  [%s] 90%% interval coverage reasonable (80-100%%)\n",
            ifelse(checks$coverage_90, "PASS", "WARN")))
cat(sprintf("  [%s] Observed mean within 95%% CI of replicated means\n",
            ifelse(checks$mean_captured, "PASS", "WARN")))
cat(sprintf("  [%s] Observed SD within 95%% CI of replicated SDs\n",
            ifelse(checks$sd_captured, "PASS", "WARN")))
cat(sprintf("  [%s] Bayesian p-values not extreme (0.05-0.95)\n",
            ifelse(checks$p_values_ok, "PASS", "WARN")))

all_pass <- all(unlist(checks))

cat(sprintf("\nOverall: %s\n",
            ifelse(all_pass, "MODEL FITS DATA WELL",
                   "REVIEW MODEL FIT - some checks show misfit")))

# ==============================================================================
# 7. RESIDUAL ANALYSIS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("RESIDUAL ANALYSIS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Compute standardized residuals
y_rep_median <- apply(y_rep, 2, median)
y_rep_sd <- apply(y_rep, 2, sd)
std_resid <- (y - y_rep_median) / y_rep_sd

cat(sprintf("Standardized residuals: mean=%.2f, sd=%.2f\n",
            mean(std_resid), sd(std_resid)))

# Check for autocorrelation in residuals
acf_resid <- acf(std_resid, lag.max = 10, plot = FALSE)
cat(sprintf("Residual autocorrelation (lag 1): %.3f\n", acf_resid$acf[2]))

# Identify outliers (|z| > 2)
outliers <- which(abs(std_resid) > 2)
cat(sprintf("Outliers (|z| > 2): %d observations (%.1f%%)\n",
            length(outliers), 100 * length(outliers) / length(std_resid)))

# Plot residuals
resid_df <- tibble(
  week = 1:N_model,
  std_resid = std_resid,
  y_obs = y,
  y_pred = y_rep_median
)

if (!is.null(metadata$dates)) {
  resid_df$date <- metadata$dates[(S + 1):N]
}

p_resid_ts <- ggplot(resid_df, aes_string(x = ifelse(!is.null(metadata$dates), "date", "week"),
                                          y = "std_resid")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "red") +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_smooth(method = "loess", se = FALSE, color = "darkblue", linewidth = 0.8) +
  labs(
    title = "Standardized Residuals Over Time",
    subtitle = "Dashed lines at +/- 2 standard deviations",
    x = ifelse(!is.null(metadata$dates), "Date", "Week"),
    y = "Standardized Residual"
  )

p_resid_hist <- ggplot(resid_df, aes(x = std_resid)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "steelblue", alpha = 0.7) +
  stat_function(fun = dnorm, color = "red", linewidth = 1) +
  labs(
    title = "Distribution of Standardized Residuals",
    subtitle = "Red line: standard normal",
    x = "Standardized Residual",
    y = "Density"
  )

p_resid_qq <- ggplot(resid_df, aes(sample = std_resid)) +
  stat_qq(color = "steelblue", alpha = 0.5) +
  stat_qq_line(color = "red") +
  labs(
    title = "Q-Q Plot of Standardized Residuals",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )

p_resid_combined <- (p_resid_ts) / (p_resid_hist | p_resid_qq)
ggsave("../results/figures/posterior_predictive_residuals.png", p_resid_combined,
       width = 12, height = 10, dpi = 150)
cat("\n  Saved: results/figures/posterior_predictive_residuals.png\n")

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("POSTERIOR PREDICTIVE CHECKS COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\nFiles saved:
  results/figures/posterior_predictive_summary.png
  results/figures/posterior_predictive_timeseries.png
  results/figures/posterior_predictive_intervals.png
  results/figures/posterior_predictive_density.png
  results/figures/posterior_predictive_stats.png
  results/figures/posterior_predictive_histograms.png
  results/figures/posterior_predictive_residuals.png
\n")
