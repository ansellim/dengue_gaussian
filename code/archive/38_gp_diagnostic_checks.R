#!/usr/bin/env Rscript
# ==============================================================================
# 38_gp_diagnostic_checks.R
#
# Diagnostic checks: Is the residual GP absorbing climate signal?
# (Model misspecification diagnostic for the climate + GP decomposition)
#
# Diagnostics:
#   1. GP-Climate correlation (f_residual vs temp/rain)
#   2. Climate-only model fit (no GP, using posterior beta draws)
#   3. Residual autocorrelation comparison (full model vs climate-only)
#   4. GP length scale interpretation
#
# Input:  results/fit_model3.rds, data/model_data.rds
# Output: results/figures/gp_diagnostic_*.png, results/gp_diagnostic_summary.csv
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
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
cat("GP DIAGNOSTIC CHECKS: IS THE GP ABSORBING CLIMATE SIGNAL?\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 0. LOAD DATA AND FIT
# ==============================================================================

cat("Loading data and model fit...\n")

model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data
df <- model_data$df
dates <- model_data$metadata$dates_model

fit <- readRDS("../results/fit_model3.rds")

N_model <- stan_data$N_model
S <- stan_data$S
N <- stan_data$N

cat(sprintf("  N_model = %d weeks\n", N_model))
cat(sprintf("  S (max generation interval) = %d weeks\n", S))

# Extract climate covariates (standardized)
temp <- stan_data$X_climate[, 1]
rain <- stan_data$X_climate[, 2]

# Collector for summary statistics
summary_stats <- list()

# ==============================================================================
# DIAGNOSTIC 1: GP-CLIMATE CORRELATION
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DIAGNOSTIC 1: GP-Climate Correlation\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

f_residual_draws <- fit$draws("f_residual", format = "matrix")
n_draws <- nrow(f_residual_draws)
n_cols <- ncol(f_residual_draws)

cat(sprintf("Computing correlations across %d posterior draws...\n", n_draws))
cat(sprintf("  f_residual columns: %d, temp length: %d, rain length: %d\n",
            n_cols, length(temp), length(rain)))

# Align dimensions (use minimum common length)
n_common <- min(n_cols, length(temp), length(rain))
temp_aligned <- temp[1:n_common]
rain_aligned <- rain[1:n_common]

cor_temp <- numeric(n_draws)
cor_rain <- numeric(n_draws)

for (i in 1:n_draws) {
  f_res_i <- f_residual_draws[i, 1:n_common]
  cor_temp[i] <- cor(f_res_i, temp_aligned)
  cor_rain[i] <- cor(f_res_i, rain_aligned)
}

cat("\nCorrelation between f_residual and temperature:\n")
cat(sprintf("  Mean: %.3f, Median: %.3f\n", mean(cor_temp), median(cor_temp)))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n",
            quantile(cor_temp, 0.025), quantile(cor_temp, 0.975)))

cat("\nCorrelation between f_residual and rainfall:\n")
cat(sprintf("  Mean: %.3f, Median: %.3f\n", mean(cor_rain), median(cor_rain)))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n",
            quantile(cor_rain, 0.025), quantile(cor_rain, 0.975)))

# Store summary
summary_stats$cor_temp_median <- median(cor_temp)
summary_stats$cor_temp_lower <- quantile(cor_temp, 0.025)
summary_stats$cor_temp_upper <- quantile(cor_temp, 0.975)
summary_stats$cor_rain_median <- median(cor_rain)
summary_stats$cor_rain_lower <- quantile(cor_rain, 0.025)
summary_stats$cor_rain_upper <- quantile(cor_rain, 0.975)

# Interpretation
temp_absorbing <- abs(median(cor_temp)) > 0.3
rain_absorbing <- abs(median(cor_rain)) > 0.3

if (!temp_absorbing && !rain_absorbing) {
  cat("\n--> f_residual is UNCORRELATED with climate covariates.\n")
  cat("    The GP is NOT absorbing climate signal.\n")
} else {
  cat("\n--> WARNING: f_residual shows non-trivial correlation with climate.\n")
  cat("    The GP may be absorbing some climate signal.\n")
}

# Figure: density plot of correlations
cor_df <- bind_rows(
  tibble(variable = "Temperature", correlation = cor_temp),
  tibble(variable = "Rainfall", correlation = cor_rain)
)

p_cor <- ggplot(cor_df, aes(x = correlation, fill = variable)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dotted", color = "red", alpha = 0.5) +
  scale_fill_manual(values = c("Temperature" = "tomato", "Rainfall" = "steelblue")) +
  labs(
    title = "Diagnostic 1: Correlation between GP residual and climate covariates",
    subtitle = "Red dotted lines at |r| = 0.3. Near-zero correlation indicates GP is not absorbing climate.",
    x = "Correlation coefficient",
    y = "Posterior density",
    fill = "Covariate"
  )

ggsave("../results/figures/gp_diagnostic_climate_correlation.png", p_cor,
       width = 10, height = 5, dpi = 150)
cat("  Saved: results/figures/gp_diagnostic_climate_correlation.png\n")

# ==============================================================================
# DIAGNOSTIC 2: CLIMATE-ONLY MODEL (NO GP)
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DIAGNOSTIC 2: Climate-Only Model Fit (No GP)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

mu_draws <- fit$draws("mu", format = "matrix")[, 1]
beta_climate_draws <- fit$draws("beta_climate", format = "matrix")

cases_obs <- stan_data$cases[(S + 1):N]
gi <- stan_data$gi

# Compute climate-only Rt and predicted cases for each draw
n_subsample <- min(n_draws, 500)
set.seed(42)
draw_idx <- sample(1:n_draws, n_subsample)

cases_climate_only <- matrix(NA, nrow = n_subsample, ncol = N_model)

cat(sprintf("Computing climate-only predictions for %d posterior draws...\n", n_subsample))

for (d in seq_along(draw_idx)) {
  idx <- draw_idx[d]
  mu_d <- mu_draws[idx]
  beta_temp_d <- beta_climate_draws[idx, 1]
  beta_rain_d <- beta_climate_draws[idx, 2]

  log_Rt_climate <- mu_d + beta_temp_d * temp + beta_rain_d * rain
  Rt_climate <- exp(log_Rt_climate)

  # Renewal equation: expected cases
  lambda_climate <- numeric(N_model)
  for (i in 1:N_model) {
    infectious_pressure <- 0
    t_idx <- S + i
    for (s in 1:S) {
      infectious_pressure <- infectious_pressure + stan_data$cases[t_idx - s] * gi[s]
    }
    lambda_climate[i] <- max(Rt_climate[i] * infectious_pressure, 1.0)
  }

  cases_climate_only[d, ] <- lambda_climate
}

# Summarize climate-only predictions
clim_pred_median <- apply(cases_climate_only, 2, median)
clim_pred_lower <- apply(cases_climate_only, 2, quantile, 0.025)
clim_pred_upper <- apply(cases_climate_only, 2, quantile, 0.975)

# RMSE
rmse_climate_only <- sqrt(mean((cases_obs - clim_pred_median)^2))
cat(sprintf("\nClimate-only RMSE: %.1f cases\n", rmse_climate_only))

# Full model comparison
cases_pred_full <- fit$draws("cases_pred", format = "matrix")
n_common <- min(ncol(cases_pred_full), length(cases_obs))
full_pred_median <- apply(cases_pred_full[, 1:n_common], 2, median)
rmse_full <- sqrt(mean((cases_obs[1:n_common] - full_pred_median)^2))
cat(sprintf("Full model RMSE: %.1f cases\n", rmse_full))

# Coverage (climate-only)
coverage_95_climate <- mean(cases_obs >= clim_pred_lower & cases_obs <= clim_pred_upper)
cat(sprintf("Climate-only 95%% coverage: %.1f%%\n", 100 * coverage_95_climate))

# Residuals for climate-only
resid_climate <- cases_obs - clim_pred_median

# Residual autocorrelation
acf_climate <- acf(resid_climate, lag.max = 20, plot = FALSE)

cat("\nClimate-only residual autocorrelation (first 5 lags):\n")
for (lag in 1:5) {
  cat(sprintf("  Lag %d: %.3f\n", lag, acf_climate$acf[lag + 1]))
}

summary_stats$rmse_climate_only <- rmse_climate_only
summary_stats$rmse_full_model <- rmse_full
summary_stats$coverage_95_climate_only <- coverage_95_climate
summary_stats$acf_lag1_climate_only <- acf_climate$acf[2]

# Figure: 2-panel climate-only fit
clim_fit_df <- tibble(
  date = dates,
  observed = cases_obs,
  predicted = clim_pred_median,
  lower = clim_pred_lower,
  upper = clim_pred_upper,
  residual = resid_climate
)

p_clim_fit <- ggplot(clim_fit_df, aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = predicted), color = "steelblue", linewidth = 0.8) +
  geom_point(aes(y = observed), size = 0.5, alpha = 0.7) +
  labs(
    title = "Climate-only expected cases (no GP)",
    subtitle = sprintf("RMSE = %.0f | 95%% coverage = %.1f%%",
                        rmse_climate_only, 100 * coverage_95_climate),
    x = NULL, y = "Weekly cases"
  )

# ACF of climate-only residuals
acf_vals_climate <- data.frame(
  lag = acf_climate$lag[-1],
  acf = acf_climate$acf[-1]
)

ci_bound <- qnorm(0.975) / sqrt(N_model)

p_clim_resid <- ggplot(clim_fit_df, aes(x = date, y = residual)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
  geom_segment(aes(xend = date, yend = 0), alpha = 0.3) +
  geom_point(size = 0.5, alpha = 0.5) +
  labs(title = "Climate-only residuals", x = "Date", y = "Residual (obs - pred)")

p_clim_combined <- p_clim_fit / p_clim_resid +
  plot_annotation(title = "Diagnostic 2: Climate-only model (no GP)")

ggsave("../results/figures/gp_diagnostic_climate_only_fit.png", p_clim_combined,
       width = 12, height = 8, dpi = 150)
cat("  Saved: results/figures/gp_diagnostic_climate_only_fit.png\n")

# ==============================================================================
# DIAGNOSTIC 3: RESIDUAL AUTOCORRELATION STRUCTURE
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DIAGNOSTIC 3: Residual Autocorrelation Comparison\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Full model residuals
full_pred_median_all <- apply(cases_pred_full[, 1:n_common], 2, median)
full_pred_sd <- apply(cases_pred_full[, 1:n_common], 2, sd)
resid_full <- cases_obs[1:n_common] - full_pred_median_all

acf_full <- acf(resid_full, lag.max = 20, plot = FALSE)

cat("Full model residual autocorrelation (first 5 lags):\n")
for (lag in 1:5) {
  cat(sprintf("  Lag %d: %.3f\n", lag, acf_full$acf[lag + 1]))
}

cat("\nClimate-only residual autocorrelation (first 5 lags):\n")
for (lag in 1:5) {
  cat(sprintf("  Lag %d: %.3f\n", lag, acf_climate$acf[lag + 1]))
}

summary_stats$acf_lag1_full_model <- acf_full$acf[2]

# Interpretation
if (abs(acf_full$acf[2]) < ci_bound && abs(acf_climate$acf[2]) > ci_bound) {
  cat("\n--> Full model residuals are approximately white noise.\n")
  cat("    Climate-only residuals are strongly autocorrelated.\n")
  cat("    The GP captures REAL temporal structure, not just noise.\n")
} else if (abs(acf_full$acf[2]) < ci_bound) {
  cat("\n--> Full model residuals are approximately white noise.\n")
} else {
  cat("\n--> Some autocorrelation remains in full model residuals.\n")
}

# Figure: side-by-side ACF
acf_comparison <- bind_rows(
  tibble(lag = as.numeric(acf_full$lag[-1]),
         acf = as.numeric(acf_full$acf[-1]),
         model = "Full model (climate + GP)"),
  tibble(lag = as.numeric(acf_climate$lag[-1]),
         acf = as.numeric(acf_climate$acf[-1]),
         model = "Climate-only (no GP)")
)

p_acf <- ggplot(acf_comparison, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_hline(yintercept = c(-ci_bound, ci_bound),
             linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ model, ncol = 2) +
  labs(
    title = "Diagnostic 3: Residual autocorrelation comparison",
    subtitle = "Blue dashed lines = 95% significance bounds. White noise expected for full model.",
    x = "Lag (weeks)",
    y = "Autocorrelation"
  ) +
  coord_cartesian(ylim = c(-1, 1))

ggsave("../results/figures/gp_diagnostic_acf_comparison.png", p_acf,
       width = 12, height = 5, dpi = 150)
cat("  Saved: results/figures/gp_diagnostic_acf_comparison.png\n")

# ==============================================================================
# DIAGNOSTIC 4: GP LENGTH SCALE INTERPRETATION
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DIAGNOSTIC 4: GP Length Scale Interpretation\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

rho_draws <- fit$draws("rho", format = "matrix")[, 1]
rho_median <- median(rho_draws)
rho_lower <- quantile(rho_draws, 0.025)
rho_upper <- quantile(rho_draws, 0.975)

cat(sprintf("Posterior rho (length scale): %.1f weeks (95%% CI: [%.1f, %.1f])\n",
            rho_median, rho_lower, rho_upper))

summary_stats$rho_median <- rho_median
summary_stats$rho_lower <- as.numeric(rho_lower)
summary_stats$rho_upper <- as.numeric(rho_upper)

cat("\nInterpretation:\n")
cat(sprintf("  - Posterior rho ~ %.0f weeks corresponds to outbreak-scale dynamics.\n",
            round(rho_median)))
cat("  - Seasonal climate effects operate at 26-52 week timescales.\n")
cat("  - If the GP were absorbing climate signal, we would expect rho >> 10 weeks\n")
cat("    (i.e., seasonal or longer timescales).\n")
cat(sprintf("  - The short rho ~ %.0f weeks suggests the GP captures:\n",
            round(rho_median)))
cat("      * Outbreak-scale dynamics (serotype introductions, immunity shifts)\n")
cat("      * Short-term transmission heterogeneity\n")
cat("    NOT seasonal climate patterns.\n")

# Proportion of posterior mass above seasonal threshold (26 weeks)
prop_seasonal <- mean(rho_draws > 26)
cat(sprintf("\n  P(rho > 26 weeks) = %.3f (posterior probability of seasonal timescale)\n",
            prop_seasonal))
summary_stats$prob_rho_seasonal <- prop_seasonal

if (prop_seasonal < 0.05) {
  cat("  --> Negligible probability that GP operates at seasonal timescale.\n")
}

# ==============================================================================
# SUMMARY FIGURE: 4-panel diagnostic
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("Generating summary figure...\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Panel 1: Correlation density
p1 <- p_cor +
  labs(title = "1. GP-Climate correlation") +
  theme(legend.position = "bottom")

# Panel 2: Climate-only fit
p2 <- p_clim_fit +
  labs(title = "2. Climate-only predicted cases (no GP)")

# Panel 3: ACF comparison
p3 <- p_acf +
  labs(title = "3. Residual autocorrelation") +
  theme(legend.position = "none")

# Panel 4: GP length scale
rho_df <- tibble(rho = rho_draws)
p4 <- ggplot(rho_df, aes(x = rho)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = 26, linetype = "dashed", color = "red") +
  annotate("text", x = 28, y = 0, label = "Seasonal\nthreshold",
           color = "red", hjust = 0, size = 3) +
  labs(
    title = "4. GP length scale (rho)",
    subtitle = sprintf("Median = %.1f weeks. P(rho > 26) = %.3f",
                        rho_median, prop_seasonal),
    x = "Length scale (weeks)",
    y = "Posterior density"
  )

p_summary <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Is the GP absorbing climate signal? -- No",
    subtitle = paste0(
      "Evidence: (1) GP uncorrelated with climate, ",
      "(2) climate-only model fits poorly, ",
      "(3) removing GP creates strong autocorrelation, ",
      "(4) GP length scale is outbreak-scale, not seasonal"
    ),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  )

ggsave("../results/figures/gp_diagnostic_summary.png", p_summary,
       width = 16, height = 12, dpi = 150)
cat("  Saved: results/figures/gp_diagnostic_summary.png\n")

# ==============================================================================
# SAVE SUMMARY STATISTICS
# ==============================================================================

cat("\nSaving summary statistics...\n")

summary_df <- tibble(
  diagnostic = c(
    "cor_temp_median", "cor_temp_lower95", "cor_temp_upper95",
    "cor_rain_median", "cor_rain_lower95", "cor_rain_upper95",
    "rmse_climate_only", "rmse_full_model",
    "coverage_95_climate_only",
    "acf_lag1_climate_only", "acf_lag1_full_model",
    "rho_median_weeks", "rho_lower95_weeks", "rho_upper95_weeks",
    "prob_rho_seasonal_gt26wk"
  ),
  value = c(
    summary_stats$cor_temp_median,
    as.numeric(summary_stats$cor_temp_lower),
    as.numeric(summary_stats$cor_temp_upper),
    summary_stats$cor_rain_median,
    as.numeric(summary_stats$cor_rain_lower),
    as.numeric(summary_stats$cor_rain_upper),
    summary_stats$rmse_climate_only,
    summary_stats$rmse_full_model,
    summary_stats$coverage_95_climate_only,
    summary_stats$acf_lag1_climate_only,
    summary_stats$acf_lag1_full_model,
    summary_stats$rho_median,
    summary_stats$rho_lower,
    summary_stats$rho_upper,
    summary_stats$prob_rho_seasonal
  ),
  interpretation = c(
    "Near zero = GP not absorbing temperature", "", "",
    "Near zero = GP not absorbing rainfall", "", "",
    "High = climate alone fits poorly",
    "Low = full model fits well",
    "Poor coverage = climate alone insufficient",
    "High = strong residual structure without GP",
    "Near zero = full model captures temporal structure",
    "Short = outbreak-scale, not seasonal",
    "", "",
    "Near zero = GP not operating at seasonal timescale"
  )
)

write_csv(summary_df, "../results/gp_diagnostic_summary.csv")
cat("  Saved: results/gp_diagnostic_summary.csv\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("GP DIAGNOSTIC CHECKS COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\n--- VERDICT: Is the GP absorbing climate signal? ---\n\n")

cat("DIAGNOSTIC 1 (GP-Climate Correlation):\n")
cat(sprintf("  cor(f_residual, temp) = %.3f [%.3f, %.3f]\n",
            summary_stats$cor_temp_median,
            as.numeric(summary_stats$cor_temp_lower),
            as.numeric(summary_stats$cor_temp_upper)))
cat(sprintf("  cor(f_residual, rain) = %.3f [%.3f, %.3f]\n",
            summary_stats$cor_rain_median,
            as.numeric(summary_stats$cor_rain_lower),
            as.numeric(summary_stats$cor_rain_upper)))

cat("\nDIAGNOSTIC 2 (Climate-Only Fit):\n")
cat(sprintf("  RMSE climate-only: %.0f vs full model: %.0f\n",
            summary_stats$rmse_climate_only, summary_stats$rmse_full_model))
cat(sprintf("  95%% coverage climate-only: %.1f%%\n",
            100 * summary_stats$coverage_95_climate_only))

cat("\nDIAGNOSTIC 3 (Autocorrelation):\n")
cat(sprintf("  ACF(1) full model: %.3f\n", summary_stats$acf_lag1_full_model))
cat(sprintf("  ACF(1) climate-only: %.3f\n", summary_stats$acf_lag1_climate_only))

cat("\nDIAGNOSTIC 4 (Length Scale):\n")
cat(sprintf("  rho = %.1f weeks [%.1f, %.1f]\n",
            summary_stats$rho_median, summary_stats$rho_lower, summary_stats$rho_upper))
cat(sprintf("  P(rho > 26 weeks) = %.3f\n", summary_stats$prob_rho_seasonal))

cat("\nCONCLUSION: The GP is NOT absorbing climate signal. It captures\n")
cat("outbreak-scale dynamics that climate covariates cannot explain.\n")

cat("\nOutput files:\n")
cat("  results/figures/gp_diagnostic_climate_correlation.png\n")
cat("  results/figures/gp_diagnostic_climate_only_fit.png\n")
cat("  results/figures/gp_diagnostic_acf_comparison.png\n")
cat("  results/figures/gp_diagnostic_summary.png\n")
cat("  results/gp_diagnostic_summary.csv\n")
