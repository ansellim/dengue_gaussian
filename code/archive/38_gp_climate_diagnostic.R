#!/usr/bin/env Rscript
# ==============================================================================
# 38_gp_climate_diagnostic.R
#
# Diagnostic: Does the GP residual absorb climate signal?
#
# If the GP is simply soaking up climate information that the linear covariates
# cannot capture, f_residual should be correlated with temperature and rainfall.
# This script tests that hypothesis via correlation analysis and CCF.
#
# Input:  results/fit_model3.rds, data/model_data.rds
# Output: results/figures/gp_diagnostic_scatter.png
#         results/figures/gp_diagnostic_ccf.png
#         results/gp_climate_diagnostic.csv
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
cat("GP RESIDUAL - CLIMATE ABSORPTION DIAGNOSTIC\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA AND FIT
# ==============================================================================

cat("Loading data and model fit...\n")

model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data

fit <- readRDS("../results/fit_model3.rds")

N_model <- stan_data$N_model
dates <- model_data$metadata$dates_model

cat(sprintf("  N_model = %d weeks\n", N_model))

# ==============================================================================
# 2. EXTRACT f_residual AND CLIMATE COVARIATES
# ==============================================================================

cat("\nExtracting posteriors...\n")

# f_residual: full posterior draws (n_draws x N_model)
f_residual_draws <- fit$draws("f_residual", format = "matrix")
n_draws <- nrow(f_residual_draws)

# Posterior median of f_residual
f_residual_median <- apply(f_residual_draws, 2, median)

# Climate covariates (standardized, from stan_data)
temp <- stan_data$X_climate[, 1]
rain <- stan_data$X_climate[, 2]

cat(sprintf("  Posterior draws: %d\n", n_draws))
cat(sprintf("  f_residual range: [%.3f, %.3f]\n",
            min(f_residual_median), max(f_residual_median)))

# ==============================================================================
# 3. POINT ESTIMATE CORRELATIONS
# ==============================================================================

cat("\n--- Point Estimate Correlations (posterior median f_residual) ---\n\n")

cor_temp_pearson  <- cor(f_residual_median, temp, method = "pearson")
cor_temp_spearman <- cor(f_residual_median, temp, method = "spearman")
cor_rain_pearson  <- cor(f_residual_median, rain, method = "pearson")
cor_rain_spearman <- cor(f_residual_median, rain, method = "spearman")

cat(sprintf("  f_residual vs Temperature:  Pearson = %.4f, Spearman = %.4f\n",
            cor_temp_pearson, cor_temp_spearman))
cat(sprintf("  f_residual vs Rainfall:     Pearson = %.4f, Spearman = %.4f\n",
            cor_rain_pearson, cor_rain_spearman))

# ==============================================================================
# 4. POSTERIOR UNCERTAINTY ON CORRELATIONS (200 subsampled draws)
# ==============================================================================

cat("\nComputing correlation uncertainty (200 posterior draws)...\n")

set.seed(123)
n_sub <- min(200, n_draws)
draw_idx <- sample(n_draws, n_sub, replace = FALSE)

cor_temp_pearson_draws  <- numeric(n_sub)
cor_temp_spearman_draws <- numeric(n_sub)
cor_rain_pearson_draws  <- numeric(n_sub)
cor_rain_spearman_draws <- numeric(n_sub)

for (i in seq_len(n_sub)) {
  f_res_i <- f_residual_draws[draw_idx[i], ]
  cor_temp_pearson_draws[i]  <- cor(f_res_i, temp, method = "pearson")
  cor_temp_spearman_draws[i] <- cor(f_res_i, temp, method = "spearman")
  cor_rain_pearson_draws[i]  <- cor(f_res_i, rain, method = "pearson")
  cor_rain_spearman_draws[i] <- cor(f_res_i, rain, method = "spearman")
}

cat(sprintf("  Temp Spearman:  median = %.4f, 95%% CI = [%.4f, %.4f]\n",
            median(cor_temp_spearman_draws),
            quantile(cor_temp_spearman_draws, 0.025),
            quantile(cor_temp_spearman_draws, 0.975)))
cat(sprintf("  Rain Spearman:  median = %.4f, 95%% CI = [%.4f, %.4f]\n",
            median(cor_rain_spearman_draws),
            quantile(cor_rain_spearman_draws, 0.025),
            quantile(cor_rain_spearman_draws, 0.975)))

# ==============================================================================
# 5. CROSS-CORRELATION (CCF) AT LAGS 0-12 WEEKS
# ==============================================================================

cat("\nComputing cross-correlations at lags 0-12...\n")

max_lag <- 12

ccf_temp <- ccf(temp, f_residual_median, lag.max = max_lag, plot = FALSE)
ccf_rain <- ccf(rain, f_residual_median, lag.max = max_lag, plot = FALSE)

# Extract non-negative lags only (climate leading residual)
lag_idx <- which(ccf_temp$lag >= 0)
ccf_df <- bind_rows(
  tibble(lag = as.numeric(ccf_temp$lag[lag_idx]),
         ccf = as.numeric(ccf_temp$acf[lag_idx]),
         variable = "Temperature"),
  tibble(lag = as.numeric(ccf_rain$lag[lag_idx]),
         ccf = as.numeric(ccf_rain$acf[lag_idx]),
         variable = "Rainfall")
)

cat("  CCF at lag 0:\n")
cat(sprintf("    Temperature: %.4f\n", ccf_df$ccf[ccf_df$variable == "Temperature" & ccf_df$lag == 0]))
cat(sprintf("    Rainfall:    %.4f\n", ccf_df$ccf[ccf_df$variable == "Rainfall" & ccf_df$lag == 0]))

# ==============================================================================
# 6. SCATTER PLOT: f_residual vs climate
# ==============================================================================

cat("\nGenerating scatter plot...\n")

scatter_df <- tibble(
  f_residual = f_residual_median,
  temperature = temp,
  rainfall = rain
)

p_temp <- ggplot(scatter_df, aes(x = temperature, y = f_residual)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("rho == %.3f", cor_temp_spearman),
           parse = TRUE, hjust = 1.1, vjust = 1.5, size = 4.5) +
  labs(x = "Temperature (standardized)",
       y = expression(f[residual] ~ "(posterior median)"),
       title = "GP residual vs Temperature")

p_rain <- ggplot(scatter_df, aes(x = rainfall, y = f_residual)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("rho == %.3f", cor_rain_spearman),
           parse = TRUE, hjust = 1.1, vjust = 1.5, size = 4.5) +
  labs(x = "Rainfall (standardized)",
       y = expression(f[residual] ~ "(posterior median)"),
       title = "GP residual vs Rainfall")

p_scatter <- p_temp + p_rain +
  plot_annotation(
    title = "Diagnostic: Is the GP Residual Absorbing Climate Signal?",
    subtitle = "Spearman rho annotated; regression line shown for visual reference"
  )

ggsave("../results/figures/gp_diagnostic_scatter.png", p_scatter,
       width = 12, height = 5, dpi = 150)
cat("  Saved: results/figures/gp_diagnostic_scatter.png\n")

# ==============================================================================
# 7. CCF PLOT
# ==============================================================================

cat("Generating CCF plot...\n")

# Approximate 95% significance threshold
sig_threshold <- 2 / sqrt(N_model)

p_ccf <- ggplot(ccf_df, aes(x = lag, y = ccf, fill = variable)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(-sig_threshold, sig_threshold),
             linetype = "dashed", color = "red", alpha = 0.6) +
  scale_x_continuous(breaks = 0:max_lag) +
  scale_fill_manual(values = c("Temperature" = "coral", "Rainfall" = "steelblue")) +
  labs(
    title = "Cross-Correlation: Climate Variables Leading GP Residual",
    subtitle = sprintf("Dashed red = 95%% significance threshold (%.3f)", sig_threshold),
    x = "Lag (weeks, climate leading)",
    y = "Cross-correlation",
    fill = "Climate Variable"
  )

ggsave("../results/figures/gp_diagnostic_ccf.png", p_ccf,
       width = 10, height = 5, dpi = 150)
cat("  Saved: results/figures/gp_diagnostic_ccf.png\n")

# ==============================================================================
# 8. INTERPRETATION
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("INTERPRETATION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

max_abs_rho <- max(abs(cor_temp_spearman), abs(cor_rain_spearman))

if (max_abs_rho < 0.1) {
  cat("CONCLUSION: GP residual is NOT absorbing climate signal.\n")
  cat(sprintf("  Max |Spearman rho| = %.4f (< 0.1 threshold)\n", max_abs_rho))
  cat("  The GP captures dynamics ORTHOGONAL to climate covariates.\n")
  conclusion <- "GP residual is NOT absorbing climate signal"
} else if (max_abs_rho < 0.3) {
  cat("CONCLUSION: Weak correlation detected between GP and climate.\n")
  cat(sprintf("  Max |Spearman rho| = %.4f\n", max_abs_rho))
  cat("  Some climate signal may leak into the GP, but the effect is modest.\n")
  conclusion <- "Weak correlation: modest climate leakage into GP"
} else {
  cat("WARNING: GP residual shows substantial correlation with climate.\n")
  cat(sprintf("  Max |Spearman rho| = %.4f (>= 0.3)\n", max_abs_rho))
  cat("  The GP may be absorbing nonlinear climate effects.\n")
  conclusion <- "Substantial correlation: GP may absorb climate signal"
}

# ==============================================================================
# 9. SAVE SUMMARY CSV
# ==============================================================================

cat("\nSaving summary...\n")

summary_df <- tibble(
  metric = c(
    "cor_temp_pearson", "cor_temp_spearman",
    "cor_rain_pearson", "cor_rain_spearman",
    "cor_temp_spearman_median_posterior", "cor_temp_spearman_lower95",
    "cor_temp_spearman_upper95",
    "cor_rain_spearman_median_posterior", "cor_rain_spearman_lower95",
    "cor_rain_spearman_upper95",
    "ccf_temp_lag0", "ccf_rain_lag0",
    "max_abs_spearman_rho", "conclusion"
  ),
  value = c(
    sprintf("%.4f", cor_temp_pearson),
    sprintf("%.4f", cor_temp_spearman),
    sprintf("%.4f", cor_rain_pearson),
    sprintf("%.4f", cor_rain_spearman),
    sprintf("%.4f", median(cor_temp_spearman_draws)),
    sprintf("%.4f", quantile(cor_temp_spearman_draws, 0.025)),
    sprintf("%.4f", quantile(cor_temp_spearman_draws, 0.975)),
    sprintf("%.4f", median(cor_rain_spearman_draws)),
    sprintf("%.4f", quantile(cor_rain_spearman_draws, 0.025)),
    sprintf("%.4f", quantile(cor_rain_spearman_draws, 0.975)),
    sprintf("%.4f", ccf_df$ccf[ccf_df$variable == "Temperature" & ccf_df$lag == 0]),
    sprintf("%.4f", ccf_df$ccf[ccf_df$variable == "Rainfall" & ccf_df$lag == 0]),
    sprintf("%.4f", max_abs_rho),
    conclusion
  )
)

write_csv(summary_df, "../results/gp_climate_diagnostic.csv")
cat("  Saved: results/gp_climate_diagnostic.csv\n")

# ==============================================================================
# 10. DONE
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("GP CLIMATE DIAGNOSTIC COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Output files:\n")
cat("  results/figures/gp_diagnostic_scatter.png  - Scatter: f_residual vs climate\n")
cat("  results/figures/gp_diagnostic_ccf.png      - Cross-correlation at lags 0-12\n")
cat("  results/gp_climate_diagnostic.csv          - Summary statistics\n")
