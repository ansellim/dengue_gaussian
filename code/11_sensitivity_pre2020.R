#!/usr/bin/env Rscript
# ==============================================================================
# 11_sensitivity_pre2020.R
#
# Sensitivity Analysis: Fit model on pre-2020 data only (WITHOUT NPI)
#
# Purpose: Estimate climate/Wolbachia effects without NPI-DENV3 confounding
#
# The 2020+ period has two coincident changes:
#   1. COVID-19 NPIs (lockdowns, movement restrictions)
#   2. DENV-3 emergence (new serotype → susceptible population)
#
# By fitting on 2012-2019 data only, we can:
#   - Estimate effects without this confounding
#   - See how much of 2020+ is "out of sample" (unexplained)
#
# NOTE: Uses model WITHOUT NPI covariate (05_model2_no_npi.stan)
#       since NPI=0 for all pre-2020 observations
#
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

dir.create("../results", showWarnings = FALSE)
dir.create("../results/figures", showWarnings = FALSE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SENSITIVITY ANALYSIS: PRE-2020 MODEL\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD AND SUBSET DATA
# ==============================================================================

cat("Loading and subsetting data to pre-2020...\n")

model_data <- readRDS("../data/model_data.rds")
stan_data_full <- model_data$stan_data
metadata <- model_data$metadata

# Get dates for the modeled period
dates_model <- metadata$dates_model
S <- stan_data_full$S

# Find cutoff index (last week before 2020)
cutoff_date <- as.Date("2020-01-01")
pre2020_idx <- which(dates_model < cutoff_date)
n_pre2020 <- length(pre2020_idx)

cat(sprintf("  Full model period: %s to %s (%d weeks)\n",
            min(dates_model), max(dates_model), length(dates_model)))
cat(sprintf("  Pre-2020 period: %s to %s (%d weeks)\n",
            min(dates_model), max(dates_model[pre2020_idx]), n_pre2020))
cat(sprintf("  Excluded 2020+ period: %d weeks\n",
            length(dates_model) - n_pre2020))

# Create pre-2020 Stan data (WITHOUT NPI covariate)
# Need to include S weeks before the modeled period for the renewal equation
# Only use first 3 covariates: temp, rain, wolbachia (exclude NPI column 4)
X_no_npi <- stan_data_full$X_full[1:n_pre2020, 1:3, drop = FALSE]

stan_data_pre2020 <- list(
  N = S + n_pre2020,  # Total weeks including burn-in
  N_model = n_pre2020,
  S = S,
  M = stan_data_full$M,
  cases = stan_data_full$cases[1:(S + n_pre2020)],
  t = stan_data_full$t[1:n_pre2020],
  L = stan_data_full$L,  # Keep same boundary factor
  gi = stan_data_full$gi,
  K = 3,  # Only 3 covariates: temp, rain, wolbachia (no NPI)
  X = X_no_npi
)

cat(sprintf("\n  Pre-2020 Stan data (no NPI):\n"))
cat(sprintf("    N = %d, N_model = %d, S = %d\n",
            stan_data_pre2020$N, stan_data_pre2020$N_model, stan_data_pre2020$S))
cat(sprintf("    K = %d covariates (temp, rain, wolbachia)\n", stan_data_pre2020$K))
cat(sprintf("    Cases range: %d - %d\n",
            min(stan_data_pre2020$cases), max(stan_data_pre2020$cases)))

# ==============================================================================
# 2. FIT PRE-2020 MODEL
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("FITTING PRE-2020 MODEL (No NPI, Short GP)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

model <- cmdstan_model("05_model2_no_npi.stan")

fit_pre2020 <- model$sample(
  data = stan_data_pre2020,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  seed = 42,
  refresh = 200
)

fit_pre2020$save_object("../results/fit_model2_pre2020.rds")
cat("\nSaved: results/fit_model2_pre2020.rds\n")

# ==============================================================================
# 3. DIAGNOSTICS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DIAGNOSTICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

diag <- fit_pre2020$diagnostic_summary()
cat(sprintf("Divergent transitions: %d\n", sum(diag$num_divergent)))
cat(sprintf("Max treedepth exceeded: %d\n", sum(diag$num_max_treedepth)))

# ==============================================================================
# 4. COMPARE PARAMETER ESTIMATES
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PARAMETER COMPARISON: PRE-2020 vs FULL MODEL\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Load full model for comparison
fit_full <- readRDS("../results/fit_model2_short_gp.rds")

# Extract key parameters (no NPI for pre-2020 model)
params_pre2020 <- c("mu", "alpha", "rho", "phi",
                    "temp_effect", "rain_effect", "wolbachia_effect")
params_full <- c("mu", "alpha", "rho", "phi",
                 "temp_effect", "rain_effect", "wolbachia_effect", "npi_effect")

summ_pre2020 <- fit_pre2020$summary(variables = params_pre2020)
summ_full <- fit_full$summary(variables = params_full)

cat("Parameter estimates:\n\n")
cat(sprintf("%-20s %15s %15s %15s\n",
            "Parameter", "Pre-2020", "Full (2012-22)", "Difference"))
cat(sprintf("%s\n", paste(rep("-", 70), collapse = "")))

# Parameters in both models
for (p in params_pre2020) {
  pre <- summ_pre2020 |> filter(variable == p)
  full <- summ_full |> filter(variable == p)

  diff <- pre$median - full$median

  cat(sprintf("%-20s %7.3f [%5.2f-%5.2f] %7.3f [%5.2f-%5.2f] %+7.3f\n",
              p,
              pre$median, pre$q5, pre$q95,
              full$median, full$q5, full$q95,
              diff))
}

# NPI effect (only in full model)
npi_full <- summ_full |> filter(variable == "npi_effect")
cat(sprintf("%-20s %7s %15s %7.3f [%5.2f-%5.2f] %7s\n",
            "npi_effect", "N/A", "",
            npi_full$median, npi_full$q5, npi_full$q95, "N/A"))

# ==============================================================================
# 5. DETAILED COMPARISON OF KEY EFFECTS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DETAILED EFFECT COMPARISON\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Extract draws
wolb_pre <- fit_pre2020$draws("wolbachia_effect", format = "matrix")[,1]
wolb_full <- fit_full$draws("wolbachia_effect", format = "matrix")[,1]

npi_full <- fit_full$draws("npi_effect", format = "matrix")[,1]

cat("WOLBACHIA EFFECT:\n")
cat(sprintf("  Pre-2020 (no NPI): %.3f [%.3f, %.3f]\n",
            median(wolb_pre), quantile(wolb_pre, 0.025), quantile(wolb_pre, 0.975)))
cat(sprintf("  Full model:        %.3f [%.3f, %.3f]\n",
            median(wolb_full), quantile(wolb_full, 0.025), quantile(wolb_full, 0.975)))
cat(sprintf("  Difference: %+.3f\n", median(wolb_pre) - median(wolb_full)))
cat(sprintf("  P(reduces Rt) Pre-2020: %.1f%%\n", 100 * mean(wolb_pre < 1)))
cat(sprintf("  P(reduces Rt) Full:     %.1f%%\n", 100 * mean(wolb_full < 1)))

cat("\nNPI EFFECT (full model only):\n")
cat(sprintf("  Full model: %.3f [%.3f, %.3f]\n",
            median(npi_full), quantile(npi_full, 0.025), quantile(npi_full, 0.975)))
cat(sprintf("  P(increases Rt): %.1f%%\n", 100 * mean(npi_full > 1)))
cat("  NOTE: Pre-2020 model does not include NPI (was 0 for all obs)\n")
cat("        NPI effect in full model is confounded with DENV-3 emergence\n")

# ==============================================================================
# 6. CLIMATE EFFECTS COMPARISON
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("CLIMATE EFFECTS (stable across periods)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

temp_pre <- fit_pre2020$draws("temp_effect", format = "matrix")[,1]
temp_full <- fit_full$draws("temp_effect", format = "matrix")[,1]
rain_pre <- fit_pre2020$draws("rain_effect", format = "matrix")[,1]
rain_full <- fit_full$draws("rain_effect", format = "matrix")[,1]

cat("TEMPERATURE EFFECT:\n")
cat(sprintf("  Pre-2020:  median=%.4f, 95%% CI [%.4f, %.4f]\n",
            median(temp_pre), quantile(temp_pre, 0.025), quantile(temp_pre, 0.975)))
cat(sprintf("  Full:      median=%.4f, 95%% CI [%.4f, %.4f]\n",
            median(temp_full), quantile(temp_full, 0.025), quantile(temp_full, 0.975)))

cat("\nRAINFALL EFFECT:\n")
cat(sprintf("  Pre-2020:  median=%.4f, 95%% CI [%.4f, %.4f]\n",
            median(rain_pre), quantile(rain_pre, 0.025), quantile(rain_pre, 0.975)))
cat(sprintf("  Full:      median=%.4f, 95%% CI [%.4f, %.4f]\n",
            median(rain_full), quantile(rain_full, 0.025), quantile(rain_full, 0.975)))

# ==============================================================================
# 6b. VARIANCE DECOMPOSITION COMPARISON
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("VARIANCE DECOMPOSITION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Pre-2020 model (no NPI)
prop_climate_pre <- fit_pre2020$draws("prop_climate", format = "matrix")[,1]
prop_wolb_pre <- fit_pre2020$draws("prop_wolbachia", format = "matrix")[,1]
prop_resid_pre <- fit_pre2020$draws("prop_residual", format = "matrix")[,1]

# Full model
prop_climate_full <- fit_full$draws("prop_climate", format = "matrix")[,1]
prop_wolb_full <- fit_full$draws("prop_wolbachia", format = "matrix")[,1]
prop_npi_full <- fit_full$draws("prop_npi", format = "matrix")[,1]
prop_resid_full <- fit_full$draws("prop_residual", format = "matrix")[,1]

cat("PRE-2020 MODEL (no NPI):\n")
cat(sprintf("  Climate:      %.1f%% [%.1f%%, %.1f%%]\n",
            100*median(prop_climate_pre), 100*quantile(prop_climate_pre, 0.025),
            100*quantile(prop_climate_pre, 0.975)))
cat(sprintf("  Wolbachia:    %.1f%% [%.1f%%, %.1f%%]\n",
            100*median(prop_wolb_pre), 100*quantile(prop_wolb_pre, 0.025),
            100*quantile(prop_wolb_pre, 0.975)))
cat(sprintf("  Residual GP:  %.1f%% [%.1f%%, %.1f%%]\n",
            100*median(prop_resid_pre), 100*quantile(prop_resid_pre, 0.025),
            100*quantile(prop_resid_pre, 0.975)))

cat("\nFULL MODEL:\n")
cat(sprintf("  Climate:      %.1f%% [%.1f%%, %.1f%%]\n",
            100*median(prop_climate_full), 100*quantile(prop_climate_full, 0.025),
            100*quantile(prop_climate_full, 0.975)))
cat(sprintf("  Wolbachia:    %.1f%% [%.1f%%, %.1f%%]\n",
            100*median(prop_wolb_full), 100*quantile(prop_wolb_full, 0.025),
            100*quantile(prop_wolb_full, 0.975)))
cat(sprintf("  NPI:          %.1f%% [%.1f%%, %.1f%%]\n",
            100*median(prop_npi_full), 100*quantile(prop_npi_full, 0.025),
            100*quantile(prop_npi_full, 0.975)))
cat(sprintf("  Residual GP:  %.1f%% [%.1f%%, %.1f%%]\n",
            100*median(prop_resid_full), 100*quantile(prop_resid_full, 0.025),
            100*quantile(prop_resid_full, 0.975)))

# ==============================================================================
# 7. POSTERIOR PREDICTIVE FOR 2020+ (OUT OF SAMPLE)
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("OUT-OF-SAMPLE ANALYSIS: 2020+ PERIOD\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Get observed cases for 2020+
post2020_idx <- which(dates_model >= cutoff_date)
cases_2020plus <- stan_data_full$cases[(S + post2020_idx)]

cat(sprintf("2020+ observed cases: n=%d weeks\n", length(cases_2020plus)))
cat(sprintf("  Mean: %.0f, SD: %.0f\n", mean(cases_2020plus), sd(cases_2020plus)))
cat(sprintf("  Range: %d - %d\n", min(cases_2020plus), max(cases_2020plus)))

# Compare with pre-2020 period
cases_pre2020 <- stan_data_pre2020$cases[(S + 1):stan_data_pre2020$N]
cat(sprintf("\nPre-2020 observed cases: n=%d weeks\n", length(cases_pre2020)))
cat(sprintf("  Mean: %.0f, SD: %.0f\n", mean(cases_pre2020), sd(cases_pre2020)))
cat(sprintf("  Range: %d - %d\n", min(cases_pre2020), max(cases_pre2020)))

cat(sprintf("\n2020+ vs Pre-2020 mean ratio: %.2f\n",
            mean(cases_2020plus) / mean(cases_pre2020)))

# Extract Rt estimates from both models for the pre-2020 period
Rt_pre2020_model <- fit_pre2020$draws("Rt", format = "matrix")
Rt_full_model <- fit_full$draws("Rt", format = "matrix")[, 1:n_pre2020]

cat("\nRt estimates for pre-2020 period:\n")
cat(sprintf("  Pre-2020 model: median=%.2f (mean across time)\n",
            median(rowMeans(Rt_pre2020_model))))
cat(sprintf("  Full model:     median=%.2f (mean across time)\n",
            median(rowMeans(Rt_full_model))))

# ==============================================================================
# 8. VISUALIZATION
# ==============================================================================

cat("\nGenerating comparison plots...\n")

theme_set(theme_minimal(base_size = 12))

# --- Plot 1: Parameter comparison ---
# Only compare parameters in both models (no NPI in pre-2020)
param_comparison <- tibble(
  parameter = rep(c("Wolbachia", "Temperature", "Rainfall"), each = 2),
  model = rep(c("Pre-2020 (no NPI)", "Full"), 3),
  median = c(median(wolb_pre), median(wolb_full),
             median(temp_pre), median(temp_full),
             median(rain_pre), median(rain_full)),
  lower = c(quantile(wolb_pre, 0.025), quantile(wolb_full, 0.025),
            quantile(temp_pre, 0.025), quantile(temp_full, 0.025),
            quantile(rain_pre, 0.025), quantile(rain_full, 0.025)),
  upper = c(quantile(wolb_pre, 0.975), quantile(wolb_full, 0.975),
            quantile(temp_pre, 0.975), quantile(temp_full, 0.975),
            quantile(rain_pre, 0.975), quantile(rain_full, 0.975))
)

p_params <- ggplot(param_comparison,
                   aes(x = parameter, y = median, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Pre-2020 (no NPI)" = "steelblue", "Full" = "darkred")) +
  labs(
    title = "Effect Estimates: Pre-2020 (no NPI) vs Full Model",
    subtitle = "Dashed line = no effect (multiplicative scale)",
    x = NULL,
    y = "Multiplicative effect on Rt",
    color = "Model"
  ) +
  theme(legend.position = "bottom")

ggsave("../results/figures/sensitivity_pre2020_params.png", p_params,
       width = 10, height = 6, dpi = 150)
cat("  Saved: results/figures/sensitivity_pre2020_params.png\n")

# --- Plot 2: Rt trajectory comparison ---
dates_pre2020 <- dates_model[1:n_pre2020]

Rt_pre2020_summary <- tibble(
  date = dates_pre2020,
  median_pre = apply(Rt_pre2020_model, 2, median),
  lower_pre = apply(Rt_pre2020_model, 2, quantile, 0.1),
  upper_pre = apply(Rt_pre2020_model, 2, quantile, 0.9),
  median_full = apply(Rt_full_model, 2, median),
  lower_full = apply(Rt_full_model, 2, quantile, 0.1),
  upper_full = apply(Rt_full_model, 2, quantile, 0.9)
)

p_Rt <- ggplot(Rt_pre2020_summary, aes(x = date)) +
  geom_ribbon(aes(ymin = lower_pre, ymax = upper_pre),
              fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = median_pre, color = "Pre-2020 (no NPI)"), linewidth = 0.8) +
  geom_line(aes(y = median_full, color = "Full model"),
            linewidth = 0.8, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = c("Pre-2020 (no NPI)" = "steelblue",
                                "Full model" = "darkred")) +
  labs(
    title = "Rt Estimates: Pre-2020 Period",
    subtitle = "Pre-2020 model (no NPI) vs full model",
    x = "Date",
    y = "Rt",
    color = NULL
  ) +
  theme(legend.position = "bottom")

ggsave("../results/figures/sensitivity_pre2020_Rt.png", p_Rt,
       width = 12, height = 6, dpi = 150)
cat("  Saved: results/figures/sensitivity_pre2020_Rt.png\n")

# --- Plot 3: Wolbachia effect density comparison ---
wolb_df <- tibble(
  value = c(wolb_pre, wolb_full),
  model = rep(c("Pre-2020 (no NPI)", "Full"), each = length(wolb_pre))
)

p_wolb <- ggplot(wolb_df, aes(x = value, fill = model)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Pre-2020 (no NPI)" = "steelblue", "Full" = "darkred")) +
  labs(
    title = "Wolbachia Effect Posterior: Pre-2020 (no NPI) vs Full Model",
    subtitle = "Red line = no effect",
    x = "Multiplicative effect on Rt",
    y = "Density",
    fill = "Model"
  ) +
  coord_cartesian(xlim = c(0, 2)) +
  theme(legend.position = "bottom")

ggsave("../results/figures/sensitivity_pre2020_wolbachia.png", p_wolb,
       width = 8, height = 6, dpi = 150)
cat("  Saved: results/figures/sensitivity_pre2020_wolbachia.png\n")

# ==============================================================================
# 9. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. WOLBACHIA EFFECT:\n")
wolb_diff <- median(wolb_pre) - median(wolb_full)
cat(sprintf("   Pre-2020 (no NPI): %.3f vs Full: %.3f (diff: %+.3f)\n",
            median(wolb_pre), median(wolb_full), wolb_diff))
if (abs(wolb_diff) < 0.1) {
  cat("   → STABLE across models (robust to 2020+ confounding)\n")
} else {
  cat("   → CHANGED between models (affected by 2020+ confounding)\n")
}

cat("\n2. NPI EFFECT:\n")
cat("   Pre-2020: NOT INCLUDED (NPI=0 for all pre-2020 observations)\n")
cat(sprintf("   Full model: %.3f (confounded with DENV-3 emergence)\n",
            median(npi_full)))
cat("   → Cannot separate NPI from serotype effect in full model\n")
cat("   → Pre-2020 model correctly excludes unidentifiable parameter\n")

cat("\n3. CLIMATE EFFECTS:\n")
temp_diff <- median(temp_pre) - median(temp_full)
rain_diff <- median(rain_pre) - median(rain_full)
cat(sprintf("   Temperature: Pre-2020=%.4f, Full=%.4f (diff: %+.4f)\n",
            median(temp_pre), median(temp_full), temp_diff))
cat(sprintf("   Rainfall: Pre-2020=%.4f, Full=%.4f (diff: %+.4f)\n",
            median(rain_pre), median(rain_full), rain_diff))

cat("\n4. INTERPRETATION:\n")
cat("   The pre-2020 model excludes NPI (which was 0 for all observations)\n")
cat("   and estimates Wolbachia/climate effects without NPI-DENV3 confounding.\n")
cat("   If Wolbachia effect is similar between models, it suggests this\n")
cat("   effect is robust to the 2020+ confounding period.\n")
cat(sprintf("   The NPI effect in the full model (%.2f) cannot be causally\n",
            median(npi_full)))
cat("   attributed to lockdowns vs DENV-3 emergence - both occurred in 2020.\n")

cat("\nFiles saved:\n")
cat("  results/fit_model2_pre2020.rds\n")
cat("  results/figures/sensitivity_pre2020_params.png\n")
cat("  results/figures/sensitivity_pre2020_Rt.png\n")
cat("  results/figures/sensitivity_pre2020_wolbachia.png\n")
