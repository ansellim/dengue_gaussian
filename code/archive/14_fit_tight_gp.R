#!/usr/bin/env Rscript
# ==============================================================================
# 14_fit_tight_gp.R
#
# Sensitivity analysis: fit Model 3 with tighter GP amplitude prior
# Compare variance decomposition with standard Model 3
#
# Key change: log_alpha ~ Normal(-1.6, 0.3) vs Normal(-1.2, 0.5)
#   Median alpha: 0.20 (vs 0.30) — constrains GP to absorb less variance
#
# Input:  data/model_data.rds, code/05_model3_climate_only_tight_gp.stan
# Output: results/fit_model3_tight_gp.rds
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)

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
cat("SENSITIVITY ANALYSIS: TIGHTER GP AMPLITUDE PRIOR\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data

# ==============================================================================
# 2. COMPILE AND FIT TIGHT GP MODEL
# ==============================================================================

cat("\nCompiling tight GP model...\n")
cat("  log_alpha ~ Normal(-1.6, 0.3): median alpha = 0.20\n")
cat("  vs standard Normal(-1.2, 0.5): median alpha = 0.30\n\n")

model_tight <- cmdstan_model("05_model3_climate_only_tight_gp.stan")

cat("Fitting model...\n")
fit_tight <- model_tight$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  seed = 42,
  refresh = 200
)

fit_tight$save_object("../results/fit_model3_tight_gp.rds")
cat("\nSaved: results/fit_model3_tight_gp.rds\n")

# ==============================================================================
# 3. DIAGNOSTICS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DIAGNOSTICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

diag <- fit_tight$diagnostic_summary()
cat(sprintf("Divergent transitions: %d\n", sum(diag$num_divergent)))
cat(sprintf("Max treedepth exceeded: %d\n", sum(diag$num_max_treedepth)))

# ==============================================================================
# 4. COMPARE WITH STANDARD MODEL 3
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMPARISON: STANDARD vs TIGHT GP\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

if (file.exists("../results/fit_model3.rds")) {
  fit_std <- readRDS("../results/fit_model3.rds")

  # --- Parameter comparison ---
  params <- c("mu", "alpha", "rho", "phi", "temp_effect", "rain_effect")
  summ_std <- fit_std$summary(variables = params)
  summ_tight <- fit_tight$summary(variables = params)

  cat("Parameter estimates:\n")
  cat(sprintf("  %-18s  %12s  %12s\n", "Parameter", "Standard", "Tight GP"))
  cat(sprintf("  %-18s  %12s  %12s\n", "---------", "--------", "--------"))
  for (i in seq_len(nrow(summ_std))) {
    cat(sprintf("  %-18s  %5.3f [%5.3f, %5.3f]  %5.3f [%5.3f, %5.3f]\n",
                summ_std$variable[i],
                summ_std$median[i], summ_std$q5[i], summ_std$q95[i],
                summ_tight$median[i], summ_tight$q5[i], summ_tight$q95[i]))
  }

  # --- Variance decomposition comparison ---
  cat("\nVariance decomposition:\n")

  pc_std <- fit_std$summary("prop_climate")
  pr_std <- fit_std$summary("prop_residual")
  pc_tight <- fit_tight$summary("prop_climate")
  pr_tight <- fit_tight$summary("prop_residual")

  cat(sprintf("  Climate:  Standard = %.1f%%, Tight GP = %.1f%%\n",
              pc_std$median * 100, pc_tight$median * 100))
  cat(sprintf("  Residual: Standard = %.1f%%, Tight GP = %.1f%%\n",
              pr_std$median * 100, pr_tight$median * 100))

  # --- Residual autocorrelation comparison ---
  S <- stan_data$S
  N <- stan_data$N
  y <- stan_data$cases[(S + 1):N]

  compute_acf <- function(fit, y) {
    y_rep <- fit$draws("cases_pred", format = "matrix")
    n <- min(ncol(y_rep), length(y))
    y_rep <- y_rep[, 1:n]
    y_sub <- y[1:n]
    y_rep_med <- apply(y_rep, 2, median)
    y_rep_sd <- apply(y_rep, 2, sd)
    resid <- (y_sub - y_rep_med) / y_rep_sd
    acf(resid, lag.max = 5, plot = FALSE)
  }

  acf_std <- compute_acf(fit_std, y)
  acf_tight <- compute_acf(fit_tight, y)

  cat("\nResidual autocorrelation (lag 1):\n")
  cat(sprintf("  Standard: %.3f\n", acf_std$acf[2]))
  cat(sprintf("  Tight GP: %.3f\n", acf_tight$acf[2]))

  # --- LOO-CV comparison ---
  cat("\nLOO-CV comparison:\n")
  lik_std <- fit_std$draws("log_lik", format = "matrix")
  lik_tight <- fit_tight$draws("log_lik", format = "matrix")

  loo_std <- loo(lik_std)
  loo_tight <- loo(lik_tight)

  comparison <- loo_compare(loo_std, loo_tight)
  print(comparison)

  # --- Coverage comparison ---
  compute_coverage <- function(fit, y) {
    y_rep <- fit$draws("cases_pred", format = "matrix")
    n <- min(ncol(y_rep), length(y))
    y_rep <- y_rep[, 1:n]
    y_sub <- y[1:n]
    q05 <- apply(y_rep, 2, quantile, 0.05)
    q95 <- apply(y_rep, 2, quantile, 0.95)
    mean(y_sub >= q05 & y_sub <= q95)
  }

  cov_std <- compute_coverage(fit_std, y)
  cov_tight <- compute_coverage(fit_tight, y)

  cat(sprintf("\n90%% interval coverage:\n"))
  cat(sprintf("  Standard: %.1f%%\n", 100 * cov_std))
  cat(sprintf("  Tight GP: %.1f%%\n", 100 * cov_tight))

} else {
  cat("Standard Model 3 not found. Run 13_fit_model3.R first.\n")
}

# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Constraining the GP amplitude prior reduces how much variance the\n")
cat("residual GP can absorb. If the variance decomposition shifts toward\n")
cat("climate, it suggests the standard GP was absorbing some climate signal.\n")
cat("If it remains similar, climate effects are genuinely small.\n\n")

cat("Files saved:\n")
cat("  results/fit_model3_tight_gp.rds\n")
