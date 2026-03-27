#!/usr/bin/env Rscript
# ==============================================================================
# 06_fit_models.R
#
# Fit all Stan models for dengue Rt estimation
# Uses cmdstanr for Stan interface
#
# Input: data/model_data.rds
# Output: results/fit_model*.rds, results/diagnostics.txt
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

# Set up paths
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(commandArgs(trailingOnly = FALSE) |>
    grep("--file=", x = _, value = TRUE) |>
    sub("--file=", "", x = _) |>
    normalizePath() |>
    dirname())
}

# Create results directory
dir.create("../results", showWarnings = FALSE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DENGUE RT ESTIMATION - MODEL FITTING\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data

cat(sprintf("  N_model = %d weeks\n", stan_data$N_model))
cat(sprintf("  M = %d basis functions\n", stan_data$M))

# ==============================================================================
# 2. COMPILE MODELS
# ==============================================================================

cat("\nCompiling Stan models...\n")

# Model 0: Baseline
cat("  Compiling Model 0 (baseline)...\n")
model0 <- cmdstan_model("03_model0_baseline.stan")

# Model 1: Climate
cat("  Compiling Model 1 (climate)...\n")
model1 <- cmdstan_model("04_model1_climate.stan")

# Model 2: Full
cat("  Compiling Model 2 (full)...\n")
model2 <- cmdstan_model("05_model2_full.stan")

cat("  All models compiled successfully.\n")

# ==============================================================================
# 3. MCMC SETTINGS
# ==============================================================================

# MCMC configuration
MCMC_CONFIG <- list(
  chains = 5,
  parallel_chains = 5,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  seed = 12345
)

cat(sprintf("\nMCMC configuration:
  Chains: %d

  Warmup: %d
  Sampling: %d
  Adapt delta: %.2f
  Max treedepth: %d
\n",
  MCMC_CONFIG$chains,
  MCMC_CONFIG$iter_warmup,
  MCMC_CONFIG$iter_sampling,
  MCMC_CONFIG$adapt_delta,
  MCMC_CONFIG$max_treedepth
))

# ==============================================================================
# 4. FIT MODEL 0 (BASELINE)
# ==============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("Fitting Model 0 (Baseline: single GP)...\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

fit0 <- model0$sample(
  data = stan_data,
  chains = MCMC_CONFIG$chains,
  parallel_chains = MCMC_CONFIG$parallel_chains,
  iter_warmup = MCMC_CONFIG$iter_warmup,
  iter_sampling = MCMC_CONFIG$iter_sampling,
  adapt_delta = MCMC_CONFIG$adapt_delta,
  max_treedepth = MCMC_CONFIG$max_treedepth,
  seed = MCMC_CONFIG$seed,
  refresh = 200
)

# Save fit
fit0$save_object("../results/fit_model0.rds")
cat("  Saved: results/fit_model0.rds\n")

# Quick diagnostics
cat("\nModel 0 diagnostics:\n")
print(fit0$diagnostic_summary())
print(fit0$summary(variables = c("mu", "alpha", "rho", "phi")))

# ==============================================================================
# 5. FIT MODEL 1 (CLIMATE)
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("Fitting Model 1 (Climate + residual GP)...\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

fit1 <- model1$sample(
  data = stan_data,
  chains = MCMC_CONFIG$chains,
  parallel_chains = MCMC_CONFIG$parallel_chains,
  iter_warmup = MCMC_CONFIG$iter_warmup,
  iter_sampling = MCMC_CONFIG$iter_sampling,
  adapt_delta = MCMC_CONFIG$adapt_delta,
  max_treedepth = MCMC_CONFIG$max_treedepth,
  seed = MCMC_CONFIG$seed,
  refresh = 200
)

# Save fit
fit1$save_object("../results/fit_model1.rds")
cat("  Saved: results/fit_model1.rds\n")

# Quick diagnostics
cat("\nModel 1 diagnostics:\n")
print(fit1$diagnostic_summary())
print(fit1$summary(variables = c("mu", "beta_climate", "alpha", "rho", "phi",
                                   "prop_climate", "residual_amplitude")))

# ==============================================================================
# 6. FIT MODEL 2 (FULL)
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("Fitting Model 2 (Climate + Wolbachia + NPI + residual GP)...\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

fit2 <- model2$sample(
  data = stan_data,
  chains = MCMC_CONFIG$chains,
  parallel_chains = MCMC_CONFIG$parallel_chains,
  iter_warmup = MCMC_CONFIG$iter_warmup,
  iter_sampling = MCMC_CONFIG$iter_sampling,
  adapt_delta = MCMC_CONFIG$adapt_delta,
  max_treedepth = MCMC_CONFIG$max_treedepth,
  seed = MCMC_CONFIG$seed,
  refresh = 200
)

# Save fit
fit2$save_object("../results/fit_model2.rds")
cat("  Saved: results/fit_model2.rds\n")

# Quick diagnostics
cat("\nModel 2 diagnostics:\n")
print(fit2$diagnostic_summary())
print(fit2$summary(variables = c("mu", "beta", "alpha", "rho", "phi",
                                   "temp_effect", "rain_effect",
                                   "wolbachia_effect", "npi_effect",
                                   "prop_climate", "prop_wolbachia",
                                   "prop_npi", "prop_residual")))

# ==============================================================================
# 7. COMPREHENSIVE DIAGNOSTICS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("Running comprehensive diagnostics...\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

# Function to check diagnostics
check_diagnostics <- function(fit, model_name) {
  diag <- fit$diagnostic_summary()

  cat(sprintf("\n%s:\n", model_name))
  cat(sprintf("  Divergent transitions: %d\n", sum(diag$num_divergent)))
  cat(sprintf("  Max treedepth exceeded: %d\n", sum(diag$num_max_treedepth)))

  # Get summary for key parameters
  summ <- fit$summary()

  # Check Rhat
  rhat_problems <- summ |>
    filter(rhat > 1.01) |>
    nrow()
  cat(sprintf("  Parameters with Rhat > 1.01: %d\n", rhat_problems))

  # Check ESS
  ess_bulk_min <- min(summ$ess_bulk, na.rm = TRUE)
  ess_tail_min <- min(summ$ess_tail, na.rm = TRUE)
  cat(sprintf("  Min ESS bulk: %.0f\n", ess_bulk_min))
  cat(sprintf("  Min ESS tail: %.0f\n", ess_tail_min))

  # Return diagnostics summary
  list(
    model = model_name,
    divergent = sum(diag$num_divergent),
    max_treedepth = sum(diag$num_max_treedepth),
    rhat_problems = rhat_problems,
    ess_bulk_min = ess_bulk_min,
    ess_tail_min = ess_tail_min,
    converged = sum(diag$num_divergent) == 0 &&
                rhat_problems == 0 &&
                ess_bulk_min > 400 &&
                ess_tail_min > 400
  )
}

diag0 <- check_diagnostics(fit0, "Model 0 (Baseline)")
diag1 <- check_diagnostics(fit1, "Model 1 (Climate)")
diag2 <- check_diagnostics(fit2, "Model 2 (Full)")

# Save diagnostics to file
sink("../results/diagnostics.txt")
cat("DENGUE RT MODEL DIAGNOSTICS\n")
cat("=" |> rep(50) |> paste(collapse = ""), "\n")
cat(sprintf("Generated: %s\n\n", Sys.time()))

for (d in list(diag0, diag1, diag2)) {
  cat(sprintf("%s:\n", d$model))
  cat(sprintf("  Divergent transitions: %d\n", d$divergent))
  cat(sprintf("  Max treedepth exceeded: %d\n", d$max_treedepth))
  cat(sprintf("  Parameters with Rhat > 1.01: %d\n", d$rhat_problems))
  cat(sprintf("  Min ESS bulk: %.0f\n", d$ess_bulk_min))
  cat(sprintf("  Min ESS tail: %.0f\n", d$ess_tail_min))
  cat(sprintf("  CONVERGED: %s\n\n", ifelse(d$converged, "YES", "NO")))
}
sink()

cat("\nDiagnostics saved to: results/diagnostics.txt\n")

# ==============================================================================
# 8. MODEL COMPARISON (LOO-CV approximation)
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("Computing LOO-CV for model comparison...\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

library(loo)

# Extract log-likelihoods
log_lik0 <- fit0$draws("log_lik", format = "matrix")
log_lik1 <- fit1$draws("log_lik", format = "matrix")
log_lik2 <- fit2$draws("log_lik", format = "matrix")

# Compute LOO
loo0 <- loo(log_lik0)
loo1 <- loo(log_lik1)
loo2 <- loo(log_lik2)

# Compare
cat("\nLOO-CV Comparison:\n")
comparison <- loo_compare(loo0, loo1, loo2)
print(comparison)

# Save comparison
saveRDS(list(loo0 = loo0, loo1 = loo1, loo2 = loo2, comparison = comparison),
        "../results/loo_comparison.rds")
cat("\nSaved: results/loo_comparison.rds\n")

# ==============================================================================
# 9. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MODEL FITTING COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\nFiles saved:
  results/fit_model0.rds - Baseline model fit
  results/fit_model1.rds - Climate model fit
  results/fit_model2.rds - Full model fit
  results/diagnostics.txt - MCMC diagnostics
  results/loo_comparison.rds - LOO-CV comparison
\n")

cat("\nNext step: Run 07_postprocess.R for posterior summaries and figures\n")
