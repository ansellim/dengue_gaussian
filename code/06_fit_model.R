#!/usr/bin/env Rscript
# ==============================================================================
# 06_fit_model.R
#
# Fit the climate-only model: covariates (temperature + rainfall) with short GP
# log(Rt) = mu + beta_temp*temp + beta_rain*rain + f_residual(t)
# Focuses on climate drivers and residual dynamics for serotype switching analysis.
#
# Input:  data/model_data.rds, code/05_model_climate_only.stan
# Output: results/fit_model.rds
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

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
cat("FITTING CLIMATE-ONLY MODEL WITH SHORT GP\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data

# Baseline GP hyperprior means (sigma fixed at 0.5 inside the Stan model).
# The Stan model takes these as data so a single compiled binary serves both
# this baseline fit and the prior sensitivity sweep in 15_gp_prior_sensitivity.R.
stan_data$log_alpha_mu <- -1.2     # alpha median ~ 0.30
stan_data$log_rho_mu   <- log(6)   # rho median = 6 weeks

cat(sprintf("  N_model = %d weeks\n", stan_data$N_model))
cat(sprintf("  Cases range: %d - %d\n", min(stan_data$cases), max(stan_data$cases)))
cat(sprintf("  Covariates: temperature + rainfall (K_climate = %d)\n", stan_data$K_climate))

# ==============================================================================
# 2. COMPILE AND FIT MODEL
# ==============================================================================

cat("\nCompiling climate-only model...\n")
cat("  log(Rt) = mu + beta_temp * temp + beta_rain * rain + f_residual(t)\n")
cat("  GP length scale prior: rho ~ LogNormal(log(6), 0.5)\n\n")

model <- cmdstan_model("05_model_climate_only.stan")

cat("Fitting model...\n")
fit <- model$sample(
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

# Save fit
fit$save_object("../results/fit_model.rds")
cat("\nSaved: results/fit_model.rds\n")

# ==============================================================================
# 3. DIAGNOSTICS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DIAGNOSTICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

diag <- fit$diagnostic_summary()
cat(sprintf("Divergent transitions: %d\n", sum(diag$num_divergent)))
cat(sprintf("Max treedepth exceeded: %d\n", sum(diag$num_max_treedepth)))

# ==============================================================================
# 4. PARAMETER SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PARAMETER SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

params <- c("mu", "alpha", "rho", "phi",
            "temp_effect", "rain_effect",
            "prop_climate", "prop_residual")
summ <- fit$summary(variables = params)
print(summ)

# ==============================================================================
# 5. POSTERIOR PREDICTIVE CHECK - RESIDUAL AUTOCORRELATION
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("RESIDUAL AUTOCORRELATION CHECK\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Extract posterior predictive samples
y_rep <- fit$draws("cases_pred", format = "matrix")
S <- stan_data$S
N <- stan_data$N
y <- stan_data$cases[(S + 1):N]

# Align dimensions
n_common <- min(ncol(y_rep), length(y))
y_rep <- y_rep[, 1:n_common]
y <- y[1:n_common]

# Compute residuals
y_rep_median <- apply(y_rep, 2, median)
y_rep_sd <- apply(y_rep, 2, sd)
std_resid <- (y - y_rep_median) / y_rep_sd

# Autocorrelation
acf_resid <- acf(std_resid, lag.max = 5, plot = FALSE)
cat("Residual autocorrelation:\n")
for (i in 1:5) {
  cat(sprintf("  Lag %d: %.3f\n", i, acf_resid$acf[i + 1]))
}

# Coverage
q05 <- apply(y_rep, 2, quantile, 0.05)
q95 <- apply(y_rep, 2, quantile, 0.95)
coverage_90 <- mean(y >= q05 & y <= q95)
cat(sprintf("\n90%% interval coverage: %.1f%% (nominal: 90%%)\n", 100 * coverage_90))

# ==============================================================================
# 6. VARIANCE DECOMPOSITION
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("VARIANCE DECOMPOSITION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

prop_climate <- fit$summary("prop_climate")
prop_residual <- fit$summary("prop_residual")

cat(sprintf("  Climate:  %.1f%% [%.1f%%, %.1f%%]\n",
            prop_climate$median * 100, prop_climate$q5 * 100, prop_climate$q95 * 100))
cat(sprintf("  Residual: %.1f%% [%.1f%%, %.1f%%]\n",
            prop_residual$median * 100, prop_residual$q5 * 100, prop_residual$q95 * 100))

# ==============================================================================
# 7. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("The model uses climate covariates (temperature + rainfall) and a residual\n")
cat("GP. The residual GP captures all non-climate variation, including\n")
cat("serotype dynamics and immunity shifts.\n\n")

cat("Next step: Run 07_posterior_predictive.R, then 08_postprocess.R\n")
cat("           Run 10_serotype_analysis.R for early warning analysis\n\n")

cat("Files saved:\n")
cat("  results/fit_model.rds\n")
