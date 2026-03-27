#!/usr/bin/env Rscript
# ==============================================================================
# 13_fit_model3.R
#
# Fit Model 3: Climate-only covariates (temperature + rainfall) with short GP
# Removes Wolbachia and NPI covariates to focus on climate drivers and
# residual dynamics for serotype switching analysis.
#
# Input:  data/model_data.rds, code/05_model3_climate_only.stan
# Output: results/fit_model3.rds
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
cat("FITTING MODEL 3: CLIMATE-ONLY WITH SHORT GP\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data

cat(sprintf("  N_model = %d weeks\n", stan_data$N_model))
cat(sprintf("  Cases range: %d - %d\n", min(stan_data$cases), max(stan_data$cases)))
cat(sprintf("  Covariates: temperature + rainfall (K_climate = %d)\n", stan_data$K_climate))

# ==============================================================================
# 2. COMPILE AND FIT MODEL
# ==============================================================================

cat("\nCompiling Model 3 (climate-only)...\n")
cat("  log(Rt) = mu + beta_temp * temp + beta_rain * rain + f_residual(t)\n")
cat("  GP length scale prior: rho ~ LogNormal(log(6), 0.5)\n\n")

model3 <- cmdstan_model("05_model3_climate_only.stan")

cat("Fitting model...\n")
fit3 <- model3$sample(
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
fit3$save_object("../results/fit_model3.rds")
cat("\nSaved: results/fit_model3.rds\n")

# ==============================================================================
# 3. DIAGNOSTICS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DIAGNOSTICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

diag <- fit3$diagnostic_summary()
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
summ <- fit3$summary(variables = params)
print(summ)

# ==============================================================================
# 5. POSTERIOR PREDICTIVE CHECK - RESIDUAL AUTOCORRELATION
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("RESIDUAL AUTOCORRELATION CHECK\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Extract posterior predictive samples
y_rep <- fit3$draws("cases_pred", format = "matrix")
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

prop_climate <- fit3$summary("prop_climate")
prop_residual <- fit3$summary("prop_residual")

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

cat("Model 3 removes Wolbachia and NPI covariates, keeping only climate\n")
cat("(temperature + rainfall) and a residual GP. The residual GP now captures\n")
cat("all non-climate variation, including serotype dynamics, immunity shifts,\n")
cat("and intervention effects.\n\n")

cat("Next step: Run 08_serotype_analysis.R for early warning analysis\n")
cat("           Run 14_fit_tight_gp.R for GP sensitivity analysis\n\n")

cat("Files saved:\n")
cat("  results/fit_model3.rds\n")
