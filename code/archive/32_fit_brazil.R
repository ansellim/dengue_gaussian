#!/usr/bin/env Rscript
# ==============================================================================
# 32_fit_brazil.R
#
# Fit the multi-kernel GP model (Matérn 3/2) to Brazil dengue data.
# Uses the same Stan model and MCMC settings as the Singapore analysis
# for direct comparison.
#
# Input:  data/brazil_model_data.rds, code/15_model3_multikernel.stan
# Output: results/fit_brazil_matern32.rds
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
cat("FITTING BRAZIL DENGUE MODEL (MATÉRN 3/2 KERNEL)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading Brazil data...\n")
model_data <- readRDS("../data/brazil_model_data.rds")
stan_data <- model_data$stan_data

cat(sprintf("  City: %s\n", model_data$metadata$city))
cat(sprintf("  N_model = %d weeks\n", stan_data$N_model))
cat(sprintf("  M = %d basis functions\n", stan_data$M))
cat(sprintf("  L = %.1f (boundary)\n", stan_data$L))
cat(sprintf("  Cases range: %d - %d\n", min(stan_data$cases), max(stan_data$cases)))
cat(sprintf("  Covariates: temperature + rainfall/humidity (K_climate = %d)\n\n",
            stan_data$K_climate))

# ==============================================================================
# 2. COMPILE AND FIT MODEL
# ==============================================================================

cat("Compiling multi-kernel Stan model...\n")
model <- cmdstan_model("15_model3_multikernel.stan")
cat("  Compilation successful.\n\n")

# Add kernel_type = 2 (Matérn 3/2) — same as Singapore default
stan_data_fit <- c(stan_data, list(kernel_type = 2))

# Check if fit already exists
fit_file <- "../results/fit_brazil_matern32.rds"
if (file.exists(fit_file)) {
  cat(sprintf("Fit already exists: %s\n", fit_file))
  cat("Loading existing fit...\n")
  fit <- readRDS(fit_file)
} else {
  cat("Fitting model (Matérn 3/2 kernel)...\n")
  cat("  log(Rt) = mu + beta_temp * temp + beta_rain * rain + f_residual(t)\n")
  cat("  Same MCMC settings as Singapore: 4 chains, 1000/1000, adapt_delta=0.95\n\n")

  t_start <- Sys.time()

  fit <- model$sample(
    data = stan_data_fit,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95,
    max_treedepth = 12,
    seed = 42,
    refresh = 200
  )

  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
  cat(sprintf("\n  Fitting completed in %.1f minutes\n", as.numeric(t_elapsed)))

  # Save fit
  fit$save_object(fit_file)
  cat(sprintf("  Saved: %s\n", fit_file))
}

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
cat("PARAMETER SUMMARY (BRAZIL)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

params <- c("mu", "alpha", "rho", "phi",
            "temp_effect", "rain_effect",
            "prop_climate", "prop_residual")

summ_brazil <- fit$summary(variables = params)
print(summ_brazil)

# ==============================================================================
# 5. POSTERIOR PREDICTIVE CHECK
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("POSTERIOR PREDICTIVE CHECK\n")
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
# 7. COMPARISON WITH SINGAPORE
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMPARISON WITH SINGAPORE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Try to load Singapore results for comparison
sg_fit_file <- "../results/fit_kernel_matern32.rds"
if (file.exists(sg_fit_file)) {
  cat("Loading Singapore Matérn 3/2 fit for comparison...\n\n")
  fit_sg <- readRDS(sg_fit_file)

  summ_sg <- fit_sg$summary(variables = params)

  # Side-by-side comparison
  cat(sprintf("%-20s %12s %12s\n", "Parameter", "Singapore", "Brazil"))
  cat(sprintf("%-20s %12s %12s\n", "---------", "---------", "------"))

  for (p in params) {
    sg_row <- summ_sg |> filter(variable == p)
    br_row <- summ_brazil |> filter(variable == p)

    if (nrow(sg_row) > 0 && nrow(br_row) > 0) {
      cat(sprintf("%-20s %5.3f [%5.3f] %5.3f [%5.3f]\n",
                  p,
                  sg_row$median, sg_row$mad,
                  br_row$median, br_row$mad))
    }
  }

  cat("\n  Key differences to examine:\n")
  cat("  - alpha (GP amplitude): higher = more non-climate variation\n")
  cat("  - rho (GP length scale): shorter = faster Rt fluctuations\n")
  cat("  - prop_climate: fraction of log(Rt) variance explained by climate\n")
  cat("  - temp_effect / rain_effect: climate sensitivity differences\n")

} else {
  cat("Singapore fit not found at: results/fit_kernel_matern32.rds\n")
  cat("Run 15_fit_kernel_comparison.R first, then re-run this script\n")
  cat("to see the cross-country comparison.\n")
}

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Fitted the Matérn 3/2 GP model to São Paulo dengue data (2012-2022).\n")
cat("The model structure is identical to Singapore:\n")
cat("  log(Rt) = mu + beta_temp * temp + beta_rain * rain + f_GP(t)\n\n")
cat("Brazil lacks Wolbachia and NPI covariates (as in Singapore Model 3),\n")
cat("so the GP captures all non-climate Rt variation.\n\n")

cat("Files saved:\n")
cat("  results/fit_brazil_matern32.rds\n")
cat("\nNext steps:\n")
cat("  - Compare Rt trajectories between Singapore and Brazil\n")
cat("  - Examine whether climate explains more/less variance in tropical Brazil\n")
cat("  - Check if GP length scale differs (serotype dynamics vs. immunity patterns)\n")
