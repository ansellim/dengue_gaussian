#!/usr/bin/env Rscript
# ==============================================================================
# 15_fit_kernel_comparison.R
#
# Cross-kernel comparison study: fit Model 3 with four different GP kernels
# using a single Stan file (15_model3_multikernel.stan) and the kernel_type
# data variable.
#
# Kernels:
#   1 = Matern 1/2 (Exponential) — rough, non-differentiable
#   2 = Matern 3/2 — once-differentiable (original default)
#   3 = Matern 5/2 — twice-differentiable
#   4 = Squared Exponential (RBF) — infinitely smooth
#
# Input:  data/model_data.rds, code/15_model3_multikernel.stan
# Output: results/fit_kernel_matern12.rds
#         results/fit_kernel_matern32.rds
#         results/fit_kernel_matern52.rds
#         results/fit_kernel_sqexp.rds
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)

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

# ==============================================================================
# Configuration
# ==============================================================================

kernel_config <- tibble(
  kernel_type = 1:4,
  name        = c("matern12", "matern32", "matern52", "sqexp"),
  label       = c("Matérn 1/2", "Matérn 3/2", "Matérn 5/2", "Squared Exp."),
  description = c(
    "Rough (continuous, non-differentiable) — Ornstein-Uhlenbeck process",
    "Once-differentiable — EpiNow2 default, balances smoothness and flexibility",
    "Twice-differentiable — smoother than 3/2, less prone to overfitting noise",
    "Infinitely smooth — risk of oversmoothing abrupt Rt changes"
  )
)

# MCMC settings — identical across all kernels for fair comparison
mcmc_settings <- list(
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  adapt_delta     = 0.95,
  max_treedepth   = 12,
  seed            = 42,
  refresh         = 200
)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("CROSS-KERNEL COMPARISON: FITTING 4 GP KERNELS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data  <- model_data$stan_data

cat(sprintf("  N_model = %d weeks\n", stan_data$N_model))
cat(sprintf("  M = %d basis functions\n", stan_data$M))
cat(sprintf("  L = %.1f (boundary factor)\n", stan_data$L))
cat(sprintf("  Cases range: %d - %d\n\n", min(stan_data$cases), max(stan_data$cases)))

# ==============================================================================
# 2. COMPILE MODEL (once)
# ==============================================================================

cat("Compiling multi-kernel Stan model...\n")
model <- cmdstan_model("15_model3_multikernel.stan")
cat("  Compilation successful.\n\n")

# ==============================================================================
# 3. FIT ALL KERNELS
# ==============================================================================

fits <- list()

for (i in seq_len(nrow(kernel_config))) {
  ktype <- kernel_config$kernel_type[i]
  kname <- kernel_config$name[i]
  klabel <- kernel_config$label[i]
  kdesc  <- kernel_config$description[i]

  cat("=" |> rep(70) |> paste(collapse = ""), "\n")
  cat(sprintf("KERNEL %d/4: %s\n", i, klabel))
  cat(sprintf("  %s\n", kdesc))
  cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

  # Add kernel_type to stan_data
  stan_data_k <- c(stan_data, list(kernel_type = ktype))

  # Check if fit already exists
  fit_file <- sprintf("../results/fit_kernel_%s.rds", kname)
  if (file.exists(fit_file)) {
    cat(sprintf("  Fit already exists: %s — skipping.\n", fit_file))
    cat("  (Delete the file to refit.)\n\n")
    fits[[kname]] <- readRDS(fit_file)
    next
  }

  # Fit
  t_start <- Sys.time()

  fit <- model$sample(
    data            = stan_data_k,
    chains          = mcmc_settings$chains,
    parallel_chains = mcmc_settings$parallel_chains,
    iter_warmup     = mcmc_settings$iter_warmup,
    iter_sampling   = mcmc_settings$iter_sampling,
    adapt_delta     = mcmc_settings$adapt_delta,
    max_treedepth   = mcmc_settings$max_treedepth,
    seed            = mcmc_settings$seed,
    refresh         = mcmc_settings$refresh
  )

  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")

  # Save
  fit$save_object(fit_file)
  fits[[kname]] <- fit
  cat(sprintf("\n  Saved: %s (%.1f min)\n\n", fit_file, as.numeric(t_elapsed)))

  # Quick diagnostics
  diag <- fit$diagnostic_summary()
  cat(sprintf("  Divergences: %d\n", sum(diag$num_divergent)))
  cat(sprintf("  Max treedepth exceeded: %d\n\n", sum(diag$num_max_treedepth)))
}

# ==============================================================================
# 4. QUICK SUMMARY TABLE
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("QUICK COMPARISON SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

params_compare <- c("alpha", "rho", "phi", "prop_climate", "prop_residual")

for (i in seq_len(nrow(kernel_config))) {
  kname  <- kernel_config$name[i]
  klabel <- kernel_config$label[i]

  cat(sprintf("--- %s ---\n", klabel))

  if (!is.null(fits[[kname]])) {
    summ <- fits[[kname]]$summary(variables = params_compare)
    print(summ, n = Inf)
  } else {
    cat("  (not loaded)\n")
  }
  cat("\n")
}

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("FITTING COMPLETE\n")
cat("Next step: Run 16_kernel_comparison_postprocess.R\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
