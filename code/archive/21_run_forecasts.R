#!/usr/bin/env Rscript
# ==============================================================================
# 21_run_forecasts.R
#
# Rolling-origin out-of-sample forecast evaluation.
# Fits the forecast Stan model at multiple training cutoffs, with and without
# the Shannon entropy covariate.
#
# Configuration:
#   - Kernel: Matern 3/2 (kernel_type = 2)
#   - Entropy: with (use_entropy=1) and without (use_entropy=0)
#   - Holdout: K_HOLDOUT = 26 weeks at the end of the series
#   - Rolling step: 4 weeks between origins
#   - Forecast horizons evaluated post-hoc: 4, 8, 13 weeks
#
# Input:
#   data/model_data.rds
#   data/forecast_entropy.rds
#   code/20_model3_forecast.stan
#
# Output:
#   results/forecasts/forecast_entropy{0,1}_origin{N}.rds
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

dir.create("../results/forecasts", showWarnings = FALSE, recursive = TRUE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ROLLING-ORIGIN OUT-OF-SAMPLE FORECASTING\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

KERNEL_TYPE   <- 2       # Matern 3/2
K_HOLDOUT     <- 26      # Weeks held out at end of series
ROLLING_STEP  <- 4       # Weeks between forecast origins
HORIZONS      <- c(4, 8, 13)  # Evaluated post-hoc

# MCMC settings
N_CHAINS      <- 4
N_WARMUP      <- 500
N_SAMPLING    <- 500
ADAPT_DELTA   <- 0.95
MAX_TREEDEPTH <- 12
SEED_BASE     <- 2026

cat(sprintf("  Kernel: Matern 3/2 (type %d)\n", KERNEL_TYPE))
cat(sprintf("  Holdout: %d weeks\n", K_HOLDOUT))
cat(sprintf("  Rolling step: %d weeks\n", ROLLING_STEP))
cat(sprintf("  Horizons: %s\n", paste(HORIZONS, collapse = ", ")))
cat(sprintf("  MCMC: %d chains x %d samples (warmup %d)\n",
            N_CHAINS, N_SAMPLING, N_WARMUP))

# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

cat("\nLoading data...\n")

model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data

# Load entropy covariate
entropy_data <- readRDS("../data/forecast_entropy.rds")
entropy_covariate <- entropy_data$entropy_covariate

N_model <- stan_data$N_model
dates_model <- model_data$metadata$dates_model

cat(sprintf("  N_model = %d weeks\n", N_model))
cat(sprintf("  Dates: %s to %s\n", min(dates_model), max(dates_model)))

# ==============================================================================
# 3. DEFINE FORECAST ORIGINS
# ==============================================================================

# The last K_HOLDOUT weeks are reserved for evaluation
# Earliest origin: N_model - K_HOLDOUT (forecast starts at N_model - K_HOLDOUT + 1)
# Latest origin: N_model - max(HORIZONS) (need at least max horizon for eval)
# We step backwards by ROLLING_STEP

max_horizon <- max(HORIZONS)
latest_origin <- N_model - max_horizon
earliest_origin <- N_model - K_HOLDOUT

origins <- seq(from = latest_origin, to = earliest_origin, by = -ROLLING_STEP)
origins <- sort(origins)

cat(sprintf("\n  Forecast origins: %d\n", length(origins)))
for (o in origins) {
  cat(sprintf("    N_train=%d -> forecast weeks %d-%d (%s to %s)\n",
              o, o + 1, min(o + max_horizon, N_model),
              dates_model[o], dates_model[min(o + max_horizon, N_model)]))
}

# ==============================================================================
# 4. COMPILE MODEL
# ==============================================================================

cat("\nCompiling forecast model...\n")
forecast_model <- cmdstan_model("20_model3_forecast.stan")
cat("  Compilation successful.\n")

# ==============================================================================
# 5. ROLLING FORECAST LOOP
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("RUNNING FORECASTS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

entropy_options <- c(0, 1)
total_fits <- length(entropy_options) * length(origins)
fit_count <- 0

for (use_ent in entropy_options) {
  for (n_train in origins) {
    fit_count <- fit_count + 1
    n_forecast <- N_model - n_train

    out_file <- sprintf("../results/forecasts/forecast_entropy%d_origin%d.rds",
                        use_ent, n_train)

    # Skip if already exists (checkpointing)
    if (file.exists(out_file)) {
      cat(sprintf("[%d/%d] SKIP entropy=%d, N_train=%d (file exists)\n",
                  fit_count, total_fits, use_ent, n_train))
      next
    }

    cat(sprintf("[%d/%d] Fitting entropy=%d, N_train=%d, N_forecast=%d...\n",
                fit_count, total_fits, use_ent, n_train, n_forecast))

    # Prepare Stan data for this origin
    fit_data <- stan_data
    fit_data$kernel_type <- KERNEL_TYPE
    fit_data$N_train <- n_train
    fit_data$N_forecast <- n_forecast
    fit_data$use_entropy <- use_ent
    fit_data$entropy_covariate <- entropy_covariate

    # Fit model
    fit <- tryCatch({
      forecast_model$sample(
        data = fit_data,
        chains = N_CHAINS,
        parallel_chains = N_CHAINS,
        iter_warmup = N_WARMUP,
        iter_sampling = N_SAMPLING,
        adapt_delta = ADAPT_DELTA,
        max_treedepth = MAX_TREEDEPTH,
        seed = SEED_BASE + n_train,
        refresh = 0
      )
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      NULL
    })

    if (is.null(fit)) {
      cat("  Skipping due to error.\n\n")
      next
    }

    # Extract diagnostics
    diag <- fit$diagnostic_summary(quiet = TRUE)
    n_div <- sum(diag$num_divergent)
    n_tree <- sum(diag$num_max_treedepth)

    cat(sprintf("  Divergences: %d, Max treedepth: %d\n", n_div, n_tree))

    # Extract draws for forecast period
    cases_forward_draws <- fit$draws("cases_forward", format = "matrix")
    rt_draws <- fit$draws("Rt", format = "matrix")
    log_lik_draws <- fit$draws("log_lik", format = "matrix")

    # Also extract key parameters
    param_draws <- fit$draws(
      variables = c("mu", "alpha", "rho", "phi", "beta_entropy",
                    "beta_climate[1]", "beta_climate[2]"),
      format = "draws_matrix"
    )

    # Observed cases for the forecast window
    S <- stan_data$S
    obs_cases_forecast <- stan_data$cases[(S + n_train + 1):(S + N_model)]

    # Save compact result
    result <- list(
      # Forecast draws (full N_model, caller extracts forecast columns)
      cases_forward = cases_forward_draws,
      Rt = rt_draws,
      log_lik = log_lik_draws,

      # Key parameters
      params = param_draws,

      # Metadata
      N_train = n_train,
      N_forecast = n_forecast,
      N_model = N_model,
      use_entropy = use_ent,
      kernel_type = KERNEL_TYPE,
      dates_model = dates_model,
      observed_cases = stan_data$cases,
      gi = stan_data$gi,
      S = S,

      # Diagnostics
      n_divergences = n_div,
      n_max_treedepth = n_tree,

      # Observed forecast-period cases for easy evaluation
      obs_cases_forecast = obs_cases_forecast
    )

    saveRDS(result, out_file)
    cat(sprintf("  Saved: %s\n\n", basename(out_file)))
  }
}

# ==============================================================================
# 6. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("FORECAST RUNS COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# List saved files
saved_files <- list.files("../results/forecasts", pattern = "^forecast_.*\\.rds$",
                          full.names = FALSE)
cat(sprintf("  Saved %d forecast files in results/forecasts/\n", length(saved_files)))
for (f in saved_files) cat(sprintf("    %s\n", f))

cat("\nNext step: Run 22_evaluate_forecasts.R\n")
