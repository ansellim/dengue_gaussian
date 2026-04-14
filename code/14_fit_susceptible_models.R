#!/usr/bin/env Rscript
# ==============================================================================
# 14_fit_susceptible_models.R
#
# Methodological audit pipeline (M0 / M1 / M2): fit M1 (climate + log S_pop + GP)
# and M2 (climate + log S_dom + GP) across an S0 sensitivity sweep, holding all
# other priors and MCMC settings identical to 06_fit_model.R.
#
# Configurations (loaded from data/susceptible_covariates.rds, which is built
# by 14_prepare_susceptible_covariates.R):
#
#   uniform_S0_050  uniform_S0_075  uniform_S0_095  historical
#
# We fit:
#   - M1 across the three uniform S0 values (S_pop is identical between
#     historical and uniform_S0_050, so historical M1 is not refit separately).
#   - M2 across all four configurations (the historical configuration uses
#     per-serotype S0 derived from observed early-period proportions).
#
# That gives 7 new fits. M0 is the existing fit_model.rds and is not refit.
#
# Output filenames (results/):
#   fit_M1_susceptible_pop_S050.rds
#   fit_M1_susceptible_pop_S075.rds
#   fit_M1_susceptible_pop_S095.rds
#   fit_M2_susceptible_dom_S050.rds
#   fit_M2_susceptible_dom_S075.rds
#   fit_M2_susceptible_dom_S095.rds
#   fit_M2_susceptible_dom_historical.rds
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)

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

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("FIT M1 + M2 SUSCEPTIBLE MODELS (S0 sweep)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. Load data
# ==============================================================================

cat("Loading data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data
stan_data$log_alpha_mu <- -1.2
stan_data$log_rho_mu   <- log(6)

susc <- readRDS("../data/susceptible_covariates.rds")

cat(sprintf("  N_model = %d weeks; K_climate = %d\n",
            stan_data$N_model, stan_data$K_climate))
cat(sprintf("  Configurations available: %s\n",
            paste(names(susc$configs), collapse = ", ")))

# ==============================================================================
# 2. Compile shared Stan model (one binary serves M1 and M2 across all configs)
# ==============================================================================

cat("\nCompiling 14_model_climate_susceptible.stan ...\n")
model <- cmdstan_model("14_model_climate_susceptible.stan")

# ==============================================================================
# 3. Fitting helper
# ==============================================================================

fit_one <- function(label, model_letter, s_vec, output_filename) {
  cat("\n")
  cat("=" |> rep(70) |> paste(collapse = ""), "\n")
  cat(sprintf("%s [%s]\n", model_letter, label))
  cat("=" |> rep(70) |> paste(collapse = ""), "\n")

  sd <- stan_data
  sd$K_S <- 1L
  sd$X_S <- matrix(s_vec, ncol = 1)

  t_start <- Sys.time()
  fit <- model$sample(
    data = sd,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95,
    max_treedepth = 12,
    seed = 42,
    refresh = 200
  )
  cat(sprintf("\n%s elapsed: %s\n", model_letter,
              format(round(difftime(Sys.time(), t_start, units = "mins"), 2))))

  fit$save_object(file.path("../results", output_filename))
  cat(sprintf("Saved: results/%s\n", output_filename))

  diag <- fit$diagnostic_summary()
  cat(sprintf("  Divergent: %d, max-treedepth: %d\n",
              sum(diag$num_divergent), sum(diag$num_max_treedepth)))

  print(fit$summary(variables = c("alpha", "rho", "beta_S",
                                  "prop_susceptible", "prop_residual")))
  invisible(fit)
}

# ==============================================================================
# 4. Run M1 sweep (S_pop, three uniform S0 values; historical S_pop = uniform_S0_050)
# ==============================================================================

m1_sweep <- list(
  list(label = "uniform_S0_050",
       s_vec = susc$configs$uniform_S0_050$log_S_pop_std,
       file  = "fit_M1_susceptible_pop_S050.rds"),
  list(label = "uniform_S0_075",
       s_vec = susc$configs$uniform_S0_075$log_S_pop_std,
       file  = "fit_M1_susceptible_pop_S075.rds"),
  list(label = "uniform_S0_095",
       s_vec = susc$configs$uniform_S0_095$log_S_pop_std,
       file  = "fit_M1_susceptible_pop_S095.rds")
)

for (cfg in m1_sweep) {
  fit_one(cfg$label, "M1: climate + log(S_pop) + GP",
          cfg$s_vec, cfg$file)
}

# ==============================================================================
# 5. Run M2 sweep (S_dom, three uniform + historical)
# ==============================================================================

m2_sweep <- list(
  list(label = "uniform_S0_050",
       s_vec = susc$configs$uniform_S0_050$log_S_dom_std,
       file  = "fit_M2_susceptible_dom_S050.rds"),
  list(label = "uniform_S0_075",
       s_vec = susc$configs$uniform_S0_075$log_S_dom_std,
       file  = "fit_M2_susceptible_dom_S075.rds"),
  list(label = "uniform_S0_095",
       s_vec = susc$configs$uniform_S0_095$log_S_dom_std,
       file  = "fit_M2_susceptible_dom_S095.rds"),
  list(label = "historical (per-serotype S0 from early proportions)",
       s_vec = susc$configs$historical$log_S_dom_std,
       file  = "fit_M2_susceptible_dom_historical.rds")
)

for (cfg in m2_sweep) {
  fit_one(cfg$label, "M2: climate + log(S_dom) + GP",
          cfg$s_vec, cfg$file)
}

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("FIT COMPLETE. Next: Rscript 14_compare_susceptible_models.R\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
