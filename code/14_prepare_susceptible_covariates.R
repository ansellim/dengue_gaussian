#!/usr/bin/env Rscript
# ==============================================================================
# 14_prepare_susceptible_covariates.R
#
# Methodological audit pipeline (M0 / M1 / M2): prepare susceptible-fraction
# covariates for the model-comparison test of whether explicit susceptible
# structure improves on the climate-only + GP baseline.
#
# This is a SEPARATE methodological audit. The headline analysis (climate-only
# model + post-hoc residual GP) deliberately holds serotype information out so
# the residual GP can serve as an independent post-hoc observable. The audit
# breaks that principle in a contained side analysis.
#
# Output: data/susceptible_covariates.rds with several configurations:
#
#   - uniform_S0_050  : S_pop and S_dom with S0 = 0.50 across all pools
#   - uniform_S0_075  : S_pop and S_dom with S0 = 0.75 across all pools
#   - uniform_S0_095  : S_pop and S_dom with S0 = 0.95 across all pools
#   - historical      : S_pop with S0 = overall_seroprev; S_dom with per-serotype
#                       S0_j derived from observed early-period proportions
#                       (this is the headline configuration for M2)
#
# Per-serotype historical S0_j formula:
#   S_j(0) = 1 - (overall_seroprev * mean(prop_j over early window))
# where "early window" is 2013-2014 (the two earliest years of serotype data;
# proxy for pre-period dominance) and overall_seroprev is set from Singapore
# adult dengue serology literature (~50 percent).
#
# This makes the per-serotype starting fractions data-driven rather than a
# uniform guess. The motivation: serotypes that were dominant in the recent
# past have more cumulative immunity (lower S), serotypes that were absent
# have more naive susceptibles (higher S).
#
# All configurations are then standardized to mean 0, sd 1 before being
# written to disk -- so the absolute S0 level barely affects the standardized
# covariate that goes into Stan, and the audit conclusion should be invariant
# across the sweep. The sweep exists to demonstrate that invariance.
#
# Input:
#   data/model_data.rds
#   data/serotype_props_2013_2016.csv
#   data/monthly_sero_type_props_all_data.csv
#
# Output:
#   data/susceptible_covariates.rds
# ==============================================================================

library(tidyverse)
library(lubridate)

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- grep("--file=", args, value = TRUE)
  if (length(script_path) > 0) {
    setwd(dirname(normalizePath(sub("--file=", "", script_path))))
  }
}

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PREPARE SUSCEPTIBLE COVARIATES (M0/M1/M2 AUDIT, S0 SWEEP)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# Fixed parameters (the same across all configurations)
# ==============================================================================

EF <- 10                       # under-reporting expansion factor (Tan et al.)
N_POP <- 5500000               # Singapore population
BIRTHS_PER_YEAR <- 33000
BIRTH_RATE_WEEKLY <- (BIRTHS_PER_YEAR / 52) / N_POP
WANING_PER_YEAR <- 0.003       # Tan et al. 2019
WANING_PER_WEEK <- WANING_PER_YEAR / 52

# Singapore adult dengue overall seroprevalence (any DENV exposure).
# Literature gives ~45-55%; we use 0.50 as a central rough value.
OVERALL_SEROPREV <- 0.50

# S0 sweep values for the uniform configurations
UNIFORM_S0_VALUES <- c(0.50, 0.75, 0.95)

# Time window used to derive per-serotype historical S0
EARLY_WINDOW_YEARS <- 2013:2014

cat(sprintf("  EF = %d, N_pop = %s\n", EF, format(N_POP, big.mark = ",")))
cat(sprintf("  Births = %d/yr, waning = %.4f/yr\n",
            BIRTHS_PER_YEAR, WANING_PER_YEAR))
cat(sprintf("  Overall seroprev assumed = %.2f\n", OVERALL_SEROPREV))
cat(sprintf("  Uniform S0 sweep: {%s}\n",
            paste(sprintf("%.2f", UNIFORM_S0_VALUES), collapse = ", ")))
cat(sprintf("  Historical window: %s\n",
            paste(range(EARLY_WINDOW_YEARS), collapse = "-")))

# ==============================================================================
# 1. Load model data
# ==============================================================================

cat("\nLoading model data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data
dates_model <- as.Date(model_data$metadata$dates_model)
N_model <- length(dates_model)
S_gi <- stan_data$S
N <- stan_data$N
cases_all <- stan_data$cases
cases_model <- cases_all[(S_gi + 1):N]

cat(sprintf("  N_model = %d weeks (%s to %s)\n",
            N_model, min(dates_model), max(dates_model)))

# ==============================================================================
# 2. Load and align serotype data
# ==============================================================================

cat("\nLoading serotype data...\n")

sero_late <- read_csv("../data/monthly_sero_type_props_all_data.csv",
                      show_col_types = FALSE)
sero_early <- read_csv("../data/serotype_props_2013_2016.csv",
                       show_col_types = FALSE)

if ("Month" %in% names(sero_late) && !("month" %in% names(sero_late))) {
  sero_late <- rename(sero_late, month = Month)
}

sero_late_monthly <- sero_late |>
  mutate(month = as.character(month)) |>
  group_by(month) |>
  summarize(D1 = mean(D1_prop, na.rm = TRUE),
            D2 = mean(D2_prop, na.rm = TRUE),
            D3 = mean(D3_prop, na.rm = TRUE),
            D4 = mean(D4_prop, na.rm = TRUE),
            .groups = "drop") |>
  mutate(date = as.Date(paste0(month, "-15"))) |>
  select(date, D1, D2, D3, D4)

month_col <- if ("Month" %in% names(sero_early)) "Month" else "month"
sero_early_monthly <- sero_early |>
  mutate(date = as.Date(paste0(.data[[month_col]], "-15"))) |>
  rename(D1 = D1_prop, D2 = D2_prop, D3 = D3_prop, D4 = D4_prop) |>
  select(date, D1, D2, D3, D4)

sero_monthly <- bind_rows(sero_early_monthly, sero_late_monthly) |>
  distinct(date, .keep_all = TRUE) |>
  arrange(date) |>
  mutate(total = D1 + D2 + D3 + D4,
         D1 = D1 / total, D2 = D2 / total,
         D3 = D3 / total, D4 = D4 / total) |>
  select(-total) |>
  filter(!is.na(D1))

cat(sprintf("  %d months of serotype data (%s to %s)\n",
            nrow(sero_monthly),
            min(sero_monthly$date), max(sero_monthly$date)))

# Interpolate to weekly, aligned to dates_model
sero_weekly <- tibble(date = dates_model)
for (sero in c("D1", "D2", "D3", "D4")) {
  vals <- sero_monthly[[sero]]
  monthly_dates_num <- as.numeric(sero_monthly$date)
  weekly_dates_num <- as.numeric(dates_model)
  in_range <- weekly_dates_num >= min(monthly_dates_num) &
              weekly_dates_num <= max(monthly_dates_num)
  interp <- rep(NA_real_, N_model)
  if (sum(in_range) > 0) {
    interp[in_range] <- approx(monthly_dates_num, vals,
                               xout = weekly_dates_num[in_range],
                               rule = 2)$y
  }
  sero_weekly[[sero]] <- interp
}
sero_weekly <- sero_weekly |>
  mutate(total = D1 + D2 + D3 + D4,
         D1 = ifelse(is.na(total), 0.25, D1 / total),
         D2 = ifelse(is.na(total), 0.25, D2 / total),
         D3 = ifelse(is.na(total), 0.25, D3 / total),
         D4 = ifelse(is.na(total), 0.25, D4 / total)) |>
  select(-total)

is_pre_serotype <- year(dates_model) < min(EARLY_WINDOW_YEARS)

# ==============================================================================
# 3. Compute per-serotype historical S0 from observed early proportions
# ==============================================================================

cat("\nDeriving historical per-serotype S0 from observed early proportions...\n")

# Mean monthly proportions over the early window
early_props <- sero_monthly |>
  filter(year(date) %in% EARLY_WINDOW_YEARS) |>
  summarise(D1 = mean(D1), D2 = mean(D2),
            D3 = mean(D3), D4 = mean(D4))
early_props_vec <- as.numeric(early_props[1, ])
names(early_props_vec) <- c("D1", "D2", "D3", "D4")

cat(sprintf("  Mean serotype proportions over %s:\n",
            paste(range(EARLY_WINDOW_YEARS), collapse = "-")))
for (j in 1:4) {
  cat(sprintf("    %s: %.4f\n", names(early_props_vec)[j], early_props_vec[j]))
}

historical_S0_per_serotype <- 1 - OVERALL_SEROPREV * early_props_vec
cat(sprintf("\n  Per-serotype historical S0 = 1 - %.2f * early_prop_j:\n",
            OVERALL_SEROPREV))
for (j in 1:4) {
  cat(sprintf("    S_%s(0) = %.4f\n",
              names(historical_S0_per_serotype)[j],
              historical_S0_per_serotype[j]))
}

# For S_pop in the historical configuration, use the overall seroprevalence
# as the starting fraction (consistent with the per-serotype derivation)
historical_S0_pop <- 1 - OVERALL_SEROPREV
cat(sprintf("\n  Historical S0_pop = 1 - %.2f = %.4f\n",
            OVERALL_SEROPREV, historical_S0_pop))

# ==============================================================================
# 4. Reconstruction helpers
# ==============================================================================

# Reconstruct S_pop given a starting fraction
reconstruct_S_pop <- function(s0_pop) {
  S_pop <- numeric(N_model)
  S_pop[1] <- s0_pop
  for (t in 2:N_model) {
    infections <- cases_model[t - 1] * EF
    S_pop[t] <- S_pop[t - 1] - infections / N_POP
    S_pop[t] <- S_pop[t] + BIRTH_RATE_WEEKLY
    S_pop[t] <- S_pop[t] + WANING_PER_WEEK * (1 - S_pop[t])
    S_pop[t] <- max(0.01, min(1.0, S_pop[t]))
  }
  S_pop
}

# Reconstruct S_dom given a vector of per-serotype starting fractions
reconstruct_S_dom <- function(s0_vec) {
  stopifnot(length(s0_vec) == 4)
  S_sero <- matrix(0, nrow = N_model, ncol = 4)
  for (j in 1:4) S_sero[1, j] <- s0_vec[j]
  colnames(S_sero) <- c("D1", "D2", "D3", "D4")

  cases_sero <- matrix(0, nrow = N_model, ncol = 4)
  for (j in 1:4) {
    sero_col <- c("D1", "D2", "D3", "D4")[j]
    cases_sero[, j] <- round(cases_model * sero_weekly[[sero_col]])
  }

  for (t in 2:N_model) {
    for (j in 1:4) {
      infections_j <- cases_sero[t - 1, j] * EF
      S_sero[t, j] <- S_sero[t - 1, j] - infections_j / N_POP
      S_sero[t, j] <- S_sero[t, j] + BIRTH_RATE_WEEKLY
      S_sero[t, j] <- S_sero[t, j] + WANING_PER_WEEK * (1 - S_sero[t, j])
    }
    S_sero[t, ] <- pmax(0.01, pmin(1.0, S_sero[t, ]))
  }

  S_dom <- numeric(N_model)
  for (t in 1:N_model) {
    props <- c(sero_weekly$D1[t], sero_weekly$D2[t],
               sero_weekly$D3[t], sero_weekly$D4[t])
    dom_idx <- which.max(props)
    S_dom[t] <- S_sero[t, dom_idx]
  }

  # Pre-2013 fill (no serotype data) with 2013+ mean
  if (any(is_pre_serotype)) {
    mean_post <- mean(S_dom[!is_pre_serotype])
    S_dom[is_pre_serotype] <- mean_post
  }
  list(S_dom = S_dom, S_sero = S_sero)
}

standardize <- function(x) {
  lx <- log(x)
  list(values = lx, mean = mean(lx), sd = sd(lx),
       std = (lx - mean(lx)) / sd(lx))
}

# ==============================================================================
# 5. Build all configurations
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("BUILDING CONFIGURATIONS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

build_config <- function(label, s0_pop, s0_dom_vec) {
  cat(sprintf("\n[%s]\n", label))
  cat(sprintf("  S0_pop = %.4f\n", s0_pop))
  cat(sprintf("  S0_dom_vec = (%s)\n",
              paste(sprintf("%.4f", s0_dom_vec), collapse = ", ")))

  S_pop <- reconstruct_S_pop(s0_pop)
  dom <- reconstruct_S_dom(s0_dom_vec)
  S_dom <- dom$S_dom

  log_S_pop <- standardize(S_pop)
  log_S_dom <- standardize(S_dom)

  cat(sprintf("  S_pop range = [%.4f, %.4f], log range = %.4f, sd = %.4f\n",
              min(S_pop), max(S_pop),
              diff(range(log(S_pop))), log_S_pop$sd))
  cat(sprintf("  S_dom range = [%.4f, %.4f], log range = %.4f, sd = %.4f\n",
              min(S_dom), max(S_dom),
              diff(range(log(S_dom))), log_S_dom$sd))

  list(
    label = label,
    s0_pop = s0_pop,
    s0_dom_vec = s0_dom_vec,
    log_S_pop_std = log_S_pop$std,
    log_S_dom_std = log_S_dom$std,
    S_pop = S_pop,
    S_dom = S_dom,
    S_sero = dom$S_sero,
    metadata = list(
      log_S_pop_mean = log_S_pop$mean,
      log_S_pop_sd = log_S_pop$sd,
      log_S_dom_mean = log_S_dom$mean,
      log_S_dom_sd = log_S_dom$sd
    )
  )
}

configs <- list()
for (s0 in UNIFORM_S0_VALUES) {
  label <- sprintf("uniform_S0_%03d", round(100 * s0))
  configs[[label]] <- build_config(label, s0_pop = s0,
                                   s0_dom_vec = rep(s0, 4))
}
configs[["historical"]] <- build_config(
  "historical",
  s0_pop = historical_S0_pop,
  s0_dom_vec = unname(historical_S0_per_serotype)
)

# ==============================================================================
# 6. Save
# ==============================================================================

output <- list(
  configs = configs,
  dates = dates_model,
  prop_weekly = as.matrix(sero_weekly[, c("D1", "D2", "D3", "D4")]),
  metadata = list(
    EF = EF,
    N_pop = N_POP,
    births_per_year = BIRTHS_PER_YEAR,
    waning_per_year = WANING_PER_YEAR,
    overall_seroprev = OVERALL_SEROPREV,
    uniform_s0_values = UNIFORM_S0_VALUES,
    early_window_years = EARLY_WINDOW_YEARS,
    early_serotype_props = early_props_vec,
    historical_s0_per_serotype = historical_S0_per_serotype,
    historical_s0_pop = historical_S0_pop,
    n_pre_serotype_filled = sum(is_pre_serotype),
    created = Sys.time()
  )
)

saveRDS(output, "../data/susceptible_covariates.rds")
cat("\nSaved: data/susceptible_covariates.rds\n")
cat(sprintf("  Configurations: %s\n", paste(names(configs), collapse = ", ")))

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PREP COMPLETE. Next: Rscript 14_fit_susceptible_models.R\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
