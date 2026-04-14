#!/usr/bin/env Rscript
# ==============================================================================
# 41_prepare_susceptible_covariate.R
#
# Compute effective susceptible fraction S_eff[t] from serotype + case data
# for use as a covariate in the GP renewal equation model.
#
# S_eff[t] = sum_i( prop_i[t] * S_i[t] )
# where S_i[t] tracks depletion of serotype-i susceptibles over time.
#
# Input:
#   data/model_data.rds
#   data/serotype_props_2013_2016.csv
#   data/monthly_sero_type_props_all_data.csv
#
# Output:
#   data/susceptible_covariate.rds
# ==============================================================================

library(tidyverse)
library(lubridate)
library(mgcv)

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

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PREPARE SUSCEPTIBLE FRACTION COVARIATE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD MODEL DATA
# ==============================================================================

cat("Loading model data...\n")
model_data <- readRDS("../data/model_data.rds")
df <- model_data$df
dates_model <- model_data$metadata$dates_model
N_model <- length(dates_model)
S <- model_data$stan_data$S
N <- model_data$stan_data$N

cat(sprintf("  N_model = %d weeks, date range: %s to %s\n",
            N_model, min(dates_model), max(dates_model)))

# Weekly cases for the full dataset (including burn-in)
cases_all <- df$cases
dates_all <- df$date

# ==============================================================================
# 2. LOAD AND COMBINE SEROTYPE DATA
# ==============================================================================

cat("\nLoading serotype data...\n")

SEROTYPE_DIR <- file.path(dirname(getwd()), "data")
if (!file.exists(file.path(SEROTYPE_DIR, "serotype_props_2013_2016.csv"))) {
  SEROTYPE_DIR <- file.path(dirname(dirname(getwd())), "DengueFever", "data", "nea")
}
stopifnot(
  "Serotype data directory not found." = dir.exists(SEROTYPE_DIR)
)

# 2013-2016: monthly proportions
sero_early <- read_csv(
  file.path(SEROTYPE_DIR, "serotype_props_2013_2016.csv"),
  show_col_types = FALSE
) |>
  mutate(year_month = as.Date(paste0(Month, "-01"))) |>
  select(year_month, D1_prop, D2_prop, D3_prop, D4_prop)

# 2017-2022: weekly rows with monthly resolution -- take one row per month
sero_late <- read_csv(
  file.path(SEROTYPE_DIR, "monthly_sero_type_props_all_data.csv"),
  show_col_types = FALSE
) |>
  mutate(year_month = as.Date(paste0(Month, "-01"))) |>
  group_by(year_month) |>
  slice(1) |>
  ungroup() |>
  select(year_month, D1_prop, D2_prop, D3_prop, D4_prop)

serotype_monthly <- bind_rows(sero_early, sero_late) |>
  arrange(year_month)

cat(sprintf("  Serotype data: %d months (%s to %s)\n",
            nrow(serotype_monthly),
            min(serotype_monthly$year_month),
            max(serotype_monthly$year_month)))

# ==============================================================================
# 3. GAM-SMOOTH SEROTYPE PROPORTIONS AND INTERPOLATE TO WEEKLY
# ==============================================================================

cat("\nSmoothing serotype proportions via GAM and interpolating to weekly...\n")

sero_names <- c("D1_prop", "D2_prop", "D3_prop", "D4_prop")

serotype_monthly <- serotype_monthly |>
  mutate(t_numeric = as.numeric(year_month - min(year_month)) / 30.44)

# Fit GAM for each serotype
gam_models <- list()
for (j in 1:4) {
  p <- pmin(pmax(serotype_monthly[[sero_names[j]]], 0.001), 0.999)
  gam_models[[j]] <- gam(p ~ s(t_numeric, k = 20, bs = "tp"),
                          family = quasibinomial(), data = serotype_monthly)
  cat(sprintf("  %s GAM deviance explained: %.1f%%\n",
              sero_names[j], summary(gam_models[[j]])$dev.expl * 100))
}

# Create weekly time grid aligned to model dates (full dataset dates)
weekly_df <- tibble(date = dates_all) |>
  mutate(
    # Map each week to the serotype monthly time scale
    t_numeric = as.numeric(date - min(serotype_monthly$year_month)) / 30.44
  )

# Predict smoothed proportions at weekly resolution
# For dates outside serotype range (2012), extrapolation will be handled below
sero_range_start <- min(serotype_monthly$year_month)
sero_range_end <- max(serotype_monthly$year_month) + months(1) - days(1)

prop_weekly <- matrix(NA_real_, nrow = nrow(weekly_df), ncol = 4)
for (j in 1:4) {
  pred <- predict(gam_models[[j]], newdata = weekly_df, type = "response")
  prop_weekly[, j] <- pred
}

# Renormalize rows to sum to 1
row_sums <- rowSums(prop_weekly)
prop_weekly <- prop_weekly / row_sums

colnames(prop_weekly) <- c("D1", "D2", "D3", "D4")

cat(sprintf("  Interpolated to %d weekly time points\n", nrow(prop_weekly)))

# ==============================================================================
# 4. COMPUTE SEROTYPE-SPECIFIC SUSCEPTIBLE FRACTIONS
# ==============================================================================

cat("\nComputing serotype-specific susceptible fractions...\n")

# Parameters
expansion_factor <- 8    # Tan et al. central estimate
N_pop <- 5500000         # Singapore population
initial_S_i <- 0.75      # ~50% overall seroprevalence -> ~75% per serotype

cat(sprintf("  Expansion factor: %d\n", expansion_factor))
cat(sprintf("  Population: %s\n", format(N_pop, big.mark = ",")))
cat(sprintf("  Initial S per serotype: %.2f\n", initial_S_i))

# For each serotype: cumulative infections and susceptible depletion
N_weeks <- nrow(weekly_df)
S_i <- matrix(NA_real_, nrow = N_weeks, ncol = 4)
colnames(S_i) <- c("S_D1", "S_D2", "S_D3", "S_D4")

for (j in 1:4) {
  # Estimated infections of serotype j each week
  infections_j <- cases_all * prop_weekly[, j] * expansion_factor

  # Cumulative infections

  cum_infections_j <- cumsum(infections_j)

  # Susceptible fraction for serotype j
  S_i[, j] <- pmax(0.01, initial_S_i - cum_infections_j / N_pop)
}

cat("  Serotype susceptible fractions at start/end:\n")
for (j in 1:4) {
  cat(sprintf("    %s: %.3f -> %.3f\n",
              colnames(S_i)[j], S_i[1, j], S_i[N_weeks, j]))
}

# ==============================================================================
# 5. COMPUTE EFFECTIVE SUSCEPTIBLE FRACTION
# ==============================================================================

cat("\nComputing effective susceptible fraction S_eff...\n")

# S_eff[t] = sum_i(prop_i[t] * S_i[t])
S_eff_all <- rowSums(prop_weekly * S_i)

cat(sprintf("  S_eff range: %.3f to %.3f\n", min(S_eff_all), max(S_eff_all)))

# ==============================================================================
# 6. EXTRACT FOR MODEL PERIOD AND HANDLE 2012 (NO SEROTYPE DATA)
# ==============================================================================

cat("\nAligning to model period and handling 2012...\n")

# Model period starts at index S+1 (after burn-in)
S_eff_model <- S_eff_all[(S + 1):N]
S_i_model <- S_i[(S + 1):N, ]
prop_weekly_model <- prop_weekly[(S + 1):N, ]

# Identify weeks in 2012 (before serotype data starts in 2013-01)
is_2012 <- year(dates_model) < 2013

# For 2012 weeks, S_eff is based on GAM extrapolation which may be unreliable.
# Replace with mean S_eff from the serotype-data period (2013+).
S_eff_2013_plus <- S_eff_model[!is_2012]
mean_S_eff <- mean(S_eff_2013_plus)

cat(sprintf("  Weeks in 2012 (pre-serotype data): %d\n", sum(is_2012)))
cat(sprintf("  Mean S_eff from 2013+: %.4f\n", mean_S_eff))

S_eff_model[is_2012] <- mean_S_eff

# ==============================================================================
# 7. COMPUTE LOG AND STANDARDIZE
# ==============================================================================

cat("\nComputing log(S_eff) and standardizing...\n")

log_S_eff <- log(S_eff_model)
log_S_eff_mean <- mean(log_S_eff)
log_S_eff_sd <- sd(log_S_eff)
log_S_eff_std <- (log_S_eff - log_S_eff_mean) / log_S_eff_sd

cat(sprintf("  log(S_eff): mean = %.4f, sd = %.4f\n", log_S_eff_mean, log_S_eff_sd))
cat(sprintf("  Standardized range: [%.2f, %.2f]\n", min(log_S_eff_std), max(log_S_eff_std)))

# ==============================================================================
# 8. SAVE OUTPUT
# ==============================================================================

cat("\nSaving output...\n")

output <- list(
  log_S_eff_std = log_S_eff_std,
  S_eff_raw = S_eff_model,
  S_i = S_i_model,
  prop_weekly = prop_weekly_model,
  dates = dates_model,
  metadata = list(
    expansion_factor = expansion_factor,
    N_pop = N_pop,
    initial_S = initial_S_i,
    log_S_eff_mean = log_S_eff_mean,
    log_S_eff_sd = log_S_eff_sd,
    n_2012_filled = sum(is_2012),
    mean_S_eff_fill = mean_S_eff,
    created = Sys.time()
  )
)

saveRDS(output, "../data/susceptible_covariate.rds")
cat("  Saved: data/susceptible_covariate.rds\n")

# Quick diagnostic plot
cat("\n  S_eff summary by year:\n")
tibble(date = dates_model, S_eff = S_eff_model) |>
  mutate(year = year(date)) |>
  group_by(year) |>
  summarize(mean_S = mean(S_eff), min_S = min(S_eff), max_S = max(S_eff),
            .groups = "drop") |>
  print()

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUSCEPTIBLE COVARIATE PREPARATION COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("\nNext step: Run 41_fit_and_compare.R to fit the model with S_eff covariate\n")
