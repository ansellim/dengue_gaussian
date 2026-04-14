#!/usr/bin/env Rscript
# ==============================================================================
# 20_prepare_forecast_data.R
#
# Prepare Shannon entropy covariate for the forecast model.
# Computes weekly lagged entropy from monthly serotype proportions.
#
# Input:
#   data/model_data.rds
#   data/serotype_props_2013_2016.csv (monthly, 2013-2016)
#   data/monthly_sero_type_props_all_data.csv (monthly via weekly rows, 2017-2022)
#
# Output:
#   data/forecast_entropy.rds
# ==============================================================================

library(tidyverse)
library(lubridate)

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
cat("PREPARE FORECAST ENTROPY COVARIATE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD MODEL DATA
# ==============================================================================

cat("Loading model data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data
dates_model <- model_data$metadata$dates_model

N_model <- stan_data$N_model
S <- stan_data$S

cat(sprintf("  N_model = %d weeks\n", N_model))
cat(sprintf("  Date range: %s to %s\n", min(dates_model), max(dates_model)))

# ==============================================================================
# 2. LOAD AND COMBINE SEROTYPE DATA
# ==============================================================================

cat("\nLoading serotype data...\n")

# Locate serotype files
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

serotype <- bind_rows(sero_early, sero_late) |>
  arrange(year_month) |>
  distinct(year_month, .keep_all = TRUE)

cat(sprintf("  Serotype data: %d months (%s to %s)\n",
            nrow(serotype), min(serotype$year_month), max(serotype$year_month)))

# ==============================================================================
# 3. COMPUTE MONTHLY SHANNON ENTROPY
# ==============================================================================

cat("\nComputing Shannon entropy...\n")

# Shannon entropy: H = -sum(p_i * log(p_i)), treating 0*log(0) = 0
serotype <- serotype |>
  mutate(
    entropy = -(
      ifelse(D1_prop > 0, D1_prop * log(D1_prop), 0) +
      ifelse(D2_prop > 0, D2_prop * log(D2_prop), 0) +
      ifelse(D3_prop > 0, D3_prop * log(D3_prop), 0) +
      ifelse(D4_prop > 0, D4_prop * log(D4_prop), 0)
    )
  )

cat(sprintf("  Entropy range: %.3f to %.3f\n",
            min(serotype$entropy), max(serotype$entropy)))
cat(sprintf("  Max possible (4 serotypes): %.3f\n", log(4)))

# ==============================================================================
# 4. INTERPOLATE MONTHLY ENTROPY TO WEEKLY
# ==============================================================================

cat("\nInterpolating to weekly resolution...\n")

# Convert monthly dates to numeric for interpolation
sero_numeric <- as.numeric(serotype$year_month)
dates_numeric <- as.numeric(dates_model)

# Use spline interpolation for smooth weekly values
interp <- spline(
  x = sero_numeric,
  y = serotype$entropy,
  xout = dates_numeric,
  method = "natural"
)

entropy_weekly <- interp$y

# Clamp to valid range [0, log(4)]
entropy_weekly <- pmax(entropy_weekly, 0)
entropy_weekly <- pmin(entropy_weekly, log(4))

# Fill any NAs at edges with nearest non-NA value
if (any(is.na(entropy_weekly))) {
  cat("  Warning: NAs in interpolated entropy, filling with nearest values\n")
  # Forward fill then backward fill
  for (i in 2:length(entropy_weekly)) {
    if (is.na(entropy_weekly[i])) entropy_weekly[i] <- entropy_weekly[i - 1]
  }
  for (i in (length(entropy_weekly) - 1):1) {
    if (is.na(entropy_weekly[i])) entropy_weekly[i] <- entropy_weekly[i + 1]
  }
}

cat(sprintf("  Weekly entropy: %d values\n", length(entropy_weekly)))

# ==============================================================================
# 5. CREATE LAGGED ENTROPY
# ==============================================================================

cat("\nCreating lagged entropy...\n")

# Lag of 43 weeks (~10 months), matching CCF peak from serotype analysis
ENTROPY_LAG <- 43

# Shift entropy forward by ENTROPY_LAG weeks
# entropy_lagged[i] = entropy_weekly[i - ENTROPY_LAG]
entropy_lagged <- rep(NA_real_, N_model)

for (i in 1:N_model) {
  if (i > ENTROPY_LAG) {
    entropy_lagged[i] <- entropy_weekly[i - ENTROPY_LAG]
  }
}

# Fill first ENTROPY_LAG weeks with mean of available lagged values
mean_entropy <- mean(entropy_lagged, na.rm = TRUE)
entropy_lagged[is.na(entropy_lagged)] <- mean_entropy

cat(sprintf("  Lag: %d weeks (~%.0f months)\n", ENTROPY_LAG, ENTROPY_LAG / 4.33))
cat(sprintf("  Fill value for first %d weeks: %.3f (mean of lagged values)\n",
            ENTROPY_LAG, mean_entropy))

# ==============================================================================
# 6. STANDARDIZE USING FULL-SERIES STATISTICS
# ==============================================================================

cat("\nStandardizing lagged entropy...\n")

# Store raw statistics for reference
entropy_mean <- mean(entropy_lagged)
entropy_sd <- sd(entropy_lagged)

entropy_std <- (entropy_lagged - entropy_mean) / entropy_sd

cat(sprintf("  Mean: %.3f, SD: %.3f\n", entropy_mean, entropy_sd))
cat(sprintf("  Standardized range: [%.2f, %.2f]\n",
            min(entropy_std), max(entropy_std)))

# ==============================================================================
# 7. SAVE OUTPUT
# ==============================================================================

cat("\nSaving output...\n")

output <- list(
  entropy_covariate = entropy_std,
  metadata = list(
    lag_weeks = ENTROPY_LAG,
    entropy_mean = entropy_mean,
    entropy_sd = entropy_sd,
    n_weeks = N_model,
    dates = dates_model,
    description = "Weekly Shannon entropy of serotype proportions, lagged by 43 weeks, standardized",
    created = Sys.time()
  )
)

saveRDS(output, "../data/forecast_entropy.rds")
cat("  Saved: data/forecast_entropy.rds\n")

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DONE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")
cat("Next step: Run 21_run_forecasts.R\n")
