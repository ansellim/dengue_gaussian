#!/usr/bin/env Rscript
# ==============================================================================
# 31_prepare_brazil_data.R
#
# Merge, lag, standardize, and format Brazil dengue data for Stan models.
# Mirrors the exact stan_data structure from 02_prepare_model_data.R so the
# same Stan model (15_model3_multikernel.stan) can be used without changes.
#
# Input:  data/brazil_raw_cases.csv, data/brazil_raw_weather.csv
# Output: data/brazil_model_data.rds
# ==============================================================================

library(tidyverse)
library(lubridate)

# Set working directory to script location
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
cat("BRAZIL DENGUE - DATA PREPARATION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD RAW DATA
# ==============================================================================

cat("Loading raw data files...\n")

cases <- read_csv("../data/brazil_raw_cases.csv", show_col_types = FALSE)
weather <- read_csv("../data/brazil_raw_weather.csv", show_col_types = FALSE)

cat(sprintf("  Cases: %d weeks (%s to %s)\n",
            nrow(cases), min(cases$date), max(cases$date)))
cat(sprintf("  Weather: %d weeks\n", nrow(weather)))

# ==============================================================================
# 2. MERGE DATASETS
# ==============================================================================

cat("\nMerging datasets...\n")

# Ensure dates are Date type
cases <- cases |> mutate(date = as.Date(date))
weather <- weather |> mutate(date = as.Date(date))

# Merge cases and weather on date
df <- cases |>
  left_join(weather |> select(-any_of("cases")), by = "date")

# Filter to complete cases and study period (2012-2022)
df <- df |>
  filter(year(date) >= 2012, year(date) <= 2022) |>
  filter(!is.na(cases), !is.na(temp_mean), !is.na(rainfall_total))

cat(sprintf("  Merged dataset: %d weeks\n", nrow(df)))

# ==============================================================================
# 3. CREATE LAGGED CLIMATE COVARIATES
# ==============================================================================

cat("\nCreating lagged climate covariates...\n")

# Same lag as Singapore (4 weeks) for comparability
CLIMATE_LAG <- 4

df <- df |>
  arrange(date) |>
  mutate(
    temp_lag = lag(temp_mean, CLIMATE_LAG),
    rain_lag = lag(rainfall_total, CLIMATE_LAG)
  )

# Remove rows with missing lagged values
df <- df |> filter(!is.na(temp_lag), !is.na(rain_lag))

cat(sprintf("  After applying lag: %d weeks\n", nrow(df)))

# ==============================================================================
# 4. STANDARDIZE COVARIATES
# ==============================================================================

cat("\nStandardizing covariates...\n")

# Standardize climate variables (z-score)
temp_mean_val <- mean(df$temp_lag, na.rm = TRUE)
temp_sd_val <- sd(df$temp_lag, na.rm = TRUE)
rain_mean_val <- mean(df$rain_lag, na.rm = TRUE)
rain_sd_val <- sd(df$rain_lag, na.rm = TRUE)

df <- df |>
  mutate(
    temp_std = (temp_lag - temp_mean_val) / temp_sd_val,
    rain_std = (rain_lag - rain_mean_val) / rain_sd_val
  )

cat(sprintf("  Temperature: mean=%.2f, sd=%.2f\n", temp_mean_val, temp_sd_val))
cat(sprintf("  Rainfall/humidity: mean=%.1f, sd=%.1f\n", rain_mean_val, rain_sd_val))

# ==============================================================================
# 5. COMPUTE GENERATION INTERVAL PMF
# ==============================================================================

cat("\nComputing generation interval from component distributions...\n")
cat("  (Same dengue biology as Singapore; using baseline ~27C EIP)\n")

# Identical function from 02_prepare_model_data.R
# The generation interval is a property of dengue biology and does not change
# by country. The EIP may differ with temperature, but we use the same baseline
# for comparability and let the climate covariates capture temperature effects.

compute_gi_from_components <- function(
    iip_mean = 5.9,
    iip_sd = 1.5,
    viremia_duration = 5,
    eip_mean = 11,
    eip_sd = 4,
    bite_mean = 2,
    max_weeks = 6,
    n_sim = 200000,
    seed = 42
) {
  set.seed(seed)

  # Component 1: Intrinsic incubation period (Gamma)
  iip_shape <- (iip_mean / iip_sd)^2
  iip_rate <- iip_mean / iip_sd^2
  t_iip <- rgamma(n_sim, shape = iip_shape, rate = iip_rate)

  # Component 2: Timing within viremic period (Uniform)
  t_viremia <- runif(n_sim, 0, viremia_duration)

  # Component 3: Extrinsic incubation period (Gamma)
  eip_shape <- (eip_mean / eip_sd)^2
  eip_rate <- eip_mean / eip_sd^2
  t_eip <- rgamma(n_sim, shape = eip_shape, rate = eip_rate)

  # Component 4: Mosquito biting delay (Exponential)
  t_bite <- rexp(n_sim, rate = 1 / bite_mean)

  # Total generation interval (days)
  gi_days <- t_iip + t_viremia + t_eip + t_bite

  # Discretize to weekly PMF
  gi_weeks <- gi_days / 7
  pmf <- numeric(max_weeks)
  for (k in 1:max_weeks) {
    pmf[k] <- mean(gi_weeks > (k - 1) & gi_weeks <= k)
  }
  pmf <- pmf / sum(pmf)

  list(
    pmf = pmf,
    mean_days = mean(gi_days),
    sd_days = sd(gi_days),
    mean_weeks = mean(gi_days) / 7,
    sd_weeks = sd(gi_days) / 7,
    q025 = unname(quantile(gi_days, 0.025)),
    q975 = unname(quantile(gi_days, 0.975))
  )
}

# Baseline GI: same as Singapore (~27C)
gi_base <- compute_gi_from_components(eip_mean = 11, eip_sd = 4)
gi_baseline <- gi_base$pmf

cat(sprintf("  GI mean: %.1f days (%.2f weeks)\n", gi_base$mean_days, gi_base$mean_weeks))
cat(sprintf("  Weekly PMF: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n",
            gi_baseline[1], gi_baseline[2], gi_baseline[3],
            gi_baseline[4], gi_baseline[5], gi_baseline[6]))

# ==============================================================================
# 6. PREPARE DATA FOR STAN
# ==============================================================================

cat("\nPreparing Stan data structure...\n")

# Number of time points
N <- nrow(df)

# Maximum generation interval (burn-in period) — same as Singapore
S <- 6

# Number of modeled time points (after burn-in)
N_model <- N - S

# HSGP parameters — same as Singapore
M <- 200  # Number of HSGP basis functions
L_factor <- 1.5  # Boundary factor

# Time indices for GP (in weeks, centered)
t_model <- (1:N_model)
t_center <- mean(t_model)
t_centered <- t_model - t_center

# Boundary for HSGP (in weeks)
L <- L_factor * (max(t_model) - min(t_model)) / 2

# Case counts
cases_vec <- df$cases

# Covariates (for modeled time points)
X_climate <- cbind(
  temp = df$temp_std[(S + 1):N],
  rain = df$rain_std[(S + 1):N]
)

# Stan data list — identical structure to Singapore model_data
stan_data <- list(
  # Dimensions
  N = N,
  N_model = N_model,
  S = S,
  M = M,

  # Observations
  cases = cases_vec,

  # Time (centered, in weeks)
  t = t_centered,
  L = L,

  # Generation interval
  gi = gi_baseline,

  # Covariates (climate only: temperature + rainfall/humidity)
  K_climate = 2,
  X_climate = X_climate
)

# ==============================================================================
# 7. SAVE OUTPUT
# ==============================================================================

cat("\nSaving output...\n")

output <- list(
  df = df,
  stan_data = stan_data,
  metadata = list(
    city = "São Paulo",
    geocode = 3550308,
    climate_lag = CLIMATE_LAG,
    temp_mean = temp_mean_val,
    temp_sd = temp_sd_val,
    rain_mean = rain_mean_val,
    rain_sd = rain_sd_val,
    dates_model = df$date[(S + 1):N],
    gi_summary = list(
      mean_days = gi_base$mean_days,
      sd_days = gi_base$sd_days,
      q025 = gi_base$q025,
      q975 = gi_base$q975
    ),
    created = Sys.time()
  )
)

saveRDS(output, "../data/brazil_model_data.rds")
cat("  Saved: data/brazil_model_data.rds\n")

# Also save as CSV for inspection
df |>
  select(date, cases, temp_std, rain_std) |>
  write_csv("../data/brazil_model_data.csv")
cat("  Saved: data/brazil_model_data.csv\n")

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DATA PREPARATION COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat(sprintf("\nDataset summary:
  City: São Paulo, Brazil
  Total weeks: %d
  Modeled weeks (after burn-in): %d
  Date range: %s to %s

  Case counts:
    Min: %.0f
    Median: %.0f
    Max: %.0f

  Covariates (standardized):
    Temperature (lag %d weeks): mean=0, sd=1
    Rainfall/humidity (lag %d weeks): mean=0, sd=1

  HSGP configuration:
    Basis functions (M): %d
    Boundary factor: %.1f
    L = %.1f weeks

  Generation interval (baseline ~27C):
    Mean: %.1f days (%.2f weeks)
    Week 1: %.1f%%
    Week 2: %.1f%%
    Week 3: %.1f%%
    Week 4: %.1f%%
    Week 5: %.1f%%
    Week 6: %.1f%%
\n",
            N, N_model,
            min(df$date), max(df$date),
            min(cases_vec), median(cases_vec), max(cases_vec),
            CLIMATE_LAG, CLIMATE_LAG,
            M, L_factor, L,
            gi_base$mean_days, gi_base$mean_weeks,
            gi_baseline[1] * 100, gi_baseline[2] * 100,
            gi_baseline[3] * 100, gi_baseline[4] * 100,
            gi_baseline[5] * 100, gi_baseline[6] * 100))

cat("\nNext step: Run 32_fit_brazil.R to fit the model\n")
