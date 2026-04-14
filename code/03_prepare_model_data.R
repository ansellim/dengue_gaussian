#!/usr/bin/env Rscript
# ==============================================================================
# 03_prepare_model_data.R
#
# Merge, lag, standardize, and format data for Stan models
# Input: data/raw_*.csv files from 01_acquire_data.py
# Output: data/model_data.rds (list ready for Stan)
# ==============================================================================

library(tidyverse)
library(lubridate)

# Set working directory to script location
# Works both in RStudio and command line
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # Command line: assume we're in the project directory
}

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DENGUE RT ESTIMATION - DATA PREPARATION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD RAW DATA
# ==============================================================================

cat("Loading raw data files...\n")

dengue <- read_csv("../data/raw_dengue_cases.csv", show_col_types = FALSE)
weather <- read_csv("../data/raw_weather.csv", show_col_types = FALSE)

cat(sprintf("  Dengue: %d weeks (%s to %s)\n",
            nrow(dengue), min(dengue$date), max(dengue$date)))
cat(sprintf("  Weather: %d weeks\n", nrow(weather)))

# ==============================================================================
# 2. MERGE DATASETS
# ==============================================================================

cat("\nMerging datasets...\n")

# Ensure all dates are Date type
dengue <- dengue |> mutate(date = as.Date(date))
weather <- weather |> mutate(date = as.Date(date))

# Create year_week for joining (handles date alignment issues)
dengue <- dengue |> mutate(year_week = paste0(year(date), "-", sprintf("%02d", isoweek(date))))
weather <- weather |> mutate(year_week = paste0(year(date), "-", sprintf("%02d", isoweek(date))))

# Merge datasets using year_week
df <- dengue |>
  left_join(weather |> select(-date), by = "year_week")

# Filter to complete cases and study period (2012-2022)
df <- df |>
  filter(year(date) >= 2012, year(date) <= 2022) |>
  filter(!is.na(cases), !is.na(temp_mean), !is.na(rainfall_total))

cat(sprintf("  Merged dataset: %d weeks\n", nrow(df)))

# ==============================================================================
# 3. CREATE LAGGED CLIMATE COVARIATES
# ==============================================================================

cat("\nCreating lagged climate covariates...\n")

# Configuration for lag (baseline = 4 weeks)
CLIMATE_LAG <- 4

df <- df |>
  arrange(date) |>
  mutate(
    # Lagged climate variables
    temp_lag = lag(temp_mean, CLIMATE_LAG),
    rain_lag = lag(rainfall_total, CLIMATE_LAG),

    # Also create alternative lags for sensitivity analysis
    temp_lag2 = lag(temp_mean, 2),
    temp_lag3 = lag(temp_mean, 3),
    temp_lag6 = lag(temp_mean, 6),
    rain_lag2 = lag(rainfall_total, 2),
    rain_lag3 = lag(rainfall_total, 3),
    rain_lag6 = lag(rainfall_total, 6)
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
    rain_std = (rain_lag - rain_mean_val) / rain_sd_val,

    # Alternative lags, standardized with same parameters
    temp_std_lag2 = (temp_lag2 - temp_mean_val) / temp_sd_val,
    temp_std_lag3 = (temp_lag3 - temp_mean_val) / temp_sd_val,
    temp_std_lag6 = (temp_lag6 - temp_mean_val) / temp_sd_val,
    rain_std_lag2 = (rain_lag2 - rain_mean_val) / rain_sd_val,
    rain_std_lag3 = (rain_lag3 - rain_mean_val) / rain_sd_val,
    rain_std_lag6 = (rain_lag6 - rain_mean_val) / rain_sd_val
  )

# ASSUMPTIONS:
# 1. 100% case ascertainment rate. Singapore has mandatory dengue notification
#    and active laboratory surveillance, so under-reporting is expected to be
#    minimal. The model does not include an ascertainment fraction parameter.
#    Note: asymptomatic infections (~75% globally; Bhatt et al. 2013) are not
#    captured. Rt reflects symptomatic-case-to-symptomatic-case transmission.
# 2. Constant reporting delay (~5-7 days onset-to-notification) absorbed into
#    the generation interval and negative binomial overdispersion. This is
#    standard for weekly dengue Rt estimation (Cori et al. 2013; Lau et al. 2022).

cat(sprintf("  Temperature: mean=%.2f, sd=%.2f\n", temp_mean_val, temp_sd_val))
cat(sprintf("  Rainfall: mean=%.1f, sd=%.1f\n", rain_mean_val, rain_sd_val))

# ==============================================================================
# 5. COMPUTE GENERATION INTERVAL PMF
# ==============================================================================

cat("\nComputing generation interval from component distributions...\n")
cat("  (Monte Carlo convolution of sub-process distributions)\n")

# The dengue generation interval (infection-to-infection) is the convolution of:
#
#   GI = T_iip + T_viremia + T_eip + T_bite
#
# Component 1 - Intrinsic incubation period (T_iip):
#   Time from human infection to symptom onset / viremia onset.
#   Gamma; mean 5.9 days, SD 1.5 days.
#   Source: Chan & Johansson (2012, PLOS ONE) - systematic review (n=153):
#     "The mean IIP estimate was 5.9 days, with 95% expected between days 3
#     and 10." Gamma fit: shape=16, rate=1.78, 95% range [3.4, 9.0] days.
#
# Component 2 - Viremia bite timing (T_viremia):
#   Time within the viremic period until a mosquito acquires the virus.
#   Uniform(0, D) where D = viremia duration (~5 days).
#   Sources:
#     CDC dengue training: "People are infectious to mosquitoes from 2 days
#       before, to 5 days after illness onset."
#     PMC reviews: "Most people are viremic for about 4-5 days."
#     Chan & Johansson (2012): "most individuals become infectious within a
#       day before or after the onset of disease."
#   Note: viremia begins ~1d before symptom onset, so T_iip + T_viremia
#   slightly overestimates by ~0.5d. Negligible at weekly resolution.
#
# Component 3 - Extrinsic incubation period (T_eip):
#   Time for DENV to replicate in the mosquito and reach salivary glands.
#   Gamma; mean depends on temperature (the dominant source of GI variation).
#   Source: Chan & Johansson (2012, PLOS ONE) - meta-analysis (n=146), log-
#     normal model with temperature covariate (betaT = -0.08 per degree C):
#     "At 25C: mean 15 days (95% CI: 10, 20)" with 95% range [5, 33] days.
#     "At 30C: mean 6.5 days (95% CI: 4.8, 8.8)" with 95% range [2.4, 15] days.
#     Interpolating at 27C: mean ~11 days.
#   Corroborated by Rohani et al. (2009): "The virus was first detected on
#     Day 9 at 26C and 28C and on Day 5 at 30C for both dengue 2 and 4."
#     (first-detection times; means are longer)
#   SD set to give CV ~0.35 (between single-population variation and the
#     meta-analytic CV of ~0.48 which includes inter-study heterogeneity).
#
# Component 4 - Mosquito biting delay (T_bite):
#   Time from mosquito becoming infectious to biting a susceptible human.
#   Exponential; mean ~2 days.
#   Reasoning: Gonotrophic cycle of Ae. aegypti at 27C is ~7 days (Delatte
#     et al. 2009: "mean duration of 7.10 (+/-3.38) days at 27C"). Ae. aegypti
#     feeds multiple times per cycle (Scott et al. 2000, J Med Entomol: 42%
#     double meals, 5% triple meals). Effective blood-feeding interval is thus
#     ~4-5 days. Time from EIP completion to next bite averages roughly half
#     this interval: ~2 days.
#
# This replaces the previous approach of fitting a single gamma to published
# serial interval estimates. Advantages:
#   1. Uses infection-to-infection timing (what the renewal equation requires),
#      not symptom-to-symptom serial intervals
#   2. Each component draws on a well-characterized literature distribution
#   3. Sensitivity analysis maps to biologically meaningful temperature
#      variation via the EIP (the dominant source of GI uncertainty)

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

  # Discretize to weekly PMF: P(week k) = P((k-1)*7 < GI <= k*7)
  gi_weeks <- gi_days / 7
  pmf <- numeric(max_weeks)
  for (k in 1:max_weeks) {
    pmf[k] <- mean(gi_weeks > (k - 1) & gi_weeks <= k)
  }
  tail_beyond <- mean(gi_weeks > max_weeks)
  pmf <- pmf / sum(pmf)

  list(
    pmf = pmf,
    mean_days = mean(gi_days),
    sd_days = sd(gi_days),
    mean_weeks = mean(gi_days) / 7,
    sd_weeks = sd(gi_days) / 7,
    q025 = unname(quantile(gi_days, 0.025)),
    q975 = unname(quantile(gi_days, 0.975)),
    tail_beyond_max = tail_beyond,
    components = data.frame(
      component = c("IIP", "Viremia timing", "EIP", "Bite delay"),
      mean_days = c(mean(t_iip), mean(t_viremia), mean(t_eip), mean(t_bite)),
      sd_days = c(sd(t_iip), sd(t_viremia), sd(t_eip), sd(t_bite))
    )
  )
}

# --- Baseline: empirical mean weekly temperature across the study period ---
# Chan & Johansson (2012) report EIP ~15d at 25C and ~6.5d at 30C; we linearly
# interpolate (slope = -1.7 d/C) at the actual study-period mean rather than at
# an idealised round-number reference, so the GI matches the data we actually
# have. Rohani et al. (2009) corroborate: first detection at day 9 at 26C/28C
# (means are longer than first detection).
empirical_mean_temp_C <- mean(df$temp_mean, na.rm = TRUE)
eip_slope <- (6.5 - 15) / (30 - 25)        # days per degree C
eip_at_mean <- 15 + eip_slope * (empirical_mean_temp_C - 25)

cat(sprintf("\n  Empirical mean weekly temperature: %.2f C\n", empirical_mean_temp_C))
cat(sprintf("  Interpolated EIP at empirical mean: %.2f days (vs idealised 11.0 d at 27 C)\n",
            eip_at_mean))

gi_base <- compute_gi_from_components(eip_mean = eip_at_mean, eip_sd = 4)
gi_baseline <- gi_base$pmf

cat(sprintf("\n  Component breakdown (baseline, %.1f C):\n", empirical_mean_temp_C))
for (i in 1:nrow(gi_base$components)) {
  cat(sprintf("    %-16s mean = %4.1f days,  SD = %3.1f days\n",
              gi_base$components$component[i],
              gi_base$components$mean_days[i],
              gi_base$components$sd_days[i]))
}
cat(sprintf("\n  Total GI: mean = %.1f days (%.2f wk), SD = %.1f days (%.2f wk)\n",
            gi_base$mean_days, gi_base$mean_weeks,
            gi_base$sd_days, gi_base$sd_weeks))
cat(sprintf("  95%% interval: [%.1f, %.1f] days\n", gi_base$q025, gi_base$q975))
cat(sprintf("  Weekly PMF: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n",
            gi_baseline[1], gi_baseline[2], gi_baseline[3],
            gi_baseline[4], gi_baseline[5], gi_baseline[6]))
cat(sprintf("  Probability beyond %d weeks: %.2f%%\n",
            length(gi_baseline), gi_base$tail_beyond_max * 100))

# --- Sensitivity: Short GI (warmer periods ~30C -> shorter EIP) ---
# Chan & Johansson (2012): "At 30C: mean 6.5 days (95% CI: 4.8, 8.8)"
gi_short_res <- compute_gi_from_components(eip_mean = 6.5, eip_sd = 2.5)
gi_short <- gi_short_res$pmf
cat(sprintf("\n  Short GI (EIP=6.5d, ~30C): mean = %.1f days (%.2f wk)\n",
            gi_short_res$mean_days, gi_short_res$mean_weeks))
cat(sprintf("    PMF: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n",
            gi_short[1], gi_short[2], gi_short[3],
            gi_short[4], gi_short[5], gi_short[6]))

# --- Sensitivity: Long GI (cooler periods ~25C -> longer EIP) ---
# Chan & Johansson (2012): "At 25C: mean 15 days (95% CI: 10, 20)"
gi_long_res <- compute_gi_from_components(eip_mean = 15, eip_sd = 5)
gi_long <- gi_long_res$pmf
cat(sprintf("  Long GI (EIP=15d, ~25C): mean = %.1f days (%.2f wk)\n",
            gi_long_res$mean_days, gi_long_res$mean_weeks))
cat(sprintf("    PMF: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n",
            gi_long[1], gi_long[2], gi_long[3],
            gi_long[4], gi_long[5], gi_long[6]))

# ==============================================================================
# 6. PREPARE DATA FOR STAN
# ==============================================================================

cat("\nPreparing Stan data structure...\n")

# Number of time points
N <- nrow(df)

# Maximum generation interval (burn-in period)
# S=6 to accommodate the long-GI sensitivity scenario (EIP=15d at ~25C,
# total GI mean ~25d with 97.5th percentile ~37d)
S <- 6

# Number of modeled time points (after burn-in)
N_model <- N - S

# HSGP parameters
M <- 200  # Number of HSGP basis functions
          # Required: M >= ceil(sqrt(3)/rho_min * 2*L/pi)
          # At 5th percentile of rho prior with median 6 weeks
          # (rho_min ~ 2.6 weeks): M ~ 183. Use 200 for safety.
L_factor <- 1.5  # Boundary factor

# Time indices for GP (in weeks, centered but not rescaled)
t_model <- (1:N_model)
t_center <- mean(t_model)
t_centered <- t_model - t_center  # Centered at 0, in units of weeks

# Boundary for HSGP (in weeks)
L <- L_factor * (max(t_model) - min(t_model)) / 2

# Case counts
cases <- df$cases

# Covariates (for modeled time points)
X_climate <- cbind(
  temp = df$temp_std[(S + 1):N],
  rain = df$rain_std[(S + 1):N]
)

# Stan data list (climate-only model)
stan_data <- list(
  # Dimensions
  N = N,
  N_model = N_model,
  S = S,
  M = M,

  # Observations
  cases = cases,

  # Time (centered, in weeks)
  t = t_centered,
  L = L,

  # Generation interval
  gi = gi_baseline,

  # Covariates (climate only: temperature + rainfall)
  K_climate = 2,
  X_climate = X_climate
)

# Also prepare sensitivity analysis data
stan_data_sens <- list(
  gi_short = gi_short,
  gi_long = gi_long,
  X_climate_lag2 = cbind(
    temp = df$temp_std_lag2[(S + 1):N],
    rain = df$rain_std_lag2[(S + 1):N]
  ),
  X_climate_lag3 = cbind(
    temp = df$temp_std_lag3[(S + 1):N],
    rain = df$rain_std_lag3[(S + 1):N]
  ),
  X_climate_lag6 = cbind(
    temp = df$temp_std_lag6[(S + 1):N],
    rain = df$rain_std_lag6[(S + 1):N]
  )
)

# ==============================================================================
# 7. SAVE OUTPUT
# ==============================================================================

cat("\nSaving output...\n")

# Save main data
output <- list(
  df = df,
  stan_data = stan_data,
  stan_data_sens = stan_data_sens,
  metadata = list(
    climate_lag = CLIMATE_LAG,
    temp_mean = temp_mean_val,
    temp_sd = temp_sd_val,
    rain_mean = rain_mean_val,
    rain_sd = rain_sd_val,
    dates_model = df$date[(S + 1):N],
    gi_components = gi_base$components,
    gi_summary = list(
      mean_days = gi_base$mean_days,
      sd_days = gi_base$sd_days,
      q025 = gi_base$q025,
      q975 = gi_base$q975
    ),
    created = Sys.time()
  )
)

saveRDS(output, "../data/model_data.rds")
cat("  Saved: data/model_data.rds\n")

# Also save as CSV for inspection
df |>
  select(date, cases, temp_std, rain_std) |>
  write_csv("../data/model_data.csv")
cat("  Saved: data/model_data.csv\n")

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DATA PREPARATION COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat(sprintf("\nDataset summary:
  Total weeks: %d
  Modeled weeks (after burn-in): %d
  Date range: %s to %s

  Case counts:
    Min: %.0f
    Median: %.0f
    Max: %.0f

  Covariates (climate only, standardized):
    Temperature (lag %d weeks): mean=0, sd=1
    Rainfall (lag %d weeks): mean=0, sd=1

  Assumptions:
    100%% case ascertainment (mandatory notification)
    Reporting delay (~5-7d) absorbed at weekly resolution

  HSGP configuration:
    Basis functions (M): %d
    Boundary factor: %.1f

  Generation interval (component-based, baseline ~27C):
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
            min(cases), median(cases), max(cases),
            CLIMATE_LAG, CLIMATE_LAG,
            M, L_factor,
            gi_base$mean_days, gi_base$mean_weeks,
            gi_baseline[1] * 100, gi_baseline[2] * 100,
            gi_baseline[3] * 100, gi_baseline[4] * 100,
            gi_baseline[5] * 100, gi_baseline[6] * 100))

cat("\nNext step: Fit Stan models using 06_fit_model.R\n")
