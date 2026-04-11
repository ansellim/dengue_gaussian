#!/usr/bin/env Rscript
# ==============================================================================
# 26_spectral_analysis.R
#
# Spectral analysis of the GP residual (f_residual) to identify dominant
# periodicities in non-climate-driven dengue transmission dynamics.
#
# Input:  results/fit_model3.rds, data/model_data.rds,
#         results/serotype_switch_timing.csv
# Output: results/figures/spectral_periodogram.png
#         results/figures/spectral_with_context.png
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)
library(lubridate)

# Set up paths
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(commandArgs(trailingOnly = FALSE) |>
    grep("--file=", x = _, value = TRUE) |>
    sub("--file=", "", x = _) |>
    normalizePath() |>
    dirname())
}

dir.create("../results/figures", showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank()))

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SPECTRAL ANALYSIS OF GP RESIDUAL\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA AND FIT
# ==============================================================================

cat("Loading data and model fit...\n")

model_data <- readRDS("../data/model_data.rds")
dates <- model_data$metadata$dates_model
N_model <- length(dates)

fit <- readRDS("../results/fit_model3.rds")
cat(sprintf("  %d modeled weeks, %s to %s\n", N_model, min(dates), max(dates)))

# ==============================================================================
# 2. EXTRACT f_residual POSTERIOR DRAWS
# ==============================================================================

cat("\nExtracting f_residual posterior draws...\n")

f_residual_draws <- fit$draws("f_residual", format = "matrix")
n_draws_total <- nrow(f_residual_draws)
cat(sprintf("  %d posterior draws x %d time points\n", n_draws_total, ncol(f_residual_draws)))

# Posterior median f_residual time series
f_residual_median <- apply(f_residual_draws, 2, median)
f_residual_lower <- apply(f_residual_draws, 2, quantile, 0.025)
f_residual_upper <- apply(f_residual_draws, 2, quantile, 0.975)

# ==============================================================================
# 3. COMPUTE POWER SPECTRUM OF MEDIAN f_residual
# ==============================================================================

cat("\nComputing power spectrum...\n")

# Hanning window function
hanning_window <- function(n) {
  0.5 * (1 - cos(2 * pi * seq(0, n - 1) / (n - 1)))
}

# Compute power spectrum for a time series using FFT with Hanning window
compute_spectrum <- function(x) {
  n <- length(x)
  # Apply Hanning window
  w <- hanning_window(n)
  x_windowed <- (x - mean(x)) * w
  # FFT
  X <- fft(x_windowed)
  # One-sided power spectral density (normalise by window energy)
  psd <- (Mod(X)^2) / (sum(w^2) * n)
  # Keep positive frequencies only (exclude DC)
  n_freq <- floor(n / 2)
  freqs <- (1:n_freq) / n  # cycles per week
  power <- 2 * psd[2:(n_freq + 1)]  # one-sided: multiply by 2
  periods <- 1 / freqs  # period in weeks
  list(freq = freqs, period = periods, power = power)
}

# Spectrum of the posterior median
spec_median <- compute_spectrum(f_residual_median)

cat(sprintf("  Frequency resolution: 1/%d = %.4f cycles/week\n", N_model, 1 / N_model))
cat(sprintf("  Maximum resolvable period: %d weeks (%.1f years)\n",
            N_model, N_model / 52))
cat(sprintf("  Minimum resolvable period: 2 weeks\n"))

# ==============================================================================
# 4. SPECTRAL UNCERTAINTY FROM POSTERIOR DRAWS
# ==============================================================================

cat("\nComputing spectral uncertainty from 200 posterior draws...\n")

set.seed(42)
draw_idx <- sample(n_draws_total, min(200, n_draws_total))
n_freq <- length(spec_median$freq)

# Matrix to store power spectra from each draw
spec_matrix <- matrix(NA_real_, nrow = length(draw_idx), ncol = n_freq)

for (i in seq_along(draw_idx)) {
  spec_i <- compute_spectrum(f_residual_draws[draw_idx[i], ])
  spec_matrix[i, ] <- spec_i$power
}

spec_lower <- apply(spec_matrix, 2, quantile, 0.025)
spec_upper <- apply(spec_matrix, 2, quantile, 0.975)

cat("  Done.\n")

# ==============================================================================
# 5. IDENTIFY DOMINANT PERIODS
# ==============================================================================

cat("\nIdentifying dominant periods...\n")

# Find peaks: periods where power exceeds 2x the median spectral power
median_power <- median(spec_median$power)
peak_mask <- spec_median$power > 2 * median_power

# Filter to plausible range (4-260 weeks)
in_range <- spec_median$period >= 4 & spec_median$period <= 260
peak_idx <- which(peak_mask & in_range)

# Find local maxima among peaks
dominant_periods <- c()
if (length(peak_idx) > 1) {
  for (j in peak_idx) {
    is_local_max <- TRUE
    if (j > 1 && spec_median$power[j] < spec_median$power[j - 1]) is_local_max <- FALSE
    if (j < n_freq && spec_median$power[j] < spec_median$power[j + 1]) is_local_max <- FALSE
    if (is_local_max) dominant_periods <- c(dominant_periods, j)
  }
}

# If no strict local maxima, take top peaks above threshold
if (length(dominant_periods) == 0 && length(peak_idx) > 0) {
  top_idx <- peak_idx[order(spec_median$power[peak_idx], decreasing = TRUE)]
  dominant_periods <- top_idx[1:min(3, length(top_idx))]
}

cat(sprintf("  Median spectral power: %.4e\n", median_power))
cat(sprintf("  Peaks above 2x median power: %d\n", length(peak_idx)))

if (length(dominant_periods) > 0) {
  cat("  Dominant periods:\n")
  for (dp in dominant_periods) {
    period_wk <- spec_median$period[dp]
    period_mo <- period_wk / (52 / 12)
    ratio <- spec_median$power[dp] / median_power
    cat(sprintf("    %.0f weeks (%.1f months) -- %.1fx median power\n",
                period_wk, period_mo, ratio))
  }
} else {
  cat("  No dominant peaks found above 2x median power\n")
}

# Reference periods
ref_periods <- tibble(
  period = c(52, 104, 156),
  label = c("1 year\n(52 wk)", "2 years\n(104 wk)", "3 years\n(156 wk)")
)

# ==============================================================================
# 6. FIGURE 1: SPECTRAL PERIODOGRAM
# ==============================================================================

cat("\nGenerating spectral periodogram plot...\n")

spec_df <- tibble(
  period = spec_median$period,
  power = spec_median$power,
  power_lower = spec_lower,
  power_upper = spec_upper
) |>
  filter(period >= 4, period <= 260)

p_periodogram <- ggplot(spec_df, aes(x = period)) +
  geom_ribbon(aes(ymin = power_lower, ymax = power_upper),
              fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = power), color = "steelblue", linewidth = 0.8) +
  geom_vline(data = ref_periods, aes(xintercept = period),
             linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_text(data = ref_periods, aes(x = period, y = max(spec_df$power) * 0.85,
                                     label = label),
            hjust = -0.1, size = 3, color = "gray30") +
  scale_x_log10(
    breaks = c(4, 8, 13, 26, 52, 104, 156, 260),
    labels = c("4", "8", "13", "26", "52", "104", "156", "260"),
    limits = c(4, 260)
  ) +
  scale_y_log10() +
  annotation_logticks(sides = "lb") +
  labs(
    title = "Power Spectrum of GP Residual (f_residual)",
    subtitle = "Hanning-windowed FFT with 95% posterior credible interval",
    x = "Period (weeks, log scale)",
    y = "Spectral density (log scale)"
  )

# Annotate dominant peaks if they exist
if (length(dominant_periods) > 0) {
  peak_df <- tibble(
    period = spec_median$period[dominant_periods],
    power = spec_median$power[dominant_periods]
  ) |>
    filter(period >= 4, period <= 260)

  if (nrow(peak_df) > 0) {
    peak_df <- peak_df |>
      mutate(label = sprintf("%.0f wk", period))

    p_periodogram <- p_periodogram +
      geom_point(data = peak_df, aes(x = period, y = power),
                 color = "red", size = 3) +
      geom_text(data = peak_df, aes(x = period, y = power, label = label),
                vjust = -1, color = "red", size = 3.5, fontface = "bold")
  }
}

ggsave("../results/figures/spectral_periodogram.png", p_periodogram,
       width = 10, height = 6, dpi = 150)
cat("  Saved: results/figures/spectral_periodogram.png\n")

# ==============================================================================
# 7. FIGURE 2: SPECTRAL WITH CONTEXT (2-panel)
# ==============================================================================

cat("\nGenerating spectral context plot...\n")

# Load serotype switches for annotation
serotype_switches <- NULL
if (file.exists("../results/serotype_switch_timing.csv")) {
  serotype_switches <- read_csv("../results/serotype_switch_timing.csv",
                                 show_col_types = FALSE)
  cat(sprintf("  Loaded %d serotype switch events\n", nrow(serotype_switches)))
}

# Top panel: f_residual time series
ts_df <- tibble(
  date = dates,
  median = f_residual_median,
  lower = f_residual_lower,
  upper = f_residual_upper
)

p_timeseries <- ggplot(ts_df, aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = median), color = "steelblue", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")

# Add serotype switch lines if available
if (!is.null(serotype_switches) && "switch_date" %in% colnames(serotype_switches)) {
  switch_dates <- as.Date(serotype_switches$switch_date)
  switch_labels <- serotype_switches$switch_event
  switch_df <- tibble(date = switch_dates, label = switch_labels) |>
    filter(date >= min(dates), date <= max(dates))

  if (nrow(switch_df) > 0) {
    p_timeseries <- p_timeseries +
      geom_vline(data = switch_df, aes(xintercept = date),
                 linetype = "dotted", color = "red", linewidth = 0.6) +
      geom_text(data = switch_df,
                aes(x = date, y = max(ts_df$upper) * 0.9, label = label),
                angle = 90, vjust = -0.3, size = 2.5, color = "red")
  }
}

p_timeseries <- p_timeseries +
  labs(
    title = "GP Residual Time Series",
    subtitle = "Posterior median and 95% CrI; red lines = serotype switches",
    x = "Date",
    y = expression(f[residual])
  )

# Bottom panel: periodogram (reuse from above)
p_periodogram_bottom <- p_periodogram +
  labs(title = "Power Spectrum of GP Residual")

# Combine with patchwork
p_combined <- p_timeseries / p_periodogram_bottom +
  plot_annotation(
    title = "Spectral Analysis of Non-Climate Dengue Dynamics",
    subtitle = "Singapore, 2012-2022",
    tag_levels = "A"
  )

ggsave("../results/figures/spectral_with_context.png", p_combined,
       width = 12, height = 10, dpi = 150)
cat("  Saved: results/figures/spectral_with_context.png\n")

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SPECTRAL ANALYSIS SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

if (length(dominant_periods) > 0) {
  cat("Dominant periods (power > 2x median spectral density):\n")
  for (dp in dominant_periods) {
    period_wk <- spec_median$period[dp]
    period_yr <- period_wk / 52
    ratio <- spec_median$power[dp] / median_power
    cat(sprintf("  %.0f weeks (%.1f years) -- %.1fx median power\n",
                period_wk, period_yr, ratio))
  }

  # Interpret peaks
  cat("\nInterpretation:\n")
  for (dp in dominant_periods) {
    period_wk <- spec_median$period[dp]
    if (period_wk >= 45 && period_wk <= 60) {
      cat(sprintf("  ~%.0f weeks: consistent with annual seasonal cycling\n", period_wk))
    } else if (period_wk >= 90 && period_wk <= 130) {
      cat(sprintf("  ~%.0f weeks: consistent with 2-year serotype cycling\n", period_wk))
    } else if (period_wk >= 130 && period_wk <= 180) {
      cat(sprintf("  ~%.0f weeks: consistent with 2-3 year multi-annual cycling\n", period_wk))
    } else if (period_wk >= 180) {
      cat(sprintf("  ~%.0f weeks: low-frequency trend (>3 years)\n", period_wk))
    } else {
      cat(sprintf("  ~%.0f weeks: sub-annual periodicity\n", period_wk))
    }
  }
} else {
  cat("No dominant peaks found above 2x median spectral density.\n")
  cat("The residual GP may exhibit broad-spectrum variability\n")
  cat("rather than periodic oscillations.\n")
}

cat("\nKnown dengue cycling timescales for reference:\n")
cat("  ~52 weeks (1 year): seasonal transmission cycle\n")
cat("  ~104-156 weeks (2-3 years): serotype replacement cycling\n")
cat("  ~208 weeks (4 years): four-serotype rotation period\n")

cat("\n")
cat("Output files:\n")
cat("  results/figures/spectral_periodogram.png\n")
cat("  results/figures/spectral_with_context.png\n")
