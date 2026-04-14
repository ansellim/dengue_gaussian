#!/usr/bin/env Rscript
# ==============================================================================
# 37_susceptible_reconstruction.R
#
# Reconstruct the effective susceptible fraction from the GP residual and
# compare against seroprevalence data (Tan et al. 2019, Am J Epidemiol).
#
# Hypothesis: f_residual tracks the effective susceptible fraction S/N.
# Biological link: log(Rt_eff) = log(Rt_int) + log(S/N), so if climate +
# intrinsic transmissibility are captured by mu + f_climate, then
# f_residual ~ log(S/N) + noise.
#
# Input:
#   results/fit_model3.rds (GP model fit)
#   data/model_data.rds (dates, cases)
#   Serotype data (via same paths as 08_serotype_analysis.R)
#
# Output:
#   results/figures/susceptible_index_trajectory.png
#   results/figures/susceptible_vs_seroprevalence.png
#   results/figures/susceptible_depletion_model.png
#   results/figures/susceptible_by_serotype_period.png
#   results/susceptible_reconstruction_summary.csv
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)
library(lubridate)
library(mgcv)

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
cat("SUSCEPTIBLE FRACTION RECONSTRUCTION FROM GP RESIDUAL\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA AND MODEL FIT
# ==============================================================================

cat("Loading data and model fit...\n")

model_data <- readRDS("../data/model_data.rds")
df <- model_data$df
dates <- model_data$metadata$dates_model
N_model <- length(dates)

if (file.exists("../results/fit_model3.rds")) {
  fit <- readRDS("../results/fit_model3.rds")
  cat("  Using Model 3 (climate-only) fit\n")
  model_label <- "Model 3 (climate-only)"
} else {
  stop("No model fit found. Run 13_fit_model3.R first.")
}

cat(sprintf("  %d modeled weeks (%s to %s)\n", N_model, min(dates), max(dates)))

# ==============================================================================
# 2. EXTRACT f_residual POSTERIOR
# ==============================================================================

cat("\nExtracting f_residual posterior...\n")

f_residual_draws <- fit$draws("f_residual", format = "matrix")
n_draws <- nrow(f_residual_draws)
cat(sprintf("  %d posterior draws x %d weeks\n", n_draws, N_model))

# Posterior summary per week
f_residual_summary <- tibble(
  date      = dates,
  median    = apply(f_residual_draws, 2, median),
  lower_95  = apply(f_residual_draws, 2, quantile, 0.025),
  upper_95  = apply(f_residual_draws, 2, quantile, 0.975),
  lower_80  = apply(f_residual_draws, 2, quantile, 0.10),
  upper_80  = apply(f_residual_draws, 2, quantile, 0.90)
)

# ==============================================================================
# 3. DERIVE EFFECTIVE SUSCEPTIBLE INDEX
# ==============================================================================

cat("\nDeriving susceptible index S_index = exp(f_residual)...\n")

# S_index draws: exponentiate each draw
S_index_draws <- exp(f_residual_draws)

S_index_summary <- tibble(
  date       = dates,
  S_median   = apply(S_index_draws, 2, median),
  S_lower_95 = apply(S_index_draws, 2, quantile, 0.025),
  S_upper_95 = apply(S_index_draws, 2, quantile, 0.975),
  S_lower_80 = apply(S_index_draws, 2, quantile, 0.10),
  S_upper_80 = apply(S_index_draws, 2, quantile, 0.90)
)

# Normalize: S_normalized = S_index / mean(S_index) per draw, then summarize
S_index_mean_per_draw <- rowMeans(S_index_draws)
S_norm_draws <- S_index_draws / S_index_mean_per_draw

S_norm_summary <- tibble(
  date          = dates,
  S_norm_median = apply(S_norm_draws, 2, median),
  S_norm_lower  = apply(S_norm_draws, 2, quantile, 0.025),
  S_norm_upper  = apply(S_norm_draws, 2, quantile, 0.975)
)

cat(sprintf("  S_index range: %.2f - %.2f (median across weeks)\n",
            min(S_index_summary$S_median), max(S_index_summary$S_median)))

# ==============================================================================
# 4. AGGREGATE TO ANNUAL
# ==============================================================================

cat("\nAggregating to annual means...\n")

annual_df <- S_index_summary |>
  mutate(year = year(date)) |>
  filter(year >= 2012, year <= 2022) |>
  group_by(year) |>
  summarize(
    S_annual_mean   = mean(S_median),
    S_annual_lower  = mean(S_lower_95),
    S_annual_upper  = mean(S_upper_95),
    n_weeks         = n(),
    .groups = "drop"
  )

cat("  Annual mean S_index:\n")
for (i in seq_len(nrow(annual_df))) {
  cat(sprintf("    %d: %.3f (%.3f - %.3f), n=%d weeks\n",
              annual_df$year[i], annual_df$S_annual_mean[i],
              annual_df$S_annual_lower[i], annual_df$S_annual_upper[i],
              annual_df$n_weeks[i]))
}

# ==============================================================================
# 5. TAN ET AL. SEROPREVALENCE DATA
# ==============================================================================

cat("\nSetting up Tan et al. seroprevalence reference data...\n")

sero_ref <- tibble(
  year          = c(2004, 2009, 2013, 2017),
  seroprevalence = c(0.45, 0.508, 0.498, 0.486),
  susceptible    = 1 - c(0.45, 0.508, 0.498, 0.486),
  source        = c("Prior serosurvey", "Tan et al.", "Tan et al.", "Tan et al.")
)

cat("  Seroprevalence reference points:\n")
for (i in seq_len(nrow(sero_ref))) {
  cat(sprintf("    %d: seroprevalence %.1f%%, susceptible %.1f%% (%s)\n",
              sero_ref$year[i], sero_ref$seroprevalence[i] * 100,
              sero_ref$susceptible[i] * 100, sero_ref$source[i]))
}

# ==============================================================================
# 6. LOAD SEROTYPE DATA AND IDENTIFY SWITCHES
# ==============================================================================

cat("\nLoading serotype data for switch identification...\n")

SEROTYPE_DIR <- file.path(dirname(getwd()), "data")
if (!file.exists(file.path(SEROTYPE_DIR, "serotype_props_2013_2016.csv"))) {
  SEROTYPE_DIR <- file.path(dirname(dirname(getwd())), "DengueFever", "data", "nea")
}
stopifnot(
  "Serotype data directory not found. Check SEROTYPE_DIR path." =
    dir.exists(SEROTYPE_DIR)
)

# 2013-2016
sero_early <- read_csv(
  file.path(SEROTYPE_DIR, "serotype_props_2013_2016.csv"),
  show_col_types = FALSE
) |>
  mutate(year_month = as.Date(paste0(Month, "-01"))) |>
  select(year_month, D1_prop, D2_prop, D3_prop, D4_prop)

# 2017-2022
sero_late <- read_csv(
  file.path(SEROTYPE_DIR, "monthly_sero_type_props_all_data.csv"),
  show_col_types = FALSE
) |>
  mutate(year_month = as.Date(paste0(Month, "-01"))) |>
  group_by(year_month) |>
  slice(1) |>
  ungroup() |>
  select(year_month, D1_prop, D2_prop, D3_prop, D4_prop)

serotype <- bind_rows(sero_early, sero_late) |> arrange(year_month)

# Smooth serotype proportions via GAM (same approach as 08_serotype_analysis.R)
sero_names <- c("D1_prop", "D2_prop", "D3_prop", "D4_prop")
serotype <- serotype |>
  mutate(t_numeric = as.numeric(year_month - min(year_month)) / 30.44)

smooth_mat <- matrix(NA_real_, nrow = nrow(serotype), ncol = 4)
for (j in 1:4) {
  p <- pmin(pmax(serotype[[sero_names[j]]], 0.001), 0.999)
  fit_j <- gam(p ~ s(t_numeric, k = 20, bs = "tp"),
               family = quasibinomial(), data = serotype)
  smooth_mat[, j] <- predict(fit_j, type = "response")
}
row_sums <- rowSums(smooth_mat)
smooth_mat <- smooth_mat / row_sums

serotype <- serotype |>
  mutate(
    D1_smooth = smooth_mat[, 1], D2_smooth = smooth_mat[, 2],
    D3_smooth = smooth_mat[, 3], D4_smooth = smooth_mat[, 4],
    dominant = case_when(
      D1_smooth >= D2_smooth & D1_smooth >= D3_smooth & D1_smooth >= D4_smooth ~ "DENV-1",
      D2_smooth >= D1_smooth & D2_smooth >= D3_smooth & D2_smooth >= D4_smooth ~ "DENV-2",
      D3_smooth >= D1_smooth & D3_smooth >= D2_smooth & D3_smooth >= D4_smooth ~ "DENV-3",
      TRUE ~ "DENV-4"
    ),
    prev_dominant = lag(dominant)
  )

persistent_switches <- serotype |>
  filter(dominant != prev_dominant, !is.na(prev_dominant)) |>
  select(year_month, from = prev_dominant, to = dominant)

cat(sprintf("  Serotype switches:\n"))
for (i in seq_len(nrow(persistent_switches))) {
  s <- persistent_switches[i, ]
  cat(sprintf("    %s: %s -> %s\n", s$year_month, s$from, s$to))
}

# Build dominance periods table
switch_dates_vec <- persistent_switches$year_month
all_months <- sort(unique(serotype$year_month))

# Create periods from switch dates
period_starts <- c(min(all_months), switch_dates_vec)
period_ends   <- c(switch_dates_vec - months(1), max(all_months))
period_sero   <- character(length(period_starts))

for (k in seq_along(period_starts)) {
  # Find dominant serotype in the middle of the period
  mid_date <- period_starts[k] + floor((period_ends[k] - period_starts[k]) / 2)
  closest <- which.min(abs(serotype$year_month - mid_date))
  period_sero[k] <- serotype$dominant[closest]
}

dominance_periods <- tibble(
  period     = seq_along(period_starts),
  start      = period_starts,
  end        = period_ends,
  serotype   = period_sero
)

cat("\n  Dominance periods:\n")
for (i in seq_len(nrow(dominance_periods))) {
  d <- dominance_periods[i, ]
  cat(sprintf("    Period %d: %s to %s (%s)\n", d$period, d$start, d$end, d$serotype))
}

# ==============================================================================
# 7. CUMULATIVE DEPLETION MODEL
# ==============================================================================

cat("\nBuilding cumulative depletion model...\n")

# Get observed cases aligned with model dates
cases_vec <- df$cases[(length(df$cases) - N_model + 1):length(df$cases)]

# Singapore population ~5.5 million, expansion factor ~6 (Tan et al.)
N_pop <- 5500000
expansion_factor <- 6

# Assign each week to a dominance period
week_period <- rep(NA_integer_, N_model)
for (i in seq_len(nrow(dominance_periods))) {
  in_period <- which(dates >= dominance_periods$start[i] &
                     dates <= (dominance_periods$end[i] + days(15)))
  week_period[in_period] <- i
}
# Fill NAs with nearest period
for (i in which(is.na(week_period))) {
  dists <- abs(as.numeric(dates[i] - dominance_periods$start))
  week_period[i] <- which.min(dists)
}

# Cumulative depletion within each period (reset at serotype switches)
S_depletion <- numeric(N_model)
for (p in unique(week_period)) {
  idx <- which(week_period == p)
  cum_infections <- cumsum(cases_vec[idx] * expansion_factor)
  S_depletion[idx] <- 1 - cum_infections / N_pop
}

depletion_df <- tibble(
  date        = dates,
  cases       = cases_vec,
  S_depletion = S_depletion,
  period      = week_period
)

cat(sprintf("  S_depletion range: %.3f - %.3f\n",
            min(S_depletion), max(S_depletion)))

# ==============================================================================
# 8. FIGURE 1: SUSCEPTIBLE INDEX TRAJECTORY (3-PANEL)
# ==============================================================================

cat("\nGenerating Figure 1: susceptible_index_trajectory.png...\n")

sero_colors <- c(
  "DENV-1" = "#E41A1C", "DENV-2" = "#377EB8",
  "DENV-3" = "#4DAF4A", "DENV-4" = "#984EA3"
)

# Seroprevalence annotation points (only 2013, 2017 within model range)
sero_annot <- sero_ref |>
  filter(year >= year(min(dates)), year <= year(max(dates))) |>
  mutate(date = as.Date(paste0(year, "-07-01")))

# Panel A: f_residual
p_fresid <- ggplot(f_residual_summary, aes(x = date)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
              fill = "steelblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80),
              fill = "steelblue", alpha = 0.4) +
  geom_line(aes(y = median), color = "steelblue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(data = persistent_switches,
             aes(xintercept = year_month),
             linetype = "dotted", color = "red", linewidth = 0.7) +
  labs(
    title = "A. Residual GP (f_residual)",
    y = expression(f[residual] ~ "(log scale)"),
    x = NULL
  ) +
  theme(axis.text.x = element_blank())

# Panel B: S_index = exp(f_residual)
p_sindex <- ggplot(S_index_summary, aes(x = date)) +
  geom_ribbon(aes(ymin = S_lower_95, ymax = S_upper_95),
              fill = "darkorange", alpha = 0.2) +
  geom_ribbon(aes(ymin = S_lower_80, ymax = S_upper_80),
              fill = "darkorange", alpha = 0.4) +
  geom_line(aes(y = S_median), color = "darkorange", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_vline(data = persistent_switches,
             aes(xintercept = year_month),
             linetype = "dotted", color = "red", linewidth = 0.7) +
  # Annotate Tan et al. seroprevalence
  {if (nrow(sero_annot) > 0)
    geom_point(data = sero_annot,
               aes(x = date, y = susceptible / mean(sero_ref$susceptible)),
               shape = 17, size = 3, color = "black")
  } +
  {if (nrow(sero_annot) > 0)
    geom_text(data = sero_annot,
              aes(x = date,
                  y = susceptible / mean(sero_ref$susceptible),
                  label = paste0("Tan et al.\n", year, ": ",
                                 round(susceptible * 100, 1), "% susc.")),
              vjust = -0.8, size = 3)
  } +
  labs(
    title = "B. Susceptible Index (S_index = exp(f_residual))",
    y = "S_index (relative)",
    x = NULL
  ) +
  theme(axis.text.x = element_blank())

# Panel C: observed cases
cases_df <- tibble(date = dates, cases = cases_vec)

p_cases <- ggplot(cases_df, aes(x = date, y = cases)) +
  geom_bar(stat = "identity", fill = "gray60", alpha = 0.7, width = 7) +
  geom_vline(data = persistent_switches,
             aes(xintercept = year_month),
             linetype = "dotted", color = "red", linewidth = 0.7) +
  labs(
    title = "C. Observed Weekly Cases",
    y = "Cases",
    x = "Date"
  )

p_fig1 <- p_fresid / p_sindex / p_cases +
  plot_annotation(
    title = "Effective Susceptible Fraction Reconstruction from GP Residual",
    subtitle = "Red dotted lines = serotype switches; triangles = Tan et al. seroprevalence estimates"
  )

ggsave("../results/figures/susceptible_index_trajectory.png", p_fig1,
       width = 13, height = 10, dpi = 150)
cat("  Saved: results/figures/susceptible_index_trajectory.png\n")

# ==============================================================================
# 9. FIGURE 2: S_INDEX VS SEROPREVALENCE
# ==============================================================================

cat("Generating Figure 2: susceptible_vs_seroprevalence.png...\n")

# Merge annual S_index with seroprevalence where available
annual_with_sero <- annual_df |>
  left_join(sero_ref, by = "year")

p_fig2 <- ggplot() +
  # Annual S_index trajectory
  geom_line(data = annual_df, aes(x = year, y = S_annual_mean),
            color = "darkorange", linewidth = 1) +
  geom_ribbon(data = annual_df,
              aes(x = year, ymin = S_annual_lower, ymax = S_annual_upper),
              fill = "darkorange", alpha = 0.2) +
  geom_point(data = annual_df, aes(x = year, y = S_annual_mean),
             color = "darkorange", size = 3) +
  # Seroprevalence points (rescaled: susceptible fraction / mean to match S_index scale)
  {
    sero_in_range <- sero_ref |> filter(year >= 2012, year <= 2022)
    if (nrow(sero_in_range) > 0) {
      # Scale susceptible fraction to S_index scale using mean ratio
      matched <- annual_with_sero |> filter(!is.na(susceptible))
      if (nrow(matched) > 0) {
        scale_factor <- mean(matched$S_annual_mean) / mean(matched$susceptible)
      } else {
        scale_factor <- mean(annual_df$S_annual_mean) / mean(sero_ref$susceptible)
      }
      list(
        geom_point(data = sero_in_range,
                   aes(x = year, y = susceptible * scale_factor),
                   shape = 17, size = 4, color = "black"),
        geom_text(data = sero_in_range,
                  aes(x = year, y = susceptible * scale_factor,
                      label = paste0(round(susceptible * 100, 1), "%")),
                  vjust = -1.2, size = 3.5)
      )
    }
  } +
  scale_x_continuous(breaks = 2012:2022) +
  labs(
    title = "Annual Susceptible Index vs Tan et al. Seroprevalence Estimates",
    subtitle = "Orange = GP-derived S_index (annual mean); triangles = seroprevalence-based susceptible fraction (rescaled)",
    x = "Year",
    y = "S_index (annual mean)"
  )

ggsave("../results/figures/susceptible_vs_seroprevalence.png", p_fig2,
       width = 10, height = 6, dpi = 150)
cat("  Saved: results/figures/susceptible_vs_seroprevalence.png\n")

# ==============================================================================
# 10. FIGURE 3: SUSCEPTIBLE DEPLETION MODEL (2-PANEL)
# ==============================================================================

cat("Generating Figure 3: susceptible_depletion_model.png...\n")

# Normalize both S_index and S_depletion to [0, 1] range for comparison
S_med <- S_index_summary$S_median
S_med_norm <- (S_med - min(S_med)) / (max(S_med) - min(S_med))
S_dep_norm <- (S_depletion - min(S_depletion)) / (max(S_depletion) - min(S_depletion))

comparison_df <- tibble(
  date       = dates,
  S_gp_norm  = S_med_norm,
  S_dep_norm = S_dep_norm,
  period     = week_period
)

# Panel A: overlay
p_overlay <- ggplot(comparison_df, aes(x = date)) +
  geom_line(aes(y = S_gp_norm, color = "GP-derived S_index"),
            linewidth = 0.8) +
  geom_line(aes(y = S_dep_norm, color = "Cumulative depletion model"),
            linewidth = 0.8) +
  geom_vline(data = persistent_switches,
             aes(xintercept = year_month),
             linetype = "dotted", color = "red", linewidth = 0.7) +
  scale_color_manual(values = c("GP-derived S_index" = "darkorange",
                                "Cumulative depletion model" = "steelblue"),
                     name = NULL) +
  labs(
    title = "A. Normalized susceptible indices: GP-derived vs cumulative case depletion",
    y = "Normalized index (0-1)",
    x = NULL
  ) +
  theme(legend.position = "top",
        axis.text.x = element_blank())

# Panel B: Spearman correlation within each dominance period
period_corrs <- comparison_df |>
  group_by(period) |>
  summarize(
    n_weeks   = n(),
    rho       = if (n() >= 5) cor(S_gp_norm, S_dep_norm, method = "spearman") else NA_real_,
    p_value   = if (n() >= 5) cor.test(S_gp_norm, S_dep_norm, method = "spearman")$p.value else NA_real_,
    start     = min(date),
    end       = max(date),
    .groups   = "drop"
  ) |>
  left_join(dominance_periods |> select(period, serotype), by = "period")

cat("\n  Spearman correlations by dominance period:\n")
for (i in seq_len(nrow(period_corrs))) {
  pc <- period_corrs[i, ]
  cat(sprintf("    Period %d (%s, %s to %s): rho = %.3f, p = %.4f, n = %d\n",
              pc$period, pc$serotype, pc$start, pc$end,
              ifelse(is.na(pc$rho), NA, pc$rho),
              ifelse(is.na(pc$p_value), NA, pc$p_value),
              pc$n_weeks))
}

p_corr <- ggplot(period_corrs |> filter(!is.na(rho)),
                 aes(x = factor(period), y = rho, fill = serotype)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_text(aes(label = sprintf("r=%.2f\n(p=%.3f)", rho, p_value)),
            vjust = ifelse(period_corrs$rho[!is.na(period_corrs$rho)] >= 0, -0.2, 1.2),
            size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = sero_colors, name = "Dominant serotype") +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    title = "B. Spearman correlation: S_index vs cumulative depletion (within each period)",
    x = "Dominance period",
    y = "Spearman rho"
  )

p_fig3 <- p_overlay / p_corr +
  plot_annotation(
    title = "Does cumulative case-driven depletion track the GP-derived susceptible index?",
    subtitle = "Red dotted lines = serotype switches (depletion resets at each switch)"
  )

ggsave("../results/figures/susceptible_depletion_model.png", p_fig3,
       width = 13, height = 9, dpi = 150)
cat("  Saved: results/figures/susceptible_depletion_model.png\n")

# ==============================================================================
# 11. FIGURE 4: S_INDEX BY SEROTYPE PERIOD (START VS END)
# ==============================================================================

cat("Generating Figure 4: susceptible_by_serotype_period.png...\n")

# Compute S_index at start and end of each dominance period
# Use the first and last 8 weeks (approximately 2 months) for stability
window <- 8

period_se <- tibble()
for (p in seq_len(nrow(dominance_periods))) {
  idx <- which(week_period == p)
  if (length(idx) < 2 * window) {
    half <- floor(length(idx) / 2)
    start_idx <- idx[1:half]
    end_idx   <- idx[(half + 1):length(idx)]
  } else {
    start_idx <- idx[1:window]
    end_idx   <- idx[(length(idx) - window + 1):length(idx)]
  }

  # Posterior summaries at start and end of period
  S_start_draws <- rowMeans(S_index_draws[, start_idx, drop = FALSE])
  S_end_draws   <- rowMeans(S_index_draws[, end_idx, drop = FALSE])

  period_se <- bind_rows(period_se, tibble(
    period   = p,
    serotype = dominance_periods$serotype[p],
    phase    = c("Start", "End"),
    S_median = c(median(S_start_draws), median(S_end_draws)),
    S_lower  = c(quantile(S_start_draws, 0.025), quantile(S_end_draws, 0.025)),
    S_upper  = c(quantile(S_start_draws, 0.975), quantile(S_end_draws, 0.975)),
    label    = paste0("P", p, " (", dominance_periods$serotype[p], ")")
  ))
}

period_se$phase <- factor(period_se$phase, levels = c("Start", "End"))

p_fig4 <- ggplot(period_se, aes(x = label, y = S_median, fill = phase)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = S_lower, ymax = S_upper),
                position = position_dodge(width = 0.7), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("Start" = "#2ca02c", "End" = "#d62728"),
                    name = "Phase within period") +
  labs(
    title = "Susceptible Index at Start vs End of Each Serotype Dominance Period",
    subtitle = "Prediction: S_index should be higher at the start (new susceptible pool) and lower at the end (depletion)",
    x = "Dominance period",
    y = "S_index (median with 95% CrI)"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("../results/figures/susceptible_by_serotype_period.png", p_fig4,
       width = 10, height = 6, dpi = 150)
cat("  Saved: results/figures/susceptible_by_serotype_period.png\n")

# ==============================================================================
# 12. SUMMARY STATISTICS
# ==============================================================================

cat("\nCompiling summary statistics...\n")

# Overall correlation between S_index and depletion
overall_rho <- cor(S_med, S_depletion, method = "spearman")
overall_p   <- cor.test(S_med, S_depletion, method = "spearman")$p.value

# Combine summaries
summary_rows <- list()

# Annual mean S_index
for (i in seq_len(nrow(annual_df))) {
  summary_rows <- c(summary_rows, list(tibble(
    metric   = "annual_mean_S_index",
    year     = annual_df$year[i],
    period   = NA_integer_,
    phase    = NA_character_,
    value    = annual_df$S_annual_mean[i],
    lower    = annual_df$S_annual_lower[i],
    upper    = annual_df$S_annual_upper[i],
    note     = ""
  )))
}

# Overall correlation
summary_rows <- c(summary_rows, list(tibble(
  metric = "overall_spearman_rho",
  year   = NA_integer_,
  period = NA_integer_,
  phase  = NA_character_,
  value  = overall_rho,
  lower  = NA_real_,
  upper  = NA_real_,
  note   = sprintf("p = %.4f", overall_p)
)))

# Period-level correlations
for (i in seq_len(nrow(period_corrs))) {
  pc <- period_corrs[i, ]
  summary_rows <- c(summary_rows, list(tibble(
    metric = "period_spearman_rho",
    year   = NA_integer_,
    period = pc$period,
    phase  = NA_character_,
    value  = ifelse(is.na(pc$rho), NA_real_, pc$rho),
    lower  = NA_real_,
    upper  = NA_real_,
    note   = sprintf("%s, p=%.4f, n=%d",
                     pc$serotype,
                     ifelse(is.na(pc$p_value), NA, pc$p_value),
                     pc$n_weeks)
  )))
}

# Start vs end of each period
for (i in seq_len(nrow(period_se))) {
  ps <- period_se[i, ]
  summary_rows <- c(summary_rows, list(tibble(
    metric = "period_start_end_S_index",
    year   = NA_integer_,
    period = ps$period,
    phase  = as.character(ps$phase),
    value  = ps$S_median,
    lower  = ps$S_lower,
    upper  = ps$S_upper,
    note   = ps$serotype
  )))
}

summary_df <- bind_rows(summary_rows)
write_csv(summary_df, "../results/susceptible_reconstruction_summary.csv")
cat("  Saved: results/susceptible_reconstruction_summary.csv\n")

# ==============================================================================
# 13. INTERPRETATION
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("INTERPRETATION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Check if S_index tracks serotype switches
start_end_diff <- period_se |>
  select(period, serotype, phase, S_median) |>
  pivot_wider(names_from = phase, values_from = S_median) |>
  mutate(diff = Start - End)

n_consistent <- sum(start_end_diff$diff > 0, na.rm = TRUE)
n_periods    <- nrow(start_end_diff)

cat(sprintf("1. SEROTYPE SWITCH TRACKING:\n"))
cat(sprintf("   In %d of %d dominance periods, S_index was higher at the start\n",
            n_consistent, n_periods))
cat(sprintf("   than at the end (consistent with susceptible depletion).\n"))
if (n_consistent > n_periods / 2) {
  cat("   -> CONSISTENT with the hypothesis that f_residual tracks immunity dynamics.\n")
} else {
  cat("   -> INCONSISTENT: S_index does not systematically decline within periods.\n")
}

cat(sprintf("\n2. CUMULATIVE DEPLETION CORRELATION:\n"))
cat(sprintf("   Overall Spearman rho = %.3f (p = %.4f)\n", overall_rho, overall_p))
if (abs(overall_rho) > 0.3 && overall_p < 0.05) {
  cat("   -> MODERATE TO STRONG correlation: cumulative cases track the GP residual.\n")
  cat("      f_residual has a plausible mechanistic explanation via susceptible depletion.\n")
} else if (abs(overall_rho) > 0.1) {
  cat("   -> WEAK correlation: some alignment, but other drivers likely contribute.\n")
} else {
  cat("   -> NO meaningful correlation: f_residual likely reflects other drivers\n")
  cat("      (e.g., spatial heterogeneity, behavioral changes, vector dynamics).\n")
}

cat(sprintf("\n3. SEROPREVALENCE COMPARISON:\n"))
cat("   Tan et al. estimate overall susceptible fraction ~50% (stable 2013-2017).\n")
cat("   The GP-derived S_index captures RELATIVE variation around a baseline.\n")
cat("   Serotype-specific susceptibility (relevant for dominant serotype) could be\n")
cat("   much higher than overall 50%, especially after a serotype switch.\n")

cat("\n4. CAVEATS:\n")
cat("   - S_index is a RELATIVE measure, not an absolute susceptible fraction.\n")
cat("   - Expansion factor (6x) and population size are approximate.\n")
cat("   - Cumulative depletion model ignores births, deaths, and waning immunity.\n")
cat("   - Only 2 seroprevalence reference points overlap the model period.\n")

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUSCEPTIBLE RECONSTRUCTION COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\nOutput files:\n")
cat("  results/figures/susceptible_index_trajectory.png     - 3-panel trajectory\n")
cat("  results/figures/susceptible_vs_seroprevalence.png    - Annual vs Tan et al.\n")
cat("  results/figures/susceptible_depletion_model.png      - GP vs depletion model\n")
cat("  results/figures/susceptible_by_serotype_period.png   - Start vs end by period\n")
cat("  results/susceptible_reconstruction_summary.csv       - Summary statistics\n")
