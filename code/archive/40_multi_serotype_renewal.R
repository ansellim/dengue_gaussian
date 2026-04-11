#!/usr/bin/env Rscript
# ==============================================================================
# 40_multi_serotype_renewal.R
#
# Multi-serotype renewal equation model for dengue in Singapore (2013-2022)
#
# Mechanistic alternative to GP + renewal equation:
#   - Splits total cases into serotype-specific streams
#   - Tracks susceptible depletion per serotype
#   - Estimates serotype-specific Rt via Cori method
#   - Predicts Rt from immunity dynamics (R0 * S_i)
#   - Compares aggregate mechanistic Rt with GP-estimated Rt
#
# Sensitivity analysis over expansion factor (EF = 6, 8, 10, 14)
#
# Input:
#   data/model_data.rds
#   data/serotype_props_2013_2016.csv
#   data/monthly_sero_type_props_all_data.csv
#   results/fit_model3.rds
#
# Output:
#   results/figures/multi_serotype_rt_trajectories.png
#   results/figures/multi_serotype_aggregate_comparison.png
#   results/figures/multi_serotype_sensitivity.png
#   results/figures/multi_serotype_susceptible_trajectories.png
#   results/multi_serotype_summary.csv
# ==============================================================================

library(tidyverse)
library(lubridate)
library(mgcv)
library(cmdstanr)
library(posterior)
library(patchwork)
library(zoo)

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

# Serotype colours used throughout
sero_cols <- c("DENV-1" = "#E41A1C", "DENV-2" = "#377EB8",
               "DENV-3" = "#4DAF4A", "DENV-4" = "#984EA3")

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MULTI-SEROTYPE RENEWAL EQUATION MODEL\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")

model_data <- readRDS("../data/model_data.rds")
stan_data  <- model_data$stan_data
df_all     <- model_data$df
gi         <- stan_data$gi   # PMF of length S = 6
S_gi       <- stan_data$S    # = 6

dates_all  <- df_all$date
cases_all  <- df_all$cases
N_all      <- length(cases_all)

cat(sprintf("  Total weekly series: %d weeks (%s to %s)\n",
            N_all, min(dates_all), max(dates_all)))
cat(sprintf("  Generation interval PMF (S=%d): [%s]\n",
            S_gi, paste(round(gi, 3), collapse = ", ")))

# ==============================================================================
# 2. LOAD AND SMOOTH SEROTYPE PROPORTIONS
# ==============================================================================

cat("\nLoading serotype data...\n")

SEROTYPE_DIR <- file.path(dirname(getwd()), "data")

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
  arrange(year_month)

cat(sprintf("  Serotype data: %d months (%s to %s)\n",
            nrow(serotype), min(serotype$year_month), max(serotype$year_month)))

# GAM-smooth each serotype proportion independently, then renormalize
cat("  Smoothing serotype proportions via GAM...\n")

sero_names <- c("D1_prop", "D2_prop", "D3_prop", "D4_prop")

serotype <- serotype |>
  mutate(t_numeric = as.numeric(year_month - min(year_month)) / 30.44)

smooth_mat <- matrix(NA_real_, nrow = nrow(serotype), ncol = 4)
for (j in 1:4) {
  p <- pmin(pmax(serotype[[sero_names[j]]], 0.001), 0.999)
  fit_j <- gam(p ~ s(t_numeric, k = 20, bs = "tp"),
               family = quasibinomial(), data = serotype)
  smooth_mat[, j] <- predict(fit_j, type = "response")
  cat(sprintf("    %s GAM deviance explained: %.1f%%\n",
              sero_names[j], summary(fit_j)$dev.expl * 100))
}

# Renormalize rows to sum to 1
row_sums <- rowSums(smooth_mat)
smooth_mat <- smooth_mat / row_sums

serotype$D1_smooth <- smooth_mat[, 1]
serotype$D2_smooth <- smooth_mat[, 2]
serotype$D3_smooth <- smooth_mat[, 3]
serotype$D4_smooth <- smooth_mat[, 4]

# ==============================================================================
# 3. INTERPOLATE SEROTYPE PROPORTIONS TO WEEKLY
# ==============================================================================

cat("\nInterpolating serotype proportions to weekly resolution...\n")

# Restrict to 2013-2022 period (serotype data availability)
sero_start <- as.Date("2013-01-01")
sero_end   <- as.Date("2022-12-31")

idx_weekly <- which(dates_all >= sero_start & dates_all <= sero_end)
dates_weekly <- dates_all[idx_weekly]
cases_weekly <- cases_all[idx_weekly]
N_weekly <- length(dates_weekly)

cat(sprintf("  Analysis period: %d weeks (%s to %s)\n",
            N_weekly, min(dates_weekly), max(dates_weekly)))

# For each serotype, interpolate smoothed monthly proportions to weekly dates
# using linear interpolation from the monthly mid-points
sero_month_dates <- serotype$year_month + 14  # approx mid-month

prop_weekly <- matrix(NA_real_, nrow = N_weekly, ncol = 4)
for (j in 1:4) {
  prop_weekly[, j] <- approx(
    x = as.numeric(sero_month_dates),
    y = smooth_mat[, j],
    xout = as.numeric(dates_weekly),
    rule = 2  # extend constant beyond boundaries
  )$y
}

# Renormalize weekly proportions
row_sums_w <- rowSums(prop_weekly)
prop_weekly <- prop_weekly / row_sums_w

colnames(prop_weekly) <- c("D1", "D2", "D3", "D4")

# ==============================================================================
# 4. LOAD GP Rt FOR COMPARISON
# ==============================================================================

cat("\nLoading GP-estimated Rt from Model 3...\n")

fit <- readRDS("../results/fit_model3.rds")
dates_model <- model_data$metadata$dates_model

# Extract Rt posterior
rt_draws <- fit$draws("Rt", format = "matrix")

rt_gp <- tibble(
  date = dates_model,
  rt_median = apply(rt_draws, 2, median),
  rt_lower  = apply(rt_draws, 2, quantile, 0.025),
  rt_upper  = apply(rt_draws, 2, quantile, 0.975)
)

# Restrict to 2013-2022
rt_gp <- rt_gp |> filter(date >= sero_start, date <= sero_end)

cat(sprintf("  GP Rt: %d weeks in analysis window\n", nrow(rt_gp)))

# Match dates between weekly analysis and GP Rt
common_dates <- intersect(dates_weekly, rt_gp$date)
idx_common_weekly <- match(common_dates, dates_weekly)
idx_common_gp     <- match(common_dates, rt_gp$date)

cat(sprintf("  Matched %d common weeks for comparison\n", length(common_dates)))

# ==============================================================================
# 5. DEFINE CORE ANALYSIS FUNCTION
# ==============================================================================

run_multi_serotype_analysis <- function(
    cases_weekly,
    prop_weekly,
    dates_weekly,
    gi,
    expansion_factor,
    N_pop = 5500000,
    S_init = 0.75,  # initial susceptible fraction per serotype
    window = 4,     # Cori sliding window width
    rt_gp_matched = NULL,  # GP Rt values at matched dates
    common_idx = NULL      # indices into weekly arrays for common dates
) {

  N <- length(cases_weekly)
  S_gi <- length(gi)
  n_sero <- 4

  # --- Step 1: Serotype-specific weekly cases ---
  cases_sero <- matrix(0, nrow = N, ncol = n_sero)
  for (j in 1:n_sero) {
    cases_sero[, j] <- round(cases_weekly * prop_weekly[, j])
  }

  # --- Step 2: Track susceptible depletion per serotype ---
  infections_sero <- cases_sero * expansion_factor
  cum_inf_sero    <- apply(infections_sero, 2, cumsum)

  S_sero <- matrix(NA_real_, nrow = N, ncol = n_sero)
  for (j in 1:n_sero) {
    S_sero[, j] <- pmax(0.01, S_init - cum_inf_sero[, j] / N_pop)
  }

  # --- Step 3: Cori Rt estimation per serotype (4-week sliding window) ---
  # Infectiousness potential: IP_i[t] = sum(cases_i[t-s] * gi[s], s=1..S)
  IP_sero <- matrix(0, nrow = N, ncol = n_sero)
  for (j in 1:n_sero) {
    for (t in (S_gi + 1):N) {
      for (s in 1:S_gi) {
        IP_sero[t, j] <- IP_sero[t, j] + cases_sero[t - s, j] * gi[s]
      }
    }
  }

  # Sliding window Cori estimate with Gamma posterior
  # Gamma(1 + sum_cases_window, 0.2 + sum_IP_window)
  rt_cori <- matrix(NA_real_, nrow = N, ncol = n_sero)
  rt_cori_lower <- matrix(NA_real_, nrow = N, ncol = n_sero)
  rt_cori_upper <- matrix(NA_real_, nrow = N, ncol = n_sero)

  start_t <- max(S_gi + 1, window)
  for (j in 1:n_sero) {
    for (t in start_t:N) {
      t_start <- max(S_gi + 1, t - window + 1)
      sum_cases <- sum(cases_sero[t_start:t, j])
      sum_IP    <- sum(IP_sero[t_start:t, j])

      alpha_post <- 1 + sum_cases
      beta_post  <- 0.2 + sum_IP

      if (sum_IP > 0.01) {
        rt_cori[t, j]       <- alpha_post / beta_post  # posterior mean
        rt_cori_lower[t, j] <- qgamma(0.025, shape = alpha_post, rate = beta_post)
        rt_cori_upper[t, j] <- qgamma(0.975, shape = alpha_post, rate = beta_post)
      }
    }
  }

  # --- Step 4: Estimate R0 per serotype ---
  # R0_i = max smoothed Rt_i during dominant periods (prop > 10%)
  R0_est <- numeric(n_sero)
  for (j in 1:n_sero) {
    dominant_mask <- prop_weekly[, j] > 0.10 & !is.na(rt_cori[, j])
    if (sum(dominant_mask) > 10) {
      # Smooth Rt for finding R0
      rt_smooth <- rollmean(rt_cori[dominant_mask, j], k = min(8, sum(dominant_mask)),
                            fill = NA, align = "center")
      R0_est[j] <- max(rt_smooth, na.rm = TRUE)
    } else {
      R0_est[j] <- NA
    }
  }

  # Use common R0 as fallback for any serotype with too little data
  R0_common <- median(R0_est, na.rm = TRUE)
  R0_est[is.na(R0_est)] <- R0_common

  cat(sprintf("    EF=%d: R0 estimates = [%.2f, %.2f, %.2f, %.2f]\n",
              expansion_factor, R0_est[1], R0_est[2], R0_est[3], R0_est[4]))

  # --- Step 5: Predicted Rt from immunity ---
  rt_predicted_sero <- matrix(NA_real_, nrow = N, ncol = n_sero)
  for (j in 1:n_sero) {
    rt_predicted_sero[, j] <- R0_est[j] * S_sero[, j]
  }

  # --- Step 6: Aggregate predicted Rt ---
  rt_predicted_agg <- rowSums(prop_weekly * rt_predicted_sero)

  # --- Step 7: Effective susceptible fraction (weighted) ---
  S_effective <- rowSums(prop_weekly * S_sero)

  # --- Step 8: Compare with GP Rt ---
  correlation <- NA
  rmse <- NA

  if (!is.null(rt_gp_matched) && !is.null(common_idx)) {
    rt_pred_common <- rt_predicted_agg[common_idx]
    valid <- !is.na(rt_pred_common) & !is.na(rt_gp_matched)
    if (sum(valid) > 10) {
      correlation <- cor(rt_pred_common[valid], rt_gp_matched[valid])
      rmse <- sqrt(mean((rt_pred_common[valid] - rt_gp_matched[valid])^2))
    }
  }

  list(
    cases_sero = cases_sero,
    S_sero = S_sero,
    S_effective = S_effective,
    IP_sero = IP_sero,
    rt_cori = rt_cori,
    rt_cori_lower = rt_cori_lower,
    rt_cori_upper = rt_cori_upper,
    rt_predicted_sero = rt_predicted_sero,
    rt_predicted_agg = rt_predicted_agg,
    R0_est = R0_est,
    correlation = correlation,
    rmse = rmse,
    expansion_factor = expansion_factor
  )
}

# ==============================================================================
# 6. RUN ANALYSIS FOR MULTIPLE EXPANSION FACTORS
# ==============================================================================

cat("\nRunning multi-serotype analysis...\n")

EF_values <- c(6, 8, 10, 14)
results_list <- list()

for (ef in EF_values) {
  results_list[[as.character(ef)]] <- run_multi_serotype_analysis(
    cases_weekly = cases_weekly,
    prop_weekly  = prop_weekly,
    dates_weekly = dates_weekly,
    gi           = gi,
    expansion_factor = ef,
    rt_gp_matched = rt_gp$rt_median[idx_common_gp],
    common_idx    = idx_common_weekly
  )
}

# Use EF=10 as the reference scenario
res <- results_list[["10"]]

# ==============================================================================
# 7. FIGURE 1: SEROTYPE-SPECIFIC Rt TRAJECTORIES (4 panels)
# ==============================================================================

cat("\nGenerating Figure 1: Serotype-specific Rt trajectories...\n")

sero_labels <- c("DENV-1", "DENV-2", "DENV-3", "DENV-4")

# Build long data frame for plotting
plot_data_sero <- list()
for (j in 1:4) {
  mask <- prop_weekly[, j] > 0.05  # only where serotype > 5%
  plot_data_sero[[j]] <- tibble(
    date = dates_weekly,
    serotype = sero_labels[j],
    rt_cori = res$rt_cori[, j],
    rt_cori_lower = res$rt_cori_lower[, j],
    rt_cori_upper = res$rt_cori_upper[, j],
    rt_predicted = res$rt_predicted_sero[, j],
    prop = prop_weekly[, j],
    visible = mask
  )
}
plot_df_sero <- bind_rows(plot_data_sero) |>
  mutate(
    rt_cori       = ifelse(visible, rt_cori, NA),
    rt_cori_lower = ifelse(visible, rt_cori_lower, NA),
    rt_cori_upper = ifelse(visible, rt_cori_upper, NA),
    rt_predicted  = ifelse(visible, rt_predicted, NA)
  )

fig1 <- ggplot(plot_df_sero, aes(x = date)) +
  geom_ribbon(aes(ymin = rt_cori_lower, ymax = rt_cori_upper),
              fill = "grey70", alpha = 0.4, na.rm = TRUE) +
  geom_line(aes(y = rt_cori, colour = "Cori estimate"), linewidth = 0.5, na.rm = TRUE) +
  geom_line(aes(y = rt_predicted, colour = "R0 x S(t)"),
            linewidth = 0.7, linetype = "dashed", na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey40") +
  facet_wrap(~serotype, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("Cori estimate" = "black",
                                  "R0 x S(t)" = "#D95F02")) +
  coord_cartesian(ylim = c(0, 4)) +
  labs(
    title = "Serotype-specific Rt trajectories (EF = 10)",
    subtitle = "Shown only during periods where serotype proportion > 5%",
    x = NULL, y = expression(R[t]),
    colour = NULL
  ) +
  theme(legend.position = "bottom")

ggsave("../results/figures/multi_serotype_rt_trajectories.png",
       fig1, width = 10, height = 10, dpi = 200, bg = "white")
cat("  Saved: results/figures/multi_serotype_rt_trajectories.png\n")

# ==============================================================================
# 8. FIGURE 2: AGGREGATE COMPARISON WITH GP Rt (2 panels)
# ==============================================================================

cat("Generating Figure 2: Aggregate comparison with GP Rt...\n")

# Top panel: GP Rt vs mechanistic Rt_predicted_agg
agg_df <- tibble(
  date = dates_weekly,
  rt_mech = res$rt_predicted_agg,
  S_eff = res$S_effective
)

# Merge with GP Rt
agg_df <- agg_df |>
  left_join(rt_gp |> select(date, rt_median, rt_lower, rt_upper), by = "date")

# Detect serotype switch points (dominant serotype changes)
dominant_sero <- apply(prop_weekly, 1, which.max)
switch_idx <- which(diff(dominant_sero) != 0)
switch_dates <- dates_weekly[switch_idx]

p_top <- ggplot(agg_df, aes(x = date)) +
  geom_ribbon(aes(ymin = rt_lower, ymax = rt_upper),
              fill = "#1B9E77", alpha = 0.2, na.rm = TRUE) +
  geom_line(aes(y = rt_median, colour = "GP estimate (Model 3)"),
            linewidth = 0.6, na.rm = TRUE) +
  geom_line(aes(y = rt_mech, colour = "Mechanistic (R0 x S)"),
            linewidth = 0.7, na.rm = TRUE) +
  geom_vline(xintercept = switch_dates, linetype = "dashed",
             colour = "grey50", alpha = 0.6) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey40") +
  scale_colour_manual(values = c("GP estimate (Model 3)" = "#1B9E77",
                                  "Mechanistic (R0 x S)" = "#D95F02")) +
  coord_cartesian(ylim = c(0, 3.5)) +
  labs(
    title = sprintf("Aggregate Rt comparison (EF = 10, r = %.2f, RMSE = %.2f)",
                    res$correlation, res$rmse),
    x = NULL, y = expression(R[t]),
    colour = NULL
  ) +
  theme(legend.position = "bottom")

p_bottom <- ggplot(agg_df, aes(x = date, y = S_eff)) +
  geom_line(colour = "#7570B3", linewidth = 0.8) +
  geom_vline(xintercept = switch_dates, linetype = "dashed",
             colour = "grey50", alpha = 0.6) +
  labs(
    subtitle = "Weighted effective susceptible fraction",
    x = NULL, y = expression(S[eff](t))
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme()

fig2 <- p_top / p_bottom + plot_layout(heights = c(2, 1))

ggsave("../results/figures/multi_serotype_aggregate_comparison.png",
       fig2, width = 10, height = 7, dpi = 200, bg = "white")
cat("  Saved: results/figures/multi_serotype_aggregate_comparison.png\n")

# ==============================================================================
# 9. FIGURE 3: SENSITIVITY ANALYSIS OVER EXPANSION FACTOR
# ==============================================================================

cat("Generating Figure 3: Sensitivity to expansion factor...\n")

sens_df <- tibble(
  EF = EF_values,
  correlation = sapply(results_list, function(r) r$correlation),
  RMSE = sapply(results_list, function(r) r$rmse)
)

p_corr <- ggplot(sens_df, aes(x = factor(EF), y = correlation)) +
  geom_col(fill = "#1B9E77", alpha = 0.7, width = 0.5) +
  geom_text(aes(label = sprintf("%.3f", correlation)), vjust = -0.5, size = 3.5) +
  labs(
    title = "Sensitivity of mechanistic Rt to expansion factor",
    x = "Expansion factor (EF)",
    y = "Correlation with GP Rt"
  ) +
  coord_cartesian(ylim = c(0, max(sens_df$correlation, na.rm = TRUE) * 1.15))

p_rmse <- ggplot(sens_df, aes(x = factor(EF), y = RMSE)) +
  geom_col(fill = "#D95F02", alpha = 0.7, width = 0.5) +
  geom_text(aes(label = sprintf("%.3f", RMSE)), vjust = -0.5, size = 3.5) +
  labs(
    x = "Expansion factor (EF)",
    y = "RMSE vs GP Rt"
  ) +
  coord_cartesian(ylim = c(0, max(sens_df$RMSE, na.rm = TRUE) * 1.15))

fig3 <- p_corr / p_rmse

ggsave("../results/figures/multi_serotype_sensitivity.png",
       fig3, width = 7, height = 7, dpi = 200, bg = "white")
cat("  Saved: results/figures/multi_serotype_sensitivity.png\n")

# ==============================================================================
# 10. FIGURE 4: SUSCEPTIBLE TRAJECTORIES BY SEROTYPE
# ==============================================================================

cat("Generating Figure 4: Susceptible trajectories by serotype...\n")

susc_df <- tibble(
  date = rep(dates_weekly, 4),
  serotype = rep(sero_labels, each = N_weekly),
  susceptible = c(res$S_sero[, 1], res$S_sero[, 2],
                  res$S_sero[, 3], res$S_sero[, 4])
)

fig4 <- ggplot(susc_df, aes(x = date, y = susceptible, colour = serotype)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = sero_cols) +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_vline(xintercept = switch_dates, linetype = "dashed",
             colour = "grey50", alpha = 0.5) +
  labs(
    title = "Serotype-specific susceptible fraction (EF = 10)",
    subtitle = sprintf("Initial S = 0.75 per serotype; population = 5.5M; dashed = serotype switch"),
    x = NULL, y = expression(S[i](t)),
    colour = NULL
  ) +
  theme(legend.position = "bottom")

ggsave("../results/figures/multi_serotype_susceptible_trajectories.png",
       fig4, width = 10, height = 5, dpi = 200, bg = "white")
cat("  Saved: results/figures/multi_serotype_susceptible_trajectories.png\n")

# ==============================================================================
# 11. SUMMARY TABLE
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# R0 estimates per serotype (reference EF = 10)
cat("R0 estimates per serotype (EF = 10):\n")
for (j in 1:4) {
  cat(sprintf("  %s: R0 = %.2f\n", sero_labels[j], res$R0_est[j]))
}

cat("\nCorrelation and RMSE vs GP Rt by expansion factor:\n")
for (ef_char in names(results_list)) {
  r <- results_list[[ef_char]]
  cat(sprintf("  EF = %2s: r = %.3f, RMSE = %.3f\n",
              ef_char, r$correlation, r$rmse))
}

# Save summary CSV
summary_rows <- list()
for (ef_char in names(results_list)) {
  r <- results_list[[ef_char]]
  for (j in 1:4) {
    summary_rows[[length(summary_rows) + 1]] <- tibble(
      expansion_factor = as.integer(ef_char),
      serotype = sero_labels[j],
      R0_estimate = r$R0_est[j],
      correlation_with_GP = r$correlation,
      RMSE_with_GP = r$rmse
    )
  }
}

summary_csv <- bind_rows(summary_rows)
write_csv(summary_csv, "../results/multi_serotype_summary.csv")
cat("\nSaved: results/multi_serotype_summary.csv\n")

cat("\nDone.\n")
