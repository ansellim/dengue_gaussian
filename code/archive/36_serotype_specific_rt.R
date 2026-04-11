#!/usr/bin/env Rscript
# ==============================================================================
# 36_serotype_specific_rt.R
#
# Serotype-specific Rt estimation using the Cori et al. (2013) method
#
# Biological question: Does DENV-3 transmit faster than DENV-1 or DENV-2
# in Singapore? DENV-3 emerged in 2020 after decades of low circulation,
# entering a largely susceptible population.
#
# Approach:
#   1. Approximate serotype-specific weekly cases by multiplying total cases
#      by GAM-smoothed serotype proportions (interpolated to weekly)
#   2. Estimate instantaneous Rt per serotype via the Cori renewal equation
#      with a 4-week sliding window and Bayesian Gamma conjugate posterior
#   3. Compare Rt distributions across serotypes and dominance periods
#
# Input:
#   data/model_data.rds
#   data/serotype_props_2013_2016.csv
#   data/monthly_sero_type_props_all_data.csv
#
# Output:
#   results/figures/serotype_specific_rt_trajectories.png
#   results/figures/serotype_specific_rt_boxplot.png
#   results/figures/serotype_specific_rt_by_period.png
#   results/figures/serotype_specific_cases.png
#   results/serotype_specific_rt_summary.csv
# ==============================================================================

library(tidyverse)
library(lubridate)
library(mgcv)
library(patchwork)

# Install zoo if not available (for rolling means)
if (!requireNamespace("zoo", quietly = TRUE)) {
  install.packages("zoo", repos = "https://cloud.r-project.org")
}
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

# Serotype color palette (consistent throughout)
sero_colors <- c("DENV-1" = "#E41A1C", "DENV-2" = "#377EB8",
                 "DENV-3" = "#4DAF4A", "DENV-4" = "#984EA3")

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SEROTYPE-SPECIFIC RT ESTIMATION (CORI METHOD)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")

model_data <- readRDS("../data/model_data.rds")
df <- model_data$df
gi <- model_data$stan_data$gi   # Weekly generation interval PMF
S  <- model_data$stan_data$S    # Max GI lag (6 weeks)

cat(sprintf("  Cases: %d weeks (%s to %s)\n",
            nrow(df), min(df$date), max(df$date)))
cat(sprintf("  GI PMF (S=%d): [%s]\n", S,
            paste(round(gi, 3), collapse = ", ")))

# --- Load serotype data (same approach as 08_serotype_analysis.R) ---

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

# ==============================================================================
# 2. GAM-SMOOTH SEROTYPE PROPORTIONS AND INTERPOLATE TO WEEKLY
# ==============================================================================

cat("\nSmoothing serotype proportions via GAM...\n")

# Add numeric time variable for GAM
serotype <- serotype |>
  mutate(t_numeric = as.numeric(year_month - min(year_month)) / 30.44)

sero_names <- c("D1_prop", "D2_prop", "D3_prop", "D4_prop")
smooth_mat <- matrix(NA_real_, nrow = nrow(serotype), ncol = 4)

for (j in 1:4) {
  p <- pmin(pmax(serotype[[sero_names[j]]], 0.001), 0.999)
  fit_j <- gam(p ~ s(t_numeric, k = 20, bs = "tp"),
               family = quasibinomial(), data = serotype)
  smooth_mat[, j] <- predict(fit_j, type = "response")
  cat(sprintf("  %s GAM deviance explained: %.1f%%\n",
              sero_names[j], summary(fit_j)$dev.expl * 100))
}

# Renormalize rows to sum to 1
smooth_mat <- smooth_mat / rowSums(smooth_mat)

serotype <- serotype |>
  mutate(D1_smooth = smooth_mat[, 1], D2_smooth = smooth_mat[, 2],
         D3_smooth = smooth_mat[, 3], D4_smooth = smooth_mat[, 4])

# --- Interpolate monthly smoothed proportions to weekly resolution ---
cat("\nInterpolating to weekly resolution...\n")

# Restrict cases to serotype data availability (2013-2022)
df_sero <- df |>
  filter(year(date) >= 2013, year(date) <= 2022)

# Monthly dates and smoothed proportions as numeric for interpolation
month_dates_num <- as.numeric(serotype$year_month)
week_dates_num  <- as.numeric(df_sero$date)

# Interpolate each serotype proportion to weekly dates
weekly_props <- tibble(date = df_sero$date)
for (j in 1:4) {
  sero_label <- c("D1", "D2", "D3", "D4")[j]
  smooth_col <- c("D1_smooth", "D2_smooth", "D3_smooth", "D4_smooth")[j]
  interp <- approx(x = month_dates_num,
                   y = serotype[[smooth_col]],
                   xout = week_dates_num,
                   rule = 2)  # extend edge values
  weekly_props[[paste0(sero_label, "_prop")]] <- pmax(interp$y, 0)
}

# Renormalize weekly proportions to sum to 1
prop_cols <- c("D1_prop", "D2_prop", "D3_prop", "D4_prop")
row_totals <- rowSums(weekly_props[, prop_cols])
for (col in prop_cols) {
  weekly_props[[col]] <- weekly_props[[col]] / row_totals
}

# Compute serotype-specific cases
weekly_props <- weekly_props |>
  mutate(
    total_cases = df_sero$cases,
    cases_DENV1 = round(total_cases * D1_prop),
    cases_DENV2 = round(total_cases * D2_prop),
    cases_DENV3 = round(total_cases * D3_prop),
    cases_DENV4 = round(total_cases * D4_prop)
  )

cat(sprintf("  Weekly serotype-specific cases: %d weeks\n", nrow(weekly_props)))

# ==============================================================================
# 3. ESTIMATE RT PER SEROTYPE USING CORI METHOD
# ==============================================================================

cat("\nEstimating serotype-specific Rt (Cori method, 4-week sliding window)...\n")

# Cori method with sliding window w:
#   Rt[t] = sum(cases[t-w+1:t]) / sum(IP[t-w+1:t])
# where IP[t] = sum_{s=1}^{S} cases[t-s] * gi[s]
#
# Bayesian conjugate: cases[t] ~ Poisson(Rt * IP[t])
#   Prior: Rt ~ Gamma(a_prior, b_prior)
#   Posterior: Gamma(a_prior + sum_cases, b_prior + sum_IP)

cori_estimate <- function(cases_vec, gi, S, window = 4,
                          a_prior = 1, b_prior = 0.2) {
  n <- length(cases_vec)

  # Step 1: compute infectious pressure at each time point
  IP <- rep(0, n)
  for (t in (S + 1):n) {
    for (s in 1:S) {
      IP[t] <- IP[t] + cases_vec[t - s] * gi[s]
    }
  }

  # Step 2: sliding window sums and Bayesian posterior
  results <- tibble(
    t = 1:n,
    cases = cases_vec,
    IP = IP,
    Rt_median = NA_real_,
    Rt_mean = NA_real_,
    Rt_lower = NA_real_,
    Rt_upper = NA_real_
  )

  start_t <- S + window  # need S weeks burn-in + window weeks
  for (t in start_t:n) {
    idx <- (t - window + 1):t
    sum_cases <- sum(cases_vec[idx])
    sum_IP    <- sum(IP[idx])

    if (sum_IP < 0.01) {
      # Not enough infectious pressure -- Rt undefined
      next
    }

    # Gamma posterior parameters
    a_post <- a_prior + sum_cases
    b_post <- b_prior + sum_IP

    results$Rt_mean[t]   <- a_post / b_post
    results$Rt_median[t] <- qgamma(0.5, shape = a_post, rate = b_post)
    results$Rt_lower[t]  <- qgamma(0.025, shape = a_post, rate = b_post)
    results$Rt_upper[t]  <- qgamma(0.975, shape = a_post, rate = b_post)
  }

  results
}

# Estimate Rt for each serotype
serotype_labels <- c("DENV-1", "DENV-2", "DENV-3", "DENV-4")
case_cols <- c("cases_DENV1", "cases_DENV2", "cases_DENV3", "cases_DENV4")

rt_all <- list()
for (j in 1:4) {
  cat(sprintf("  Estimating Rt for %s...\n", serotype_labels[j]))
  est <- cori_estimate(weekly_props[[case_cols[j]]], gi, S, window = 4)
  est$date <- weekly_props$date
  est$serotype <- serotype_labels[j]
  # Add the proportion for filtering noisy periods
  est$prop <- weekly_props[[prop_cols[j]]]
  rt_all[[j]] <- est
}

rt_df <- bind_rows(rt_all)

cat(sprintf("  Total Rt estimates: %d rows (%d per serotype)\n",
            nrow(rt_df), nrow(weekly_props)))

# ==============================================================================
# 4. COMPARE SEROTYPE-SPECIFIC RT
# ==============================================================================

cat("\nComparing serotype-specific Rt...\n")

# Define dominance periods based on known Singapore serotype history
# DENV-1 dominant ~2013-2015, DENV-2 dominant ~2016-2019, DENV-3 emerged 2020+
dominance_periods <- tribble(
  ~serotype,  ~period_label,       ~start,        ~end,
  "DENV-1",  "2013-2015",  "2013-01-01",  "2015-12-31",
  "DENV-2",  "2016-2019",  "2016-01-01",  "2019-12-31",
  "DENV-3",  "2020-2022",  "2020-01-01",  "2022-12-31"
) |>
  mutate(start = as.Date(start), end = as.Date(end))

# Compute summary statistics per serotype during dominant periods
rt_dominant <- rt_df |>
  filter(!is.na(Rt_median), prop > 0.05) |>
  inner_join(dominance_periods, by = "serotype") |>
  filter(date >= start, date <= end)

summary_by_period <- rt_dominant |>
  group_by(serotype, period_label) |>
  summarize(
    n_weeks = n(),
    Rt_median_of_medians = median(Rt_median, na.rm = TRUE),
    Rt_mean_of_medians = mean(Rt_median, na.rm = TRUE),
    Rt_q25 = quantile(Rt_median, 0.25, na.rm = TRUE),
    Rt_q75 = quantile(Rt_median, 0.75, na.rm = TRUE),
    Rt_sd = sd(Rt_median, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Rt summary by serotype during dominant period ---\n")
for (i in seq_len(nrow(summary_by_period))) {
  r <- summary_by_period[i, ]
  cat(sprintf("  %s (%s): median Rt = %.3f (IQR: %.3f-%.3f), n=%d weeks\n",
              r$serotype, r$period_label,
              r$Rt_median_of_medians, r$Rt_q25, r$Rt_q75, r$n_weeks))
}

# Overall summary (all weeks with prop > 5%)
summary_overall <- rt_df |>
  filter(!is.na(Rt_median), prop > 0.05) |>
  group_by(serotype) |>
  summarize(
    n_weeks = n(),
    Rt_median_of_medians = median(Rt_median, na.rm = TRUE),
    Rt_mean_of_medians = mean(Rt_median, na.rm = TRUE),
    Rt_q25 = quantile(Rt_median, 0.25, na.rm = TRUE),
    Rt_q75 = quantile(Rt_median, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Rt summary by serotype (overall, prop > 5%) ---\n")
for (i in seq_len(nrow(summary_overall))) {
  r <- summary_overall[i, ]
  cat(sprintf("  %s: median Rt = %.3f (IQR: %.3f-%.3f), n=%d weeks\n",
              r$serotype, r$Rt_median_of_medians, r$Rt_q25, r$Rt_q75, r$n_weeks))
}

# --- P(Rt_DENV3 > Rt_DENV2) from posterior overlap ---
# For each week in the DENV-3 dominant period, draw from Gamma posterior
# and compare against DENV-2 posterior at the same time

cat("\nComputing P(Rt_DENV3 > Rt_DENV2)...\n")

# Get overlapping weeks where both have sufficient data
rt_wide <- rt_df |>
  filter(!is.na(Rt_median)) |>
  select(date, serotype, cases, IP) |>
  pivot_wider(names_from = serotype, values_from = c(cases, IP),
              names_sep = "_")

# Filter to DENV-3 dominant period where both have infectious pressure
rt_compare <- rt_wide |>
  filter(date >= as.Date("2020-01-01"), date <= as.Date("2022-12-31")) |>
  filter(!is.na(`IP_DENV-3`), !is.na(`IP_DENV-2`),
         `IP_DENV-3` > 0.01)

n_sim <- 10000
set.seed(42)
prob_d3_gt_d2 <- numeric(nrow(rt_compare))

a_prior <- 1; b_prior <- 0.2

for (i in seq_len(nrow(rt_compare))) {
  # DENV-3 posterior
  a3 <- a_prior + rt_compare$`cases_DENV-3`[i]
  b3 <- b_prior + rt_compare$`IP_DENV-3`[i]
  draws3 <- rgamma(n_sim, shape = a3, rate = b3)

  # DENV-2 posterior
  a2 <- a_prior + rt_compare$`cases_DENV-2`[i]
  b2 <- b_prior + rt_compare$`IP_DENV-2`[i]
  draws2 <- rgamma(n_sim, shape = a2, rate = b2)

  prob_d3_gt_d2[i] <- mean(draws3 > draws2)
}

cat(sprintf("  P(Rt_DENV3 > Rt_DENV2) per week (2020-2022):\n"))
cat(sprintf("    Mean: %.3f\n", mean(prob_d3_gt_d2)))
cat(sprintf("    Median: %.3f\n", median(prob_d3_gt_d2)))
cat(sprintf("    Range: [%.3f, %.3f]\n", min(prob_d3_gt_d2), max(prob_d3_gt_d2)))

# ==============================================================================
# 5. FIGURES
# ==============================================================================

cat("\nGenerating figures...\n")

# --- Figure 1: Rt trajectories ---

# Only show Rt when serotype proportion > 5% (otherwise too noisy)
rt_plot_df <- rt_df |>
  filter(!is.na(Rt_median), prop > 0.05) |>
  # Cap Rt at 5 for plotting (extreme values from low-case periods)
  mutate(
    Rt_median = pmin(Rt_median, 5),
    Rt_lower  = pmin(Rt_lower, 5),
    Rt_upper  = pmin(Rt_upper, 5)
  )

# Serotype switch dates (approximate)
switch_dates <- as.Date(c("2015-07-01", "2020-01-01"))

p_rt <- ggplot(rt_plot_df, aes(x = date, color = serotype, fill = serotype)) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), alpha = 0.15,
              color = NA) +
  geom_line(aes(y = Rt_median), linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = switch_dates, linetype = "dotted",
             color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = sero_colors) +
  scale_fill_manual(values = sero_colors) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 0.5)) +
  labs(
    title = "Serotype-Specific Rt Estimates (Cori Method, 4-week Window)",
    subtitle = "Singapore 2013-2022 | shown only when serotype proportion > 5%",
    x = "Date", y = expression(R[t]),
    color = "Serotype", fill = "Serotype"
  ) +
  annotate("text", x = as.Date("2014-03-01"), y = 3.8,
           label = "DENV-1\ndominant", size = 3, color = sero_colors["DENV-1"]) +
  annotate("text", x = as.Date("2017-10-01"), y = 3.8,
           label = "DENV-2\ndominant", size = 3, color = sero_colors["DENV-2"]) +
  annotate("text", x = as.Date("2021-06-01"), y = 3.8,
           label = "DENV-3\ndominant", size = 3, color = sero_colors["DENV-3"])

ggsave("../results/figures/serotype_specific_rt_trajectories.png", p_rt,
       width = 14, height = 6, dpi = 150)
cat("  Saved: results/figures/serotype_specific_rt_trajectories.png\n")

# --- Figure 2: Boxplot of Rt by serotype during dominant periods ---

p_box <- ggplot(rt_dominant, aes(x = serotype, y = Rt_median, fill = serotype)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = sero_colors) +
  coord_cartesian(ylim = c(0, 4)) +
  labs(
    title = "Distribution of Weekly Rt by Serotype (During Dominant Periods)",
    subtitle = paste0(
      "DENV-1 (2013-2015), DENV-2 (2016-2019), DENV-3 (2020-2022) | ",
      sprintf("P(Rt_DENV3 > Rt_DENV2) = %.2f", mean(prob_d3_gt_d2))
    ),
    x = "Serotype", y = expression(R[t]),
    fill = "Serotype"
  ) +
  theme(legend.position = "none")

ggsave("../results/figures/serotype_specific_rt_boxplot.png", p_box,
       width = 7, height = 6, dpi = 150)
cat("  Saved: results/figures/serotype_specific_rt_boxplot.png\n")

# --- Figure 3: Mean Rt per serotype per dominance period (bar chart) ---

p_bar <- ggplot(summary_by_period,
                aes(x = serotype, y = Rt_mean_of_medians, fill = serotype)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = Rt_q25, ymax = Rt_q75), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_text(aes(label = period_label), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = sero_colors) +
  coord_cartesian(ylim = c(0, max(summary_by_period$Rt_q75, na.rm = TRUE) * 1.3)) +
  labs(
    title = "Mean Rt per Serotype During Dominance Period",
    subtitle = "Error bars show IQR of weekly Rt estimates",
    x = "Serotype", y = expression(paste("Mean ", R[t])),
    fill = "Serotype"
  ) +
  theme(legend.position = "none")

ggsave("../results/figures/serotype_specific_rt_by_period.png", p_bar,
       width = 7, height = 6, dpi = 150)
cat("  Saved: results/figures/serotype_specific_rt_by_period.png\n")

# --- Figure 4: Stacked area chart of serotype-specific cases ---

cases_long <- weekly_props |>
  select(date, cases_DENV1, cases_DENV2, cases_DENV3, cases_DENV4) |>
  pivot_longer(cols = starts_with("cases_"),
               names_to = "serotype", values_to = "cases") |>
  mutate(serotype = case_when(
    serotype == "cases_DENV1" ~ "DENV-1",
    serotype == "cases_DENV2" ~ "DENV-2",
    serotype == "cases_DENV3" ~ "DENV-3",
    serotype == "cases_DENV4" ~ "DENV-4"
  )) |>
  mutate(serotype = factor(serotype, levels = c("DENV-4", "DENV-3", "DENV-2", "DENV-1")))

p_cases <- ggplot(cases_long, aes(x = date, y = cases, fill = serotype)) +
  geom_area(alpha = 0.8) +
  scale_fill_manual(values = sero_colors,
                    breaks = c("DENV-1", "DENV-2", "DENV-3", "DENV-4")) +
  geom_vline(xintercept = switch_dates, linetype = "dotted",
             color = "gray30", linewidth = 0.5) +
  labs(
    title = "Estimated Serotype-Specific Weekly Dengue Cases",
    subtitle = "Singapore 2013-2022 | Total cases x GAM-smoothed serotype proportions",
    x = "Date", y = "Weekly Cases",
    fill = "Serotype"
  )

ggsave("../results/figures/serotype_specific_cases.png", p_cases,
       width = 14, height = 5, dpi = 150)
cat("  Saved: results/figures/serotype_specific_cases.png\n")

# ==============================================================================
# 6. SAVE SUMMARY STATISTICS
# ==============================================================================

cat("\nSaving summary statistics...\n")

# Combine overall and period-specific summaries
summary_out <- bind_rows(
  summary_overall |> mutate(period = "overall (prop > 5%)"),
  summary_by_period |>
    rename(period = period_label) |>
    select(serotype, period, n_weeks, Rt_median_of_medians,
           Rt_mean_of_medians, Rt_q25, Rt_q75)
) |>
  mutate(
    P_Rt_DENV3_gt_DENV2 = ifelse(
      serotype == "DENV-3" & period %in% c("2020-2022", "overall (prop > 5%)"),
      mean(prob_d3_gt_d2), NA_real_
    )
  )

write_csv(summary_out, "../results/serotype_specific_rt_summary.csv")
cat("  Saved: results/serotype_specific_rt_summary.csv\n")

# ==============================================================================
# 7. FINAL SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SEROTYPE-SPECIFIC RT ESTIMATION COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\n--- KEY FINDINGS ---\n\n")

cat("Rt by serotype during dominant period:\n")
for (i in seq_len(nrow(summary_by_period))) {
  r <- summary_by_period[i, ]
  cat(sprintf("  %s (%s): median Rt = %.3f, mean Rt = %.3f (n=%d weeks)\n",
              r$serotype, r$period_label,
              r$Rt_median_of_medians, r$Rt_mean_of_medians, r$n_weeks))
}

cat(sprintf("\nP(Rt_DENV3 > Rt_DENV2) during 2020-2022: %.3f (mean across weeks)\n",
            mean(prob_d3_gt_d2)))

cat("\nInterpretation:\n")
if (mean(prob_d3_gt_d2) > 0.75) {
  cat("  DENV-3 shows evidence of higher transmissibility than DENV-2.\n")
  cat("  This is consistent with DENV-3 entering a largely naive population\n")
  cat("  after decades of low circulation in Singapore.\n")
} else if (mean(prob_d3_gt_d2) > 0.5) {
  cat("  Weak evidence that DENV-3 Rt exceeds DENV-2 Rt.\n")
  cat("  Differences may reflect population immunity rather than\n")
  cat("  intrinsic transmissibility.\n")
} else {
  cat("  No evidence that DENV-3 transmits faster than DENV-2.\n")
}

cat("\nOutput files:\n")
cat("  results/figures/serotype_specific_rt_trajectories.png\n")
cat("  results/figures/serotype_specific_rt_boxplot.png\n")
cat("  results/figures/serotype_specific_rt_by_period.png\n")
cat("  results/figures/serotype_specific_cases.png\n")
cat("  results/serotype_specific_rt_summary.csv\n")
