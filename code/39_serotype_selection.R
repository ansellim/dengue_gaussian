#!/usr/bin/env Rscript
# ==============================================================================
# 39_serotype_selection.R
#
# Estimate the serotype selection coefficient for the DENV-2 to DENV-3
# transition in Singapore (Oct 2019 - Sep 2021).
#
# Inspired by Volz et al. (2021, Nature) who estimated B.1.1.7's transmission
# advantage using the log-ratio of variant frequencies.
#
# Three approaches:
#   Option 1: Descriptive log-ratio plot (raw + GAM-smoothed)
#   Option 2: Piecewise linear selection coefficient (pre/during/post COVID)
#   Option 3: Renewal equation decomposition with serotype-specific Rt
#
# Input:
#   data/serotype_props_2013_2016.csv
#   data/monthly_sero_type_props_all_data.csv
#   results/fit_model3.rds, data/model_data.rds
#
# Output:
#   results/figures/serotype_selection_logratio.png
#   results/figures/serotype_selection_piecewise.png
#   results/figures/serotype_selection_advantage.png
#   results/figures/serotype_selection_summary.png
#   results/serotype_selection_summary.csv
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

SEROTYPE_DIR <- file.path(dirname(getwd()), "data")
if (!file.exists(file.path(SEROTYPE_DIR, "serotype_props_2013_2016.csv"))) {
  SEROTYPE_DIR <- file.path(dirname(dirname(getwd())), "DengueFever", "data", "nea")
}
stopifnot(
  "Serotype data directory not found. Check SEROTYPE_DIR path." =
    dir.exists(SEROTYPE_DIR)
)

dir.create("../results/figures", showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank()))

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SEROTYPE SELECTION COEFFICIENT ESTIMATION\n")
cat("DENV-2 to DENV-3 transition, Singapore\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD AND COMBINE SEROTYPE DATA
# ==============================================================================

cat("Loading serotype data...\n")

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

cat(sprintf("  Combined serotype data: %d months (%s to %s)\n",
            nrow(serotype), min(serotype$year_month), max(serotype$year_month)))

# ==============================================================================
# 2. FOCUS ON TRANSITION PERIOD AND COMPUTE LOG-RATIOS
# ==============================================================================

cat("\nFocusing on DENV-2 to DENV-3 transition period (Oct 2019 - Sep 2021)...\n")

transition <- serotype |>
  filter(year_month >= as.Date("2019-10-01"),
         year_month <= as.Date("2021-09-01")) |>
  mutate(
    t_numeric = as.numeric(year_month - min(year_month)) / 30.44,  # months since start
    # Raw log-ratio with epsilon for zero-handling
    D2_safe = pmax(D2_prop, 0.001),
    D3_safe = pmax(D3_prop, 0.001),
    log_ratio_raw = log(D3_safe / D2_safe)
  )

cat(sprintf("  Transition period: %d months\n", nrow(transition)))
cat(sprintf("  D3 range: %.3f to %.3f\n", min(transition$D3_prop), max(transition$D3_prop)))
cat(sprintf("  D2 range: %.3f to %.3f\n", min(transition$D2_prop), max(transition$D2_prop)))

# GAM-smoothed log-ratio
# Smooth D2 and D3 proportions separately, then compute log-ratio
fit_d2 <- gam(pmin(pmax(D2_prop, 0.001), 0.999) ~ s(t_numeric, k = 10, bs = "tp"),
              family = quasibinomial(), data = transition)
fit_d3 <- gam(pmin(pmax(D3_prop, 0.001), 0.999) ~ s(t_numeric, k = 10, bs = "tp"),
              family = quasibinomial(), data = transition)

transition <- transition |>
  mutate(
    D2_smooth = predict(fit_d2, type = "response"),
    D3_smooth = predict(fit_d3, type = "response"),
    log_ratio_smooth = log(D3_smooth / D2_smooth)
  )

cat(sprintf("  D2 GAM deviance explained: %.1f%%\n", summary(fit_d2)$dev.expl * 100))
cat(sprintf("  D3 GAM deviance explained: %.1f%%\n", summary(fit_d3)$dev.expl * 100))

# ==============================================================================
# 3. OPTION 1: DESCRIPTIVE LOG-RATIO PLOT
# ==============================================================================

cat("\n--- OPTION 1: Descriptive Log-Ratio Plot ---\n")

# COVID phase date ranges
covid_phases <- tibble(
  phase = c("Circuit Breaker", "Phase 1", "Phase 2"),
  start = as.Date(c("2020-04-07", "2020-06-02", "2020-06-19")),
  end   = as.Date(c("2020-06-01", "2020-06-18", "2020-12-28")),
  fill  = c("gray30", "gray50", "gray75")
)

p_logratio <- ggplot(transition, aes(x = year_month)) +
  # COVID phase shading

  annotate("rect",
           xmin = covid_phases$start, xmax = covid_phases$end,
           ymin = -Inf, ymax = Inf,
           fill = covid_phases$fill, alpha = 0.3) +
  annotate("text",
           x = covid_phases$start + (covid_phases$end - covid_phases$start) / 2,
           y = max(transition$log_ratio_raw, na.rm = TRUE) * 0.95,
           label = covid_phases$phase,
           size = 2.5, angle = 90, vjust = 0, color = "gray30") +
  # Parity line

  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  # Raw data points
  geom_point(aes(y = log_ratio_raw), size = 2.5, color = "gray30", alpha = 0.7) +
  # Smoothed line
  geom_line(aes(y = log_ratio_smooth), color = "steelblue", linewidth = 1.2) +
  # Annotation
  annotate("text",
           x = as.Date("2020-01-01"),
           y = min(transition$log_ratio_raw, na.rm = TRUE) * 0.9,
           label = "If constant fitness advantage,\nthis should be linear",
           size = 3, fontface = "italic", color = "gray40", hjust = 0) +
  labs(
    title = "Serotype replacement dynamics: DENV-3 vs DENV-2",
    subtitle = "Singapore, Oct 2019 - Sep 2021",
    x = "Date",
    y = "log(DENV-3 / DENV-2)"
  )

ggsave("../results/figures/serotype_selection_logratio.png", p_logratio,
       width = 10, height = 6, dpi = 150)
cat("  Saved: results/figures/serotype_selection_logratio.png\n")

# ==============================================================================
# 4. OPTION 2: PIECEWISE LINEAR SELECTION COEFFICIENT
# ==============================================================================

cat("\n--- OPTION 2: Piecewise Linear Selection Coefficient ---\n")

# Define periods
transition <- transition |>
  mutate(
    period = case_when(
      year_month >= as.Date("2019-10-01") & year_month <= as.Date("2020-03-01") ~ "Pre-COVID",
      year_month >= as.Date("2020-04-01") & year_month <= as.Date("2020-12-01") ~ "During COVID",
      year_month >= as.Date("2021-01-01") & year_month <= as.Date("2021-09-01") ~ "Post-COVID",
      TRUE ~ NA_character_
    ),
    period = factor(period, levels = c("Pre-COVID", "During COVID", "Post-COVID"))
  )

# Fit piecewise linear regressions
period_results <- list()

for (p in levels(transition$period)) {
  sub_df <- transition |> filter(period == p)
  if (nrow(sub_df) < 3) {
    cat(sprintf("  Skipping %s: only %d data points\n", p, nrow(sub_df)))
    next
  }

  # Time in months relative to start of this period
  sub_df <- sub_df |>
    mutate(t_period = as.numeric(year_month - min(year_month)) / 30.44)

  lm_fit <- lm(log_ratio_raw ~ t_period, data = sub_df)
  s <- summary(lm_fit)
  ci <- confint(lm_fit, "t_period", level = 0.95)

  period_results[[p]] <- tibble(
    period = p,
    slope = coef(lm_fit)["t_period"],
    slope_lower = ci[1],
    slope_upper = ci[2],
    intercept = coef(lm_fit)["(Intercept)"],
    r_squared = s$r.squared,
    p_value = coef(s)["t_period", "Pr(>|t|)"],
    n_months = nrow(sub_df)
  )

  cat(sprintf("  %s: s = %.3f (95%% CI: %.3f to %.3f), R2 = %.3f, p = %.4f, n = %d\n",
              p, coef(lm_fit)["t_period"], ci[1], ci[2],
              s$r.squared, coef(s)["t_period", "Pr(>|t|)"], nrow(sub_df)))
}

piecewise_summary <- bind_rows(period_results)

# Overall linear model (ignoring COVID disruption)
lm_overall <- lm(log_ratio_raw ~ t_numeric, data = transition)
s_overall <- summary(lm_overall)
ci_overall <- confint(lm_overall, "t_numeric", level = 0.95)

cat(sprintf("\n  Overall (ignoring COVID): s = %.3f (95%% CI: %.3f to %.3f), R2 = %.3f, p = %.4f\n",
            coef(lm_overall)["t_numeric"], ci_overall[1], ci_overall[2],
            s_overall$r.squared, coef(s_overall)["t_numeric", "Pr(>|t|)"]))

# Generate predictions for each period's fitted line
pred_lines <- list()
for (p in levels(transition$period)) {
  sub_df <- transition |> filter(period == p)
  if (nrow(sub_df) < 3) next
  sub_df <- sub_df |>
    mutate(t_period = as.numeric(year_month - min(year_month)) / 30.44)
  lm_fit <- lm(log_ratio_raw ~ t_period, data = sub_df)
  sub_df$fitted <- predict(lm_fit)
  pred_lines[[p]] <- sub_df |> select(year_month, period, fitted)
}
pred_df <- bind_rows(pred_lines)

# Period colors
period_colors <- c("Pre-COVID" = "#2166AC", "During COVID" = "#B2182B", "Post-COVID" = "#1B7837")

# Build annotation labels
anno_labels <- piecewise_summary |>
  mutate(label = sprintf("s = %.3f [%.3f, %.3f]", slope, slope_lower, slope_upper))

# Compute label positions: midpoint of each period's date range, with the fitted value
anno_pos <- pred_df |>
  group_by(period) |>
  summarize(
    x = median(year_month),
    y = median(fitted),
    .groups = "drop"
  ) |>
  left_join(anno_labels |> select(period, label), by = "period")

p_piecewise <- ggplot(transition |> filter(!is.na(period)), aes(x = year_month)) +
  # COVID phase shading

  annotate("rect",
           xmin = covid_phases$start, xmax = covid_phases$end,
           ymin = -Inf, ymax = Inf,
           fill = covid_phases$fill, alpha = 0.3) +
  # Parity line
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  # Raw data points colored by period
  geom_point(aes(y = log_ratio_raw, color = period), size = 2.5, alpha = 0.7) +
  # Fitted regression lines
  geom_line(data = pred_df, aes(y = fitted, color = period),
            linewidth = 1.2) +
  # Slope annotations
  geom_label(data = anno_pos,
             aes(x = x, y = y + 0.4, label = label, color = period),
             size = 3, fill = "white", alpha = 0.8, label.size = 0.3,
             show.legend = FALSE) +
  scale_color_manual(values = period_colors, name = "Period") +
  labs(
    title = "Piecewise selection coefficients: DENV-3 vs DENV-2",
    subtitle = "Slope = selection coefficient s (per month); s > 0 means DENV-3 growing",
    x = "Date",
    y = "log(DENV-3 / DENV-2)"
  ) +
  theme(legend.position = "bottom")

ggsave("../results/figures/serotype_selection_piecewise.png", p_piecewise,
       width = 10, height = 6, dpi = 150)
cat("  Saved: results/figures/serotype_selection_piecewise.png\n")

# ==============================================================================
# 4b. BOOTSTRAP CIs AND PERMUTATION TESTS FOR OPTION 2
# ==============================================================================

cat("\n--- Bootstrap CIs and Permutation Tests for Piecewise Slopes ---\n")

set.seed(42)
n_boot <- 10000
n_perm <- 10000

boot_results <- list()
perm_results <- list()

for (p in levels(transition$period)) {
  sub_df <- transition |> filter(period == p)
  if (nrow(sub_df) < 3) next

  sub_df <- sub_df |>
    mutate(t_period = as.numeric(year_month - min(year_month)) / 30.44)

  n_obs <- nrow(sub_df)
  observed_slope <- coef(lm(log_ratio_raw ~ t_period, data = sub_df))["t_period"]

  # --- Bootstrap ---
  boot_slopes <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(n_obs, n_obs, replace = TRUE)
    boot_df <- sub_df[idx, ]
    boot_slopes[b] <- coef(lm(log_ratio_raw ~ t_period, data = boot_df))["t_period"]
  }
  boot_ci <- quantile(boot_slopes, c(0.025, 0.975))

  boot_results[[p]] <- tibble(
    period = p,
    boot_lower = boot_ci[1],
    boot_upper = boot_ci[2]
  )

  # --- Permutation test ---
  perm_slopes <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    perm_df <- sub_df
    perm_df$log_ratio_raw <- sample(perm_df$log_ratio_raw)
    perm_slopes[b] <- coef(lm(log_ratio_raw ~ t_period, data = perm_df))["t_period"]
  }

  # Two-sided p-value
  perm_p_two <- mean(abs(perm_slopes) >= abs(observed_slope))
  # One-sided p-value (slope >= observed)
  perm_p_one <- mean(perm_slopes >= observed_slope)

  perm_results[[p]] <- tibble(
    period = p,
    perm_p_twosided = perm_p_two,
    perm_p_onesided = perm_p_one
  )

  cat(sprintf("  %s: bootstrap CI [%.3f, %.3f], permutation p (two-sided) = %.4f\n",
              p, boot_ci[1], boot_ci[2], perm_p_two))
}

boot_summary <- bind_rows(boot_results)
perm_summary <- bind_rows(perm_results)

# Merge bootstrap and permutation results into piecewise_summary
piecewise_summary <- piecewise_summary |>
  left_join(boot_summary, by = "period") |>
  left_join(perm_summary, by = "period")

# Print formatted summary
cat("\nSelection coefficient estimates with uncertainty:\n")
for (i in seq_len(nrow(piecewise_summary))) {
  row <- piecewise_summary[i, ]
  cat(sprintf("  %-14s s = %.3f, parametric CI [%.3f, %.3f], bootstrap CI [%.3f, %.3f], permutation p = %.4f\n",
              paste0(row$period, ":"),
              row$slope, row$slope_lower, row$slope_upper,
              row$boot_lower, row$boot_upper,
              row$perm_p_twosided))
}

cat("\n  Comparison: If bootstrap CIs agree with parametric CIs, the result is robust.\n")
cat("  Permutation p-values test whether the trend could arise by chance.\n")

# --- Update piecewise figure with bootstrap CIs ---

# Recompute annotation labels with both CI types
anno_labels_updated <- piecewise_summary |>
  mutate(label = sprintf("s = %.3f\nParametric [%.3f, %.3f]\nBootstrap [%.3f, %.3f]\nPerm p = %.4f",
                          slope, slope_lower, slope_upper,
                          boot_lower, boot_upper, perm_p_twosided))

anno_pos_updated <- pred_df |>
  group_by(period) |>
  summarize(
    x = median(year_month),
    y = median(fitted),
    .groups = "drop"
  ) |>
  left_join(anno_labels_updated |> select(period, label, boot_lower, boot_upper, slope),
            by = "period")

# Build error bar data: one point per period at the midpoint date
# Show both parametric and bootstrap CIs as error bars
ci_errorbar_data <- piecewise_summary |>
  left_join(
    pred_df |>
      group_by(period) |>
      summarize(x = median(year_month), y_fitted = median(fitted), .groups = "drop"),
    by = "period"
  )

p_piecewise_boot <- ggplot(transition |> filter(!is.na(period)), aes(x = year_month)) +
  # COVID phase shading
  annotate("rect",
           xmin = covid_phases$start, xmax = covid_phases$end,
           ymin = -Inf, ymax = Inf,
           fill = covid_phases$fill, alpha = 0.3) +
  # Parity line
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  # Raw data points colored by period
  geom_point(aes(y = log_ratio_raw, color = period), size = 2.5, alpha = 0.7) +
  # Fitted regression lines
  geom_line(data = pred_df, aes(y = fitted, color = period),
            linewidth = 1.2) +
  # Slope annotations with bootstrap CIs
  geom_label(data = anno_pos_updated,
             aes(x = x, y = y + 0.6, label = label, color = period),
             size = 2.5, fill = "white", alpha = 0.8, label.size = 0.3,
             show.legend = FALSE, lineheight = 0.9) +
  scale_color_manual(values = period_colors, name = "Period") +
  labs(
    title = "Piecewise selection coefficients with bootstrap CIs and permutation tests",
    subtitle = "Slope = selection coefficient s (per month); parametric + bootstrap CIs; permutation p-values",
    x = "Date",
    y = "log(DENV-3 / DENV-2)"
  ) +
  theme(legend.position = "bottom")

ggsave("../results/figures/serotype_selection_piecewise.png", p_piecewise_boot,
       width = 11, height = 7, dpi = 150)
cat("  Updated: results/figures/serotype_selection_piecewise.png (with bootstrap CIs)\n")

# Update p_piecewise reference for the summary figure
p_piecewise <- p_piecewise_boot

# ==============================================================================
# 5. OPTION 3: RENEWAL EQUATION WITH SEROTYPE-SPECIFIC Rt
# ==============================================================================

cat("\n--- OPTION 3: Renewal Equation Decomposition ---\n")

# Load model fit and extract Rt posterior
if (file.exists("../results/fit_model3.rds")) {
  fit <- readRDS("../results/fit_model3.rds")
  cat("  Using Model 3 (climate-only) fit\n")
} else if (file.exists("../results/fit_model2.rds")) {
  fit <- readRDS("../results/fit_model2.rds")
  cat("  Using Model 2 (full) fit as fallback\n")
} else {
  stop("No model fit found. Run 13_fit_model3.R first.")
}

model_data <- readRDS("../data/model_data.rds")
dates_model <- model_data$metadata$dates_model

# Extract Rt draws
Rt_draws <- fit$draws("Rt", format = "matrix")
n_draws <- nrow(Rt_draws)
cat(sprintf("  %d posterior draws, %d modeled weeks\n", n_draws, ncol(Rt_draws)))

# Map model weeks to calendar months, restrict to transition period
week_df <- tibble(
  week_idx = seq_along(dates_model),
  date = dates_model,
  year_month = floor_date(dates_model, "month")
) |>
  filter(year_month >= as.Date("2019-10-01"),
         year_month <= as.Date("2021-09-01"))

# Average Rt within each month
month_groups <- week_df |>
  group_by(year_month) |>
  summarize(week_indices = list(week_idx), .groups = "drop")

Rt_monthly <- matrix(NA_real_, nrow = n_draws, ncol = nrow(month_groups))
for (m in seq_len(nrow(month_groups))) {
  idx <- month_groups$week_indices[[m]]
  Rt_monthly[, m] <- rowMeans(Rt_draws[, idx, drop = FALSE])
}

# Join with serotype proportions
trans_rt <- month_groups |>
  select(year_month) |>
  inner_join(
    transition |> select(year_month, D2_prop, D3_prop, D2_smooth, D3_smooth),
    by = "year_month"
  )

if (nrow(trans_rt) == 0) {
  cat("  WARNING: No overlapping months between Rt and serotype data.\n")
  cat("  Skipping Option 3.\n")
  p_advantage <- ggplot() + annotate("text", x = 0.5, y = 0.5,
    label = "Option 3 not available:\nNo overlap between Rt and serotype data") +
    theme_void()
  delta_summary <- tibble(
    delta_median = NA_real_, delta_lower = NA_real_, delta_upper = NA_real_
  )
} else {
  n_overlap <- nrow(trans_rt)
  cat(sprintf("  Overlapping months for decomposition: %d\n", n_overlap))

  # Align Rt_monthly columns to overlapping months
  overlap_idx <- which(month_groups$year_month %in% trans_rt$year_month)
  Rt_overlap <- Rt_monthly[, overlap_idx, drop = FALSE]

  # Use smoothed proportions for stability; renormalize D2+D3 to sum to 1
  # (ignore D1, D4 which are negligible in this period)
  prop_D2 <- trans_rt$D2_smooth / (trans_rt$D2_smooth + trans_rt$D3_smooth)
  prop_D3 <- trans_rt$D3_smooth / (trans_rt$D2_smooth + trans_rt$D3_smooth)

  # Grid search over delta: Rt_D3 = (1 + delta) * Rt_D2
  # Rt_agg = Rt_D2 * (prop_D2 + (1+delta) * prop_D3)
  # => Rt_D2 = Rt_agg / (prop_D2 + (1+delta) * prop_D3)
  #
  # To find optimal delta, we use the log-ratio approach:
  # From the proportions, the expected change in log(D3/D2) over one generation

  # is log(1+delta). We fit delta to match the observed log-ratio trajectory.
  #
  # For each posterior draw of Rt_agg, find delta that minimizes sum of squared
  # differences between implied and observed serotype proportions.
  #
  # Simplified approach: for each draw, compute delta from the relationship
  # between Rt and proportion changes.

  delta_grid <- seq(0, 1, by = 0.01)
  n_delta <- length(delta_grid)

  # For each draw, find delta that minimizes discrepancy
  # Use a subsample for computational tractability
  n_sub <- min(n_draws, 500)
  set.seed(42)
  draw_idx <- sample(n_draws, n_sub)

  delta_posterior <- numeric(n_sub)

  cat("  Estimating delta via grid search over posterior draws...\n")

  for (i in seq_len(n_sub)) {
    rt_agg <- Rt_overlap[draw_idx[i], ]

    # For each delta, compute implied Rt_D2 and Rt_D3, then compute
    # implied next-month proportion shift and compare to observed
    sse <- numeric(n_delta)

    for (d in seq_len(n_delta)) {
      delta <- delta_grid[d]
      denom <- prop_D2 + (1 + delta) * prop_D3
      Rt_D2 <- rt_agg / denom
      Rt_D3 <- (1 + delta) * Rt_D2

      # Implied proportion change: next month's D3 proportion
      # D3_next_prop ~ D3_prop * Rt_D3 / (D2_prop * Rt_D2 + D3_prop * Rt_D3)
      implied_D3 <- (prop_D3 * Rt_D3) / (prop_D2 * Rt_D2 + prop_D3 * Rt_D3)

      # Compare with actual next-month proportion (shifted by 1)
      if (n_overlap > 1) {
        observed_next <- prop_D3[2:n_overlap]
        implied_next  <- implied_D3[1:(n_overlap - 1)]
        sse[d] <- sum((observed_next - implied_next)^2)
      } else {
        sse[d] <- Inf
      }
    }

    delta_posterior[i] <- delta_grid[which.min(sse)]
  }

  delta_summary <- tibble(
    delta_median = median(delta_posterior),
    delta_lower  = quantile(delta_posterior, 0.025),
    delta_upper  = quantile(delta_posterior, 0.975)
  )

  cat(sprintf("  Estimated delta: %.3f (95%% CrI: %.3f to %.3f)\n",
              delta_summary$delta_median, delta_summary$delta_lower, delta_summary$delta_upper))
  cat(sprintf("  Interpretation: DENV-3 Rt is ~%.0f%% (%.0f%% to %.0f%%) higher than DENV-2\n",
              delta_summary$delta_median * 100, delta_summary$delta_lower * 100,
              delta_summary$delta_upper * 100))

  # Compute Rt_D2 and Rt_D3 at MAP delta for plotting
  delta_map <- delta_summary$delta_median
  # Use posterior median Rt for the trajectory plot
  rt_agg_median <- apply(Rt_overlap, 2, median)
  denom_map <- prop_D2 + (1 + delta_map) * prop_D3
  Rt_D2_map <- rt_agg_median / denom_map
  Rt_D3_map <- (1 + delta_map) * Rt_D2_map

  rt_sero_df <- tibble(
    year_month = trans_rt$year_month,
    Rt_aggregate = rt_agg_median,
    Rt_DENV2 = Rt_D2_map,
    Rt_DENV3 = Rt_D3_map
  ) |>
    pivot_longer(cols = starts_with("Rt_"), names_to = "series", values_to = "Rt") |>
    mutate(series = str_replace(series, "Rt_", "") |>
             str_replace("DENV", "DENV-") |>
             str_replace("aggregate", "Aggregate"))

  # Top panel: Rt trajectories
  p_rt_sero <- ggplot(rt_sero_df, aes(x = year_month, y = Rt, color = series)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1) +
    scale_color_manual(
      values = c("Aggregate" = "gray40", "DENV-2" = "#D95F02", "DENV-3" = "#7570B3"),
      name = ""
    ) +
    labs(
      title = sprintf("Serotype-specific Rt (delta = %.3f)", delta_map),
      x = "Date",
      y = expression(R[t])
    ) +
    theme(legend.position = "bottom")

  # Bottom panel: posterior of delta
  p_delta_post <- ggplot(tibble(delta = delta_posterior), aes(x = delta)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30,
                   fill = "steelblue", color = "white", alpha = 0.7) +
    geom_vline(xintercept = delta_summary$delta_median, color = "red",
               linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "gray50", linetype = "dotted") +
    annotate("text", x = delta_summary$delta_median + 0.02,
             y = Inf, vjust = 1.5,
             label = sprintf("Median = %.3f\n95%% CrI: [%.3f, %.3f]",
                             delta_summary$delta_median,
                             delta_summary$delta_lower,
                             delta_summary$delta_upper),
             size = 3.5, hjust = 0) +
    annotate("text", x = 0.02, y = Inf, vjust = 3,
             label = "delta = 0: no advantage",
             size = 2.8, color = "gray50", hjust = 0, fontface = "italic") +
    labs(
      title = "Posterior distribution of relative fitness advantage (delta)",
      subtitle = "Rt(DENV-3) = (1 + delta) * Rt(DENV-2)",
      x = expression(delta),
      y = "Density"
    )

  p_advantage <- p_rt_sero / p_delta_post +
    plot_layout(heights = c(1, 1))

  ggsave("../results/figures/serotype_selection_advantage.png", p_advantage,
         width = 10, height = 9, dpi = 150)
  cat("  Saved: results/figures/serotype_selection_advantage.png\n")
}

# ==============================================================================
# 6. SUMMARY FIGURE (3-panel patchwork)
# ==============================================================================

cat("\n--- Generating Summary Figure ---\n")

# Recreate panel A and B with compact titles for patchwork
pA <- p_logratio + labs(title = "A) Log-ratio with COVID phases")
pB <- p_piecewise + labs(title = "B) Piecewise selection coefficients")

if (inherits(p_advantage, "patchwork")) {
  # Extract just the Rt panel for the summary
  pC <- p_advantage[[1]] + labs(title = "C) Serotype-specific Rt decomposition")
  p_summary <- (pA / pB / pC) +
    plot_annotation(
      title = "Serotype selection coefficient: DENV-2 to DENV-3 transition",
      subtitle = "Singapore, Oct 2019 - Sep 2021",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
} else {
  pC <- p_advantage + labs(title = "C) Serotype-specific Rt decomposition")
  p_summary <- (pA / pB / pC) +
    plot_annotation(
      title = "Serotype selection coefficient: DENV-2 to DENV-3 transition",
      subtitle = "Singapore, Oct 2019 - Sep 2021",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
}

ggsave("../results/figures/serotype_selection_summary.png", p_summary,
       width = 11, height = 16, dpi = 150)
cat("  Saved: results/figures/serotype_selection_summary.png\n")

# ==============================================================================
# 7. SAVE SUMMARY TABLE
# ==============================================================================

cat("\n--- Saving summary CSV ---\n")

# Combine all results
overall_row <- tibble(
  period = "Overall (ignoring COVID)",
  slope = coef(lm_overall)["t_numeric"],
  slope_lower = ci_overall[1],
  slope_upper = ci_overall[2],
  intercept = coef(lm_overall)["(Intercept)"],
  r_squared = s_overall$r.squared,
  p_value = coef(s_overall)["t_numeric", "Pr(>|t|)"],
  n_months = nrow(transition),
  boot_lower = NA_real_,
  boot_upper = NA_real_,
  perm_p_twosided = NA_real_,
  perm_p_onesided = NA_real_
)

selection_summary <- bind_rows(piecewise_summary, overall_row) |>
  mutate(
    delta_median = delta_summary$delta_median,
    delta_lower_95 = delta_summary$delta_lower,
    delta_upper_95 = delta_summary$delta_upper
  )

write_csv(selection_summary, "../results/serotype_selection_summary.csv")
cat("  Saved: results/serotype_selection_summary.csv\n")

# ==============================================================================
# 8. INTERPRETATION
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("INTERPRETATION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("1. Linearity of the log-ratio:\n")
cat("   The log(DENV-3/DENV-2) trajectory is clearly NON-LINEAR. A constant\n")
cat("   fitness advantage would produce a straight line (Volz et al. 2021),\n")
cat("   but the observed pattern shows two distinct waves separated by a\n")
cat("   reversal during COVID restrictions. This indicates the replacement\n")
cat("   dynamics are CONTEXT-DEPENDENT, not driven by intrinsic fitness alone.\n\n")

cat("2. COVID disruption:\n")
cat(sprintf("   Pre-COVID selection coefficient:    s = %.3f per month\n",
            piecewise_summary$slope[piecewise_summary$period == "Pre-COVID"]))
cat(sprintf("   During COVID selection coefficient: s = %.3f per month\n",
            piecewise_summary$slope[piecewise_summary$period == "During COVID"]))
cat(sprintf("   Post-COVID selection coefficient:   s = %.3f per month\n",
            piecewise_summary$slope[piecewise_summary$period == "Post-COVID"]))
cat("   COVID measures disrupted the replacement process. The reversal during\n")
cat("   the Circuit Breaker suggests that transmission reduction differentially\n")
cat("   affected DENV-3 (the rising serotype), likely because it depended on\n")
cat("   sustained human mobility to access the large susceptible pool.\n\n")

cat("3. Fitness advantage (delta):\n")
if (!is.na(delta_summary$delta_median)) {
  cat(sprintf("   Estimated delta = %.3f (95%% CrI: %.3f to %.3f)\n",
              delta_summary$delta_median, delta_summary$delta_lower, delta_summary$delta_upper))
  if (delta_summary$delta_lower > 0) {
    cat("   Delta is significantly different from 0.\n\n")
  } else {
    cat("   The 95% CrI includes 0, so the evidence is not conclusive.\n\n")
  }
} else {
  cat("   Could not estimate delta (no overlapping Rt/serotype data).\n\n")
}

cat("4. CAVEAT:\n")
cat("   The estimated selection coefficient conflates intrinsic viral fitness\n")
cat("   with population-level immunity effects. DENV-3's apparent advantage\n")
cat("   likely reflects the large susceptible pool (DENV-3 had not circulated\n")
cat("   widely since ~2008) rather than intrinsic transmissibility. This is\n")
cat("   an ecological fitness advantage, not a virological one.\n\n")

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SEROTYPE SELECTION ANALYSIS COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("\nOutput files:\n")
cat("  results/figures/serotype_selection_logratio.png   - Option 1: log-ratio plot\n")
cat("  results/figures/serotype_selection_piecewise.png  - Option 2: piecewise regression\n")
cat("  results/figures/serotype_selection_advantage.png  - Option 3: Rt decomposition\n")
cat("  results/figures/serotype_selection_summary.png    - 3-panel summary\n")
cat("  results/serotype_selection_summary.csv            - Summary table\n")
