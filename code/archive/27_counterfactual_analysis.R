#!/usr/bin/env Rscript
# ==============================================================================
# 27_counterfactual_analysis.R
#
# Counterfactual analysis: excess dengue cases attributable to non-climate
# drivers, estimated by comparing the full model Rt to a climate-only
# counterfactual Rt (with the residual GP zeroed out).
#
# log(Rt_actual) = mu + f_climate + f_residual
# log(Rt_cf)     = mu + f_climate   (counterfactual: no residual)
#
# Input:  results/fit_model3.rds, data/model_data.rds,
#         results/serotype_switch_timing.csv
# Output: results/figures/counterfactual_rt.png
#         results/figures/counterfactual_excess_cases.png
#         results/figures/counterfactual_switch_attribution.png
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
cat("COUNTERFACTUAL ANALYSIS - EXCESS CASES FROM NON-CLIMATE DRIVERS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA AND FIT
# ==============================================================================

cat("Loading data and model fit...\n")

model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data
dates <- model_data$metadata$dates_model
N_model <- length(dates)

S <- stan_data$S
N <- stan_data$N
cases_all <- stan_data$cases  # full case vector (length N)
gi <- stan_data$gi            # generation interval PMF (length S)

fit <- readRDS("../results/fit_model3.rds")
cat(sprintf("  %d modeled weeks, S = %d, N = %d\n", N_model, S, N))
cat(sprintf("  Date range: %s to %s\n", min(dates), max(dates)))

# ==============================================================================
# 2. EXTRACT POSTERIOR DRAWS (thinned to 500)
# ==============================================================================

cat("\nExtracting posterior draws...\n")

mu_draws <- fit$draws("mu", format = "matrix")            # n_draws x 1
f_climate_draws <- fit$draws("f_climate", format = "matrix")   # n_draws x N_model
f_residual_draws <- fit$draws("f_residual", format = "matrix") # n_draws x N_model
Rt_draws <- fit$draws("Rt", format = "matrix")             # n_draws x N_model
phi_draws <- fit$draws("phi", format = "matrix")           # n_draws x 1

n_draws_total <- nrow(mu_draws)
cat(sprintf("  Total draws: %d\n", n_draws_total))

# Thin to 500 for tractability
set.seed(42)
n_thin <- min(500, n_draws_total)
thin_idx <- sort(sample(n_draws_total, n_thin))

mu_draws <- mu_draws[thin_idx, , drop = FALSE]
f_climate_draws <- f_climate_draws[thin_idx, , drop = FALSE]
f_residual_draws <- f_residual_draws[thin_idx, , drop = FALSE]
Rt_draws <- Rt_draws[thin_idx, , drop = FALSE]
phi_draws <- phi_draws[thin_idx, , drop = FALSE]

cat(sprintf("  Thinned to %d draws\n", n_thin))

# Also extract predicted cases for the observed fit
cases_pred_draws <- fit$draws("cases_pred", format = "matrix")[thin_idx, , drop = FALSE]

# ==============================================================================
# 3. COMPUTE COUNTERFACTUAL Rt
# ==============================================================================

cat("\nComputing counterfactual (climate-only) Rt...\n")

# Counterfactual: log(Rt_cf) = mu + f_climate (no residual)
log_Rt_cf_draws <- sweep(f_climate_draws, 1, mu_draws[, 1], "+")
Rt_cf_draws <- exp(log_Rt_cf_draws)

cat("  Done.\n")

# ==============================================================================
# 4. COMPUTE INFECTIOUS PRESSURE AND EXPECTED CASES
# ==============================================================================

cat("\nReconstructing infectious pressure and expected cases...\n")

# Infectious pressure: ip[t] = sum_{s=1}^{S} cases[S+t-s] * gi[s]
# This uses observed cases, identical for actual and counterfactual
ip <- numeric(N_model)
for (i in 1:N_model) {
  t_idx <- S + i
  for (s in 1:S) {
    ip[i] <- ip[i] + cases_all[t_idx - s] * gi[s]
  }
}

cat(sprintf("  Infectious pressure range: %.1f - %.1f\n", min(ip), max(ip)))

# Actual expected cases: lambda[t] = Rt[t] * ip[t], floored at 1
# Counterfactual expected cases: lambda_cf[t] = Rt_cf[t] * ip[t], floored at 1
lambda_draws <- sweep(Rt_draws, 2, ip, "*")
lambda_draws <- pmax(lambda_draws, 1.0)

lambda_cf_draws <- sweep(Rt_cf_draws, 2, ip, "*")
lambda_cf_draws <- pmax(lambda_cf_draws, 1.0)

# Excess expected cases per week
excess_draws <- lambda_draws - lambda_cf_draws

cat("  Expected cases and excess computed for all draws.\n")

# ==============================================================================
# 5. SUMMARIZE
# ==============================================================================

cat("\nSummarizing posterior quantities...\n")

# Observed cases (model period)
cases_obs <- cases_all[(S + 1):N]

# Predicted cases summary
pred_summary <- tibble(
  date = dates,
  observed = cases_obs,
  pred_median = apply(cases_pred_draws, 2, median),
  pred_lower = apply(cases_pred_draws, 2, quantile, 0.025),
  pred_upper = apply(cases_pred_draws, 2, quantile, 0.975)
)

# Rt summary (actual and counterfactual)
rt_summary <- tibble(
  date = dates,
  Rt_median = apply(Rt_draws, 2, median),
  Rt_lower = apply(Rt_draws, 2, quantile, 0.025),
  Rt_upper = apply(Rt_draws, 2, quantile, 0.975),
  Rt_cf_median = apply(Rt_cf_draws, 2, median),
  Rt_cf_lower = apply(Rt_cf_draws, 2, quantile, 0.025),
  Rt_cf_upper = apply(Rt_cf_draws, 2, quantile, 0.975)
)

# Excess cases summary
excess_summary <- tibble(
  date = dates,
  excess_median = apply(excess_draws, 2, median),
  excess_lower = apply(excess_draws, 2, quantile, 0.025),
  excess_upper = apply(excess_draws, 2, quantile, 0.975),
  excess_mean = colMeans(excess_draws)
)

# Cumulative excess (per draw, then summarize)
cum_excess_draws <- t(apply(excess_draws, 1, cumsum))
excess_summary <- excess_summary |>
  mutate(
    cum_median = apply(cum_excess_draws, 2, median),
    cum_lower = apply(cum_excess_draws, 2, quantile, 0.025),
    cum_upper = apply(cum_excess_draws, 2, quantile, 0.975)
  )

# Annual excess
excess_summary <- excess_summary |>
  mutate(year = year(date))

annual_excess <- excess_summary |>
  group_by(year) |>
  summarize(
    annual_excess = sum(excess_mean),
    .groups = "drop"
  )

# ==============================================================================
# 6. LOAD SEROTYPE SWITCHES
# ==============================================================================

cat("\nLoading serotype switch data...\n")

serotype_switches <- NULL
if (file.exists("../results/serotype_switch_timing.csv")) {
  serotype_switches <- read_csv("../results/serotype_switch_timing.csv",
                                 show_col_types = FALSE)
  cat(sprintf("  Loaded %d serotype switch events\n", nrow(serotype_switches)))
} else {
  cat("  Warning: serotype_switch_timing.csv not found. Skipping switch attribution.\n")
}

# Switch dates for vertical lines
switch_vlines <- NULL
if (!is.null(serotype_switches) && "switch_date" %in% colnames(serotype_switches)) {
  switch_vlines <- tibble(
    date = as.Date(serotype_switches$switch_date),
    label = serotype_switches$switch_event
  ) |>
    filter(date >= min(dates), date <= max(dates))
}

# ==============================================================================
# 7. FIGURE 1: COUNTERFACTUAL Rt (2-panel)
# ==============================================================================

cat("\nGenerating counterfactual Rt plot...\n")

# Top panel: observed cases with model predicted
p_cases <- ggplot(pred_summary, aes(x = date)) +
  geom_bar(aes(y = observed), stat = "identity", fill = "gray75", alpha = 0.7) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
              fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = pred_median), color = "steelblue", linewidth = 0.7) +
  labs(
    title = "Observed vs Model-Predicted Weekly Dengue Cases",
    x = NULL,
    y = "Weekly cases"
  )

if (!is.null(switch_vlines) && nrow(switch_vlines) > 0) {
  p_cases <- p_cases +
    geom_vline(data = switch_vlines, aes(xintercept = date),
               linetype = "dotted", color = "red", linewidth = 0.5)
}

# Bottom panel: actual vs counterfactual Rt with shaded gap
# Prepare data for gap shading
gap_df <- rt_summary |>
  mutate(
    fill_color = ifelse(Rt_median > Rt_cf_median, "Excess (Rt > Rt_cf)", "Deficit (Rt < Rt_cf)")
  )

p_rt <- ggplot(rt_summary, aes(x = date)) +
  # Shade the gap between actual and counterfactual
  geom_ribbon(aes(ymin = pmin(Rt_median, Rt_cf_median),
                  ymax = pmax(Rt_median, Rt_cf_median),
                  fill = ifelse(Rt_median > Rt_cf_median, "above", "below")),
              alpha = 0.3, show.legend = TRUE) +
  # Counterfactual Rt (climate-only)
  geom_ribbon(aes(ymin = Rt_cf_lower, ymax = Rt_cf_upper),
              fill = "darkorange", alpha = 0.15) +
  geom_line(aes(y = Rt_cf_median, color = "Climate-only (counterfactual)"),
            linewidth = 0.7) +
  # Actual Rt
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper),
              fill = "steelblue", alpha = 0.15) +
  geom_line(aes(y = Rt_median, color = "Actual (full model)"),
            linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_fill_manual(
    values = c("above" = "red", "below" = "blue"),
    labels = c("above" = "Rt > counterfactual", "below" = "Rt < counterfactual"),
    name = "Gap"
  ) +
  scale_color_manual(
    values = c("Actual (full model)" = "steelblue",
               "Climate-only (counterfactual)" = "darkorange"),
    name = ""
  ) +
  labs(
    title = "Actual vs Counterfactual Rt",
    subtitle = "Red shading: non-climate drivers increase Rt; Blue: decrease Rt",
    x = "Date",
    y = expression(R[t])
  ) +
  theme(legend.position = "bottom")

if (!is.null(switch_vlines) && nrow(switch_vlines) > 0) {
  p_rt <- p_rt +
    geom_vline(data = switch_vlines, aes(xintercept = date),
               linetype = "dotted", color = "red", linewidth = 0.5)
}

p_cf_rt <- p_cases / p_rt +
  plot_annotation(
    title = "Counterfactual Analysis: Climate-Only vs Full Model",
    subtitle = "Singapore Dengue, 2012-2022",
    tag_levels = "A"
  )

ggsave("../results/figures/counterfactual_rt.png", p_cf_rt,
       width = 13, height = 9, dpi = 150)
cat("  Saved: results/figures/counterfactual_rt.png\n")

# ==============================================================================
# 8. FIGURE 2: EXCESS CASES TIME SERIES
# ==============================================================================

cat("\nGenerating excess cases plot...\n")

# Color by sign of excess
excess_pos <- excess_summary |> mutate(
  excess_pos = pmax(excess_median, 0),
  excess_neg = pmin(excess_median, 0)
)

p_excess <- ggplot(excess_summary, aes(x = date)) +
  # Weekly excess ribbon
  geom_ribbon(aes(ymin = excess_lower, ymax = excess_upper),
              fill = "gray70", alpha = 0.3) +
  # Color excess by sign
  geom_bar(data = excess_pos |> filter(excess_median >= 0),
           aes(x = date, y = excess_median),
           stat = "identity", fill = "red", alpha = 0.5, width = 7) +
  geom_bar(data = excess_pos |> filter(excess_median < 0),
           aes(x = date, y = excess_median),
           stat = "identity", fill = "blue", alpha = 0.5, width = 7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray30") +
  # Cumulative excess as secondary line
  geom_line(aes(y = cum_median / max(abs(cum_median)) *
                  max(abs(c(excess_summary$excess_lower, excess_summary$excess_upper))) * 0.5),
            color = "black", linewidth = 0.8, linetype = "solid") +
  labs(
    title = "Weekly Excess Cases Attributable to Non-Climate Drivers",
    subtitle = "Red = excess (non-climate increases transmission); Blue = deficit; Black line = cumulative (scaled)",
    x = "Date",
    y = "Excess cases / week (median)"
  )

if (!is.null(switch_vlines) && nrow(switch_vlines) > 0) {
  p_excess <- p_excess +
    geom_vline(data = switch_vlines, aes(xintercept = date),
               linetype = "dotted", color = "red", linewidth = 0.5)
}

# Annotate annual excess
annual_labels <- annual_excess |>
  mutate(
    date = as.Date(sprintf("%d-07-01", year)),
    label = sprintf("%+.0f", annual_excess)
  ) |>
  filter(date >= min(dates), date <= max(dates))

if (nrow(annual_labels) > 0) {
  y_pos <- max(excess_summary$excess_upper, na.rm = TRUE) * 0.95
  p_excess <- p_excess +
    geom_text(data = annual_labels, aes(x = date, y = y_pos, label = label),
              size = 2.8, color = "gray20", fontface = "italic")
}

ggsave("../results/figures/counterfactual_excess_cases.png", p_excess,
       width = 13, height = 6, dpi = 150)
cat("  Saved: results/figures/counterfactual_excess_cases.png\n")

# ==============================================================================
# 9. FIGURE 3: SWITCH ATTRIBUTION BAR CHART
# ==============================================================================

if (!is.null(serotype_switches) && nrow(serotype_switches) > 0) {
  cat("\nComputing excess cases per serotype switch (0-6 month post-switch)...\n")

  switch_excess <- list()

  for (i in seq_len(nrow(serotype_switches))) {
    sw_date <- as.Date(serotype_switches$switch_date[i])
    sw_label <- serotype_switches$switch_event[i]

    # 0-6 month window post-switch
    window_start <- sw_date
    window_end <- sw_date %m+% months(6)

    # Find weeks in this window
    wk_idx <- which(dates >= window_start & dates < window_end)

    if (length(wk_idx) == 0) {
      cat(sprintf("  %s: no model weeks in post-switch window\n", sw_label))
      next
    }

    # Sum excess over window for each draw
    excess_per_draw <- rowSums(excess_draws[, wk_idx, drop = FALSE])

    switch_excess[[i]] <- tibble(
      switch_event = sw_label,
      switch_date = sw_date,
      excess_median = median(excess_per_draw),
      excess_lower = quantile(excess_per_draw, 0.025),
      excess_upper = quantile(excess_per_draw, 0.975),
      excess_mean = mean(excess_per_draw),
      n_weeks = length(wk_idx)
    )

    cat(sprintf("  %s: %.0f excess cases (95%% CrI: %.0f to %.0f) over %d weeks\n",
                sw_label, median(excess_per_draw),
                quantile(excess_per_draw, 0.025),
                quantile(excess_per_draw, 0.975),
                length(wk_idx)))
  }

  switch_excess_df <- bind_rows(switch_excess)

  if (nrow(switch_excess_df) > 0) {
    # Order by date
    switch_excess_df <- switch_excess_df |>
      mutate(switch_event = fct_reorder(switch_event, switch_date))

    p_switch <- ggplot(switch_excess_df,
                       aes(x = switch_event, y = excess_median)) +
      geom_col(aes(fill = excess_median > 0), width = 0.7, show.legend = FALSE) +
      geom_errorbar(aes(ymin = excess_lower, ymax = excess_upper),
                    width = 0.25, linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "solid", color = "gray30") +
      scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
      coord_flip() +
      labs(
        title = "Excess Cases per Serotype Switch (0-6 Months Post-Switch)",
        subtitle = "Median with 95% CrI; red = excess, blue = deficit",
        x = NULL,
        y = "Cumulative excess cases"
      )

    ggsave("../results/figures/counterfactual_switch_attribution.png", p_switch,
           width = 10, height = max(4, 1.5 + 0.8 * nrow(switch_excess_df)), dpi = 150)
    cat("  Saved: results/figures/counterfactual_switch_attribution.png\n")
  }
} else {
  cat("\nSkipping switch attribution (no serotype switch data).\n")
  switch_excess_df <- tibble()
}

# ==============================================================================
# 10. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("COUNTERFACTUAL ANALYSIS SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Total excess over full period
total_excess_per_draw <- rowSums(excess_draws)
cat(sprintf("Total estimated excess cases (2012-2022) from non-climate drivers:\n"))
cat(sprintf("  Median: %.0f cases\n", median(total_excess_per_draw)))
cat(sprintf("  95%% CrI: %.0f to %.0f\n",
            quantile(total_excess_per_draw, 0.025),
            quantile(total_excess_per_draw, 0.975)))

cat("\nAnnual breakdown:\n")
for (i in seq_len(nrow(annual_excess))) {
  cat(sprintf("  %d: %+.0f excess cases\n",
              annual_excess$year[i], annual_excess$annual_excess[i]))
}

if (!is.null(serotype_switches) && exists("switch_excess_df") && nrow(switch_excess_df) > 0) {
  cat("\nExcess cases per serotype switch (0-6 months post-switch):\n")
  for (i in seq_len(nrow(switch_excess_df))) {
    cat(sprintf("  %s: %.0f (95%% CrI: %.0f to %.0f)\n",
                switch_excess_df$switch_event[i],
                switch_excess_df$excess_median[i],
                switch_excess_df$excess_lower[i],
                switch_excess_df$excess_upper[i]))
  }

  # Which switch produced the most excess?
  max_idx <- which.max(switch_excess_df$excess_median)
  cat(sprintf("\nLargest excess: %s with %.0f excess cases\n",
              switch_excess_df$switch_event[max_idx],
              switch_excess_df$excess_median[max_idx]))
}

cat("\nOutput files:\n")
cat("  results/figures/counterfactual_rt.png\n")
cat("  results/figures/counterfactual_excess_cases.png\n")
if (exists("switch_excess_df") && nrow(switch_excess_df) > 0) {
  cat("  results/figures/counterfactual_switch_attribution.png\n")
}
