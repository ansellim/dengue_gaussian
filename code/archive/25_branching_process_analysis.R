#!/usr/bin/env Rscript
# ==============================================================================
# 25_branching_process_analysis.R
#
# Connects GP-estimated Rt posterior to Galton-Watson branching process theory.
# Computes time-varying extinction probabilities, superspreading metrics,
# and analyses around serotype switch events.
#
# Theory (SPH6102 Week 2):
#   NegBin offspring distribution with mean Rt and overdispersion k.
#   PGF: G(s) = (1 + (Rt/k)(1-s))^(-k)
#   Extinction probability q is the smallest root of q = G(q) in [0,1].
#   When Rt <= 1, q = 1 (certain extinction).
#   When Rt > 1, q < 1 (positive probability of sustained outbreak).
#
# Input:
#   results/fit_model3.rds, data/model_data.rds
#   results/serotype_switch_timing.csv
#
# Output:
#   results/figures/branching_extinction_probability.png
#   results/figures/branching_serotype_switches.png
#   results/figures/branching_superspreading.png
#   results/figures/branching_offspring_distribution.png
#   results/branching_process_summary.csv
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
cat("BRANCHING PROCESS ANALYSIS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. EXTRACT Rt POSTERIOR AND OVERDISPERSION
# ==============================================================================

cat("Loading data and model fit...\n")

model_data <- readRDS("../data/model_data.rds")
dates <- model_data$metadata$dates_model
N_model <- length(dates)

# Load best available fit
fit_file <- NULL
for (f in c("../results/fit_model3.rds",
            "../results/fit_model2_short_gp.rds",
            "../results/fit_model2.rds")) {
  if (file.exists(f)) {
    fit_file <- f
    break
  }
}

if (is.null(fit_file)) {
  stop("No fitted model found. Run 13_fit_model3.R first.")
}

cat(sprintf("  Loading: %s\n", fit_file))
fit <- readRDS(fit_file)

# Rt draws: matrix (n_draws x N_model)
Rt_draws <- fit$draws("Rt", format = "matrix")
n_draws <- nrow(Rt_draws)

# phi (NegBin overdispersion) draws: vector of length n_draws
# Stan NegBin2: Var = mu + mu^2/phi, so phi is the dispersion parameter k
phi_draws <- as.vector(fit$draws("phi", format = "matrix"))

cat(sprintf("  %d posterior draws, %d modeled weeks\n", n_draws, N_model))
cat(sprintf("  Date range: %s to %s\n", min(dates), max(dates)))
cat(sprintf("  phi (k) posterior median: %.2f (95%% CrI: %.2f - %.2f)\n",
            median(phi_draws),
            quantile(phi_draws, 0.025),
            quantile(phi_draws, 0.975)))

# ==============================================================================
# 2. TIME-VARYING EXTINCTION PROBABILITY
# ==============================================================================

cat("\nComputing time-varying extinction probabilities...\n")

#' Compute extinction probability for NegBin(mu = Rt, k) offspring.
#'
#' PGF: G(s) = (1 + (Rt/k)(1-s))^(-k)
#' When Rt <= 1: q = 1.
#' When Rt > 1: solve q = G(q) via uniroot on h(s) = G(s) - s in (0, 1).
#'
#' @param Rt Reproduction number (scalar)
#' @param k  Overdispersion parameter (scalar)
#' @return   Extinction probability q in [0, 1]
extinction_prob <- function(Rt, k) {
  if (Rt <= 1) return(1.0)
  pgf <- function(s) (1 + (Rt / k) * (1 - s))^(-k)
  h <- function(s) pgf(s) - s
  # h(0) = (1 + Rt/k)^(-k) > 0; h(1) = 0 (trivial root)
  # For Rt > 1 there is a non-trivial root in (0, 1)
  root <- tryCatch(
    uniroot(h, lower = 1e-12, upper = 1 - 1e-12, tol = 1e-10)$root,
    error = function(e) {
      # Fallback: fixed-point iteration
      q <- 0.5
      for (i in 1:500) {
        q_new <- pgf(q)
        if (abs(q_new - q) < 1e-12) break
        q <- q_new
      }
      q
    }
  )
  return(root)
}

# Vectorised over a single draw (all weeks)
extinction_prob_vec <- Vectorize(extinction_prob)

# Compute q for each draw and each week
# Use a thinned set if n_draws is very large to keep runtime manageable
max_draws <- min(n_draws, 1000)
draw_idx <- sort(sample.int(n_draws, max_draws))

cat(sprintf("  Using %d posterior draws (thinned from %d)...\n",
            max_draws, n_draws))

q_matrix <- matrix(NA_real_, nrow = max_draws, ncol = N_model)
pb <- txtProgressBar(min = 0, max = max_draws, style = 3)
for (i in seq_along(draw_idx)) {
  d <- draw_idx[i]
  q_matrix[i, ] <- extinction_prob_vec(Rt_draws[d, ], phi_draws[d])
  setTxtProgressBar(pb, i)
}
close(pb)

# Summarise extinction probability
q_summary <- tibble(
  date        = dates,
  q_median    = apply(q_matrix, 2, median),
  q_lower     = apply(q_matrix, 2, quantile, 0.025),
  q_upper     = apply(q_matrix, 2, quantile, 0.975),
  p_sustain_median = 1 - apply(q_matrix, 2, median),
  p_sustain_lower  = 1 - apply(q_matrix, 2, quantile, 0.975),
  p_sustain_upper  = 1 - apply(q_matrix, 2, quantile, 0.025)
)

# Rt summary (recompute for plotting)
rt_summary <- tibble(
  date      = dates,
  Rt_median = apply(Rt_draws, 2, median),
  Rt_lower  = apply(Rt_draws, 2, quantile, 0.025),
  Rt_upper  = apply(Rt_draws, 2, quantile, 0.975)
)

cat("\n  Extinction probability computed.\n")

# ==============================================================================
# 3. SEROTYPE SWITCH ANALYSIS
# ==============================================================================

cat("\nAnalysing branching dynamics around serotype switches...\n")

switch_file <- "../results/serotype_switch_timing.csv"
if (file.exists(switch_file)) {
  switches <- read_csv(switch_file, show_col_types = FALSE) |>
    mutate(switch_date = as.Date(switch_date))

  cat(sprintf("  Loaded %d serotype switch events\n", nrow(switches)))

  # For each switch, compare P(outbreak sustains) in 3 months before vs after
  window_days <- 91  # ~3 months

  switch_comparison <- switches |>
    rowwise() |>
    mutate(
      before_start = switch_date - window_days,
      after_end    = switch_date + window_days,
      idx_before   = list(which(dates >= before_start & dates < switch_date)),
      idx_after    = list(which(dates >= switch_date & dates <= after_end)),
      p_sustain_before = if (length(idx_before) > 0)
        mean(q_summary$p_sustain_median[idx_before]) else NA_real_,
      p_sustain_after  = if (length(idx_after) > 0)
        mean(q_summary$p_sustain_median[idx_after]) else NA_real_,
      delta_p_sustain  = p_sustain_after - p_sustain_before
    ) |>
    ungroup() |>
    select(switch_event, switch_date, from, to,
           p_sustain_before, p_sustain_after, delta_p_sustain)

  cat("\n  Before vs after switch (mean P(outbreak sustains)):\n")
  for (i in 1:nrow(switch_comparison)) {
    cat(sprintf("    %s: %.3f -> %.3f (delta = %+.3f)\n",
                switch_comparison$switch_event[i],
                switch_comparison$p_sustain_before[i],
                switch_comparison$p_sustain_after[i],
                switch_comparison$delta_p_sustain[i]))
  }
} else {
  cat("  serotype_switch_timing.csv not found; skipping switch analysis.\n")
  switches <- NULL
  switch_comparison <- NULL
}

# ==============================================================================
# 4. SUPERSPREADING ANALYSIS
# ==============================================================================

cat("\nComputing superspreading metrics...\n")

# For each week, compute:
#   (a) P(offspring > 20) under NegBin(Rt, k)
#   (b) Expected proportion of transmission from top 20% of cases

# Use posterior medians for a clean time series
Rt_med <- apply(Rt_draws, 2, median)
k_med  <- median(phi_draws)

# (a) P(offspring > 20)
p_super <- 1 - pnbinom(20, size = k_med, mu = Rt_med)

# (b) Proportion of transmission from top 20% (Lorenz/Gini approach)
#     For NegBin(mu, k), compute via Monte Carlo for representative values
compute_top20_share <- function(Rt, k, n_mc = 50000) {
  if (Rt < 0.01) return(0)
  offspring <- rnbinom(n_mc, size = k, mu = Rt)
  offspring_sorted <- sort(offspring)
  # Cumulative share of transmission
  cum_share <- cumsum(offspring_sorted) / sum(offspring_sorted)
  # Top 20% means bottom 80% cutoff
  cutoff_idx <- floor(0.8 * n_mc)
  if (cutoff_idx < 1) return(1)
  bottom80_share <- cum_share[cutoff_idx]
  top20_share <- 1 - bottom80_share
  return(top20_share)
}

# Compute for a grid of Rt values and interpolate for speed
set.seed(42)
Rt_grid <- seq(0, max(Rt_med) + 0.5, by = 0.05)
top20_grid <- sapply(Rt_grid, compute_top20_share, k = k_med)
top20_share <- approx(Rt_grid, top20_grid, xout = Rt_med, rule = 2)$y

superspreading_df <- tibble(
  date         = dates,
  Rt_median    = Rt_med,
  p_offspring_gt20 = p_super,
  top20_share  = top20_share
)

cat(sprintf("  Median P(offspring > 20): %.4f\n", median(p_super)))
cat(sprintf("  Median top-20%% transmission share: %.1f%%\n",
            100 * median(top20_share)))

# ==============================================================================
# 5. FIGURES
# ==============================================================================

cat("\nGenerating figures...\n")

switch_dates <- if (!is.null(switches)) switches$switch_date else NULL

# --- Figure 1: Extinction probability time series ---

p_rt <- ggplot(rt_summary, aes(x = date)) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper),
              fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = Rt_median), color = "steelblue", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  { if (!is.null(switch_dates))
    geom_vline(xintercept = switch_dates, linetype = "dotted",
               color = "gray40", alpha = 0.7) } +
  labs(x = NULL, y = expression(R[t]),
       title = "Branching Process Analysis of Dengue Rt",
       subtitle = "Singapore, 2012-2022") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")

# Identify sustained-outbreak periods (median P(sustain) > 0.5)
sustain_periods <- q_summary |>
  mutate(sustained = p_sustain_median > 0.5)

p_q <- ggplot(q_summary, aes(x = date)) +
  geom_ribbon(aes(ymin = q_lower, ymax = q_upper),
              fill = "darkorange", alpha = 0.3) +
  geom_line(aes(y = q_median), color = "darkorange", linewidth = 0.8) +
  # Highlight periods where P(sustain) > 0.5, i.e. q < 0.5
  geom_rug(data = filter(sustain_periods, sustained),
           aes(x = date), sides = "b", color = "red3", alpha = 0.5,
           length = unit(0.05, "npc")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  { if (!is.null(switch_dates))
    geom_vline(xintercept = switch_dates, linetype = "dotted",
               color = "gray40", alpha = 0.7) } +
  labs(x = "Date",
       y = "Extinction probability (q)",
       caption = "Red ticks: P(outbreak sustains) > 0.5. Dotted lines: serotype switches.") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(limits = c(0, 1))

p_fig1 <- p_rt / p_q + plot_layout(heights = c(1, 1))
ggsave("../results/figures/branching_extinction_probability.png", p_fig1,
       width = 12, height = 8, dpi = 150)
cat("  Saved: results/figures/branching_extinction_probability.png\n")

# --- Figure 2: Serotype switch comparison ---

if (!is.null(switch_comparison)) {
  switch_long <- switch_comparison |>
    select(switch_event, p_sustain_before, p_sustain_after) |>
    pivot_longer(cols = c(p_sustain_before, p_sustain_after),
                 names_to = "period", values_to = "p_sustain") |>
    mutate(
      period = recode(period,
                      "p_sustain_before" = "3 months before",
                      "p_sustain_after"  = "3 months after"),
      period = factor(period, levels = c("3 months before", "3 months after"))
    )

  p_fig2 <- ggplot(switch_long,
                    aes(x = switch_event, y = p_sustain, fill = period)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("3 months before" = "steelblue",
                                 "3 months after"  = "tomato")) +
    labs(
      title = "Outbreak Sustainability Around Serotype Switches",
      subtitle = "Mean P(outbreak sustains) in 3-month windows",
      x = NULL, y = "P(outbreak sustains)", fill = NULL
    ) +
    coord_flip() +
    theme(legend.position = "bottom")

  ggsave("../results/figures/branching_serotype_switches.png", p_fig2,
         width = 10, height = 5, dpi = 150)
  cat("  Saved: results/figures/branching_serotype_switches.png\n")
}

# --- Figure 3: Superspreading time series ---

p_super_ts <- ggplot(superspreading_df, aes(x = date)) +
  geom_line(aes(y = p_offspring_gt20), color = "purple4", linewidth = 0.7) +
  { if (!is.null(switch_dates))
    geom_vline(xintercept = switch_dates, linetype = "dotted",
               color = "gray40", alpha = 0.7) } +
  labs(x = NULL, y = "P(offspring > 20)",
       title = "Superspreading Potential Over Time") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")

p_top20_ts <- ggplot(superspreading_df, aes(x = date)) +
  geom_line(aes(y = top20_share), color = "darkgreen", linewidth = 0.7) +
  { if (!is.null(switch_dates))
    geom_vline(xintercept = switch_dates, linetype = "dotted",
               color = "gray40", alpha = 0.7) } +
  labs(x = "Date", y = "Top-20% transmission share",
       caption = "Proportion of total transmission attributed to the top 20% of cases") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = scales::percent_format())

p_fig3 <- p_super_ts / p_top20_ts
ggsave("../results/figures/branching_superspreading.png", p_fig3,
       width = 12, height = 7, dpi = 150)
cat("  Saved: results/figures/branching_superspreading.png\n")

# --- Figure 4: Offspring distributions for selected weeks ---

# Select representative weeks: low Rt, Rt ~ 1, high Rt, serotype switch
week_candidates <- tibble(date = dates, Rt = Rt_med) |>
  mutate(week_idx = row_number())

idx_low  <- week_candidates |> filter(Rt < 0.6) |> slice_sample(n = 1)
idx_one  <- week_candidates |> filter(abs(Rt - 1) < 0.1) |> slice_sample(n = 1)
idx_high <- week_candidates |> slice_max(Rt, n = 1)

# Near a serotype switch
idx_switch <- NULL
if (!is.null(switches) && nrow(switches) > 0) {
  target_date <- switches$switch_date[1]
  idx_switch <- week_candidates |>
    mutate(dist = abs(as.numeric(date - target_date))) |>
    slice_min(dist, n = 1) |>
    select(-dist)
}

selected_weeks <- bind_rows(
  idx_low  |> mutate(label = sprintf("Low Rt (%.2f, %s)", Rt, date)),
  idx_one  |> mutate(label = sprintf("Rt ~ 1 (%.2f, %s)", Rt, date)),
  idx_high |> mutate(label = sprintf("High Rt (%.2f, %s)", Rt, date)),
  idx_switch |> mutate(label = sprintf("Switch (%.2f, %s)", Rt, date))
)

offspring_x <- 0:40
offspring_dists <- selected_weeks |>
  rowwise() |>
  do({
    row <- .
    probs <- dnbinom(offspring_x, size = k_med, mu = row$Rt)
    tibble(label = row$label, x = offspring_x, prob = probs)
  }) |>
  ungroup() |>
  mutate(label = factor(label, levels = selected_weeks$label))

p_fig4 <- ggplot(offspring_dists, aes(x = x, y = prob)) +
  geom_col(fill = "steelblue", alpha = 0.8, width = 0.8) +
  facet_wrap(~ label, scales = "free_y", ncol = 2) +
  labs(
    title = "Offspring Distributions Under NegBin Branching Model",
    subtitle = sprintf("Using posterior median k = %.2f", k_med),
    x = "Number of secondary cases",
    y = "Probability"
  )

ggsave("../results/figures/branching_offspring_distribution.png", p_fig4,
       width = 10, height = 7, dpi = 150)
cat("  Saved: results/figures/branching_offspring_distribution.png\n")

# ==============================================================================
# 6. SUMMARY STATISTICS
# ==============================================================================

cat("\nComputing summary statistics...\n")

# Mean extinction probability by year
yearly_summary <- q_summary |>
  mutate(year = year(date)) |>
  group_by(year) |>
  summarise(
    mean_q           = mean(q_median),
    mean_p_sustain   = mean(p_sustain_median),
    weeks_sustained  = sum(p_sustain_median > 0.5),
    mean_Rt          = mean(rt_summary$Rt_median[rt_summary$date %in% date]),
    .groups = "drop"
  )

cat("\n  Mean extinction probability and sustainability by year:\n")
for (i in 1:nrow(yearly_summary)) {
  cat(sprintf("    %d: q = %.3f, P(sustain) = %.3f, weeks with P(sustain)>0.5: %d\n",
              yearly_summary$year[i],
              yearly_summary$mean_q[i],
              yearly_summary$mean_p_sustain[i],
              yearly_summary$weeks_sustained[i]))
}

# Correlation between Rt and q
cor_rt_q <- cor(rt_summary$Rt_median, q_summary$q_median, method = "spearman")
cat(sprintf("\n  Spearman correlation (Rt vs q): %.3f\n", cor_rt_q))

# Serotype switch significance (simple permutation test on delta)
switch_signif <- NULL
if (!is.null(switch_comparison)) {
  mean_delta <- mean(switch_comparison$delta_p_sustain, na.rm = TRUE)
  n_positive <- sum(switch_comparison$delta_p_sustain > 0, na.rm = TRUE)
  n_total <- sum(!is.na(switch_comparison$delta_p_sustain))
  cat(sprintf("\n  Serotype switches: %d/%d show increased P(sustain) after switch\n",
              n_positive, n_total))
  cat(sprintf("  Mean delta P(sustain): %+.3f\n", mean_delta))
  switch_signif <- tibble(
    metric = c("n_switches_increase", "n_switches_total", "mean_delta_p_sustain"),
    value = c(n_positive, n_total, round(mean_delta, 4))
  )
}

# Assemble and save summary CSV
summary_out <- bind_rows(
  yearly_summary |>
    mutate(category = "yearly") |>
    pivot_longer(cols = c(mean_q, mean_p_sustain, weeks_sustained, mean_Rt),
                 names_to = "metric", values_to = "value") |>
    mutate(label = paste0(year, "_", metric)) |>
    select(category, label, metric, value),
  tibble(category = "correlation", label = "Rt_vs_q_spearman",
         metric = "spearman_rho", value = round(cor_rt_q, 4)),
  if (!is.null(switch_signif))
    switch_signif |> mutate(category = "serotype_switch", label = metric)
)

write_csv(summary_out, "../results/branching_process_summary.csv")
cat("  Saved: results/branching_process_summary.csv\n")

# ==============================================================================
# DONE
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("BRANCHING PROCESS ANALYSIS COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\nOutput files:\n")
cat("  results/figures/branching_extinction_probability.png\n")
cat("  results/figures/branching_serotype_switches.png\n")
cat("  results/figures/branching_superspreading.png\n")
cat("  results/figures/branching_offspring_distribution.png\n")
cat("  results/branching_process_summary.csv\n")
