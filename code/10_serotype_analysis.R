#!/usr/bin/env Rscript
# ==============================================================================
# 10_serotype_analysis.R
#
# Post-hoc analysis: residual GP vs serotype dynamics
# Focus: Early warning system for Rt increases using serotype switching
#
# Questions addressed:
#   1. Are increases in Rt contemporaneous with, preceded by, or followed by
#      serotype switching?
#   2. Can entropy rise/fall signal upcoming Rt surges?
#   3. How soon before a serotype switch can we observe an Rt rise?
#
# Input:
#   results/fit_model.rds, data/model_data.rds
#   Serotype data:
#     - serotype_props_2013_2016.csv (monthly, 2013-2016)
#     - monthly_sero_type_props_all_data.csv (monthly via weekly rows, 2017-2022)
#
# Output:
#   results/figures/serotype_panel.png
#   results/figures/serotype_switch_profiles.png
#   results/figures/serotype_ccf.png
#   results/figures/serotype_entropy_profiles.png
#   results/figures/serotype_rt_entropy_overlay.png
#   results/serotype_residual_monthly.csv
#   results/serotype_ccf_results.csv
#   results/serotype_switch_timing.csv
#   results/serotype_post_switch_dynamics.csv
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)
library(lubridate)
library(mgcv)       # GAM multinomial logistic regression for serotype smoothing

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

# Path to serotype data (in separate DengueFever project, same Dropbox folder)
# Look in project data/ first, fall back to DengueFever sibling project
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
cat("POST-HOC SEROTYPE ANALYSIS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD MODEL FIT AND EXTRACT f_residual
# ==============================================================================

cat("Loading model fit and data...\n")

model_data <- readRDS("../data/model_data.rds")

fit_file <- "../results/fit_model.rds"
if (!file.exists(fit_file)) {
  stop("No model fit found. Run 06_fit_model.R first.")
}
fit <- readRDS(fit_file)
model_label <- "climate-only model"
cat("  Using climate-only model fit\n")

dates <- model_data$metadata$dates_model
N_model <- length(dates)

# f_residual posterior draws: matrix (n_draws x N_model)
# f_residual(t) = log(Rt) - mu - X_t * beta
# Lives on the log(Rt) scale: >0 means Rt higher than covariates predict
f_residual_draws <- fit$draws("f_residual", format = "matrix")
n_draws <- nrow(f_residual_draws)

cat(sprintf("  %d posterior draws, %d modeled weeks\n", n_draws, N_model))
cat(sprintf("  Date range: %s to %s\n", min(dates), max(dates)))

# ==============================================================================
# 2. AGGREGATE f_residual TO MONTHLY
# ==============================================================================

cat("\nAggregating f_residual to monthly resolution...\n")

# Map each week to its calendar month
week_month <- tibble(
  week_idx = 1:N_model,
  date = dates,
  year_month = floor_date(dates, "month")
)

month_groups <- week_month |>
  group_by(year_month) |>
  summarize(week_indices = list(week_idx), .groups = "drop")

# Average f_residual within each month for each posterior draw
n_months <- nrow(month_groups)
f_residual_monthly <- matrix(NA_real_, nrow = n_draws, ncol = n_months)

for (m in seq_len(n_months)) {
  idx <- month_groups$week_indices[[m]]
  f_residual_monthly[, m] <- rowMeans(f_residual_draws[, idx, drop = FALSE])
}

monthly_summary <- tibble(
  year_month = month_groups$year_month,
  median    = apply(f_residual_monthly, 2, median),
  lower_95  = apply(f_residual_monthly, 2, quantile, 0.025),
  upper_95  = apply(f_residual_monthly, 2, quantile, 0.975),
  lower_80  = apply(f_residual_monthly, 2, quantile, 0.10),
  upper_80  = apply(f_residual_monthly, 2, quantile, 0.90)
)

cat(sprintf("  Aggregated to %d months\n", n_months))

# ==============================================================================
# 3. LOAD AND COMBINE SEROTYPE DATA
# ==============================================================================

cat("\nLoading serotype data...\n")

# 2013-2016: monthly proportions
sero_early <- read_csv(
  file.path(SEROTYPE_DIR, "serotype_props_2013_2016.csv"),
  show_col_types = FALSE
) |>
  mutate(year_month = as.Date(paste0(Month, "-01"))) |>
  select(year_month, D1_prop, D2_prop, D3_prop, D4_prop)

# 2017-2022: weekly rows with monthly resolution — take one row per month
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
# 4. SMOOTH SEROTYPE PROPORTIONS AND COMPUTE METRICS
# ==============================================================================
#
# Raw monthly proportions are noisy because only a subset of cases are serotyped
# each month. Following Finch et al. (2025, Nat Comms), we smooth proportions
# using a GAM multinomial logistic regression before computing dominance.
# This filters sampling noise and avoids spurious month-to-month flips in the
# dominant serotype identity.
#
# Approach: fit a GAM with multinomial logistic link where each serotype's
# log-odds (relative to a reference) is a smooth function of time. We use
# mgcv::multinom with thin plate regression splines. Only data up to each
# time point is used in principle (Finch et al. use a forward-looking constraint),
# but for a post-hoc descriptive analysis we smooth over the full series.
# ==============================================================================

cat("\nSmoothing serotype proportions via GAM multinomial logistic regression...\n")

# Prepare data for mgcv::multinom: need counts or a matrix of proportions
# We have proportions — convert to a response matrix for multinom GAM
# mgcv expects a matrix of category counts (or proportions with weights)
n_sero_months <- nrow(serotype)
serotype <- serotype |>
  mutate(t_numeric = as.numeric(year_month - min(year_month)) / 30.44)  # months since start

# Smooth each serotype proportion independently via GAM with logit link,
# then renormalize so the four smoothed proportions sum to 1.
# This avoids issues with mgcv::multinom which expects indicator responses,
# and is equivalent in effect: we get a smooth trend for each serotype that
# filters monthly sampling noise while preserving the compositional constraint.
sero_names <- c("D1_prop", "D2_prop", "D3_prop", "D4_prop")
smooth_names <- c("D1_smooth", "D2_smooth", "D3_smooth", "D4_smooth")

smooth_mat <- matrix(NA_real_, nrow = nrow(serotype), ncol = 4)
for (j in 1:4) {
  # Clamp proportions away from 0/1 for logit transform stability
  p <- pmin(pmax(serotype[[sero_names[j]]], 0.001), 0.999)
  fit_j <- gam(p ~ s(t_numeric, k = 20, bs = "tp"),
               family = quasibinomial(), data = serotype)
  smooth_mat[, j] <- predict(fit_j, type = "response")
  cat(sprintf("  %s GAM deviance explained: %.1f%%\n",
              sero_names[j], summary(fit_j)$dev.expl * 100))
}

# Renormalize rows to sum to 1
row_sums <- rowSums(smooth_mat)
smooth_mat <- smooth_mat / row_sums

serotype <- serotype |>
  mutate(D1_smooth = smooth_mat[,1], D2_smooth = smooth_mat[,2],
         D3_smooth = smooth_mat[,3], D4_smooth = smooth_mat[,4])

# Dominant serotype from SMOOTHED proportions (argmax)
serotype <- serotype |>
  mutate(
    dominant_raw = case_when(
      D1_prop >= D2_prop & D1_prop >= D3_prop & D1_prop >= D4_prop ~ "DENV-1",
      D2_prop >= D1_prop & D2_prop >= D3_prop & D2_prop >= D4_prop ~ "DENV-2",
      D3_prop >= D1_prop & D3_prop >= D2_prop & D3_prop >= D4_prop ~ "DENV-3",
      TRUE ~ "DENV-4"
    ),
    dominant = case_when(
      D1_smooth >= D2_smooth & D1_smooth >= D3_smooth & D1_smooth >= D4_smooth ~ "DENV-1",
      D2_smooth >= D1_smooth & D2_smooth >= D3_smooth & D2_smooth >= D4_smooth ~ "DENV-2",
      D3_smooth >= D1_smooth & D3_smooth >= D2_smooth & D3_smooth >= D4_smooth ~ "DENV-3",
      TRUE ~ "DENV-4"
    ),
    # Shannon entropy computed on RAW proportions (reflects actual diversity observed)
    entropy = -(
      ifelse(D1_prop > 0, D1_prop * log(D1_prop), 0) +
      ifelse(D2_prop > 0, D2_prop * log(D2_prop), 0) +
      ifelse(D3_prop > 0, D3_prop * log(D3_prop), 0) +
      ifelse(D4_prop > 0, D4_prop * log(D4_prop), 0)
    ),
    # Also compute smoothed entropy for plotting
    entropy_smooth = -(
      ifelse(D1_smooth > 0, D1_smooth * log(D1_smooth), 0) +
      ifelse(D2_smooth > 0, D2_smooth * log(D2_smooth), 0) +
      ifelse(D3_smooth > 0, D3_smooth * log(D3_smooth), 0) +
      ifelse(D4_smooth > 0, D4_smooth * log(D4_smooth), 0)
    )
  )

# Identify serotype switches from smoothed dominant serotype
# No persistence filter needed — smoothing eliminates spurious flips
serotype <- serotype |>
  mutate(prev_dominant = lag(dominant))

persistent_switches <- serotype |>
  filter(dominant != prev_dominant, !is.na(prev_dominant)) |>
  select(year_month, from = prev_dominant, to = dominant)

cat(sprintf("\n  Serotype switches (from GAM-smoothed proportions):\n"))
for (i in seq_len(nrow(persistent_switches))) {
  s <- persistent_switches[i, ]
  cat(sprintf("    %s: %s -> %s\n", s$year_month, s$from, s$to))
}

# Report how many spurious flips the smoothing removed
n_raw_switches <- sum(
  serotype$dominant_raw != lag(serotype$dominant_raw),
  na.rm = TRUE
)
cat(sprintf("\n  Raw data had %d month-to-month dominant serotype changes\n",
            n_raw_switches))
cat(sprintf("  GAM smoothing reduced this to %d switch events\n",
            nrow(persistent_switches)))

# ==============================================================================
# 5. MERGE AND RESTRICT TO OVERLAPPING PERIOD
# ==============================================================================

analysis_df <- monthly_summary |>
  inner_join(
    serotype |> select(-prev_dominant, -t_numeric),
    by = "year_month"
  )

cat(sprintf("\n  Joined analysis data: %d months\n", nrow(analysis_df)))

write_csv(analysis_df, "../results/serotype_residual_monthly.csv")
cat("  Saved: results/serotype_residual_monthly.csv\n")

# ==============================================================================
# 6. PANEL PLOT
# ==============================================================================

cat("\nGenerating panel plot...\n")

sero_colors <- c(
  "DENV-1" = "#E41A1C",
  "DENV-2" = "#377EB8",
  "DENV-3" = "#4DAF4A",
  "DENV-4" = "#984EA3"
)

# Panel 1: f_residual with credible bands
p_residual <- ggplot(analysis_df, aes(x = year_month)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
              fill = "steelblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80),
              fill = "steelblue", alpha = 0.4) +
  geom_line(aes(y = median), color = "steelblue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(
    data = persistent_switches,
    aes(xintercept = year_month),
    linetype = "dotted", color = "red", linewidth = 0.7
  ) +
  labs(
    title = paste("Residual GP from", model_label),
    subtitle = expression(
      f[residual](t) *
      ": unexplained variation in log(" * R[t] *
      ") after controlling for climate"
    ),
    y = expression(f[residual] ~ "(log scale)"),
    x = NULL
  ) +
  theme(axis.text.x = element_blank())

# Panel 2: GAM-smoothed serotype proportions (stacked area)
smooth_long <- analysis_df |>
  select(year_month, D1_smooth, D2_smooth, D3_smooth, D4_smooth) |>
  pivot_longer(cols = -year_month, names_to = "serotype", values_to = "proportion") |>
  mutate(serotype = recode(serotype,
    D1_smooth = "DENV-1", D2_smooth = "DENV-2",
    D3_smooth = "DENV-3", D4_smooth = "DENV-4"
  ))

# Also prepare raw proportions as points for overlay
raw_long <- analysis_df |>
  select(year_month, D1_prop, D2_prop, D3_prop, D4_prop) |>
  pivot_longer(cols = -year_month, names_to = "serotype", values_to = "proportion") |>
  mutate(serotype = recode(serotype,
    D1_prop = "DENV-1", D2_prop = "DENV-2",
    D3_prop = "DENV-3", D4_prop = "DENV-4"
  ))

p_proportions <- ggplot() +
  geom_area(data = smooth_long,
            aes(x = year_month, y = proportion, fill = serotype),
            alpha = 0.6, position = "stack") +
  scale_fill_manual(values = sero_colors, name = "Serotype") +
  geom_vline(
    data = persistent_switches,
    aes(xintercept = year_month),
    linetype = "dotted", color = "red", linewidth = 0.7
  ) +
  labs(
    title = "GAM-smoothed serotype proportions",
    y = "Proportion",
    x = NULL
  ) +
  theme(axis.text.x = element_blank())

# Panel 3: dominant serotype bar (from smoothed proportions)
p_dominant <- ggplot(analysis_df, aes(x = year_month, y = 1, fill = dominant)) +
  geom_tile(height = 1) +
  scale_fill_manual(values = sero_colors, name = "Dominant\nserotype") +
  geom_vline(
    data = persistent_switches,
    aes(xintercept = year_month),
    linetype = "dotted", color = "red", linewidth = 0.7
  ) +
  labs(y = NULL, x = NULL) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    axis.text.x  = element_blank()
  )

# Panel 4: Shannon entropy (raw + smoothed)
p_entropy <- ggplot(analysis_df, aes(x = year_month)) +
  geom_point(aes(y = entropy), color = "darkorange", size = 1, alpha = 0.5) +
  geom_line(aes(y = entropy_smooth), color = "darkorange", linewidth = 0.8) +
  geom_vline(
    data = persistent_switches,
    aes(xintercept = year_month),
    linetype = "dotted", color = "red", linewidth = 0.7
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    y = "Shannon entropy",
    x = "Date",
    caption = paste(
      "Red dotted lines: serotype switches from GAM-smoothed proportions",
      "(Finch et al. 2025 approach). Points: raw entropy; line: smoothed."
    )
  )

# Assemble
p_panel <- p_residual / p_proportions / p_dominant / p_entropy +
  plot_layout(heights = c(2, 2.5, 0.6, 2)) +
  plot_annotation(
    title = "Post-hoc Serotype Analysis: Residual GP vs Serotype Dynamics",
    subtitle = "Singapore, 2013-2022"
  )

ggsave("../results/figures/serotype_panel.png", p_panel,
       width = 14, height = 9, dpi = 150)
cat("  Saved: results/figures/serotype_panel.png\n")

# ==============================================================================
# 7. TEMPORAL PROFILES AROUND SWITCH EVENTS
# ==============================================================================

cat("\nGenerating temporal profiles around serotype switches...\n")

WINDOW_BEFORE <- 6   # months before switch
WINDOW_AFTER  <- 12  # months after switch

switch_profiles <- list()

for (i in seq_len(nrow(persistent_switches))) {
  switch_date <- persistent_switches$year_month[i]
  switch_label <- sprintf(
    "%s -> %s (%s)",
    persistent_switches$from[i],
    persistent_switches$to[i],
    format(switch_date, "%b %Y")
  )

  # Window bounds
  window_start <- switch_date %m-% months(WINDOW_BEFORE)
  window_end   <- switch_date %m+% months(WINDOW_AFTER)

  # Find months in window that have posterior data
  month_idx <- which(
    month_groups$year_month >= window_start &
    month_groups$year_month <= window_end
  )
  if (length(month_idx) == 0) next

  window_months <- month_groups$year_month[month_idx]
  window_draws  <- f_residual_monthly[, month_idx, drop = FALSE]

  # Months relative to switch (0 = switch month)
  months_relative <- round(
    as.numeric(difftime(window_months, switch_date, units = "days")) / 30.44
  )

  switch_profiles[[i]] <- tibble(
    switch_event    = switch_label,
    switch_date     = switch_date,
    year_month      = window_months,
    months_relative = months_relative,
    median    = apply(window_draws, 2, median),
    lower_95  = apply(window_draws, 2, quantile, 0.025),
    upper_95  = apply(window_draws, 2, quantile, 0.975),
    lower_80  = apply(window_draws, 2, quantile, 0.10),
    upper_80  = apply(window_draws, 2, quantile, 0.90)
  )
}

profiles_df <- bind_rows(switch_profiles)

# Split switches chronologically into two groups so each fits on its own slide
# without the facet grid being chopped off. Threshold 2019-01-01 separates the
# 2013-2016 DENV-1/2 cycling period (4 switches) from the 2020-2021 DENV-2/3
# transition era (3 switches).
GROUP_THRESHOLD <- as.Date("2019-01-01")

make_profiles_plot <- function(df, subtitle_extra) {
  ggplot(df, aes(x = months_relative)) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
                fill = "steelblue", alpha = 0.2) +
    geom_ribbon(aes(ymin = lower_80, ymax = upper_80),
                fill = "steelblue", alpha = 0.4) +
    geom_line(aes(y = median), color = "steelblue", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "red",
               linewidth = 0.7) +
    facet_wrap(~ switch_event, ncol = 2, scales = "free_y") +
    labs(
      title = expression(
        "Temporal Profile of" ~ f[residual] ~ "Around Serotype Switches"
      ),
      subtitle = subtitle_extra,
      x = "Months relative to switch",
      y = expression(f[residual] ~ "(log scale)"),
      caption = paste(
        "Red line = month of switch.",
        "Shaded: 80% and 95% credible intervals of monthly-averaged residual GP"
      )
    )
}

if (nrow(profiles_df) > 0) {
  profiles_early <- profiles_df |> filter(switch_date < GROUP_THRESHOLD)
  profiles_late  <- profiles_df |> filter(switch_date >= GROUP_THRESHOLD)

  n_early <- length(unique(profiles_early$switch_event))
  n_late  <- length(unique(profiles_late$switch_event))

  # Figure 1: 2013-2016 DENV-1/2 cycling period
  if (nrow(profiles_early) > 0) {
    p1 <- make_profiles_plot(
      profiles_early,
      sprintf("Part 1: 2013--2016 DENV-1 / DENV-2 cycling period (%d events)", n_early)
    )
    h1 <- max(4, 3.5 * ceiling(n_early / 2))
    ggsave("../results/figures/serotype_switch_profiles_part1.png", p1,
           width = 12, height = h1, dpi = 150)
    cat("  Saved: results/figures/serotype_switch_profiles_part1.png\n")
  }

  # Figure 2: 2020-2021 DENV-2 <-> DENV-3 transition era
  if (nrow(profiles_late) > 0) {
    p2 <- make_profiles_plot(
      profiles_late,
      sprintf("Part 2: 2020--2021 DENV-2 / DENV-3 transition era (%d events)", n_late)
    )
    h2 <- max(4, 3.5 * ceiling(n_late / 2))
    ggsave("../results/figures/serotype_switch_profiles_part2.png", p2,
           width = 12, height = h2, dpi = 150)
    cat("  Saved: results/figures/serotype_switch_profiles_part2.png\n")
  }

  # Keep the combined figure for the report (one cohesive multi-panel)
  p_combined <- make_profiles_plot(profiles_df,
    "Red line = month of switch. Predicted: rise post-switch, gradual decline")
  fig_height <- max(4, 3.5 * ceiling(nrow(persistent_switches) / 2))
  ggsave("../results/figures/serotype_switch_profiles.png", p_combined,
         width = 12, height = fig_height, dpi = 150)
  cat("  Saved: results/figures/serotype_switch_profiles.png\n")
} else {
  cat("  WARNING: No persistent serotype switches found.\n")
}

# ==============================================================================
# 8. CROSS-CORRELATION ANALYSIS (CCF)
# ==============================================================================
#
# Question: Does Rt increase precede, coincide with, or follow serotype switches?
# Method: CCF between monthly median f_residual and serotype indicators.
#
# Sign convention (R's ccf(x, y)):
#   rho(k) = cor(x[t+k], y[t])   with x = f_residual, y = serotype indicator
#   Positive lag (+k) => f_residual at t+k correlated with indicator at t
#                     => indicator LEADS f_residual (early-warning direction)
#   Negative lag (-k) => f_residual LEADS indicator
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("CROSS-CORRELATION ANALYSIS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Compute serotype turnover indicators
#
# "Dominant serotype proportion" = proportion of whichever serotype is most
# prevalent in each month (from GAM-smoothed proportions).
#
# "Serotype turnover rate" = |delta_dominant|: absolute month-to-month change
# in dominant proportion. Captures the *speed* of serotype landscape disruption
# without the directional ambiguity of raw first differences. High values
# indicate the dominant serotype is rapidly losing or gaining ground.

analysis_df <- analysis_df |>
  mutate(
    # Proportion of the dominant serotype (smoothed)
    dominant_prop = pmap_dbl(
      list(D1_smooth, D2_smooth, D3_smooth, D4_smooth),
      ~ max(c(..1, ..2, ..3, ..4))
    ),
    # Absolute rate of change — serotype turnover rate
    abs_delta_dominant = abs(dominant_prop - lag(dominant_prop)),
    # 3-month rolling change in entropy
    delta_entropy_3m = entropy_smooth - lag(entropy_smooth, 3)
  )

# --- CCF 1: f_residual median vs smoothed entropy ---
# Hypothesis: entropy rises (co-circulation) -> new serotype meets susceptible
# population -> Rt increases -> cases surge.
# If entropy is an early warning for Rt, peak CCF at POSITIVE lag
# (entropy leads, f_residual follows).
ccf_entropy <- ccf(
  analysis_df$median,
  analysis_df$entropy_smooth,
  lag.max = 12,
  na.action = na.pass,
  plot = FALSE
)

# --- CCF 2: f_residual median vs serotype turnover rate ---
ccf_turnover <- ccf(
  analysis_df$median,
  analysis_df$abs_delta_dominant,
  lag.max = 12,
  na.action = na.pass,
  plot = FALSE
)

# Save CCF results
ccf_results <- tibble(
  lag = ccf_entropy$lag[, 1, 1],
  ccf_vs_entropy = ccf_entropy$acf[, 1, 1],
  ccf_vs_turnover_rate = ccf_turnover$acf[, 1, 1]
)
write_csv(ccf_results, "../results/serotype_ccf_results.csv")
cat("  Saved: results/serotype_ccf_results.csv\n")

# Print key lags
# Sign convention (see block header): positive lag => indicator LEADS f_residual.
cat("\nCross-correlation: f_residual vs Shannon entropy:\n")
cat("  (Positive lag = entropy leads Rt; Negative lag = Rt leads entropy)\n")
for (lag_val in c(-6, -3, -1, 0, 1, 3, 6, 10, 11, 12)) {
  idx <- which(ccf_results$lag == lag_val)
  if (length(idx) == 1) {
    cat(sprintf("  Lag %+3d months: r = %+.3f%s\n",
                lag_val, ccf_results$ccf_vs_entropy[idx],
                ifelse(lag_val > 0, "  (entropy leads Rt)", "")))
  }
}

cat("\nCross-correlation: f_residual vs serotype turnover rate:\n")
cat("  (Positive lag = turnover leads Rt; Negative lag = Rt leads turnover)\n")
for (lag_val in c(-6, -3, -1, 0, 1, 3, 6, 10, 11, 12)) {
  idx <- which(ccf_results$lag == lag_val)
  if (length(idx) == 1) {
    cat(sprintf("  Lag %+3d months: r = %+.3f%s\n",
                lag_val, ccf_results$ccf_vs_turnover_rate[idx],
                ifelse(lag_val > 0, "  (turnover leads Rt)", "")))
  }
}

# Identify peak CCF lags
peak_entropy_lag <- ccf_results$lag[which.max(abs(ccf_results$ccf_vs_entropy))]
peak_turnover_lag <- ccf_results$lag[which.max(abs(ccf_results$ccf_vs_turnover_rate))]
cat(sprintf("\nPeak |CCF| with entropy at lag %+d months\n", peak_entropy_lag))
cat(sprintf("Peak |CCF| with turnover rate at lag %+d months\n", peak_turnover_lag))

# Plot CCF
n_obs <- sum(!is.na(analysis_df$median) & !is.na(analysis_df$entropy_smooth))
ci_bound <- 1.96 / sqrt(n_obs)

ccf_plot_df <- ccf_results |>
  pivot_longer(-lag, names_to = "indicator", values_to = "ccf") |>
  mutate(
    indicator = recode(indicator,
      ccf_vs_entropy = "Shannon entropy (smoothed)",
      ccf_vs_turnover_rate = "Serotype turnover rate (|delta dominant|)"
    ),
    highlight = lag %in% c(8, 9)
  )

p_ccf <- ggplot(ccf_plot_df, aes(x = lag, y = ccf, fill = highlight)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_hline(yintercept = c(-ci_bound, ci_bound),
             linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  scale_fill_manual(
    values = c(`FALSE` = "steelblue", `TRUE` = "firebrick"),
    guide = "none"
  ) +
  facet_wrap(~ indicator, ncol = 1) +
  scale_x_continuous(breaks = seq(-12, 12, by = 3)) +
  labs(
    title = expression("Cross-Correlation:" ~ f[residual] ~ "vs Serotype Indicators"),
    subtitle = paste(
      "Bar at lag k = cor(f_residual at month t+k, indicator at month t).",
      "Dashed lines: 95% noise band. Red bars: lag +8/+9 (turnover peak)."
    ),
    x = "\u2190 Rt leads serotype indicator     |     Lag (months)     |     serotype indicator leads Rt \u2192",
    y = "Cross-correlation"
  )

ggsave("../results/figures/serotype_ccf.png", p_ccf,
       width = 10, height = 8, dpi = 150)
cat("  Saved: results/figures/serotype_ccf.png\n")

# ==============================================================================
# 9. PRE-SWITCH Rt DYNAMICS
# ==============================================================================
#
# For each switch: is f_residual rising in the months before the switch?
# Quantify lead time: how many months before the switch can we detect Rt rise?
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PRE-SWITCH Rt DYNAMICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

switch_timing <- list()

for (i in seq_len(nrow(persistent_switches))) {
  switch_date <- persistent_switches$year_month[i]
  switch_label <- sprintf("%s -> %s (%s)",
    persistent_switches$from[i], persistent_switches$to[i],
    format(switch_date, "%b %Y"))

  # Find month indices relative to this switch
  early_window <- switch_date %m-% months(6)  # 6-3 months before
  late_window  <- switch_date %m-% months(3)  # 3-0 months before

  early_idx <- which(
    month_groups$year_month >= early_window &
    month_groups$year_month < late_window
  )
  late_idx <- which(
    month_groups$year_month >= late_window &
    month_groups$year_month < switch_date
  )

  if (length(early_idx) < 2 || length(late_idx) < 2) {
    cat(sprintf("  %s: insufficient data for pre-switch analysis\n", switch_label))
    next
  }

  # Posterior probability that f_residual is rising (late > early)
  early_mean <- rowMeans(f_residual_monthly[, early_idx, drop = FALSE])
  late_mean  <- rowMeans(f_residual_monthly[, late_idx, drop = FALSE])
  p_rising <- mean(late_mean > early_mean)

  # Lead time: earliest month (going backward from switch) where median > 0
  pre_switch_months <- which(
    month_groups$year_month >= (switch_date %m-% months(12)) &
    month_groups$year_month < switch_date
  )
  pre_switch_medians <- apply(f_residual_monthly[, pre_switch_months, drop = FALSE],
                              2, median)
  pre_switch_dates <- month_groups$year_month[pre_switch_months]

  # Find first month (chronologically) where median > 0 in pre-switch window
  positive_idx <- which(pre_switch_medians > 0)
  if (length(positive_idx) > 0) {
    first_positive_date <- pre_switch_dates[min(positive_idx)]
    lead_months <- round(as.numeric(difftime(switch_date, first_positive_date,
                                             units = "days")) / 30.44)
  } else {
    first_positive_date <- NA
    lead_months <- NA
  }

  cat(sprintf("  %s:\n", switch_label))
  cat(sprintf("    P(f_residual rising before switch) = %.1f%%\n", p_rising * 100))
  cat(sprintf("    Lead time (median > 0): %s months\n",
              ifelse(is.na(lead_months), "N/A", as.character(lead_months))))

  switch_timing[[i]] <- tibble(
    switch_event = switch_label,
    switch_date = switch_date,
    from = persistent_switches$from[i],
    to = persistent_switches$to[i],
    p_rising = p_rising,
    lead_months = lead_months,
    first_positive_date = first_positive_date
  )
}

switch_timing_df <- bind_rows(switch_timing)
write_csv(switch_timing_df, "../results/serotype_switch_timing.csv")
cat("\n  Saved: results/serotype_switch_timing.csv\n")

# ==============================================================================
# 9b. POST-SWITCH Rt DYNAMICS
# ==============================================================================
#
# If serotype replacement DRIVES Rt increases, we expect f_residual to rise
# in the months AFTER a switch (new serotype meets susceptible population).
# This complements Section 9: pre-switch looks for anticipatory signals,
# post-switch confirms the causal mechanism.
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("POST-SWITCH Rt DYNAMICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

post_switch_results <- list()

for (i in seq_len(nrow(persistent_switches))) {
  switch_date <- persistent_switches$year_month[i]
  switch_label <- sprintf("%s -> %s (%s)",
    persistent_switches$from[i], persistent_switches$to[i],
    format(switch_date, "%b %Y"))

  # Post-switch windows: 0-3 months and 3-6 months after switch
  early_post <- switch_date %m+% months(0)
  mid_post   <- switch_date %m+% months(3)
  late_post  <- switch_date %m+% months(6)

  idx_0_3 <- which(
    month_groups$year_month >= early_post &
    month_groups$year_month < mid_post
  )
  idx_3_6 <- which(
    month_groups$year_month >= mid_post &
    month_groups$year_month < late_post
  )

  # Also get pre-switch baseline (-6 to -3 months before)
  baseline_start <- switch_date %m-% months(6)
  baseline_end   <- switch_date %m-% months(3)
  idx_baseline <- which(
    month_groups$year_month >= baseline_start &
    month_groups$year_month < baseline_end
  )

  if (length(idx_0_3) < 2 || length(idx_baseline) < 2) {
    cat(sprintf("  %s: insufficient data for post-switch analysis\n", switch_label))
    next
  }

  # Posterior probability that f_residual is elevated post-switch vs baseline
  baseline_mean <- rowMeans(f_residual_monthly[, idx_baseline, drop = FALSE])
  post_0_3_mean <- rowMeans(f_residual_monthly[, idx_0_3, drop = FALSE])

  p_elevated_0_3 <- mean(post_0_3_mean > baseline_mean)
  median_change_0_3 <- median(post_0_3_mean - baseline_mean)

  if (length(idx_3_6) >= 2) {
    post_3_6_mean <- rowMeans(f_residual_monthly[, idx_3_6, drop = FALSE])
    p_elevated_3_6 <- mean(post_3_6_mean > baseline_mean)
    median_change_3_6 <- median(post_3_6_mean - baseline_mean)
  } else {
    p_elevated_3_6 <- NA
    median_change_3_6 <- NA
  }

  # Posterior probability that f_residual > 0 in 0-3 months post-switch
  post_median_0_3 <- apply(f_residual_monthly[, idx_0_3, drop = FALSE], 2, median)
  p_positive_post <- mean(post_median_0_3 > 0)

  cat(sprintf("  %s:\n", switch_label))
  cat(sprintf("    P(Rt elevated 0-3m post-switch vs baseline) = %.1f%%\n",
              p_elevated_0_3 * 100))
  cat(sprintf("    Median change in f_residual (0-3m): %+.3f\n", median_change_0_3))
  if (!is.na(p_elevated_3_6)) {
    cat(sprintf("    P(Rt elevated 3-6m post-switch vs baseline) = %.1f%%\n",
                p_elevated_3_6 * 100))
  }

  post_switch_results[[i]] <- tibble(
    switch_event = switch_label,
    switch_date = switch_date,
    from = persistent_switches$from[i],
    to = persistent_switches$to[i],
    p_elevated_0_3m = p_elevated_0_3,
    median_change_0_3m = median_change_0_3,
    p_elevated_3_6m = p_elevated_3_6,
    median_change_3_6m = median_change_3_6,
    p_positive_post = p_positive_post
  )
}

post_switch_df <- bind_rows(post_switch_results)
write_csv(post_switch_df, "../results/serotype_post_switch_dynamics.csv")
cat("\n  Saved: results/serotype_post_switch_dynamics.csv\n")

# ==============================================================================
# 10. ENTROPY PROFILES AROUND CHANGEPOINTS
# ==============================================================================
#
# For each switch: plot entropy rise/fall in a window around the changepoint.
# Overlay f_residual to see temporal alignment.
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ENTROPY PROFILES AROUND CHANGEPOINTS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

entropy_profiles <- list()

for (i in seq_len(nrow(persistent_switches))) {
  switch_date <- persistent_switches$year_month[i]
  switch_label <- sprintf("%s -> %s (%s)",
    persistent_switches$from[i], persistent_switches$to[i],
    format(switch_date, "%b %Y"))

  window_start <- switch_date %m-% months(WINDOW_BEFORE)
  window_end   <- switch_date %m+% months(WINDOW_AFTER)

  window_data <- analysis_df |>
    filter(year_month >= window_start, year_month <= window_end) |>
    mutate(
      months_relative = round(
        as.numeric(difftime(year_month, switch_date, units = "days")) / 30.44
      ),
      switch_event = switch_label
    )

  if (nrow(window_data) > 0) {
    entropy_profiles[[i]] <- window_data
  }
}

entropy_profiles_df <- bind_rows(entropy_profiles)

if (nrow(entropy_profiles_df) > 0) {
  # Plot entropy profiles with f_residual overlay
  p_entropy_profiles <- ggplot(entropy_profiles_df, aes(x = months_relative)) +
    # Entropy (primary y-axis)
    geom_point(aes(y = entropy), color = "darkorange", size = 1.5, alpha = 0.6) +
    geom_line(aes(y = entropy_smooth), color = "darkorange", linewidth = 0.9) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "red", linewidth = 0.7) +
    facet_wrap(~ switch_event, ncol = 2, scales = "free_y") +
    labs(
      title = "Shannon Entropy Around Serotype Switch Events",
      subtitle = paste(
        "Red line = month of switch.",
        "Predicted: entropy rises before switch (co-circulation),",
        "then falls (new dominant emerges)"
      ),
      x = "Months relative to switch",
      y = "Shannon entropy",
      caption = "Points: raw entropy; line: GAM-smoothed entropy"
    )

  fig_height <- max(4, 3.5 * ceiling(nrow(persistent_switches) / 2))
  ggsave("../results/figures/serotype_entropy_profiles.png", p_entropy_profiles,
         width = 12, height = fig_height, dpi = 150)
  cat("  Saved: results/figures/serotype_entropy_profiles.png\n")
}

# ==============================================================================
# 11. COMBINED Rt + ENTROPY OVERLAY
# ==============================================================================
#
# For each switch: dual-panel plot showing f_residual (with CrI) and entropy
# on aligned time axes. This directly shows whether Rt and entropy co-move.
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMBINED Rt + ENTROPY OVERLAY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

if (nrow(entropy_profiles_df) > 0) {
  overlay_plots <- list()

  for (sw in unique(entropy_profiles_df$switch_event)) {
    sw_data <- entropy_profiles_df |> filter(switch_event == sw)

    # Top: f_residual with credible bands
    p_top <- ggplot(sw_data, aes(x = months_relative)) +
      geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
                  fill = "steelblue", alpha = 0.2) +
      geom_ribbon(aes(ymin = lower_80, ymax = upper_80),
                  fill = "steelblue", alpha = 0.4) +
      geom_line(aes(y = median), color = "steelblue", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "red",
                 linewidth = 0.7) +
      labs(
        title = sw,
        y = expression(f[residual] ~ "(log scale)"),
        x = NULL
      ) +
      theme(axis.text.x = element_blank())

    # Bottom: entropy
    p_bottom <- ggplot(sw_data, aes(x = months_relative)) +
      geom_point(aes(y = entropy), color = "darkorange", size = 1.5, alpha = 0.6) +
      geom_line(aes(y = entropy_smooth), color = "darkorange", linewidth = 0.9) +
      geom_vline(xintercept = 0, linetype = "dotted", color = "red",
                 linewidth = 0.7) +
      labs(
        y = "Shannon entropy",
        x = "Months relative to switch"
      )

    overlay_plots[[sw]] <- p_top / p_bottom + plot_layout(heights = c(2, 1))
  }

  # Combine all switch events
  n_switches <- length(overlay_plots)
  p_overlay_all <- wrap_plots(overlay_plots, ncol = 1) +
    plot_annotation(
      title = expression(
        "Combined" ~ f[residual] ~ "and Shannon Entropy Around Serotype Switches"
      ),
      subtitle = paste(
        "Top: residual GP (blue, with 80%/95% CrI).",
        "Bottom: Shannon entropy (orange).",
        "Red line: switch month."
      )
    )

  fig_height <- max(6, 5 * n_switches)
  ggsave("../results/figures/serotype_rt_entropy_overlay.png", p_overlay_all,
         width = 12, height = fig_height, dpi = 150)
  cat("  Saved: results/figures/serotype_rt_entropy_overlay.png\n")

  # --- Print assessment ---
  cat("\n--- TEMPORAL RELATIONSHIP ASSESSMENT ---\n\n")

  for (i in seq_len(nrow(persistent_switches))) {
    sw <- persistent_switches[i, ]
    sw_label <- sprintf("%s -> %s (%s)", sw$from, sw$to,
                        format(sw$year_month, "%b %Y"))
    sw_data <- entropy_profiles_df |>
      filter(switch_event == sw_label)

    if (nrow(sw_data) == 0) next

    # Check if f_residual is positive before the switch
    pre_data <- sw_data |> filter(months_relative < 0, months_relative >= -6)
    post_data <- sw_data |> filter(months_relative > 0, months_relative <= 6)

    pre_median <- if (nrow(pre_data) > 0) mean(pre_data$median) else NA
    post_median <- if (nrow(post_data) > 0) mean(post_data$median) else NA

    pre_entropy <- if (nrow(pre_data) > 0) mean(pre_data$entropy_smooth) else NA
    post_entropy <- if (nrow(post_data) > 0) mean(post_data$entropy_smooth) else NA

    cat(sprintf("  %s:\n", sw_label))
    if (!is.na(pre_median) && !is.na(post_median)) {
      direction <- ifelse(post_median > pre_median, "RISES", "FALLS")
      cat(sprintf("    f_residual: %s after switch (pre=%.3f, post=%.3f)\n",
                  direction, pre_median, post_median))
    }
    if (!is.na(pre_entropy) && !is.na(post_entropy)) {
      e_direction <- ifelse(post_entropy < pre_entropy, "FALLS", "RISES")
      cat(sprintf("    Entropy: %s after switch (pre=%.3f, post=%.3f)\n",
                  e_direction, pre_entropy, post_entropy))
    }
    cat("\n")
  }
}

# ==============================================================================
# 12. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SEROTYPE EARLY WARNING ANALYSIS COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\n--- INTERPRETATION GUIDE ---\n\n")
cat("CAUSAL HYPOTHESIS:\n")
cat("  Entropy rises (serotype co-circulation)\n")
cat("    -> new serotype encounters susceptible population\n")
cat("    -> Rt increases (f_residual > 0)\n")
cat("    -> cases surge\n")
cat("\n")
cat("If this causal chain holds, the CCF should show:\n")
cat("  - Negative CCF lag: entropy/turnover LEADS Rt (serotype change precedes Rt rise)\n")
cat("  - Post-switch P(elevated) >> 50%: Rt rises AFTER serotype switches\n")
cat("\n")
cat("EARLY WARNING SIGNALS:\n")
cat("  1. ENTROPY as early warning for Rt: if entropy rises before Rt,\n")
cat("     monitoring serotype diversity could predict Rt surges.\n")
cat("  2. Rt as early warning for CASES: since Rt increases precede case\n")
cat("     surges by definition (renewal equation), detecting Rt anomalies\n")
cat("     (f_residual > 0) provides lead time for case increases.\n")
cat("\n")
cat("The two-stage early warning chain:\n")
cat("  Entropy/serotype monitoring -> Rt anomaly detection -> case surge prediction\n")
cat("\nThis analysis is descriptive and hypothesis-generating. With ~3-4 major\n")
cat("switches in 10 years, there is no sample size for formal testing.\n")

cat("\nOutput files:\n")
cat("  results/figures/serotype_panel.png              - Four-panel overview\n")
cat("  results/figures/serotype_switch_profiles.png     - Event-centered f_residual\n")
cat("  results/figures/serotype_ccf.png                 - Cross-correlation analysis\n")
cat("  results/figures/serotype_entropy_profiles.png    - Entropy around switches\n")
cat("  results/figures/serotype_rt_entropy_overlay.png  - Combined Rt + entropy\n")
cat("  results/serotype_residual_monthly.csv            - Monthly data\n")
cat("  results/serotype_ccf_results.csv                 - CCF values\n")
cat("  results/serotype_switch_timing.csv               - Lead time estimates\n")
