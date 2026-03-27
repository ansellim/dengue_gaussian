#!/usr/bin/env Rscript
# ==============================================================================
# 08_serotype_analysis.R
#
# Post-hoc analysis: Model 2 residual GP vs serotype dynamics
#
# Hypothesis: Serotype replacement events increase Rt because when a new
# serotype becomes dominant, the population has low immunity to it, creating
# a large susceptible pool. f_residual should show: a rise near serotype
# switches, followed by gradual decline during stable dominance.
#
# Input:
#   results/fit_model2.rds, data/model_data.rds
#   Serotype data from DengueFever project:
#     - serotype_props_2013_2016.csv (monthly, 2013-2016)
#     - monthly_sero_type_props_all_data.csv (monthly via weekly rows, 2017-2022)
#
# Output:
#   results/figures/serotype_panel.png
#   results/figures/serotype_switch_profiles.png
#   results/serotype_residual_monthly.csv
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

cat("Loading Model 2 fit and data...\n")

model_data <- readRDS("../data/model_data.rds")
fit2 <- readRDS("../results/fit_model2.rds")

dates <- model_data$metadata$dates_model
N_model <- length(dates)

# f_residual posterior draws: matrix (n_draws x N_model)
# f_residual(t) = log(Rt) - mu - X_t * beta
# Lives on the log(Rt) scale: >0 means Rt higher than covariates predict
f_residual_draws <- fit2$draws("f_residual", format = "matrix")
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
    title = "Residual GP from Model 2",
    subtitle = expression(
      f[residual](t) *
      ": unexplained variation in log(" * R[t] *
      ") after controlling for climate, Wolbachia, NPIs"
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
  plot_layout(heights = c(4, 2.5, 0.6, 2)) +
  plot_annotation(
    title = "Post-hoc Serotype Analysis: Residual GP vs Serotype Dynamics",
    subtitle = "Singapore, 2013-2022"
  )

ggsave("../results/figures/serotype_panel.png", p_panel,
       width = 14, height = 12, dpi = 150)
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

if (nrow(profiles_df) > 0) {
  p_profiles <- ggplot(profiles_df, aes(x = months_relative)) +
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
        "Temporal Profile of" ~ f[residual] ~ "Around Serotype Switch Events"
      ),
      subtitle = paste(
        "Red line = month of switch.",
        "Predicted: rise after switch (fresh susceptible pool),",
        "gradual decline (immunity accumulates)"
      ),
      x = "Months relative to switch",
      y = expression(f[residual] ~ "(log scale)"),
      caption = paste(
        "Shaded: 80% and 95% credible intervals",
        "of monthly-averaged residual GP"
      )
    )

  fig_height <- max(4, 3.5 * ceiling(nrow(persistent_switches) / 2))
  ggsave("../results/figures/serotype_switch_profiles.png", p_profiles,
         width = 12, height = fig_height, dpi = 150)
  cat("  Saved: results/figures/serotype_switch_profiles.png\n")
} else {
  cat("  WARNING: No persistent serotype switches found.\n")
}

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SEROTYPE ANALYSIS COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\n--- INTERPRETATION GUIDE ---\n\n")
cat("The serotype-immunity hypothesis predicts that f_residual should:\n")
cat("  1. RISE around serotype switches (new dominant serotype encounters\n")
cat("     a largely susceptible population -> Rt exceeds what covariates predict)\n")
cat("  2. DECLINE during stable dominance (herd immunity to the dominant\n")
cat("     serotype accumulates -> Rt falls back toward covariate prediction)\n")
cat("\nExamine the switch profile plots:\n")
cat("  - Clear separation of pre/post credible intervals -> supports hypothesis\n")
cat("  - Heavy overlap of intervals -> data cannot distinguish effect from noise\n")
cat("  - Flat or opposite pattern -> hypothesis not supported at national level\n")
cat("\nThis analysis is descriptive and hypothesis-generating. With ~3-4 major\n")
cat("switches in 10 years, there is no sample size for formal testing. Evidence\n")
cat("comes from whether f_residual's temporal profile matches the predicted\n")
cat("pattern in direction and timing across individual events.\n")

cat("\nOutput files:\n")
cat("  results/figures/serotype_panel.png          - Three-panel overview\n")
cat("  results/figures/serotype_switch_profiles.png - Event-centered profiles\n")
cat("  results/serotype_residual_monthly.csv        - Monthly data for further analysis\n")
