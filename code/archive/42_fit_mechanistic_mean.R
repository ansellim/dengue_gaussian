#!/usr/bin/env Rscript
# ==============================================================================
# 42_fit_mechanistic_mean.R
#
# Fit Model E: Mechanistic mean function (immunity dynamics) with GP residual.
#
# Key idea:
#   Model A (05/13): log(Rt) = mu + f_climate + f_residual (flat mean)
#   Model E (this):  log(Rt) = log(R0 * S_eff) + f_climate + f_residual
#
# The GP residual now captures deviations FROM the immunity prediction.
# If alpha shrinks -> immunity explains most of what the GP was absorbing.
#
# Input:
#   data/model_data.rds
#   data/susceptible_covariate.rds (if available, from 41_prepare_susceptible_covariate.R)
#   code/42_model_mechanistic_mean.stan
#   results/fit_model3.rds (Model A baseline for comparison)
#
# Output:
#   results/fit_mechanistic_mean.rds
#   results/figures/mechanistic_mean_comparison.png
#   results/figures/mechanistic_mean_variance.png
#   results/figures/mechanistic_mean_gp_diagnostic.png
#   results/figures/mechanistic_mean_R0.png
#   results/mechanistic_mean_summary.csv
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(patchwork)
library(loo)
library(lubridate)
library(mgcv)

# Set up paths
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- grep("--file=", args, value = TRUE)
  if (length(script_path) > 0) {
    setwd(dirname(normalizePath(sub("--file=", "", script_path))))
  }
}

dir.create("../results", showWarnings = FALSE)
dir.create("../results/figures", showWarnings = FALSE)

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank()))

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("FITTING MODEL E: MECHANISTIC MEAN (IMMUNITY DYNAMICS) + GP RESIDUAL\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data
df <- model_data$df
dates <- model_data$metadata$dates_model
N_model <- stan_data$N_model

cat(sprintf("  N_model = %d weeks\n", N_model))
cat(sprintf("  Cases range: %d - %d\n", min(stan_data$cases), max(stan_data$cases)))

# ==============================================================================
# 2. PREPARE log_S_eff (EFFECTIVE SUSCEPTIBLE FRACTION)
# ==============================================================================

cat("\nPreparing log(S_eff) mechanistic mean...\n")

if (file.exists("../data/susceptible_covariate.rds")) {
  cat("  Loading pre-computed susceptible_covariate.rds...\n")
  susc_data <- readRDS("../data/susceptible_covariate.rds")
  S_eff <- susc_data$S_eff
  log_S_eff <- log(S_eff)
  cat(sprintf("  S_eff range: %.4f - %.4f\n", min(S_eff), max(S_eff)))
} else {
  cat("  susceptible_covariate.rds not found, computing inline...\n")

  # --- Load serotype data ---
  SEROTYPE_DIR <- file.path(dirname(getwd()), "data")

  sero_early <- read_csv(
    file.path(SEROTYPE_DIR, "serotype_props_2013_2016.csv"),
    show_col_types = FALSE
  ) |>
    mutate(year_month = as.Date(paste0(Month, "-01"))) |>
    select(year_month, D1_prop, D2_prop, D3_prop, D4_prop)

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

  # --- GAM smooth serotype proportions ---
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

  # --- Identify serotype switch dates ---
  persistent_switches <- serotype |>
    filter(dominant != prev_dominant, !is.na(prev_dominant)) |>
    select(year_month, from = prev_dominant, to = dominant)

  cat("  Serotype switches:\n")
  for (i in seq_len(nrow(persistent_switches))) {
    s <- persistent_switches[i, ]
    cat(sprintf("    %s: %s -> %s\n", s$year_month, s$from, s$to))
  }

  # --- Interpolate serotype dominance to weekly dates ---
  switch_dates_vec <- persistent_switches$year_month
  all_months <- sort(unique(serotype$year_month))

  period_starts <- c(min(all_months), switch_dates_vec)
  period_ends   <- c(switch_dates_vec - months(1), max(all_months))
  period_sero   <- character(length(period_starts))

  for (k in seq_along(period_starts)) {
    mid_date <- period_starts[k] + floor((period_ends[k] - period_starts[k]) / 2)
    closest <- which.min(abs(serotype$year_month - mid_date))
    period_sero[k] <- serotype$dominant[closest]
  }

  dominance_periods <- tibble(
    period  = seq_along(period_starts),
    start   = period_starts,
    end     = period_ends,
    serotype = period_sero
  )

  # --- Assign each modeled week to a dominance period ---
  week_period <- rep(NA_integer_, N_model)
  for (i in seq_len(nrow(dominance_periods))) {
    in_period <- which(dates >= dominance_periods$start[i] &
                       dates <= (dominance_periods$end[i] + days(15)))
    week_period[in_period] <- i
  }
  for (i in which(is.na(week_period))) {
    dists <- abs(as.numeric(dates[i] - dominance_periods$start))
    week_period[i] <- which.min(dists)
  }

  # --- Cumulative depletion model ---
  N_pop <- 5500000
  expansion_factor <- 8
  initial_S <- 0.75

  cases_vec <- df$cases[(length(df$cases) - N_model + 1):length(df$cases)]

  S_eff <- numeric(N_model)
  for (p in unique(week_period)) {
    idx <- which(week_period == p)
    cum_infections <- cumsum(cases_vec[idx] * expansion_factor)
    S_eff[idx] <- initial_S * (1 - cum_infections / N_pop)
  }

  # Floor to avoid log(0)
  S_eff <- pmax(S_eff, 0.01)
  log_S_eff <- log(S_eff)

  cat(sprintf("  S_eff range: %.4f - %.4f\n", min(S_eff), max(S_eff)))
  cat(sprintf("  log_S_eff range: %.4f - %.4f\n", min(log_S_eff), max(log_S_eff)))
}

# ==============================================================================
# 3. COMPILE AND FIT MODEL
# ==============================================================================

cat("\nCompiling Model E (mechanistic mean)...\n")
cat("  log(Rt) = log_R0 + log_S_eff + beta_temp * temp + beta_rain * rain + f_residual(t)\n")
cat("  Prior: log_R0 ~ N(log(1.3), 0.3) => R0 median ~1.3\n\n")

model_e <- cmdstan_model("42_model_mechanistic_mean.stan")

# Add log_S_eff to stan_data
stan_data_e <- stan_data
stan_data_e$log_S_eff <- log_S_eff

cat("Fitting model...\n")
fit_e <- model_e$sample(
  data = stan_data_e,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  seed = 42,
  refresh = 200
)

# Save fit
fit_e$save_object("../results/fit_mechanistic_mean.rds")
cat("\nSaved: results/fit_mechanistic_mean.rds\n")

# ==============================================================================
# 4. DIAGNOSTICS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DIAGNOSTICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

diag <- fit_e$diagnostic_summary()
cat(sprintf("Divergent transitions: %d\n", sum(diag$num_divergent)))
cat(sprintf("Max treedepth exceeded: %d\n", sum(diag$num_max_treedepth)))

# ==============================================================================
# 5. PARAMETER SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PARAMETER SUMMARY (MODEL E)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

params_e <- c("log_R0", "R0", "alpha", "rho", "phi",
              "temp_effect", "rain_effect",
              "prop_mechanistic", "prop_climate", "prop_residual")
summ_e <- fit_e$summary(variables = params_e)
print(summ_e)

# ==============================================================================
# 6. LOAD MODEL A (BASELINE) FOR COMPARISON
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMPARISON: MODEL A (FLAT MEAN) VS MODEL E (MECHANISTIC MEAN)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

if (file.exists("../results/fit_model3.rds")) {
  fit_a <- readRDS("../results/fit_model3.rds")
  cat("  Loaded Model A (fit_model3.rds) for comparison\n\n")
  has_model_a <- TRUE
} else {
  cat("  WARNING: results/fit_model3.rds not found. Skipping comparison.\n\n")
  has_model_a <- FALSE
}

# ==============================================================================
# 7. LOO COMPARISON
# ==============================================================================

if (has_model_a) {
  cat("Computing LOO-CV...\n")

  loo_e <- fit_e$loo()
  loo_a <- fit_a$loo()

  loo_comp <- loo_compare(list("Model_A_flat_mean" = loo_a,
                               "Model_E_mechanistic" = loo_e))
  cat("\nLOO comparison (positive elpd_diff favors the first model):\n")
  print(loo_comp)
  cat("\n")
}

# ==============================================================================
# 8. GP HYPERPARAMETER COMPARISON
# ==============================================================================

if (has_model_a) {
  cat("GP hyperparameter comparison:\n")

  alpha_a <- fit_a$draws("alpha", format = "draws_matrix")
  alpha_e <- fit_e$draws("alpha", format = "draws_matrix")
  rho_a   <- fit_a$draws("rho", format = "draws_matrix")
  rho_e   <- fit_e$draws("rho", format = "draws_matrix")

  cat(sprintf("  alpha (GP amplitude):\n"))
  cat(sprintf("    Model A: median = %.4f [%.4f, %.4f]\n",
              median(alpha_a), quantile(alpha_a, 0.025), quantile(alpha_a, 0.975)))
  cat(sprintf("    Model E: median = %.4f [%.4f, %.4f]\n",
              median(alpha_e), quantile(alpha_e, 0.025), quantile(alpha_e, 0.975)))
  cat(sprintf("    Ratio (E/A): %.3f\n", median(alpha_e) / median(alpha_a)))

  cat(sprintf("  rho (GP length scale):\n"))
  cat(sprintf("    Model A: median = %.2f [%.2f, %.2f]\n",
              median(rho_a), quantile(rho_a, 0.025), quantile(rho_a, 0.975)))
  cat(sprintf("    Model E: median = %.2f [%.2f, %.2f]\n",
              median(rho_e), quantile(rho_e, 0.025), quantile(rho_e, 0.975)))
  cat("\n")
}

# ==============================================================================
# 9. VARIANCE DECOMPOSITION
# ==============================================================================

cat("Variance decomposition (Model E):\n")
prop_mech_e <- fit_e$summary("prop_mechanistic")
prop_clim_e <- fit_e$summary("prop_climate")
prop_resid_e <- fit_e$summary("prop_residual")

cat(sprintf("  Mechanistic (immunity): %.1f%% [%.1f%%, %.1f%%]\n",
            prop_mech_e$median * 100, prop_mech_e$q5 * 100, prop_mech_e$q95 * 100))
cat(sprintf("  Climate:                %.1f%% [%.1f%%, %.1f%%]\n",
            prop_clim_e$median * 100, prop_clim_e$q5 * 100, prop_clim_e$q95 * 100))
cat(sprintf("  Residual (GP):          %.1f%% [%.1f%%, %.1f%%]\n",
            prop_resid_e$median * 100, prop_resid_e$q5 * 100, prop_resid_e$q95 * 100))

if (has_model_a) {
  cat("\nVariance decomposition (Model A for reference):\n")
  prop_clim_a <- fit_a$summary("prop_climate")
  prop_resid_a <- fit_a$summary("prop_residual")
  cat(sprintf("  Climate:  %.1f%%\n", prop_clim_a$median * 100))
  cat(sprintf("  Residual: %.1f%%\n", prop_resid_a$median * 100))
}
cat("\n")

# ==============================================================================
# 10. R0 ESTIMATE
# ==============================================================================

cat("R0 estimate:\n")
R0_summ <- fit_e$summary("R0")
cat(sprintf("  R0 = %.3f [%.3f, %.3f] (95%% CrI)\n",
            R0_summ$median, R0_summ$q5, R0_summ$q95))
cat("  Tan et al. (2019) reference: R0 ~ 1.28 - 1.30\n\n")

# ==============================================================================
# 11. FIGURES
# ==============================================================================

cat("Generating figures...\n")

# --- Extract posteriors ---
Rt_e_draws <- fit_e$draws("Rt", format = "matrix")
f_resid_e_draws <- fit_e$draws("f_residual", format = "matrix")
mean_func_draws <- fit_e$draws("mean_function", format = "matrix")

Rt_e_summary <- tibble(
  date   = dates,
  median = apply(Rt_e_draws, 2, median),
  lower  = apply(Rt_e_draws, 2, quantile, 0.10),
  upper  = apply(Rt_e_draws, 2, quantile, 0.90),
  model  = "Model E (mechanistic)"
)

f_resid_e_summary <- tibble(
  date   = dates,
  median = apply(f_resid_e_draws, 2, median),
  lower  = apply(f_resid_e_draws, 2, quantile, 0.10),
  upper  = apply(f_resid_e_draws, 2, quantile, 0.90),
  model  = "Model E (mechanistic)"
)

mean_func_summary <- tibble(
  date   = dates,
  median = apply(mean_func_draws, 2, median),
  lower  = apply(mean_func_draws, 2, quantile, 0.10),
  upper  = apply(mean_func_draws, 2, quantile, 0.90)
)

if (has_model_a) {
  Rt_a_draws <- fit_a$draws("Rt", format = "matrix")
  f_resid_a_draws <- fit_a$draws("f_residual", format = "matrix")

  Rt_a_summary <- tibble(
    date   = dates,
    median = apply(Rt_a_draws, 2, median),
    lower  = apply(Rt_a_draws, 2, quantile, 0.10),
    upper  = apply(Rt_a_draws, 2, quantile, 0.90),
    model  = "Model A (flat mean)"
  )

  f_resid_a_summary <- tibble(
    date   = dates,
    median = apply(f_resid_a_draws, 2, median),
    lower  = apply(f_resid_a_draws, 2, quantile, 0.10),
    upper  = apply(f_resid_a_draws, 2, quantile, 0.90),
    model  = "Model A (flat mean)"
  )
}

# --------------------------------------------------------------------------
# Figure 1: mechanistic_mean_comparison.png (3-panel)
# --------------------------------------------------------------------------

cat("  Figure 1: mechanistic_mean_comparison.png\n")

model_colors <- c("Model A (flat mean)" = "steelblue",
                   "Model E (mechanistic)" = "darkorange")

# Panel A: Rt from both models
if (has_model_a) {
  Rt_both <- bind_rows(Rt_a_summary, Rt_e_summary)
} else {
  Rt_both <- Rt_e_summary
}

p1_top <- ggplot(Rt_both, aes(x = date, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(aes(y = median), linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = model_colors, name = NULL) +
  scale_fill_manual(values = model_colors, name = NULL) +
  labs(title = "A. Estimated Rt (median + 80% CrI)",
       y = "Rt", x = NULL) +
  theme(legend.position = "top",
        axis.text.x = element_blank())

# Panel B: Mechanistic mean function
p1_mid <- ggplot(mean_func_summary, aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "darkorange", alpha = 0.2) +
  geom_line(aes(y = median), color = "darkorange", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "B. Mechanistic mean: log(R0 * S_eff)",
       y = "log(R0 * S_eff)", x = NULL) +
  theme(axis.text.x = element_blank())

# Panel C: f_residual from both models
if (has_model_a) {
  resid_both <- bind_rows(f_resid_a_summary, f_resid_e_summary)
} else {
  resid_both <- f_resid_e_summary
}

p1_bot <- ggplot(resid_both, aes(x = date, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(aes(y = median), linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = model_colors, name = NULL) +
  scale_fill_manual(values = model_colors, name = NULL) +
  labs(title = "C. GP residual f_residual (has Model E's residual shrunk?)",
       y = "f_residual", x = "Date") +
  theme(legend.position = "top")

p_fig1 <- p1_top / p1_mid / p1_bot +
  plot_annotation(
    title = "Model A (flat mean) vs Model E (mechanistic immunity mean)",
    subtitle = "If immunity explains Rt variation, Model E's GP residual should be smaller"
  )

ggsave("../results/figures/mechanistic_mean_comparison.png", p_fig1,
       width = 13, height = 10, dpi = 150)
cat("    Saved: results/figures/mechanistic_mean_comparison.png\n")

# --------------------------------------------------------------------------
# Figure 2: mechanistic_mean_variance.png
# --------------------------------------------------------------------------

cat("  Figure 2: mechanistic_mean_variance.png\n")

var_e <- tibble(
  component = c("Mechanistic\n(immunity)", "Climate", "Residual\n(GP)"),
  proportion = c(prop_mech_e$median, prop_clim_e$median, prop_resid_e$median),
  model = "Model E (mechanistic mean)"
)

if (has_model_a) {
  var_a <- tibble(
    component = c("Mechanistic\n(immunity)", "Climate", "Residual\n(GP)"),
    proportion = c(0, prop_clim_a$median, prop_resid_a$median),
    model = "Model A (flat mean)"
  )
  var_both <- bind_rows(var_a, var_e)
} else {
  var_both <- var_e
}

var_both$component <- factor(var_both$component,
                             levels = c("Mechanistic\n(immunity)", "Climate", "Residual\n(GP)"))

p_fig2 <- ggplot(var_both, aes(x = component, y = proportion * 100, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%", proportion * 100)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Model A (flat mean)" = "steelblue",
                                "Model E (mechanistic mean)" = "darkorange"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 105)) +
  labs(
    title = "Variance Decomposition: Model A vs Model E",
    subtitle = "How much variance does the mechanistic immunity mean explain?",
    x = NULL,
    y = "Proportion of log(Rt) variance (%)"
  ) +
  theme(legend.position = "top")

ggsave("../results/figures/mechanistic_mean_variance.png", p_fig2,
       width = 8, height = 6, dpi = 150)
cat("    Saved: results/figures/mechanistic_mean_variance.png\n")

# --------------------------------------------------------------------------
# Figure 3: mechanistic_mean_gp_diagnostic.png (2-panel)
# --------------------------------------------------------------------------

cat("  Figure 3: mechanistic_mean_gp_diagnostic.png\n")

if (has_model_a) {
  alpha_df <- bind_rows(
    tibble(value = as.numeric(alpha_a), model = "Model A (flat mean)"),
    tibble(value = as.numeric(alpha_e), model = "Model E (mechanistic)")
  )

  rho_df <- bind_rows(
    tibble(value = as.numeric(rho_a), model = "Model A (flat mean)"),
    tibble(value = as.numeric(rho_e), model = "Model E (mechanistic)")
  )
} else {
  alpha_df <- tibble(value = as.numeric(alpha_e), model = "Model E (mechanistic)")
  rho_df   <- tibble(value = as.numeric(rho_e), model = "Model E (mechanistic)")
}

p3_left <- ggplot(alpha_df, aes(x = value, fill = model, color = model)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values = model_colors, name = NULL) +
  scale_color_manual(values = model_colors, name = NULL) +
  labs(title = "A. GP amplitude (alpha)",
       subtitle = "Shrinkage = immunity absorbs variance",
       x = "alpha", y = "Density") +
  theme(legend.position = "top")

p3_right <- ggplot(rho_df, aes(x = value, fill = model, color = model)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values = model_colors, name = NULL) +
  scale_color_manual(values = model_colors, name = NULL) +
  labs(title = "B. GP length scale (rho)",
       subtitle = "Change indicates different residual structure",
       x = "rho (weeks)", y = "Density") +
  theme(legend.position = "top")

p_fig3 <- p3_left | p3_right

ggsave("../results/figures/mechanistic_mean_gp_diagnostic.png", p_fig3,
       width = 12, height = 5, dpi = 150)
cat("    Saved: results/figures/mechanistic_mean_gp_diagnostic.png\n")

# --------------------------------------------------------------------------
# Figure 4: mechanistic_mean_R0.png
# --------------------------------------------------------------------------

cat("  Figure 4: mechanistic_mean_R0.png\n")

R0_draws <- as.numeric(fit_e$draws("R0", format = "draws_matrix"))

p_fig4 <- ggplot(tibble(R0 = R0_draws), aes(x = R0)) +
  geom_density(fill = "darkorange", alpha = 0.4, color = "darkorange", linewidth = 0.8) +
  geom_vline(xintercept = 1.3, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 1.3, y = Inf, label = "Tan et al. (2019)\nR0 ~ 1.3",
           vjust = 1.5, hjust = -0.1, color = "red", size = 3.5) +
  labs(
    title = "Posterior of R0 (basic reproduction number)",
    subtitle = sprintf("Median = %.3f [%.3f, %.3f] (95%% CrI)",
                       R0_summ$median, R0_summ$q5, R0_summ$q95),
    x = "R0",
    y = "Density"
  )

ggsave("../results/figures/mechanistic_mean_R0.png", p_fig4,
       width = 7, height = 5, dpi = 150)
cat("    Saved: results/figures/mechanistic_mean_R0.png\n")

# ==============================================================================
# 12. PRINT SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

if (has_model_a) {
  cat("LOO comparison:\n")
  print(loo_comp)
  cat("\n")
}

cat(sprintf("R0 estimate: %.3f [%.3f, %.3f]\n",
            R0_summ$median, R0_summ$q5, R0_summ$q95))
cat(sprintf("  Tan et al. reference: 1.28-1.30\n\n"))

cat("Variance decomposition (Model E):\n")
cat(sprintf("  Mechanistic (immunity): %.1f%%\n", prop_mech_e$median * 100))
cat(sprintf("  Climate:                %.1f%%\n", prop_clim_e$median * 100))
cat(sprintf("  Residual (GP):          %.1f%%\n", prop_resid_e$median * 100))

if (has_model_a) {
  cat(sprintf("\nGP amplitude alpha:\n"))
  cat(sprintf("  Model A: %.4f\n", median(as.numeric(alpha_a))))
  cat(sprintf("  Model E: %.4f\n", median(as.numeric(alpha_e))))
  cat(sprintf("  Ratio (E/A): %.3f\n", median(as.numeric(alpha_e)) / median(as.numeric(alpha_a))))
}

# ==============================================================================
# 13. SAVE SUMMARY CSV
# ==============================================================================

cat("\nSaving summary CSV...\n")

summary_rows <- list()

# R0
summary_rows[[1]] <- tibble(
  metric = "R0", estimate = R0_summ$median,
  lower = R0_summ$q5, upper = R0_summ$q95,
  note = "Tan et al. reference: 1.28-1.30"
)

# Variance decomposition (Model E)
summary_rows[[2]] <- tibble(
  metric = "prop_mechanistic_E", estimate = prop_mech_e$median,
  lower = prop_mech_e$q5, upper = prop_mech_e$q95, note = "Model E"
)
summary_rows[[3]] <- tibble(
  metric = "prop_climate_E", estimate = prop_clim_e$median,
  lower = prop_clim_e$q5, upper = prop_clim_e$q95, note = "Model E"
)
summary_rows[[4]] <- tibble(
  metric = "prop_residual_E", estimate = prop_resid_e$median,
  lower = prop_resid_e$q5, upper = prop_resid_e$q95, note = "Model E"
)

# GP amplitude
summary_rows[[5]] <- tibble(
  metric = "alpha_E", estimate = median(as.numeric(fit_e$draws("alpha", format = "draws_matrix"))),
  lower = quantile(as.numeric(fit_e$draws("alpha", format = "draws_matrix")), 0.025),
  upper = quantile(as.numeric(fit_e$draws("alpha", format = "draws_matrix")), 0.975),
  note = "Model E GP amplitude"
)

if (has_model_a) {
  summary_rows[[6]] <- tibble(
    metric = "alpha_A", estimate = median(as.numeric(alpha_a)),
    lower = quantile(as.numeric(alpha_a), 0.025),
    upper = quantile(as.numeric(alpha_a), 0.975),
    note = "Model A GP amplitude"
  )

  summary_rows[[7]] <- tibble(
    metric = "alpha_ratio_E_over_A",
    estimate = median(as.numeric(alpha_e)) / median(as.numeric(alpha_a)),
    lower = NA_real_, upper = NA_real_,
    note = "Ratio < 1 means immunity absorbed GP variance"
  )

  # LOO
  loo_diff <- loo_comp[2, "elpd_diff"]
  loo_se   <- loo_comp[2, "se_diff"]
  summary_rows[[8]] <- tibble(
    metric = "loo_elpd_diff", estimate = loo_diff,
    lower = loo_diff - 1.96 * loo_se, upper = loo_diff + 1.96 * loo_se,
    note = "Positive favors first model in loo_compare"
  )

  # Model A variance decomposition
  summary_rows[[9]] <- tibble(
    metric = "prop_climate_A", estimate = prop_clim_a$median,
    lower = prop_clim_a$q5, upper = prop_clim_a$q95, note = "Model A"
  )
  summary_rows[[10]] <- tibble(
    metric = "prop_residual_A", estimate = prop_resid_a$median,
    lower = prop_resid_a$q5, upper = prop_resid_a$q95, note = "Model A"
  )
}

summary_df <- bind_rows(summary_rows)
write_csv(summary_df, "../results/mechanistic_mean_summary.csv")
cat("  Saved: results/mechanistic_mean_summary.csv\n")

# ==============================================================================
# 14. DONE
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MECHANISTIC MEAN MODEL FITTING COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\nOutput files:\n")
cat("  results/fit_mechanistic_mean.rds                        - Model E fit object\n")
cat("  results/figures/mechanistic_mean_comparison.png          - 3-panel Rt comparison\n")
cat("  results/figures/mechanistic_mean_variance.png            - Variance decomposition\n")
cat("  results/figures/mechanistic_mean_gp_diagnostic.png       - GP alpha/rho posteriors\n")
cat("  results/figures/mechanistic_mean_R0.png                  - R0 posterior\n")
cat("  results/mechanistic_mean_summary.csv                     - Summary statistics\n")
