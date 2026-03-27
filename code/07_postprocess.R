#!/usr/bin/env Rscript
# ==============================================================================
# 07_postprocess.R
#
# Posterior summaries, decomposition plots, and diagnostics for Model 3
# (climate-only: temperature + rainfall + residual GP)
#
# Input: results/fit_model3.rds, data/model_data.rds
# Output: results/figures/*.png, results/*.csv
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(patchwork)

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

# Create figures directory
dir.create("../results/figures", showWarnings = FALSE, recursive = TRUE)

# Set theme for plots
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank()))

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DENGUE RT ESTIMATION - POST-PROCESSING\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA AND FIT
# ==============================================================================

cat("Loading data and model fit...\n")

model_data <- readRDS("../data/model_data.rds")
df <- model_data$df
dates <- model_data$metadata$dates_model

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

N_model <- length(dates)
cat(sprintf("  Loaded %d modeled weeks\n", N_model))

# ==============================================================================
# 2. EXTRACT POSTERIORS
# ==============================================================================

cat("\nExtracting posterior summaries...\n")

# Rt trajectory summary
Rt_draws <- fit$draws("Rt", format = "matrix")

rt_summary <- tibble(
  date = dates,
  Rt_mean = colMeans(Rt_draws),
  Rt_median = apply(Rt_draws, 2, median),
  Rt_lower = apply(Rt_draws, 2, quantile, 0.025),
  Rt_upper = apply(Rt_draws, 2, quantile, 0.975),
  Rt_lower50 = apply(Rt_draws, 2, quantile, 0.25),
  Rt_upper50 = apply(Rt_draws, 2, quantile, 0.75)
)

# ==============================================================================
# 3. PARAMETER SUMMARIES
# ==============================================================================

cat("\nGenerating parameter summaries...\n")

params <- fit$summary(variables = c("mu", "alpha", "rho", "phi",
                                      "temp_effect", "rain_effect",
                                      "prop_climate", "prop_residual"))
write_csv(params, "../results/parameter_summary.csv")
cat("  Saved: results/parameter_summary.csv\n")

# Print key results
cat("\n--- Key Parameter Estimates ---\n")

# Effect sizes
effects <- fit$summary(variables = c("temp_effect", "rain_effect"))
cat("\nMultiplicative effects on Rt (per 1 SD change):\n")
for (i in 1:nrow(effects)) {
  cat(sprintf("  %s: %.3f (90%% CI: %.3f - %.3f)\n",
              effects$variable[i],
              effects$median[i],
              effects$q5[i],
              effects$q95[i]))
}

# Variance decomposition
decomp <- fit$summary(variables = c("prop_climate", "prop_residual"))
cat("\nVariance decomposition (proportion of log(Rt) variance):\n")
for (i in 1:nrow(decomp)) {
  cat(sprintf("  %s: %.1f%% (90%% CI: %.1f%% - %.1f%%)\n",
              decomp$variable[i],
              decomp$median[i] * 100,
              decomp$q5[i] * 100,
              decomp$q95[i] * 100))
}

# ==============================================================================
# 4. RT TRAJECTORY PLOT
# ==============================================================================

cat("\nGenerating Rt trajectory plots...\n")

# Observed cases for context
cases_df <- tibble(
  date = dates,
  cases = df$cases[(length(df$cases) - N_model + 1):length(df$cases)]
)

# Rt with cases (dual axis)
p_rt_cases <- ggplot() +
  geom_bar(data = cases_df, aes(x = date, y = cases / 500),
           stat = "identity", fill = "gray80", alpha = 0.7) +
  geom_ribbon(data = rt_summary,
              aes(x = date, ymin = Rt_lower, ymax = Rt_upper),
              fill = "steelblue", alpha = 0.3) +
  geom_ribbon(data = rt_summary,
              aes(x = date, ymin = Rt_lower50, ymax = Rt_upper50),
              fill = "steelblue", alpha = 0.5) +
  geom_line(data = rt_summary,
            aes(x = date, y = Rt_median), color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(
    name = expression(R[t]),
    sec.axis = sec_axis(~ . * 500, name = "Weekly Cases")
  ) +
  labs(
    title = "Dengue Rt Estimates (Climate-Only Model) with Weekly Cases",
    subtitle = "Singapore, 2012-2022",
    x = "Date"
  )

ggsave("../results/figures/rt_with_cases.png", p_rt_cases,
       width = 12, height = 6, dpi = 150)
cat("  Saved: results/figures/rt_with_cases.png\n")

# ==============================================================================
# 5. DECOMPOSITION PLOTS
# ==============================================================================

cat("\nGenerating decomposition plots...\n")

# Extract component posteriors
f_climate <- fit$draws("f_climate", format = "matrix")
f_residual <- fit$draws("f_residual", format = "matrix")

# Summarize each component
summarize_component <- function(draws, name) {
  tibble(
    date = dates,
    component = name,
    mean = colMeans(draws),
    median = apply(draws, 2, median),
    lower = apply(draws, 2, quantile, 0.025),
    upper = apply(draws, 2, quantile, 0.975)
  )
}

components <- bind_rows(
  summarize_component(f_climate, "Climate"),
  summarize_component(f_residual, "Residual GP")
)

# Decomposition plot
p_decomp <- ggplot(components, aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = median), color = "steelblue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ component, ncol = 1, scales = "free_y") +
  labs(
    title = "Decomposition of log(Rt) into Components",
    subtitle = "Climate (temp + rain) and Residual GP",
    x = "Date",
    y = "Effect on log(Rt)"
  )

ggsave("../results/figures/decomposition.png", p_decomp,
       width = 12, height = 7, dpi = 150)
cat("  Saved: results/figures/decomposition.png\n")

# ==============================================================================
# 6. POSTERIOR PREDICTIVE CHECK
# ==============================================================================

cat("\nGenerating posterior predictive checks...\n")

cases_pred <- fit$draws("cases_pred", format = "matrix")
cases_obs <- df$cases[(length(df$cases) - N_model + 1):length(df$cases)]

# Align dimensions
n_common <- min(ncol(cases_pred), length(cases_obs))
cases_pred <- cases_pred[, 1:n_common]
cases_obs_aligned <- cases_obs[1:n_common]

ppc_summary <- tibble(
  date = dates[1:n_common],
  observed = cases_obs_aligned,
  pred_median = apply(cases_pred, 2, median),
  pred_lower = apply(cases_pred, 2, quantile, 0.025),
  pred_upper = apply(cases_pred, 2, quantile, 0.975),
  pred_lower80 = apply(cases_pred, 2, quantile, 0.1),
  pred_upper80 = apply(cases_pred, 2, quantile, 0.9)
)

# Coverage
coverage_95 <- mean(cases_obs_aligned >= ppc_summary$pred_lower &
                    cases_obs_aligned <= ppc_summary$pred_upper)
coverage_80 <- mean(cases_obs_aligned >= ppc_summary$pred_lower80 &
                    cases_obs_aligned <= ppc_summary$pred_upper80)

cat(sprintf("\nPosterior predictive coverage:\n"))
cat(sprintf("  80%% interval: %.1f%%\n", 100 * coverage_80))
cat(sprintf("  95%% interval: %.1f%%\n", 100 * coverage_95))

# PPC plot
p_ppc <- ggplot(ppc_summary, aes(x = date)) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), fill = "steelblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = pred_lower80, ymax = pred_upper80), fill = "steelblue", alpha = 0.4) +
  geom_line(aes(y = pred_median), color = "steelblue", linewidth = 0.5) +
  geom_point(aes(y = observed), size = 0.5, alpha = 0.7) +
  labs(
    title = "Posterior Predictive Check (Climate-Only Model)",
    subtitle = sprintf("80%% coverage: %.1f%%, 95%% coverage: %.1f%%",
                       100 * coverage_80, 100 * coverage_95),
    x = "Date",
    y = "Weekly Cases"
  )

ggsave("../results/figures/ppc_model2.png", p_ppc,
       width = 12, height = 5, dpi = 150)
cat("  Saved: results/figures/ppc_model2.png\n")

# ==============================================================================
# 7. EFFECT SIZE PLOT
# ==============================================================================

cat("\nGenerating effect size plot...\n")

p_effects <- mcmc_areas(fit$draws(c("temp_effect", "rain_effect")),
                         prob = 0.8, prob_outer = 0.95) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Climate Effects on Rt (Multiplicative Scale)",
    subtitle = "Values < 1 indicate Rt reduction; > 1 indicate Rt increase",
    x = "Effect on Rt"
  )

ggsave("../results/figures/effect_sizes.png", p_effects,
       width = 8, height = 4, dpi = 150)
cat("  Saved: results/figures/effect_sizes.png\n")

# ==============================================================================
# 8. SAVE ALL SUMMARIES
# ==============================================================================

cat("\nSaving summary data...\n")

write_csv(rt_summary, "../results/rt_trajectories.csv")
cat("  Saved: results/rt_trajectories.csv\n")

write_csv(components, "../results/decomposition.csv")
cat("  Saved: results/decomposition.csv\n")

write_csv(ppc_summary, "../results/ppc_summary.csv")
cat("  Saved: results/ppc_summary.csv\n")

# ==============================================================================
# 9. FINAL SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("POST-PROCESSING COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\n--- KEY FINDINGS ---\n\n")

cat("Climate effects on Rt:\n")
for (i in 1:nrow(effects)) {
  cat(sprintf("  %s: %.3f (90%% CI: %.3f - %.3f)\n",
              effects$variable[i],
              effects$median[i],
              effects$q5[i],
              effects$q95[i]))
}

cat("\nVariance decomposition:\n")
for (i in 1:nrow(decomp)) {
  cat(sprintf("  %s: %.1f%%\n", decomp$variable[i], decomp$median[i] * 100))
}

cat("\n\nOutput files:\n")
cat("  results/figures/rt_with_cases.png - Rt with observed cases\n")
cat("  results/figures/decomposition.png - Component decomposition\n")
cat("  results/figures/ppc_model2.png - Posterior predictive check\n")
cat("  results/figures/effect_sizes.png - Climate effect sizes\n")
cat("  results/*.csv - Summary tables\n")
