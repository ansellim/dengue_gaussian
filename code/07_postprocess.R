#!/usr/bin/env Rscript
# ==============================================================================
# 07_postprocess.R
#
# Posterior summaries, decomposition plots, and model comparison
#
# Input: results/fit_model*.rds, data/model_data.rds
# Output: results/figures/*.png, results/summary_tables.csv
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(patchwork)
library(loo)

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
# 1. LOAD DATA AND FITS
# ==============================================================================

cat("Loading data and model fits...\n")

model_data <- readRDS("../data/model_data.rds")
df <- model_data$df
dates <- model_data$metadata$dates_model

fit0 <- readRDS("../results/fit_model0.rds")
fit1 <- readRDS("../results/fit_model1.rds")
fit2 <- readRDS("../results/fit_model2.rds")

N_model <- length(dates)
cat(sprintf("  Loaded %d modeled weeks\n", N_model))

# ==============================================================================
# 2. EXTRACT POSTERIORS
# ==============================================================================

cat("\nExtracting posterior summaries...\n")

# Function to summarize Rt trajectory
summarize_rt <- function(fit, model_name) {
  Rt_draws <- fit$draws("Rt", format = "matrix")

  tibble(
    date = dates,
    model = model_name,
    Rt_mean = colMeans(Rt_draws),
    Rt_median = apply(Rt_draws, 2, median),
    Rt_lower = apply(Rt_draws, 2, quantile, 0.025),
    Rt_upper = apply(Rt_draws, 2, quantile, 0.975),
    Rt_lower50 = apply(Rt_draws, 2, quantile, 0.25),
    Rt_upper50 = apply(Rt_draws, 2, quantile, 0.75)
  )
}

rt_summary <- bind_rows(
  summarize_rt(fit0, "Model 0: Baseline"),
  summarize_rt(fit1, "Model 1: Climate"),
  summarize_rt(fit2, "Model 2: Full")
)

# ==============================================================================
# 3. PARAMETER SUMMARIES
# ==============================================================================

cat("\nGenerating parameter summaries...\n")

# Model 0 parameters
params0 <- fit0$summary(variables = c("mu", "alpha", "rho", "phi")) |>
  mutate(model = "Model 0")

# Model 1 parameters
params1 <- fit1$summary(variables = c("mu", "beta_climate", "alpha", "rho", "phi",
                                        "prop_climate", "residual_amplitude")) |>
  mutate(model = "Model 1")

# Model 2 parameters
params2 <- fit2$summary(variables = c("mu", "beta", "alpha", "rho", "phi",
                                        "temp_effect", "rain_effect",
                                        "wolbachia_effect", "npi_effect",
                                        "prop_climate", "prop_wolbachia",
                                        "prop_npi", "prop_residual")) |>
  mutate(model = "Model 2")

# Combine and save
all_params <- bind_rows(params0, params1, params2)
write_csv(all_params, "../results/parameter_summary.csv")
cat("  Saved: results/parameter_summary.csv\n")

# Print key results
cat("\n--- Key Parameter Estimates (Model 2) ---\n")

# Effect sizes
effects <- fit2$summary(variables = c("temp_effect", "rain_effect",
                                        "wolbachia_effect", "npi_effect"))
cat("\nMultiplicative effects on Rt:\n")
for (i in 1:nrow(effects)) {
  cat(sprintf("  %s: %.3f (95%% CI: %.3f - %.3f)\n",
              effects$variable[i],
              effects$mean[i],
              effects$q5[i],
              effects$q95[i]))
}

# Variance decomposition
decomp <- fit2$summary(variables = c("prop_climate", "prop_wolbachia",
                                       "prop_npi", "prop_residual"))
cat("\nVariance decomposition (proportion of log(Rt) variance):\n")
for (i in 1:nrow(decomp)) {
  cat(sprintf("  %s: %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
              decomp$variable[i],
              decomp$mean[i] * 100,
              decomp$q5[i] * 100,
              decomp$q95[i] * 100))
}

# ==============================================================================
# 4. RT TRAJECTORY PLOTS
# ==============================================================================

cat("\nGenerating Rt trajectory plots...\n")

# Observed cases for context
cases_df <- tibble(
  date = dates,
  cases = df$cases[(length(df$cases) - N_model + 1):length(df$cases)]
)

# Plot 1: All models comparison
p_rt_compare <- ggplot(rt_summary, aes(x = date)) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper, fill = model), alpha = 0.2) +
  geom_ribbon(aes(ymin = Rt_lower50, ymax = Rt_upper50, fill = model), alpha = 0.4) +
  geom_line(aes(y = Rt_median, color = model), linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  facet_wrap(~ model, ncol = 1) +
  scale_y_continuous(limits = c(0, 3)) +
  labs(
    title = "Estimated Effective Reproduction Number (Rt) - Model Comparison",
    subtitle = "Shaded regions: 50% and 95% credible intervals",
    x = "Date",
    y = expression(R[t]),
    color = "Model",
    fill = "Model"
  ) +
  theme(legend.position = "none")

ggsave("../results/figures/rt_comparison.png", p_rt_compare,
       width = 12, height = 10, dpi = 150)
cat("  Saved: results/figures/rt_comparison.png\n")

# Plot 2: Model 2 Rt with cases
p_rt_cases <- ggplot() +
  # Cases (secondary axis)
  geom_bar(data = cases_df, aes(x = date, y = cases / 500),
           stat = "identity", fill = "gray80", alpha = 0.7) +
  # Rt
  geom_ribbon(data = filter(rt_summary, model == "Model 2: Full"),
              aes(x = date, ymin = Rt_lower, ymax = Rt_upper),
              fill = "steelblue", alpha = 0.3) +
  geom_ribbon(data = filter(rt_summary, model == "Model 2: Full"),
              aes(x = date, ymin = Rt_lower50, ymax = Rt_upper50),
              fill = "steelblue", alpha = 0.5) +
  geom_line(data = filter(rt_summary, model == "Model 2: Full"),
            aes(x = date, y = Rt_median), color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(
    name = expression(R[t]),
    sec.axis = sec_axis(~ . * 500, name = "Weekly Cases")
  ) +
  labs(
    title = "Dengue Rt Estimates (Model 2: Full) with Weekly Cases",
    subtitle = "Singapore, 2012-2022",
    x = "Date"
  )

ggsave("../results/figures/rt_with_cases.png", p_rt_cases,
       width = 12, height = 6, dpi = 150)
cat("  Saved: results/figures/rt_with_cases.png\n")

# ==============================================================================
# 5. DECOMPOSITION PLOTS (MODEL 2)
# ==============================================================================

cat("\nGenerating decomposition plots...\n")

# Extract component posteriors
f_climate <- fit2$draws("f_climate", format = "matrix")
f_wolbachia <- fit2$draws("f_wolbachia", format = "matrix")
f_npi <- fit2$draws("f_npi", format = "matrix")
f_residual <- fit2$draws("f_residual", format = "matrix")

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
  summarize_component(f_wolbachia, "Wolbachia"),
  summarize_component(f_npi, "COVID-19 NPI"),
  summarize_component(f_residual, "Residual GP")
)

# Decomposition plot
p_decomp <- ggplot(components, aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = median), color = "steelblue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ component, ncol = 1, scales = "free_y") +
  labs(
    title = "Decomposition of log(Rt) into Components (Model 2)",
    subtitle = "Effect on log scale; positive = increases Rt, negative = decreases Rt",
    x = "Date",
    y = "Effect on log(Rt)"
  )

ggsave("../results/figures/decomposition.png", p_decomp,
       width = 12, height = 10, dpi = 150)
cat("  Saved: results/figures/decomposition.png\n")

# ==============================================================================
# 6. POSTERIOR PREDICTIVE CHECK
# ==============================================================================

cat("\nGenerating posterior predictive checks...\n")

# Extract predicted cases
cases_pred0 <- fit0$draws("cases_pred", format = "matrix")
cases_pred1 <- fit1$draws("cases_pred", format = "matrix")
cases_pred2 <- fit2$draws("cases_pred", format = "matrix")

# Observed cases
cases_obs <- df$cases[(length(df$cases) - N_model + 1):length(df$cases)]

# Summarize predictions
summarize_ppc <- function(draws, model_name) {
  tibble(
    date = dates,
    model = model_name,
    observed = cases_obs,
    pred_mean = colMeans(draws),
    pred_median = apply(draws, 2, median),
    pred_lower = apply(draws, 2, quantile, 0.025),
    pred_upper = apply(draws, 2, quantile, 0.975),
    pred_lower80 = apply(draws, 2, quantile, 0.1),
    pred_upper80 = apply(draws, 2, quantile, 0.9)
  )
}

ppc_summary <- bind_rows(
  summarize_ppc(cases_pred0, "Model 0"),
  summarize_ppc(cases_pred1, "Model 1"),
  summarize_ppc(cases_pred2, "Model 2")
)

# Compute coverage
compute_coverage <- function(df) {
  tibble(
    model = unique(df$model),
    coverage_95 = mean(df$observed >= df$pred_lower & df$observed <= df$pred_upper),
    coverage_80 = mean(df$observed >= df$pred_lower80 & df$observed <= df$pred_upper80)
  )
}

coverage <- ppc_summary |>
  group_by(model) |>
  group_modify(~ compute_coverage(.x)) |>
  ungroup()

cat("\nPosterior predictive coverage:\n")
print(coverage)

# PPC plot for Model 2
p_ppc <- ggplot(filter(ppc_summary, model == "Model 2"), aes(x = date)) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), fill = "steelblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = pred_lower80, ymax = pred_upper80), fill = "steelblue", alpha = 0.4) +
  geom_line(aes(y = pred_median), color = "steelblue", linewidth = 0.5) +
  geom_point(aes(y = observed), size = 0.5, alpha = 0.7) +
  labs(
    title = "Posterior Predictive Check (Model 2: Full)",
    subtitle = sprintf("80%% coverage: %.1f%%, 95%% coverage: %.1f%%",
                       filter(coverage, model == "Model 2")$coverage_80 * 100,
                       filter(coverage, model == "Model 2")$coverage_95 * 100),
    x = "Date",
    y = "Weekly Cases"
  )

ggsave("../results/figures/ppc_model2.png", p_ppc,
       width = 12, height = 5, dpi = 150)
cat("  Saved: results/figures/ppc_model2.png\n")

# ==============================================================================
# 7. GP HYPERPARAMETER COMPARISON
# ==============================================================================

cat("\nComparing GP hyperparameters across models...\n")

# Extract GP hyperparameters
gp_params <- bind_rows(
  fit0$summary(variables = c("alpha", "rho")) |> mutate(model = "Model 0"),
  fit1$summary(variables = c("alpha", "rho")) |> mutate(model = "Model 1"),
  fit2$summary(variables = c("alpha", "rho")) |> mutate(model = "Model 2")
)

cat("\nGP hyperparameter comparison:\n")
print(gp_params |> select(model, variable, mean, q5, q95))

# Does adding covariates reduce residual GP amplitude?
cat("\n--- Interpretation ---\n")
cat("If covariates explain Rt variation, the residual GP should have:\n")
cat("  - Smaller amplitude (alpha) in Models 1 and 2\n")
cat("  - Possibly shorter length scale (rho) if covariates absorb slow trends\n")

# ==============================================================================
# 8. INTERVENTION EFFECT VISUALIZATION
# ==============================================================================

cat("\nGenerating intervention effect plots...\n")

# Wolbachia coverage trajectory
wolbachia_df <- tibble(
  date = dates,
  coverage = df$wolbachia_coverage[(length(df$wolbachia_coverage) - N_model + 1):length(df$wolbachia_coverage)]
)

# NPI trajectory
npi_df <- tibble(
  date = dates,
  npi = df$npi_intensity[(length(df$npi_intensity) - N_model + 1):length(df$npi_intensity)]
)

# Extract wolbachia effect posterior
wolbachia_effect <- fit2$draws("wolbachia_effect", format = "draws_array")
npi_effect <- fit2$draws("npi_effect", format = "draws_array")

# Effect size plot
p_effects <- mcmc_areas(fit2$draws(c("temp_effect", "rain_effect",
                                       "wolbachia_effect", "npi_effect")),
                         prob = 0.8, prob_outer = 0.95) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Covariate Effects on Rt (Multiplicative Scale)",
    subtitle = "Values < 1 indicate Rt reduction; > 1 indicate Rt increase",
    x = "Effect on Rt"
  )

ggsave("../results/figures/effect_sizes.png", p_effects,
       width = 8, height = 5, dpi = 150)
cat("  Saved: results/figures/effect_sizes.png\n")

# ==============================================================================
# 9. SAVE ALL SUMMARIES
# ==============================================================================

cat("\nSaving summary data...\n")

# Save Rt trajectories
write_csv(rt_summary, "../results/rt_trajectories.csv")
cat("  Saved: results/rt_trajectories.csv\n")

# Save component decomposition
write_csv(components, "../results/decomposition.csv")
cat("  Saved: results/decomposition.csv\n")

# Save PPC summaries
write_csv(ppc_summary, "../results/ppc_summary.csv")
cat("  Saved: results/ppc_summary.csv\n")

# Save coverage statistics
write_csv(coverage, "../results/ppc_coverage.csv")
cat("  Saved: results/ppc_coverage.csv\n")

# ==============================================================================
# 10. FINAL SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("POST-PROCESSING COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\n--- KEY FINDINGS ---\n\n")

# Model comparison
loo_results <- readRDS("../results/loo_comparison.rds")
cat("Model Comparison (LOO-CV):\n")
print(loo_results$comparison)

cat("\n\nInterpretation:\n")
cat("  - Positive elpd_diff indicates better predictive performance\n")
cat("  - If Model 2 outperforms Model 0, covariates improve prediction\n")
cat("  - Examine variance decomposition to understand relative importance\n")

# Wolbachia effect
wolb_eff <- fit2$summary("wolbachia_effect")
cat(sprintf("\n\nWolbachia Effect:\n"))
cat(sprintf("  At full coverage, Rt multiplied by: %.3f (95%% CI: %.3f - %.3f)\n",
            wolb_eff$mean, wolb_eff$q5, wolb_eff$q95))
cat(sprintf("  Percentage reduction: %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
            (1 - wolb_eff$mean) * 100,
            (1 - wolb_eff$q95) * 100,
            (1 - wolb_eff$q5) * 100))

cat("\n\nOutput files:\n")
cat("  results/figures/rt_comparison.png - Rt estimates across models\n")
cat("  results/figures/rt_with_cases.png - Rt with observed cases\n")
cat("  results/figures/decomposition.png - Component decomposition\n")
cat("  results/figures/ppc_model2.png - Posterior predictive check\n")
cat("  results/figures/effect_sizes.png - Covariate effect sizes\n")
cat("  results/*.csv - Summary tables\n")
