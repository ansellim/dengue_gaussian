#!/usr/bin/env Rscript
# ==============================================================================
# 41_fit_and_compare.R
#
# Fit Model B (climate + susceptible fraction + GP) and compare against
# Model A (climate + GP only) to test whether immunity dynamics reduce
# the GP residual's role.
#
# Model A: log(Rt) = mu + f_climate + f_GP                (existing)
# Model B: log(Rt) = mu + f_climate + beta_S * log_S_eff + f_GP  (new)
#
# Input:
#   data/model_data.rds
#   data/susceptible_covariate.rds (from 41_prepare_susceptible_covariate.R)
#   code/41_model_with_susceptible.stan
#   results/fit_model3.rds (Model A, existing)
#
# Output:
#   results/fit_model_susceptible.rds
#   results/figures/susceptible_model_comparison.png
#   results/figures/susceptible_variance_decomposition.png
#   results/figures/susceptible_gp_amplitude_comparison.png
#   results/figures/susceptible_beta_posterior.png
#   results/susceptible_comparison_summary.csv
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
library(patchwork)
library(lubridate)

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
dir.create("../results/figures", showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank()))

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("FIT MODEL B (CLIMATE + SUSCEPTIBLE + GP) AND COMPARE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")

model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data
dates_model <- model_data$metadata$dates_model

susceptible <- readRDS("../data/susceptible_covariate.rds")

cat(sprintf("  N_model = %d weeks\n", stan_data$N_model))
cat(sprintf("  log_S_eff_std length = %d\n", length(susceptible$log_S_eff_std)))

# Add susceptible covariate to Stan data
stan_data_B <- stan_data
stan_data_B$log_S_eff <- susceptible$log_S_eff_std

# ==============================================================================
# 2. COMPILE AND FIT MODEL B
# ==============================================================================

cat("\nCompiling Model B (climate + susceptible + GP)...\n")
cat("  log(Rt) = mu + beta_temp*temp + beta_rain*rain + beta_S*log_S_eff + f_GP\n\n")

model_B <- cmdstan_model("41_model_with_susceptible.stan")

cat("Fitting Model B...\n")
fit_B <- model_B$sample(
  data = stan_data_B,
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
fit_B$save_object("../results/fit_model_susceptible.rds")
cat("\nSaved: results/fit_model_susceptible.rds\n")

# ==============================================================================
# 3. LOAD MODEL A (EXISTING)
# ==============================================================================

cat("\nLoading Model A (climate-only) fit...\n")

if (!file.exists("../results/fit_model3.rds")) {
  stop("Model A fit not found. Run 13_fit_model3.R first.")
}

fit_A <- readRDS("../results/fit_model3.rds")
cat("  Loaded: results/fit_model3.rds\n")

# ==============================================================================
# 4. DIAGNOSTICS FOR MODEL B
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MODEL B DIAGNOSTICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

diag_B <- fit_B$diagnostic_summary()
cat(sprintf("Divergent transitions: %d\n", sum(diag_B$num_divergent)))
cat(sprintf("Max treedepth exceeded: %d\n", sum(diag_B$num_max_treedepth)))

# ==============================================================================
# 5. PARAMETER SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MODEL B PARAMETER SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

params_B <- c("mu", "alpha", "rho", "phi", "beta_S",
              "temp_effect", "rain_effect", "susceptible_effect",
              "prop_climate", "prop_susceptible", "prop_residual")
summ_B <- fit_B$summary(variables = params_B)
print(summ_B)

# ==============================================================================
# 6. LOO-CV COMPARISON
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("LOO-CV COMPARISON: MODEL A vs MODEL B\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

log_lik_A <- fit_A$draws("log_lik", format = "matrix")
log_lik_B <- fit_B$draws("log_lik", format = "matrix")

loo_A <- loo(log_lik_A)
loo_B <- loo(log_lik_B)

cat("Model A (climate + GP):\n")
print(loo_A)
cat("\nModel B (climate + susceptible + GP):\n")
print(loo_B)

loo_comp <- loo_compare(list("Model_A_climate_GP" = loo_A,
                              "Model_B_climate_susc_GP" = loo_B))
cat("\nLOO comparison (positive = better):\n")
print(loo_comp)

# ==============================================================================
# 7. GP AMPLITUDE COMPARISON
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("GP AMPLITUDE COMPARISON\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

alpha_A <- fit_A$draws("alpha", format = "matrix")[, 1]
alpha_B <- fit_B$draws("alpha", format = "matrix")[, 1]

cat(sprintf("  Model A alpha: median = %.4f [%.4f, %.4f]\n",
            median(alpha_A), quantile(alpha_A, 0.025), quantile(alpha_A, 0.975)))
cat(sprintf("  Model B alpha: median = %.4f [%.4f, %.4f]\n",
            median(alpha_B), quantile(alpha_B, 0.025), quantile(alpha_B, 0.975)))
cat(sprintf("  Ratio (B/A): %.3f\n", median(alpha_B) / median(alpha_A)))

rho_A <- fit_A$draws("rho", format = "matrix")[, 1]
rho_B <- fit_B$draws("rho", format = "matrix")[, 1]

cat(sprintf("\n  Model A rho: median = %.2f [%.2f, %.2f]\n",
            median(rho_A), quantile(rho_A, 0.025), quantile(rho_A, 0.975)))
cat(sprintf("  Model B rho: median = %.2f [%.2f, %.2f]\n",
            median(rho_B), quantile(rho_B, 0.025), quantile(rho_B, 0.975)))

# ==============================================================================
# 8. BETA_S POSTERIOR
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUSCEPTIBLE FRACTION EFFECT (beta_S)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

beta_S <- fit_B$draws("beta_S", format = "matrix")[, 1]
susc_effect <- fit_B$draws("susceptible_effect", format = "matrix")[, 1]

cat(sprintf("  beta_S: median = %.4f [%.4f, %.4f]\n",
            median(beta_S), quantile(beta_S, 0.025), quantile(beta_S, 0.975)))
cat(sprintf("  P(beta_S > 0) = %.3f\n", mean(beta_S > 0)))
cat(sprintf("  Multiplicative effect (exp(beta_S)): %.3f [%.3f, %.3f]\n",
            median(susc_effect),
            quantile(susc_effect, 0.025),
            quantile(susc_effect, 0.975)))

# ==============================================================================
# 9. VARIANCE DECOMPOSITION COMPARISON
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("VARIANCE DECOMPOSITION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Model A
prop_clim_A <- fit_A$draws("prop_climate", format = "matrix")[, 1]
prop_resid_A <- fit_A$draws("prop_residual", format = "matrix")[, 1]

cat("Model A (climate + GP):\n")
cat(sprintf("  Climate:  %.2f%% [%.2f%%, %.2f%%]\n",
            median(prop_clim_A) * 100,
            quantile(prop_clim_A, 0.025) * 100,
            quantile(prop_clim_A, 0.975) * 100))
cat(sprintf("  Residual: %.2f%% [%.2f%%, %.2f%%]\n",
            median(prop_resid_A) * 100,
            quantile(prop_resid_A, 0.025) * 100,
            quantile(prop_resid_A, 0.975) * 100))

# Model B
prop_clim_B <- fit_B$draws("prop_climate", format = "matrix")[, 1]
prop_susc_B <- fit_B$draws("prop_susceptible", format = "matrix")[, 1]
prop_resid_B <- fit_B$draws("prop_residual", format = "matrix")[, 1]

cat("\nModel B (climate + susceptible + GP):\n")
cat(sprintf("  Climate:      %.2f%% [%.2f%%, %.2f%%]\n",
            median(prop_clim_B) * 100,
            quantile(prop_clim_B, 0.025) * 100,
            quantile(prop_clim_B, 0.975) * 100))
cat(sprintf("  Susceptible:  %.2f%% [%.2f%%, %.2f%%]\n",
            median(prop_susc_B) * 100,
            quantile(prop_susc_B, 0.025) * 100,
            quantile(prop_susc_B, 0.975) * 100))
cat(sprintf("  Residual:     %.2f%% [%.2f%%, %.2f%%]\n",
            median(prop_resid_B) * 100,
            quantile(prop_resid_B, 0.025) * 100,
            quantile(prop_resid_B, 0.975) * 100))

# ==============================================================================
# 10. FIGURES
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("GENERATING FIGURES\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Serotype switch dates (approximate, for reference lines)
switch_dates <- as.Date(c("2013-07-01", "2016-01-01", "2020-01-01"))

# ---------- Figure 1: Model comparison (3-panel) ----------

cat("  Figure 1: Model comparison...\n")

# Extract Rt and f_residual posteriors
Rt_A <- fit_A$draws("Rt", format = "matrix")
Rt_B <- fit_B$draws("Rt", format = "matrix")
f_resid_A <- fit_A$draws("f_residual", format = "matrix")
f_resid_B <- fit_B$draws("f_residual", format = "matrix")

N_model <- length(dates_model)

comparison_df <- tibble(
  date = dates_model,
  # Model A Rt
  Rt_A_median = apply(Rt_A, 2, median),
  Rt_A_lo = apply(Rt_A, 2, quantile, 0.025),
  Rt_A_hi = apply(Rt_A, 2, quantile, 0.975),
  # Model B Rt
  Rt_B_median = apply(Rt_B, 2, median),
  Rt_B_lo = apply(Rt_B, 2, quantile, 0.025),
  Rt_B_hi = apply(Rt_B, 2, quantile, 0.975),
  # Model A f_residual
  fA_median = apply(f_resid_A, 2, median),
  fA_lo = apply(f_resid_A, 2, quantile, 0.025),
  fA_hi = apply(f_resid_A, 2, quantile, 0.975),
  # Model B f_residual
  fB_median = apply(f_resid_B, 2, median),
  fB_lo = apply(f_resid_B, 2, quantile, 0.025),
  fB_hi = apply(f_resid_B, 2, quantile, 0.975),
  # S_eff
  S_eff = susceptible$S_eff_raw
)

p1_top <- ggplot(comparison_df, aes(x = date)) +
  geom_ribbon(aes(ymin = Rt_A_lo, ymax = Rt_A_hi), fill = "steelblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = Rt_B_lo, ymax = Rt_B_hi), fill = "firebrick", alpha = 0.2) +
  geom_line(aes(y = Rt_A_median, color = "Model A (climate + GP)"), linewidth = 0.6) +
  geom_line(aes(y = Rt_B_median, color = "Model B (+ susceptible)"), linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = switch_dates, linetype = "dotted", color = "grey50", alpha = 0.6) +
  scale_color_manual(values = c("Model A (climate + GP)" = "steelblue",
                                 "Model B (+ susceptible)" = "firebrick")) +
  labs(y = expression(R[t]), color = NULL,
       title = "Reproduction number: Model A vs Model B") +
  theme(legend.position = "top")

p1_mid <- ggplot(comparison_df, aes(x = date)) +
  geom_ribbon(aes(ymin = fA_lo, ymax = fA_hi), fill = "steelblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = fB_lo, ymax = fB_hi), fill = "firebrick", alpha = 0.2) +
  geom_line(aes(y = fA_median, color = "Model A"), linewidth = 0.6) +
  geom_line(aes(y = fB_median, color = "Model B"), linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = switch_dates, linetype = "dotted", color = "grey50", alpha = 0.6) +
  scale_color_manual(values = c("Model A" = "steelblue", "Model B" = "firebrick")) +
  labs(y = "f_residual (GP)", color = NULL,
       title = "GP residual: has it shrunk in Model B?") +
  theme(legend.position = "top")

p1_bot <- ggplot(comparison_df, aes(x = date, y = S_eff)) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  geom_vline(xintercept = switch_dates, linetype = "dotted", color = "grey50", alpha = 0.6) +
  labs(y = expression(S[eff]), title = "Effective susceptible fraction") +
  theme()

p1 <- p1_top / p1_mid / p1_bot +
  plot_annotation(tag_levels = "A")

ggsave("../results/figures/susceptible_model_comparison.png", p1,
       width = 12, height = 10, dpi = 300, bg = "white")
cat("    Saved: susceptible_model_comparison.png\n")

# ---------- Figure 2: Variance decomposition ----------

cat("  Figure 2: Variance decomposition...\n")

var_df <- tibble(
  model = rep(c("Model A\n(climate + GP)",
                "Model B\n(climate + susc + GP)"), each = 3),
  component = rep(c("Climate", "Susceptible", "Residual GP"), 2),
  median = c(
    median(prop_clim_A) * 100, 0, median(prop_resid_A) * 100,
    median(prop_clim_B) * 100, median(prop_susc_B) * 100, median(prop_resid_B) * 100
  ),
  lo = c(
    quantile(prop_clim_A, 0.025) * 100, 0, quantile(prop_resid_A, 0.025) * 100,
    quantile(prop_clim_B, 0.025) * 100, quantile(prop_susc_B, 0.025) * 100,
    quantile(prop_resid_B, 0.025) * 100
  ),
  hi = c(
    quantile(prop_clim_A, 0.975) * 100, 0, quantile(prop_resid_A, 0.975) * 100,
    quantile(prop_clim_B, 0.975) * 100, quantile(prop_susc_B, 0.975) * 100,
    quantile(prop_resid_B, 0.975) * 100
  )
)

var_df$component <- factor(var_df$component,
                            levels = c("Climate", "Susceptible", "Residual GP"))

p2 <- ggplot(var_df, aes(x = component, y = median, fill = component)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", median)),
            vjust = -0.5, size = 3.5) +
  facet_wrap(~model) +
  scale_fill_manual(values = c("Climate" = "steelblue",
                                "Susceptible" = "darkgreen",
                                "Residual GP" = "grey60")) +
  labs(y = "Variance explained (%)",
       title = "Variance decomposition: did susceptible fraction reduce GP dominance?") +
  ylim(0, 105) +
  theme(legend.position = "none")

ggsave("../results/figures/susceptible_variance_decomposition.png", p2,
       width = 10, height = 5, dpi = 300, bg = "white")
cat("    Saved: susceptible_variance_decomposition.png\n")

# ---------- Figure 3: GP amplitude comparison ----------

cat("  Figure 3: GP amplitude comparison...\n")

alpha_df <- tibble(
  alpha = c(alpha_A, alpha_B),
  model = rep(c("Model A (climate + GP)", "Model B (+ susceptible)"),
              c(length(alpha_A), length(alpha_B)))
)

p3 <- ggplot(alpha_df, aes(x = alpha, fill = model)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = median(alpha_A), linetype = "dashed", color = "steelblue") +
  geom_vline(xintercept = median(alpha_B), linetype = "dashed", color = "firebrick") +
  scale_fill_manual(values = c("Model A (climate + GP)" = "steelblue",
                                "Model B (+ susceptible)" = "firebrick")) +
  labs(x = expression(alpha ~ "(GP amplitude)"),
       y = "Posterior density",
       fill = NULL,
       title = "GP amplitude: lower in Model B = susceptible absorbed GP variance") +
  theme(legend.position = "top")

ggsave("../results/figures/susceptible_gp_amplitude_comparison.png", p3,
       width = 8, height = 5, dpi = 300, bg = "white")
cat("    Saved: susceptible_gp_amplitude_comparison.png\n")

# ---------- Figure 4: beta_S posterior ----------

cat("  Figure 4: beta_S posterior...\n")

beta_df <- tibble(beta_S = beta_S)

p4 <- ggplot(beta_df, aes(x = beta_S)) +
  geom_density(fill = "darkgreen", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = median(beta_S), linetype = "solid", color = "darkgreen") +
  annotate("text",
           x = median(beta_S), y = 0,
           label = sprintf("median = %.3f\n95%% CrI: [%.3f, %.3f]\nP(>0) = %.3f",
                           median(beta_S),
                           quantile(beta_S, 0.025),
                           quantile(beta_S, 0.975),
                           mean(beta_S > 0)),
           hjust = -0.1, vjust = -0.5, size = 3.5) +
  labs(x = expression(beta[S] ~ "(effect of log susceptible fraction on log " * R[t] * ")"),
       y = "Posterior density",
       title = expression("Posterior of " * beta[S] * ": susceptible fraction effect")) +
  theme()

ggsave("../results/figures/susceptible_beta_posterior.png", p4,
       width = 8, height = 5, dpi = 300, bg = "white")
cat("    Saved: susceptible_beta_posterior.png\n")

# ==============================================================================
# 11. SUMMARY TABLE
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Build summary
summary_df <- tibble(
  metric = c(
    "LOO_elpd_A", "LOO_elpd_B", "LOO_elpd_diff", "LOO_se_diff",
    "alpha_A_median", "alpha_B_median", "alpha_ratio_B_over_A",
    "rho_A_median", "rho_B_median",
    "beta_S_median", "beta_S_lo95", "beta_S_hi95", "P_beta_S_positive",
    "susceptible_effect_median",
    "prop_climate_A", "prop_residual_A",
    "prop_climate_B", "prop_susceptible_B", "prop_residual_B"
  ),
  value = c(
    loo_A$estimates["elpd_loo", "Estimate"],
    loo_B$estimates["elpd_loo", "Estimate"],
    loo_comp[2, "elpd_diff"],
    loo_comp[2, "se_diff"],
    median(alpha_A), median(alpha_B),
    median(alpha_B) / median(alpha_A),
    median(rho_A), median(rho_B),
    median(beta_S), quantile(beta_S, 0.025), quantile(beta_S, 0.975),
    mean(beta_S > 0),
    median(susc_effect),
    median(prop_clim_A) * 100, median(prop_resid_A) * 100,
    median(prop_clim_B) * 100, median(prop_susc_B) * 100, median(prop_resid_B) * 100
  )
)

print(summary_df, n = Inf)

write_csv(summary_df, "../results/susceptible_comparison_summary.csv")
cat("\nSaved: results/susceptible_comparison_summary.csv\n")

# Print interpretation
cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("INTERPRETATION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

elpd_diff <- loo_comp[2, "elpd_diff"]
se_diff <- loo_comp[2, "se_diff"]
alpha_ratio <- median(alpha_B) / median(alpha_A)

cat(sprintf("LOO comparison: elpd_diff = %.1f (SE = %.1f)\n", elpd_diff, se_diff))
if (abs(elpd_diff) < 2 * se_diff) {
  cat("  -> Models are indistinguishable by LOO (difference < 2 SE)\n")
} else if (elpd_diff > 0) {
  cat("  -> Model B preferred (adding S_eff improves prediction)\n")
} else {
  cat("  -> Model A preferred (S_eff does not help prediction)\n")
}

cat(sprintf("\nGP amplitude ratio (B/A): %.3f\n", alpha_ratio))
if (alpha_ratio < 0.9) {
  cat("  -> GP amplitude shrunk >10%%: S_eff absorbed some GP variance\n")
} else {
  cat("  -> GP amplitude did not shrink much: S_eff may be redundant with GP\n")
}

cat(sprintf("\nbeta_S: %.3f [%.3f, %.3f], P(>0) = %.3f\n",
            median(beta_S), quantile(beta_S, 0.025), quantile(beta_S, 0.975),
            mean(beta_S > 0)))
if (quantile(beta_S, 0.025) > 0) {
  cat("  -> Significantly positive: more susceptibles -> higher Rt\n")
} else if (mean(beta_S > 0) > 0.9) {
  cat("  -> Positive but 95%% CrI includes zero: suggestive but not conclusive\n")
} else {
  cat("  -> Not clearly positive: susceptible fraction may not drive Rt\n")
}

cat(sprintf("\nVariance explained by susceptible fraction: %.1f%%\n",
            median(prop_susc_B) * 100))

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Files saved:\n")
cat("  results/fit_model_susceptible.rds\n")
cat("  results/figures/susceptible_model_comparison.png\n")
cat("  results/figures/susceptible_variance_decomposition.png\n")
cat("  results/figures/susceptible_gp_amplitude_comparison.png\n")
cat("  results/figures/susceptible_beta_posterior.png\n")
cat("  results/susceptible_comparison_summary.csv\n")
