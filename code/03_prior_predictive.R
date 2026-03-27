#!/usr/bin/env Rscript
# ==============================================================================
# 03_prior_predictive.R
#
# Prior predictive checks for dengue Rt estimation model
# Verifies that priors generate plausible case counts and Rt trajectories
#
# Input: data/model_data.rds
# Output: results/figures/prior_predictive_*.png
# ==============================================================================

library(tidyverse)
library(posterior)
library(bayesplot)
library(patchwork)

# Set up paths
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # Get script directory when run from command line
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- grep("--file=", args, value = TRUE)
  if (length(script_path) > 0) {
    setwd(dirname(normalizePath(sub("--file=", "", script_path))))
  }
}

# Create directories
dir.create("../results", showWarnings = FALSE)
dir.create("../results/figures", showWarnings = FALSE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PRIOR PREDICTIVE CHECKS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data
metadata <- model_data$metadata

observed_cases <- stan_data$cases[(stan_data$S + 1):stan_data$N]
N_model <- stan_data$N_model

cat(sprintf("  N_model = %d weeks\n", N_model))
cat(sprintf("  Observed cases: min=%.0f, median=%.0f, max=%.0f\n",
            min(observed_cases), median(observed_cases), max(observed_cases)))

# ==============================================================================
# 2. SIMULATE FROM PRIORS
# ==============================================================================

cat("\nSimulating from priors...\n")

set.seed(42)
n_sims <- 500  # Number of prior simulations

# Storage for prior predictive samples
prior_Rt <- matrix(NA, nrow = n_sims, ncol = N_model)
prior_cases <- matrix(NA, nrow = n_sims, ncol = N_model)

# Prior specifications (from 05_model3_climate_only.stan)
# mu ~ Normal(0, 0.5)
# beta_climate[1] (temp) ~ Normal(0, 0.5)
# beta_climate[2] (rain) ~ Normal(0, 0.5)
# log_alpha ~ Normal(-1.2, 0.5)
# log_rho ~ Normal(log(6), 0.5)
# phi ~ HalfNormal(0, 10) -- we'll use truncated normal

cat("  Prior specifications:\n")
cat("    mu ~ Normal(0, 0.5)\n")
cat("    beta_temp ~ Normal(0, 0.5)\n")
cat("    beta_rain ~ Normal(0, 0.5)\n")
cat("    log_alpha ~ Normal(-1.2, 0.5)\n")
cat("    log_rho ~ Normal(log(6), 0.5)\n")
cat("    phi ~ HalfNormal(0, 10)\n\n")

# Extract covariates
X_climate <- stan_data$X_climate
gi <- stan_data$gi
S <- stan_data$S
M <- stan_data$M
L <- stan_data$L
t_vec <- stan_data$t

# Precompute HSGP basis matrix (same computation as in Stan)
sqrt_eigenvalues <- (1:M) * pi / (2 * L)
PHI <- matrix(0, N_model, M)
for (ii in 1:N_model) {
  for (jj in 1:M) {
    PHI[ii, jj] <- sin(sqrt_eigenvalues[jj] * (t_vec[ii] + L)) / sqrt(L)
  }
}

# Matern 3/2 spectral density
spd_matern32 <- function(omega, alpha, rho) {
  s3r <- sqrt(3) / rho
  alpha^2 * 4 * s3r^3 / (s3r^2 + omega^2)^2
}

for (i in 1:n_sims) {
  # Sample from priors
  mu <- rnorm(1, 0, 0.5)
  beta_climate <- c(
    rnorm(1, 0, 0.5),      # temp
    rnorm(1, 0, 0.5)       # rain
  )
  log_alpha <- rnorm(1, -1.2, 0.5)
  log_rho <- rnorm(1, log(6), 0.5)
  phi <- abs(rnorm(1, 0, 10))  # Half-normal

  alpha <- exp(log_alpha)
  rho <- exp(log_rho)

  # Compute covariate effects (climate only)
  f_climate <- X_climate[, 1] * beta_climate[1] + X_climate[, 2] * beta_climate[2]

  # GP realization via HSGP spectral decomposition
  spd_weights <- sqrt(sapply(sqrt_eigenvalues, spd_matern32,
                              alpha = alpha, rho = rho))
  beta_gp <- rnorm(M)
  f_residual <- as.numeric(PHI %*% (spd_weights * beta_gp))

  # log(Rt)
  log_Rt <- mu + f_climate + f_residual
  prior_Rt[i, ] <- exp(log_Rt)

  # Simulate cases using renewal equation
  cases_sim <- stan_data$cases  # Start with observed for infectious pressure

  for (t in 1:N_model) {
    t_idx <- S + t
    infectious_pressure <- sum(cases_sim[(t_idx - 1):(t_idx - S)] * rev(gi))
    lambda <- max(exp(log_Rt[t]) * infectious_pressure, 1)

    # Negative binomial draw
    if (is.finite(lambda) && lambda > 0 && is.finite(phi) && phi > 0) {
      cases_sim[t_idx] <- rnbinom(1, size = phi, mu = lambda)
    } else {
      cases_sim[t_idx] <- NA
    }
  }

  prior_cases[i, ] <- cases_sim[(S + 1):stan_data$N]
}

cat(sprintf("  Completed %d simulations\n", n_sims))

# ==============================================================================
# 3. PRIOR PREDICTIVE CHECKS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PRIOR PREDICTIVE SUMMARIES\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Check Rt ranges
Rt_median <- apply(prior_Rt, 1, median, na.rm = TRUE)
Rt_max <- apply(prior_Rt, 1, max, na.rm = TRUE)
Rt_min <- apply(prior_Rt, 1, min, na.rm = TRUE)

cat("Rt from prior:\n")
cat(sprintf("  Median Rt (across simulations): %.2f (95%% interval: %.2f - %.2f)\n",
            median(Rt_median, na.rm = TRUE),
            quantile(Rt_median, 0.025, na.rm = TRUE),
            quantile(Rt_median, 0.975, na.rm = TRUE)))
cat(sprintf("  Max Rt observed in any simulation: %.2f\n", max(Rt_max, na.rm = TRUE)))
cat(sprintf("  Min Rt observed in any simulation: %.4f\n", min(Rt_min, na.rm = TRUE)))

# Check case ranges
cases_median <- apply(prior_cases, 1, median, na.rm = TRUE)
cases_max <- apply(prior_cases, 1, max, na.rm = TRUE)

cat("\nCases from prior:\n")
cat(sprintf("  Median cases (across simulations): %.0f (95%% interval: %.0f - %.0f)\n",
            median(cases_median, na.rm = TRUE),
            quantile(cases_median, 0.025, na.rm = TRUE),
            quantile(cases_median, 0.975, na.rm = TRUE)))
cat(sprintf("  Max cases in any simulation: %s\n",
            format(max(cases_max, na.rm = TRUE), big.mark = ",")))

# Proportion of simulations with extreme values
prop_Rt_over_10 <- mean(Rt_max > 10, na.rm = TRUE)
prop_cases_over_10k <- mean(cases_max > 10000, na.rm = TRUE)

cat("\nExtreme value checks:\n")
cat(sprintf("  Simulations with Rt > 10: %.1f%%\n", 100 * prop_Rt_over_10))
cat(sprintf("  Simulations with cases > 10,000/week: %.1f%%\n", 100 * prop_cases_over_10k))

# ==============================================================================
# 4. VISUALIZATIONS
# ==============================================================================

cat("\nGenerating visualizations...\n")

# Theme
theme_set(theme_minimal(base_size = 12))

# --- Plot 1: Prior Rt trajectories ---
# Sample 50 trajectories for visualization
set.seed(123)
sample_idx <- sample(1:n_sims, min(50, n_sims))

Rt_df <- as_tibble(t(prior_Rt[sample_idx, ])) |>
  mutate(week = 1:N_model) |>
  pivot_longer(-week, names_to = "sim", values_to = "Rt")

p_Rt_prior <- ggplot(Rt_df, aes(x = week, y = Rt, group = sim)) +
  geom_line(alpha = 0.2, color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0, 5)) +
  labs(
    title = "Prior Predictive: Rt Trajectories",
    subtitle = sprintf("%d random draws from prior (clipped at Rt=5)", length(sample_idx)),
    x = "Week",
    y = "Rt"
  )

# --- Plot 2: Prior case distribution vs observed ---
prior_cases_flat <- as.vector(prior_cases)
prior_cases_flat <- prior_cases_flat[!is.na(prior_cases_flat) & prior_cases_flat < 1e6]

cases_df <- tibble(
  cases = c(prior_cases_flat, observed_cases),
  type = c(rep("Prior Predictive", length(prior_cases_flat)),
           rep("Observed", length(observed_cases)))
)

p_cases_dist <- ggplot(cases_df, aes(x = cases, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 50, alpha = 0.6, position = "identity") +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values = c("Prior Predictive" = "steelblue",
                               "Observed" = "darkred")) +
  labs(
    title = "Prior Predictive: Case Distribution",
    subtitle = "Log scale; comparing prior samples to observed data",
    x = "Weekly Cases",
    y = "Density",
    fill = ""
  ) +
  theme(legend.position = "bottom")

# --- Plot 3: Prior predictive intervals for Rt ---
Rt_summary <- tibble(
  week = 1:N_model,
  median = apply(prior_Rt, 2, median, na.rm = TRUE),
  q05 = apply(prior_Rt, 2, quantile, 0.05, na.rm = TRUE),
  q25 = apply(prior_Rt, 2, quantile, 0.25, na.rm = TRUE),
  q75 = apply(prior_Rt, 2, quantile, 0.75, na.rm = TRUE),
  q95 = apply(prior_Rt, 2, quantile, 0.95, na.rm = TRUE)
)

p_Rt_intervals <- ggplot(Rt_summary, aes(x = week)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "steelblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "steelblue", alpha = 0.4) +
  geom_line(aes(y = median), color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0, 4)) +
  labs(
    title = "Prior Predictive: Rt Intervals",
    subtitle = "Median with 50% and 90% intervals",
    x = "Week",
    y = "Rt"
  )

# --- Plot 4: Prior distributions for key parameters ---
param_samples <- tibble(
  mu = rnorm(4000, 0, 0.5),
  alpha = exp(rnorm(4000, -1.2, 0.5)),
  rho = exp(rnorm(4000, log(6), 0.5))
) |>
  mutate(
    baseline_Rt = exp(mu)
  )

p_prior_baseline <- ggplot(param_samples, aes(x = baseline_Rt)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  coord_cartesian(xlim = c(0, 5)) +
  labs(
    title = "Prior: Baseline Rt",
    subtitle = "Normal(0, 0.5) on log scale",
    x = "Baseline Rt",
    y = "Density"
  )

p_prior_gp <- ggplot(param_samples, aes(x = alpha)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "steelblue", alpha = 0.7) +
  coord_cartesian(xlim = c(0, 2)) +
  labs(
    title = "Prior: GP Amplitude",
    subtitle = "LogNormal(-1.2, 0.5)",
    x = "alpha (marginal SD)",
    y = "Density"
  )

# --- Combine and save ---
p_combined <- (p_Rt_prior | p_Rt_intervals) / (p_cases_dist | p_prior_baseline)

ggsave("../results/figures/prior_predictive_summary.png", p_combined,
       width = 14, height = 10, dpi = 150)
cat("  Saved: results/figures/prior_predictive_summary.png\n")

# Individual plots
ggsave("../results/figures/prior_predictive_Rt_trajectories.png", p_Rt_prior,
       width = 10, height = 6, dpi = 150)
ggsave("../results/figures/prior_predictive_cases.png", p_cases_dist,
       width = 8, height = 6, dpi = 150)

# ==============================================================================
# 5. ASSESSMENT
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PRIOR PREDICTIVE ASSESSMENT\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Check criteria
checks <- list(
  Rt_reasonable = median(Rt_median, na.rm = TRUE) > 0.5 &&
                  median(Rt_median, na.rm = TRUE) < 3,
  Rt_not_extreme = prop_Rt_over_10 < 0.2,
  cases_overlap = any(prior_cases_flat > min(observed_cases) &
                      prior_cases_flat < max(observed_cases)),
  cases_not_extreme = prop_cases_over_10k < 0.3
)

cat("Checks:\n")
cat(sprintf("  [%s] Median Rt in reasonable range (0.5-3)\n",
            ifelse(checks$Rt_reasonable, "PASS", "WARN")))
cat(sprintf("  [%s] <20%% simulations with Rt > 10\n",
            ifelse(checks$Rt_not_extreme, "PASS", "WARN")))
cat(sprintf("  [%s] Prior cases overlap with observed range\n",
            ifelse(checks$cases_overlap, "PASS", "WARN")))
cat(sprintf("  [%s] <30%% simulations with extreme cases (>10k/week)\n",
            ifelse(checks$cases_not_extreme, "PASS", "WARN")))

all_pass <- all(unlist(checks))

cat(sprintf("\nOverall: %s\n",
            ifelse(all_pass, "PRIORS APPEAR REASONABLE",
                   "REVIEW PRIORS - some checks failed")))

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PRIOR PREDICTIVE CHECKS COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("\nNext step: Fit the model, then run 09_posterior_predictive.R\n")
