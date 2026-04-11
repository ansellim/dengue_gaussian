#!/usr/bin/env Rscript
# ==============================================================================
# 10_fit_short_gp.R
#
# Fit model with shorter GP length scale to reduce residual autocorrelation
# Compare with standard model
#
# Key change: rho prior LogNormal(log(6), 0.5) vs LogNormal(log(15), 0.5)
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

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

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("FITTING MODEL WITH SHORTER GP LENGTH SCALE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data <- model_data$stan_data

cat(sprintf("  N_model = %d weeks\n", stan_data$N_model))
cat(sprintf("  Cases range: %d - %d\n", min(stan_data$cases), max(stan_data$cases)))

# ==============================================================================
# 2. COMPILE AND FIT MODEL
# ==============================================================================

cat("\nCompiling model with shorter GP length scale...\n")
cat("  Key change: rho ~ LogNormal(log(6), 0.5) instead of log(15)\n")
cat("  Median length scale: 6 weeks (vs 15 weeks)\n\n")

model_short_gp <- cmdstan_model("05_model2_full_short_gp.stan")

cat("Fitting model...\n")
fit_short_gp <- model_short_gp$sample(
  data = stan_data,
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
fit_short_gp$save_object("../results/fit_model2_short_gp.rds")
cat("\nSaved: results/fit_model2_short_gp.rds\n")

# ==============================================================================
# 3. DIAGNOSTICS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DIAGNOSTICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

diag <- fit_short_gp$diagnostic_summary()
cat(sprintf("Divergent transitions: %d\n", sum(diag$num_divergent)))
cat(sprintf("Max treedepth exceeded: %d\n", sum(diag$num_max_treedepth)))

# ==============================================================================
# 4. PARAMETER SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PARAMETER SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Key parameters
params <- c("mu", "alpha", "rho", "phi",
            "wolbachia_effect", "npi_effect", "temp_effect", "rain_effect")
summ <- fit_short_gp$summary(variables = params)
print(summ)

# ==============================================================================
# 5. POSTERIOR PREDICTIVE CHECK - RESIDUAL AUTOCORRELATION
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("RESIDUAL AUTOCORRELATION CHECK\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Extract posterior predictive samples
y_rep <- fit_short_gp$draws("cases_pred", format = "matrix")
S <- stan_data$S
N <- stan_data$N
y <- stan_data$cases[(S + 1):N]

# Align dimensions
n_common <- min(ncol(y_rep), length(y))
y_rep <- y_rep[, 1:n_common]
y <- y[1:n_common]

# Compute residuals
y_rep_median <- apply(y_rep, 2, median)
y_rep_sd <- apply(y_rep, 2, sd)
std_resid <- (y - y_rep_median) / y_rep_sd

# Autocorrelation
acf_resid <- acf(std_resid, lag.max = 5, plot = FALSE)
cat("Residual autocorrelation:\n")
for (i in 1:5) {
  cat(sprintf("  Lag %d: %.3f\n", i, acf_resid$acf[i + 1]))
}

# Coverage
q05 <- apply(y_rep, 2, quantile, 0.05)
q95 <- apply(y_rep, 2, quantile, 0.95)
coverage_90 <- mean(y >= q05 & y <= q95)
cat(sprintf("\n90%% interval coverage: %.1f%% (nominal: 90%%)\n", 100 * coverage_90))

# ==============================================================================
# 6. COMPARE WITH STANDARD MODEL
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMPARISON WITH STANDARD MODEL (rho ~ 15 weeks)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

if (file.exists("../results/fit_model2.rds")) {
  fit_standard <- readRDS("../results/fit_model2.rds")

  # Standard model residuals
  y_rep_std <- fit_standard$draws("cases_pred", format = "matrix")
  y_rep_std <- y_rep_std[, 1:n_common]
  y_rep_median_std <- apply(y_rep_std, 2, median)
  y_rep_sd_std <- apply(y_rep_std, 2, sd)
  std_resid_std <- (y - y_rep_median_std) / y_rep_sd_std

  acf_std <- acf(std_resid_std, lag.max = 1, plot = FALSE)

  cat("Residual autocorrelation (lag 1):\n")
  cat(sprintf("  Standard model (rho ~ 15 weeks): %.3f\n", acf_std$acf[2]))
  cat(sprintf("  Short GP model (rho ~ 6 weeks):  %.3f\n", acf_resid$acf[2]))

  # Coverage comparison
  q05_std <- apply(y_rep_std, 2, quantile, 0.05)
  q95_std <- apply(y_rep_std, 2, quantile, 0.95)
  coverage_90_std <- mean(y >= q05_std & y <= q95_std)

  cat(sprintf("\n90%% interval coverage:\n"))
  cat(sprintf("  Standard model: %.1f%%\n", 100 * coverage_90_std))
  cat(sprintf("  Short GP model: %.1f%%\n", 100 * coverage_90))

  # Effect estimates comparison
  cat("\nIntervention effect estimates:\n")

  wolb_std <- fit_standard$summary("wolbachia_effect")
  wolb_new <- fit_short_gp$summary("wolbachia_effect")
  cat(sprintf("  Wolbachia effect:\n"))
  cat(sprintf("    Standard: %.3f (95%% CI: %.3f - %.3f)\n",
              wolb_std$median, wolb_std$q5, wolb_std$q95))
  cat(sprintf("    Short GP: %.3f (95%% CI: %.3f - %.3f)\n",
              wolb_new$median, wolb_new$q5, wolb_new$q95))

  npi_std <- fit_standard$summary("npi_effect")
  npi_new <- fit_short_gp$summary("npi_effect")
  cat(sprintf("  NPI effect:\n"))
  cat(sprintf("    Standard: %.3f (95%% CI: %.3f - %.3f)\n",
              npi_std$median, npi_std$q5, npi_std$q95))
  cat(sprintf("    Short GP: %.3f (95%% CI: %.3f - %.3f)\n",
              npi_new$median, npi_new$q5, npi_new$q95))

  # Length scale comparison
  rho_std <- fit_standard$summary("rho")
  rho_new <- fit_short_gp$summary("rho")
  cat(sprintf("\n  GP length scale (rho):\n"))
  cat(sprintf("    Standard: %.1f weeks (95%% CI: %.1f - %.1f)\n",
              rho_std$median, rho_std$q5, rho_std$q95))
  cat(sprintf("    Short GP: %.1f weeks (95%% CI: %.1f - %.1f)\n",
              rho_new$median, rho_new$q5, rho_new$q95))

} else {
  cat("Standard model not found. Run 06_fit_models.R first for comparison.\n")
}

# ==============================================================================
# 7. VISUALIZATION
# ==============================================================================

cat("\nGenerating comparison plots...\n")

# Residual time series
resid_df <- tibble(
  week = 1:length(std_resid),
  short_gp = std_resid
)

if (exists("std_resid_std")) {
  resid_df$standard <- std_resid_std

  resid_long <- resid_df |>
    pivot_longer(-week, names_to = "model", values_to = "residual") |>
    mutate(model = ifelse(model == "short_gp", "Short GP (rho~6 wks)",
                          "Standard (rho~15 wks)"))

  p_resid_compare <- ggplot(resid_long, aes(x = week, y = residual, color = model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(alpha = 0.7) +
    facet_wrap(~model, ncol = 1) +
    labs(
      title = "Standardized Residuals: Standard vs Short GP",
      x = "Week",
      y = "Standardized Residual"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  ggsave("../results/figures/residual_comparison_short_gp.png", p_resid_compare,
         width = 10, height = 8, dpi = 150)
  cat("  Saved: results/figures/residual_comparison_short_gp.png\n")
}

# ACF comparison
if (exists("acf_std")) {
  acf_df <- tibble(
    lag = 0:5,
    standard = acf_std$acf[1:6],
    short_gp = acf_resid$acf[1:6]
  ) |>
    pivot_longer(-lag, names_to = "model", values_to = "acf") |>
    mutate(model = ifelse(model == "short_gp", "Short GP (rho~6 wks)",
                          "Standard (rho~15 wks)"))

  p_acf <- ggplot(acf_df |> filter(lag > 0), aes(x = lag, y = acf, fill = model)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "blue", alpha = 0.5) +
    labs(
      title = "Residual Autocorrelation Comparison",
      subtitle = "Dashed lines show approximate 95% bounds for white noise",
      x = "Lag (weeks)",
      y = "Autocorrelation",
      fill = "Model"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave("../results/figures/acf_comparison_short_gp.png", p_acf,
         width = 8, height = 6, dpi = 150)
  cat("  Saved: results/figures/acf_comparison_short_gp.png\n")
}

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("The shorter GP length scale allows the model to capture faster epidemic\n")
cat("fluctuations, which should reduce residual autocorrelation.\n\n")

cat("Trade-off: A more flexible GP may absorb some signal that would otherwise\n")
cat("be attributed to covariates. Check if intervention effects changed.\n\n")

cat("Files saved:\n")
cat("  results/fit_model2_short_gp.rds\n")
cat("  results/figures/residual_comparison_short_gp.png\n")
cat("  results/figures/acf_comparison_short_gp.png\n")
