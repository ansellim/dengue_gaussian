#!/usr/bin/env Rscript
# ==============================================================================
# 35_fit_nonlinear_climate.R
#
# Fit nonlinear climate-Rt model: replace linear beta*temp + beta*rain with
# separate 1D Gaussian Processes in climate space that learn the nonlinear
# mapping from temperature/rainfall to log(Rt).
#
# Hypothesis: the real climate-Rt relationship is nonlinear (threshold/hump-
# shaped), but a linear model on Singapore's narrow temperature range (25-31C)
# cannot detect it.
#
# Input:  data/model_data.rds, code/35_model_nonlinear_climate.stan
# Output: results/fit_nonlinear_climate.rds
#         results/figures/nonlinear_climate_response_temp.png
#         results/figures/nonlinear_climate_response_rain.png
#         results/figures/nonlinear_climate_variance_decomposition.png
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
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- grep("--file=", args, value = TRUE)
  if (length(script_path) > 0) {
    setwd(dirname(normalizePath(sub("--file=", "", script_path))))
  }
}

dir.create("../results", showWarnings = FALSE)
dir.create("../results/figures", showWarnings = FALSE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("FITTING NONLINEAR CLIMATE MODEL (GP IN CLIMATE SPACE)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data_base <- model_data$stan_data
meta <- model_data$metadata

cat(sprintf("  N_model = %d weeks\n", stan_data_base$N_model))
cat(sprintf("  Cases range: %d - %d\n",
            min(stan_data_base$cases), max(stan_data_base$cases)))

# ==============================================================================
# 2. PREPARE CLIMATE GP DATA
# ==============================================================================

cat("\nPreparing climate GP inputs...\n")

# Extract standardized temperature and rainfall from X_climate
# Column 1 = standardized temperature, Column 2 = standardized rainfall
temp_std <- stan_data_base$X_climate[, 1]
rain_std <- stan_data_base$X_climate[, 2]

cat(sprintf("  Temperature (standardized): range [%.2f, %.2f]\n",
            min(temp_std), max(temp_std)))
cat(sprintf("  Rainfall (standardized): range [%.2f, %.2f]\n",
            min(rain_std), max(rain_std)))

# Compute HSGP boundaries for climate GPs
# L = 1.5 * half-range of the input variable
L_temp <- 1.5 * (max(temp_std) - min(temp_std)) / 2
L_rain <- 1.5 * (max(rain_std) - min(rain_std)) / 2

cat(sprintf("  L_temp = %.3f\n", L_temp))
cat(sprintf("  L_rain = %.3f\n", L_rain))

# Number of basis functions for climate GPs
M_climate <- 20
cat(sprintf("  M_climate = %d basis functions\n", M_climate))

# Use Matern 3/2 for temporal GP (the default from Model 3)
kernel_type <- 2

# Build Stan data list
stan_data <- c(
  stan_data_base,
  list(
    kernel_type = kernel_type,
    M_climate   = M_climate,
    temp_std    = as.vector(temp_std),
    rain_std    = as.vector(rain_std),
    L_temp      = L_temp,
    L_rain      = L_rain
  )
)

# ==============================================================================
# 3. COMPILE AND FIT MODEL
# ==============================================================================

cat("\nCompiling nonlinear climate model...\n")
cat("  log(Rt) = mu + f_temp(temp) + f_rain(rain) + f_residual(t)\n")
cat("  f_temp ~ GP_SE in temperature space\n")
cat("  f_rain ~ GP_SE in rainfall space\n")
cat("  f_residual ~ GP_Matern32 in time\n\n")

model <- cmdstan_model("35_model_nonlinear_climate.stan")

cat("Fitting model (4 chains, 1000 warmup + 1000 sampling)...\n")
fit <- model$sample(
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
fit$save_object("../results/fit_nonlinear_climate.rds")
cat("\nSaved: results/fit_nonlinear_climate.rds\n")

# ==============================================================================
# 4. DIAGNOSTICS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DIAGNOSTICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

diag <- fit$diagnostic_summary()
cat(sprintf("Divergent transitions: %d\n", sum(diag$num_divergent)))
cat(sprintf("Max treedepth exceeded: %d\n", sum(diag$num_max_treedepth)))

# ==============================================================================
# 5. PARAMETER SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PARAMETER SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

params <- c("mu", "alpha_temp", "rho_temp", "alpha_rain", "rho_rain",
            "alpha", "rho", "phi",
            "prop_temp", "prop_rain", "prop_residual", "prop_climate")
summ <- fit$summary(variables = params)
print(summ)

# ==============================================================================
# 6. EXTRACT CLIMATE GP FUNCTIONS FOR PLOTTING
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("EXTRACTING CLIMATE RESPONSE CURVES\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Extract f_temp and f_rain posterior draws
f_temp_draws <- fit$draws("f_temp", format = "matrix")  # ndraws x N_model
f_rain_draws <- fit$draws("f_rain", format = "matrix")

# Recover original-scale temperature and rainfall
# temp_std = (temp_lag - temp_mean) / temp_sd => temp_lag = temp_std * temp_sd + temp_mean
temp_original <- temp_std * meta$temp_sd + meta$temp_mean
rain_original <- rain_std * meta$rain_sd + meta$rain_mean

cat(sprintf("  Temperature range (original): %.1f - %.1f C\n",
            min(temp_original), max(temp_original)))
cat(sprintf("  Rainfall range (original): %.0f - %.0f mm\n",
            min(rain_original), max(rain_original)))

# For plotting, we need to sort by the climate variable
# and compute posterior summaries

# --- Temperature response ---
temp_order <- order(temp_original)
temp_sorted <- temp_original[temp_order]
f_temp_sorted <- f_temp_draws[, temp_order]

# Compute posterior summaries
f_temp_median <- apply(f_temp_sorted, 2, median)
f_temp_q025 <- apply(f_temp_sorted, 2, quantile, 0.025)
f_temp_q975 <- apply(f_temp_sorted, 2, quantile, 0.975)
f_temp_q10 <- apply(f_temp_sorted, 2, quantile, 0.10)
f_temp_q90 <- apply(f_temp_sorted, 2, quantile, 0.90)

temp_plot_df <- tibble(
  temp = temp_sorted,
  median = f_temp_median,
  q025 = f_temp_q025,
  q975 = f_temp_q975,
  q10 = f_temp_q10,
  q90 = f_temp_q90
)

# --- Rainfall response ---
rain_order <- order(rain_original)
rain_sorted <- rain_original[rain_order]
f_rain_sorted <- f_rain_draws[, rain_order]

f_rain_median <- apply(f_rain_sorted, 2, median)
f_rain_q025 <- apply(f_rain_sorted, 2, quantile, 0.025)
f_rain_q975 <- apply(f_rain_sorted, 2, quantile, 0.975)
f_rain_q10 <- apply(f_rain_sorted, 2, quantile, 0.10)
f_rain_q90 <- apply(f_rain_sorted, 2, quantile, 0.90)

rain_plot_df <- tibble(
  rain = rain_sorted,
  median = f_rain_median,
  q025 = f_rain_q025,
  q975 = f_rain_q975,
  q10 = f_rain_q10,
  q90 = f_rain_q90
)

# ==============================================================================
# 7. FIGURE 1: TEMPERATURE RESPONSE CURVE
# ==============================================================================

cat("\nGenerating figures...\n")

# Try to load linear model for comparison
linear_fit <- NULL
linear_beta_temp <- NULL
linear_beta_rain <- NULL
tryCatch({
  linear_fit <- readRDS("../results/fit_model3.rds")
  linear_beta_temp <- median(linear_fit$draws("beta_climate[1]", format = "matrix"))
  linear_beta_rain <- median(linear_fit$draws("beta_climate[2]", format = "matrix"))
  cat("  Loaded linear model for comparison\n")
  cat(sprintf("    Linear beta_temp = %.4f (on standardized scale)\n", linear_beta_temp))
  cat(sprintf("    Linear beta_rain = %.4f (on standardized scale)\n", linear_beta_rain))
}, error = function(e) {
  cat("  Could not load linear model (results/fit_model3.rds) for comparison\n")
})

# Temperature response plot
p_temp <- ggplot(temp_plot_df, aes(x = temp)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#E63946", alpha = 0.15) +
  geom_ribbon(aes(ymin = q10, ymax = q90), fill = "#E63946", alpha = 0.25) +
  geom_line(aes(y = median), color = "#E63946", linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_rug(data = tibble(temp = temp_original), aes(x = temp),
           sides = "b", alpha = 0.1, color = "grey30") +
  labs(
    title = "Nonlinear Temperature-Rt Response (GP in Climate Space)",
    subtitle = "f_temp(temperature): learned effect on log(Rt)",
    x = expression("Temperature (" * degree * "C, lag-4 weeks)"),
    y = expression(f[temp] ~ "(effect on log " * R[t] * ")")
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# Add linear comparison if available
if (!is.null(linear_beta_temp)) {
  # Linear prediction: beta_temp * temp_std = beta_temp * (temp - mean) / sd
  temp_seq <- seq(min(temp_original), max(temp_original), length.out = 200)
  temp_std_seq <- (temp_seq - meta$temp_mean) / meta$temp_sd
  linear_pred <- linear_beta_temp * temp_std_seq

  p_temp <- p_temp +
    geom_line(data = tibble(temp = temp_seq, y = linear_pred),
              aes(x = temp, y = y),
              linetype = "dotted", color = "grey40", linewidth = 0.8) +
    annotate("text", x = min(temp_original) + 0.5, y = max(f_temp_q975) * 0.9,
             label = "Dotted = linear model", color = "grey40", hjust = 0, size = 3.5)
}

ggsave("../results/figures/nonlinear_climate_response_temp.png",
       p_temp, width = 10, height = 6, dpi = 300)
cat("  Saved: results/figures/nonlinear_climate_response_temp.png\n")

# ==============================================================================
# 8. FIGURE 2: RAINFALL RESPONSE CURVE
# ==============================================================================

p_rain <- ggplot(rain_plot_df, aes(x = rain)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#457B9D", alpha = 0.15) +
  geom_ribbon(aes(ymin = q10, ymax = q90), fill = "#457B9D", alpha = 0.25) +
  geom_line(aes(y = median), color = "#457B9D", linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_rug(data = tibble(rain = rain_original), aes(x = rain),
           sides = "b", alpha = 0.1, color = "grey30") +
  labs(
    title = "Nonlinear Rainfall-Rt Response (GP in Climate Space)",
    subtitle = "f_rain(rainfall): learned effect on log(Rt)",
    x = "Rainfall (mm/week, lag-4 weeks)",
    y = expression(f[rain] ~ "(effect on log " * R[t] * ")")
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# Add linear comparison if available
if (!is.null(linear_beta_rain)) {
  rain_seq <- seq(min(rain_original), max(rain_original), length.out = 200)
  rain_std_seq <- (rain_seq - meta$rain_mean) / meta$rain_sd
  linear_pred_rain <- linear_beta_rain * rain_std_seq

  p_rain <- p_rain +
    geom_line(data = tibble(rain = rain_seq, y = linear_pred_rain),
              aes(x = rain, y = y),
              linetype = "dotted", color = "grey40", linewidth = 0.8) +
    annotate("text", x = min(rain_original) + 10, y = max(f_rain_q975) * 0.9,
             label = "Dotted = linear model", color = "grey40", hjust = 0, size = 3.5)
}

ggsave("../results/figures/nonlinear_climate_response_rain.png",
       p_rain, width = 10, height = 6, dpi = 300)
cat("  Saved: results/figures/nonlinear_climate_response_rain.png\n")

# ==============================================================================
# 9. FIGURE 3: VARIANCE DECOMPOSITION COMPARISON
# ==============================================================================

# Nonlinear model variance proportions
prop_temp_draws <- fit$draws("prop_temp", format = "matrix")
prop_rain_draws <- fit$draws("prop_rain", format = "matrix")
prop_residual_draws <- fit$draws("prop_residual", format = "matrix")
prop_climate_draws <- fit$draws("prop_climate", format = "matrix")

nonlinear_decomp <- tibble(
  component = c("Temperature GP", "Rainfall GP", "Climate (total)", "Residual GP"),
  model = "Nonlinear (GP)",
  median = c(median(prop_temp_draws), median(prop_rain_draws),
             median(prop_climate_draws), median(prop_residual_draws)),
  q025 = c(quantile(prop_temp_draws, 0.025), quantile(prop_rain_draws, 0.025),
           quantile(prop_climate_draws, 0.025), quantile(prop_residual_draws, 0.025)),
  q975 = c(quantile(prop_temp_draws, 0.975), quantile(prop_rain_draws, 0.975),
           quantile(prop_climate_draws, 0.975), quantile(prop_residual_draws, 0.975))
)

# Try to get linear model decomposition for comparison
linear_decomp <- NULL
if (!is.null(linear_fit)) {
  tryCatch({
    prop_clim_lin <- linear_fit$draws("prop_climate", format = "matrix")
    prop_resid_lin <- linear_fit$draws("prop_residual", format = "matrix")

    linear_decomp <- tibble(
      component = c("Climate (total)", "Residual GP"),
      model = "Linear",
      median = c(median(prop_clim_lin), median(prop_resid_lin)),
      q025 = c(quantile(prop_clim_lin, 0.025), quantile(prop_resid_lin, 0.025)),
      q975 = c(quantile(prop_clim_lin, 0.975), quantile(prop_resid_lin, 0.975))
    )
  }, error = function(e) {
    cat("  Could not extract linear model variance decomposition\n")
  })
}

# Build comparison data frame
if (!is.null(linear_decomp)) {
  # For the comparison plot, focus on climate vs residual
  compare_df <- bind_rows(
    linear_decomp,
    nonlinear_decomp |> filter(component %in% c("Climate (total)", "Residual GP"))
  ) |>
    mutate(
      component = factor(component, levels = c("Climate (total)", "Residual GP")),
      model = factor(model, levels = c("Linear", "Nonlinear (GP)"))
    )
} else {
  compare_df <- nonlinear_decomp |>
    filter(component %in% c("Climate (total)", "Residual GP")) |>
    mutate(
      component = factor(component, levels = c("Climate (total)", "Residual GP")),
      model = "Nonlinear (GP)"
    )
}

p_var <- ggplot(compare_df, aes(x = component, y = median * 100, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = q025 * 100, ymax = q975 * 100),
                position = position_dodge(width = 0.7), width = 0.2) +
  scale_fill_manual(values = c("Linear" = "#A8DADC", "Nonlinear (GP)" = "#E63946")) +
  labs(
    title = "Variance Decomposition: Linear vs Nonlinear Climate Model",
    subtitle = "Proportion of log(Rt) variance explained by each component",
    x = NULL,
    y = "% of total variance",
    fill = "Model"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

# Also add the temperature/rainfall breakdown for the nonlinear model
nonlinear_detail <- nonlinear_decomp |>
  filter(component %in% c("Temperature GP", "Rainfall GP")) |>
  mutate(component = factor(component))

p_detail <- ggplot(nonlinear_detail, aes(x = component, y = median * 100)) +
  geom_col(fill = "#E63946", width = 0.5) +
  geom_errorbar(aes(ymin = q025 * 100, ymax = q975 * 100), width = 0.15) +
  labs(
    title = "Nonlinear Model: Climate Breakdown",
    subtitle = "Temperature vs rainfall contributions",
    x = NULL,
    y = "% of total variance"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

p_combined <- p_var + p_detail + plot_layout(widths = c(2, 1))

ggsave("../results/figures/nonlinear_climate_variance_decomposition.png",
       p_combined, width = 14, height = 6, dpi = 300)
cat("  Saved: results/figures/nonlinear_climate_variance_decomposition.png\n")

# ==============================================================================
# 10. SUMMARY COMPARISON
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("RESULTS SUMMARY: LINEAR vs NONLINEAR CLIMATE MODEL\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("--- Nonlinear model (GP in climate space) ---\n")
cat(sprintf("  Climate variance (total):  %.2f%% [%.2f%%, %.2f%%]\n",
            median(prop_climate_draws) * 100,
            quantile(prop_climate_draws, 0.025) * 100,
            quantile(prop_climate_draws, 0.975) * 100))
cat(sprintf("    Temperature GP:          %.2f%% [%.2f%%, %.2f%%]\n",
            median(prop_temp_draws) * 100,
            quantile(prop_temp_draws, 0.025) * 100,
            quantile(prop_temp_draws, 0.975) * 100))
cat(sprintf("    Rainfall GP:             %.2f%% [%.2f%%, %.2f%%]\n",
            median(prop_rain_draws) * 100,
            quantile(prop_rain_draws, 0.025) * 100,
            quantile(prop_rain_draws, 0.975) * 100))
cat(sprintf("  Residual variance:         %.2f%% [%.2f%%, %.2f%%]\n",
            median(prop_residual_draws) * 100,
            quantile(prop_residual_draws, 0.025) * 100,
            quantile(prop_residual_draws, 0.975) * 100))

if (!is.null(linear_decomp)) {
  cat("\n--- Linear model (beta * x) ---\n")
  cat(sprintf("  Climate variance:          %.2f%% [%.2f%%, %.2f%%]\n",
              linear_decomp$median[1] * 100,
              linear_decomp$q025[1] * 100,
              linear_decomp$q975[1] * 100))
  cat(sprintf("  Residual variance:         %.2f%% [%.2f%%, %.2f%%]\n",
              linear_decomp$median[2] * 100,
              linear_decomp$q025[2] * 100,
              linear_decomp$q975[2] * 100))

  ratio <- median(prop_climate_draws) / linear_decomp$median[1]
  cat(sprintf("\n  Variance ratio (nonlinear / linear climate): %.1fx\n", ratio))
}

cat("\n--- GP hyperparameters ---\n")
cat(sprintf("  Temperature GP: alpha = %.4f, rho = %.3f (in std units)\n",
            median(fit$draws("alpha_temp", format = "matrix")),
            median(fit$draws("rho_temp", format = "matrix"))))
cat(sprintf("  Rainfall GP:    alpha = %.4f, rho = %.3f (in std units)\n",
            median(fit$draws("alpha_rain", format = "matrix")),
            median(fit$draws("rho_rain", format = "matrix"))))
cat(sprintf("  Temporal GP:    alpha = %.4f, rho = %.3f weeks\n",
            median(fit$draws("alpha", format = "matrix")),
            median(fit$draws("rho", format = "matrix"))))

cat("\n--- Shape assessment ---\n")
# Check for nonlinearity by comparing GP range to posterior uncertainty
f_temp_range <- max(f_temp_median) - min(f_temp_median)
f_rain_range <- max(f_rain_median) - min(f_rain_median)
avg_temp_width <- mean(f_temp_q975 - f_temp_q025)
avg_rain_width <- mean(f_rain_q975 - f_rain_q025)

cat(sprintf("  Temperature GP: range = %.4f, avg 95%% CrI width = %.4f\n",
            f_temp_range, avg_temp_width))
if (f_temp_range > avg_temp_width) {
  cat("    => Signal detected: GP range exceeds posterior uncertainty\n")
} else {
  cat("    => Weak/no signal: GP range within posterior uncertainty\n")
}

cat(sprintf("  Rainfall GP:    range = %.4f, avg 95%% CrI width = %.4f\n",
            f_rain_range, avg_rain_width))
if (f_rain_range > avg_rain_width) {
  cat("    => Signal detected: GP range exceeds posterior uncertainty\n")
} else {
  cat("    => Weak/no signal: GP range within posterior uncertainty\n")
}

cat("\n")
cat("Files saved:\n")
cat("  results/fit_nonlinear_climate.rds\n")
cat("  results/figures/nonlinear_climate_response_temp.png\n")
cat("  results/figures/nonlinear_climate_response_rain.png\n")
cat("  results/figures/nonlinear_climate_variance_decomposition.png\n")
