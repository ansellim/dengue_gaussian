#!/usr/bin/env Rscript
# ==============================================================================
# 15_gp_prior_sensitivity.R
#
# Sweep the GP hyperprior MEANS for both alpha (amplitude) and rho (length
# scale). Tests whether the variance decomposition is robust across the
# prior choices for the GP, beyond just the existing tight-alpha sensitivity.
#
# Grid:
#   alpha median in {0.15, 0.30, 0.60}  -> log_alpha_mu in {log(0.15), ...}
#   rho   median in {3,    6,    15}    -> log_rho_mu   in {log(3), ...}
#
# Sigma is held fixed at 0.5 on the log scale (matches baseline).
#
# Output:
#   results/gp_prior_sensitivity_summary.csv
#   results/figures/gp_prior_sensitivity.png
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)

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

cat("======================================================================\n")
cat("GP PRIOR SENSITIVITY SWEEP\n")
cat("======================================================================\n\n")

model_data <- readRDS("../data/model_data.rds")
stan_data_base <- model_data$stan_data

cat("Compiling model...\n")
model <- cmdstan_model("05_model_climate_only.stan")

# ------------------------------------------------------------------
# Grid
# ------------------------------------------------------------------
alpha_medians <- c(0.15, 0.30, 0.60)
rho_medians   <- c(3, 6, 15)

grid <- expand_grid(
  alpha_median = alpha_medians,
  rho_median   = rho_medians
) |>
  mutate(
    log_alpha_mu = log(alpha_median),
    log_rho_mu   = log(rho_median),
    label        = sprintf("alpha~%.2f, rho~%dwks", alpha_median, rho_median)
  )

cat(sprintf("Sweeping %d (alpha, rho) combinations...\n\n", nrow(grid)))

# ------------------------------------------------------------------
# Fit each cell
# ------------------------------------------------------------------
results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  cat(sprintf("[%d/%d] %s\n", i, nrow(grid), grid$label[i]))

  stan_data_i <- c(
    stan_data_base,
    list(log_alpha_mu = grid$log_alpha_mu[i],
         log_rho_mu   = grid$log_rho_mu[i])
  )

  fit_i <- model$sample(
    data = stan_data_i,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95,
    max_treedepth = 12,
    seed = 42 + i,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )

  diag <- fit_i$diagnostic_summary()
  ndiv <- sum(diag$num_divergent)

  draws_i <- fit_i$draws(
    variables = c("alpha", "rho", "temp_effect", "rain_effect",
                  "prop_climate", "prop_residual", "phi"),
    format = "draws_matrix"
  )

  results[[i]] <- tibble(
    variable = colnames(draws_i),
    mean     = apply(draws_i, 2, mean),
    median   = apply(draws_i, 2, median),
    sd       = apply(draws_i, 2, sd),
    q2.5     = apply(draws_i, 2, quantile, 0.025),
    q97.5    = apply(draws_i, 2, quantile, 0.975)
  ) |>
    mutate(
      alpha_median_prior = grid$alpha_median[i],
      rho_median_prior   = grid$rho_median[i],
      label              = grid$label[i],
      divergences        = ndiv
    )

  summ <- results[[i]]  # for the print below

  cat(sprintf("    -> alpha=%.3f, rho=%.2f, prop_climate=%.3f, divergences=%d\n",
              summ |> filter(variable == "alpha") |> pull(median),
              summ |> filter(variable == "rho") |> pull(median),
              summ |> filter(variable == "prop_climate") |> pull(median),
              ndiv))
}

results_df <- bind_rows(results)
write_csv(results_df, "../results/gp_prior_sensitivity_summary.csv")
cat("\nSaved: results/gp_prior_sensitivity_summary.csv\n")

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------
plot_df <- results_df |>
  filter(variable %in% c("alpha", "rho", "prop_climate", "prop_residual", "temp_effect")) |>
  mutate(
    alpha_label = sprintf("alpha prior median: %.2f", alpha_median_prior),
    rho_label   = sprintf("rho prior median: %d wks", rho_median_prior),
    panel = case_when(
      variable == "alpha"         ~ "alpha (posterior)",
      variable == "rho"           ~ "rho (posterior, wks)",
      variable == "prop_climate"  ~ "prop_climate",
      variable == "prop_residual" ~ "prop_residual",
      variable == "temp_effect"   ~ "temp effect (xRt)"
    ),
    panel = factor(panel, levels = c(
      "prop_climate", "prop_residual", "alpha (posterior)", "rho (posterior, wks)", "temp effect (xRt)"
    ))
  )

p <- ggplot(plot_df, aes(x = factor(rho_median_prior),
                         y = median,
                         color = factor(alpha_median_prior),
                         group = factor(alpha_median_prior))) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5),
                position = position_dodge(width = 0.4), width = 0.2) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  facet_wrap(~ panel, scales = "free_y", nrow = 1) +
  scale_color_brewer(palette = "Set1", name = "alpha prior median") +
  labs(
    title = "GP prior sensitivity: posterior summaries across a 3x3 grid of prior means",
    subtitle = sprintf("alpha medians {0.15, 0.30, 0.60} x rho medians {3, 6, 15} weeks; baseline = (0.30, 6)"),
    x = "rho prior median (weeks)",
    y = "Posterior median (95% CrI)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("../results/figures/gp_prior_sensitivity.png", p,
       width = 13, height = 5, dpi = 150)
cat("Saved: results/figures/gp_prior_sensitivity.png\n")

# ------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------
cat("\n=== Variance decomposition across prior grid ===\n")
results_df |>
  filter(variable == "prop_climate") |>
  transmute(alpha_median_prior, rho_median_prior,
            prop_climate_pct = sprintf("%.1f%% [%.1f, %.1f]",
                                       100 * median,
                                       100 * q2.5,
                                       100 * q97.5),
            divergences) |>
  print()

cat("\n=== Range of prop_climate across the grid ===\n")
prop_clim <- results_df |> filter(variable == "prop_climate")
cat(sprintf("  min median: %.3f%%\n", 100 * min(prop_clim$median)))
cat(sprintf("  max median: %.3f%%\n", 100 * max(prop_clim$median)))
cat(sprintf("  range:      %.3f pp\n", 100 * (max(prop_clim$median) - min(prop_clim$median))))

cat("\nDone.\n")
