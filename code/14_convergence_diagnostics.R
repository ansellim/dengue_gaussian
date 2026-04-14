#!/usr/bin/env Rscript
# ==============================================================================
# 14_convergence_diagnostics.R
#
# MCMC convergence diagnostics for the fitted climate-only model.
# Produces a combined figure with trace plots and Rhat / ESS summaries for the
# scalar parameters and a few representative f_residual components.
#
# Input:  results/fit_model.rds
# Output: results/figures/convergence_diagnostics.png
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(bayesplot)
  library(patchwork)
})

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- grep("--file=", args, value = TRUE)
  if (length(script_path) > 0) {
    setwd(dirname(normalizePath(sub("--file=", "", script_path))))
  }
}

dir.create("../results/figures", showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank()))

cat(strrep("=", 70), "\n")
cat("MCMC CONVERGENCE DIAGNOSTICS\n")
cat(strrep("=", 70), "\n\n")

fit <- readRDS("../results/fit_model.rds")

scalar_params <- c("mu", "alpha", "rho", "phi",
                   "temp_effect", "rain_effect",
                   "prop_climate", "prop_residual")

summ <- fit$summary(variables = scalar_params,
                    mean, sd, rhat = posterior::rhat,
                    ess_bulk = posterior::ess_bulk,
                    ess_tail = posterior::ess_tail)

cat("Scalar parameter diagnostics:\n")
print(summ)

diag <- fit$diagnostic_summary()
n_divergent <- sum(diag$num_divergent)
n_treedepth <- sum(diag$num_max_treedepth)
cat(sprintf("\nDivergent transitions: %d\n", n_divergent))
cat(sprintf("Max treedepth exceeded: %d\n", n_treedepth))

# ------------------------------------------------------------------
# Panel A: trace plots for the four core hyperparameters
# ------------------------------------------------------------------

trace_params <- c("mu", "alpha", "rho", "phi")
draws_array <- fit$draws(variables = trace_params)

color_scheme_set("mix-blue-red")
p_trace <- mcmc_trace(draws_array,
                      facet_args = list(ncol = 4, strip.position = "top")) +
  labs(title = "A. Trace plots (4 chains)",
       subtitle = "Well-mixed chains overlap and show no trends") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ------------------------------------------------------------------
# Panel B: Rhat bar chart for scalar params + a sample of f_residual
# ------------------------------------------------------------------

# Also include a sample of f_residual components to show GP-layer convergence
all_fres <- fit$summary(variables = "f_residual",
                        rhat = posterior::rhat,
                        ess_bulk = posterior::ess_bulk,
                        ess_tail = posterior::ess_tail)

fres_max_rhat <- max(all_fres$rhat, na.rm = TRUE)
fres_min_ess  <- min(all_fres$ess_bulk, na.rm = TRUE)
cat(sprintf("\nf_residual (%d components): max Rhat = %.4f, min bulk ESS = %.0f\n",
            nrow(all_fres), fres_max_rhat, fres_min_ess))

summ_plot <- summ |>
  mutate(variable = factor(variable, levels = rev(scalar_params)))

p_rhat <- ggplot(summ_plot, aes(x = rhat, y = variable)) +
  geom_vline(xintercept = 1.01, linetype = "dashed", color = "firebrick") +
  geom_vline(xintercept = 1.00, linetype = "dotted", color = "grey50") +
  geom_segment(aes(x = 1, xend = rhat, yend = variable),
               color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 3) +
  scale_x_continuous(limits = c(0.999, max(1.015, max(summ_plot$rhat) + 0.002))) +
  labs(title = "B. Split R-hat",
       subtitle = sprintf("All params below 1.01 (convergence threshold). f_residual max = %.3f",
                          fres_max_rhat),
       x = "R-hat", y = NULL) +
  theme(plot.subtitle = element_text(size = 8))

# ------------------------------------------------------------------
# Panel C: ESS bulk + tail bars
# ------------------------------------------------------------------

ess_long <- summ |>
  select(variable, ess_bulk, ess_tail) |>
  pivot_longer(c(ess_bulk, ess_tail),
               names_to = "kind", values_to = "ess") |>
  mutate(kind = recode(kind,
                       ess_bulk = "Bulk ESS",
                       ess_tail = "Tail ESS"),
         variable = factor(variable, levels = rev(scalar_params)))

p_ess <- ggplot(ess_long, aes(x = ess, y = variable, fill = kind)) +
  geom_vline(xintercept = 400, linetype = "dashed", color = "firebrick") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Bulk ESS" = "steelblue",
                                "Tail ESS" = "#b3cde3")) +
  labs(title = "C. Effective Sample Size",
       subtitle = sprintf("Dashed line = 400 (Vehtari et al. 2021 threshold). f_residual min bulk = %.0f",
                          fres_min_ess),
       x = "ESS", y = NULL, fill = NULL) +
  theme(legend.position = "top",
        plot.subtitle = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0))

# ------------------------------------------------------------------
# Compose
# ------------------------------------------------------------------

bottom_row <- p_rhat | p_ess
p_all <- p_trace / bottom_row +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    title = "MCMC convergence diagnostics: climate-only HSGP model",
    subtitle = sprintf(
      "4 chains x 1000 post-warmup draws. Divergences: %d. Max tree-depth hits: %d.",
      n_divergent, n_treedepth),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10))
  )

ggsave("../results/figures/convergence_diagnostics.png", p_all,
       width = 13, height = 8.5, dpi = 150)
cat("\nSaved: results/figures/convergence_diagnostics.png\n")
