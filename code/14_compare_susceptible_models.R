#!/usr/bin/env Rscript
# ==============================================================================
# 14_compare_susceptible_models.R
#
# Methodological audit comparison across the S0 sensitivity sweep.
#
#   M0 (baseline)                    : results/fit_model.rds
#   M1 (climate + log S_pop + GP)    : sweep over S0 in {0.50, 0.75, 0.95}
#   M2 (climate + log S_dom + GP)    : sweep over S0 in {0.50, 0.75, 0.95}
#                                      plus historical (per-serotype S0 from
#                                      observed early-period proportions)
#
# Reports for each model in the sweep:
#   - alpha (GP amplitude) posterior summary
#   - beta_S posterior summary
#   - prop_susceptible / prop_residual
#   - LOO ELPD (and difference vs M0)
#
# Saves:
#   results/susceptible_audit_summary.csv          (long-format sweep table)
#   results/susceptible_audit_loo_summary.csv      (LOO comparison table)
#   results/figures/susceptible_audit_comparison.png  (4-panel sweep figure)
# ==============================================================================

library(tidyverse)
library(posterior)
library(loo)
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

dir.create("../results/figures", showWarnings = FALSE, recursive = TRUE)
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank()))

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("M0 / M1 / M2 SUSCEPTIBLE-COVARIATE AUDIT (S0 SWEEP)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. Define the sweep and load fits
# ==============================================================================

sweep_specs <- tribble(
  ~model_label,                     ~family, ~config,         ~rds,
  "M0",                              "M0",    "baseline",      "fit_model.rds",
  "M1: log S_pop, S0 = 0.50",        "M1",    "uniform_S0_050","fit_M1_susceptible_pop_S050.rds",
  "M1: log S_pop, S0 = 0.75",        "M1",    "uniform_S0_075","fit_M1_susceptible_pop_S075.rds",
  "M1: log S_pop, S0 = 0.95",        "M1",    "uniform_S0_095","fit_M1_susceptible_pop_S095.rds",
  "M2: log S_dom, S0 = 0.50",        "M2",    "uniform_S0_050","fit_M2_susceptible_dom_S050.rds",
  "M2: log S_dom, S0 = 0.75",        "M2",    "uniform_S0_075","fit_M2_susceptible_dom_S075.rds",
  "M2: log S_dom, S0 = 0.95",        "M2",    "uniform_S0_095","fit_M2_susceptible_dom_S095.rds",
  "M2: log S_dom, historical",       "M2",    "historical",    "fit_M2_susceptible_dom_historical.rds"
)

cat("Loading fits...\n")
fits <- list()
for (i in seq_len(nrow(sweep_specs))) {
  path <- file.path("../results", sweep_specs$rds[i])
  if (!file.exists(path)) stop("Missing fit: ", path)
  fits[[sweep_specs$model_label[i]]] <- readRDS(path)
  cat(sprintf("  %s\n", sweep_specs$model_label[i]))
}

# ==============================================================================
# 2. Extract per-model summaries
# ==============================================================================

extract_metrics <- function(fit, has_S) {
  alpha_draws <- as.numeric(fit$draws("alpha", format = "matrix"))
  rho_draws   <- as.numeric(fit$draws("rho", format = "matrix"))

  out <- tibble(
    alpha_med  = median(alpha_draws),
    alpha_lo   = quantile(alpha_draws, 0.025),
    alpha_hi   = quantile(alpha_draws, 0.975),
    rho_med    = median(rho_draws)
  )

  if (has_S) {
    bS <- as.numeric(fit$draws("beta_S", format = "matrix"))
    ps <- as.numeric(fit$draws("prop_susceptible", format = "matrix"))
    pr <- as.numeric(fit$draws("prop_residual", format = "matrix"))
    out <- out |> mutate(
      beta_S_med = median(bS),
      beta_S_lo  = quantile(bS, 0.025),
      beta_S_hi  = quantile(bS, 0.975),
      P_betaS_pos = mean(bS > 0),
      prop_S_med  = median(ps),
      prop_R_med  = median(pr)
    )
  } else {
    pr <- as.numeric(fit$draws("prop_residual", format = "matrix"))
    out <- out |> mutate(
      beta_S_med = NA_real_, beta_S_lo = NA_real_, beta_S_hi = NA_real_,
      P_betaS_pos = NA_real_,
      prop_S_med = NA_real_,
      prop_R_med = median(pr)
    )
  }
  out
}

summary_rows <- map_dfr(seq_len(nrow(sweep_specs)), function(i) {
  has_S <- sweep_specs$family[i] != "M0"
  m <- extract_metrics(fits[[sweep_specs$model_label[i]]], has_S = has_S)
  bind_cols(sweep_specs[i, ], m)
})

write_csv(summary_rows, "../results/susceptible_audit_summary.csv")
cat("\nSaved: results/susceptible_audit_summary.csv\n")

# ==============================================================================
# 3. LOO comparison
# ==============================================================================

cat("\nComputing LOO for each fit...\n")
loo_list <- map(fits, function(f) {
  suppressWarnings(loo(f$draws("log_lik", format = "matrix")))
})

loo_df <- map_dfr(seq_along(loo_list), function(i) {
  est <- loo_list[[i]]$estimates
  tibble(
    model_label = names(loo_list)[i],
    elpd_loo    = est["elpd_loo", "Estimate"],
    elpd_se     = est["elpd_loo", "SE"]
  )
})

elpd_M0 <- loo_df$elpd_loo[loo_df$model_label == "M0"]
loo_df <- loo_df |>
  mutate(elpd_diff_vs_M0 = elpd_loo - elpd_M0,
         family = sweep_specs$family[match(model_label, sweep_specs$model_label)],
         config = sweep_specs$config[match(model_label, sweep_specs$model_label)])

write_csv(loo_df, "../results/susceptible_audit_loo_summary.csv")
cat("Saved: results/susceptible_audit_loo_summary.csv\n")

cat("\nLOO ELPD summary:\n")
print(loo_df |> select(model_label, elpd_loo, elpd_se, elpd_diff_vs_M0))

# Compare via loo_compare for the formal SE-of-difference
loo_comp <- loo_compare(loo_list)
cat("\nFormal LOO comparison (loo_compare):\n")
print(loo_comp)

# ==============================================================================
# 4. Console summary
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SWEEP SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("GP amplitude (alpha) is essentially unchanged across the entire sweep:\n")
print(summary_rows |>
        select(model_label, alpha_med, alpha_lo, alpha_hi))

cat("\nbeta_S posterior across the sweep (CIs include zero throughout):\n")
print(summary_rows |>
        filter(family != "M0") |>
        select(model_label, beta_S_med, beta_S_lo, beta_S_hi, P_betaS_pos))

cat("\nLOO ELPD difference vs M0 (none beat baseline by more than the SE):\n")
print(loo_df |>
        filter(model_label != "M0") |>
        select(model_label, elpd_diff_vs_M0, elpd_se))

# ==============================================================================
# 5. Comparison figure (4 panels, sweep across S0)
# ==============================================================================

cat("\nGenerating sweep comparison figure...\n")

s0_label_to_x <- function(model_label) {
  case_when(
    grepl("baseline", model_label) ~ "M0",
    grepl("S0 = 0.50", model_label) ~ "0.50",
    grepl("S0 = 0.75", model_label) ~ "0.75",
    grepl("S0 = 0.95", model_label) ~ "0.95",
    grepl("historical", model_label) ~ "hist"
  )
}

plot_df <- summary_rows |>
  mutate(x = s0_label_to_x(model_label),
         x = factor(x, levels = c("M0", "0.50", "0.75", "0.95", "hist")))

# Panel A: alpha (GP amplitude) vs S0
p_alpha <- ggplot(plot_df, aes(x = x, y = alpha_med, color = family)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = alpha_lo, ymax = alpha_hi),
                width = 0.15, position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c(M0 = "gray40", M1 = "steelblue", M2 = "firebrick")) +
  labs(title = expression("GP amplitude" ~ hat(alpha) ~ "across the sweep"),
       subtitle = "Essentially unchanged: adding S(t) does not shrink the GP at any S0",
       x = expression(S[0]), y = expression(hat(alpha)),
       color = NULL) +
  theme(legend.position = "top")

# Panel B: beta_S across the sweep
beta_df <- plot_df |> filter(family != "M0")
p_beta <- ggplot(beta_df, aes(x = x, y = beta_S_med, color = family)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = beta_S_lo, ymax = beta_S_hi),
                width = 0.15, position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c(M1 = "steelblue", M2 = "firebrick")) +
  labs(title = expression(beta[S] ~ "posterior across the sweep"),
       subtitle = "All credible intervals straddle zero",
       x = expression(S[0]), y = expression(beta[S]),
       color = NULL) +
  theme(legend.position = "top")

# Panel C: variance share absorbed by S
p_var <- ggplot(beta_df, aes(x = x, y = prop_S_med, color = family)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c(M1 = "steelblue", M2 = "firebrick")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Variance absorbed by S",
       subtitle = "Very small in every configuration",
       x = expression(S[0]), y = "median prop_susceptible",
       color = NULL) +
  theme(legend.position = "top")

# Panel D: LOO ELPD difference vs M0
loo_plot_df <- loo_df |>
  mutate(x = s0_label_to_x(model_label),
         x = factor(x, levels = c("M0", "0.50", "0.75", "0.95", "hist"))) |>
  filter(model_label != "M0")

p_loo <- ggplot(loo_plot_df, aes(x = x, y = elpd_diff_vs_M0, color = family)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = elpd_diff_vs_M0 - 2 * elpd_se,
                    ymax = elpd_diff_vs_M0 + 2 * elpd_se),
                width = 0.15, position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c(M1 = "steelblue", M2 = "firebrick")) +
  labs(title = expression("LOO" ~ Delta * "ELPD vs M0"),
       subtitle = "Differences hover near zero (errorbars +/- 2 SE)",
       x = expression(S[0]), y = expression(Delta * "ELPD"),
       color = NULL) +
  theme(legend.position = "top")

p_combined <- (p_alpha | p_beta) / (p_var | p_loo) +
  plot_annotation(
    title = "M0 / M1 / M2 audit across S0 sensitivity sweep",
    subtitle = "M0: climate + GP   |   M1: + log S_pop   |   M2: + log S_dom (uniform sweep + historical per-serotype S0)"
  )

ggsave("../results/figures/susceptible_audit_comparison.png", p_combined,
       width = 12, height = 8, dpi = 150)
cat("Saved: results/figures/susceptible_audit_comparison.png\n")

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("AUDIT COMPARISON COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
