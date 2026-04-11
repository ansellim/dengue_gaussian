#!/usr/bin/env Rscript
# ==============================================================================
# 16_kernel_comparison_postprocess.R
#
# Cross-kernel comparison: LOO-CV, WAIC, posterior predictive checks,
# variance decomposition, Rt trajectories, GP hyperparameters, and
# serotype signal detectability across 4 GP kernels.
#
# Input:  results/fit_kernel_{matern12,matern32,matern52,sqexp}.rds
#         data/model_data.rds
#         results/serotype_switch_timing.csv (if available)
# Output: results/figures/kernel_*.png
#         results/kernel_comparison_summary.csv
#         results/kernel_loo_comparison.csv
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
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- grep("--file=", args, value = TRUE)
  if (length(script_path) > 0) {
    setwd(dirname(normalizePath(sub("--file=", "", script_path))))
  }
}

dir.create("../results/figures", showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank()))

# Kernel metadata
kernel_config <- tibble(
  name  = c("matern12", "matern32", "matern52", "sqexp"),
  label = c("Matérn 1/2", "Matérn 3/2", "Matérn 5/2", "Squared Exp."),
  color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
  smoothness = c("C⁰ (continuous)", "C¹ (once diff.)",
                 "C² (twice diff.)", "C∞ (infinitely smooth)")
)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("CROSS-KERNEL COMPARISON: POST-PROCESSING\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA AND FITS
# ==============================================================================

cat("Loading data and model fits...\n")

model_data <- readRDS("../data/model_data.rds")
stan_data  <- model_data$stan_data
dates      <- model_data$metadata$dates_model
S          <- stan_data$S
N          <- stan_data$N
y_obs      <- stan_data$cases[(S + 1):N]
N_model    <- length(y_obs)

fits <- list()
for (kname in kernel_config$name) {
  f <- sprintf("../results/fit_kernel_%s.rds", kname)
  if (file.exists(f)) {
    fits[[kname]] <- readRDS(f)
    cat(sprintf("  Loaded: %s\n", f))
  } else {
    cat(sprintf("  MISSING: %s — skipping this kernel\n", f))
  }
}

if (length(fits) < 2) {
  stop("Need at least 2 kernel fits for comparison. Run 15_fit_kernel_comparison.R first.")
}

available_kernels <- kernel_config |> filter(name %in% names(fits))
cat(sprintf("\n  %d kernels loaded for comparison.\n\n", nrow(available_kernels)))

# ==============================================================================
# 2. MCMC DIAGNOSTICS TABLE
# ==============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MCMC DIAGNOSTICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

diag_table <- map_dfr(available_kernels$name, function(kname) {
  fit <- fits[[kname]]
  diag <- fit$diagnostic_summary()
  summ_all <- fit$summary()
  tibble(
    kernel       = kname,
    divergences  = sum(diag$num_divergent),
    max_treedepth_exceeded = sum(diag$num_max_treedepth),
    max_rhat     = max(summ_all$rhat, na.rm = TRUE),
    min_ess_bulk = min(summ_all$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(summ_all$ess_tail, na.rm = TRUE)
  )
})

print(diag_table, n = Inf, width = Inf)
cat("\n")

# ==============================================================================
# 3. LOO-CV COMPARISON
# ==============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("LOO-CV COMPARISON\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

loo_results <- list()
for (kname in available_kernels$name) {
  log_lik <- fits[[kname]]$draws("log_lik", format = "matrix")
  loo_results[[kname]] <- loo(log_lik, cores = 4)
  cat(sprintf("  %s: ELPD = %.1f (SE = %.1f), p_loo = %.1f\n",
              kname,
              loo_results[[kname]]$estimates["elpd_loo", "Estimate"],
              loo_results[[kname]]$estimates["elpd_loo", "SE"],
              loo_results[[kname]]$estimates["p_loo", "Estimate"]))
}

cat("\nLOO-CV Model Comparison:\n")
loo_comp <- loo_compare(loo_results)
print(loo_comp)

# Save LOO comparison
loo_comp_df <- as.data.frame(loo_comp)
loo_comp_df$kernel <- rownames(loo_comp_df)
write_csv(loo_comp_df, "../results/kernel_loo_comparison.csv")
cat("\n  Saved: results/kernel_loo_comparison.csv\n\n")

# ==============================================================================
# 4. WAIC COMPARISON
# ==============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("WAIC COMPARISON\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

waic_results <- list()
for (kname in available_kernels$name) {
  log_lik <- fits[[kname]]$draws("log_lik", format = "matrix")
  waic_results[[kname]] <- waic(log_lik)
  cat(sprintf("  %s: WAIC ELPD = %.1f (SE = %.1f), p_waic = %.1f\n",
              kname,
              waic_results[[kname]]$estimates["elpd_waic", "Estimate"],
              waic_results[[kname]]$estimates["elpd_waic", "SE"],
              waic_results[[kname]]$estimates["p_waic", "Estimate"]))
}

cat("\nWAIC Model Comparison:\n")
waic_comp <- loo_compare(waic_results)
print(waic_comp)
cat("\n")

# ==============================================================================
# 5. GP HYPERPARAMETER POSTERIORS
# ==============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("GP HYPERPARAMETER POSTERIORS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

gp_posteriors <- map_dfr(available_kernels$name, function(kname) {
  fit <- fits[[kname]]
  summ <- fit$summary(variables = c("alpha", "rho"))
  summ$kernel <- kname
  summ$label  <- available_kernels$label[available_kernels$name == kname]
  summ
})
print(gp_posteriors |> select(kernel, variable, median, q5, q95), n = Inf)
cat("\n")

# Plot: GP hyperparameter comparison (ridge plots)
gp_draws <- map_dfr(available_kernels$name, function(kname) {
  fit <- fits[[kname]]
  alpha_d <- as.vector(fit$draws("alpha", format = "matrix"))
  rho_d   <- as.vector(fit$draws("rho", format = "matrix"))
  label <- available_kernels$label[available_kernels$name == kname]
  bind_rows(
    tibble(kernel = label, parameter = "alpha (amplitude)", value = alpha_d),
    tibble(kernel = label, parameter = "rho (length scale, weeks)", value = rho_d)
  )
}) |>
  mutate(kernel = factor(kernel, levels = rev(available_kernels$label)))

p_gp_hyper <- ggplot(gp_draws, aes(x = value, y = kernel, fill = kernel)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.5) +
  facet_wrap(~parameter, scales = "free_x") +
  scale_fill_manual(values = setNames(rev(available_kernels$color),
                                      rev(available_kernels$label))) +
  labs(title = "GP Hyperparameter Posteriors by Kernel",
       x = "Value", y = NULL) +
  theme(legend.position = "none")

ggsave("../results/figures/kernel_gp_hyperparameters.png", p_gp_hyper,
       width = 10, height = 5, dpi = 200)
cat("  Saved: results/figures/kernel_gp_hyperparameters.png\n\n")

# ==============================================================================
# 6. VARIANCE DECOMPOSITION COMPARISON
# ==============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("VARIANCE DECOMPOSITION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

vardecomp <- map_dfr(available_kernels$name, function(kname) {
  fit <- fits[[kname]]
  pc <- as.vector(fit$draws("prop_climate", format = "matrix"))
  pr <- as.vector(fit$draws("prop_residual", format = "matrix"))
  label <- available_kernels$label[available_kernels$name == kname]
  tibble(
    kernel = label,
    climate_median  = median(pc) * 100,
    climate_q5      = quantile(pc, 0.05) * 100,
    climate_q95     = quantile(pc, 0.95) * 100,
    residual_median = median(pr) * 100,
    residual_q5     = quantile(pr, 0.05) * 100,
    residual_q95    = quantile(pr, 0.95) * 100
  )
})

print(vardecomp, width = Inf)
cat("\n")

# Plot: Variance decomposition
vardecomp_draws <- map_dfr(available_kernels$name, function(kname) {
  fit <- fits[[kname]]
  pc <- as.vector(fit$draws("prop_climate", format = "matrix"))
  label <- available_kernels$label[available_kernels$name == kname]
  tibble(kernel = label, prop_climate = pc * 100)
}) |>
  mutate(kernel = factor(kernel, levels = available_kernels$label))

p_vardecomp <- ggplot(vardecomp_draws, aes(x = prop_climate, fill = kernel)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = setNames(available_kernels$color,
                                      available_kernels$label)) +
  labs(title = "Climate Proportion of Rt Variance by Kernel",
       x = "% Rt Variance Explained by Climate",
       y = "Posterior Density", fill = "Kernel") +
  coord_cartesian(xlim = c(0, 5))

ggsave("../results/figures/kernel_variance_decomposition.png", p_vardecomp,
       width = 8, height = 5, dpi = 200)
cat("  Saved: results/figures/kernel_variance_decomposition.png\n\n")

# ==============================================================================
# 7. POSTERIOR PREDICTIVE CHECKS
# ==============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("POSTERIOR PREDICTIVE CHECKS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

ppc_table <- map_dfr(available_kernels$name, function(kname) {
  fit <- fits[[kname]]
  y_rep <- fit$draws("cases_pred", format = "matrix")
  n <- min(ncol(y_rep), N_model)
  y_rep <- y_rep[, 1:n]
  y <- y_obs[1:n]

  q025 <- apply(y_rep, 2, quantile, 0.025)
  q10  <- apply(y_rep, 2, quantile, 0.10)
  q90  <- apply(y_rep, 2, quantile, 0.90)
  q975 <- apply(y_rep, 2, quantile, 0.975)

  label <- available_kernels$label[available_kernels$name == kname]
  tibble(
    kernel       = label,
    coverage_80  = mean(y >= q10 & y <= q90) * 100,
    coverage_95  = mean(y >= q025 & y <= q975) * 100
  )
})

print(ppc_table, width = Inf)
cat("\n")

# ==============================================================================
# 8. Rt TRAJECTORY COMPARISON
# ==============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("Rt TRAJECTORY COMPARISON\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

rt_comparison <- map_dfr(available_kernels$name, function(kname) {
  fit <- fits[[kname]]
  Rt_draws <- fit$draws("Rt", format = "matrix")
  n <- min(ncol(Rt_draws), length(dates))
  label <- available_kernels$label[available_kernels$name == kname]

  tibble(
    kernel = label,
    date   = dates[1:n],
    median = apply(Rt_draws[, 1:n], 2, median),
    q025   = apply(Rt_draws[, 1:n], 2, quantile, 0.025),
    q10    = apply(Rt_draws[, 1:n], 2, quantile, 0.10),
    q90    = apply(Rt_draws[, 1:n], 2, quantile, 0.90),
    q975   = apply(Rt_draws[, 1:n], 2, quantile, 0.975)
  )
}) |>
  mutate(kernel = factor(kernel, levels = available_kernels$label))

# Full Rt comparison plot
p_rt <- ggplot(rt_comparison, aes(x = date)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = kernel), alpha = 0.15) +
  geom_ribbon(aes(ymin = q10, ymax = q90, fill = kernel), alpha = 0.3) +
  geom_line(aes(y = median, color = kernel), linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = setNames(available_kernels$color,
                                       available_kernels$label)) +
  scale_fill_manual(values = setNames(available_kernels$color,
                                      available_kernels$label)) +
  labs(title = "Rt Trajectory Comparison Across GP Kernels",
       subtitle = "Shaded: 80% and 95% CrI",
       x = NULL, y = "Rt", color = "Kernel", fill = "Kernel") +
  theme(legend.position = "bottom")

ggsave("../results/figures/kernel_rt_comparison.png", p_rt,
       width = 12, height = 6, dpi = 200)
cat("  Saved: results/figures/kernel_rt_comparison.png\n")

# Faceted version for clarity
p_rt_facet <- ggplot(rt_comparison, aes(x = date)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "steelblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = q10, ymax = q90), fill = "steelblue", alpha = 0.4) +
  geom_line(aes(y = median), color = "steelblue", linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  facet_wrap(~kernel, ncol = 1) +
  labs(title = "Rt Trajectory by Kernel (Faceted)",
       subtitle = "Shaded: 80% and 95% CrI",
       x = NULL, y = "Rt")

ggsave("../results/figures/kernel_rt_faceted.png", p_rt_facet,
       width = 12, height = 10, dpi = 200)
cat("  Saved: results/figures/kernel_rt_faceted.png\n\n")

# ==============================================================================
# 9. RESIDUAL GP (f_residual) COMPARISON
# ==============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("RESIDUAL GP COMPARISON\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

fresid_comparison <- map_dfr(available_kernels$name, function(kname) {
  fit <- fits[[kname]]
  fr_draws <- fit$draws("f_residual", format = "matrix")
  n <- min(ncol(fr_draws), length(dates))
  label <- available_kernels$label[available_kernels$name == kname]

  tibble(
    kernel = label,
    date   = dates[1:n],
    median = apply(fr_draws[, 1:n], 2, median),
    q025   = apply(fr_draws[, 1:n], 2, quantile, 0.025),
    q975   = apply(fr_draws[, 1:n], 2, quantile, 0.975)
  )
}) |>
  mutate(kernel = factor(kernel, levels = available_kernels$label))

p_fresid <- ggplot(fresid_comparison, aes(x = date)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = median), color = "steelblue", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~kernel, ncol = 1) +
  labs(title = "Residual GP (f_residual) by Kernel",
       subtitle = "Shaded: 95% CrI. How smooth/rough each kernel fits the unexplained signal.",
       x = NULL, y = "f_residual (log Rt scale)")

ggsave("../results/figures/kernel_residual_gp.png", p_fresid,
       width = 12, height = 10, dpi = 200)
cat("  Saved: results/figures/kernel_residual_gp.png\n\n")

# ==============================================================================
# 10. SEROTYPE SIGNAL DETECTABILITY
# ==============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SEROTYPE SIGNAL DETECTABILITY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

switch_file <- "../results/serotype_switch_timing.csv"
if (file.exists(switch_file)) {

  switch_dates <- read_csv(switch_file, show_col_types = FALSE)
  cat(sprintf("  Loaded %d serotype switch events.\n\n", nrow(switch_dates)))

  # For each kernel and each switch, compute P(elevated f_residual) in 0-3 months post-switch
  signal_detect <- map_dfr(available_kernels$name, function(kname) {
    fit <- fits[[kname]]
    fr_draws <- fit$draws("f_residual", format = "matrix")
    n <- min(ncol(fr_draws), length(dates))
    label <- available_kernels$label[available_kernels$name == kname]

    map_dfr(seq_len(nrow(switch_dates)), function(s) {
      switch_date <- as.Date(switch_dates$switch_date[s])
      transition  <- switch_dates$switch_event[s]

      # Find weeks in 0-13 weeks (0-3 months) post-switch
      post_idx <- which(dates >= switch_date &
                        dates < switch_date + 91 &
                        seq_along(dates) <= n)

      if (length(post_idx) == 0) {
        return(tibble(kernel = label, switch_date = switch_date,
                      transition = switch_event,
                      p_elevated = NA_real_,
                      median_increase = NA_real_))
      }

      # For each draw, compute mean f_residual in post-switch window
      post_mean <- apply(fr_draws[, post_idx, drop = FALSE], 1, mean)

      tibble(
        kernel          = label,
        switch_date     = switch_date,
        transition      = transition,
        p_elevated      = mean(post_mean > 0),
        median_increase = median(post_mean)
      )
    })
  }) |>
    mutate(kernel = factor(kernel, levels = available_kernels$label))

  print(signal_detect, n = Inf, width = Inf)

  # Plot
  p_signal <- ggplot(signal_detect |> filter(!is.na(p_elevated)),
                     aes(x = transition, y = p_elevated, fill = kernel)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = setNames(available_kernels$color,
                                        available_kernels$label)) +
    geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey40") +
    labs(title = "Serotype Signal Detectability by Kernel",
         subtitle = "P(elevated f_residual) in 0-3 months post-switch. Dashed line = 0.75 threshold.",
         x = "Serotype Transition", y = "P(elevated f_residual)",
         fill = "Kernel") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom") +
    coord_cartesian(ylim = c(0, 1))

  ggsave("../results/figures/kernel_serotype_detectability.png", p_signal,
         width = 10, height = 6, dpi = 200)
  cat("\n  Saved: results/figures/kernel_serotype_detectability.png\n")

} else {
  cat("  Serotype switch file not found — skipping signal detectability analysis.\n")
  cat("  (Run 08_serotype_analysis.R first to generate it.)\n")
}

# ==============================================================================
# 11. SUMMARY TABLE
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMPREHENSIVE SUMMARY TABLE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

summary_table <- diag_table |>
  left_join(ppc_table |>
              mutate(name = available_kernels$name[match(kernel, available_kernels$label)]),
            by = c("kernel" = "name")) |>
  left_join(vardecomp |>
              mutate(name = available_kernels$name[match(kernel, available_kernels$label)]),
            by = c("kernel" = "name")) |>
  left_join(
    map_dfr(names(loo_results), function(kname) {
      tibble(
        kernel   = kname,
        elpd_loo = loo_results[[kname]]$estimates["elpd_loo", "Estimate"],
        se_loo   = loo_results[[kname]]$estimates["elpd_loo", "SE"],
        p_loo    = loo_results[[kname]]$estimates["p_loo", "Estimate"]
      )
    }),
    by = "kernel"
  )

write_csv(summary_table, "../results/kernel_comparison_summary.csv")
cat("  Saved: results/kernel_comparison_summary.csv\n\n")

# Print nicely
cat("Key comparison metrics:\n\n")
for (i in seq_len(nrow(summary_table))) {
  r <- summary_table[i, ]
  lbl <- available_kernels$label[available_kernels$name == r$kernel]
  cat(sprintf("  %s:\n", lbl))
  cat(sprintf("    ELPD-LOO: %.1f (SE %.1f) | p_loo: %.1f\n",
              r$elpd_loo, r$se_loo, r$p_loo))
  cat(sprintf("    PPC: 80%% cov = %.1f%%, 95%% cov = %.1f%%\n",
              r$coverage_80, r$coverage_95))
  cat(sprintf("    Variance: climate = %.1f%%, residual = %.1f%%\n",
              r$climate_median, r$residual_median))
  cat(sprintf("    Divergences: %d, Max Rhat: %.4f\n\n",
              r$divergences, r$max_rhat))
}

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("CROSS-KERNEL COMPARISON COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Figures saved:\n")
cat("  results/figures/kernel_gp_hyperparameters.png\n")
cat("  results/figures/kernel_variance_decomposition.png\n")
cat("  results/figures/kernel_rt_comparison.png\n")
cat("  results/figures/kernel_rt_faceted.png\n")
cat("  results/figures/kernel_residual_gp.png\n")
cat("  results/figures/kernel_serotype_detectability.png\n\n")
cat("Tables saved:\n")
cat("  results/kernel_comparison_summary.csv\n")
cat("  results/kernel_loo_comparison.csv\n")
