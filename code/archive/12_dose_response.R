#!/usr/bin/env Rscript
# ==============================================================================
# 12_dose_response.R
#
# Dose-Response Analysis: Wolbachia Coverage vs Rt Reduction
#
# This script computes and visualizes the relationship between Wolbachia
# coverage and the estimated effect on Rt.
#
# IMPORTANT CAVEAT: This analysis assumes a LINEAR dose-response relationship,
# which cannot be verified with the available data (coverage only observed
# increasing monotonically from 0% to ~30%). Extrapolation beyond observed
# coverage levels is speculative.
#
# ==============================================================================

library(tidyverse)
library(posterior)
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

dir.create("../results/figures", showWarnings = FALSE, recursive = TRUE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DOSE-RESPONSE ANALYSIS: WOLBACHIA COVERAGE vs Rt REDUCTION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD MODEL FIT AND DATA
# ==============================================================================

cat("Loading model fit and data...\n")

fit <- readRDS("../results/fit_model2_short_gp.rds")
model_data <- readRDS("../data/model_data.rds")

# Get observed Wolbachia coverage range
wolbachia_coverage <- model_data$stan_data$X_full[, 3]
max_observed_coverage <- max(wolbachia_coverage)
min_observed_coverage <- min(wolbachia_coverage)

cat(sprintf("  Observed coverage range: %.1f%% to %.1f%%\n",
            100 * min_observed_coverage, 100 * max_observed_coverage))

# ==============================================================================
# 2. EXTRACT POSTERIOR DRAWS FOR WOLBACHIA COEFFICIENT
# ==============================================================================

cat("\nExtracting posterior draws...\n")

# beta[3] is the Wolbachia coefficient
# The model: f_wolbachia = coverage * beta[3]
# So effect at coverage c is: exp(beta[3] * c)

# Extract beta[3] draws
beta_wolbachia_draws <- fit$draws("beta[3]", format = "matrix")[, 1]
n_draws <- length(beta_wolbachia_draws)

cat(sprintf("  Number of posterior draws: %d\n", n_draws))
cat(sprintf("  beta[3] posterior median: %.3f\n", median(beta_wolbachia_draws)))
cat(sprintf("  beta[3] posterior 95%% CI: [%.3f, %.3f]\n",
            quantile(beta_wolbachia_draws, 0.025),
            quantile(beta_wolbachia_draws, 0.975)))

# ==============================================================================
# 3. COMPUTE DOSE-RESPONSE CURVE
# ==============================================================================

cat("\nComputing dose-response curve...\n")

# Coverage levels to evaluate (0% to 100%)
# Include key policy-relevant levels explicitly
coverage_levels <- sort(unique(c(seq(0, 1, by = 0.05), 0.30, 0.70)))

# Compute effect (Rt multiplier) at each coverage level for each draw
# Effect at coverage c = exp(beta[3] * c)
dose_response_matrix <- sapply(coverage_levels, function(c) {
  exp(beta_wolbachia_draws * c)
})

# Compute summary statistics
dose_response_summary <- tibble(
  coverage = coverage_levels,
  coverage_pct = coverage_levels * 100,
  median = apply(dose_response_matrix, 2, median),
  mean = apply(dose_response_matrix, 2, mean),
  q025 = apply(dose_response_matrix, 2, quantile, 0.025),
  q975 = apply(dose_response_matrix, 2, quantile, 0.975),
  q10 = apply(dose_response_matrix, 2, quantile, 0.10),
  q90 = apply(dose_response_matrix, 2, quantile, 0.90),
  # Rt reduction (1 - multiplier)
  reduction_median = 1 - apply(dose_response_matrix, 2, median),
  reduction_q025 = 1 - apply(dose_response_matrix, 2, quantile, 0.975),  # Note: flipped

  reduction_q975 = 1 - apply(dose_response_matrix, 2, quantile, 0.025),
  # Is this extrapolation?
  extrapolated = coverage_levels > max_observed_coverage
)

# ==============================================================================
# 4. KEY RESULTS TABLE
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DOSE-RESPONSE ESTIMATES\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

key_coverages <- c(0.10, 0.20, 0.30, 0.50, 0.70, 1.00)
key_results <- dose_response_summary |>
  filter(coverage %in% key_coverages)

cat(sprintf("%-12s %-20s %-20s %-12s\n",
            "Coverage", "Rt Multiplier", "Rt Reduction", "Data Status"))
cat(sprintf("%s\n", paste(rep("-", 70), collapse = "")))

for (i in 1:nrow(key_results)) {
  row <- key_results[i, ]
  status <- ifelse(row$extrapolated, "EXTRAPOLATED", "Observed")
  cat(sprintf("%-12s %-20s %-20s %-12s\n",
              sprintf("%.0f%%", row$coverage_pct),
              sprintf("%.2f [%.2f, %.2f]", row$median, row$q025, row$q975),
              sprintf("%.0f%% [%.0f%%, %.0f%%]",
                      100 * row$reduction_median,
                      100 * row$reduction_q025,
                      100 * row$reduction_q975),
              status))
}

# Highlight 100% coverage result
cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PROJECTED EFFECT AT 100% COVERAGE (EXTRAPOLATION)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

result_100 <- dose_response_summary |> filter(coverage == 1)
cat(sprintf("Rt multiplier at 100%% coverage: %.2f [95%% CI: %.2f, %.2f]\n",
            result_100$median, result_100$q025, result_100$q975))
cat(sprintf("Rt reduction at 100%% coverage:  %.0f%% [95%% CI: %.0f%%, %.0f%%]\n",
            100 * result_100$reduction_median,
            100 * result_100$reduction_q025,
            100 * result_100$reduction_q975))

# Probability of protective effect at different coverage levels
cat("\nProbability of protective effect (Rt multiplier < 1):\n")
for (c in c(0.30, 0.50, 0.70, 1.00)) {
  idx <- which(coverage_levels == c)
  prob_protective <- mean(dose_response_matrix[, idx] < 1)
  cat(sprintf("  At %.0f%% coverage: %.1f%%\n", c * 100, prob_protective * 100))
}

# ==============================================================================
# 5. COMPARISON WITH AWED TRIAL
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMPARISON WITH AWED TRIAL (External Reference)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("AWED Trial (Utarini et al., 2021 NEJM):\n")
cat("  - Cluster-randomized trial in Yogyakarta, Indonesia\n")
cat("  - Wolbachia-treated areas had 77% protective efficacy\n")
cat("  - This corresponds to Rt multiplier of ~0.23\n\n")

cat("Our extrapolation to 100% coverage:\n")
cat(sprintf("  - Rt multiplier: %.2f [%.2f, %.2f]\n",
            result_100$median, result_100$q025, result_100$q975))
cat(sprintf("  - Rt reduction: %.0f%% [%.0f%%, %.0f%%]\n",
            100 * result_100$reduction_median,
            100 * result_100$reduction_q025,
            100 * result_100$reduction_q975))

awed_in_ci <- result_100$q025 <= 0.23 & result_100$q975 >= 0.23
cat(sprintf("\nAWED trial estimate (0.23) %s our 95%% CI\n",
            ifelse(awed_in_ci, "falls WITHIN", "falls OUTSIDE")))

# ==============================================================================
# 6. VISUALIZATION
# ==============================================================================

cat("\nGenerating dose-response visualization...\n")

theme_set(theme_minimal(base_size = 12))

# Main dose-response plot
p_dose_response <- ggplot(dose_response_summary, aes(x = coverage_pct)) +
  # 95% CI ribbon - different colors for observed vs extrapolated
  geom_ribbon(data = dose_response_summary |> filter(!extrapolated),
              aes(ymin = q025, ymax = q975),
              fill = "steelblue", alpha = 0.3) +
  geom_ribbon(data = dose_response_summary |> filter(extrapolated),
              aes(ymin = q025, ymax = q975),
              fill = "darkorange", alpha = 0.2) +
  # 80% CI ribbon
  geom_ribbon(data = dose_response_summary |> filter(!extrapolated),
              aes(ymin = q10, ymax = q90),
              fill = "steelblue", alpha = 0.3) +
  geom_ribbon(data = dose_response_summary |> filter(extrapolated),
              aes(ymin = q10, ymax = q90),
              fill = "darkorange", alpha = 0.2) +
  # Median line
  geom_line(data = dose_response_summary |> filter(!extrapolated),
            aes(y = median), color = "steelblue", linewidth = 1.2) +
  geom_line(data = dose_response_summary |> filter(extrapolated),
            aes(y = median), color = "darkorange", linewidth = 1.2, linetype = "dashed") +
  # Reference line at 1 (no effect)
  geom_hline(yintercept = 1, linetype = "dotted", color = "gray40") +
  # Vertical line at max observed coverage
  geom_vline(xintercept = max_observed_coverage * 100,
             linetype = "dashed", color = "gray40", linewidth = 0.8) +
  # AWED reference point
  annotate("point", x = 100, y = 0.23, shape = 18, size = 4, color = "darkred") +
  annotate("text", x = 95, y = 0.23, label = "AWED trial\n(77% efficacy)",
           hjust = 1, size = 3, color = "darkred") +
  # Annotation for boundary
  annotate("text", x = max_observed_coverage * 100 + 2, y = 0.95,
           label = "Observed data\nends here", hjust = 0, size = 3, color = "gray30") +
  # Labels
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.25), limits = c(0, 1.5)) +
  labs(
    title = "Dose-Response: Wolbachia Coverage vs Effect on Rt",
    subtitle = "Blue = observed range (0-30%); Orange dashed = extrapolation (assumes linearity)",
    x = "Wolbachia Coverage (%)",
    y = "Rt Multiplier\n(1 = no effect, <1 = protective)"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Rt reduction version (more intuitive for policy)
p_reduction <- ggplot(dose_response_summary, aes(x = coverage_pct)) +
  # 95% CI ribbon
  geom_ribbon(data = dose_response_summary |> filter(!extrapolated),
              aes(ymin = reduction_q025 * 100, ymax = reduction_q975 * 100),
              fill = "steelblue", alpha = 0.3) +
  geom_ribbon(data = dose_response_summary |> filter(extrapolated),
              aes(ymin = reduction_q025 * 100, ymax = reduction_q975 * 100),
              fill = "darkorange", alpha = 0.2) +
  # Median line
  geom_line(data = dose_response_summary |> filter(!extrapolated),
            aes(y = reduction_median * 100), color = "steelblue", linewidth = 1.2) +
  geom_line(data = dose_response_summary |> filter(extrapolated),
            aes(y = reduction_median * 100), color = "darkorange", linewidth = 1.2,
            linetype = "dashed") +
  # Reference line at 0 (no effect)
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
  # Vertical line at max observed coverage
  geom_vline(xintercept = max_observed_coverage * 100,
             linetype = "dashed", color = "gray40", linewidth = 0.8) +
  # AWED reference
  annotate("point", x = 100, y = 77, shape = 18, size = 4, color = "darkred") +
  annotate("text", x = 95, y = 77, label = "AWED trial", hjust = 1, size = 3, color = "darkred") +
  # Labels
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(breaks = seq(-50, 100, by = 25), limits = c(-50, 100)) +
  labs(
    title = "Projected Rt Reduction by Wolbachia Coverage",
    subtitle = "Negative values indicate Rt increase (within uncertainty bounds)",
    x = "Wolbachia Coverage (%)",
    y = "Rt Reduction (%)"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Combined figure
p_combined <- p_dose_response / p_reduction +
  plot_annotation(
    caption = str_wrap(
      "CAVEAT: Extrapolation beyond 30% coverage (orange region) assumes a linear dose-response relationship, which cannot be verified with available data. Coverage has only been observed increasing monotonically within a narrow range. Alternative dose-response shapes (threshold, saturation, or accelerating effects) are plausible but not distinguishable.",
      width = 100
    ),
    theme = theme(plot.caption = element_text(size = 9, color = "gray30", hjust = 0))
  )

ggsave("../results/figures/dose_response_wolbachia.png", p_combined,
       width = 10, height = 12, dpi = 150)
cat("  Saved: results/figures/dose_response_wolbachia.png\n")

# ==============================================================================
# 7. SUMMARY TABLE FOR REPORTING
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUMMARY FOR REPORTING\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Create a clean table for reporting
report_table <- dose_response_summary |>
  filter(coverage %in% c(0.10, 0.20, 0.30, 0.50, 0.70, 1.00)) |>
  mutate(
    coverage_label = sprintf("%.0f%%", coverage_pct),
    rt_multiplier = sprintf("%.2f (%.2f–%.2f)", median, q025, q975),
    rt_reduction = sprintf("%.0f%% (%.0f%%–%.0f%%)",
                           reduction_median * 100,
                           reduction_q025 * 100,
                           reduction_q975 * 100),
    status = ifelse(extrapolated, "Extrapolated*", "Observed")
  ) |>
  select(Coverage = coverage_label,
         `Rt Multiplier (95% CI)` = rt_multiplier,
         `Rt Reduction (95% CI)` = rt_reduction,
         Status = status)

print(report_table)

cat("\n* Extrapolation assumes linear dose-response; cannot be verified with data.\n")

# Save table
write_csv(report_table, "../results/dose_response_table.csv")
cat("\nSaved: results/dose_response_table.csv\n")

# ==============================================================================
# 8. CAVEATS
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("IMPORTANT CAVEATS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("1. LINEARITY ASSUMPTION: The dose-response curve assumes each additional\n")
cat("   percentage point of coverage yields the same proportional Rt reduction.\n")
cat("   This is a modeling assumption, not a verified biological relationship.\n\n")

cat("2. EXTRAPOLATION BEYOND DATA: Coverage has only been observed in the range\n")
cat(sprintf("   %.0f%% to %.0f%%. Projections beyond this range are speculative.\n\n",
            min_observed_coverage * 100, max_observed_coverage * 100))

cat("3. ALTERNATIVE DOSE-RESPONSE SHAPES that are equally plausible:\n")
cat("   - Threshold: No effect until critical coverage reached\n")
cat("   - Saturation: Diminishing returns at higher coverage\n")
cat("   - Accelerating: Herd effects amplify protection at high coverage\n\n")

cat("4. COVERAGE DEFINITION: 'Coverage' here means proportion of households in\n")
cat("   Wolbachia release areas, NOT the fraction of mosquitoes carrying\n")
cat("   Wolbachia or the fraction of transmission events blocked.\n\n")

cat("5. ECOLOGICAL STUDY: This is a national-level time series analysis, not\n")
cat("   a randomized trial. Confounding with time trends cannot be ruled out.\n")
