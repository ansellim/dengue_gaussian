#!/usr/bin/env Rscript
# ==============================================================================
# 17_export_for_vae.R
#
# Export GP posterior draws and case data for the Python VAE comparison script.
#
# Input:
#   results/fit_model3.rds, data/model_data.rds
#   results/serotype_switch_timing.csv, results/serotype_residual_monthly.csv
#
# Output:
#   results/vae_export/f_residual_draws.csv   (200 draws x N_model columns)
#   results/vae_export/Rt_draws.csv           (200 draws x N_model columns)
#   results/vae_export/f_residual_summary.csv (date, median, q025, q975)
#   results/vae_export/cases_weekly.csv       (date, cases, week_index)
#   results/vae_export/serotype_info.csv      (merged serotype switch + entropy)
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)

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

dir.create("../results/vae_export", showWarnings = FALSE, recursive = TRUE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("EXPORTING DATA FOR VAE COMPARISON\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA AND FIT
# ==============================================================================

cat("Loading data and model fit...\n")

model_data <- readRDS("../data/model_data.rds")
dates_model <- model_data$metadata$dates_model
N_model <- length(dates_model)
cat(sprintf("  N_model = %d weeks\n", N_model))

fit <- readRDS("../results/fit_model3.rds")
cat("  Loaded fit_model3.rds\n")

# ==============================================================================
# 2. EXTRACT AND SUBSAMPLE POSTERIOR DRAWS
# ==============================================================================

cat("\nExtracting posterior draws...\n")

set.seed(42)

# f_residual draws
f_draws_all <- fit$draws("f_residual", format = "matrix")
cat(sprintf("  f_residual: %d draws x %d weeks\n", nrow(f_draws_all), ncol(f_draws_all)))
idx <- sample(nrow(f_draws_all), size = min(200, nrow(f_draws_all)))
f_draws_sub <- f_draws_all[idx, ]

# Rt draws
Rt_draws_all <- fit$draws("Rt", format = "matrix")
cat(sprintf("  Rt: %d draws x %d weeks\n", nrow(Rt_draws_all), ncol(Rt_draws_all)))
Rt_draws_sub <- Rt_draws_all[idx, ]

# ==============================================================================
# 3. SAVE DRAWS (matrix CSVs, no row/column names from Stan)
# ==============================================================================

cat("\nSaving posterior draw matrices...\n")

# Strip Stan parameter names so columns are just V1, V2, ...
colnames(f_draws_sub) <- paste0("V", seq_len(ncol(f_draws_sub)))
colnames(Rt_draws_sub) <- paste0("V", seq_len(ncol(Rt_draws_sub)))

write_csv(as_tibble(f_draws_sub), "../results/vae_export/f_residual_draws.csv")
cat("  Saved: results/vae_export/f_residual_draws.csv\n")

write_csv(as_tibble(Rt_draws_sub), "../results/vae_export/Rt_draws.csv")
cat("  Saved: results/vae_export/Rt_draws.csv\n")

# ==============================================================================
# 4. f_residual SUMMARY
# ==============================================================================

cat("\nComputing f_residual summary...\n")

f_summary <- tibble(
  date   = dates_model,
  median = apply(f_draws_all, 2, median),
  q025   = apply(f_draws_all, 2, quantile, 0.025),
  q975   = apply(f_draws_all, 2, quantile, 0.975)
)

write_csv(f_summary, "../results/vae_export/f_residual_summary.csv")
cat("  Saved: results/vae_export/f_residual_summary.csv\n")

# ==============================================================================
# 5. CASES (full time series including burn-in)
# ==============================================================================

cat("\nExporting weekly cases...\n")

cases_raw <- read_csv("../data/raw_dengue_cases.csv", show_col_types = FALSE)

cases_out <- cases_raw |>
  mutate(week_index = row_number()) |>
  select(date, cases, week_index)

write_csv(cases_out, "../results/vae_export/cases_weekly.csv")
cat(sprintf("  Saved: results/vae_export/cases_weekly.csv (%d rows)\n", nrow(cases_out)))

# ==============================================================================
# 6. SEROTYPE INFO (merge switch timing + monthly residual/entropy)
# ==============================================================================

cat("\nExporting serotype info...\n")

serotype_monthly <- read_csv("../results/serotype_residual_monthly.csv",
                             show_col_types = FALSE)

# Keep columns useful for the VAE comparison
serotype_info <- serotype_monthly |>
  select(year_month, median, lower_95, upper_95,
         dominant, entropy, entropy_smooth,
         D1_prop, D2_prop, D3_prop, D4_prop)

write_csv(serotype_info, "../results/vae_export/serotype_info.csv")
cat(sprintf("  Saved: results/vae_export/serotype_info.csv (%d rows)\n",
            nrow(serotype_info)))

# Also copy switch timing for convenience
file.copy("../results/serotype_switch_timing.csv",
          "../results/vae_export/serotype_switch_timing.csv",
          overwrite = TRUE)
cat("  Copied: results/vae_export/serotype_switch_timing.csv\n")

# ==============================================================================
# DONE
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("EXPORT COMPLETE\n")
cat(sprintf("  Model dates: %s to %s (%d weeks)\n",
            min(dates_model), max(dates_model), N_model))
cat(sprintf("  Total case weeks: %d (including %d burn-in weeks)\n",
            nrow(cases_out), nrow(cases_out) - N_model))
cat("  Output directory: results/vae_export/\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
