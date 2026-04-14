#!/usr/bin/env Rscript
# ==============================================================================
# 12_lag_sensitivity.R
#
# Sensitivity analysis for the climate covariate lag. Refits the climate-only
# model with temperature and rainfall lagged by 2, 3, 4 (baseline), and 6 weeks,
# and compares climate effect sizes and variance decomposition across lags.
#
# Input:  data/model_data.rds, code/05_model_climate_only.stan
# Output: results/fit_model_lag{2,3,4,6}.rds
#         results/lag_sensitivity_summary.csv
#         results/figures/lag_sensitivity.png
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

cat(strrep("=", 70), "\n")
cat("CLIMATE LAG SENSITIVITY ANALYSIS\n")
cat(strrep("=", 70), "\n\n")

# ==============================================================================
# 1. LOAD DATA AND MODEL
# ==============================================================================

cat("Loading data and compiling model...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data_base <- model_data$stan_data
stan_data_sens <- model_data$stan_data_sens

# Baseline GP hyperprior means (must be set because the Stan model now reads
# them as data; same values as 06_fit_model.R).
stan_data_base$log_alpha_mu <- -1.2
stan_data_base$log_rho_mu   <- log(6)

model <- cmdstan_model("05_model_climate_only.stan")

# Build per-lag X_climate matrices
lag_x <- list(
  `2` = stan_data_sens$X_climate_lag2,
  `3` = stan_data_sens$X_climate_lag3,
  `4` = stan_data_base$X_climate,
  `6` = stan_data_sens$X_climate_lag6
)

stopifnot(all(sapply(lag_x, function(x) !is.null(x))))
stopifnot(all(sapply(lag_x, function(x) !anyNA(x))))

cat(sprintf("  Lags to fit: %s\n", paste(names(lag_x), collapse = ", ")))
cat(sprintf("  N_model = %d weeks per fit\n", stan_data_base$N_model))

# ==============================================================================
# 2. FIT FOR EACH LAG
# ==============================================================================

fit_one_lag <- function(lag_weeks, X_climate) {
  cat(strrep("-", 70), "\n")
  cat(sprintf("Fitting climate-only model with climate lag = %s weeks\n", lag_weeks))
  cat(strrep("-", 70), "\n")

  stan_data_k <- stan_data_base
  stan_data_k$X_climate <- X_climate

  fit <- model$sample(
    data = stan_data_k,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95,
    max_treedepth = 12,
    seed = 42,
    refresh = 500
  )

  fit_path <- sprintf("../results/fit_model_lag%s.rds", lag_weeks)
  fit$save_object(fit_path)
  cat(sprintf("  Saved: %s\n", fit_path))

  diag <- fit$diagnostic_summary()
  cat(sprintf("  Divergences: %d | Max treedepth: %d\n",
              sum(diag$num_divergent), sum(diag$num_max_treedepth)))

  summ <- fit$summary(
    variables = c("mu", "temp_effect", "rain_effect",
                  "prop_climate", "prop_residual", "alpha", "rho")
  )

  extract_row <- function(var) {
    s <- summ[summ$variable == var, ]
    tibble(
      variable = var,
      median = s$median,
      q5 = s$q5,
      q95 = s$q95
    )
  }

  summary_long <- bind_rows(
    extract_row("temp_effect"),
    extract_row("rain_effect"),
    extract_row("prop_climate"),
    extract_row("prop_residual"),
    extract_row("alpha"),
    extract_row("rho"),
    extract_row("mu")
  ) |>
    mutate(lag_weeks = as.integer(lag_weeks),
           divergences = sum(diag$num_divergent))

  summary_long
}

all_results <- list()
for (lag_k in names(lag_x)) {
  all_results[[lag_k]] <- fit_one_lag(lag_k, lag_x[[lag_k]])
}

results_df <- bind_rows(all_results)

# ==============================================================================
# 3. SAVE SUMMARY TABLE
# ==============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("LAG SENSITIVITY SUMMARY\n")
cat(strrep("=", 70), "\n\n")

results_wide <- results_df |>
  mutate(
    value = sprintf("%.3f [%.3f, %.3f]", median, q5, q95)
  ) |>
  select(variable, lag_weeks, value) |>
  pivot_wider(names_from = lag_weeks, values_from = value,
              names_prefix = "lag_") |>
  arrange(variable)

print(results_wide, n = Inf)

write_csv(results_df, "../results/lag_sensitivity_summary.csv")
cat("\nSaved: results/lag_sensitivity_summary.csv\n")

# ==============================================================================
# 4. PLOT
# ==============================================================================

cat("\nGenerating lag sensitivity figure...\n")

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank()))

plot_panel <- function(var_name, y_label, null_line = NULL) {
  dat <- results_df |>
    filter(variable == var_name) |>
    mutate(lag_weeks = factor(lag_weeks, levels = c(2, 3, 4, 6)),
           is_baseline = lag_weeks == "4")

  p <- ggplot(dat, aes(x = lag_weeks, y = median)) +
    geom_errorbar(aes(ymin = q5, ymax = q95), width = 0.2, linewidth = 0.6) +
    geom_point(aes(color = is_baseline), size = 3) +
    scale_color_manual(
      values = c(`TRUE` = "#d7191c", `FALSE` = "#2b2b2b"),
      guide = "none"
    ) +
    labs(
      x = "Climate lag (weeks)",
      y = y_label,
      title = var_name
    )

  if (!is.null(null_line)) {
    p <- p + geom_hline(yintercept = null_line,
                         linetype = "dashed", color = "gray50")
  }

  p
}

p1 <- plot_panel("temp_effect",
                 expression("Temperature effect " * (exp(beta[temp]))),
                 null_line = 1) +
  labs(title = "Temperature effect on Rt")

p2 <- plot_panel("rain_effect",
                 expression("Rainfall effect " * (exp(beta[rain]))),
                 null_line = 1) +
  labs(title = "Rainfall effect on Rt")

p3 <- plot_panel("prop_climate",
                 "Climate variance share",
                 null_line = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "Variance share: climate")

p4 <- plot_panel("prop_residual",
                 "Residual variance share",
                 null_line = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "Variance share: residual GP")

fig <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Lag sensitivity for the climate-only + residual GP model",
    subtitle = "Baseline = 4 weeks (red). Median with 90% credible interval.",
    theme = theme(plot.title = element_text(face = "bold"))
  )

ggsave("../results/figures/lag_sensitivity.png", fig,
       width = 11, height = 8, dpi = 150)
cat("  Saved: results/figures/lag_sensitivity.png\n")

cat("\nDone.\n")
