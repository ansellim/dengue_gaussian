#!/usr/bin/env Rscript
# ==============================================================================
# 12_fit_model3_clusters.R
#
# Fit Model 3: full covariates + active dengue clusters (vector activity proxy)
# Restricted to Jul 2015 - Nov 2020 (cluster data availability window)
#
# Also fits Model 2 on the same subperiod for direct comparison.
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(lubridate)
library(loo)

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(commandArgs(trailingOnly = FALSE) |>
    grep("--file=", x = _, value = TRUE) |>
    sub("--file=", "", x = _) |>
    normalizePath() |>
    dirname())
}

dir.create("../results", showWarnings = FALSE)
dir.create("../results/figures", showWarnings = FALSE)

cat("======================================================================\n")
cat("MODEL 3: FULL + CLUSTERS (2015-2020)\n")
cat("======================================================================\n\n")

# ==============================================================================
# 1. PREPARE SUBPERIOD DATA
# ==============================================================================

cat("Loading and subsetting data...\n")

model_data <- readRDS("../data/model_data.rds")
df_full <- model_data$df
metadata <- model_data$metadata

# Load cluster data
clusters <- read_csv("../data/raw_cluster_counts.csv", show_col_types = FALSE) |>
  mutate(date = as.Date(date))

cat(sprintf("  Cluster data: %d weeks (%s to %s)\n",
            nrow(clusters), min(clusters$date), max(clusters$date)))

# Merge clusters into main dataframe by closest date (within 7 days)
df_full <- df_full |> mutate(date = as.Date(date))
df_merged <- df_full |>
  left_join(
    clusters |> select(date, n_clusters),
    by = join_by(closest(date >= date))
  )

# Subset to period where cluster data exists (Jul 2015 - Nov 2020)
cluster_start <- min(clusters$date)
cluster_end <- max(clusters$date)

df_sub <- df_merged |>
  filter(date >= cluster_start, date <= cluster_end) |>
  filter(!is.na(n_clusters), !is.na(temp_std), !is.na(rain_std))

cat(sprintf("  Subperiod: %s to %s (%d weeks)\n",
            min(df_sub$date), max(df_sub$date), nrow(df_sub)))
cat(sprintf("  Clusters: min=%d, median=%d, max=%d\n",
            min(df_sub$n_clusters), median(df_sub$n_clusters), max(df_sub$n_clusters)))

# Standardize cluster count
cluster_mean <- mean(df_sub$n_clusters)
cluster_sd <- sd(df_sub$n_clusters)
df_sub <- df_sub |> mutate(clusters_std = (n_clusters - cluster_mean) / cluster_sd)

cat(sprintf("  Clusters standardized: mean=%.1f, sd=%.1f\n", cluster_mean, cluster_sd))

# ==============================================================================
# 2. FORMAT FOR STAN
# ==============================================================================

N_sub <- nrow(df_sub)
S <- model_data$stan_data$S  # Same generation interval
gi <- model_data$stan_data$gi

N_model_sub <- N_sub - S
M <- model_data$stan_data$M  # Same basis functions

# Time (centered, in weeks)
t_model <- 1:N_model_sub
t_center <- mean(t_model)
t_centered <- t_model - t_center
L <- 1.5 * (max(t_model) - min(t_model)) / 2

# Covariates (subperiod)
X_with_clusters <- cbind(
  temp = df_sub$temp_std[(S + 1):N_sub],
  rain = df_sub$rain_std[(S + 1):N_sub],
  wolbachia = df_sub$wolbachia_coverage[(S + 1):N_sub],
  npi = df_sub$npi_intensity[(S + 1):N_sub],
  clusters = df_sub$clusters_std[(S + 1):N_sub]
)

X_no_clusters <- X_with_clusters[, 1:4]

stan_data_sub <- list(
  N = N_sub,
  N_model = N_model_sub,
  S = S,
  M = M,
  cases = df_sub$cases,
  t = t_centered,
  L = L,
  gi = gi
)

# Model 3 data (with clusters)
stan_data_m3 <- c(stan_data_sub, list(
  K = 5,
  X = X_with_clusters
))

# Model 2 data on same subperiod (without clusters, for comparison)
stan_data_m2sub <- c(stan_data_sub, list(
  K_full = 4,
  X_full = X_no_clusters
))

cat(sprintf("\n  Stan data: N=%d, N_model=%d, S=%d, M=%d, L=%.1f\n",
            N_sub, N_model_sub, S, M, L))

dates_sub <- df_sub$date[(S + 1):N_sub]

# ==============================================================================
# 3. COMPILE AND FIT
# ==============================================================================

MCMC <- list(
  chains = 5, parallel_chains = 5,
  iter_warmup = 1000, iter_sampling = 1000,
  adapt_delta = 0.95, max_treedepth = 12, seed = 42
)

# --- Model 2 on subperiod (baseline comparison) ---
cat("\nCompiling Model 2 (subperiod)...\n")
model2 <- cmdstan_model("05_model2_full.stan")

cat("Fitting Model 2 on 2015-2020 subperiod...\n")
fit2_sub <- model2$sample(
  data = stan_data_m2sub,
  chains = MCMC$chains, parallel_chains = MCMC$parallel_chains,
  iter_warmup = MCMC$iter_warmup, iter_sampling = MCMC$iter_sampling,
  adapt_delta = MCMC$adapt_delta, max_treedepth = MCMC$max_treedepth,
  seed = MCMC$seed, refresh = 200
)
fit2_sub$save_object("../results/fit_model2_subperiod.rds")
cat("  Saved: results/fit_model2_subperiod.rds\n")

# --- Model 3 (with clusters) ---
cat("\nCompiling Model 3 (with clusters)...\n")
model3 <- cmdstan_model("05_model3_clusters.stan")

cat("Fitting Model 3 on 2015-2020 subperiod...\n")
fit3 <- model3$sample(
  data = stan_data_m3,
  chains = MCMC$chains, parallel_chains = MCMC$parallel_chains,
  iter_warmup = MCMC$iter_warmup, iter_sampling = MCMC$iter_sampling,
  adapt_delta = MCMC$adapt_delta, max_treedepth = MCMC$max_treedepth,
  seed = MCMC$seed, refresh = 200
)
fit3$save_object("../results/fit_model3_clusters.rds")
cat("  Saved: results/fit_model3_clusters.rds\n")

# ==============================================================================
# 4. DIAGNOSTICS AND COMPARISON
# ==============================================================================

cat("\n======================================================================\n")
cat("DIAGNOSTICS\n")
cat("======================================================================\n\n")

cat("--- Model 2 (subperiod, no clusters) ---\n")
cat("Divergences:", sum(fit2_sub$diagnostic_summary()$num_divergent), "\n")
cat("Max treedepth:", sum(fit2_sub$diagnostic_summary()$num_max_treedepth), "\n")
s2 <- fit2_sub$summary(variables = c("mu", "beta", "alpha", "rho", "phi",
                                       "temp_effect", "rain_effect",
                                       "wolbachia_effect", "npi_effect",
                                       "prop_climate", "prop_wolbachia",
                                       "prop_npi", "prop_residual"))
print(s2, n = 20)

cat("\n--- Model 3 (with clusters) ---\n")
cat("Divergences:", sum(fit3$diagnostic_summary()$num_divergent), "\n")
cat("Max treedepth:", sum(fit3$diagnostic_summary()$num_max_treedepth), "\n")
s3 <- fit3$summary(variables = c("mu", "beta", "alpha", "rho", "phi",
                                   "temp_effect", "rain_effect",
                                   "wolbachia_effect", "npi_effect",
                                   "cluster_effect",
                                   "prop_climate", "prop_wolbachia",
                                   "prop_npi", "prop_clusters",
                                   "prop_residual"))
print(s3, n = 20)

# ==============================================================================
# 5. LOO COMPARISON
# ==============================================================================

cat("\n--- LOO-CV Comparison ---\n")
loo2 <- fit2_sub$loo()
loo3 <- fit3$loo()
comp <- loo_compare(list(model2_sub = loo2, model3_clusters = loo3))
print(comp)

# ==============================================================================
# 6. VARIANCE DECOMPOSITION COMPARISON PLOT
# ==============================================================================

cat("\nGenerating comparison plots...\n")
theme_set(theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank()))

# Extract variance proportions
var_m2 <- fit2_sub$summary(variables = c("prop_climate", "prop_wolbachia",
                                          "prop_npi", "prop_residual"))
var_m3 <- fit3$summary(variables = c("prop_climate", "prop_wolbachia",
                                       "prop_npi", "prop_clusters",
                                       "prop_residual"))

var_df <- bind_rows(
  var_m2 |> mutate(model = "Model 2 (no clusters)"),
  var_m3 |> mutate(model = "Model 3 (with clusters)")
) |>
  mutate(
    component = str_replace(variable, "prop_", "") |>
      str_to_title() |>
      str_replace("Npi", "NPI"),
    component = factor(component, levels = c("Climate", "Wolbachia", "NPI", "Clusters", "Residual"))
  )

p_var <- var_df |>
  ggplot(aes(x = component, y = median, fill = model)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = q5, ymax = q95),
                position = position_dodge(0.8), width = 0.3) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Variance Decomposition: Model 2 vs Model 3",
       subtitle = "2015-2020 subperiod; error bars: 90% CI",
       y = "Proportion of log(Rt) variance", x = NULL, fill = NULL) +
  theme(legend.position = "top")

ggsave("../results/figures/variance_m2_vs_m3.png", p_var, width = 10, height = 6, dpi = 150)
cat("  Saved: results/figures/variance_m2_vs_m3.png\n")

# PPC comparison
cases_obs <- df_sub$cases[(S + 1):N_sub]
pred2 <- fit2_sub$draws("cases_pred", format = "matrix")
pred3 <- fit3$draws("cases_pred", format = "matrix")

coverage_fn <- function(pred_mat, obs) {
  q025 <- apply(pred_mat, 2, quantile, 0.025)
  q975 <- apply(pred_mat, 2, quantile, 0.975)
  q10 <- apply(pred_mat, 2, quantile, 0.10)
  q90 <- apply(pred_mat, 2, quantile, 0.90)
  c(coverage_95 = mean(obs >= q025 & obs <= q975),
    coverage_80 = mean(obs >= q10 & obs <= q90))
}

cov2 <- coverage_fn(pred2, cases_obs)
cov3 <- coverage_fn(pred3, cases_obs)

cat(sprintf("\nPosterior predictive coverage (2015-2020):
  Model 2: 80%%=%.1f%%, 95%%=%.1f%%
  Model 3: 80%%=%.1f%%, 95%%=%.1f%%\n",
  cov2["coverage_80"]*100, cov2["coverage_95"]*100,
  cov3["coverage_80"]*100, cov3["coverage_95"]*100))

cat("\n======================================================================\n")
cat("MODEL 3 FITTING COMPLETE\n")
cat("======================================================================\n")
