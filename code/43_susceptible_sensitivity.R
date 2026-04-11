#!/usr/bin/env Rscript
# ==============================================================================
# 43_susceptible_sensitivity.R
#
# Sensitivity analysis: does national-level susceptible depletion explain
# the GP residual's amplitude? Tests across a grid of assumptions:
#   - Expansion factor (EF): 6, 8, 10, 14, 20
#   - Initial susceptible fraction (S_init): 0.50, 0.75, 0.95
#   - Reset at serotype switches: Yes / No
#   - Cross-immunity duration: fixed at 0 (simplicity)
#   - Birth replenishment + waning immunity: both off / both on
#       Births: ~33,000/year (~635/week) -- new births susceptible to all serotypes
#       Waning: 0.003/year (Tan et al. 2019) -- immune individuals become susceptible
#
# The key diagnostic: does the range of log(S_eff) match the GP amplitude?
#   GP amplitude alpha ≈ 0.30 on log(Rt) scale
#   If range(log(S_eff)) << 0.30, susceptible depletion CAN'T explain Rt
#   If range(log(S_eff)) >= 0.30, it COULD explain Rt
#
# Input:  data/model_data.rds, serotype CSVs, results/fit_model3.rds
# Output: results/figures/susceptible_sensitivity_*.png
#         results/susceptible_sensitivity_summary.csv
# ==============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)
library(mgcv)

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

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUSCEPTIBLE DEPLETION SENSITIVITY ANALYSIS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")
model_data <- readRDS("../data/model_data.rds")
stan_data  <- model_data$stan_data
dates_model <- model_data$metadata$dates_model
N_pop <- 5500000

# Load GP amplitude for reference
fit <- readRDS("../results/fit_model3.rds")
alpha_draws <- as.vector(fit$draws("alpha", format = "matrix"))
alpha_median <- median(alpha_draws)
alpha_q95 <- quantile(alpha_draws, 0.975)
cat(sprintf("  GP amplitude (alpha): median = %.3f, 97.5%% = %.3f\n",
            alpha_median, alpha_q95))

# Weekly cases (full series)
cases_all <- stan_data$cases
N <- stan_data$N
S_gi <- stan_data$S

# ==============================================================================
# 2. LOAD AND PREPARE SEROTYPE DATA
# ==============================================================================

cat("\nLoading serotype data...\n")

SEROTYPE_DIR <- "../data"
if (!file.exists(file.path(SEROTYPE_DIR, "serotype_props_2013_2016.csv"))) {
  SEROTYPE_DIR <- file.path(dirname(getwd()), "data")
}

sero_early <- read_csv(file.path(SEROTYPE_DIR, "serotype_props_2013_2016.csv"),
                       show_col_types = FALSE)
sero_late  <- read_csv(file.path(SEROTYPE_DIR, "monthly_sero_type_props_all_data.csv"),
                       show_col_types = FALSE)

# Standardize column names
if ("Month" %in% names(sero_late) && !("month" %in% names(sero_late))) {
  sero_late <- sero_late |> rename(month = Month)
}

# Extract monthly serotype proportions
sero_late_monthly <- sero_late |>
  mutate(month = as.character(month)) |>
  group_by(month) |>
  summarize(D1 = mean(D1_prop, na.rm = TRUE),
            D2 = mean(D2_prop, na.rm = TRUE),
            D3 = mean(D3_prop, na.rm = TRUE),
            D4 = mean(D4_prop, na.rm = TRUE),
            .groups = "drop") |>
  mutate(date = as.Date(paste0(month, "-15")))

# Prepare early data
if ("D1_prop" %in% names(sero_early)) {
  # Column is "Month" (capital M) with values like "2013-01"
  month_col <- if ("Month" %in% names(sero_early)) "Month" else "month"
  sero_early_monthly <- sero_early |>
    mutate(date = as.Date(paste0(.data[[month_col]], "-15"))) |>
    rename(D1 = D1_prop, D2 = D2_prop, D3 = D3_prop, D4 = D4_prop) |>
    select(date, D1, D2, D3, D4)
  sero_late_monthly <- sero_late_monthly |> select(date, D1, D2, D3, D4)
  sero_monthly <- bind_rows(sero_early_monthly, sero_late_monthly) |>
    distinct(date, .keep_all = TRUE) |>
    arrange(date)
} else {
  sero_monthly <- sero_late_monthly |> select(date, D1, D2, D3, D4)
}

# Normalize
sero_monthly <- sero_monthly |>
  mutate(total = D1 + D2 + D3 + D4,
         D1 = D1 / total, D2 = D2 / total,
         D3 = D3 / total, D4 = D4 / total) |>
  select(-total) |>
  filter(!is.na(D1))

cat(sprintf("  %d months of serotype data: %s to %s\n",
            nrow(sero_monthly),
            min(sero_monthly$date), max(sero_monthly$date)))

# Identify dominant serotype and switch dates
sero_monthly <- sero_monthly |>
  mutate(dominant = case_when(
    D1 >= pmax(D2, D3, D4) ~ "DENV-1",
    D2 >= pmax(D1, D3, D4) ~ "DENV-2",
    D3 >= pmax(D1, D2, D4) ~ "DENV-3",
    TRUE ~ "DENV-4"
  ))

switch_dates <- sero_monthly |>
  mutate(prev_dom = lag(dominant)) |>
  filter(dominant != prev_dom, !is.na(prev_dom)) |>
  pull(date)

cat(sprintf("  %d serotype switch dates identified\n", length(switch_dates)))

# Interpolate to weekly
dates_weekly <- as.Date(dates_model)
sero_weekly <- tibble(date = dates_weekly)

for (sero in c("D1", "D2", "D3", "D4")) {
  vals <- sero_monthly[[sero]]
  monthly_dates <- as.numeric(sero_monthly$date)
  weekly_dates_num <- as.numeric(dates_weekly)

  # Only interpolate within range
  in_range <- weekly_dates_num >= min(monthly_dates) &
              weekly_dates_num <= max(monthly_dates)

  interp <- rep(NA_real_, length(dates_weekly))
  if (sum(in_range) > 0) {
    interp[in_range] <- approx(monthly_dates, vals,
                                xout = weekly_dates_num[in_range],
                                rule = 2)$y
  }
  sero_weekly[[sero]] <- interp
}

# Normalize and handle NAs
sero_weekly <- sero_weekly |>
  mutate(total = D1 + D2 + D3 + D4,
         D1 = ifelse(is.na(total), 0.25, D1 / total),
         D2 = ifelse(is.na(total), 0.25, D2 / total),
         D3 = ifelse(is.na(total), 0.25, D3 / total),
         D4 = ifelse(is.na(total), 0.25, D4 / total)) |>
  select(-total)

# Dominant serotype per week
sero_weekly <- sero_weekly |>
  mutate(dominant = case_when(
    D1 >= pmax(D2, D3, D4) ~ "DENV-1",
    D2 >= pmax(D1, D3, D4) ~ "DENV-2",
    D3 >= pmax(D1, D2, D4) ~ "DENV-3",
    TRUE ~ "DENV-4"
  ))

# ==============================================================================
# 3. SENSITIVITY GRID
# ==============================================================================

cat("\nRunning sensitivity grid...\n")

# Biological parameters for birth replenishment and waning immunity
# Singapore births: ~33,000/year => ~635/week
births_per_year <- 33000
births_per_week <- births_per_year / 52
birth_rate_weekly <- births_per_week / N_pop  # as fraction of population

# Waning immunity rate from Tan et al. (2019): 0.003/year
waning_rate_yearly <- 0.003
waning_rate_weekly <- waning_rate_yearly / 52  # ~0.0000577/week

cat(sprintf("  Birth replenishment: %d births/year = %.1f/week = %.6f/week (fraction)\n",
            births_per_year, births_per_week, birth_rate_weekly))
cat(sprintf("  Waning immunity rate: %.4f/year = %.7f/week (Tan et al. 2019)\n",
            waning_rate_yearly, waning_rate_weekly))

# Focused sensitivity grid:
#   EF: {6, 8, 10, 14, 20}
#   S_init: {0.50, 0.75, 0.95}
#   Reset: {TRUE, FALSE}
#   Cross-immunity: fixed at 0 (simplicity)
#   Births + waning: toggled together {both off, both on}
# Total: 5 x 3 x 2 x 2 = 60 configurations

EF_values <- c(6, 8, 10, 14, 20)
S_init_values <- c(0.50, 0.75, 0.95)
reset_options <- c(FALSE, TRUE)
cross_immunity_years <- c(0)  # fixed for simplicity
bio_toggle <- c(FALSE, TRUE)  # births + waning together

results <- list()
counter <- 0
total_configs <- length(EF_values) * length(S_init_values) *
                 length(reset_options) * length(bio_toggle)

cat(sprintf("  Testing %d configurations...\n", total_configs))

for (ef in EF_values) {
  for (s_init in S_init_values) {
    for (do_reset in reset_options) {
      for (include_bio in bio_toggle) {
        counter <- counter + 1
        ci_years <- 0  # fixed

        include_births <- include_bio
        include_waning <- include_bio

        # --- Compute S_i[t] for each serotype ---
        n_weeks <- nrow(sero_weekly)
        S_sero <- matrix(s_init, nrow = n_weeks, ncol = 4)
        colnames(S_sero) <- c("D1", "D2", "D3", "D4")

        # Weekly serotype-specific cases
        # Use cases from the model period (after burn-in)
        # Align: cases are indexed 1..N, model dates start at week S_gi+1
        cases_model <- cases_all[(S_gi + 1):N]
        # Trim to match sero_weekly length
        n_common <- min(length(cases_model), n_weeks)

        cases_sero <- matrix(0, nrow = n_common, ncol = 4)
        for (j in 1:4) {
          sero_col <- c("D1", "D2", "D3", "D4")[j]
          cases_sero[, j] <- round(cases_model[1:n_common] *
                                     sero_weekly[[sero_col]][1:n_common])
        }

        # Track depletion with births and waning immunity
        for (t in 2:n_common) {
          for (j in 1:4) {
            infections_j <- cases_sero[t - 1, j] * ef
            depletion_j <- infections_j / N_pop
            S_sero[t, j] <- S_sero[t - 1, j] - depletion_j

            # Cross-immunity: infection with serotype j also depletes others temporarily
            if (ci_years > 0) {
              ci_weeks <- ci_years * 52
              if (t > ci_weeks) {
                # Recover cross-immunity from ci_weeks ago
                for (k in 1:4) {
                  if (k != j) {
                    recovered <- cases_sero[t - ci_weeks, j] * ef / N_pop * 0.5
                    S_sero[t, k] <- S_sero[t, k] + recovered
                  }
                }
              }
              # Apply cross-immunity depletion
              for (k in 1:4) {
                if (k != j) {
                  S_sero[t, k] <- S_sero[t, k] - infections_j / N_pop * 0.5
                }
              }
            }
          }

          # Birth replenishment: new births are susceptible to all serotypes
          if (include_births) {
            for (j in 1:4) {
              S_sero[t, j] <- S_sero[t, j] + birth_rate_weekly
            }
          }

          # Waning immunity: a fraction of immune individuals become susceptible again
          if (include_waning) {
            for (j in 1:4) {
              S_sero[t, j] <- S_sero[t, j] + waning_rate_weekly * (1 - S_sero[t, j])
            }
          }

          # Reset at serotype switches?
          if (do_reset) {
            current_date <- dates_weekly[t]
            if (any(abs(as.numeric(current_date - switch_dates)) < 7)) {
              # Find new dominant serotype
              new_dom_idx <- which.max(c(sero_weekly$D1[t], sero_weekly$D2[t],
                                          sero_weekly$D3[t], sero_weekly$D4[t]))
              # Reset the NEW dominant serotype's susceptible to s_init
              S_sero[t, new_dom_idx] <- s_init
            }
          }

          # Floor
          S_sero[t, ] <- pmax(0.01, pmin(1.0, S_sero[t, ]))
        }

        # Effective susceptible fraction (weighted by serotype proportion)
        S_eff <- numeric(n_common)
        for (t in 1:n_common) {
          props <- c(sero_weekly$D1[t], sero_weekly$D2[t],
                     sero_weekly$D3[t], sero_weekly$D4[t])
          S_eff[t] <- sum(props * S_sero[t, 1:4])
        }

        # Also compute S for dominant serotype only
        S_dominant <- numeric(n_common)
        for (t in 1:n_common) {
          dom_idx <- which.max(c(sero_weekly$D1[t], sero_weekly$D2[t],
                                 sero_weekly$D3[t], sero_weekly$D4[t]))
          S_dominant[t] <- S_sero[t, dom_idx]
        }

        # Key metrics
        log_S_eff_range <- diff(range(log(S_eff), na.rm = TRUE))
        log_S_dom_range <- diff(range(log(S_dominant), na.rm = TRUE))
        S_eff_range <- diff(range(S_eff, na.rm = TRUE))
        S_dom_range <- diff(range(S_dominant, na.rm = TRUE))

        # Can this explain the GP amplitude?
        # GP alpha ≈ 0.30 on log scale
        explains_gp <- log_S_dom_range >= alpha_median

        results[[counter]] <- tibble(
          EF = ef,
          S_init = s_init,
          reset_at_switch = do_reset,
          cross_immunity_years = ci_years,
          include_births = include_births,
          include_waning = include_waning,
          S_eff_min = min(S_eff, na.rm = TRUE),
          S_eff_max = max(S_eff, na.rm = TRUE),
          S_eff_range = S_eff_range,
          S_dom_min = min(S_dominant, na.rm = TRUE),
          S_dom_max = max(S_dominant, na.rm = TRUE),
          S_dom_range = S_dom_range,
          log_S_eff_range = log_S_eff_range,
          log_S_dom_range = log_S_dom_range,
          gp_alpha_median = alpha_median,
          ratio_to_gp = log_S_dom_range / alpha_median,
          explains_gp = explains_gp
        )
      }
    }
  }
}

results_df <- bind_rows(results)
cat(sprintf("  %d configurations tested.\n", nrow(results_df)))

# ==============================================================================
# 4. SUMMARY
# ==============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("RESULTS SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat(sprintf("GP amplitude (alpha) median: %.3f\n\n", alpha_median))

cat("Configurations where log(S_dominant) range >= GP alpha:\n")
explains <- results_df |> filter(explains_gp)
cat(sprintf("  %d / %d (%.0f%%) configurations\n\n",
            nrow(explains), nrow(results_df),
            100 * nrow(explains) / nrow(results_df)))

if (nrow(explains) > 0) {
  cat("These configurations require:\n")
  cat(sprintf("  EF range: %d - %d\n", min(explains$EF), max(explains$EF)))
  cat(sprintf("  S_init range: %.2f - %.2f\n",
              min(explains$S_init), max(explains$S_init)))
  cat(sprintf("  Reset at switch: %s\n",
              paste(unique(explains$reset_at_switch), collapse = ", ")))
  cat(sprintf("  Births+waning: %s\n",
              paste(unique(explains$include_births), collapse = ", ")))
}

cat("\n--- By reset option ---\n")
results_df |>
  group_by(reset_at_switch) |>
  summarize(
    mean_log_S_dom_range = mean(log_S_dom_range),
    max_log_S_dom_range = max(log_S_dom_range),
    pct_explains = mean(explains_gp) * 100,
    .groups = "drop"
  ) |>
  print()

cat("\n--- By EF ---\n")
results_df |>
  group_by(EF) |>
  summarize(
    mean_log_S_dom_range = mean(log_S_dom_range),
    max_log_S_dom_range = max(log_S_dom_range),
    pct_explains = mean(explains_gp) * 100,
    .groups = "drop"
  ) |>
  print()

cat("\n--- By S_init ---\n")
results_df |>
  group_by(S_init) |>
  summarize(
    mean_log_S_dom_range = mean(log_S_dom_range),
    max_log_S_dom_range = max(log_S_dom_range),
    pct_explains = mean(explains_gp) * 100,
    .groups = "drop"
  ) |>
  print()

cat("\n--- By births + waning ---\n")
results_df |>
  group_by(include_births, include_waning) |>
  summarize(
    mean_log_S_dom_range = mean(log_S_dom_range),
    max_log_S_dom_range = max(log_S_dom_range),
    pct_explains = mean(explains_gp) * 100,
    .groups = "drop"
  ) |>
  print()

cat("\n--- By EF x births+waning ---\n")
results_df |>
  group_by(EF, include_births) |>
  summarize(
    mean_log_S_dom_range = mean(log_S_dom_range),
    max_log_S_dom_range = max(log_S_dom_range),
    pct_explains = mean(explains_gp) * 100,
    .groups = "drop"
  ) |>
  mutate(bio_label = ifelse(include_births, "Births+Waning ON", "Births+Waning OFF")) |>
  print()

# Save
write_csv(results_df, "../results/susceptible_sensitivity_summary.csv")
cat("\nSaved: results/susceptible_sensitivity_summary.csv\n")

# ==============================================================================
# 5. FIGURES
# ==============================================================================

cat("\nGenerating figures...\n")

# Figure 1: Heatmap — log(S_dom) range vs GP alpha, by EF × S_init
# Faceted by reset × cross-immunity

fig1_data <- results_df |>
  mutate(
    reset_label = ifelse(reset_at_switch, "Reset S at switch", "No reset"),
    bio_label = ifelse(include_births, "Births + Waning ON", "Births + Waning OFF"),
    panel = paste(reset_label, bio_label, sep = "\n"),
    S_init_label = sprintf("S₀ = %.2f", S_init)
  )

fig1 <- ggplot(fig1_data,
               aes(x = factor(EF), y = S_init_label, fill = ratio_to_gp)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", ratio_to_gp * 100)),
            size = 2.5) +
  scale_fill_gradient2(
    low = "#d73027", mid = "#ffffbf", high = "#1a9850",
    midpoint = 1.0,
    name = "Ratio to GP\namplitude",
    limits = c(0, max(fig1_data$ratio_to_gp, 2))
  ) +
  facet_wrap(~panel, ncol = 2) +
  labs(
    title = "Can susceptible depletion explain the GP residual?",
    subtitle = sprintf("Green (≥100%%): log(S) range matches GP amplitude (α = %.2f). Red (<100%%): too small.",
                        alpha_median),
    x = "Expansion Factor (EF)",
    y = "Initial susceptible fraction (S₀)"
  )

ggsave("../results/figures/susceptible_sensitivity_heatmap.png", fig1,
       width = 12, height = 8, dpi = 200)
cat("  Saved: results/figures/susceptible_sensitivity_heatmap.png\n")

# Figure 2: Time series of S_dominant for key configurations
# Pick 4 representative configs: (no reset, EF=8), (reset, EF=8),
#                                 (no reset, EF=20), (reset, EF=20)

key_configs <- list(
  list(ef = 8,  s_init = 0.75, reset = FALSE, ci = 0, births = FALSE, waning = FALSE,
       label = "EF=8, no reset, no bio"),
  list(ef = 8,  s_init = 0.75, reset = TRUE,  ci = 0, births = FALSE, waning = FALSE,
       label = "EF=8, reset, no bio"),
  list(ef = 14, s_init = 0.75, reset = TRUE,  ci = 0, births = FALSE, waning = FALSE,
       label = "EF=14, reset, no bio"),
  list(ef = 14, s_init = 0.75, reset = TRUE,  ci = 0, births = TRUE, waning = TRUE,
       label = "EF=14, reset, births+waning"),
  list(ef = 20, s_init = 0.95, reset = TRUE,  ci = 0, births = FALSE, waning = FALSE,
       label = "EF=20, reset, no bio"),
  list(ef = 20, s_init = 0.95, reset = TRUE,  ci = 0, births = TRUE, waning = TRUE,
       label = "EF=20, reset, births+waning")
)

ts_data <- list()
for (cfg in key_configs) {
  n_weeks <- nrow(sero_weekly)
  S_sero <- matrix(cfg$s_init, nrow = n_weeks, ncol = 4)
  n_common <- min(length(cases_model), n_weeks)

  cases_sero <- matrix(0, nrow = n_common, ncol = 4)
  for (j in 1:4) {
    sero_col <- c("D1", "D2", "D3", "D4")[j]
    cases_sero[, j] <- round(cases_model[1:n_common] *
                               sero_weekly[[sero_col]][1:n_common])
  }

  for (t in 2:n_common) {
    for (j in 1:4) {
      infections_j <- cases_sero[t - 1, j] * cfg$ef
      S_sero[t, j] <- S_sero[t - 1, j] - infections_j / N_pop

      if (cfg$ci > 0) {
        ci_weeks <- cfg$ci * 52
        if (t > ci_weeks) {
          for (k in 1:4) {
            if (k != j) {
              recovered <- cases_sero[t - ci_weeks, j] * cfg$ef / N_pop * 0.5
              S_sero[t, k] <- S_sero[t, k] + recovered
            }
          }
        }
        for (k in 1:4) {
          if (k != j) {
            S_sero[t, k] <- S_sero[t, k] - infections_j / N_pop * 0.5
          }
        }
      }
    }

    # Birth replenishment
    if (isTRUE(cfg$births)) {
      for (j in 1:4) {
        S_sero[t, j] <- S_sero[t, j] + birth_rate_weekly
      }
    }

    # Waning immunity
    if (isTRUE(cfg$waning)) {
      for (j in 1:4) {
        S_sero[t, j] <- S_sero[t, j] + waning_rate_weekly * (1 - S_sero[t, j])
      }
    }

    if (cfg$reset) {
      current_date <- dates_weekly[t]
      if (any(abs(as.numeric(current_date - switch_dates)) < 7)) {
        new_dom_idx <- which.max(c(sero_weekly$D1[t], sero_weekly$D2[t],
                                    sero_weekly$D3[t], sero_weekly$D4[t]))
        S_sero[t, new_dom_idx] <- cfg$s_init
      }
    }

    S_sero[t, ] <- pmax(0.01, pmin(1.0, S_sero[t, ]))
  }

  S_dom <- numeric(n_common)
  for (t in 1:n_common) {
    dom_idx <- which.max(c(sero_weekly$D1[t], sero_weekly$D2[t],
                           sero_weekly$D3[t], sero_weekly$D4[t]))
    S_dom[t] <- S_sero[t, dom_idx]
  }

  ts_data[[cfg$label]] <- tibble(
    date = dates_weekly[1:n_common],
    S_dominant = S_dom,
    log_S = log(S_dom),
    config = cfg$label
  )
}

ts_df <- bind_rows(ts_data)

fig2 <- ggplot(ts_df, aes(x = date, y = log_S, color = config)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = switch_dates, linetype = "dashed",
             color = "grey50", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey30") +
  annotate("rect",
           ymin = -alpha_median, ymax = alpha_median,
           xmin = min(ts_df$date), xmax = max(ts_df$date),
           fill = "blue", alpha = 0.1) +
  annotate("text", x = min(ts_df$date) + 100, y = alpha_median + 0.02,
           label = sprintf("GP amplitude band (±α = ±%.2f)", alpha_median),
           size = 3, hjust = 0, color = "blue") +
  labs(
    title = "log(S_dominant) trajectory under different assumptions",
    subtitle = "Blue band = GP amplitude. If log(S) variation reaches the band, depletion COULD explain the GP.",
    x = NULL, y = "log(S_dominant)",
    color = "Configuration"
  ) +
  theme(legend.position = "bottom",
        legend.direction = "vertical")

ggsave("../results/figures/susceptible_sensitivity_trajectories.png", fig2,
       width = 12, height = 7, dpi = 200)
cat("  Saved: results/figures/susceptible_sensitivity_trajectories.png\n")

# Figure 3: Bar chart — what % of GP amplitude can depletion explain?
fig3_data <- results_df |>
  mutate(
    reset_label = ifelse(reset_at_switch, "Reset", "No reset"),
    bio_label = ifelse(include_births, "+Bio", "No bio"),
    config_short = paste(reset_label, bio_label)
  ) |>
  group_by(EF, config_short) |>
  summarize(
    mean_ratio = mean(ratio_to_gp),
    max_ratio = max(ratio_to_gp),
    .groups = "drop"
  )

fig3 <- ggplot(fig3_data, aes(x = factor(EF), y = max_ratio * 100,
                               fill = config_short)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red",
             linewidth = 0.8) +
  annotate("text", x = 0.5, y = 105, label = "100% = fully explains GP",
           hjust = 0, size = 3, color = "red") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Maximum % of GP amplitude explained by susceptible depletion",
    subtitle = "Across all S₀ values. Dashed line = GP amplitude fully explained.",
    x = "Expansion Factor", y = "% of GP amplitude explained",
    fill = "Configuration"
  ) +
  theme(legend.position = "bottom")

ggsave("../results/figures/susceptible_sensitivity_barplot.png", fig3,
       width = 10, height = 6, dpi = 200)
cat("  Saved: results/figures/susceptible_sensitivity_barplot.png\n")

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SENSITIVITY ANALYSIS COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Key question: Can national-level susceptible depletion explain\n")
cat("the GP residual amplitude (α = ", sprintf("%.3f", alpha_median), ")?\n\n")

n_explains <- sum(results_df$explains_gp)
cat(sprintf("Answer: %d / %d configurations (%.0f%%) reach the GP amplitude.\n",
            n_explains, nrow(results_df),
            100 * n_explains / nrow(results_df)))

if (n_explains == 0) {
  cat("\n=> Under NO tested assumptions does national-level susceptible\n")
  cat("   depletion match the GP residual amplitude. The Rt fluctuations\n")
  cat("   must be driven by factors beyond aggregate immunity dynamics.\n")
  cat("   Note: birth replenishment and waning immunity were tested but\n")
  cat("   do not change this conclusion.\n")
} else if (n_explains < nrow(results_df) / 2) {
  cat("\n=> Only under aggressive assumptions (high EF, reset at switch,\n")
  cat("   high initial susceptibility) does depletion approach the GP\n")
  cat("   amplitude. The conclusion is sensitive to these assumptions.\n")
  cat("   Birth replenishment (~635/week) and waning immunity (0.003/yr)\n")
  cat("   partially offset depletion, reducing the range of log(S).\n")
} else {
  cat("\n=> Under most assumptions, susceptible depletion CAN explain\n")
  cat("   the GP residual amplitude.\n")
  cat("   Birth replenishment and waning immunity have a modest effect.\n")
}

cat("\nOutput files:\n")
cat("  results/figures/susceptible_sensitivity_heatmap.png\n")
cat("  results/figures/susceptible_sensitivity_trajectories.png\n")
cat("  results/figures/susceptible_sensitivity_barplot.png\n")
cat("  results/susceptible_sensitivity_summary.csv\n")
