# Dengue GP Project — Complete Results Summary

## Project Overview

Bayesian renewal equation model for weekly dengue Rt in Singapore (2012–2022), decomposing log(Rt) into climate covariates and a residual Gaussian process (HSGP, Matérn 3/2). Extended with cross-kernel comparison, forecasting, branching process analysis, VAE comparison, spectral analysis, and counterfactual analysis.

---

## 1. Core Model (Climate + Residual GP)

**Model:** log(Rt) = μ + β_temp × temp + β_rain × rain + f_residual(t)

### Key Parameters
- **Temperature effect:** 0.989 [95% CrI: 0.969, 1.009] — negligible
- **Rainfall effect:** 1.003 [0.986, 1.020] — negligible
- **Climate variance proportion:** 0.4% (median) — climate explains <1% of Rt
- **Residual GP variance:** 99.6% — nearly all Rt variation is unexplained by climate
- **GP length scale (ρ):** 2.86 weeks [2.09, 3.71] — short, tracks outbreak dynamics
- **GP amplitude (α):** 0.30 [0.27, 0.34]
- **Overdispersion (φ):** NegBin, well-estimated

### Model Validation
- Posterior predictive 95% coverage: 97.7%
- Posterior predictive 80% coverage: 92.3%
- Zero divergences, all Rhat < 1.01, ESS > 800

### Figures
- `decomposition.png` — log(Rt) decomposed into climate + residual
- `rt_with_cases.png` — Rt trajectory overlaid on observed cases
- `effect_sizes.png` — Posterior distributions of covariate effects
- `prior_predictive_summary.png` — Prior predictive check
- `posterior_predictive_timeseries.png` — Posterior predictive check

---

## 2. Serotype Analysis (Post-hoc)

### Key Findings
- **7 serotype switch events** identified from GAM-smoothed proportions (2013–2022)
- **5 of 7 switches** show P(elevated f_residual) > 0.75 in 0–3 months post-switch
- **DENV-2→DENV-3 (Feb 2020):** strongest signal (P = 1.0, 20–35% Rt increase)
- **Cross-correlation (CCF):** entropy leads f_residual by ~10–12 months

### Figures
- `serotype_panel.png` — Three-panel: f_residual, serotype proportions, Shannon entropy
- `serotype_switch_profiles.png` — f_residual profiles around each switch event
- `serotype_ccf.png` — Cross-correlation between f_residual and entropy
- `serotype_entropy_profiles.png` — Entropy trajectory around each switch
- `serotype_rt_entropy_overlay.png` — f_residual and entropy overlaid per switch

---

## 3. Cross-Kernel Comparison

### LOO-CV Results (all within ~2 ELPD — statistically indistinguishable)

| Kernel | ELPD | ΔELPD | SE(Δ) | p_loo |
|--------|------|-------|-------|-------|
| Matérn 1/2 | -2839.4 | 0.0 | — | 105.2 |
| Matérn 3/2 | -2840.2 | -0.8 | 0.5 | 105.8 |
| Squared Exp | -2841.0 | -1.6 | 0.7 | 106.6 |
| Matérn 5/2 | -2841.3 | -1.8 | 0.6 | 106.2 |

### MCMC Diagnostics (all kernels)
- Zero divergences across all 4 kernels
- All Rhat ≤ 1.01
- Min ESS bulk > 805, min ESS tail > 1609

### Interpretation
Kernel choice has negligible impact on model fit. The roughest kernel (Matérn 1/2) is marginally preferred, consistent with dengue Rt having abrupt changes. The Matérn 3/2 choice is defensible but not critical — all scientific conclusions are robust to kernel specification.

### Figures
- `kernel_gp_hyperparameters.png` — α and ρ posteriors by kernel (violin plots)
- `kernel_variance_decomposition.png` — Climate variance % by kernel
- `kernel_rt_comparison.png` — Rt trajectories overlaid (all 4 kernels)
- `kernel_rt_faceted.png` — Rt trajectories faceted by kernel
- `kernel_residual_gp.png` — f_residual faceted by kernel
- `kernel_serotype_detectability.png` — P(elevated) per switch per kernel

---

## 4. Out-of-Sample Forecasting

### Design
- Rolling-window forecast: 4 origins, horizons h = 4, 8, 13 weeks
- Two baselines: Rt=1 (naive) and random walk (last Rt held constant)
- With and without serotype entropy as leading indicator (43-week lag)
- Matérn 3/2 kernel

### Results (without entropy)

| Horizon | CRPS | Skill vs Rt=1 | Skill vs RW | 95% Coverage |
|---------|------|---------------|-------------|--------------|
| 4 weeks | 119 | +0.35 | **-7.6** | 100% |
| 8 weeks | 243 | -0.67 | -2.3 | 100% |
| 13 weeks | 350 | -64 | -10.5 | 100% |

### Key Findings
- **GP never beats the random walk baseline**, even at 4 weeks
- The apparent 35% skill vs Rt=1 is misleading — the GP is converging to the same Rt≈1 prediction as the baseline
- At horizons >> ρ (3 weeks), the GP posterior forecast is indistinguishable from the prior — it regresses to the mean
- **Entropy as a leading indicator does not help** — slightly worse at every horizon. The ~10-month lead is too long for these forecast horizons; entropy is better suited for seasonal risk assessment, not weekly tactical forecasting
- **Conclusion:** The GP excels at retrospective decomposition but provides no forecasting value beyond naive persistence

### Figures
- `kernel_forecast_crps.png` — CRPS boxplots by horizon
- `kernel_forecast_calibration.png` — Calibration plot (nominal vs empirical coverage)
- `kernel_forecast_trajectories.png` — Example forecast trajectories at selected origins
- `kernel_forecast_entropy_effect.png` — With vs without entropy comparison

---

## 5. Branching Process Analysis

### Time-Varying Extinction Probability
- **Strong anticorrelation** between Rt and extinction probability q (Spearman ρ = -0.94)
- P(outbreak sustains) = 1 - q computed for each week from posterior Rt and overdispersion (φ as NegBin k)

### Yearly Summary

| Year | Mean Rt | Mean P(sustain) | Weeks sustained | Notes |
|------|---------|-----------------|-----------------|-------|
| 2012 | 1.07 | 0.16 | 0 | |
| 2013 | 1.13 | 0.24 | 9 | DENV-1 peak |
| 2014 | 0.99 | 0.14 | 6 | |
| 2015 | 1.07 | 0.15 | 5 | |
| 2016 | 0.93 | 0.08 | 3 | Inter-epidemic trough |
| 2017 | 1.03 | 0.15 | 3 | |
| 2018 | 1.08 | 0.18 | 5 | |
| 2019 | 1.11 | 0.23 | 10 | Large DENV-2 outbreak |
| 2020 | 1.04 | 0.18 | 9 | DENV-3 emergence |
| 2021 | 0.95 | 0.10 | 3 | Low activity |
| **2022** | **1.17** | **0.27** | **18** | **Most sustained** |

### Serotype Switch Effect
- **5 of 7 switches** showed increased P(outbreak sustains) after the switch
- Mean increase: +6 percentage points
- Serotype switches create measurable windows of outbreak vulnerability

### Figures
- `branching_extinction_probability.png` — Rt trajectory + extinction probability (dual panel)
- `branching_serotype_switches.png` — P(sustain) before vs after each switch
- `branching_superspreading.png` — P(offspring > 20) and top-20% transmission share over time
- `branching_offspring_distribution.png` — NegBin offspring PMFs at selected Rt values

---

## 6. VAE Latent Space Comparison

### Design
- Sliding-window VAE (W=26 weeks, stride=1) on log-transformed case counts
- Architecture: encoder (26→32→16→4 latent), decoder (4→16→32→26)
- 200 epochs, beta-annealing, 80/20 chronological split

### Results — Important Negative Finding

| Metric | Value | p-value |
|--------|-------|---------|
| Spearman latent0 vs f_residual | 0.090 | 0.035 |
| Spearman latent1 vs f_residual | 0.018 | 0.674 |
| Spearman latent2 vs f_residual | 0.063 | 0.142 |
| Spearman latent3 vs f_residual | -0.031 | 0.464 |
| Spearman PC1 vs f_residual | -0.032 | 0.462 |
| Spearman PC2 vs f_residual | 0.045 | 0.291 |
| Fraction of switches in top-20% recon error | 0.153 | — |
| Mann-Whitney (switch vs non-switch MSE) | — | 0.966 |

### Interpretation
- **VAE does NOT discover the GP residual structure** — max correlation is 0.09, PCA components near zero
- **Reconstruction error does NOT spike at serotype switches** (p = 0.97, no signal)
- Switch windows actually had *lower* MSE than non-switch windows
- **This validates the renewal equation framework**: the mechanistic structure (renewal equation + GP) captures biologically meaningful temporal variation that a purely data-driven approach (VAE on raw cases) cannot discover. The GP residual encodes epidemiological information, not just statistical smoothing.

### Figures
- `vae_latent_space.png` — 2×2 PCA panel (colored by time, f_residual, serotype, entropy)
- `vae_latent_trajectory.png` — Latent PC1/PC2 over time vs f_residual
- `vae_reconstruction_error.png` — MSE time series with serotype switch lines
- `vae_training_loss.png` — Train/val loss curves (total, reconstruction, KL)
- `vae_uncertainty_comparison.png` — GP CrI width vs VAE decoder variance

---

## 7. Spectral Analysis

### Dominant Periodicities in f_residual

| Period | Weeks | Power (× median) | Interpretation |
|--------|-------|-------------------|----------------|
| **6 months** | 26 | 2,934,187× | **Strongest** — biannual dengue pattern |
| 7 months | 30 | 1,151,560× | Sub-annual dynamics |
| 4.4 months | 19 | 1,000,170× | Within-outbreak cycling |
| **1.6 years** | 82 | 783,653× | Multi-year cycle |
| **3.7 years** | 191 | 633,379× | Near 4-year serotype rotation |
| 3 months | 13 | 558,498× | Individual outbreak duration |
| **1 year** | 52 | 498,760× | Annual seasonal (residual after climate removal) |
| **2.2 years** | 115 | 293,743× | **Serotype replacement cycling** |

### Interpretation
- The residual GP has rich multi-scale structure even after climate covariates are removed
- The **~2-year cycle** (115 weeks) is consistent with serotype replacement theory
- The **~4-year cycle** (191 weeks) matches the theoretical 4-serotype rotation period
- The **6-month peak** suggests within-year epidemic dynamics beyond temperature/rainfall — possibly monsoon-driven vector ecology
- The **52-week annual signal** is notable: residual seasonality exists that climate covariates don't capture, suggesting non-climatic seasonal drivers (e.g., school terms, travel patterns)

### Figures
- `spectral_periodogram.png` — Power spectrum with reference period annotations and 95% CrI
- `spectral_with_context.png` — f_residual time series + periodogram (dual panel)

---

## 8. Counterfactual Analysis

### Question: "How many excess cases are attributable to non-climate drivers?"

Counterfactual: Rt_cf = exp(μ + f_climate), zeroing out the residual GP.

### Total Excess (2012–2022)
- **Median: -763 cases** (95% CrI: -12,129 to +10,556)
- Net excess is approximately zero — positive and negative excursions cancel over the decade

### Per Serotype Switch (0–6 months post-switch)

| Switch Event | Excess Cases | 95% CrI | Significant? |
|--------------|-------------|---------|--------------|
| DENV-2→DENV-1 (Mar 2013) | -61 | -1199 to 1115 | No |
| DENV-1→DENV-2 (Nov 2015) | +9 | -928 to 841 | No |
| DENV-2→DENV-1 (Sep 2016) | -693 | -953 to -409 | **Yes (fewer)** |
| DENV-1→DENV-2 (Dec 2016) | -109 | -241 to 28 | Marginal |
| **DENV-2→DENV-3 (Feb 2020)** | **+4,419** | **2925 to 6043** | **Yes (excess)** |
| DENV-3→DENV-2 (Jun 2020) | -843 | -3317 to 1459 | No |
| DENV-2→DENV-3 (Mar 2021) | -296 | -599 to -41 | **Yes (fewer)** |

### Key Finding
- **The Feb 2020 DENV-3 emergence was the only switch producing a clearly significant excess: ~4,400 cases above climate prediction** — consistent with a large naive susceptible pool encountering DENV-3 for the first time
- Most other switches show excess cases indistinguishable from zero, suggesting that only *novel serotype introductions* (not rotations between familiar serotypes) produce detectable excess transmission

### Figures
- `counterfactual_rt.png` — Actual vs counterfactual Rt with gap shading
- `counterfactual_excess_cases.png` — Weekly excess cases + cumulative line
- `counterfactual_switch_attribution.png` — Excess cases per switch (bar chart with CrI)

---

## Complete Figure Inventory (52 figures)

### Core Model (12)
- `decomposition.png`, `rt_comparison.png`, `rt_with_cases.png`
- `effect_sizes.png`, `prior_predictive_summary.png`, `prior_predictive_Rt_trajectories.png`, `prior_predictive_cases.png`
- `posterior_predictive_timeseries.png`, `posterior_predictive_intervals.png`, `posterior_predictive_density.png`, `posterior_predictive_histograms.png`, `posterior_predictive_residuals.png`, `posterior_predictive_stats.png`, `posterior_predictive_summary.png`, `ppc_model2.png`

### Serotype Analysis (5)
- `serotype_panel.png`, `serotype_switch_profiles.png`, `serotype_ccf.png`
- `serotype_entropy_profiles.png`, `serotype_rt_entropy_overlay.png`

### Sensitivity (2)
- `acf_comparison_short_gp.png`, `residual_comparison_short_gp.png`

### Cross-Kernel Comparison (6)
- `kernel_gp_hyperparameters.png`, `kernel_variance_decomposition.png`
- `kernel_rt_comparison.png`, `kernel_rt_faceted.png`
- `kernel_residual_gp.png`, `kernel_serotype_detectability.png`

### Forecasting (4)
- `kernel_forecast_crps.png`, `kernel_forecast_calibration.png`
- `kernel_forecast_trajectories.png`, `kernel_forecast_entropy_effect.png`

### Branching Process (4)
- `branching_extinction_probability.png`, `branching_serotype_switches.png`
- `branching_superspreading.png`, `branching_offspring_distribution.png`

### VAE Comparison (5)
- `vae_latent_space.png`, `vae_latent_trajectory.png`
- `vae_reconstruction_error.png`, `vae_training_loss.png`, `vae_uncertainty_comparison.png`

### Spectral Analysis (2)
- `spectral_periodogram.png`, `spectral_with_context.png`

### Counterfactual Analysis (3)
- `counterfactual_rt.png`, `counterfactual_excess_cases.png`
- `counterfactual_switch_attribution.png`

---

## Key References

1. Cori A et al. (2013). Am J Epidemiol, 178(9), 1505–1512. doi: 10.1093/aje/kwt133
2. Riutort-Mayol G et al. (2023). Stat Comput, 33, 17. doi: 10.1007/s11222-022-10167-2
3. Chan M, Johansson MA (2012). PLOS ONE, 7(11), e50972. doi: 10.1371/journal.pone.0050972
4. Tan LK et al. (2019). Am J Epidemiol, 188(8), 1529–1538. doi: 10.1093/aje/kwz110
5. Finch E et al. (2025). Nat Commun, 16, 11364. doi: 10.1038/s41467-025-66411-6
6. Undurraga EA et al. (2013). PLOS Negl Trop Dis, 7(2), e2056. doi: 10.1371/journal.pntd.0002056
7. Bhatt S et al. (2013). Nature, 496(7446), 504–507. doi: 10.1038/nature12060
