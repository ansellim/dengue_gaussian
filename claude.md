# Project Memory: Dengue Rt Estimation with Gaussian Processes

## Project Overview
Bayesian epidemiological model for estimating the effective reproduction number (Rt) of dengue in Singapore (2012-2022) using:
- Renewal equation framework
- Hilbert Space Gaussian Process (HSGP) approximation with Matern 3/2 kernel
- Semiparametric decomposition: log(Rt) = mu + f_climate + f_residual
- Negative binomial observation model

## Current Status (Updated March 2026)

### What Has Been Done
1. **Data Pipeline**: Complete acquisition and preparation of dengue cases, weather, NPI (OxCGRT), Wolbachia coverage, and serotype data
2. **Model Development**: Multiple model versions tested (v0-v3); current model is Model 3 (climate-only: temp + rain + residual GP)
3. **Model Simplification**: Removed Wolbachia and NPI covariates due to confounding (NPI-serotype) and limited signal (~30% coverage)
4. **Model Validation**: Prior predictive checks, posterior predictive checks, residual diagnostics
5. **Sensitivity Analysis**: Tight GP amplitude prior to test variance decomposition robustness
6. **Serotype Early Warning**: Cross-correlation, pre-switch Rt dynamics, entropy profiling

### Model Evolution
- **Model 0**: Single GP baseline (no covariates)
- **Model 1**: Climate covariates + residual GP
- **Model 2**: Full model (climate + Wolbachia + NPI + residual GP) — retired due to confounding
- **Model 3**: Climate-only (temp + rain + residual GP) — **CURRENT**

### Key Findings (from previous Model 2, retained for reference)
| Parameter | Estimate | 95% CI | Interpretation |
|-----------|----------|--------|----------------|
| Temperature | 0.99 | [0.96, 1.02] | Negligible effect at weekly resolution |
| Rainfall | 1.01 | [0.99, 1.04] | Negligible effect at weekly resolution |

### Why Wolbachia/NPI Were Removed
1. **NPI-Serotype Confounding**: Cannot separate NPI effects from DENV-3 emergence (both in 2020)
2. **Wolbachia Coverage Uncertainty**: Only ~30% coverage by 2022, effect diluted at national level
3. **Residual GP Dominates**: ~90% of variance in residual GP, suggesting unmeasured drivers
4. **Focus Shift**: Model now serves as platform for serotype early warning analysis

---

## File Structure

### Code Files

| File | Purpose | Status |
|------|---------|--------|
| `01_acquire_data.py` | Download dengue, weather, NPI, Wolbachia, serotype data | Complete |
| `02_prepare_model_data.R` | Merge, lag, standardize, format for Stan (climate covariates only) | Complete |
| `03_prior_predictive.R` | Verify priors generate plausible data | Complete |
| `05_model3_climate_only.stan` | Model 3: Climate-only (temp + rain) + residual GP — **CURRENT** | Complete |
| `05_model3_climate_only_tight_gp.stan` | Model 3 variant: Tighter GP amplitude prior | Complete |
| `07_postprocess.R` | Posterior summaries, decomposition, effect sizes | Complete |
| `08_serotype_analysis.R` | Serotype early warning: CCF, pre-switch dynamics, entropy | Complete |
| `09_posterior_predictive.R` | Model validation via posterior predictive checks | Complete |
| `13_fit_model3.R` | Fit Model 3 (climate-only) | Complete |
| `14_fit_tight_gp.R` | Sensitivity: tight GP amplitude prior | Complete |

#### Legacy Files (Model 2, retained for reference)
| File | Purpose |
|------|---------|
| `03_model0_baseline.stan` | Model 0: Single GP baseline |
| `04_model1_climate.stan` | Model 1: Climate covariates + residual GP |
| `05_model2_full.stan` | Model 2: Full model (original, rho~15 weeks) |
| `05_model2_full_short_gp.stan` | Model 2: Short GP (rho~6 weeks) |
| `05_model2_no_npi.stan` | Model 2: No NPI covariate |
| `06_fit_models.R` | Fit Models 0, 1, 2 with comparison |
| `10_fit_short_gp.R` | Fit and compare short GP model |
| `11_sensitivity_pre2020.R` | Pre-2020 sensitivity analysis |
| `12_dose_response.R` | Wolbachia dose-response extrapolation |

### Data Files

| File | Description |
|------|-------------|
| `data/model_data.rds` | Prepared Stan data with climate covariates |
| `data/raw_dengue_cases.csv` | Weekly case counts (2012-2022) |
| `data/raw_weather.csv` | Temperature and rainfall from Meteostat |
| `data/raw_npi_oxcgrt.csv` | OxCGRT Stringency Index for Singapore |
| `data/raw_wolbachia.csv` | Wolbachia coverage estimates |
| `data/monthly_sero_type_props_all_data.csv` | Monthly serotype proportions |

### Results Files

| File | Description |
|------|-------------|
| `results/fit_model3.rds` | **Current best model fit** |
| `results/fit_model3_tight_gp.rds` | Tight GP sensitivity analysis fit |
| `results/fit_model2_short_gp.rds` | Legacy Model 2 fit |

---

## Analysis Pipeline

```
┌─────────────────────────────────────────────────────────────────────┐
│ 01_acquire_data.py                                                   │
│   └── Downloads: dengue, weather, NPI, Wolbachia, serotype          │
├─────────────────────────────────────────────────────────────────────┤
│ 02_prepare_model_data.R                                              │
│   └── Creates: data/model_data.rds (climate covariates only)        │
├─────────────────────────────────────────────────────────────────────┤
│ 03_prior_predictive.R                                                │
│   └── Validates: priors generate plausible Rt and cases             │
├─────────────────────────────────────────────────────────────────────┤
│ 13_fit_model3.R                                                      │
│   └── Creates: results/fit_model3.rds                               │
├─────────────────────────────────────────────────────────────────────┤
│ 07_postprocess.R                                                     │
│   └── Creates: effect estimates, decomposition plots                │
├─────────────────────────────────────────────────────────────────────┤
│ 09_posterior_predictive.R                                            │
│   └── Validates: model fits data, checks residual autocorrelation   │
├─────────────────────────────────────────────────────────────────────┤
│ 08_serotype_analysis.R                                               │
│   └── Serotype early warning: CCF, entropy, pre-switch dynamics     │
├─────────────────────────────────────────────────────────────────────┤
│ 14_fit_tight_gp.R                                                    │
│   └── Sensitivity: tight GP amplitude prior                         │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Model Specification

### Current Model: `05_model3_climate_only.stan`

```
log(Rt) = mu + f_climate + f_residual

where:
  f_climate    = beta_temp * temp_lag + beta_rain * rain_lag
  f_residual   ~ HSGP(Matern 3/2, alpha, rho)
```

### Key Prior Choices

| Parameter | Prior | Justification |
|-----------|-------|---------------|
| mu | Normal(0, 0.5) | Baseline Rt median=1, 95% CI [0.37, 2.7] |
| beta_climate | Normal(0, 0.5) | Weakly informative for temperature and rainfall |
| log_alpha | Normal(-1.2, 0.5) | GP marginal SD ~0.3 |
| log_rho | Normal(log(6), 0.5) | Short length scale ~6 weeks for faster dynamics |
| phi | Normal(0, 10)+ | Overdispersion |

### Observation Model Assumptions
- 100% case ascertainment (Singapore mandatory notification)
- Reporting delay (~5-7 days) absorbed at weekly resolution (Cori et al. 2013; Lau et al. 2022)
- Serial interval approximated by generation interval at weekly resolution
- Asymptomatic infections (~75%, Bhatt et al. 2013) not captured; Rt reflects symptomatic transmission

---

## Serotype Early Warning Analysis

The residual GP captures non-climate variation including serotype dynamics:
- **Cross-correlation (CCF)**: f_residual vs Shannon entropy and dominant serotype change rate
- **Pre-switch Rt dynamics**: P(f_residual rising) in 3 months before serotype switch
- **Entropy profiles**: Shannon entropy trajectory in [-6, +12] month windows around switches
- **Combined overlay**: f_residual (with CrI) and entropy per switch event

---

## Critical Caveats

### 1. Climate Effects Are Negligible
Temperature and rainfall have near-zero effects on weekly Rt, consistent with climate acting through slower mechanisms (vector ecology) not captured at weekly resolution.

### 2. Residual GP Dominates
The residual GP captures the vast majority of Rt variation, suggesting important unmeasured drivers (serotypes, immunity, spatial heterogeneity).

### 3. Single Time Series
No spatial variation; limited causal inference possible from national-level data alone.

---

## Environment

- **R**: renv for package management
- **Python**: uv for data acquisition
- **Stan**: CmdStan 2.38.0
- **Key packages**: cmdstanr, tidyverse, posterior, bayesplot, patchwork

---

## References

See `specification.md` for full reference list. Key papers:
- Cori et al. (2013) - Renewal equation framework
- Riutort-Mayol et al. (2023) - HSGP approximation
- Utarini et al. (2021, NEJM) - AWED Wolbachia trial
- Chan & Johansson (2012) - Dengue generation interval
- Lau et al. (2022, PLOS Comp Biol) - Dengue Rt estimation, serial interval approach
