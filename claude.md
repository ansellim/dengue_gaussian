# Project Memory: Dengue Rt Estimation with Gaussian Processes

## Project Overview
Bayesian epidemiological model for estimating the effective reproduction number (Rt) of dengue in Singapore (2012-2022) using:
- Renewal equation framework
- Hilbert Space Gaussian Process (HSGP) approximation with Matern 3/2 kernel
- Semiparametric decomposition: log(Rt) = mu + f_climate + f_wolbachia + f_npi + f_residual
- Negative binomial observation model

## Current Status (Updated March 2026)

### What Has Been Done
1. **Data Pipeline**: Complete acquisition and preparation of dengue cases, weather, NPI (OxCGRT), Wolbachia coverage, and serotype data
2. **Model Development**: Multiple model versions tested (v1-v4), with short GP model (rho~6 weeks) as current best
3. **Model Validation**: Prior predictive checks, posterior predictive checks, residual diagnostics
4. **Sensitivity Analysis**: Pre-2020 analysis without NPI to assess robustness to confounding
5. **Post-hoc Analysis**: Serotype dynamics vs residual GP correlation
6. **Dose-Response Analysis**: Wolbachia coverage vs Rt reduction with extrapolation to 100% coverage

### Key Findings
| Parameter | Estimate | 95% CI | Interpretation |
|-----------|----------|--------|----------------|
| Wolbachia effect | 0.64 | [0.30, 1.35] | ~36% Rt reduction, 87% probability of protective effect |
| NPI effect | 1.21 | [0.86, 1.77] | CONFOUNDED with DENV-3 emergence in 2020 |
| Temperature | 0.99 | [0.96, 1.02] | Negligible effect at weekly resolution |
| Rainfall | 1.01 | [0.99, 1.04] | Negligible effect at weekly resolution |

### Known Limitations
1. **NPI-Serotype Confounding**: Cannot separate NPI effects from DENV-3 emergence (both in 2020)
2. **Wolbachia Coverage Uncertainty**: Only ~30% coverage by 2022, effect diluted at national level
3. **Residual GP Dominates**: ~90% of variance in residual GP, suggesting unmeasured drivers
4. **Single Time Series**: No spatial variation, limited causal inference

---

## File Structure

### Code Files

| File | Purpose | Status |
|------|---------|--------|
| `01_acquire_data.py` | Download dengue, weather, NPI, Wolbachia, serotype data | Complete |
| `02_prepare_model_data.R` | Merge, lag, standardize, format for Stan | Complete |
| `03_prior_predictive.R` | Verify priors generate plausible data | Complete |
| `03_model0_baseline.stan` | Model 0: Single GP baseline | Complete |
| `04_model1_climate.stan` | Model 1: Climate covariates + residual GP | Complete |
| `05_model2_full.stan` | Model 2: Full model (original, rho~15 weeks) | Complete |
| `05_model2_full_short_gp.stan` | Model 2: Short GP (rho~6 weeks) - **PREFERRED** | Complete |
| `05_model2_no_npi.stan` | Model 2: No NPI covariate (for pre-2020 analysis) | Complete |
| `06_fit_models.R` | Fit Models 0, 1, 2 with comparison | Complete |
| `07_postprocess.R` | Posterior summaries, decomposition, effect sizes | Complete |
| `08_serotype_analysis.R` | Post-hoc: residual GP vs serotype dynamics | Complete |
| `09_posterior_predictive.R` | Model validation via posterior predictive checks | Complete |
| `10_fit_short_gp.R` | Fit and compare short GP model | Complete |
| `11_sensitivity_pre2020.R` | Pre-2020 sensitivity analysis (no NPI model) | Complete |
| `12_dose_response.R` | Wolbachia dose-response extrapolation | Complete |

### Data Files

| File | Description |
|------|-------------|
| `data/model_data.rds` | Prepared Stan data with all covariates |
| `data/raw_dengue_cases.csv` | Weekly case counts (2012-2022) |
| `data/raw_weather.csv` | Temperature and rainfall from Meteostat |
| `data/raw_npi_oxcgrt.csv` | OxCGRT Stringency Index for Singapore |
| `data/raw_wolbachia.csv` | Wolbachia coverage estimates |

### Results Files

| File | Description |
|------|-------------|
| `results/fit_model2_short_gp.rds` | **Current best model fit** |
| `results/fit_model2_pre2020.rds` | Pre-2020 sensitivity analysis fit |
| `results/fit_model2_v3_oxcgrt.rds` | Original long GP model fit |

---

## Analysis Pipeline

```
┌─────────────────────────────────────────────────────────────────────┐
│ 01_acquire_data.py                                                   │
│   └── Downloads: dengue, weather, NPI, Wolbachia, serotype          │
├─────────────────────────────────────────────────────────────────────┤
│ 02_prepare_model_data.R                                              │
│   └── Creates: data/model_data.rds                                  │
├─────────────────────────────────────────────────────────────────────┤
│ 03_prior_predictive.R                                                │
│   └── Validates: priors generate plausible Rt and cases             │
├─────────────────────────────────────────────────────────────────────┤
│ 06_fit_models.R or 10_fit_short_gp.R                                │
│   └── Creates: results/fit_model2_short_gp.rds (PREFERRED)          │
├─────────────────────────────────────────────────────────────────────┤
│ 07_postprocess.R                                                     │
│   └── Creates: effect estimates, decomposition plots                │
├─────────────────────────────────────────────────────────────────────┤
│ 09_posterior_predictive.R                                            │
│   └── Validates: model fits data, checks residual autocorrelation   │
├─────────────────────────────────────────────────────────────────────┤
│ 08_serotype_analysis.R                                               │
│   └── Post-hoc: correlates residual GP with serotype dynamics       │
├─────────────────────────────────────────────────────────────────────┤
│ 11_sensitivity_pre2020.R                                             │
│   └── Sensitivity: pre-2020 model without NPI covariate             │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Model Specification

### Current Best Model: `05_model2_full_short_gp.stan`

```
log(Rt) = mu + f_climate + f_wolbachia + f_npi + f_residual

where:
  f_climate    = beta_temp * temp_lag + beta_rain * rain_lag
  f_wolbachia  = beta_wolbachia * coverage_t
  f_npi        = beta_npi * stringency_t
  f_residual   ~ HSGP(Matern 3/2, alpha, rho)
```

### Key Prior Choices

| Parameter | Prior | Justification |
|-----------|-------|---------------|
| mu | Normal(0, 0.5) | Baseline Rt median=1, 95% CI [0.37, 2.7] |
| beta_wolbachia | Normal(-0.3, 0.5) | AWED trial: 77% efficacy, attenuated by partial coverage |
| beta_npi | Normal(0, 0.5) | Agnostic: peridomestic Aedes may increase exposure at home |
| log_alpha | Normal(-1.2, 0.5) | GP marginal SD ~0.3 |
| log_rho | Normal(log(6), 0.5) | **Short length scale** ~6 weeks for faster dynamics |
| phi | Normal(0, 10)+ | Overdispersion |

### Why Short GP (rho~6 weeks)?
- Original model (rho~15 weeks) had high residual autocorrelation (0.84 at lag 1)
- Short GP reduced autocorrelation to 0.51 and improved 90% coverage to 94%
- Trade-off: more flexible GP may absorb some covariate signal

---

## Key Results Summary

### Effect Estimates (Short GP Model)

```
INTERVENTION EFFECTS (multiplicative on Rt):
  Wolbachia:  0.64 [0.30, 1.35]  P(reduces Rt) = 87%
  NPI:        1.21 [0.86, 1.77]  CONFOUNDED with DENV-3

CLIMATE EFFECTS (per 1 SD change):
  Temperature: 0.99 [0.96, 1.02]  Negligible
  Rainfall:    1.01 [0.99, 1.04]  Negligible

VARIANCE DECOMPOSITION:
  Climate:     0.9%
  Wolbachia:   4.0%
  NPI:         4.2%
  Residual GP: 89.8%  ← Most variation unexplained by covariates
```

### Sensitivity Analysis: Pre-2020 (No NPI)

```
Wolbachia effect is ROBUST to exclusion of 2020+ period:
  Pre-2020 (no NPI): 0.65 [0.26, 1.57]  P(reduces Rt) = 84%
  Full model:        0.64 [0.30, 1.35]  P(reduces Rt) = 87%
  Difference: +0.01 (negligible)
```

### Serotype Correlation

Strong correlation between residual GP and DENV-3 emergence:
- Correlation with DENV-3 proportion: r = 0.81
- DENV-3 emerged in 2020 simultaneously with COVID NPIs
- This explains why NPI effect cannot be causally interpreted

### Dose-Response Analysis (Extrapolation)

Projected Wolbachia effects at different coverage levels (assuming linear dose-response):

| Coverage | Rt Multiplier (95% CI) | Rt Reduction (95% CI) | Status |
|----------|------------------------|----------------------|--------|
| 10% | 0.96 (0.89–1.03) | 4% (-3%–11%) | Observed |
| 20% | 0.91 (0.78–1.06) | 9% (-6%–22%) | Observed |
| 30% | 0.87 (0.69–1.10) | 13% (-10%–31%) | Observed |
| 50% | 0.80 (0.55–1.16) | 20% (-16%–45%) | **Extrapolated** |
| 70% | 0.73 (0.43–1.24) | 27% (-24%–57%) | **Extrapolated** |
| 100% | 0.64 (0.30–1.35) | 36% (-35%–70%) | **Extrapolated** |

**AWED Trial comparison:** 77% efficacy (Rt multiplier ~0.23) falls outside our 95% CI at 100% coverage, suggesting either different contexts or high uncertainty in our estimate.

**Critical caveat:** Extrapolation beyond 30% assumes linear dose-response, which cannot be verified. Alternative shapes (threshold, saturation, accelerating) are plausible.

---

## Critical Caveats

### 1. NPI Effect is NOT Causal
The NPI coefficient (1.21) reflects "something changed in 2020-2022" but cannot distinguish:
- COVID lockdowns (more time at home with peridomestic Aedes)
- DENV-3 emergence (susceptible population)
- Disrupted vector control
- Climate anomalies
- Healthcare-seeking changes

### 2. Wolbachia Effect is Plausible but Wide
- 95% CI includes 1 (no effect)
- ~30% coverage limits population-level signal
- GP may absorb some of the temporal trend in coverage

### 3. Residual GP Dominates
~90% of variance in residual GP suggests:
- Important unmeasured drivers (serotypes, immunity, spatial heterogeneity)
- Model structure may not fully capture transmission dynamics

---

## Gaps and Future Work

### Incomplete Analyses
1. **Leave-future-out cross-validation**: Specified but not implemented
2. **Generation interval sensitivity**: Short/baseline/long GI not tested
3. **Climate lag sensitivity**: Only 4-week lag tested (2, 6 weeks pending)

### Potential Extensions (see below for data sources)
1. Vector abundance data (Gravitrap counts)
2. Spatial heterogeneity (district-level analysis)
3. Human mobility (Google/Apple mobility data)
4. Population immunity proxies
5. Additional vector control interventions

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
- Lim et al. (2024, Lancet Microbe) - Singapore Wolbachia efficacy
