# Project Memory: Dengue Rt Estimation with Gaussian Processes

## Project Overview
Bayesian epidemiological model for estimating the effective reproduction number (Rt) of dengue in Singapore (2012-2022) using:
- Renewal equation framework
- Hilbert Space Gaussian Process (HSGP) approximation with Matern 3/2 kernel
- Semiparametric decomposition: log(Rt) = mu + f_climate + f_residual
- Negative binomial observation model

## Project Focus

The project focus is on **serotype-driven Rt dynamics and outbreak prediction**. The two primary aims are:
1. **Serotype effects on Rt**: Characterizing how serotype replacement events affect Rt via the immunity-depletion mechanism, using the residual GP (f_residual) as a proxy for non-climate drivers
2. **Outbreak prediction**: Evaluating serotype diversity (Shannon entropy) as an early warning indicator for Rt elevation, leveraging the ~10-month leading relationship identified in CCF analysis

## Current Status

### What Has Been Done
1. **Data Pipeline**: Complete acquisition and preparation of dengue cases, weather, and serotype data
2. **Model Development**: Climate-only renewal equation model with HSGP residual (Model 3: temp + rain + residual GP)
3. **Model Validation**: Prior predictive checks, posterior predictive checks (95% coverage: 97.7%, 80% coverage: 92.3%)
4. **Sensitivity Analysis**: Tight GP amplitude prior confirms variance decomposition is robust
5. **Serotype Analysis**: Cross-correlation, pre-switch Rt dynamics, entropy profiling complete

### Key Findings (Model 3)
| Parameter | Estimate | 95% CrI | Interpretation |
|-----------|----------|---------|----------------|
| Temperature | 0.989 | [0.969, 1.009] | Negligible effect at weekly resolution |
| Rainfall | 1.003 | [0.986, 1.020] | Negligible effect at weekly resolution |
| prop_climate | 0.4% (median) | | Climate explains <1% of Rt variance |
| prop_residual | 99.6% (median) | | Residual GP captures nearly all Rt variation |
| GP rho | 2.86 weeks (median) | [2.09, 3.71] | Short length scale tracks outbreak dynamics |
| GP alpha | 0.30 (median) | [0.27, 0.34] | Moderate amplitude on log(Rt) scale |

### Serotype Analysis Key Results
- 5 of 7 switch events show P(elevated f_residual) > 0.75 in 0-3 months post-switch
- DENV-2 to DENV-3 transitions (2020-2021) show strongest signals (P=1.0, median increase 0.2-0.3 on log scale)
- CCF peaks at lag +10 to +12 months (entropy leading f_residual), suggesting actionable lead time

---

## File Structure

### Code Files

| File | Purpose | Status |
|------|---------|--------|
| `01_acquire_data.py` | Download dengue, weather, serotype data | Complete |
| `02_prepare_model_data.R` | Merge, lag, standardize, format for Stan (climate covariates only) | Complete |
| `03_prior_predictive.R` | Verify priors generate plausible data | Complete |
| `05_model3_climate_only.stan` | Climate-only (temp + rain) + residual GP -- **CURRENT** | Complete |
| `07_postprocess.R` | Posterior summaries, decomposition, effect sizes | Complete |
| `08_serotype_analysis.R` | Serotype early warning: CCF, pre-switch dynamics, entropy | Complete |
| `09_posterior_predictive.R` | Model validation via posterior predictive checks | Complete |
| `13_fit_model3.R` | Fit climate-only model | Complete |

### Data Files

| File | Description |
|------|-------------|
| `data/model_data.rds` | Prepared Stan data with climate covariates |
| `data/raw_dengue_cases.csv` | Weekly case counts (2012-2022) |
| `data/raw_weather.csv` | Temperature and rainfall from Meteostat |
| `data/monthly_sero_type_props_all_data.csv` | Monthly serotype proportions |

### Results Files

| File | Description |
|------|-------------|
| `results/fit_model3.rds` | **Current best model fit** |
| `results/fit_model3_tight_gp.rds` | Tight GP sensitivity analysis fit |
| `results/serotype_ccf_results.csv` | CCF: f_residual vs entropy and turnover rate |
| `results/serotype_post_switch_dynamics.csv` | P(elevated) and median change per switch event |
| `results/serotype_switch_timing.csv` | Switch dates, P(rising), lead months |
| `results/serotype_residual_monthly.csv` | Monthly aggregated f_residual posteriors |

---

## Analysis Pipeline

```
01_acquire_data.py
  └── Downloads: dengue, weather, serotype data
02_prepare_model_data.R
  └── Creates: data/model_data.rds (climate covariates only)
03_prior_predictive.R
  └── Validates: priors generate plausible Rt and cases
13_fit_model3.R
  └── Creates: results/fit_model3.rds
07_postprocess.R
  └── Creates: effect estimates, decomposition plots
09_posterior_predictive.R
  └── Validates: model fits data, checks residual autocorrelation
08_serotype_analysis.R
  └── Serotype early warning: CCF, entropy, pre-switch dynamics
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
- **Cross-correlation (CCF)**: f_residual vs Shannon entropy and dominant serotype change rate; entropy leads by ~10 months
- **Pre-switch Rt dynamics**: P(f_residual rising) in 3 months before serotype switch
- **Entropy profiles**: Shannon entropy trajectory in [-6, +12] month windows around switches
- **Combined overlay**: f_residual (with CrI) and entropy per switch event
- **Post-switch elevation**: P(elevated f_residual) in 0-3 and 3-6 month windows after switch

---

## Critical Caveats

### 1. Climate Effects Are Negligible
Temperature and rainfall explain <1% of Rt variance at weekly resolution, consistent with climate acting through slower mechanisms (vector ecology) not captured at this timescale.

### 2. Residual GP Captures Nearly All Variation
The residual GP accounts for ~99.5% of Rt variance, indicating important unmeasured drivers (serotypes, immunity, spatial heterogeneity).

### 3. Single Time Series
No spatial variation; limited causal inference possible from national-level data alone.

### 4. Serotype Analysis Is Descriptive
With only ~3-4 major serotype switches in 10 years, formal hypothesis testing is not feasible. Evidence is based on whether temporal profiles of f_residual match predicted immunity-depletion patterns across individual events.

---

## Environment

- **R**: renv for package management
- **Python**: uv for data acquisition
- **Stan**: CmdStan 2.38.0
- **Key packages**: cmdstanr, tidyverse, posterior, bayesplot, patchwork, mgcv

---

## References

See `specification.md` for full reference list. Key papers:
- Cori et al. (2013) - Renewal equation framework
- Riutort-Mayol et al. (2023) - HSGP approximation
- Chan & Johansson (2012) - Dengue generation interval
- Finch et al. (2025, Nat Commun) - GAM smoothing of serotype proportions
- Lau et al. (2022, PLOS Comp Biol) - Dengue Rt estimation, serial interval approach
