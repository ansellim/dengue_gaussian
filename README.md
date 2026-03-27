# Semiparametric Decomposition of Dengue Rt via Composite Gaussian Process Kernels

Bayesian estimation of the effective reproduction number (Rt) for dengue in Singapore (2012-2022) using a renewal equation framework with Hilbert Space Gaussian Process (HSGP) approximation. The model decomposes log(Rt) into additive components for climate, Wolbachia vector control, COVID-19 NPIs, and a residual GP capturing unmeasured drivers.

See `specification.md` for the full model specification, prior justification, and analysis plan.

## Pipeline

All code is in the `code/` directory. Scripts are numbered in execution order. Run each from the `code/` directory.

| Step | File | Language | Description |
|------|------|----------|-------------|
| 1 | `01_acquire_data.py` | Python | Downloads weekly dengue case counts (data.gov.sg), meteorological data (Meteostat), Wolbachia coverage, and OxCGRT NPI stringency index. Outputs `data/raw_*.csv`. |
| 2 | `02_prepare_model_data.R` | R | Merges datasets, creates lagged climate covariates, standardizes, computes the generation interval PMF from component distributions, and formats everything for Stan. Outputs `data/model_data.rds`. |
| 3 | `03_model0_baseline.stan` | Stan | Model 0: single GP on log(Rt), no covariates. Benchmark model. |
| 4 | `04_model1_climate.stan` | Stan | Model 1: climate covariates (lagged temperature and rainfall) + residual GP. |
| 5 | `05_model2_full.stan` | Stan | Model 2: climate + Wolbachia coverage + NPI intensity + residual GP. Primary model. |
| 6 | `06_fit_models.R` | R | Compiles and fits all three Stan models via CmdStan. Runs MCMC diagnostics, LOO-CV comparison. Outputs `results/fit_model*.rds`. |
| 7 | `07_postprocess.R` | R | Extracts posterior summaries, generates Rt trajectory plots, variance decomposition, posterior predictive checks, and effect size estimates. Outputs figures and CSVs to `results/`. |
| 8 | `08_serotype_analysis.R` | R | Post-hoc analysis comparing the Model 2 residual GP against serotype dynamics (2013-2022). Tests whether unexplained Rt variation aligns with serotype replacement events. |

## Requirements

- **Python**: managed via `uv` (see `pyproject.toml`). Key packages: `requests`, `meteostat`, `pandas`.
- **R**: managed via `renv` (see `renv/`). Key packages: `cmdstanr`, `tidyverse`, `posterior`, `bayesplot`, `patchwork`.
- **Stan**: CmdStan (tested with v2.38.0).

## Data

All data are stored in `data/` and are acquired programmatically by `01_acquire_data.py` from open sources. Serotype proportion data are sourced from a separate NEA dataset (see `specification.md` section 6.5).

## Results

Model outputs, diagnostic summaries, and figures are written to `results/`.
