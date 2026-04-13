# Semiparametric Decomposition of Dengue Rt with a Gaussian Process Residual

Bayesian estimation of the effective reproduction number (Rt) for dengue in Singapore (2012-2022) using a renewal equation framework with Hilbert Space Gaussian Process (HSGP) approximation. The model decomposes log(Rt) into additive components for climate covariates (lagged temperature and rainfall) and a residual GP capturing unmeasured drivers.

See `specification.md` for the full model specification, prior justification, and analysis plan.

## Pipeline

All code is in the `code/` directory. Scripts are numbered in execution order. Run each from the `code/` directory.

| Step | File | Language | Description |
|------|------|----------|-------------|
| 1 | `01_acquire_data.py` | Python | Downloads weekly dengue case counts (data.gov.sg) and meteorological data (Meteostat). Outputs `data/raw_*.csv`. |
| 2 | `02_prepare_model_data.R` | R | Merges datasets, creates lagged climate covariates, standardizes, computes the generation interval PMF from component distributions, and formats everything for Stan. Outputs `data/model_data.rds`. |
| 3 | `03_prior_predictive.R` | R | Verifies that the priors generate plausible Rt trajectories and case counts. |
| 4 | `05_model3_climate_only.stan` | Stan | Renewal equation model: log(Rt) = mu + beta_temp*temp + beta_rain*rain + f_residual(t), with f_residual an HSGP (Matern 3/2). |
| 5 | `13_fit_model3.R` | R | Compiles and fits the Stan model via CmdStan. Runs MCMC diagnostics. Outputs `results/fit_model3.rds`. |
| 6 | `07_postprocess.R` | R | Extracts posterior summaries, generates Rt trajectory plots, variance decomposition, and effect size estimates. Outputs figures and CSVs to `results/`. |
| 7 | `09_posterior_predictive.R` | R | Posterior predictive checks and residual autocorrelation diagnostics. |
| 8 | `08_serotype_analysis.R` | R | Post-hoc analysis comparing the residual GP against serotype dynamics (2013-2022). Tests whether unexplained Rt variation aligns with serotype replacement events. |

## Requirements

- **Python**: managed via `uv` (see `pyproject.toml`). Key packages: `requests`, `meteostat`, `pandas`.
- **R**: managed via `renv` (see `renv/`). Key packages: `cmdstanr`, `tidyverse`, `posterior`, `bayesplot`, `patchwork`.
- **Stan**: CmdStan (tested with v2.38.0).

## Data

All data are stored in `data/` and are acquired programmatically by `01_acquire_data.py` from open sources. Serotype proportion data are sourced from a separate NEA dataset (see `specification.md` section 6.5).

## Results

Model outputs, diagnostic summaries, and figures are written to `results/`.
