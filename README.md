# Serotype Signals in a Bayesian Rt Model of Dengue in Singapore

Bayesian estimation of the effective reproduction number (Rt) for dengue in Singapore (2012--2022) using a renewal equation framework with Hilbert Space Gaussian Process (HSGP) approximation. The model decomposes log(Rt) into climate covariates and a residual GP, then uses the residual as a lens to study how serotype dynamics **associate with** transmission --- deliberately associational rather than mechanistic, since a companion model comparison shows explicit serotype-aware susceptibility structure does not absorb the residual.

**Course project for SPH6102: Inferential Infectious Disease Modelling**, National University of Singapore, AY2025/2026 Sem 2.

## Model

```
log(Rt) = mu + f_climate + f_residual

f_climate  = beta_temp * temp_lag4w + beta_rain * rain_lag4w
f_residual ~ HSGP(Matern 3/2, alpha, rho)   [110 basis functions]
```

Cases follow a renewal equation with negative binomial observation model. The generation interval PMF is derived from dengue biology (intrinsic incubation + extrinsic incubation + mosquito lifespan; Chan & Johansson 2012).

## Key Findings

1. **Climate explains <1% of Rt variance** at weekly resolution. Temperature and rainfall multiplicative effects on Rt: 0.989 [0.969, 1.009] and 1.003 [0.986, 1.020]. This is consistent with climate acting through slower vector-ecology mechanisms not captured at this timescale.
2. **The residual GP captures ~99.5% of Rt variation**, with a short length scale (~2.9 weeks) that tracks individual outbreak dynamics.
3. **5 of 7 serotype switch events** show elevated residual Rt (P > 0.75 in 0--3 months post-switch), consistent with immunity-depletion driving transmission surges.
4. **Serotype entropy leads Rt by ~10 months** (CCF analysis), suggesting a potentially actionable early warning window for outbreak preparedness.
5. **Kernel sensitivity analysis** (Matern 1/2, 3/2, 5/2, squared exponential) confirms all findings are robust to kernel choice.
6. **Forecasting evaluation** shows the retrospective GP does not translate to useful prospective predictions (rapid regression to prior mean within ~3 weeks), honestly characterising the model's limitations.
7. **Susceptibility audit (M0/M1/M2)** adds explicit serotype-blind and serotype-aware susceptibility covariates. GP amplitude is unchanged across the full sensitivity sweep and $\beta_S$ credible intervals include zero in every configuration, showing that the residual cannot be captured by simple SIR-style susceptible accounting. Serotype--Rt associations remain real; the mechanistic "serotype-driven immunity depletion" reading is weakened.

## Project Structure

```
dengue_gaussian/
|-- code/                   # Analysis scripts (numbered in execution order)
|   |-- 01_acquire_data.py            # Download dengue, weather, serotype data
|   |-- 02_acquire_opendengue.py      # Acquire OpenDengue multi-country data
|   |-- 03_prepare_model_data.R       # Merge, lag, standardise, format for Stan
|   |-- 04_prior_predictive.R         # Prior predictive checks
|   |-- 05_model_climate_only.stan    # Renewal equation + HSGP Stan model
|   |-- 06_fit_model.R                # Compile and fit Stan model via CmdStan
|   |-- 07_posterior_predictive.R     # Posterior predictive checks
|   |-- 08_postprocess.R              # Posterior summaries, decomposition, effects
|   |-- 09_residual_decomposition.R   # Residual GP contribution accounting
|   |-- 10_serotype_analysis.R        # CCF, pre-switch dynamics, entropy profiles
|   |-- 11_serotype_selection.R       # Serotype selection advantage analysis
|   |-- 12_lag_sensitivity.R          # Climate-lag sensitivity analysis
|   |-- 13_susceptible_sensitivity.R  # Susceptible pool sensitivity analysis
|   +-- archive/                      # Exploratory/superseded scripts
|
|-- data/                   # Input data (acquired programmatically)
|   |-- raw_dengue_cases.csv        # Weekly notified cases (2012-2022), data.gov.sg
|   |-- raw_weather.csv             # Daily temperature/rainfall, Meteostat
|   |-- monthly_sero_type_props_all_data.csv  # Serotype proportions, NEA
|   |-- serotype_props_2013_2016.csv          # Early serotype data
|   |-- WeeklyInfectiousDiseaseBulletinCases.csv  # MOH bulletin data
|   |-- raw_cluster_counts.csv      # Dengue cluster counts
|   |-- dengue_clusters_archive.zip # Archived cluster data
|   +-- model_data.rds              # Prepared Stan input
|
|-- results/                # Model outputs (large .rds files gitignored)
|   |-- fit_model.rds               # Main model posterior draws (~104 MB)
|   |-- fit_model_tight_gp.rds      # Tight GP sensitivity fit (~104 MB)
|   |-- fit_kernel_*.rds            # Kernel comparison fits (4 x ~103 MB)
|   |-- figures/                    # All diagnostic and analysis figures
|   |-- forecasts/                  # Forecast evaluation RDS files
|   |-- parameter_summary.csv       # Parameter estimates with CrI
|   |-- decomposition.csv           # Climate vs residual variance decomposition
|   |-- rt_trajectories.csv         # Weekly Rt posterior summaries
|   |-- serotype_ccf_results.csv    # Cross-correlation: f_residual vs entropy
|   |-- serotype_post_switch_dynamics.csv  # P(elevated Rt) per switch event
|   +-- serotype_switch_timing.csv  # Switch dates, P(rising), lead months
|
|-- report/                 # Written report and slides
|   |-- report.qmd                  # Quarto source
|   |-- slides.qmd                  # Slides source
|   |-- slides.pdf                  # Compiled slides
|   +-- references.bib              # Bibliography
|
|-- renv.lock               # R package lockfile
|-- pyproject.toml           # Python project config (uv)
|-- uv.lock                 # Python dependency lockfile
+-- .gitignore
```

## Reproduction

### Prerequisites

- **R** (>= 4.3) with `renv`
- **Python** (>= 3.11) with `uv`
- **CmdStan** (>= 2.38.0)

### Setup

```bash
# R dependencies
Rscript -e 'renv::restore()'

# Python dependencies
uv sync

# Install CmdStan (if not already installed)
Rscript -e 'cmdstanr::install_cmdstan(version = "2.38.0")'
```

### Run the pipeline

Scripts are numbered in execution order:

```bash
# 1. Acquire data
uv run code/01_acquire_data.py

# 2. Prepare model data
Rscript code/03_prepare_model_data.R

# 3. Prior predictive checks (optional)
Rscript code/04_prior_predictive.R

# 4. Fit model (~30-60 min depending on hardware)
Rscript code/06_fit_model.R

# 5. Posterior predictive checks
Rscript code/07_posterior_predictive.R

# 6. Postprocess
Rscript code/08_postprocess.R

# 7. Residual decomposition
Rscript code/09_residual_decomposition.R

# 8. Serotype analysis
Rscript code/10_serotype_analysis.R

# 9. Render report
quarto render report/report.qmd
```

## References

- Cori A, Ferguson NM, Fraser C, Cauchemez S (2013). A new framework and software to estimate time-varying reproduction numbers during epidemics. *Am J Epidemiol*, 178(9), 1505--1512.
- Riutort-Mayol G, Burkner PC, Andersen MR et al. (2023). Practical Hilbert space approximate Bayesian Gaussian processes for probabilistic programming. *Stat Comput*, 33, 17.
- Chan M, Johansson MA (2012). The incubation periods of dengue viruses. *PLOS ONE*, 7(11), e50972.
- Finch E, Chang C, Kucharski A et al. (2025). Climate variation and serotype competition drive dengue outbreak dynamics in Singapore. *Nat Commun*, 16, 11364.
- Lau M et al. (2022). Comparing dengue transmission dynamics with and without considering serotype interactions. *PLOS Comp Biol*.
