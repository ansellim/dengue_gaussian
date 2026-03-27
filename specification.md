# Semiparametric Decomposition of Dengue Rt via Composite Gaussian Process Kernels

## Project Specification

### 1. Motivation and Scope

Singapore's dengue epidemiology is shaped by multiple interacting drivers operating at
different timescales: climate variation (temperature and rainfall affecting vector biology),
vector control interventions (particularly the phased rollout of Project Wolbachia),
non-pharmaceutical interventions during the COVID-19 pandemic, and slower-moving
processes such as serotype cycling and population immunity shifts. Standard approaches
to estimating the effective reproduction number (Rt) either use a single flexible smoother
(as in EpiNow2) that conflates all sources of variation, or purely parametric covariate
models that impose rigid functional forms.

This project builds a renewal equation-based model of dengue Rt in Singapore where
log(Rt) is decomposed into additive components, each reflecting a distinct epidemiological
driver with its own prior structure. The framework allows us to ask: conditional on
measured climate variation and known interventions, what residual temporal structure
remains in Rt, and what does that residual tell us about unmeasured drivers?

**What this project is not:** The Rt estimated here is an effective human-to-human
reproduction number derived from collapsing the human-mosquito-human transmission
cycle into a single generation interval kernel. It captures changes in both human-to-mosquito
and mosquito-to-human transmission rates, as well as vector population dynamics, without
distinguishing them. The decomposition is therefore a decomposition of *effective
transmissibility variation*, not a mechanistic attribution to specific transmission
pathway components. This is standard for renewal equation-based Rt estimation applied
to vector-borne diseases, but should be stated explicitly in any write-up.

**Scope:** A small project (1-2 weeks), implemented in Stan (with R or Python wrapper),
using only open-source data. The project progresses through a series of models of
increasing complexity.


### 2. Model Progression

**Model 0 (Baseline):** Single GP on log(Rt), renewal equation, negative binomial
observation model. No covariates. Purpose: sanity check, comparison benchmark.

**Model 1 (Climate):** Replace the single GP with:
  log(Rt) = intercept + beta_temp * temp_lag + beta_rain * rain_lag + f_residual(t)
where temp_lag and rain_lag are lagged weekly meteorological covariates, and f_residual
is a GP on time capturing unexplained variation. Purpose: quantify how much of Rt
variation is attributable to climate, and examine the residual structure.

**Model 2 (Climate + Interventions):** Extend Model 1 with:
  log(Rt) = intercept + beta_temp * temp_lag + beta_rain * rain_lag
             + beta_wolbachia * coverage_t + beta_npi * lockdown_t + f_residual(t)
where coverage_t is the approximate proportion of households under Wolbachia releases
and lockdown_t encodes COVID-19 NPI phases. Purpose: assess whether intervention
effects on Rt are detectable after accounting for climate.

**Post-hoc Serotype Analysis (after Model 2):** Rather than incorporating serotype data
as a covariate in the model, we perform a post-hoc analysis of the Model 2 residual GP
against serotype dynamics. The residual GP from Model 2 captures Rt variation not
explained by climate, Wolbachia, or NPIs. We ask: does this residual correlate with
serotype switching and diversity patterns?

Models 0-2 are fit on the full 2012-2022 period. The post-hoc serotype analysis is
restricted to 2013-2022 (the period for which serotype proportion data are available).
The residual GP posterior is averaged from weekly to monthly resolution for comparison
with the monthly serotype metrics.

Two serotype metrics are examined:
1. **Dominant serotype identity** — which serotype (D1-D4) has the highest proportion
   in each month, plotted as a categorical timeline alongside the residual GP
2. **Shannon entropy** — H = -Σ p_i log(p_i), measuring serotype diversity; higher
   entropy indicates more co-circulating serotypes and potentially larger aggregate
   susceptible pools

The analysis is descriptive: if the residual GP shows elevated transmissibility during
periods of serotype replacement (e.g., the D1→D2 switch in 2015-2016), this constitutes
evidence that serotype-driven immunity dynamics are an important unmeasured driver of
Rt variation in Singapore.


### 3. Generative Model Structure

The model has three layers.

#### 3.1 Layer 1: Latent Process for Rt

For Model 0:

    log(Rt) = mu + f(t)

where mu is an intercept and f(t) ~ GP(0, k(t, t')) is a zero-mean Gaussian process.

For Models 1-2, the GP is replaced by:

    log(Rt) = mu + X_t * beta + f_residual(t)

where X_t is a row vector of covariates at time t, beta is a vector of regression
coefficients, and f_residual(t) ~ GP(0, k(t, t')) captures residual temporal variation.

**Why not a separate periodic/seasonal kernel?** In equatorial Singapore, dengue
seasonality is largely climate-driven: monsoon patterns affect temperature and rainfall,
which in turn modulate Aedes breeding cycles and viral extrinsic incubation periods. A
separate periodic kernel would be collinear with climate covariates that already carry
the seasonal signal. Additionally, fixing the period at 52 weeks assumes a rigid annual
cycle that Singapore's dengue epidemiology does not follow cleanly -- multi-year
epidemic waves driven by serotype switching, interannual ENSO variation, and secular
trends in vector control all disrupt simple annual periodicity. Estimating the period
would add a parameter that is poorly identified given the collinearity with climate.
The cleaner solution is to let the climate covariates absorb whatever seasonal structure
exists, and let the residual GP capture any remaining temporal correlation.

#### 3.2 Layer 2: Transmission Model (Renewal Equation)

Expected cases at week t are given by the discretized renewal equation:

    lambda_t = R_t * sum_{s=1}^{S} c_{t-s} * w_s

where:
- c_{t-s} is the observed case count at week t-s
- w_s is the discretized generation interval PMF (probability that a secondary case
  occurs s weeks after the primary case)
- S is the maximum generation interval in weeks
- The sum (total infectiousness) represents the effective infectious pressure at time t

For the first S weeks, the renewal equation cannot be computed because it requires
case counts from before the observation window. We handle this by treating the first S
weeks as a burn-in period: observed cases are used as given, and Rt is only modeled
from week S+1 onward.

**Generation interval construction.** Rather than fitting a single parametric
distribution to published serial interval estimates (which conflates the generation
interval with the serial interval — see Gostic et al. 2020), we construct the
generation interval from its biological components via Monte Carlo convolution:

    GI = T_iip + T_viremia + T_eip + T_bite

where:
- T_iip ~ Gamma(mean=5.9d, SD=1.5d): intrinsic incubation period, infection to
  symptom onset / viremia onset. Source: Chan & Johansson (2012, PLOS ONE),
  systematic review (n=153): "The mean IIP estimate was 5.9 days, with 95%
  expected between days 3 and 10." Gamma fit: shape=16, 95% range [3.4, 9.0]d.
- T_viremia ~ Uniform(0, 5d): time within the ~5-day viremic period until a
  mosquito acquires the virus. Sources: CDC dengue materials ("infectious to
  mosquitoes from 2 days before, to 5 days after illness onset"); literature
  consensus of 4-5 day high-viremia window. Note: viremia begins ~1d before
  symptom onset; this slightly overestimates T_iip + T_viremia by ~0.5d,
  negligible at weekly resolution.
- T_eip ~ Gamma(mean depends on temperature, SD ~CV 0.35): extrinsic incubation
  period in the mosquito, the dominant source of GI variation. Source:
  Chan & Johansson (2012, PLOS ONE) meta-analysis (n=146), log-normal model with
  temperature covariate (βT = −0.08/°C): "At 25°C: mean 15 days (95% CI: 10, 20)"
  and "At 30°C: mean 6.5 days (95% CI: 4.8, 8.8)". Interpolation at 27°C gives
  ~11d. Corroborated by Rohani et al. (2009): "virus first detected on Day 9 at
  26°C and 28°C and on Day 5 at 30°C" (first-detection times; means are longer).
- T_bite ~ Exponential(mean=2d): time from mosquito becoming infectious to
  biting a susceptible human. Derived from: Delatte et al. (2009, PLOS ONE):
  gonotrophic cycle "mean duration of 7.10 (±3.38) days at 27°C"; Scott et al.
  (2000, J Med Entomol): Ae. aegypti takes ~1.5 blood meals per gonotrophic
  cycle (42% double meals, 5% triple). Effective feeding interval ~4-5d;
  waiting time to next bite after EIP completion averages ~2d.

This approach: (1) targets the generation interval (infection-to-infection) rather
than the serial interval (symptom-to-symptom), which the renewal equation requires;
(2) draws each component from a well-characterized literature distribution; and
(3) maps the sensitivity analysis to biologically meaningful temperature variation
via the EIP.

At weekly resolution, we discretize into a PMF over weeks 1 through 6 (S = 6).
Sensitivity analysis varies the EIP to reflect Singapore's temperature range:
- Short GI (warm, ~30°C): EIP mean 6.5d → total GI mean ~17d (~2.4 weeks)
- Baseline GI (~27°C): EIP mean 11d → total GI mean ~21d (~3.0 weeks)
- Long GI (cool, ~25°C): EIP mean 15d → total GI mean ~25d (~3.6 weeks)

#### 3.3 Layer 3: Observation Model

    c_t ~ NegBin(lambda_t, phi)

Negative binomial parameterized with mean lambda_t and overdispersion phi, where
Var(c_t) = lambda_t + lambda_t^2 / phi. The NegBin accommodates week-to-week
variance in dengue case counts that exceeds what a Poisson model allows. Dengue
surveillance data is known to exhibit overdispersion due to clustered transmission,
spatially heterogeneous reporting, and serotype-dependent case detection ratios
(Tan et al., 2019, Am J Epidemiol).


### 4. Gaussian Process Specification

#### 4.1 Kernel Choice: Matern 3/2

The Matern 3/2 kernel:

    k(t, t') = alpha^2 * (1 + sqrt(3) * |t - t'| / rho) * exp(-sqrt(3) * |t - t'| / rho)

produces once-differentiable sample paths. This is smoother than a Matern 1/2 (Ornstein-
Uhlenbeck, which produces rough, noisy paths) but less smooth than the squared
exponential (infinitely differentiable, which tends to oversmooth epidemic dynamics).
The Matern 3/2 is the default in EpiNow2 (Abbott et al., 2020) for Rt estimation
precisely because it permits Rt to change relatively quickly while maintaining continuity.

**Why not a squared exponential?** The SE kernel enforces infinite smoothness, which
can oversmooth genuine abrupt changes in Rt (e.g., at the onset of COVID lockdowns
or following a serotype switch). The Matern 3/2 allows sharper transitions while still
providing regularization.

The kernel has two hyperparameters:
- alpha (marginal standard deviation): controls the amplitude of variation in log(Rt)
- rho (length scale): controls the temporal correlation -- how quickly log(Rt) can change

#### 4.2 Hilbert Space GP Approximation

Exact GP inference requires O(N^3) computation for N time points (~520 modeled weeks
after burn-in). We use the Hilbert space basis function approximation (Riutort-Mayol
et al., 2023, Statistics and Computing; Solin and Sarkka, 2020):

    f(t) ~ sum_{j=1}^{M} beta_j * phi_j(t)

where:
- phi_j(t) = (1/sqrt(L)) * sin(j * pi * (t + L) / (2L)) are Laplacian eigenfunctions
  on the domain [-L, L]
- L = boundary_factor * (max(t) - min(t)) / 2 is the boundary (typically
  boundary_factor = 1.5)
- beta_j ~ Normal(0, sqrt(S(lambda_j))) where S is the spectral density of the
  Matern 3/2 kernel evaluated at eigenvalue lambda_j = (j * pi / (2L))^2
- M is the number of basis functions (typically 20-40)

The spectral density of the Matern 3/2 kernel in 1D is:

    S(omega) = alpha^2 * 4 * (sqrt(3)/rho)^3 / ((sqrt(3)/rho)^2 + omega^2)^2

This reduces GP inference from O(N^3) to O(M * N), making it tractable in Stan via
HMC/NUTS.

#### 4.3 Number of Basis Functions

The approximation quality improves with M but each basis function adds a parameter.
Riutort-Mayol et al. (2023) recommend choosing M such that the highest eigenfrequency
exceeds the effective bandwidth of the GP's spectral density. Since time is scaled to
[-1, 1], rho in scaled units is small (e.g., 10 weeks ≈ 0.035 in scaled time for
N_model ≈ 574), requiring more basis functions than an unscaled parameterization.
Using M >= ceil(sqrt(3)/rho_min × 2L/π) at the 5th percentile of the rho prior:
- Models 0, 1 (rho median ~10 weeks, 5th percentile ~4.4 weeks): M ≈ 110
- Models 2+ (rho median ~15 weeks, 5th percentile ~6.6 weeks): M ≈ 75
We use M = 110 for all models for simplicity.


### 5. Prior Specification

All priors are chosen to be weakly informative, encoding broad epidemiological
plausibility without being strongly prescriptive. The rationale for each:

#### 5.1 Intercept: mu ~ Normal(0, 0.5)

This places the prior median of Rt at exp(0) = 1, consistent with endemic dengue where
Rt fluctuates around 1. The SD of 0.5 on the log scale means the prior 95% interval for
the baseline Rt is approximately [exp(-1), exp(1)] = [0.37, 2.72], which comfortably
spans the plausible range for dengue in Singapore.

#### 5.2 GP length scale: log(rho) ~ Normal(log(rho_weeks), 0.5)

The GP length scale rho is in weeks (time is passed to Stan centered but unscaled).
We parameterize on the log scale for better HMC geometry.

For Models 0 and 1 (no covariates; GP must capture all temporal structure):
  log(rho) ~ Normal(log(10), 0.5), giving a prior median of 10 weeks with 95%
  interval [3.8, 26.6] weeks. The shorter length scale allows the GP to track
  epidemic dynamics at the timescale of individual outbreaks.

For Models 2+ (climate and intervention covariates absorb short-term variation):
  log(rho) ~ Normal(log(15), 0.5), giving a prior median of 15 weeks with 95%
  interval [5.6, 40.0] weeks. The longer length scale reflects that the residual
  GP needs only to capture variation not explained by the covariates (e.g.,
  serotype dynamics, immunity shifts).

#### 5.3 GP amplitude: log(alpha) ~ Normal(−1.2, 0.5)

The GP marginal standard deviation alpha is parameterized on the log scale for better
HMC geometry. The prior places the median at exp(−1.2) ≈ 0.3, with 95% interval
approximately [0.11, 0.80]. This implies typical excursions in log(Rt) of roughly
±0.3 to ±0.8, corresponding to Rt varying by factors of roughly 1.3× to 2.2× around
the covariate-predicted value. Dengue Rt in Singapore historically ranges from about
0.5 during inter-epidemic troughs to 2–3 during epidemic peaks, so this prior is
permissive but prevents implausibly large GP amplitudes.

#### 5.4 Overdispersion: phi ~ half-Normal(0, 10)

For the negative binomial parameterization with Var = mu + mu^2/phi, larger phi means
less overdispersion (approaching Poisson as phi -> infinity). A half-Normal(0, 10) is
broadly permissive, allowing for substantial overdispersion (phi near 1) or near-Poisson
behavior (phi >> 10). This follows EpiNow2's approach. Weekly dengue counts in
Singapore range from ~50 to ~2000+, and the variance structure likely varies across
this range.

#### 5.5 Covariate coefficients (Models 1-2): beta ~ Normal(0, 0.5)

Weakly informative, centered at zero (no effect), with SD chosen so that a 1-SD change
in the covariate (after standardization) shifts log(Rt) by 0.5, corresponding to a
multiplicative effect on Rt of exp(0.5) ~ 1.65. This is permissive for climate effects
(where we expect modest effects) and for Wolbachia (where published efficacy estimates
suggest substantial Rt reduction of ~40-70% in release areas, though diluted at the
national level).


### 6. Data Specification

#### 6.1 Outcome: Weekly Dengue Case Counts (2012-2022)

Source: MOH Weekly Infectious Diseases Bulletin via data.gov.sg
(dataset ID: d_ca168b2cb763640d72c4600a68f9909e). Covers January 2012 to December
2022, weekly resolution, ~570 observations. Includes dengue fever and dengue
hemorrhagic fever (aggregated to total dengue).

#### 6.2 Meteorological Covariates: Temperature and Rainfall (2012-2022)

Source: Meteostat Python package, Singapore Changi Airport station (WMO: 48698).
Daily resolution, aggregated to weekly to match case data. Variables:
- Weekly mean temperature (degrees C)
- Weekly total rainfall (mm)

Both covariates are standardized (z-scored) before entering the model. A lag structure
is applied to reflect the delay between climate conditions and case onset: temperature
and rainfall at time t affect mosquito population dynamics and viral replication, which
in turn affect human cases ~2-6 weeks later. The baseline specification uses a 4-week
lag (matching the approximate sum of the extrinsic incubation period and the
human-to-notification delay). Sensitivity analysis should explore lags of 2, 4, and
6 weeks.

Data validation: cross-reference Meteostat data against the data.gov.sg Historical
Daily Weather Records (dataset ID: d_03bb2eb67ad645d0188342fa74ad7066, covering
2009-2017) for overlapping years.

#### 6.3 Wolbachia Coverage (2016-2022)

Source: Reconstructed from NEA press releases and published literature.

Approximate proportion of Singapore households under Wolbachia releases:

| Period          | Coverage (%) | Source                                         |
|-----------------|-------------|------------------------------------------------|
| 2012 - Sep 2016 | 0           | Pre-intervention                                |
| Oct 2016 - 2018 | ~0.5        | 3 small pilot sites; negligible at national level|
| 2019             | ~5.7        | Lim et al. (2024) Lancet Microbe               |
| 2020             | ~23.3       | Lim et al. (2024) Lancet Microbe               |
| 2021             | ~24-26      | MSE parliamentary replies                       |
| Mid-2022         | ~26-31      | NEA press release (Jul 2022 expansion)          |

Entered as a continuous covariate (proportion, 0 to 1), linearly interpolated between
known data points. Not standardized (already on a natural scale).

#### 6.4 COVID-19 NPI Phases (2020)

Source: Public government announcements. Encoded as a categorical step function:

| Period                | NPI Level | Coding |
|-----------------------|-----------|--------|
| Before 7 Apr 2020     | None      | 0      |
| 7 Apr - 1 Jun 2020    | Circuit Breaker (full lockdown) | 1 |
| 2 Jun - 18 Jun 2020   | Phase 1 (partial reopening) | 0.67 |
| 19 Jun - 28 Dec 2020  | Phase 2 (most activities resumed) | 0.33 |
| 28 Dec 2020 - mid 2021| Phase 3 (near-normal) | 0.1 |
| After mid-2021        | None      | 0      |

This coding reflects the approximate intensity of mobility restrictions. An alternative
is to use Google Community Mobility Reports as a continuous covariate, though these
are only available through October 2022.


#### 6.5 Serotype Proportions (2013-2022)

Source: NEA sequencing data, stored in the DengueFever project directory. Two files
are combined to cover the full period:
- `serotype_props_2013_2016.csv`: monthly serotype proportions (D1-D4) from NEA
  sequencing, January 2013 to December 2016
- `monthly_sero_type_props_all_data.csv`: monthly serotype proportions expanded to
  weekly rows, January 2017 to December 2022

Both files contain columns D1_prop, D2_prop, D3_prop, D4_prop summing to ~1. The
2013-2016 file has actual monthly variation; the 2017-2022 file repeats monthly values
across weeks within each month.

For the post-hoc analysis, the combined series is used at monthly resolution.

**Smoothing.** Raw monthly proportions are noisy because only a subset of cases
are serotyped each month, producing sampling variation that causes spurious
month-to-month flips in the dominant serotype identity. Following Finch et al.
(2025, Nat Comms), we smooth the raw proportions using a GAM multinomial logistic
regression (R package mgcv, family = multinom) where each serotype's log-odds
(relative to a reference) is a smooth function of time (thin plate regression
spline). The smoothed proportions represent the estimated underlying trajectory
of serotype composition, filtering out sampling noise.

Two derived metrics are computed from the smoothed proportions:
- Dominant serotype: argmax of smoothed proportions per month
- Serotype switch: month where the smoothed dominant serotype identity changes
  (no persistence filter needed — the smoothing eliminates spurious flips)

Shannon entropy (H = -Σ_{i=1}^{4} p_i log(p_i), with convention 0 log 0 = 0)
is computed on both the raw proportions (reflecting observed diversity) and the
smoothed proportions (for trend visualization).

This data is NOT used in the Stan models (Models 0-2). It is used only in the
post-hoc analysis of the Model 2 residual GP.


### 7. Implementation Plan

All code in Stan with R (or Python) wrapper scripts.

#### 7.1 File Structure

```
dengue_gaussian/
  code/
    01_acquire_data.py        # Download and process data
    02_prepare_model_data.R   # Merge, lag, standardize, format for Stan
    03_model0_baseline.stan   # Model 0: single GP baseline
    04_model1_climate.stan    # Model 1: climate covariates + residual GP
    05_model2_full.stan       # Model 2: climate + Wolbachia + NPI + residual GP
    06_fit_models.R           # Fit all models, diagnostics
    07_postprocess.R          # Posterior summaries, decomposition plots, comparison
    08_serotype_analysis.R    # Post-hoc: residual GP vs serotype dynamics (2013-2022)
  data/                       # Raw and processed data
  results/                    # MCMC output, figures
  specification.md            # This file
  README.md                   # Pipeline overview
```

#### 7.2 Diagnostics Checklist

- MCMC convergence: R-hat < 1.01, bulk ESS > 400, tail ESS > 400 for all parameters
- Divergent transitions: zero (or investigate source if present)
- Prior predictive check: simulate Rt trajectories from priors, verify epidemiological
  plausibility
- Posterior predictive check: overlay observed cases on posterior predictive distribution
- Leave-future-out cross-validation: hold out last 26 weeks, compare pointwise predictive
  accuracy across models
- Sensitivity analysis: generation interval specification (3 variants), climate lag
  (2, 4, 6 weeks), prior widths on GP hyperparameters

#### 7.3 Model Comparison

- Leave-future-out predictive log-likelihood (preferred over LOO-CV for time series)
- Posterior predictive coverage: what fraction of observed weeks fall within the
  80%/95% posterior predictive intervals?
- Visual comparison of Rt trajectories across models
- For Models 1-2: examine the residual GP -- does adding covariates reduce the amplitude
  and shorten the length scale of the residual, indicating that the covariates are
  explaining structured variation?


#### 7.4 Post-hoc Serotype Analysis

**What f_residual is.** The residual GP f_residual(t) is the component of log(Rt) not
explained by the intercept and covariates in Model 2:

    f_residual(t) = log(Rt) - mu - beta_temp * temp_lag - beta_rain * rain_lag
                    - beta_wolbachia * coverage_t - beta_npi * lockdown_t

It lives on the log(Rt) scale: f_residual(t) = 0 means Rt is fully explained by
covariates; f_residual(t) > 0 means Rt is higher than covariates predict (unexplained
excess transmission); f_residual(t) < 0 means lower. It has a full posterior
distribution at each time point, not a single value.

**The hypothesis.** Serotype replacement events increase Rt because infection with one
dengue serotype confers lasting immunity to that serotype but not to others. When a
new serotype rises to dominance, the population has low immunity to it (accumulated
immunity is to the previously dominant serotype), creating a large effective susceptible
pool that elevates transmission. As the new dominant serotype persists, immunity
accumulates and Rt declines. If this mechanism is an important driver, f_residual
should show a specific temporal profile around serotype switches: a rise near the time
of the switch, followed by a gradual decline during stable dominance.

**Procedure:**

1. **Aggregate to monthly:** For each posterior draw, average the weekly f_residual
   values within each calendar month to obtain monthly f_residual. Compute posterior
   median and 80%/95% credible intervals of the monthly residual. Restrict to
   2013-2022 (the period with serotype data).

2. **Prepare serotype data:** Combine the 2013-2016 and 2017-2022 serotype proportion
   files into a single monthly series. Smooth the raw proportions using a GAM
   multinomial logistic regression (mgcv, thin plate regression splines) following
   Finch et al. (2025, Nat Comms). Compute dominant serotype identity from the
   smoothed proportions (argmax). Identify serotype switches as months where the
   smoothed dominant serotype identity changes. Compute Shannon entropy on both raw
   and smoothed proportions.

3. **Panel plot (4 panels, aligned x-axes):**
   - Top: monthly f_residual posterior median with 80%/95% credible bands
   - Second: GAM-smoothed serotype proportions (stacked area)
   - Third: dominant serotype identity as a color-coded bar (D1-D4)
   - Bottom: Shannon entropy (raw as points, smoothed as line)
   - Vertical lines marking serotype switches from GAM-smoothed proportions

4. **Temporal profile around switch events:** For each identified serotype switch,
   extract the posterior distribution of monthly f_residual in a window around the
   switch (e.g., 6 months before to 12 months after). Plot these as posterior
   ribbon plots centered on the switch date. The predicted pattern is:
   - Pre-switch months: f_residual near zero or negative (immunity to the incumbent
     serotype has accumulated, suppressing transmission below what covariates predict)
   - Post-switch months: f_residual positive (fresh susceptible pool elevates
     transmission beyond what covariates predict)
   - Gradual decline in f_residual over the months following the switch as immunity
     to the new dominant serotype builds

   The credible intervals indicate whether the shift is distinguishable from GP
   uncertainty. A clear positive shift with non-overlapping before/after credible
   intervals is strong evidence; a shift that is swallowed by wide intervals means
   the data cannot speak to the hypothesis.

5. **Interpretation.** This analysis is descriptive and hypothesis-generating, not a
   formal hypothesis test. With only ~3-4 major serotype switches in 10 years, there
   is no meaningful sample size for frequentist testing. The evidence comes from
   whether the temporal profile of f_residual matches the specific pattern predicted
   by the immunity mechanism — in both direction (positive after switch) and timing
   (aligned with observed serotype data) — across the individual events. If the
   residual GP is instead flat or shows fluctuations unrelated to serotype dynamics,
   that is evidence against this mechanism being an important driver at the national
   level.


### 8. Sensitivity Analyses

#### 8.1 Generation Interval

The GI is constructed from component distributions (see §3.2). Sensitivity analysis
varies the EIP mean to reflect Singapore's temperature range, which is the dominant
source of GI uncertainty:

| Scenario | EIP mean (days) | Approx. temp | Total GI mean | Rationale                                |
|----------|----------------|-------------|---------------|-------------------------------------------|
| Short    | 6.5            | ~30°C       | ~17d (~2.4 wk)| Chan & Johansson (2012): 6.5d at 30°C    |
| Baseline | 11             | ~27°C       | ~21d (~3.0 wk)| Chan & Johansson (2012): interpolated     |
| Long     | 15             | ~25°C       | ~25d (~3.6 wk)| Chan & Johansson (2012): 15d at 25°C     |

#### 8.2 Climate Lag

Test lags of 2, 4, and 6 weeks between meteorological covariates and case counts.
The 4-week lag is the baseline, motivated by the extrinsic incubation period (~10 days)
plus the intrinsic incubation period (~5 days) plus reporting delay (~5-7 days),
totaling approximately 3-4 weeks.

#### 8.3 GP Hyperparameter Priors

Vary the prior on rho (length scale) to assess sensitivity. For Models 2+
(baseline: log(rho) ~ Normal(log(15), 0.5), median 15 weeks):
- Shorter: log(rho) ~ Normal(log(10), 0.5), median 10 weeks
- Baseline: log(rho) ~ Normal(log(15), 0.5), median 15 weeks
- Longer: log(rho) ~ Normal(log(20), 0.5), median 20 weeks

#### 8.4 Overdispersion Structure

Compare negative binomial vs Poisson observation model (phi -> infinity) to verify
that overdispersion is warranted by the data.


### 9. Key References

**Rt estimation and renewal equation:**
- Cori A, Ferguson NM, Fraser C, Cauchemez S (2013). A new framework and software to
  estimate time-varying reproduction numbers during epidemics. Am J Epidemiol, 178(9),
  1505-1512.
- Abbott S, Hellewell J, Thompson RN et al. (2020). Estimating the time-varying
  reproduction number of SARS-CoV-2 using national and subnational case counts.
  Wellcome Open Res, 5:112.
- Gostic KM, McGough L, Baskerville EB et al. (2020). Practical considerations for
  measuring the effective reproductive number, Rt. PLoS Comput Biol, 16(12), e1008409.

**Gaussian processes in epidemiology:**
- Riutort-Mayol G, Burkner PC, Andersen MR, Solin A, Vehtari A (2023). Practical
  Hilbert space approximate Bayesian Gaussian processes for probabilistic programming.
  Stat Comput, 33, 17.
- Bhatt S, Gething PW, Brady OJ et al. (2013). The global distribution and burden
  of dengue. Nature, 496(7446), 504-507.

**Dengue generation interval / serial interval:**
- Chan M, Johansson MA (2012). The incubation periods of dengue viruses. PLoS One,
  7(11), e50972.
- Aldstadt J, Yoon IK, Tannitisupawong D et al. (2012). Space-time analysis of
  hospitalised dengue patients in rural Thailand reveals important temporal intervals
  in the pattern of dengue virus transmission. Trop Med Int Health, 17(9), 1076-1085.

**Singapore dengue epidemiology:**
- Ang LW, Cutter J, Lee V (2018). Dengue in Singapore from 2004 to 2016: cyclical
  epidemic patterns dominated by serotypes 1 and 2. Am J Trop Med Hyg, 99(1), 204-210.
- Tan LK, Low SL, Sun H et al. (2019). Force of infection and true infection rate of
  dengue in Singapore: implications for dengue control and management. Am J Epidemiol,
  188(8), 1529-1538.

**Wolbachia intervention:**
- Lim JT, Bansal S, Chong CS et al. (2024). Efficacy of Wolbachia-mediated sterility
  to reduce the incidence of dengue: a synthetic control study in Singapore. Lancet
  Microbe, 5, e422-e432.
- Ong J, Ho SH, Soh SXH et al. (2020). Gravitrap deployment for adult Aedes aegypti
  surveillance and its impact on dengue cases. PLoS Negl Trop Dis, 14(8), e0008528.

**Climate and dengue in Singapore:**
- Lim JT et al. (2025). Climate variation and serotype competition drive dengue
  outbreak dynamics in Singapore. Nat Commun (in press / 2025).

**Stan and Bayesian workflow:**
- Carpenter B, Gelman A, Hoffman MD et al. (2017). Stan: a probabilistic programming
  language. J Stat Softw, 76(1).
- Gabry J, Simpson D, Vehtari A, Betancourt M, Gelman A (2019). Visualization in
  Bayesian workflow. J R Stat Soc A, 182(2), 389-402.
