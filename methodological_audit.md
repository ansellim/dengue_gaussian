# Systematic Methodological Audit: GP Renewal Equation Model vs EpiNow2

## Overview

This document provides a systematic comparison of every methodological choice in our dengue Rt estimation model against EpiNow2 (Abbott et al. 2020), the most widely used Bayesian Rt estimation framework. For each choice, we document: (1) what EpiNow2 does, (2) what we do, (3) why we differ (or agree), and (4) the implications.

---

## 1. Core Framework

### 1.1 Renewal Equation

| Aspect | EpiNow2 | Our Model |
|--------|---------|-----------|
| **Equation** | λ_t = R_t × Σ c_{t-s} × w_s | λ_t = R_t × Σ c_{t-s} × w_s |
| **Cases used** | Latent infections (back-calculated) | Observed cases directly |
| **Infectious pressure** | Estimated infections convolved with GI | Observed cases convolved with GI |

**Key difference:** EpiNow2 operates on *latent infections* that are linked to observations through a delay distribution and day-of-week effects. Our model uses *observed cases directly* in the renewal equation.

**Justification for our choice:**
- At **weekly resolution**, the distinction between latent infections and observed cases is less critical because reporting delays (~5-7 days for dengue in Singapore) are shorter than or equal to the time step
- Singapore has mandatory dengue notification, so the case-to-notification delay is short and stable
- A constant reporting fraction cancels in the renewal equation Rt ratio (Cori et al. 2013)
- The negative binomial observation model absorbs residual reporting noise

**Potential critique & response:**
- *"Using observed cases in the infectious pressure biases Rt downward during epidemic growth."*
  Response: This is the same approach as Cori et al. (2013). The bias is small at weekly resolution because the reporting delay (5-7d) is within one time step (7d). EpiNow2's latent infection approach is more important for daily data with multi-day delays.

### 1.2 Temporal Resolution

| Aspect | EpiNow2 | Our Model |
|--------|---------|-----------|
| **Resolution** | Daily (primary use case) | Weekly |
| **Rationale** | Designed for COVID-19 real-time surveillance | Weekly dengue data from MOH |

**Justification:**
- Dengue data in Singapore are reported weekly (MOH Weekly Infectious Diseases Bulletin)
- Weekly resolution absorbs day-of-week effects and short reporting delays
- Trade-off: loses within-week dynamics but gains stability and simplifies the delay model

**Implication:** Many EpiNow2 features (day-of-week effects, explicit reporting delay distributions, right-truncation adjustment) are unnecessary at weekly resolution. This is not a limitation but a deliberate scope choice.

---

## 2. GP Prior on log(Rt)

### 2.1 Kernel Choice

| Aspect | EpiNow2 | Our Model |
|--------|---------|-----------|
| **Default kernel** | Matérn 3/2 | Matérn 3/2 |
| **Alternatives tested** | Matérn 5/2, SE (user can specify) | Matérn 1/2, 3/2, 5/2, SE (systematic comparison) |
| **Sensitivity analysis** | Not standard | Full cross-kernel comparison (Section 3.4 of report) |

**Our addition:** We go beyond EpiNow2 by systematically comparing 4 kernels and showing ΔELPD < 2 across all, confirming the choice is not critical. This is a strength.

### 2.2 GP Approximation

| Aspect | EpiNow2 | Our Model |
|--------|---------|-----------|
| **Method** | HSGP (Hilbert Space GP) since v1.4+ | HSGP |
| **Basis functions** | M = 17 (default) to ~40 | M = 110 |
| **Boundary factor** | L = 1.5 (default) | L = 1.5 |

**Key difference:** We use many more basis functions (M=110 vs EpiNow2's default ~17-40).

**Justification:**
- Our time series is longer (574 weeks vs typical COVID dashboards of ~50-200 days)
- M must satisfy M ≥ ⌈√3/ρ_min × 2L/π⌉ where ρ_min is the 5th percentile of the length-scale prior
- At ρ_min ≈ 2.2 weeks and our domain length, M_min ≈ 100; we use 110 for safety margin
- EpiNow2's lower M is appropriate for shorter time series

### 2.3 GP Hyperparameter Priors

| Parameter | EpiNow2 Default | Our Model | Comment |
|-----------|----------------|-----------|---------|
| **Length scale (ρ)** | log(ρ) ~ N(log(21), 0.7) [~21 days] | log(ρ) ~ N(log(6), 0.5) [~6 weeks = 42 days] | Different scales (days vs weeks) |
| **Amplitude (α)** | log(α) ~ N(0, 0.3) [~1.0] | log(α) ~ N(-1.2, 0.5) [~0.3] | More conservative amplitude |
| **GP coefficients** | Non-centered | Non-centered | Same |

**Effective comparison (converting to same units):**
- EpiNow2 ρ: median ~21 days (3 weeks), 95% CI [~6d, ~72d]
- Our ρ: median ~42 days (6 weeks), 95% CI [~15d, ~113d]
- Our prior is wider and centred on a longer length scale, but the *posterior* converges to ~20 days (2.86 weeks), which is close to EpiNow2's prior median

**Amplitude comparison:**
- EpiNow2 α: median ~1.0 on log(Rt) scale → implies ±1 excursions, i.e., Rt from 0.37 to 2.72
- Our α: median ~0.3 → implies ±0.3 excursions, i.e., Rt from 0.74 to 1.35
- Our prior is more conservative; the posterior (α ≈ 0.30) is well within our prior

**Justification for tighter amplitude:**
- Dengue Rt in Singapore shows moderate variation (0.5-2.0), less extreme than COVID-19
- The tighter prior prevents GP overfitting to noise in the tails
- Follows Riutort-Mayol et al. (2023) recommendation for moderately informative GP priors

### 2.4 IBM Prior (Extension)

EpiNow2 does not offer an IBM (Integrated Brownian Motion) prior alternative. Our extension (Section 2.7 of report) shows the IBM gives identical results with O(N) cost and no boundary artifacts. This provides additional robustness evidence not available in EpiNow2 analyses.

---

## 3. Observation Model

### 3.1 Likelihood

| Aspect | EpiNow2 | Our Model |
|--------|---------|-----------|
| **Distribution** | Negative binomial | Negative binomial |
| **Parameterisation** | NegBin(λ, φ) | NegBin_2(λ, φ) |
| **Mean-variance** | Var = λ + λ²/φ | Var = λ + λ²/φ |

**Same choice.** Both use the "NB2" parameterisation where φ controls overdispersion.

### 3.2 Overdispersion Prior

| Aspect | EpiNow2 | Our Model |
|--------|---------|-----------|
| **Prior** | 1/√φ ~ half-N(0, 1) (reciprocal sqrt) | φ ~ half-N(0, 10) |

**Difference:** EpiNow2 places the prior on the reciprocal square root of overdispersion, which induces a heavy-tailed prior on φ itself. Our prior on φ directly is simpler but achieves a similar effect.

**Note:** The specification.md states "This follows EpiNow2's approach" for the overdispersion prior, but the actual parameterisation differs slightly. The substantive effect is similar: both allow the data to determine overdispersion level without strong constraint.

---

## 4. Generation Interval

### 4.1 Construction Method

| Aspect | EpiNow2 | Our Model |
|--------|---------|-----------|
| **Input** | User provides parametric distribution (e.g., Gamma/LogNormal mean+SD) | Monte Carlo convolution of 4 biological components |
| **Uncertainty** | Can propagate GI uncertainty via bootstrapping or joint estimation | Fixed (no uncertainty propagation) |
| **Components** | Single distribution (serial interval or GI) | IIP + viremia + EIP + biting interval |
| **Temperature dependence** | Not built-in | EIP varies by temperature (multi-country extension) |

**Key difference:** Our GI construction is more mechanistic (built from biology), while EpiNow2 takes a published estimate as input.

**Justification:**
- The component-based approach is more transparent and allows temperature adjustment
- For dengue, the GI involves a human-mosquito-human cycle, making component decomposition natural
- The trade-off is that we do not propagate GI uncertainty into Rt estimates

**Potential critique & response:**
- *"Fixing the GI ignores an important source of uncertainty."*
  Response: Valid criticism. However, sensitivity analysis across 3 EIP scenarios (6.5d, 11d, 15d) shows qualitatively stable results. Joint estimation of GI and Rt is possible but creates identifiability challenges with weekly data (only 6 GI bins). EpiNow2's uncertainty propagation is a genuine advantage when GI is poorly known.

### 4.2 Generation Interval vs Serial Interval

| Aspect | EpiNow2 | Our Model |
|--------|---------|-----------|
| **What's used** | Either SI or GI (user choice) | GI (biological components) |
| **Distinction** | Explicit — different modelling implications | GI approximated by SI at weekly resolution |

At weekly resolution, the GI/SI distinction matters less because both are discretised to the same few bins. Our S=6 weekly bins capture the full GI mass.

---

## 5. Reporting Delay Model

### 5.1 Delay Convolution

| Aspect | EpiNow2 | Our Model |
|--------|---------|-----------|
| **Delay modelling** | Explicit convolution: cases = infections ⊗ delay_dist | No explicit delay model |
| **Delay distributions** | Incubation period + reporting delay (parametric, can be uncertain) | Absorbed into weekly resolution |
| **Right-truncation** | Adjusts for delayed reporting of recent cases | Not needed (weekly aggregation) |
| **Day-of-week effects** | Modelled via random effects | Not applicable (weekly data) |

**This is the largest structural difference from EpiNow2.**

**Justification for omitting delay convolution:**
1. **Weekly resolution absorbs delays:** Dengue notification delay in Singapore is ~5-7 days, which is within one weekly time step. When data are aggregated weekly, most cases appear in the correct week.
2. **Stable reporting process:** Singapore's mandatory notification system produces consistent delays, unlike the variable COVID-19 reporting that motivated EpiNow2's delay model.
3. **Cori et al. (2013) precedent:** The original EpiEstim approach for weekly data does not model reporting delays explicitly.
4. **Lau et al. (2022):** Singapore dengue Rt estimation at weekly resolution without explicit delays.

**When our approach would be problematic:**
- Daily data with 7+ day delays (COVID-19 scenario)
- Settings with highly variable reporting (weekend effects, backlogs)
- Real-time "nowcasting" where right-truncation matters
- These don't apply to our weekly historical dengue data

---

## 6. Additive Decomposition (Our Innovation)

### 6.1 Model Structure

| Aspect | EpiNow2 | Our Model |
|--------|---------|-----------|
| **log(Rt) structure** | GP only (single smooth function) | μ + β·X_climate + f_residual(t) |
| **Covariates** | Not standard (can be added manually) | Built into model structure |
| **Variance decomposition** | Not available | Var(f_climate)/Var(total) computed per posterior draw |
| **Causal interpretation** | None (GP conflates all drivers) | Partial (conditional on climate, what remains?) |

**This is the core methodological innovation of our approach.**

EpiNow2 estimates a single smooth Rt trajectory. Our model asks: how much of Rt variation is attributable to measured drivers (climate) versus unmeasured drivers (captured by the residual GP)?

**Strengths over EpiNow2:**
- Enables attribution of Rt variation to specific drivers
- The residual GP has interpretable structure (can be compared to serotype data post-hoc)
- Climate effects are estimated with proper uncertainty

**Limitations vs EpiNow2:**
- More complex model specification
- Requires careful prior specification to avoid confounding between climate and GP
- The GP can still absorb climate-driven variation if the length scale is short enough
- Climate explains <1% in our case, making the decomposition somewhat trivial

### 6.2 Post-hoc Serotype Analysis

EpiNow2 has no analogue of our post-hoc serotype analysis. This is a novel contribution: using the GP residual as a probe for unmeasured drivers. However, the approach is inherently descriptive (n=7 events, no formal testing).

---

## 7. Inference

### 7.1 Computational Framework

| Aspect | EpiNow2 | Our Model |
|--------|---------|-----------|
| **Backend** | Stan (CmdStan or RStan) | Stan (CmdStan) |
| **Sampler** | HMC/NUTS | HMC/NUTS |
| **Chains/iterations** | 4 chains, 1000+1000 (default) | 4 chains, 1000+1000 |
| **Parameterisation** | Non-centered GP | Non-centered GP |
| **Diagnostics** | R-hat, ESS, divergences | R-hat, ESS, divergences |

**Same approach.** Both use Stan's HMC/NUTS with standard diagnostics.

### 7.2 Model Comparison

| Aspect | EpiNow2 | Our Model |
|--------|---------|-----------|
| **Method** | Not built-in (user can compute LOO) | LOO-CV via PSIS (Vehtari et al. 2017) |
| **Forecasting evaluation** | Not built-in | CRPS, calibration, skill scores vs baselines |
| **Cross-validation** | Not standard | Cross-kernel LOO comparison |

**Our addition:** Systematic model comparison is a strength, including the cross-kernel LOO comparison and forecasting evaluation against meaningful baselines (random walk, not just Rt=1).

---

## 8. Features in EpiNow2 That We Deliberately Omit

### 8.1 Right-Truncation Adjustment
- **What it does:** Adjusts for the fact that recent cases may not yet be reported
- **Why we omit it:** Our data is historical (2012-2022), not real-time. All cases are fully reported. Even for real-time use, weekly aggregation reduces truncation effects.

### 8.2 Day-of-Week Effects
- **What it does:** Models systematic variation in reporting by day of week
- **Why we omit it:** Weekly data has no day-of-week variation by construction.

### 8.3 Latent Infection Back-Calculation
- **What it does:** Estimates a latent infection time series that is then convolved with delays to produce observed cases
- **Why we omit it:** At weekly resolution with stable short delays, the distinction between infections and cases is minimal. Adding a latent infection layer would increase model complexity without meaningful improvement.

### 8.4 Multiple Delay Distributions
- **What it does:** Chains incubation period + reporting delay distributions
- **Why we omit it:** Same as 8.3 — unnecessary at weekly resolution with Singapore's stable notification system.

### 8.5 Gaussian Process Random Walk Option
- **What it does:** EpiNow2 can use a random walk instead of GP for Rt dynamics
- **Why we don't use it:** The GP is our primary object of interest (decomposition, spectral analysis). A random walk would be less informative for our research questions.

---

## 9. Features in Our Model That EpiNow2 Lacks

1. **Additive decomposition** of log(Rt) into covariate + residual components
2. **Cross-kernel comparison** (4 kernels systematically evaluated)
3. **IBM prior alternative** (O(N) Markov state-space, no boundary artifacts)
4. **Component-based GI construction** (mechanistic, temperature-adjustable)
5. **Post-hoc serotype analysis** of GP residual
6. **Variance decomposition** (per-draw climate vs residual attribution)
7. **Spectral analysis** of GP residual for multi-scale temporal structure
8. **Branching process analysis** linking Rt to extinction/outbreak probability
9. **Counterfactual analysis** (excess cases with GP zeroed out)
10. **Multi-country application** with temperature-adjusted GI (12 countries)

---

## 10. Summary: When to Use Which Approach

### Use EpiNow2 when:
- Real-time surveillance with daily data
- Reporting delays are substantial and variable
- Right-truncation adjustment is needed for nowcasting
- A single Rt trajectory is sufficient (no decomposition needed)
- Quick, standardised Rt estimation is the goal
- Day-of-week reporting effects are present

### Use our approach when:
- Research question involves attribution of Rt to specific drivers
- Weekly data where delays are absorbed
- Post-hoc analysis of residual temporal structure is desired
- Comparison across kernel specifications is important
- Forecasting evaluation against meaningful baselines is needed
- Multi-country application with biology-based GI is desired

---

## 11. Anticipated Examiner Questions and Responses

### Q1: "Why didn't you just use EpiNow2?"

**Response:** EpiNow2 estimates a single smooth Rt trajectory, conflating all sources of variation. Our research question requires *decomposing* Rt into measured (climate) and unmeasured components. The additive structure log(Rt) = μ + β·X + f_residual is not available in EpiNow2. We share the same core framework (renewal equation, Matérn 3/2 GP, negative binomial, HSGP, HMC/NUTS) but extend it with covariate decomposition and post-hoc residual analysis.

### Q2: "Isn't your model just EpiNow2 with covariates?"

**Response:** Structurally, yes — and that's the point. We deliberately build on the same methodological foundation (renewal equation + GP) so that differences in results are attributable to the decomposition, not to different modelling frameworks. The innovations are: (1) additive decomposition, (2) systematic kernel comparison, (3) IBM prior alternative, (4) component-based GI, and (5) post-hoc serotype analysis of the residual.

### Q3: "Why don't you model reporting delays like EpiNow2?"

**Response:** EpiNow2's delay model was designed for daily COVID-19 data where reporting delays of 7-14+ days create substantial right-truncation and day-of-week effects. For weekly dengue data with ~5-7 day notification delays, the delay is within one time step and absorbed by weekly aggregation. Adding a delay model would increase complexity without improving inference. This follows the approach of Cori et al. (2013) and Lau et al. (2022) for weekly dengue Rt.

### Q4: "How do you know the GP isn't just absorbing what a delay model would capture?"

**Response:** The GP length scale (ρ ≈ 3 weeks) is too long to capture within-week reporting variation — it tracks multi-week outbreak dynamics. If reporting delays were biasing Rt, we would expect systematic autocorrelation in posterior predictive residuals at lag 1 week, which we do not observe. The 97.7% posterior predictive coverage suggests the model is well-calibrated without explicit delay modelling.

### Q5: "EpiNow2 propagates GI uncertainty. Don't you lose something by fixing it?"

**Response:** Yes, this is a genuine limitation. We partially address it with sensitivity analysis across 3 EIP scenarios (6.5d, 11d, 15d at 27°C), which shows qualitatively stable results. The component-based construction provides transparency about which biological parameter drives GI uncertainty (primarily the EIP). Full joint estimation of GI and Rt would be a valuable extension but creates identifiability challenges with only S=6 weekly bins.

### Q6: "Your amplitude prior (α ~ 0.3) is much tighter than EpiNow2's default (α ~ 1.0). Doesn't this constrain the GP?"

**Response:** The posterior α ≈ 0.30 (95% CrI: 0.27-0.34) is well within our prior, indicating the data are informative and the prior is not constraining. EpiNow2's wider default (α ~ 1.0) is appropriate for COVID-19 where Rt can change from 3 to 0.5 rapidly. For endemic dengue in Singapore, Rt fluctuates between ~0.5 and ~2.0, corresponding to log(Rt) excursions of ±0.7 — within 2-3 standard deviations of our α prior. Sensitivity analysis with a tighter GP prior confirms robustness.

### Q7: "How does your forecasting compare to EpiNow2's?"

**Response:** EpiNow2 is not primarily a forecasting tool — it's designed for retrospective Rt estimation. Neither EpiNow2 nor our model should be expected to provide good forecasts because the GP prior is a smoothing prior, not a mechanistic forward model. Our explicit forecasting evaluation (GP loses to random walk at all horizons) is actually a contribution: it quantifies a limitation that applies equally to EpiNow2 but is rarely tested. Forecasting dengue requires models with explicit vector dynamics and immunity tracking.

### Q8: "Would your results change if you used EpiNow2 instead?"

**Response:** For the single-GP Rt trajectory: unlikely to differ materially, given we use the same kernel, similar HSGP approximation, and same observation model. The Rt estimates would be nearly identical. The difference is that EpiNow2 would give us *one* Rt trajectory with no decomposition, while our approach gives us the same Rt trajectory *plus* the ability to attribute variance to climate and identify residual structure for post-hoc serotype analysis. The additional structure is the point.

---

## 12. Checklist: Every Methodological Choice Audited

| # | Choice | Same as EpiNow2? | Justified? | Evidence |
|---|--------|-------------------|------------|----------|
| 1 | Renewal equation | ✅ Yes | ✅ | Standard framework |
| 2 | Matérn 3/2 kernel | ✅ Yes (default) | ✅ | Cross-kernel comparison shows robustness |
| 3 | HSGP approximation | ✅ Yes | ✅ | M=110 appropriate for 574-week series |
| 4 | Neg binomial obs model | ✅ Yes | ✅ | Overdispersion well-calibrated (φ≈6) |
| 5 | Non-centered parameterisation | ✅ Yes | ✅ | Standard for GP in Stan |
| 6 | No delay convolution | ❌ Differs | ✅ | Weekly resolution absorbs 5-7d delays |
| 7 | No right-truncation | ❌ Differs | ✅ | Historical data, fully reported |
| 8 | No day-of-week effects | ❌ Differs (N/A) | ✅ | Weekly data, not applicable |
| 9 | Additive decomposition | ❌ Novel | ✅ | Core research question requires it |
| 10 | Fixed GI (component-based) | ❌ Differs | ⚠️ Partial | Sensitivity tested but no joint estimation |
| 11 | Weekly resolution | ❌ Differs | ✅ | Data availability; appropriate for dengue |
| 12 | Observed cases in renewal eq | ❌ Differs | ✅ | Standard for weekly data (Cori 2013) |
| 13 | GP amplitude prior (α~0.3) | ❌ Tighter | ✅ | Posterior within prior; dengue-appropriate |
| 14 | GP length scale prior (ρ~6wk) | ❌ Longer | ✅ | Posterior converges to ~3wk regardless |
| 15 | Overdispersion prior | ≈ Similar | ✅ | Different parameterisation, similar effect |
| 16 | Climate covariates in model | ❌ Novel | ✅ | Enables decomposition |
| 17 | Post-hoc serotype analysis | ❌ Novel | ✅ | Avoids identifiability issues |
| 18 | Cross-kernel comparison | ❌ Novel | ✅ | Demonstrates robustness |
| 19 | IBM prior alternative | ❌ Novel | ✅ | Computational + robustness benefit |
| 20 | Constant ascertainment | ✅ Same assumption | ⚠️ Partial | Cancels in renewal eq; serotype-varying rho is a limitation |

**Legend:** ✅ = Fully justified | ⚠️ = Justified with acknowledged limitations | ❌ = Not justified (none)

---

## 13. Key Talking Points for Presentation

1. **"We build on the EpiNow2 framework"** — Same renewal equation + GP + NegBin + HSGP foundation. This is intentional: we want to extend, not replace.

2. **"The key innovation is decomposition"** — EpiNow2 gives one smooth Rt. We give the same Rt decomposed into climate + residual, enabling attribution.

3. **"We omit delay modelling for good reason"** — Weekly dengue data with stable 5-7d delays ≠ daily COVID data with variable 7-14d delays. Cori et al. (2013) and Lau et al. (2022) do the same.

4. **"Our robustness checks go further"** — Cross-kernel (4 kernels), IBM prior, sensitivity analyses. Most EpiNow2 applications use defaults without systematic validation.

5. **"The forecasting result is honest"** — GP-based Rt estimation (EpiNow2 included) provides no forecasting skill. We quantify this explicitly; most studies don't.
