# Analysis Log

## 2026-03-24: Post-hoc Serotype Analysis

### Data

- Serotype proportions (D1–D4): monthly, 2013-01 to 2022-12 (120 months)
- Combined from `serotype_props_2013_2016.csv` (2013–2016) and `monthly_sero_type_props_all_data.csv` (2017–2022)
- GAM multinomial logistic regression used to smooth raw serotype proportions before computing dominant serotype (reduces 23 raw switches to 7 persistent switches)

### Identified Serotype Switches (GAM-smoothed)

| Date | Transition |
|------|-----------|
| Mar 2013 | DENV-2 → DENV-1 |
| Nov 2015 | DENV-1 → DENV-2 |
| Sep 2016 | DENV-2 → DENV-1 |
| Dec 2016 | DENV-1 → DENV-2 |
| Feb 2020 | DENV-2 → DENV-3 |
| Jun 2020 | DENV-3 → DENV-2 |
| Mar 2021 | DENV-2 → DENV-3 |

Major persistent switches: DENV-2→DENV-1 (2013), DENV-1→DENV-2 (2015), DENV-2→DENV-3 (2020/2021). The 2016 switches are rapid back-and-forth during a transitional period.

### Panel Plot (serotype_panel.png)

**Top panel — Residual GP**: f_residual shows clear oscillations on ~2–3 month timescales with amplitude ~±0.5 on log(Rt) scale. Notably elevated in 2014 (peak DENV-1 epidemic), 2016 (DENV-1/2 transition), and 2020 (DENV-3 emergence). The widening credible bands after 2020 coincide with less certain serotype dynamics.

**Middle panel — Serotype proportions**: DENV-1 dominated 2013–2015, DENV-2 dominated 2016–2019, DENV-3 emerged in 2020 and became dominant by 2021. DENV-4 was present at low levels throughout.

**Bottom panel — Shannon entropy**: Entropy peaks during transition periods (2015–2016 DENV-1/2 transition, 2020 DENV-2/3 transition), indicating higher serotype diversity. These entropy peaks approximately coincide with periods of elevated residual GP, consistent with the hypothesis that serotype diversity increases the effective susceptible pool.

### Switch Profile Analysis (serotype_switch_profiles.png)

The temporal profiles of f_residual around each switch event show **mixed evidence** for the serotype-immunity hypothesis:

**Supporting evidence:**
- **DENV-2 → DENV-3 (Feb 2020)**: Clear rise in f_residual around the switch, with sustained positive values post-switch. This is the strongest signal.
- **DENV-2 → DENV-3 (Mar 2021)**: f_residual rises sharply post-switch, consistent with a fresh susceptible pool driving excess transmission.
- **DENV-2 → DENV-1 (Mar 2013)**: f_residual is positive in the months following the switch.

**Ambiguous/non-supporting:**
- **DENV-1 → DENV-2 (Nov 2015)**: f_residual was already declining before the switch and continued declining after — no clear post-switch rise.
- **DENV-2 → DENV-1 (Sep 2016)**: f_residual drops sharply post-switch, opposite to prediction.
- **DENV-3 → DENV-2 (Jun 2020)**: Transient switch during a volatile period; profile is noisy.

### Interpretation

The residual GP captures substantial temporal structure that aligns partially with serotype dynamics, particularly the 2020–2021 DENV-3 emergence. However, the pattern is not consistent across all switch events. With only ~3–4 major switches over 10 years, the evidence is **suggestive but not conclusive** that serotype-driven immunity dynamics are a major driver of the residual Rt variation.

The residual GP likely reflects a combination of serotype dynamics, population immunity shifts, spatial heterogeneity in transmission, and other unmeasured factors — not serotype dynamics alone.

### Figures

- `results/figures/serotype_panel.png` — Three-panel overview (residual GP, serotype proportions, Shannon entropy)
- `results/figures/serotype_switch_profiles.png` — Event-centered temporal profiles of f_residual around each switch
- `results/serotype_residual_monthly.csv` — Monthly analysis data

---

## 2026-03-24: Model Fitting and Diagnostics

### Configuration

| Setting | Value |
|---------|-------|
| Model | Climate-only: log(Rt) = mu + beta_temp*temp + beta_rain*rain + f_residual(t) |
| Chains | 4 |
| Warmup / Sampling | 1000 / 1000 |
| adapt_delta | 0.95 |
| max_treedepth | 12 |
| M (HSGP basis functions) | 200 |
| Time input | Centered weeks (unscaled), L = 429.75 |
| Generation interval | Component-based: IIP (5.9d) + viremia (Unif 0–5d) + EIP (11d at 27°C) + bite delay (2d); mean 21.4d, S = 6 |

### Priors

| Parameter | Prior | Interpretation |
|-----------|-------|----------------|
| mu | Normal(0, 0.5) | Baseline Rt median 1.0, 95% [0.37, 2.72] |
| beta_temp | Normal(0, 0.5) | Per 1-SD temperature change |
| beta_rain | Normal(0, 0.5) | Per 1-SD rainfall change |
| log_alpha | Normal(-1.2, 0.5) | GP amplitude median 0.30, 95% [0.11, 0.80] |
| log_rho | Normal(log(6), 0.5) | GP length scale median 6 weeks, 95% [2.2, 16.0] |
| phi | HalfNormal(0, 10) | NegBin overdispersion |

### MCMC Diagnostics

Zero divergences, no treedepth issues, all Rhat < 1.01, minimum ESS > 800.

### GP Hyperparameter Posteriors

| Parameter | Posterior median | 95% CrI |
|-----------|------------------|---------|
| alpha (amplitude) | 0.30 | [0.27, 0.34] |
| rho (length scale, weeks) | 2.86 | [2.09, 3.71] |

Posterior rho (~3 weeks) is substantially shorter than the prior median (6 weeks). The data strongly prefer a short GP correlation length.

### Effect Sizes — Multiplicative on Rt

| Effect | Median | 95% CI |
|--------|--------|--------|
| Temperature (per 1 SD) | 0.989 | [0.969, 1.009] |
| Rainfall (per 1 SD) | 1.003 | [0.986, 1.020] |

Climate effects are negligible (95% CIs tightly bracket 1).

### Variance Decomposition

| Component | Proportion (median) | 95% CI |
|-----------|---------------------|--------|
| Climate | 0.4% | [0.0%, 1.4%] |
| Residual GP | 99.6% | [98.6%, 100%] |

The residual GP dominates, capturing nearly all Rt variation at weekly resolution.

### Posterior Predictive Check

The posterior predictive intervals track observed epidemic peaks well (2013–2014, 2016, 2019–2020, 2022). However, pointwise intervals under-cover weekly fluctuations — the negative binomial dispersion is too tight to absorb weekly shock structure. Flagged as an observation-layer limitation; the latent Rt decomposition is still credible.

### Prior Predictive Check

Week-level Rt distribution from 5,000 prior draws:
- Median Rt: 0.98
- 95% interval: [0.16, 5.96]

Prior generates plausible Rt values consistent with published Singapore range (Lim et al. 2022: 0.45–2.91). Case-level forward simulations show explosive growth in some draws due to renewal equation feedback without susceptible depletion — a known artifact, not a prior deficiency.

### Key Figures

- `results/figures/prior_predictive_summary.png` — Prior predictive Rt and cases
- `results/figures/rt_with_cases.png` — Rt overlaid on observed cases
- `results/figures/decomposition.png` — log(Rt) component decomposition
- `results/figures/posterior_predictive_timeseries.png` — Posterior predictive check
- `results/figures/effect_sizes.png` — Covariate effect size posteriors
