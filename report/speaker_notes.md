# Speaker Notes: Serotype-Driven Dynamics of Dengue Rt in Singapore

Total: ~15 minutes | 18 content slides + title + closing | ~1 min per slide

---

## Slide 1: Why Decompose Dengue Rt?
**[Target: 1 min | Cumulative: 1:00]**

Today I'll present our work on decomposing dengue transmission dynamics in Singapore. Dengue is the most widespread mosquito-borne disease globally, with roughly 390 million infections per year. Singapore is a particularly interesting case: despite world-class public health infrastructure, it still experiences recurrent epidemics.

The fundamental problem is that standard Rt estimation methods — Cori, EpiNow2 — give you a single smooth function. They tell you what Rt is, but not WHY it's changing. Is it climate? Serotype switching? Something else?

Our approach explicitly decomposes log Rt into climate-attributable and residual components using Gaussian Processes. We chose GPs specifically because they allow principled additive decomposition with covariates, provide Bayesian uncertainty quantification, and have interpretable hyperparameters — unlike black-box alternatives like variational autoencoders or foundation models.

---

## Slide 2: Research Questions
**[Target: 45 sec | Cumulative: 1:45]**

We address four specific research questions. First: at weekly resolution, how much of Rt variation can we actually attribute to climate? Second: does the residual GP — the part left over after accounting for climate — show patterns that align with serotype switching events? Third: are our findings sensitive to the GP kernel we choose? And fourth, practically: can this framework produce useful forecasts?

Beyond these core questions, we also extend the analysis with branching process theory, a VAE comparison, spectral analysis, and counterfactual excess case estimation. I'll touch on each of these briefly.

---

## Slide 3: Data Sources
**[Target: 45 sec | Cumulative: 2:30]**

Our data spans 11 years: January 2012 to December 2022 — that's 574 weeks. Dengue case counts come from MOH's weekly infectious diseases bulletin via data.gov.sg. Climate data — mean temperature and total rainfall — come from the Changi Airport weather station through the Meteostat API. These are standardised and lagged by 4 weeks to reflect the time climate conditions take to propagate through *Aedes* vector ecology — egg development, larval and pupal maturation, adult emergence, and biting, roughly two to three weeks — plus the human case-ascertainment delay of intrinsic incubation and onset-to-notification, roughly another two weeks.

A critical design choice: serotype surveillance data from NEA is used only in post-hoc analysis, never as a model covariate. This is deliberate — serotype proportions are derived from the same case surveillance data used as the outcome, so including them as covariates would create identifiability problems.

---

## Slide 4: Model — Semiparametric Decomposition
**[Target: 1 min | Cumulative: 3:30]**

Here's the core model. Log Rt is decomposed additively into an intercept, linear climate effects at 4-week lag, and a residual Gaussian process. The GP uses a Matérn 3/2 kernel — the same default used in EpiNow2 — because it produces once-differentiable sample paths: smooth enough to be continuous, but flexible enough to capture rapid Rt changes during outbreaks.

The key computational trick is the Hilbert Space GP approximation with 110 basis functions, which reduces complexity from cubic to linear in the number of time points. This makes HMC inference tractable in Stan.

The transmission model is a standard discretised renewal equation. Importantly, the generation interval is constructed from biological components — intrinsic incubation, viremic period, extrinsic incubation, and biting interval — rather than simply fitting a parametric distribution to observed serial intervals.

---

## Slide 5: Priors and Identifiability
**[Target: 45 sec | Cumulative: 4:15]**

Our priors are weakly informative but grounded in epidemiology. The intercept prior centres Rt at 1, appropriate for endemic dengue. The GP length scale prior has a median of 15 weeks, reflecting that the residual should capture variation on medium timescales. Climate coefficients are centred at zero — no effect. All priors were validated through prior predictive simulation to ensure they produce plausible Rt trajectories and case counts before seeing any data.

The GP amplitude prior with median 0.3 on the log-Rt scale implies typical excursions of plus-or-minus 0.3 to 0.8 — reasonable for dengue dynamics.

---

## Slide 6: Climate Explains <1% of Rt Variance
**[Target: 1 min | Cumulative: 5:15]**

This is the headline finding. Climate effects on weekly Rt are negligible. A 1-SD increase in temperature multiplies Rt by 0.989 and rainfall by 1.003 — both statistically indistinguishable from 1. The variance decomposition confirms: climate explains only 0.4 percent, with the residual GP absorbing 99.6 percent.

Now, this requires careful interpretation. This does NOT mean climate doesn't matter for dengue. Climate almost certainly acts through seasonal modulation of vector populations, but on longer timescales. At weekly resolution, faster processes dominate — and we'll argue these are serotype-driven immunity dynamics.

In the figure, you can see the GP residual in orange tracking virtually all the temporal structure, while the climate component in blue is essentially flat.

---

## Slide 7: Model Validation
**[Target: 30 sec | Cumulative: 5:45]**

Before trusting the decomposition, we validate the model. The posterior predictive check captures both the timing and magnitude of all major epidemic peaks: the 2013 DENV-1 wave, the 2019 and 2020 outbreaks, and the large 2022 epidemic. We flag that observation-layer calibration is poor — pointwise predictive intervals under-cover weekly fluctuations because the negative binomial dispersion is too tight to absorb weekly shocks. The latent Rt decomposition is credible, but the pointwise case distribution is under-specified and we note this as a limitation.

---

## Slide 8: Serotype Switches Associate with Elevated Rt
**[Target: 1 min | Cumulative: 6:45]**

Now the key question: what's driving the residual GP? This three-panel figure shows the answer. Top: the GP residual with credible intervals. Middle: smoothed serotype proportions, showing sequential dominance patterns. Bottom: Shannon entropy.

The vertical lines mark serotype switch events. Five of seven switches coincide with positive excursions in the GP residual. The strongest signal is the February 2020 DENV-2 to DENV-3 transition, with posterior probability 1.0 of elevation and a 20 to 35 percent Rt increase beyond climate prediction. This is consistent with immunity depletion: when DENV-3 re-emerged after a long absence, there was a large naive susceptible pool.

---

## Slide 9: Temporal Profiles Around Switches
**[Target: 45 sec | Cumulative: 7:30]**

Zooming in on individual switch profiles, we see a pattern consistent with immunity depletion theory: Rt rises at the switch point, then gradually declines during stable dominance as the susceptible pool is depleted by accumulating infections. This is clearest for the 2020 DENV-3 emergence but visible across most events.

I want to be transparent: with only 7 switch events over 10 years, this is necessarily descriptive. We cannot do formal hypothesis testing. But the temporal profiles match theoretical predictions in both direction — positive post-switch — and timing — aligned with observed serotype transitions.

---

## Slide 10: Entropy Leads Rt by ~10 Months
**[Target: 45 sec | Cumulative: 8:15]**

The cross-correlation analysis reveals something potentially actionable. Shannon entropy — measuring serotype diversity — leads the GP residual by approximately 10 to 12 months. This makes biological sense: rising diversity signals increasing serotype competition, which precedes susceptible pool dynamics that eventually elevate Rt.

A 10-month lead is potentially useful for seasonal risk assessment — if you see entropy rising now, you might expect elevated transmission 10 months later. However, as we'll see in the forecasting results, this doesn't translate to tactical weekly forecast improvements. The signal operates on too long a timescale.

---

## Slide 11: Kernel Comparison — Robust
**[Target: 45 sec | Cumulative: 9:00]**

A natural concern is whether our findings depend on the Matérn 3/2 kernel choice. We refitted the model with all four standard kernels: Matérn 1/2, 3/2, 5/2, and squared exponential. Left panel: Rt trajectories are nearly identical, with delta-ELPD less than 2 for all — statistically indistinguishable fits.

More importantly, on the right: all four kernels detect the same 5 of 7 serotype switches at the 0.75 probability threshold. The Matérn 1/2 is marginally preferred, consistent with dengue Rt exhibiting abrupt changes. But the scientific conclusions are completely robust to kernel specification.

---

## Slide 12: Forecasting — No Predictive Skill
**[Target: 1 min | Cumulative: 10:00]**

Here's our honest negative result. The GP never beats a random walk baseline at any forecast horizon — 4, 8, or 13 weeks. The right panel shows why: the GP posterior rapidly regresses to the prior mean of Rt approximately equals 1 within about 3 weeks, matching the estimated length scale.

This is fundamental, not an implementation issue. The GP is a smoothing prior that regularises temporal variation retrospectively. It contains no mechanistic forward model and therefore provides no information about future trajectory beyond mean reversion. For operational dengue forecasting, you would need models with explicit vector dynamics and immunity tracking. Reporting this negative result honestly is important.

---

## Slide 13: Branching Process — Vulnerability Windows
**[Target: 45 sec | Cumulative: 10:45]**

Connecting Rt estimates to outbreak risk via branching process theory, we find a strong anticorrelation between Rt and extinction probability — Spearman rho of minus 0.94. Even small Rt changes translate to large changes in outbreak sustaining probability.

Five of seven serotype switches show increased sustaining probability post-switch, with a mean increase of 6 percentage points. The 2022 outbreak is striking: 18 consecutive weeks of sustaining probability above 0.5 and mean Rt of 1.17 — the most sustained transmission in our observation period, coinciding with DENV-3 circulation.

---

## Slide 14: VAE Comparison — Meaningful Negative Result
**[Target: 45 sec | Cumulative: 11:30]**

To validate our mechanistic framework, we compared it against a variational autoencoder trained purely on case count data. The VAE completely fails to discover the GP residual structure — maximum correlation of just 0.09 with any latent dimension, and no clustering by serotype, entropy, or time period.

This is a meaningful negative result. It confirms that the renewal equation plus GP framework extracts epidemiologically meaningful information about transmission dynamics that purely data-driven pattern recognition simply cannot discover from case counts alone. The mechanistic structure of the renewal equation is doing real work here.

---

## Slide 15: Spectral Analysis & Counterfactual
**[Target: 1 min | Cumulative: 12:30]**

Two final extensions. On the left, spectral analysis of the GP residual reveals multi-scale structure matching theoretical predictions: a dominant 26-week biannual cycle from Singapore's two monsoon seasons, residual annual seasonality at 52 weeks, and importantly, peaks at 2.2 and 3.7 years matching the predicted multi-serotype cycling periods from Wearing and Rohani 2006.

On the right, counterfactual analysis quantifies excess cases per switch. Only the February 2020 DENV-3 emergence produces clearly significant excess: roughly 4,400 additional cases. Most switches produce excess indistinguishable from zero — suggesting that only truly novel serotype introductions, not familiar rotations, drive detectable excess transmission at the national level.

---

## Slide 16: Limitations
**[Target: 45 sec | Cumulative: 13:15]**

I want to be upfront about limitations. The GP dominance is partly structural: with a 3-week length scale, the GP absorbs any smooth variation — so 99.6 percent variance share doesn't prove serotype causation. Climate may be more important at seasonal timescales. The forecasting failure is fundamental, not a tuning problem.

Epidemiologically, 7 switch events give us no meaningful sample size for formal testing. Case ascertainment likely varies by serotype — estimates range from 3-to-1 to 16-to-1 depending on serotype and period. We treat Singapore as a single spatial unit, and our 11-year window provides limited power for characterising multi-year cycles.

---

## Slide 17: Conclusions & Future Directions
**[Target: 1 min | Cumulative: 14:15]**

To summarise: five key takeaways. First, climate is genuinely negligible at weekly resolution. Second, serotype switches associate with elevated Rt, with entropy providing roughly 10 months of lead time. Third, all conclusions are robust across four kernel specifications. Fourth, the GP provides no forecasting skill — an honest negative result that's important to report. Fifth, the VAE comparison validates that mechanistic modelling adds real value over purely data-driven approaches.

Looking forward, three extensions would address the most important limitations. Integrating GP decomposition with an explicit multi-serotype SIR model for causal attribution. Spatial disaggregation using dengue cluster data. And hybrid approaches combining retrospective GP decomposition with mechanistic forward models for operational forecasting. Thank you.

---

## Slide 18: Thank You / Questions
**[Target: 45 sec | Cumulative: 15:00]**

Thank you for your attention. I'm happy to take any questions. The full code repository and data pipeline are available in the submitted materials.

---

# Q&A Preparation

**Likely questions and suggested responses:**

**Q: Why not include serotype as a covariate in the model?**
A: Serotype proportions are derived from the same case surveillance data — including them creates an identifiability problem. The post-hoc approach avoids this while still allowing us to examine the association.

**Q: Isn't 99.5% GP variance decomposition just overfitting?**
A: The GP has informative priors (log-alpha, log-rho) and the posterior predictive tracks the timing and magnitude of all major peaks. The variance share is high because the GP's length scale (~3 weeks) captures dynamics that climate (seasonal timescale) cannot. But yes, the GP's flexibility means it absorbs non-climate variation whether or not it's serotype-related. We acknowledge this explicitly.

**Q: Why does climate show no effect when other studies find climate matters?**
A: Scale-dependent. Climate affects vector populations on seasonal/monthly timescales. At weekly resolution, faster processes dominate. Finch et al. 2025 found similar results — climate explains a modest fraction after serotype effects. A model at monthly or seasonal aggregation would likely give climate more weight.

**Q: Could the 10-month lag be spurious?**
A: With only 7 switch events, we can't rule it out. The lag is consistent with theoretical predictions about susceptible pool dynamics, but the sample size precludes formal testing. We present it as suggestive evidence, not a definitive finding.

**Q: Why not use EpiNow2 directly?**
A: EpiNow2 estimates Rt but doesn't decompose it into driver components. Our model adds the additive structure with explicit climate covariates plus residual GP, which enables the decomposition analysis. The GP implementation (HSGP) is similar in spirit to EpiNow2's approach.

**Q: What would you do differently?**
A: Three things: (1) jointly estimate the generation interval instead of fixing it, (2) use a multi-serotype SIR model instead of post-hoc analysis for causal attribution, and (3) incorporate spatial cluster data for neighbourhood-level analysis.

**Q: How does this compare to Finch et al. 2025?**
A: Finch et al. use a different framework — GAM-based model with explicit serotype competition terms. They find that climate explains a modest fraction after serotype effects, consistent with our finding. Our contribution is the GP decomposition approach and the specific characterisation of the residual's temporal structure and periodicities.

**Q: Why Matérn 3/2 specifically?**
A: It's the default in EpiNow2 for Rt estimation because it produces once-differentiable paths — smooth enough for continuity but rough enough for rapid changes in Rt. Our kernel comparison shows it doesn't matter much — all four kernels give essentially the same results.

**Q: Why is the generation interval fixed rather than estimated?**
A: We constructed it from well-characterised biological components (Chan & Johansson 2012) rather than fitting to data. Joint estimation would be ideal but adds complexity. Sensitivity analysis across temperature-dependent EIP scenarios shows qualitatively stable results.

**Q: How do you handle the COVID-19 period (2020-2021)?**
A: The model treats 2020-2021 the same as other periods. The large 2020 outbreak coincides with the DENV-3 emergence, not COVID-19 lockdowns per se. Lockdowns may have affected mosquito exposure patterns, but our model doesn't explicitly account for this — it's captured in the GP residual.
