---
output:
  pdf_document: default
  html_document: default
---
# Serotype-Driven Dynamics of Dengue Effective Reproduction Number in Singapore: A Semiparametric Bayesian Approach

## Background and Rationale

Singapore's dengue epidemiology is shaped by drivers operating at multiple timescales, including climate variation, serotype cycling, and population immunity shifts. Standard approaches to estimating the effective reproduction number (Rt) either use a single flexible smoother that conflates all sources of variation, or purely parametric models that impose rigid functional forms. Neither approach is well suited to disentangling serotype-driven immunity dynamics from other sources of Rt variation.

This project develops a renewal equation model of weekly dengue Rt in Singapore (2012--2022) where log(Rt) is decomposed into additive components: climate covariates (temperature and rainfall at 4-week lag) and a residual Gaussian process (GP) capturing unexplained temporal structure. The central hypothesis is that serotype replacement events elevate Rt by expanding the effective susceptible pool, and that the residual GP can serve as a signal for early detection of serotype-driven outbreaks.

## Methodology

**Transmission model.** Expected weekly cases follow a discretized renewal equation, with a generation interval constructed from biological components (intrinsic incubation period, viremia window, extrinsic incubation period, and mosquito biting interval) via Monte Carlo convolution, discretized to a 6-week PMF (Chan and Johansson 2012, PLOS ONE).

**Latent process.** The log-transformed Rt is modeled as log(Rt) = mu + beta_temp * temp + beta_rain * rain + f_residual(t), where f_residual is a zero-mean GP with Matern 3/2 kernel. The GP is implemented via Hilbert Space approximation (HSGP; Riutort-Mayol et al. 2023, Stat Comput) with 110 basis functions, reducing computation from O(N^3) to O(MN) and enabling fully Bayesian inference in Stan.

**Observation model.** Cases follow a negative binomial distribution with overdispersion parameter phi, accommodating the clustered and spatially heterogeneous nature of dengue surveillance data. Singapore's mandatory notification system supports the assumption of complete symptomatic case ascertainment.

**Priors.** Weakly informative: mu ~ Normal(0, 0.5) centers baseline Rt at 1; GP length scale prior median at 6 weeks allows tracking of outbreak-scale dynamics; GP amplitude prior median at 0.3 on the log(Rt) scale permits Rt to vary by factors of roughly 1.3x to 2.2x around covariate predictions.

**Serotype analysis.** Monthly serotype proportions (DENV-1 through DENV-4, 2013--2022) are smoothed using GAM multinomial logistic regression (following Finch et al. 2025, Nat Commun). Shannon entropy and dominant serotype identity are derived from the smoothed proportions. The posterior distribution of f_residual is aggregated to monthly resolution and examined in temporal windows around serotype switch events via cross-correlation, pre-switch trend analysis, and event-centered ribbon plots.

## Key Findings to Date

**Climate effects are negligible at weekly resolution.** Temperature and rainfall coefficients are centered near zero (multiplicative effects on Rt of 0.989 [95% CrI: 0.969, 1.009] and 1.003 [0.986, 1.020], respectively), consistent with climate acting through slower vector ecology mechanisms not resolved at weekly timescales.

**The residual GP dominates Rt variation.** The climate component accounts for only 0.5% of variance in log(Rt) (median 0.4%), with the residual GP capturing effectively all structured variation. This indicates that unmeasured drivers, plausibly serotype-immunity dynamics, are the primary determinants of Rt fluctuation.

**Serotype switches are associated with elevated residual Rt.** In 5 of 7 identified switch events, the posterior probability of elevated f_residual in the 0--3 months post-switch exceeds 0.75 (range: 0.002 to 1.0). The DENV-2 to DENV-3 transitions in 2020--2021 show particularly strong signals (P(elevated) = 1.0), with median f_residual increases of 0.2--0.3 on the log scale (corresponding to 20--35% increases in Rt beyond climate-predicted levels).

**Cross-correlation suggests a leading relationship.** The CCF between f_residual and Shannon entropy peaks at lag +10 to +12 months (entropy leading), suggesting that increasing serotype diversity precedes residual Rt elevation by approximately 10 months, a potentially actionable lead time for outbreak preparedness.

**Model validation.** Posterior predictive checks show 95% coverage of 97.7% and 80% coverage of 92.3%, indicating good calibration. The model is robust to a tighter GP amplitude prior, confirming that the variance decomposition is not an artifact of prior specification.

## Current Direction

The project is shifting from a general decomposition framework toward two specific aims: (1) formally characterizing the effect of serotype replacement on Rt, examining the temporal profile of f_residual around switch events to test whether the pattern matches the immunity-depletion mechanism (rise at switch, gradual decline during stable dominance); and (2) evaluating the feasibility of serotype-based outbreak prediction, leveraging the ~10-month lead of entropy signals over Rt elevation as a candidate early warning indicator. These aims exploit the semiparametric structure of the model, where the GP residual isolates non-climate Rt variation that can be attributed post hoc to serotype dynamics without the identifiability problems of including serotype as an in-model covariate.

## References

1. Bhatt S, Gething PW, Brady OJ et al. (2013). The global distribution and burden of dengue. *Nature*, 496(7446), 504--507.
2. Chan M, Johansson MA (2012). The incubation periods of dengue viruses. *PLOS ONE*, 7(11), e50972.
3. Cori A, Ferguson NM, Fraser C, Cauchemez S (2013). A new framework and software to estimate time-varying reproduction numbers during epidemics. *Am J Epidemiol*, 178(9), 1505--1512.
4. Finch E et al. (2025). Climate variation and serotype competition drive dengue outbreak dynamics in Singapore. *Nat Commun*.
5. Lau MAY et al. (2022). New framework for reducing dengue Rt using a serial interval distribution from genomic data. *PLOS Comput Biol*.
6. Riutort-Mayol G, Burkner PC, Andersen MR, Solin A, Vehtari A (2023). Practical Hilbert space approximate Bayesian Gaussian processes for probabilistic programming. *Stat Comput*, 33, 17.
