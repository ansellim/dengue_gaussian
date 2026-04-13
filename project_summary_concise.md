---
output:
  pdf_document:
    geometry: margin=2cm
---

# Serotype-Driven Dynamics of Dengue Rt in Singapore: A Semiparametric Bayesian Approach

## Background

We develop a Bayesian renewal equation model for weekly dengue Rt in Singapore (2012--2022), decomposing log(Rt) into climate covariates (temperature and rainfall at 4-week lag) and a residual Gaussian process (GP). The residual GP isolates non-climate Rt variation that can be examined post hoc for serotype-driven immunity dynamics, avoiding the identifiability problems of including serotype as an in-model covariate.

## Methods

Expected weekly cases follow a discretized renewal equation (Cori et al. 2013) with a generation interval derived from dengue biology (Chan & Johansson 2012). The latent process is log(Rt) = $\mu$ + $\beta_{temp}$ $\cdot$ temp + $\beta_{rain}$ $\cdot$ rain + $f_{residual}(t)$, where $f_{residual}$ is a zero-mean GP with Matérn 3/2 kernel, implemented via Hilbert Space approximation (HSGP; Riutort-Mayol et al. 2023) with 110 basis functions. Cases follow a negative binomial observation model. Monthly serotype proportions (DENV-1--4) are smoothed via GAM multinomial logistic regression (Finch et al. 2025), from which Shannon entropy and switch events are derived and compared against posterior $f_{residual}$ distributions.

## Key Findings

- **Climate effects are negligible at weekly resolution.** Temperature and rainfall multiplicative effects on Rt: 0.989 [95% CrI: 0.969, 1.009] and 1.003 [0.986, 1.020]. Climate explains <1% of Rt variance.

- **The residual GP captures ~99.5% of Rt variation**, indicating that unmeasured drivers---plausibly serotype-immunity dynamics---dominate Rt fluctuation.

- **Serotype switches associate with elevated Rt.** In 5/7 switch events, P(elevated $f_{residual}$) > 0.75 in the 0--3 months post-switch. DENV-2$\to$DENV-3 transitions (2020--2021) show the strongest signals (P = 1.0; 20--35% Rt increase beyond climate prediction).

- **Serotype diversity leads Rt elevation by ~10 months.** Cross-correlation between $f_{residual}$ and Shannon entropy peaks at lag +10--12 months (entropy leading), suggesting a potentially actionable early warning window.

- **Model validation.** Posterior predictive 95% coverage: 97.7%; 80% coverage: 92.3%. Variance decomposition is robust to tighter GP priors.

## Current Direction

Two aims: (1) characterize the temporal profile of $f_{residual}$ around serotype switches to test the immunity-depletion mechanism; (2) evaluate serotype diversity (Shannon entropy) as an early warning indicator for Rt elevation, leveraging the ~10-month leading relationship.

## References

1. Chan M, Johansson MA (2012). The incubation periods of dengue viruses. *PLOS ONE*, 7(11), e50972.
2. Cori A, Ferguson NM, Fraser C, Cauchemez S (2013). A new framework and software to estimate time-varying reproduction numbers during epidemics. *Am J Epidemiol*, 178(9), 1505--1512. doi: 10.1093/aje/kwt133.
3. Finch E, Chang C, Kucharski A et al. (2025). Climate variation and serotype competition drive dengue outbreak dynamics in Singapore. *Nat Commun*, 16, 11364. doi: 10.1038/s41467-025-66411-6.
4. Riutort-Mayol G, Bürkner PC, Andersen MR et al. (2023). Practical Hilbert space approximate Bayesian Gaussian processes for probabilistic programming. *Stat Comput*, 33, 17. doi: 10.1007/s11222-022-10167-2.
