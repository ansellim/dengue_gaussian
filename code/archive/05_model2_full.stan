// ==============================================================================
// Model 2 v3: Informative intervention priors based on literature
//
// Changes from v2:
// 1. Wolbachia prior: Normal(-0.3, 0.5) based on AWED trial evidence
//    - Utarini et al. (2021, NEJM) showed 77% protective efficacy
//    - With partial coverage (~30%), population effect attenuated
//    - Prior median: exp(-0.3) = 0.74 (26% reduction in Rt)
//    - Prior 95% CI: [0.28, 1.97] - still allows for no effect
//    - P(effect > 1) ~ 27% under prior - data can override
//
// 2. NPI prior: Normal(0, 0.5) - regularized but agnostic on direction
//    - Literature mixed: lockdowns reduce mobility but Aedes is peridomestic
//    - Some studies show decreased dengue, others increased
//    - Tightened from Normal(0, 1.0) to regularize extreme effects
//    - Prior 95% CI: [0.37, 2.7] on multiplicative scale
//
// Prior sensitivity: Run with v2 priors to assess robustness
// ==============================================================================

functions {
  // Spectral density of Matern 3/2 kernel in 1D
  real spd_matern32(real omega, real alpha, real rho) {
    real sqrt3_over_rho = sqrt(3.0) / rho;
    real numerator = square(alpha) * 4.0 * pow(sqrt3_over_rho, 3);
    real denominator = square(square(sqrt3_over_rho) + square(omega));
    return numerator / denominator;
  }
}

data {
  int<lower=1> N;            // Total number of weeks
  int<lower=1> N_model;      // Number of modeled weeks (N - S)
  int<lower=1> S;            // Maximum generation interval in weeks
  int<lower=1> M;            // Number of HSGP basis functions

  array[N] int<lower=0> cases;
  vector[N_model] t;
  real<lower=0> L;
  simplex[S] gi;

  int<lower=1> K_full;
  matrix[N_model, K_full] X_full;
}

transformed data {
  matrix[N_model, M] PHI;
  vector[M] sqrt_eigenvalues;

  for (j in 1:M) {
    sqrt_eigenvalues[j] = (j * pi()) / (2 * L);
  }

  for (i in 1:N_model) {
    for (j in 1:M) {
      PHI[i, j] = sin(sqrt_eigenvalues[j] * (t[i] + L)) / sqrt(L);
    }
  }
}

parameters {
  real mu;
  vector[K_full] beta;

  // Log-parameterized GP hyperparameters (better geometry)
  real log_alpha;
  real log_rho;

  vector[M] beta_gp_raw;  // Raw (non-centered) GP coefficients
  real<lower=0> phi;
}

transformed parameters {
  // Transform to natural scale
  real<lower=0> alpha = exp(log_alpha);
  real<lower=0> rho = exp(log_rho);

  vector[N_model] f_climate;
  vector[N_model] f_wolbachia;
  vector[N_model] f_npi;
  vector[N_model] f_residual;
  vector[N_model] log_Rt;
  vector<lower=0>[N_model] lambda;

  // Covariate effects
  f_climate = X_full[, 1] * beta[1] + X_full[, 2] * beta[2];
  f_wolbachia = X_full[, 3] * beta[3];
  f_npi = X_full[, 4] * beta[4];

  // Residual GP with non-centered parameterization
  {
    vector[M] spd_weights;
    for (j in 1:M) {
      spd_weights[j] = sqrt(spd_matern32(sqrt_eigenvalues[j], alpha, rho));
    }
    f_residual = PHI * (spd_weights .* beta_gp_raw);
  }

  // log(Rt)
  log_Rt = mu + f_climate + f_wolbachia + f_npi + f_residual;

  // Renewal equation with better numerical stability
  for (i in 1:N_model) {
    real infectious_pressure = 0;
    int t_idx = S + i;

    for (s in 1:S) {
      infectious_pressure += cases[t_idx - s] * gi[s];
    }

    // Use softplus-like floor to avoid hard cutoff
    lambda[i] = fmax(exp(log_Rt[i]) * infectious_pressure, 1.0);
  }
}

model {
  // Priors
  mu ~ normal(0, 0.5);

  // Climate effects (regularized, no strong directional prior)
  beta[1] ~ normal(0, 0.5);  // Temperature
  beta[2] ~ normal(0, 0.5);  // Rainfall

  // ============================================================================
  // INTERVENTION PRIORS - Literature-informed
  // ============================================================================

  // Wolbachia effect on Rt
  // Literature basis:
  //   - AWED trial (Utarini et al., 2021): 77% reduction in symptomatic dengue
  //   - Effect attenuated at population level with partial coverage
  //   - Biological mechanism: Wolbachia blocks DENV replication in Aedes aegypti
  // Prior: Normal(-0.3, 0.5)
  //   - Median multiplicative effect: exp(-0.3) = 0.74 (26% Rt reduction)
  //   - 95% prior interval: [exp(-1.28), exp(0.68)] = [0.28, 1.97]
  //   - P(increases Rt | prior) = P(beta[3] > 0) = 27%
  //   - Weakly informative: favors reduction but data can show no effect
  beta[3] ~ normal(-0.3, 0.5);

  // NPI effect on Rt (COVID-19 measures)
  // Literature basis:
  //   - Singapore 2020 was a RECORD dengue year (>35,000 cases) despite lockdowns
  //   - Aedes aegypti is peridomestic; more time at home may INCREASE exposure
  //   - Reduced outdoor mobility could decrease transmission, but evidence is mixed
  //   - Vector control activities may have been disrupted during lockdowns
  //   - Singapore Circuit Breaker (Apr-Jun 2020) coincided with outbreak peak
  // Prior: Normal(0, 0.5)
  //   - Centered at NO EFFECT - genuinely agnostic on direction
  //   - P(Rt increases | prior) = 50% - allows for increased transmission
  //   - P(Rt decreases | prior) = 50%
  //   - 95% prior interval: [0.37, 2.7] on multiplicative scale
  //   - Regularized vs Normal(0,1) to prevent implausibly extreme effects
  beta[4] ~ normal(0, 0.5);

  // GP hyperparameters (log-scale priors)
  log_alpha ~ normal(-1.6, 0.4);  // alpha median ~0.2, 95% ~[0.09, 0.44]
  log_rho ~ normal(log(6), 0.5);  // rho median 6 weeks, 95% ~[2.2, 16.0]

  // GP basis coefficients (non-centered)
  beta_gp_raw ~ std_normal();

  // Overdispersion
  phi ~ normal(0, 10);

  // Likelihood
  for (i in 1:N_model) {
    int t_idx = S + i;
    cases[t_idx] ~ neg_binomial_2(lambda[i], phi);
  }
}

generated quantities {
  vector[N_model] Rt = exp(log_Rt);

  array[N_model] int cases_pred;
  vector[N_model] log_lik;

  for (i in 1:N_model) {
    int t_idx = S + i;
    cases_pred[i] = neg_binomial_2_rng(lambda[i], phi);
    log_lik[i] = neg_binomial_2_lpmf(cases[t_idx] | lambda[i], phi);
  }

  // Effect sizes (multiplicative scale)
  real temp_effect = exp(beta[1]);
  real rain_effect = exp(beta[2]);
  real wolbachia_effect = exp(beta[3]);
  real npi_effect = exp(beta[4]);

  // Variance decomposition
  real var_climate = variance(f_climate);
  real var_wolbachia = variance(f_wolbachia);
  real var_npi = variance(f_npi);
  real var_residual = variance(f_residual);
  real var_total = var_climate + var_wolbachia + var_npi + var_residual;

  real prop_climate = var_climate / var_total;
  real prop_wolbachia = var_wolbachia / var_total;
  real prop_npi = var_npi / var_total;
  real prop_residual = var_residual / var_total;
}
