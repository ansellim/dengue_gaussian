// ==============================================================================
// Model 2 with shorter GP length scale
//
// Changes from standard model:
// - GP length scale prior: LogNormal(log(6), 0.5) instead of LogNormal(log(15), 0.5)
// - Median rho = 6 weeks (vs 15 weeks)
// - 95% CI: [2.2, 16.1] weeks
// - Allows GP to capture faster epidemic fluctuations
// - Should reduce residual autocorrelation
//
// Trade-off: More flexible GP may absorb some covariate signal
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

  // Intervention priors (same as standard model)
  beta[3] ~ normal(-0.3, 0.5);  // Wolbachia
  beta[4] ~ normal(0, 0.5);     // NPI

  // GP hyperparameters (log-scale priors)
  // ============================================================================
  // SHORTER LENGTH SCALE - key change for reducing residual autocorrelation
  // ============================================================================
  log_alpha ~ normal(-1.2, 0.5);  // alpha median ~0.3 (unchanged)
  log_rho ~ normal(log(6), 0.5);  // rho median 6 weeks, 95% CI ~[2.2, 16.1]
                                   // Changed from log(15) to allow faster fluctuations

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
