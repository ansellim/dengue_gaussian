// ==============================================================================
// Climate + susceptible-fraction covariate(s) + GP residual.
//
// Methodological audit pipeline (M0 / M1 / M2). Sister model to
// 05_model_climate_only.stan. The GP residual is *kept* — this model adds an
// explicit susceptible covariate alongside it, not in place of it. The audit
// asks: does adding explicit susceptibility absorb GP variance / improve LOO?
//
// log(Rt) = mu + X_climate * beta_climate + X_S * beta_S + f_residual
//   f_residual ~ HSGP(Matern 3/2, alpha, rho)
//
// X_S has K_S columns of standardized log(S(t)). Used as:
//   K_S = 1, X_S = log(S_pop)        => M1 (single national pool)
//   K_S = 1, X_S = log(S_dominant)   => M2 (currently dominant pool)
//
// Setting K_S = 0 is not supported here on purpose; the climate-only baseline
// (M0) lives in 05_model_climate_only.stan.
//
// Caveat: log(S(t)) is constructed from cumulative cases, which is also the
// model outcome. A positive ELPD improvement with this covariate has built-in
// feedback that can absorb low-frequency drift trivially. Interpret accordingly.
// ==============================================================================

functions {
  real spd_matern32(real omega, real alpha, real rho) {
    real sqrt3_over_rho = sqrt(3.0) / rho;
    real numerator = square(alpha) * 4.0 * pow(sqrt3_over_rho, 3);
    real denominator = square(square(sqrt3_over_rho) + square(omega));
    return numerator / denominator;
  }
}

data {
  int<lower=1> N;
  int<lower=1> N_model;
  int<lower=1> S;
  int<lower=1> M;

  array[N] int<lower=0> cases;
  vector[N_model] t;
  real<lower=0> L;
  simplex[S] gi;

  int<lower=1> K_climate;
  matrix[N_model, K_climate] X_climate;

  int<lower=1> K_S;
  matrix[N_model, K_S] X_S;

  real log_alpha_mu;
  real log_rho_mu;
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
  vector[K_climate] beta_climate;
  vector[K_S] beta_S;

  real log_alpha;
  real log_rho;

  vector[M] beta_gp_raw;
  real<lower=0> phi;
}

transformed parameters {
  real<lower=0> alpha = exp(log_alpha);
  real<lower=0> rho = exp(log_rho);

  vector[N_model] f_climate     = X_climate * beta_climate;
  vector[N_model] f_susceptible = X_S * beta_S;
  vector[N_model] f_residual;
  vector[N_model] log_Rt;
  vector<lower=0>[N_model] lambda;

  {
    vector[M] spd_weights;
    for (j in 1:M) {
      spd_weights[j] = sqrt(spd_matern32(sqrt_eigenvalues[j], alpha, rho));
    }
    f_residual = PHI * (spd_weights .* beta_gp_raw);
  }

  log_Rt = mu + f_climate + f_susceptible + f_residual;

  for (i in 1:N_model) {
    real infectious_pressure = 0;
    int t_idx = S + i;

    for (s in 1:S) {
      infectious_pressure += cases[t_idx - s] * gi[s];
    }

    lambda[i] = fmax(exp(log_Rt[i]) * infectious_pressure, 1.0);
  }
}

model {
  mu ~ normal(0, 0.5);

  beta_climate ~ normal(0, 0.5);

  // Neutral prior on the susceptible-covariate slope. SIR theory predicts
  // beta_S ~ +1 on raw log(S), but log(S) here is standardized (sd ~ 0.05),
  // so a unit-1 raw effect becomes a tiny standardized effect. Normal(0, 1)
  // is wide enough that the posterior is data-driven.
  beta_S ~ normal(0, 1);

  log_alpha ~ normal(log_alpha_mu, 0.5);
  log_rho   ~ normal(log_rho_mu,   0.5);

  beta_gp_raw ~ std_normal();

  phi ~ normal(0, 5);

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

  // Variance decomposition (climate + susceptible + residual)
  real var_climate     = variance(f_climate);
  real var_susceptible = variance(f_susceptible);
  real var_residual    = variance(f_residual);
  real var_total       = var_climate + var_susceptible + var_residual;

  real prop_climate     = var_climate     / var_total;
  real prop_susceptible = var_susceptible / var_total;
  real prop_residual    = var_residual    / var_total;
}
