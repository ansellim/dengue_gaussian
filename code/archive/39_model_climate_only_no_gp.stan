// ==============================================================================
// Model: Climate-only, NO GP residual
//
// Ablation model to test whether the GP absorbs climate signal.
// log(Rt) = mu + beta_temp * temp + beta_rain * rain
//
// No GP, no basis functions, no HSGP. Observation model is still NegBin
// with the renewal equation.
//
// If climate effects (beta_temp, beta_rain) remain negligible even without
// the GP competing, that confirms climate is genuinely uninformative.
// ==============================================================================

data {
  int<lower=1> N;            // Total number of weeks
  int<lower=1> N_model;      // Number of modeled weeks (N - S)
  int<lower=1> S;            // Maximum generation interval in weeks

  array[N] int<lower=0> cases;
  simplex[S] gi;

  int<lower=1> K_climate;
  matrix[N_model, K_climate] X_climate;
}

parameters {
  real mu;
  vector[K_climate] beta_climate;
  real<lower=0> phi;
}

transformed parameters {
  vector[N_model] log_Rt;
  vector<lower=0>[N_model] lambda;

  // log(Rt) = intercept + climate effects only (no GP)
  log_Rt = mu + X_climate * beta_climate;

  // Renewal equation
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
  // Priors (same as Model 3)
  mu ~ normal(0, 0.5);
  beta_climate[1] ~ normal(0, 0.5);  // Temperature
  beta_climate[2] ~ normal(0, 0.5);  // Rainfall
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
  real temp_effect = exp(beta_climate[1]);
  real rain_effect = exp(beta_climate[2]);
}
