// ==============================================================================
// Model E: Mechanistic mean function (immunity dynamics) with GP residual
//
// Key difference from Model 3 (climate-only):
//   Model 3: log(Rt) = mu + f_climate + f_residual  (GP has constant mean)
//   Model E: log(Rt) = log(R0 * S_eff) + f_climate + f_residual
//                     = log_R0 + log_S_eff + f_climate + f_residual
//
// The GP residual now captures deviations FROM the immunity prediction.
// If the GP amplitude (alpha) shrinks dramatically -> immunity explains
// most Rt variation that the GP was previously absorbing.
//
// log_S_eff is pre-computed from cumulative case depletion with serotype-
// specific resets (expansion factor ~8, N=5.5M, initial_S=0.75).
//
// GP: HSGP approximation with Matern 3/2 spectral density.
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

  int<lower=1> K_climate;
  matrix[N_model, K_climate] X_climate;

  // Mechanistic mean: pre-computed log(effective susceptible fraction)
  vector[N_model] log_S_eff;
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
  real log_R0;                // log of basic reproduction number
  vector[K_climate] beta_climate;

  // Log-parameterized GP hyperparameters (better geometry)
  real log_alpha;
  real log_rho;

  vector[M] beta_gp_raw;     // Raw (non-centered) GP coefficients
  real<lower=0> phi;
}

transformed parameters {
  // Transform to natural scale
  real<lower=0> alpha = exp(log_alpha);
  real<lower=0> rho = exp(log_rho);

  vector[N_model] mean_function;
  vector[N_model] f_climate;
  vector[N_model] f_residual;
  vector[N_model] log_Rt;
  vector<lower=0>[N_model] lambda;

  // Mechanistic mean function: log(R0 * S_eff) = log_R0 + log_S_eff
  mean_function = log_R0 + log_S_eff;

  // Climate contribution (small adjustment)
  f_climate = X_climate[, 1] * beta_climate[1] + X_climate[, 2] * beta_climate[2];

  // Residual GP with non-centered parameterization
  {
    vector[M] spd_weights;
    for (j in 1:M) {
      spd_weights[j] = sqrt(spd_matern32(sqrt_eigenvalues[j], alpha, rho));
    }
    f_residual = PHI * (spd_weights .* beta_gp_raw);
  }

  // Total: mechanistic mean + climate + GP residual
  log_Rt = mean_function + f_climate + f_residual;

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
  // Priors
  log_R0 ~ normal(log(1.3), 0.3);  // Tan et al. (2019): R0 ~ 1.28-1.30

  // Climate effects (regularized, no strong directional prior)
  beta_climate[1] ~ normal(0, 0.5);  // Temperature
  beta_climate[2] ~ normal(0, 0.5);  // Rainfall

  // GP hyperparameters (log-scale priors, same as Model 3)
  log_alpha ~ normal(-1.2, 0.5);  // alpha median ~0.3
  log_rho ~ normal(log(6), 0.5);  // rho median 6 weeks, 95% CI ~[2.2, 16.1]

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

  // R0 on natural scale
  real R0 = exp(log_R0);

  // Effect sizes (multiplicative scale)
  real temp_effect = exp(beta_climate[1]);
  real rain_effect = exp(beta_climate[2]);

  // Variance decomposition: mechanistic + climate + residual
  real mean_mean_function = mean(mean_function);
  real var_mechanistic = variance(mean_function - mean_mean_function);
  real var_climate = variance(f_climate);
  real var_residual = variance(f_residual);
  real var_total = var_mechanistic + var_climate + var_residual;

  real prop_mechanistic = var_mechanistic / var_total;
  real prop_climate = var_climate / var_total;
  real prop_residual = var_residual / var_total;
}
