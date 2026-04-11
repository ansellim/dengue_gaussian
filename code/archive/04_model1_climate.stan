// ==============================================================================
// Model 1: Climate Covariates + Residual GP
//
// log(Rt) = mu + beta_temp * temp + beta_rain * rain + f_residual(t)
// where f_residual(t) ~ GP(0, k_matern32(t, t'))
//
// Renewal equation: lambda_t = R_t * sum_{s=1}^{S} c_{t-s} * w_s
// Observation model: c_t ~ NegBin(lambda_t, phi)
//
// GP approximated via Hilbert Space basis functions (HSGP)
// ==============================================================================

functions {
  // Spectral density of Matern 3/2 kernel in 1D
  real spd_matern32(real omega, real alpha, real rho) {
    real sqrt3_over_rho = sqrt(3.0) / rho;
    real numerator = alpha^2 * 4.0 * sqrt3_over_rho^3;
    real denominator = (sqrt3_over_rho^2 + omega^2)^2;
    return numerator / denominator;
  }
}

data {
  // Dimensions
  int<lower=1> N;            // Total number of weeks
  int<lower=1> N_model;      // Number of modeled weeks (N - S)
  int<lower=1> S;            // Maximum generation interval in weeks
  int<lower=1> M;            // Number of HSGP basis functions

  // Observations
  array[N] int<lower=0> cases;  // Weekly case counts

  // Time (scaled to approximately [-1, 1])
  vector[N_model] t;
  real<lower=0> L;           // HSGP boundary

  // Generation interval PMF
  simplex[S] gi;

  // Climate covariates (standardized)
  int<lower=1> K_climate;    // Number of climate covariates (2: temp, rain)
  matrix[N_model, K_climate] X_climate;
}

transformed data {
  // Precompute HSGP basis matrix
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
  // Intercept for log(Rt)
  real mu;

  // Climate regression coefficients
  vector[K_climate] beta_climate;

  // Residual GP hyperparameters (log-parameterized for better HMC geometry)
  real log_alpha;
  real log_rho;

  // HSGP basis coefficients
  vector[M] beta_gp;

  // Observation model
  real<lower=0> phi;         // Overdispersion parameter
}

transformed parameters {
  // Transform to natural scale
  real<lower=0> alpha = exp(log_alpha);
  real<lower=0> rho = exp(log_rho);

  // Climate effect
  vector[N_model] f_climate;

  // Residual GP function values
  vector[N_model] f_residual;

  // Log reproduction number
  vector[N_model] log_Rt;

  // Expected cases (lambda)
  vector[N_model] lambda;

  // Climate effect: X * beta
  f_climate = X_climate * beta_climate;

  // Compute residual GP
  {
    vector[M] spd_weights;
    for (j in 1:M) {
      real omega_j = sqrt_eigenvalues[j];
      spd_weights[j] = sqrt(spd_matern32(omega_j, alpha, rho));
    }
    f_residual = PHI * (spd_weights .* beta_gp);
  }

  // log(Rt) = mu + climate effect + residual GP
  log_Rt = mu + f_climate + f_residual;

  // Renewal equation
  for (i in 1:N_model) {
    real infectious_pressure = 0;
    int t_idx = S + i;

    for (s in 1:S) {
      infectious_pressure += cases[t_idx - s] * gi[s];
    }

    lambda[i] = exp(log_Rt[i]) * infectious_pressure;

    if (lambda[i] < 1e-6) {
      lambda[i] = 1e-6;
    }
  }
}

model {
  // ---- Priors ----

  // Intercept
  mu ~ normal(0, 0.5);

  // Climate coefficients: weakly informative
  beta_climate ~ normal(0, 0.5);

  // Residual GP hyperparameters (log-scale priors, rho in weeks)
  log_alpha ~ normal(-1.2, 0.5);  // alpha median ~0.3
  log_rho ~ normal(log(10), 0.5); // rho median 10 weeks, 95% ~[3.8, 26.6]

  // HSGP basis coefficients
  beta_gp ~ std_normal();

  // Overdispersion
  phi ~ normal(0, 10);

  // ---- Likelihood ----
  for (i in 1:N_model) {
    int t_idx = S + i;
    cases[t_idx] ~ neg_binomial_2(lambda[i], phi);
  }
}

generated quantities {
  // Rt values
  vector[N_model] Rt = exp(log_Rt);

  // Posterior predictive cases
  array[N_model] int cases_pred;
  for (i in 1:N_model) {
    cases_pred[i] = neg_binomial_2_rng(lambda[i], phi);
  }

  // Log-likelihood for model comparison
  vector[N_model] log_lik;
  for (i in 1:N_model) {
    int t_idx = S + i;
    log_lik[i] = neg_binomial_2_lpmf(cases[t_idx] | lambda[i], phi);
  }

  // Decomposition: proportion of variance explained by climate
  real var_climate = variance(f_climate);
  real var_residual = variance(f_residual);
  real var_total = variance(log_Rt - mu);
  real prop_climate = var_climate / (var_climate + var_residual);

  // Residual GP amplitude (empirical)
  real residual_amplitude = sd(f_residual);
}
