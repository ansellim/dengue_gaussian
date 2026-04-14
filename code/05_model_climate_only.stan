// ==============================================================================
// Climate-only model: temperature and rainfall covariates with short GP length scale
//
// Focuses on climate drivers and residual dynamics (serotype switching analysis).
//
// log(Rt) = mu + f_climate + f_residual
//   f_climate  = beta_temp * temp + beta_rain * rain
//   f_residual ~ HSGP(Matern 3/2, alpha, rho)
//
// GP hyperprior means are passed in as data (log_alpha_mu, log_rho_mu) so the
// same compiled model can serve both the baseline fit and the prior sweep in
// 15_gp_prior_sensitivity.R. Baseline values: log_alpha_mu = -1.2 (alpha
// median = 0.30), log_rho_mu = log(6) (rho median = 6 weeks). Sigma is fixed
// at 0.5 on the log scale.
//
// Assumptions:
//   - 100% case ascertainment (Singapore mandatory notification)
//   - Constant reporting delay (~5-7 days) absorbed at weekly resolution
//     (standard for weekly dengue Rt; Cori et al. 2013, Lau et al. 2022)
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

  // GP hyperprior means (log scale). See header for baseline values.
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
  vector[N_model] f_residual;
  vector[N_model] log_Rt;
  vector<lower=0>[N_model] lambda;

  // Covariate effects (climate only)
  f_climate = X_climate[, 1] * beta_climate[1] + X_climate[, 2] * beta_climate[2];

  // Residual GP with non-centered parameterization
  {
    vector[M] spd_weights;
    for (j in 1:M) {
      spd_weights[j] = sqrt(spd_matern32(sqrt_eigenvalues[j], alpha, rho));
    }
    f_residual = PHI * (spd_weights .* beta_gp_raw);
  }

  // log(Rt)
  log_Rt = mu + f_climate + f_residual;

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
  mu ~ normal(0, 0.5);

  // Climate effects (regularized, no strong directional prior)
  beta_climate[1] ~ normal(0, 0.5);  // Temperature
  beta_climate[2] ~ normal(0, 0.5);  // Rainfall

  // GP hyperparameters (log-scale priors; means passed in as data)
  log_alpha ~ normal(log_alpha_mu, 0.5);
  log_rho   ~ normal(log_rho_mu,   0.5);

  // GP basis coefficients (non-centered)
  beta_gp_raw ~ std_normal();

  // Overdispersion — tightened from Normal(0, 10)+ to push phi toward the
  // biologically sensible 5–15 range for weekly dengue counts. Previous prior
  // let the data pull phi to ~47 (near-Poisson), producing PPC intervals too
  // narrow to cover weekly fluctuations.
  phi ~ normal(0, 5);

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

  // Variance decomposition (climate + residual only)
  real var_climate = variance(f_climate);
  real var_residual = variance(f_residual);
  real var_total = var_climate + var_residual;

  real prop_climate = var_climate / var_total;
  real prop_residual = var_residual / var_total;
}
