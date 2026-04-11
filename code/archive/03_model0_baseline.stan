// ==============================================================================
// Model 0: Baseline - Single GP on log(Rt)
//
// log(Rt) = mu + f(t)
// where f(t) ~ GP(0, k_matern32(t, t'))
//
// Renewal equation: lambda_t = R_t * sum_{s=1}^{S} c_{t-s} * w_s
// Observation model: c_t ~ NegBin(lambda_t, phi)
//
// GP approximated via Hilbert Space basis functions (HSGP)
// ==============================================================================

functions {
  // Spectral density of Matern 3/2 kernel in 1D
  // S(omega) = alpha^2 * 4 * (sqrt(3)/rho)^3 / ((sqrt(3)/rho)^2 + omega^2)^2
  real spd_matern32(real omega, real alpha, real rho) {
    real sqrt3_over_rho = sqrt(3.0) / rho;
    real numerator = alpha^2 * 4.0 * sqrt3_over_rho^3;
    real denominator = (sqrt3_over_rho^2 + omega^2)^2;
    return numerator / denominator;
  }

  // HSGP basis functions: phi_j(t) = sin(j * pi * (t + L) / (2L)) / sqrt(L)
  // Input t is assumed to be in [-1, 1] (scaled)
  // L is the boundary factor
  vector hsgp_basis(vector t, int M, real L) {
    int N = num_elements(t);
    matrix[N, M] PHI;

    for (j in 1:M) {
      real lambda_j = (j * pi()) / (2 * L);
      for (i in 1:N) {
        PHI[i, j] = sin(lambda_j * (t[i] + L)) / sqrt(L);
      }
    }

    return PHI * rep_vector(1.0, M);  // Placeholder - actual computation in model
  }

  // Compute eigenvalues for HSGP
  vector hsgp_eigenvalues(int M, real L) {
    vector[M] lambda;
    for (j in 1:M) {
      lambda[j] = (j * pi() / (2 * L))^2;
    }
    return lambda;
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

  // Generation interval PMF (sums to 1)
  simplex[S] gi;
}

transformed data {
  // Precompute HSGP basis matrix
  matrix[N_model, M] PHI;
  vector[M] sqrt_eigenvalues;

  // Eigenvalues: lambda_j = (j * pi / (2L))^2
  for (j in 1:M) {
    sqrt_eigenvalues[j] = (j * pi()) / (2 * L);
  }

  // Basis functions: phi_j(t) = sin(sqrt(lambda_j) * (t + L)) / sqrt(L)
  for (i in 1:N_model) {
    for (j in 1:M) {
      PHI[i, j] = sin(sqrt_eigenvalues[j] * (t[i] + L)) / sqrt(L);
    }
  }
}

parameters {
  // Intercept for log(Rt)
  real mu;

  // GP hyperparameters (log-parameterized for better HMC geometry)
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

  // GP function values at modeled time points
  vector[N_model] f_gp;

  // Log reproduction number
  vector[N_model] log_Rt;

  // Expected cases (lambda)
  vector[N_model] lambda;

  // Compute spectral weights for HSGP
  {
    vector[M] spd_weights;
    for (j in 1:M) {
      real omega_j = sqrt_eigenvalues[j];
      spd_weights[j] = sqrt(spd_matern32(omega_j, alpha, rho));
    }

    // GP values: f = PHI * (spd_weights .* beta_gp)
    f_gp = PHI * (spd_weights .* beta_gp);
  }

  // log(Rt) = mu + f(t)
  log_Rt = mu + f_gp;

  // Renewal equation: lambda_t = R_t * sum_{s=1}^{S} c_{t-s} * w_s
  for (i in 1:N_model) {
    real infectious_pressure = 0;
    int t_idx = S + i;  // Index in full cases array

    for (s in 1:S) {
      infectious_pressure += cases[t_idx - s] * gi[s];
    }

    lambda[i] = exp(log_Rt[i]) * infectious_pressure;

    // Ensure lambda is positive (numerical stability)
    if (lambda[i] < 1e-6) {
      lambda[i] = 1e-6;
    }
  }
}

model {
  // ---- Priors ----

  // Intercept: weakly informative, centered at Rt = 1
  mu ~ normal(0, 0.5);

  // GP hyperparameters (log-scale priors, rho in weeks)
  log_alpha ~ normal(-1.2, 0.5);  // alpha median ~0.3
  log_rho ~ normal(log(10), 0.5); // rho median 10 weeks, 95% ~[3.8, 26.6]

  // HSGP basis coefficients: standard normal (scaling handled by spd_weights)
  beta_gp ~ std_normal();

  // Overdispersion: half-Normal(0, 10)
  phi ~ normal(0, 10);

  // ---- Likelihood ----

  // Negative binomial: cases ~ NegBin(lambda, phi)
  // Parameterization: Var = mu + mu^2/phi
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
}
