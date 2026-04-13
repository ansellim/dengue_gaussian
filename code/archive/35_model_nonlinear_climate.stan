// ==============================================================================
// Model: Nonlinear Climate-Rt Relationship via GP in Climate Space
//
// Replaces linear climate covariates with separate 1D Gaussian Processes
// that learn nonlinear mappings from temperature/rainfall to log(Rt).
//
// log(Rt) = mu + f_temp(temp[t]) + f_rain(rain[t]) + f_residual(t)
//
//   f_temp  ~ GP(0, k_SE(alpha_temp, rho_temp))   -- GP in temperature space
//   f_rain  ~ GP(0, k_SE(alpha_rain, rho_rain))   -- GP in rainfall space
//   f_residual ~ GP(0, k(alpha, rho))              -- GP in time (as before)
//
// Climate GPs use Squared Exponential kernel (smooth biological response).
// Temporal GP uses the data-supplied kernel_type for consistency with
// the kernel comparison framework.
//
// HSGP approximation:
//   Climate GPs: M_climate basis functions (e.g., 20), SE kernel always
//   Temporal GP: M basis functions (e.g., 200), selectable kernel
// ==============================================================================

functions {
  // Unified spectral density function for 1D stationary kernels
  real spd_kernel(real omega, real alpha, real rho, int kernel_type) {
    real result;

    if (kernel_type == 1) {
      // Matern 1/2 (Exponential)
      real inv_rho = 1.0 / rho;
      result = square(alpha) * 2.0 * inv_rho / (square(inv_rho) + square(omega));
    }
    else if (kernel_type == 2) {
      // Matern 3/2
      real sqrt3_over_rho = sqrt(3.0) / rho;
      result = square(alpha) * 4.0 * pow(sqrt3_over_rho, 3)
               / square(square(sqrt3_over_rho) + square(omega));
    }
    else if (kernel_type == 3) {
      // Matern 5/2
      real sqrt5_over_rho = sqrt(5.0) / rho;
      result = square(alpha) * (16.0 / 3.0) * pow(sqrt5_over_rho, 5)
               / pow(square(sqrt5_over_rho) + square(omega), 3);
    }
    else {
      // Squared Exponential (RBF)
      result = square(alpha) * sqrt(2.0 * pi()) * rho
               * exp(-0.5 * square(rho) * square(omega));
    }

    return result;
  }

  // SE-only spectral density (for climate GPs, always SE kernel)
  real spd_se(real omega, real alpha, real rho) {
    return square(alpha) * sqrt(2.0 * pi()) * rho
           * exp(-0.5 * square(rho) * square(omega));
  }
}

data {
  int<lower=1> N;            // Total number of weeks
  int<lower=1> N_model;      // Number of modeled weeks (N - S)
  int<lower=1> S;            // Maximum generation interval in weeks
  int<lower=1> M;            // Number of HSGP basis functions (temporal GP)

  array[N] int<lower=0> cases;
  vector[N_model] t;         // Centered time indices for temporal GP
  real<lower=0> L;           // Boundary for temporal HSGP
  simplex[S] gi;             // Generation interval PMF

  int<lower=1> K_climate;
  matrix[N_model, K_climate] X_climate;  // Standardized climate (for reference)

  int<lower=1, upper=4> kernel_type;  // Kernel for temporal GP

  // --- Nonlinear climate GP data ---
  int<lower=1> M_climate;           // Basis functions for climate GPs (e.g., 20)
  vector[N_model] temp_std;         // Standardized temperature (HSGP input)
  vector[N_model] rain_std;         // Standardized rainfall (HSGP input)
  real<lower=0> L_temp;             // Boundary for temperature HSGP
  real<lower=0> L_rain;             // Boundary for rainfall HSGP
}

transformed data {
  // --- Temporal GP basis functions (same as base model) ---
  matrix[N_model, M] PHI_time;
  vector[M] sqrt_eigenvalues_time;

  for (j in 1:M) {
    sqrt_eigenvalues_time[j] = (j * pi()) / (2 * L);
  }
  for (i in 1:N_model) {
    for (j in 1:M) {
      PHI_time[i, j] = sin(sqrt_eigenvalues_time[j] * (t[i] + L)) / sqrt(L);
    }
  }

  // --- Temperature GP basis functions (in climate space) ---
  matrix[N_model, M_climate] PHI_temp;
  vector[M_climate] sqrt_eigenvalues_temp;

  for (j in 1:M_climate) {
    sqrt_eigenvalues_temp[j] = (j * pi()) / (2 * L_temp);
  }
  for (i in 1:N_model) {
    for (j in 1:M_climate) {
      PHI_temp[i, j] = sin(sqrt_eigenvalues_temp[j] * (temp_std[i] + L_temp)) / sqrt(L_temp);
    }
  }

  // --- Rainfall GP basis functions (in climate space) ---
  matrix[N_model, M_climate] PHI_rain;
  vector[M_climate] sqrt_eigenvalues_rain;

  for (j in 1:M_climate) {
    sqrt_eigenvalues_rain[j] = (j * pi()) / (2 * L_rain);
  }
  for (i in 1:N_model) {
    for (j in 1:M_climate) {
      PHI_rain[i, j] = sin(sqrt_eigenvalues_rain[j] * (rain_std[i] + L_rain)) / sqrt(L_rain);
    }
  }
}

parameters {
  real mu;

  // --- Temperature GP hyperparameters ---
  real log_alpha_temp;
  real log_rho_temp;
  vector[M_climate] beta_temp_raw;

  // --- Rainfall GP hyperparameters ---
  real log_alpha_rain;
  real log_rho_rain;
  vector[M_climate] beta_rain_raw;

  // --- Temporal residual GP hyperparameters ---
  real log_alpha;
  real log_rho;
  vector[M] beta_gp_raw;

  // Overdispersion
  real<lower=0> phi;
}

transformed parameters {
  // Transform to natural scale
  real<lower=0> alpha_temp = exp(log_alpha_temp);
  real<lower=0> rho_temp = exp(log_rho_temp);
  real<lower=0> alpha_rain = exp(log_alpha_rain);
  real<lower=0> rho_rain = exp(log_rho_rain);
  real<lower=0> alpha = exp(log_alpha);
  real<lower=0> rho = exp(log_rho);

  vector[N_model] f_temp;
  vector[N_model] f_rain;
  vector[N_model] f_residual;
  vector[N_model] log_Rt;
  vector<lower=0>[N_model] lambda;

  // --- Temperature GP (SE kernel in temperature space) ---
  {
    vector[M_climate] spd_weights_temp;
    for (j in 1:M_climate) {
      spd_weights_temp[j] = sqrt(spd_se(sqrt_eigenvalues_temp[j], alpha_temp, rho_temp));
    }
    f_temp = PHI_temp * (spd_weights_temp .* beta_temp_raw);
  }

  // --- Rainfall GP (SE kernel in rainfall space) ---
  {
    vector[M_climate] spd_weights_rain;
    for (j in 1:M_climate) {
      spd_weights_rain[j] = sqrt(spd_se(sqrt_eigenvalues_rain[j], alpha_rain, rho_rain));
    }
    f_rain = PHI_rain * (spd_weights_rain .* beta_rain_raw);
  }

  // --- Temporal residual GP (selectable kernel) ---
  {
    vector[M] spd_weights_time;
    for (j in 1:M) {
      spd_weights_time[j] = sqrt(spd_kernel(sqrt_eigenvalues_time[j], alpha, rho, kernel_type));
    }
    f_residual = PHI_time * (spd_weights_time .* beta_gp_raw);
  }

  // log(Rt) = mu + f_temp + f_rain + f_residual
  log_Rt = mu + f_temp + f_rain + f_residual;

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
  // --- Priors ---
  mu ~ normal(0, 0.5);

  // Climate GP hyperparameters (small amplitude expected for climate effects)
  log_alpha_temp ~ normal(-2, 0.5);   // alpha_temp median ~0.14
  log_rho_temp ~ normal(0, 0.5);      // rho_temp median 1 (in standardized units)
  log_alpha_rain ~ normal(-2, 0.5);
  log_rho_rain ~ normal(0, 0.5);

  // Climate GP basis coefficients (non-centered)
  beta_temp_raw ~ std_normal();
  beta_rain_raw ~ std_normal();

  // Temporal GP hyperparameters (same as base model)
  log_alpha ~ normal(-1.2, 0.5);      // alpha median ~0.3
  log_rho ~ normal(log(6), 0.5);      // rho median 6 weeks
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

  // --- Variance decomposition (3 components) ---
  real var_temp = variance(f_temp);
  real var_rain = variance(f_rain);
  real var_residual = variance(f_residual);
  real var_total = var_temp + var_rain + var_residual;

  real prop_temp = var_temp / var_total;
  real prop_rain = var_rain / var_total;
  real prop_residual = var_residual / var_total;

  // Combined climate proportion for comparison with linear model
  real var_climate = var_temp + var_rain;
  real prop_climate = var_climate / var_total;
}
