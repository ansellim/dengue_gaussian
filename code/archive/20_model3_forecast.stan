// ==============================================================================
// Model 3 Forecast: Climate + optional entropy with selectable GP kernel
//
// Extension of 15_model3_multikernel.stan for out-of-sample forecasting.
//
// log(Rt) = mu + f_climate + f_residual + beta_entropy * entropy (if enabled)
//   f_climate  = beta_temp * temp + beta_rain * rain
//   f_residual ~ HSGP(kernel, alpha, rho)
//
// Key differences from 15_model3_multikernel.stan:
//   - Likelihood restricted to training period (1:N_train)
//   - Optional entropy covariate (use_entropy flag)
//   - Forward simulation in generated quantities for forecast period
//   - log_lik computed for all N_model time points
// ==============================================================================

functions {
  // Unified spectral density function for 1D stationary kernels
  // Returns S(omega | alpha, rho) for the selected kernel
  real spd_kernel(real omega, real alpha, real rho, int kernel_type) {
    real result;

    if (kernel_type == 1) {
      // ---- Matern 1/2 (Exponential) ----
      real inv_rho = 1.0 / rho;
      result = square(alpha) * 2.0 * inv_rho / (square(inv_rho) + square(omega));
    }
    else if (kernel_type == 2) {
      // ---- Matern 3/2 ----
      real sqrt3_over_rho = sqrt(3.0) / rho;
      result = square(alpha) * 4.0 * pow(sqrt3_over_rho, 3)
               / square(square(sqrt3_over_rho) + square(omega));
    }
    else if (kernel_type == 3) {
      // ---- Matern 5/2 ----
      real sqrt5_over_rho = sqrt(5.0) / rho;
      result = square(alpha) * (16.0 / 3.0) * pow(sqrt5_over_rho, 5)
               / pow(square(sqrt5_over_rho) + square(omega), 3);
    }
    else {
      // ---- Squared Exponential (RBF) ----
      result = square(alpha) * sqrt(2.0 * pi()) * rho
               * exp(-0.5 * square(rho) * square(omega));
    }

    return result;
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

  int<lower=1, upper=4> kernel_type;  // 1=M12, 2=M32, 3=M52, 4=SE

  // Forecast-specific data
  int<lower=1> N_train;               // Number of training time points (<= N_model)
  int<lower=0> N_forecast;            // N_model - N_train
  int<lower=0, upper=1> use_entropy;  // Whether to include entropy covariate
  vector[N_model] entropy_covariate;  // Lagged, standardized Shannon entropy
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

  // Entropy effect
  real beta_entropy;
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
      spd_weights[j] = sqrt(spd_kernel(sqrt_eigenvalues[j], alpha, rho, kernel_type));
    }
    f_residual = PHI * (spd_weights .* beta_gp_raw);
  }

  // log(Rt) = mu + climate + residual + entropy (if enabled)
  log_Rt = mu + f_climate + f_residual;
  if (use_entropy) {
    log_Rt += beta_entropy * entropy_covariate;
  }

  // Renewal equation (using observed cases for infectious pressure)
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
  // Priors — identical across all kernels for fair comparison
  mu ~ normal(0, 0.5);

  // Climate effects
  beta_climate[1] ~ normal(0, 0.5);  // Temperature
  beta_climate[2] ~ normal(0, 0.5);  // Rainfall

  // GP hyperparameters (log-scale priors)
  log_alpha ~ normal(-1.2, 0.5);  // alpha median ~0.3
  log_rho ~ normal(log(6), 0.5);  // rho median 6 weeks

  // GP basis coefficients (non-centered)
  beta_gp_raw ~ std_normal();

  // Overdispersion
  phi ~ normal(0, 10);

  // Entropy effect
  beta_entropy ~ normal(0, 0.5);

  // Likelihood — training period only
  for (i in 1:N_train) {
    int t_idx = S + i;
    cases[t_idx] ~ neg_binomial_2(lambda[i], phi);
  }
}

generated quantities {
  vector[N_model] Rt = exp(log_Rt);

  array[N_model] int cases_pred;
  vector[N_model] log_lik;

  // Forward simulation: use observed cases in training, simulated in forecast
  array[N_model] int cases_forward;

  // --- Training period: use observed cases ---
  for (i in 1:N_train) {
    int t_idx = S + i;
    cases_pred[i] = neg_binomial_2_rng(lambda[i], phi);
    log_lik[i] = neg_binomial_2_lpmf(cases[t_idx] | lambda[i], phi);
    cases_forward[i] = cases[t_idx];  // observed
  }

  // --- Forecast period: sequential forward simulation ---
  for (i in (N_train + 1):N_model) {
    real ip_fwd = 0;

    for (s in 1:S) {
      int lb = S + i - s;  // absolute index into cases[] or cases_forward[]

      if (lb <= S) {
        // Within burn-in period: use original observed cases
        ip_fwd += cases[lb] * gi[s];
      } else if (lb - S <= N_train) {
        // Training period: use observed cases
        ip_fwd += cases[lb] * gi[s];
      } else {
        // Forecast period: use forward-simulated cases
        ip_fwd += cases_forward[lb - S] * gi[s];
      }
    }

    real lambda_fwd = fmax(exp(log_Rt[i]) * ip_fwd, 1.0);
    cases_forward[i] = neg_binomial_2_rng(lambda_fwd, phi);
    cases_pred[i] = cases_forward[i];

    // Log-likelihood for forecast period (using forward-simulated lambda)
    {
      int t_idx = S + i;
      log_lik[i] = neg_binomial_2_lpmf(cases[t_idx] | lambda_fwd, phi);
    }
  }

  // Effect sizes (multiplicative scale)
  real temp_effect = exp(beta_climate[1]);
  real rain_effect = exp(beta_climate[2]);
  real entropy_effect = exp(beta_entropy);

  // Variance decomposition
  real var_climate = variance(f_climate);
  real var_residual = variance(f_residual);
  real var_total = var_climate + var_residual;

  real prop_climate = var_climate / var_total;
  real prop_residual = var_residual / var_total;
}
