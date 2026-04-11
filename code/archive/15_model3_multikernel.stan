// ==============================================================================
// Model 3 Multi-Kernel: Climate-only covariates with selectable GP kernel
//
// Identical to 05_model3_climate_only.stan but with a kernel_type data variable
// that selects the spectral density function for the HSGP approximation.
//
// log(Rt) = mu + f_climate + f_residual
//   f_climate  = beta_temp * temp + beta_rain * rain
//   f_residual ~ HSGP(kernel, alpha, rho)
//
// Kernel types:
//   1 = Matern 1/2 (Exponential / Ornstein-Uhlenbeck) — rough, non-differentiable
//   2 = Matern 3/2 — once-differentiable (current default)
//   3 = Matern 5/2 — twice-differentiable
//   4 = Squared Exponential (RBF) — infinitely smooth
//
// Ascertainment note:
//   Under constant case ascertainment rho, the renewal equation on cases yields
//   the same Rt as on infections (rho cancels; Cori et al. 2013). We therefore
//   work directly with reported cases. See Discussion for caveats regarding
//   time-varying rho (Tan et al. 2019: 1:14 in 2005-2009 to 1:6 in 2014-2017).
//
// All other model components (priors, renewal equation, observation model)
// are identical across kernels to ensure a fair comparison.
// ==============================================================================

functions {
  // Unified spectral density function for 1D stationary kernels
  // Returns S(omega | alpha, rho) for the selected kernel
  real spd_kernel(real omega, real alpha, real rho, int kernel_type) {
    real result;

    if (kernel_type == 1) {
      // ---- Matern 1/2 (Exponential) ----
      // S(w) = alpha^2 * 2 * (1/rho) / ((1/rho)^2 + w^2)
      real inv_rho = 1.0 / rho;
      result = square(alpha) * 2.0 * inv_rho / (square(inv_rho) + square(omega));
    }
    else if (kernel_type == 2) {
      // ---- Matern 3/2 ----
      // S(w) = alpha^2 * 4 * (sqrt(3)/rho)^3 / ((sqrt(3)/rho)^2 + w^2)^2
      real sqrt3_over_rho = sqrt(3.0) / rho;
      result = square(alpha) * 4.0 * pow(sqrt3_over_rho, 3)
               / square(square(sqrt3_over_rho) + square(omega));
    }
    else if (kernel_type == 3) {
      // ---- Matern 5/2 ----
      // S(w) = alpha^2 * (16/3) * (sqrt(5)/rho)^5 / ((sqrt(5)/rho)^2 + w^2)^3
      real sqrt5_over_rho = sqrt(5.0) / rho;
      result = square(alpha) * (16.0 / 3.0) * pow(sqrt5_over_rho, 5)
               / pow(square(sqrt5_over_rho) + square(omega), 3);
    }
    else {
      // ---- Squared Exponential (RBF) ----
      // S(w) = alpha^2 * sqrt(2*pi) * rho * exp(-rho^2 * w^2 / 2)
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
      spd_weights[j] = sqrt(spd_kernel(sqrt_eigenvalues[j], alpha, rho, kernel_type));
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
  // Priors — identical across all kernels for fair comparison
  mu ~ normal(0, 0.5);

  // Climate effects (regularized, no strong directional prior)
  beta_climate[1] ~ normal(0, 0.5);  // Temperature
  beta_climate[2] ~ normal(0, 0.5);  // Rainfall

  // GP hyperparameters (log-scale priors)
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
