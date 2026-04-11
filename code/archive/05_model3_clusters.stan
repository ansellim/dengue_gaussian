// ==============================================================================
// Model 3: Full model + active dengue clusters as vector activity proxy
//
// log(Rt) = mu + beta_temp*temp + beta_rain*rain + beta_wolbachia*coverage
//           + beta_npi*stringency + beta_clusters*clusters + f_residual(t)
//
// Fit on 2015-2020 subperiod (where cluster data is available)
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

  int<lower=1> K;
  matrix[N_model, K] X;  // All covariates including clusters
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
  vector[K] beta;

  real log_alpha;
  real log_rho;

  vector[M] beta_gp_raw;
  real<lower=0> phi;
}

transformed parameters {
  real<lower=0> alpha = exp(log_alpha);
  real<lower=0> rho = exp(log_rho);

  vector[N_model] f_climate;
  vector[N_model] f_wolbachia;
  vector[N_model] f_npi;
  vector[N_model] f_clusters;
  vector[N_model] f_residual;
  vector[N_model] log_Rt;
  vector<lower=0>[N_model] lambda;

  f_climate = X[, 1] * beta[1] + X[, 2] * beta[2];
  f_wolbachia = X[, 3] * beta[3];
  f_npi = X[, 4] * beta[4];
  f_clusters = X[, 5] * beta[5];

  {
    vector[M] spd_weights;
    for (j in 1:M) {
      spd_weights[j] = sqrt(spd_matern32(sqrt_eigenvalues[j], alpha, rho));
    }
    f_residual = PHI * (spd_weights .* beta_gp_raw);
  }

  log_Rt = mu + f_climate + f_wolbachia + f_npi + f_clusters + f_residual;

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

  beta[1] ~ normal(0, 0.5);    // Temperature
  beta[2] ~ normal(0, 0.5);    // Rainfall
  beta[3] ~ normal(-0.3, 0.5); // Wolbachia
  beta[4] ~ normal(0, 0.5);    // NPI
  beta[5] ~ normal(0, 0.5);    // Clusters (standardized; no directional prior)

  log_alpha ~ normal(-1.6, 0.4);
  log_rho ~ normal(log(6), 0.5);

  beta_gp_raw ~ std_normal();
  phi ~ normal(0, 10);

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

  real temp_effect = exp(beta[1]);
  real rain_effect = exp(beta[2]);
  real wolbachia_effect = exp(beta[3]);
  real npi_effect = exp(beta[4]);
  real cluster_effect = exp(beta[5]);

  real var_climate = variance(f_climate);
  real var_wolbachia = variance(f_wolbachia);
  real var_npi = variance(f_npi);
  real var_clusters = variance(f_clusters);
  real var_residual = variance(f_residual);
  real var_total = var_climate + var_wolbachia + var_npi + var_clusters + var_residual;

  real prop_climate = var_climate / var_total;
  real prop_wolbachia = var_wolbachia / var_total;
  real prop_npi = var_npi / var_total;
  real prop_clusters = var_clusters / var_total;
  real prop_residual = var_residual / var_total;
}
