data {
  int<lower=1> N;              // number of observations
  array[N] real x;             // univariate covariate (e.g. time)
  vector[N] y;                 // target variable
  int<lower=1> J;              // number of subjects/groups
  array[N] int<lower=1, upper=J> subject; // subject index for each obs
  vector[N] z;                 // continuous scaling covariate per obs
                               // (e.g. disease age / exposure; 0 = baseline)
}

transformed data {
  real delta = 1e-9;
}

parameters {
  real intercept;

  // Shared (population-level) GP
  real<lower=0> length_scale_k1;
  real<lower=0> variance_k1;

  // Individual-level GP (lgpr-style: scaled by z)
  real<lower=0> length_scale_k2;
  real<lower=0> variance_k2;

  real<lower=0> sigma;

  vector[N] eta;       // non-centered parameterisation for shared GP
  vector[N] eta_ind;   // non-centered parameterisation for individual GP
}

transformed parameters {
  vector[N] f;
  vector[N] f_1;    // shared GP component
  vector[N] f_2;    // individual (scaled) GP component

  {
    matrix[N, N] K_1;
    matrix[N, N] L_K_1;
    matrix[N, N] K_2;
    matrix[N, N] L_K_2;

    // --- Shared GP kernel (OU, unchanged) ---
    for (i in 1:N) {
      for (j in i:N) {
        real d = fabs(x[i] - x[j]);
        K_1[i, j] = square(variance_k1) * exp(-d / length_scale_k1);
        K_1[j, i] = K_1[i, j];
      }
    }
    for (n in 1:N) K_1[n, n] = K_1[n, n] + delta;

    L_K_1 = cholesky_decompose(K_1);
    f_1 = L_K_1 * eta;

    // --- Individual GP kernel (lgpr approach) ---
    // The kernel is masked to zero between different subjects,
    // and scaled by z on both sides: k2(i,j) = z_i * OU(x_i,x_j) * z_j
    // This is the lgpr "heterogeneous" component.
    for (i in 1:N) {
      for (j in i:N) {
        if (subject[i] != subject[j]) {
          // Different subjects: kernel is exactly zero (no cross-subject covariance)
          K_2[i, j] = 0.0;
          K_2[j, i] = 0.0;
        } else {
          real d = fabs(x[i] - x[j]);
          // OU base kernel, scaled by z_i and z_j (lgpr heterogeneous scaling)
          K_2[i, j] = z[i] * square(variance_k2) * exp(-d / length_scale_k2) * z[j];
          K_2[j, i] = K_2[i, j];
        }
      }
    }
    for (n in 1:N) K_2[n, n] = K_2[n, n] + delta;

    L_K_2 = cholesky_decompose(K_2);
    f_2 = L_K_2 * eta_ind;

    f = f_1 + f_2;
  }
}

model {
  // Priors (unchanged for shared GP)
  intercept        ~ normal(0, 5);
  length_scale_k1  ~ inv_gamma(5, 5);
  variance_k1      ~ inv_gamma(5, 5);

  // Priors for individual GP
  length_scale_k2  ~ inv_gamma(5, 5);
  variance_k2      ~ inv_gamma(5, 5);

  sigma     ~ std_normal();
  eta       ~ std_normal();
  eta_ind   ~ std_normal();

  // Likelihood
  y ~ normal(intercept + f, sigma);
}

generated quantities {
  vector[N] fit   = intercept + f;
  vector[N] eff_1 = f_1;   // shared GP effect
  vector[N] eff_2 = f_2;   // individual (scaled) GP effect
}
