data {
  int<lower=1> N;
  array[N] real x;
  vector[N] y;
  int<lower=1> J;
  array[N] int<lower=1, upper=J> subject;
  vector[N] z;
}

transformed data {
  real delta = 1e-9;
  real loc1 = -24.0;   // first changepoint:  anticipatory -> central
  real loc2 =  24.0;   // second changepoint: central -> recovery
}

parameters {
  real intercept_s1;   // anticipatory baseline
  real intercept_s2;   // central (around event)
  real intercept_s3;   // recovery baseline

  // Shared GP — EQ kernels for each of the three sections
  real<lower=0> length_scale_k1_s1;   // anticipatory
  real<lower=0> length_scale_k1_s2;   // central
  real<lower=0> length_scale_k1_s3;   // recovery
  real<lower=0> variance_k1_s1;
  real<lower=0> variance_k1_s2;       // will get higher-variance prior
  real<lower=0> variance_k1_s3;
  real<lower=0> steepness1;           // steeper  (anticipatory -> central)
  real<lower=0> steepness2;           // shallower (central -> recovery)

  // Individual GP — same three-section EQ changepoint kernel
  real<lower=0> length_scale_k2_s1;
  real<lower=0> length_scale_k2_s2;
  real<lower=0> length_scale_k2_s3;
  real<lower=0> variance_k2_s1;
  real<lower=0> variance_k2_s2;
  real<lower=0> variance_k2_s3;

  real<lower=0> sigma;
  vector[N] eta;
  vector[N] eta_ind;
}

transformed parameters {
  vector[N] f;
  vector[N] f_1;
  vector[N] f_2;
  vector[N] mu;
  
  {
    matrix[N, N] K_1;
    matrix[N, N] L_K_1;
    matrix[N, N] K_2;
    matrix[N, N] L_K_2;

    // Precompute sigmoid values for each observation
    // sig1[i] = sigmoid(steepness1 * (x[i] - loc1))
    // sig2[i] = sigmoid(steepness2 * (x[i] - loc2))
    vector[N] sig1;
    vector[N] sig2;
    for (n in 1:N) {
      sig1[n] = inv_logit(steepness1 * (x[n] - loc1));
      sig2[n] = inv_logit(steepness2 * (x[n] - loc2));
    }

// Blended intercept per observation

    for (n in 1:N) {
      real s1 = sig1[n];
      real s2 = sig2[n];
      mu[n] =   intercept_s1 * (1 - s1)
              + intercept_s2 * s1 * (1 - s2)
              + intercept_s3 * s1 * s2;
    }

    // --- Shared GP: three-section EQ changepoint kernel ---
    for (i in 1:N) {
      for (j in i:N) {
        real d2 = square(x[i] - x[j]);

        // Three EQ base kernels
        real k_s1 = square(variance_k1_s1) * exp(-0.5 * d2 / square(length_scale_k1_s1));
        real k_s2 = square(variance_k1_s2) * exp(-0.5 * d2 / square(length_scale_k1_s2));
        real k_s3 = square(variance_k1_s3) * exp(-0.5 * d2 / square(length_scale_k1_s3));

        // Sigmoid blending weights (GPflow changepoint formula)
        real s1i = sig1[i]; real s1j = sig1[j];
        real s2i = sig2[i]; real s2j = sig2[j];

        real w1 = (1 - s1i) * (1 - s1j);           // before loc1
        real w2 = s1i * s1j * (1 - s2i) * (1 - s2j); // between loc1 and loc2
        real w3 = s1i * s1j * s2i * s2j;            // after loc2

        K_1[i, j] = k_s1 * w1 + k_s2 * w2 + k_s3 * w3;
        K_1[j, i] = K_1[i, j];
      }
    }
    for (n in 1:N) K_1[n, n] = K_1[n, n] + delta;

    L_K_1 = cholesky_decompose(K_1);
    f_1 = L_K_1 * eta;

    // --- Individual GP: same changepoint kernel, subject-masked, z-scaled ---
    for (i in 1:N) {
      for (j in i:N) {
        if (subject[i] != subject[j]) {
          K_2[i, j] = 0.0;
          K_2[j, i] = 0.0;
        } else {
          real d2 = square(x[i] - x[j]);

          real k_s1 = square(variance_k2_s1) * exp(-0.5 * d2 / square(length_scale_k2_s1));
          real k_s2 = square(variance_k2_s2) * exp(-0.5 * d2 / square(length_scale_k2_s2));
          real k_s3 = square(variance_k2_s3) * exp(-0.5 * d2 / square(length_scale_k2_s3));

          real s1i = sig1[i]; real s1j = sig1[j];
          real s2i = sig2[i]; real s2j = sig2[j];

          real w1 = (1 - s1i) * (1 - s1j);
          real w2 = s1i * s1j * (1 - s2i) * (1 - s2j);
          real w3 = s1i * s1j * s2i * s2j;

          real k_base = k_s1 * w1 + k_s2 * w2 + k_s3 * w3;

          // lgpr z-scaling
          K_2[i, j] = z[i] * k_base * z[j];
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
  intercept_s1 ~ normal(0, 5);
  intercept_s2 ~ normal(0, 5);
  intercept_s3 ~ normal(0, 5);

  // Shared GP priors
  length_scale_k1_s1 ~ inv_gamma(1, 15);
  length_scale_k1_s2 ~ inv_gamma(2, 11);
  length_scale_k1_s3 ~ inv_gamma(1, 13);
  variance_k1_s1     ~ inv_gamma(4, 1);
  variance_k1_s2     ~ inv_gamma(5, 15);  // higher variance for central section
  variance_k1_s3     ~ inv_gamma(2, 1);
  steepness1         ~ gamma(7, 0.5);    // prior centred around higher steepness
  steepness2         ~ gamma(3, 0.5);    // prior centred around lower steepness

  // Individual GP priors
  length_scale_k2_s1 ~ inv_gamma(1, 11);
  length_scale_k2_s2 ~ inv_gamma(3, 11);
  length_scale_k2_s3 ~ inv_gamma(1, 9);
  variance_k2_s1     ~ inv_gamma(4, 3);
  variance_k2_s2     ~ inv_gamma(3, 15);
  variance_k2_s3     ~ inv_gamma(2, 3);

  sigma    ~ std_normal();
  eta      ~ std_normal();
  eta_ind  ~ std_normal();

  y ~ normal(mu + f, sigma);
}

generated quantities {
  vector[N] fit   = mu + f;
  vector[N] eff_1 = f_1;
  vector[N] eff_2 = f_2;
}
