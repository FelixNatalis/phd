data {
  int<lower=1> N1;
  array[N1] real x1;
  vector[N1] y1;
  int<lower=1> N2;
  array[N2] real x2;
}
transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
  array[N] real x;
  for (n1 in 1:N1) {
    x[n1] = x1[n1];
  }
  for (n2 in 1:N2) {
    x[N1 + n2] = x2[n2];
  }
}
parameters {
  real intercept;
  real<lower=0> length_scale;
  real<lower=0> variance;
  real<lower=0> sigma;
  vector[N] eta;
}
transformed parameters {
  vector[N] f;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = gp_exp_quad_cov(x, variance, length_scale);

    // diagonal elements
    for (n in 1:N) {
      K[n, n] = K[n, n] + delta;
    }

    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }
}
model {
  intercept ~ normal(0, 5);
  length_scale ~ inv_gamma(5, 1);
  variance ~ normal(100, 50);
  sigma ~ normal(1000, 50);
  eta ~ normal(100, 10);

  y1 ~ normal(intercept + f[1:N1], sigma);
}
generated quantities {
  vector[N1] fit = intercept + f[1:N1];
  vector[N2] y2;
  for (n2 in 1:N2) {
    y2[n2] = normal_rng(intercept + f[N1 + n2], sigma);
  }
}
