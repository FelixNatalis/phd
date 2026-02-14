functions {
vector diagSPD_EQ(real alpha, real rho, real L, int M) {
  return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
}
vector diagSPD_Matern32(real alpha, real rho, real L, int M) {
   return 2*alpha * (sqrt(3)/rho)^1.5 * inv((sqrt(3)/rho)^2 + ((pi()/2/L) * linspaced_vector(M, 1, M))^2);
}
vector diagSPD_periodic(real alpha, real rho, int M) {
  real a = 1/rho^2;
  vector[M] q = exp(log(alpha) + 0.5 * (log(2) - a + to_vector(log_modified_bessel_first_kind(linspaced_int_array(M, 1, M), a))));
  return append_row(q,q);
}
matrix PHI(int N, int M, real L, vector x) {
  return sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), linspaced_vector(M, 1, M)))/sqrt(L);
}
matrix PHI_periodic(int N, int M, real w0, vector x) {
  matrix[N,M] mw0x = diag_post_multiply(rep_matrix(w0*x, M), linspaced_vector(M, 1, M));
  return append_col(cos(mw0x), sin(mw0x));
}


}
data {
  int<lower=1> N;      // number of observations
  vector[N] x;         // univariate covariate
  vector[N] y;         // target variable
        
  real<lower=0> c_f1;  // factor c to determine the boundary value L
  int<lower=1> M_f1;   // number of basis functions for smooth function
}
transformed data {
  // Normalize data
  real xmean = mean(x);
  real ymean = mean(y);
  real xsd = sd(x);
  real ysd = sd(y);
  vector[N] xn = (x - xmean)/xsd;
  vector[N] yn = (y - ymean)/ysd;
  // Basis functions for f1
  real L_f1 = c_f1*max(xn);
  matrix[N,M_f1] PHI_f1 = PHI(N, M_f1, L_f1, xn);
}
parameters {
  vector[M_f1] beta_f1;         // the basis functions coefficients
  real<lower=0> lengthscale_f1; // lengthscale of f1
  real<lower=0> sigma_f1;       // scale of f1
  real<lower=0> sigma;          // residual scale
}
model {
  // spectral densities for f1
  vector[M_f1] diagSPD_f1 = diagSPD_EQ(sigma_f1, lengthscale_f1, L_f1, M_f1);
  // priors
  beta_f1 ~ normal(0, 1);
  lengthscale_f1 ~ lognormal(log(700/xsd), 1);
  sigma_f1 ~ normal(0, 1);
  sigma ~ normal(0, .5);
  // model
  yn ~ normal_id_glm(PHI_f1, 0.0, diagSPD_f1 .* beta_f1, sigma); 
}
generated quantities {
  vector[N] f;
  vector[N] log_lik;
  {
    // spectral densities for f1
    vector[M_f1] diagSPD_f1 = diagSPD_EQ(sigma_f1, lengthscale_f1, L_f1, M_f1);
    // function scaled back to original scale
    f = (0.0 + PHI_f1 * (diagSPD_f1 .* beta_f1))*ysd + ymean;
    // log_liks for loo
    for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | f[n], sigma*ysd);
  }
}
