data{
   int<lower=1> N_1;      // number of observations
   array[N_1] real x_1;         // univariate covariate
   vector[N_1] y_1;         // target variable 
   int<lower=1> N_2; // number of prediction points
   array[N_2] real x_2; // x for prediction
}

transformed data{
   real delta = 1e-9;
   real period_k_2 = 7;
   real period_k_3 = 30.4;  
   
   int<lower=1> N = N_1 + N_2;1`
   array[N] real x;
   for (n_1 in 1:N_1) {
     x[n_1] = x_1[n_1];
   }
   for (n_2 in 1:N_2) {
     x[N_1 + n_2] = x_2[n_2];
   }
}

parameters{
   real intercept;
   real<lower=0> length_scale_k_1;
   real<lower=0> variance_k_1;
   real<lower=0> length_scale_k_2;
   real<lower=0> variance_k_2;
   real<lower=0> length_scale_k_3;
   real<lower=0> variance_k_3;
   real<lower=0> length_scale_k_4;
   real<lower=0> variance_k_4;
   real<lower=0> sigma;
   vector[N] eta;
}

transformed parameters{
  
   vector[N] f; 
   vector[N] f_1; // ou daily, shaky
   vector[N] f_2; // periodic weekly
   vector[N] f_3; // periodic monthly
   vector[N] f_4; // eq long and smooth
    
   {
    matrix[N, N] K_1;
    matrix[N, N] K_2 = gp_periodic_cov(x, variance_k_2, length_scale_k_2, period_k_2);
    matrix[N, N] K_3 = gp_periodic_cov(x, variance_k_3, length_scale_k_3, period_k_3);
    matrix[N, N] K_4 = gp_exp_quad_cov(x, variance_k_4, length_scale_k_4);
    matrix[N, N] L_K_1;
    matrix[N, N] L_K_2;
    matrix[N, N] L_K_3;
    matrix[N, N] L_K_4;

  // OU kernel
    for (i in 1:N) {
     for (j in i:N) {
       real d = fabs(x[i] - x[j]);
       K_1[i,j] = square(variance_k_1) * exp(-d/ length_scale_k_1);
       K_1[j,i] = K_1[i,j];
     }
    }

    // diagonal elements
    for (n in 1:N) {
      K_1[n, n] = K_1[n, n] + delta;
      K_2[n, n] = K_2[n, n] + delta;
      K_3[n, n] = K_3[n, n] + delta;
      K_4[n, n] = K_4[n, n] + delta;
    }

    L_K_1 = cholesky_decompose(K_1);
    L_K_2 = cholesky_decompose(K_2);
    L_K_3 = cholesky_decompose(K_3);
    L_K_4 = cholesky_decompose(K_4);

    f_1 = L_K_1 * eta;
    f_2 = L_K_2 * eta;
    f_3 = L_K_3 * eta;
    f_4 = L_K_4 * eta;

    f = f_1 + f_2 + f_3 + f_4;
   }
}

model{
   // priors 
  // length_scale_k_1 ~ inv_gamma(3, 2); // ou, ls nonsquared, var squared
  // variance_k_1 ~ inv_gamma(7, 15);
//
  // length_scale_k_2 ~ inv_gamma(5, 15);
  // variance_k_2 ~ inv_gamma(3, 2);
//
  // length_scale_k_3 ~ inv_gamma(3, 9);
  // variance_k_3 ~ inv_gamma(3, 2);
//
  // length_scale_k_4 ~ inv_gamma(1, 6);
  // variance_k_4 ~ inv_gamma(4, 3);
//
  // intercept ~ normal(0, 2);
  // sigma ~ std_normal();
  // eta ~ std_normal();

   // model 
   y_1 ~ normal(intercept + f[1:N_1], sigma);
  
}


generated quantities{
  vector[N_1] fit = intercept + f[1:N_1];
  vector[N_2] y_2;
  for (n_2 in 1:N_2) {
    y_2[n_2] = normal_rng(intercept + f[N_1 + n_2], sigma);
  }
  vector[N] eff_1 = f_1;
  vector[N] eff_2 = f_2;
  vector[N] eff_3 = f_3;
  vector[N] eff_4 = f_4;
}
