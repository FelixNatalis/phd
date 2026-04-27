data{
   int<lower=1> N;      // number of observations
   array[N] real x;         // univariate covariate
   vector[N] y;         // target variable 
  
}

transformed data{
   real delta = 1e-9;
   real period_k2 = 7;
   real period_k3 = 30.4;  
}

parameters{
   real intercept;
   real<lower=0> length_scale_k1;
   real<lower=0> variance_k1;
   real<lower=0> length_scale_k2;
   real<lower=0> variance_k2;
   real<lower=0> length_scale_k3;
   real<lower=0> variance_k3;
   real<lower=0> sigma;
   vector[N] eta;
}

transformed parameters{
  
   vector[N] f;
   vector[N] f_1;
   vector[N] f_2;
   vector[N] f_3;
    
   {
    matrix[N, N] K_1;
    matrix[N, N] K_2 = gp_periodic_cov(x, variance_k2, length_scale_k2, period_k2);
    matrix[N, N] K_3 = gp_periodic_cov(x, variance_k3, length_scale_k3, period_k3);
    matrix[N, N] L_K_1;
    matrix[N, N] L_K_2;
    matrix[N, N] L_K_3;

  // OU kernel
    for (i in 1:N) {
     for (j in i:N) {
       real d = fabs(x[i] - x[j]);
       K_1[i,j] = square(variance_k1) * exp(-d/ length_scale_k1);
       K_1[j,i] = K_1[i,j];
     }
    }

    // diagonal elements
    for (n in 1:N) {
      K_1[n, n] = K_1[n, n] + delta;
      K_2[n, n] = K_2[n, n] + delta;
      K_3[n, n] = K_3[n, n] + delta;
    }

    L_K_1 = cholesky_decompose(K_1);
    L_K_2 = cholesky_decompose(K_2);
    L_K_3 = cholesky_decompose(K_3);

    f_1 = L_K_1 * eta;
    f_2 = L_K_2 * eta;
    f_3 = L_K_3 * eta;

    f = f_1 + f_2 + f_3;
   }
}

model{
   // priors 
   intercept ~ normal(0, 5);
   length_scale_k1 ~ inv_gamma(5, 5);
   variance_k1 ~ inv_gamma(5, 5);
   length_scale_k2 ~ inv_gamma(5, 5);
   variance_k2 ~ inv_gamma(5, 5);
   length_scale_k3 ~ inv_gamma(5, 5);
   variance_k3 ~ inv_gamma(5, 5);
   sigma ~ std_normal();
   eta ~ std_normal();

   // model 
   y ~ normal(intercept + f, sigma);
  
}


generated quantities{
  vector[N] fit = intercept + f;
  vector[N] eff_1 = f_1;
  vector[N] eff_2 = f_2;
  vector[N] eff_3 = f_3;
}
