data{
   int<lower=1> N;      // number of observations
   array[N] real x;         // univariate covariate
   vector[N] y;         // target variable 
  
}

transformed data{
   real delta = 1e-9;
   real period_k2 = 7;
  
}

parameters{
   real intercept;
   real<lower=0> length_scale_k1;
   real<lower=0> variance_k1;
   real<lower=0> length_scale_k2;
   real<lower=0> variance_k2;
   real<lower=0> sigma;
   vector[N] eta;
}

transformed parameters{
  
   vector[N] f;
    
   {
    matrix[N, N] L_K;
    matrix[N, N] K = gp_exp_quad_cov(x, variance_k1, length_scale_k1) + gp_periodic_cov(x, variance_k2, length_scale_k2, period_k2);

    // diagonal elements
    for (n in 1:N) {
      K[n, n] = K[n, n] + delta;
    }

    L_K = cholesky_decompose(K);
    f = L_K * eta;
   }
}

model{
   // priors 
   intercept ~ normal(0, 5);
   length_scale_k1 ~ inv_gamma(5, 5);
   variance_k1 ~ inv_gamma(5, 5);
   length_scale_k2 ~ inv_gamma(5, 5);
   variance_k2 ~ inv_gamma(5, 5);
   sigma ~ std_normal();
   eta ~ std_normal();

   // model 
   y ~ normal(intercept + f, sigma);
  
}


generated quantities{
  vector[N] fit = intercept + f;
}
