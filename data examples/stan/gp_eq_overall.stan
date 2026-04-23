data{
   int<lower=1> N;      // number of observations
   array[N] real x;         // univariate covariate
   vector[N] y;         // target variable 
  
}

transformed data{
   real delta = 1e-9;
  
  
}

parameters{
   real<lower=0> rho;
   real<lower=0> alpha;
   real<lower=0> sigma;
   vector[N] eta;
}

transformed parameters{
  
   vector[N] f;
    
   {
    matrix[N, N] L_K;
    matrix[N, N] K = gp_exp_quad_cov(x, alpha, rho);

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
   rho ~ inv_gamma(5, 5);
   alpha ~ std_normal();
   sigma ~ std_normal();
   eta ~ std_normal();

   // model 
   y ~ normal(f, sigma);
  
}


generated quantities{
  vector[N] fit = f;
}
