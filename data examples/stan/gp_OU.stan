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
}

model{

   {
    matrix[N, N] L_K;
    
    vector diagSPD_Matern12(real alpha, real rho, real L, int M) {
  vector[M] indices = linspaced_vector(M, 1, M);
  real factor = 2;
  vector[M] denom = rho * ((1 / rho)^2 + (pi() / 2 / L) * indices);
  return alpha * sqrt(factor * inv(denom));
}
    
    matrix[N, N] K = gp_exp_quad_cov(x, alpha, rho);

    // diagonal elements
    for (n in 1:N) {
      K[n, n] = K[n, n] + delta;
    }

    L_K = cholesky_decompose(K);
    f = L_K * eta;
   }
  
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
