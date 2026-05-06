data{
   int<lower=1> N;      // number of observations
   array[N] real x;         // univariate covariate
   vector[N] y;         // target variable 
  
}

transformed data{
   real delta = 1e-9;
  
  
}

parameters{
   real intercept;
   real<lower=0> length_scale;
   real<lower=0> variance;
   real<lower=0> sigma;
   vector[N] eta;
}

transformed parameters{
  
   vector[N] f;
    
   {
    matrix[N, N] K;
    matrix[N, N] L_K;

  // OU kernel
    for (i in 1:N) {
     for (j in i:N) {
       K[i,j] = square(variance) * exp(-square(d)/ square(length_scale));
       K[j,i] = K[i,j];
     }
    }

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
   length_scale ~ inv_gamma(5, 5);
   variance ~ inv_gamma(5, 5);
   sigma ~ std_normal();
   eta ~ std_normal();

   // model 
   y ~ normal(intercept + f, sigma);
  
}


generated quantities{
  vector[N] fit = intercept + f;
}
