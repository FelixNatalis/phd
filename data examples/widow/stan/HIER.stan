data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1,upper=S> subj[N];

  vector[N] y;
  vector[N] x;
}

parameters {
  real<lower=0> variance_k_1;
  real<lower=0> length_scale_k_1;
  real<lower=0> variance_k_2;
  real<lower=0> length_scale_k_2;
  real<lower=0> variance_k_3;
  real<lower=0> length_scale_k_3;
  real<lower=0> sigma;

  // intercepts
  real<lower=0> sigma_alpha;
  vector[S] alpha_raw;

  // changepoint parameters
  real<lower=0> steepness;
  real x_0;

  // GP latent
  vector[N] eta;
}

transformed parameters {

  vector[N] mu;
  vector[S] alpha;



  // -------------------------------------------------
  // 2. INTERCEPTS
  // -------------------------------------------------
  alpha = sigma_alpha * alpha_raw;

  // -------------------------------------------------
  // 3. GP PER SUBJECT
  // -------------------------------------------------

  mu = rep_vector(0, N);

  for (s in 1:S) {

    int Ns = 0;

    for (n in 1:N)
      if (subj[n] == s) Ns += 1;

    if (Ns > 0) {

      matrix[Ns, Ns] K;
      matrix[Ns, Ns] L;
      vector[Ns] eta_s;
      vector[Ns] x_s;
      vector[Ns] mu_s;

      int k = 1;

      for (n in 1:N) {
        if (subj[n] == s) {
          x_s[k] = x[n];   
          eta_s[k] = eta[n];
          k += 1;
        }
      }

      for (i in 1:Ns) {
        for (j in i:Ns) {
          real d = fabs(x_s[i] - x_s[j]);
          
          K[i,j] = square(variance_k_1) * exp(-d / length_scale_k_1);
          K[j,i] = K[i,j];
        }
      }

      for (i in 1:Ns)
        K[i,i] += 1e-6;

      L = cholesky_decompose(K);
      mu_s = L * eta_s;

      k = 1;
      for (n in 1:N) {
        if (subj[n] == s) {
          mu[n] = mu_s[k] + alpha[s];
          k += 1;
        }
      }
    }
  }
}

model {

  eta ~ normal(0,1);
  sigma ~ normal(0.5,0.2);
  
  alpha_raw ~ normal(0,1);
  sigma_alpha ~ normal(1,0.5);

  variance_k_1 ~ normal(1,0.5);
  length_scale_k_1 ~ lognormal(log(10), 0.5);
  
  variance_k_2 ~ normal(1,0.5);
  length_scale_k_2 ~ lognormal(log(10), 0.5);
  
  variance_k_3 ~ normal(1,0.5);
  length_scale_k_3 ~ lognormal(log(10), 0.5);

  steepness ~ lognormal(log(0.05), 0.7);
  x_0 ~ normal(0,1);
  
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N] fit = mu;
}