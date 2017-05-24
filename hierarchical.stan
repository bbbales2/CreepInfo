data {
  int<lower=1> L;
  int<lower=1> T;
  int<lower=1> labels[L];
  vector<lower=0.0>[L] stress;
  vector<lower=0.0>[L] thickness;
  vector[L] mus;
}

parameters {
  real<lower=0.0> sigma[T];
  real a_mu;
  real<lower = 0.0> a_sigma;
  real b_mu;
  real<lower = 0.0> b_sigma;
  real a[T];
  real b[T];
  real c;
}

model {
  a_sigma ~ normal(0.0, 10.01);
  b_sigma ~ normal(0.0, 10.01);
  a ~ normal(a_mu, a_sigma);
  b ~ normal(b_mu, b_sigma);
  sigma ~ normal(0.0, 5.0);
  
  for(l in 1:L) {
    mus[l] ~ normal(a[labels[l]] * log(stress[l]) + b[labels[l]] * log(1.0 / thickness[l]) + c, sigma[labels[l]]);
  }
}

generated quantities {
  vector[L] muhat;

  for(l in 1:L) {
    muhat[l] = normal_rng(a[labels[l]] * log(stress[l]) + b[labels[l]] * log(1.0 / thickness[l]) + c, sigma[labels[l]]);
  }
}
