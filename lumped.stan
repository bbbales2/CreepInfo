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
  real a;
  real b;
  real c;
}

model {
  sigma ~ normal(0.0, 5.0);
  
  for(l in 1:L) {
    mus[l] ~ normal(a * log(stress[l]) + b * log(1.0 / thickness[l]) + c, sigma[labels[l]]);
  }
}

generated quantities {
  vector[L] muhat;
  
  for(l in 1:L) {
    muhat[l] = normal_rng(a * log(stress[l]) + b * log(1.0 / thickness[l]) + c, sigma[labels[l]]);
  }
}
