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
  real<lower=0.0> lambda;
  real mu;
  vector[max(labels)] n;
  real c;
}

#transformed parameters {
#  vector[max(labels)] n = z * lambda;
#}

model {
  sigma ~ normal(0.0, 5.0);
  mu ~ normal(0.0, 10.0);
  lambda ~ normal(0.0, 1.0);
  n ~ normal(mu, lambda);
  
  for(l in 1:L) {
    mus[l] ~ normal(n[labels[l]] * log(stress[l]) + c, sigma[labels[l]]);
  }
}

generated quantities {
  vector[L] muhat;
  vector[L] mumu;
  
  for(l in 1:L) {
    mumu[l] = n[labels[l]] * log(stress[l]) + c;
    muhat[l] = normal_rng(mumu[l], sigma[labels[l]]);
  }
}
