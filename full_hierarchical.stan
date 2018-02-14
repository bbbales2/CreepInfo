data {
  int<lower=1> L;
  int<lower=1> T;
  int<lower=1> labels[L];
  vector<lower=0.0>[L] stress;
  vector<lower=0.0>[L] thickness;
  vector<lower=0.0>[T] thicknesses;
  vector[L] mus;
}

parameters {
  real<lower=0.0> sigma[T];
  real<lower=0.0> lambda_n;
  real<lower=0.0> lambda_p;
  vector[T] zn;
  vector[T] zd;
  real mu;
  real p;
  real c;
}

transformed parameters {
  vector[T] n = zn * lambda_n + mu;
  vector[T] d;
  
  for(t in 1:T) {
    d[t] = zd[t] * lambda_p + p * log(1.0 / thicknesses[t]) + c;
  }
}

model {
  sigma ~ normal(0.0, 1.0);
  mu ~ normal(0.0, 10.0);
  lambda_n ~ normal(0.0, 1.0);
  lambda_p ~ normal(0.0, 1.0);
  zn ~ normal(0.0, 1.0);
  zd ~ normal(0.0, 1.0);
  p ~ normal(0.0, 5.0);
  
  for(l in 1:L) {
    mus[l] ~ normal(n[labels[l]] * log(stress[l]) + d[labels[l]], sigma[labels[l]]);
  }
}

