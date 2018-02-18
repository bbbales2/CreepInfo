data {
  int<lower=1> L;
  int<lower=1> T;
  int<lower=1> labels[L];
  vector[L] log_stress;
  vector[T] log_inv_thicknesses;
  vector[L] mus;
}

parameters {
  real<lower=0.0> sigma[T];
  real<lower=0.0> lambda_n;
  real<lower=0.0> lambda_d;
  vector[T] zn;
  vector[T] zd;
  //real n;
  real p_n;
  real p_d;
  real c_n;
  real c_d;
}

transformed parameters {
  vector[T] n;
  vector[T] d;
  
  for(t in 1:T) {
    n[t] = zn[t] * lambda_n + p_n * log_inv_thicknesses[t] + c_n;
    d[t] = zd[t] * lambda_d + p_d * log_inv_thicknesses[t] + c_d;
  }
}

model {
  sigma ~ normal(0.0, 1.0);
  lambda_n ~ normal(0.0, 10.0);
  lambda_d ~ normal(0.0, 1.0);
  zn ~ normal(0.0, 1.0);
  zd ~ normal(0.0, 1.0);
  p_n ~ normal(0.0, 5.0);
  p_d ~ normal(0.0, 5.0);
  c_n ~ normal(0.0, 50.0);
  c_d ~ normal(0.0, 50.0);
  
  for(l in 1:L) {
    mus[l] ~ normal(n[labels[l]] * log_stress[l] + d[labels[l]], sigma[labels[l]]);
  }
}

generated quantities {
  vector[L] muhat;
  vector[L] mumu;
  
  for(l in 1:L) {
    mumu[l] = n[labels[l]] * log_stress[l] + d[labels[l]];
    muhat[l] = normal_rng(mumu[l], sigma[labels[l]]);
  }
}
