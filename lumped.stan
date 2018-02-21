data {
  int<lower=1> L;
  int<lower=1> T;
  int<lower=1> labels[L];
  vector[L] log_stress;
  vector[L] log_inv_thickness;
  vector[L] mus;
}

parameters {
  real<lower=0.0> sigma;//[T];
  real n;
  real p;
  real c;
}

model {
  sigma ~ normal(0.0, 5.0);
  
  for(l in 1:L) {
    mus[l] ~ normal(n * log_stress[l] + p * log_inv_thickness[l] + c, sigma);//[labels[l]]);
  }
}

generated quantities {
  vector[L] muhat;
  vector[L] mumu;
  vector[L] uncertainty;
  vector[L] sdo;
  vector[L] lso;
  vector[L] sdo_hat;
  vector[L] lso_hat;
  
  for(l in 1:L) {
    mumu[l] = n * log_stress[l] + p * log_inv_thickness[l] + c;
    muhat[l] = normal_rng(mumu[l], sigma);//[labels[l]]);
    uncertainty[l] = mus[l] - muhat[l];
    sdo_hat[l] = muhat[l] - p * log_inv_thickness[l];
    lso_hat[l] = muhat[l] - n * log_stress[l];
    sdo[l] = mus[l] - p * log_inv_thickness[l];
    lso[l] = mus[l] - n * log_stress[l];
  }
}
