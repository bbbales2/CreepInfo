data {
  int<lower=1> P;
  int<lower=1> S;
  int<lower=1> L;
  int<lower=1> T;
  int<lower=1> labels[L];
  vector[L] log_stress;
  vector[L] log_inv_thickness;
  vector[S] log_predict_stress;
  vector[P] log_inv_predict_thickness;
  vector[L] mus;
}

parameters {
  //real mu;
  //real<lower=0.0> lambda;
  real<lower=0.0> sigma;//[T];
  real n;
  real p;
  real c;
}

model {
  //mu ~ normal(0.0, 5.0);
  //lambda ~ normal(0.0, 5.0);
  //for(t in 1:T) {
  //  sigma[t] ~ lognormal(mu, lambda);
  //}
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
  matrix[S, P] muhat2;
  
  for(s in 1:S) {
    for(i in 1:P) {
      muhat2[s, i] = normal_rng(n * log_predict_stress[s] + p * log_inv_predict_thickness[i] + c, sigma);
    }
  }
  
  for(l in 1:L) {
    mumu[l] = n * log_stress[l] + p * log_inv_thickness[l] + c;
    muhat[l] = normal_rng(mumu[l], sigma);//[labels[l]]);
    uncertainty[l] = mus[l] - muhat[l];
    sdo_hat[l] = muhat[l] - p * log_inv_thickness[l];
    lso_hat[l] = muhat[l] - n * log_stress[l];
    sdo[l] = mus[l] - p * log_inv_thickness[l];
    lso[l] = mus[l] - n * log_stress[l];
    //sigma_t = lognormal_rng(mu, lambda);
  }
}
