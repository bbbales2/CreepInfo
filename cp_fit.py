#%%
import numpy
import matplotlib.pyplot as plt
import scipy
import pystan

filename = '/home/bbales2/CreepInfo/200nm-200A-200B/creep_200B_275MPa.mat'

data = scipy.io.loadmat(filename)
#variable                                      creepinfo_200B_275MPa
#iminf                                                           0.5
#imaxf                                                             1
#thickness                                                       200
#stress                                                          275
#heat_treatment                                                    0
#treated                                                       False
#Name: 9, dtype: object
x, y = data['creepinfo_200B_275MPa'][:200:1].T

x -= x.min()
x /= x.max()

y -= y.min()
y /= y.max()

y = y# + numpy.random.randn(len(y)) * 0.025
y -= y.min()

plt.plot(x, y)
#%%
zs = []
for n in range(1, 100):
    zs.append(sum((y[n::n] - y[:-n:n])**2))
    print n, zs[-1]

plt.show()
plt.plot(zs)
plt.show()

#%%

model_code = """
data {
  int<lower=1> N; // Number of single samples
  vector<lower=0.0>[N] t;
  vector<lower=0.0>[N] y;
}

parameters {
  real<lower=0.0> sigma;
  real<lower = 0.0> C;
  real<lower = 0.0> p;
  real<lower = 0.0> em;
}

model {
  C ~ normal(0.0, 1.0);
  p ~ normal(0.0, 1.0);
  em ~ normal(0.0, 1.0);

  y ~ normal(C * p * t ./ (1 + p * t) + em * t, sigma);
}

generated quantities {
  vector[N] yhat;
  vector[N] yhatn;

  for(n in 1:N) {
    yhat[n] <- C * p * t[n] ./ (1 + p * t[n]) + em * t[n];
    yhatn[n] <- normal_rng(C * p * t[n] ./ (1 + p * t[n]) + em * t[n], sigma);
  }
}
"""

sm = pystan.StanModel(model_code = model_code)

#%%
