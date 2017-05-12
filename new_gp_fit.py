#%%
import numpy
import matplotlib.pyplot as plt
import scipy
import pystan

filename = '/home/bbales2/CreepInfo/200nm-200A-200B/creep_200B_275MPa.mat'

data = scipy.io.loadmat(filename)
x, y = data['creepinfo_200B_275MPa'][:100:1].T

#x -= x.min()
#x /= x.max()

#y -= y.min()
#y /= y.max()

#y = y + numpy.random.randn(len(y))
x = numpy.arange(len(y))

plt.plot(x, y)

template = """
N <- {0}
x <- c({1})
y <- c({2})
"""

def arr_to_str(x):
    return ", ".join(str(xx) for xx in x)

with open('/home/bbales2/cmdstan-2.15/gp.dat', 'w') as f:
    f.write(template.format(len(x), arr_to_str(x), arr_to_str(y)))
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
  int<lower=1> N;
  real x[N];
  vector[N] y;
}

transformed data {
  real delta = 1e-9;
}

parameters {
  real<lower=0> alpha;
  real<lower=0> rho;
  real<lower=0> sigma;
  vector[N] eta;
}

model {
  vector[N] f;

  {
    matrix[N, N] Sigma = cov_exp_quad(x, alpha, rho);
    matrix[N, N] L_Sigma;
    // diagonal elements
    for (k in 1:N)
      Sigma[k, k] = Sigma[k, k] + delta;

    L_Sigma = cholesky_decompose(Sigma);
    f = L_Sigma * eta;
  }

  alpha ~ cauchy(0.0, 2.0);//gamma(0.25, 0.25);
  rho ~ cauchy(0.0, 2.0);//gamma(0.25, 0.25);
  sigma ~ cauchy(0, 1.0);
  eta ~ normal(0, 1);

  y ~ normal(f, sigma);
}

"""

sm = pystan.StanModel(model_code = model_code)
#%%
fit = sm.sampling(data = {
    'x' : x,
    'y' : y,
    #'eta_sq' : 0.1,
    #'l' : 0.01,
    'N' : len(x)
    })
#%%
modelcode = """
data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
}

transformed data {
  vector[N] mu;
  mu = rep_vector(0, N);
}

parameters {
  real<lower=0> alpha;
  real<lower=0> rho;
  //real<lower=0> sigma;
}

model {
  matrix[N, N] Sigma = cov_exp_quad(x, alpha, rho);
  matrix[N, N] L_Sigma;
  //real sq_sigma = square(sigma);
  // diagonal elements
  for (k in 1:N)
    Sigma[k, k] = Sigma[k, k] + 1e-9;//sq_sigma;

  L_Sigma = cholesky_decompose(Sigma);
  alpha ~ normal(0, 1);
  rho ~ gamma(10, 10);
  //sigma ~ cauchy(0, 5);

  y ~ multi_normal_cholesky(mu, L_Sigma);
}"""
sm = pystan.StanModel(model_code = modelcode)
#%%
fit = sm.sampling(data = {
    'x' : x,
    'y' : y,
    #'eta_sq' : 0.1,
    #'l' : 0.01,
    'N' : len(x)
    })
#%%
print fit
#%%

#%%
print fit

a = fit.extract()['a']
b = fit.extract()['b']
yhat = fit.extract()['yhat']

plt.plot(a, b, '*')
plt.show()

plt.plot(x, y)
for i in range(10):
    plt.plot(x, yhat[i, :], '*')
plt.show()

#%%
model_code = """
data {
  int<lower=1> N; // Number of single samples
  vector<lower=0.0>[N] x;
  vector<lower=0.0>[N] y;
}

transformed data {
  vector[N] z;

  for(i in 1:N)
    z[i] = 0.0;
}

parameters {
  real<lower=0> eta_sq;
  real<lower=0> inv_rho_sq;
  real<lower=0> sigma_sq;

  real a;
  real b;

  vector[N] mu;
}

transformed parameters {
  real<lower=0> rho_sq;
  rho_sq = inv(inv_rho_sq);
}

model {
  matrix[N, N] Sigma;

  // off-diagonal elements
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i, j] = eta_sq * exp(-rho_sq * pow(x[i] - x[j], 2));
      Sigma[j, i] = Sigma[i, j];
    }
  }

  // diagonal elements
  for (k in 1:N)
    Sigma[k, k] = eta_sq + 1e-6; // + jitter

  mu ~ multi_normal(z, Sigma);
  eta_sq ~ normal(0, 0.1);
  inv_rho_sq ~ normal(0, 1.0);
  sigma_sq ~ normal(0, 1.0e-4);

  y ~ normal(a * x + b + mu, sigma_sq);
}

generated quantities {
  vector[N] yhat;
  vector[N] muhat;
  vector[N] lhat;

  for(i in 1:N) {
    yhat[i] = normal_rng(a * x[i] + b + mu[i], sigma_sq);
    muhat[i] = normal_rng(mu[i], sigma_sq);
    lhat[i] = normal_rng(a * x[i] + b, sigma_sq);
  }
}
"""
"""

  vector[N] gp;

  {
    matrix[N, N] Sigma;

    // off-diagonal elements
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        Sigma[i, j] = eta_sq * exp(-rho_sq * pow(x[i] - x[j],2));
        Sigma[j, i] = Sigma[i, j];
      }
    }

    // diagonal elements
    for (k in 1:N)
      Sigma[k, k] = eta_sq + sigma_sq; // + jitter

    gp = multi_normal_rng(mu, Sigma);
    yhat = a * x + b;
  }"""
#generated quantities {
#  vector[N] yhat;
#
#  for(n in 1:N) {
#    yhat[n] <- normal_rng(a * x[n] + b, sigma);
#  }
#}
#"""

smg = pystan.StanModel(model_code = model_code)

#%%
model_code = """
data {
  int<lower=1> N; // Number of single samples
  vector<lower=0.0>[N] x;
  vector<lower=0.0>[N] y;
  real<lower=0.0> eta_sq;
  real<lower=0.0> l;
}

transformed data {
  vector[N] zero;
  matrix[N, N] L;
  cov_matrix[N] Sigma;

  for(i in 1:N)
    zero[i] = 0.0;

  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i, j] = eta_sq * exp(-pow(x[i] - x[j], 2) / (l * l));
      Sigma[j, i] = Sigma[i, j];
    }
  }

  for (k in 1:N)
    Sigma[k, k] = eta_sq + 1e-6; // + jitter

  L = cholesky_decompose(Sigma);
}

parameters {
  real<lower=0> sigma_sq;

  real a;
  real p;
  real b;

  vector[N] mu;
  vector[N] t;
}

model {
  mu ~ multi_normal_cholesky(zero, L);
  sigma_sq ~ normal(0, 1.0);

  a ~ normal(0.0, 2.0);
  b ~ normal(0.0, 1.0);

  t ~

  y ~ normal(a * x + b + mu, sigma_sq);
}

generated quantities {
  vector[N] yhat;
  vector[N] muhat;
  vector[N] lhat;

  for(i in 1:N) {
    yhat[i] = normal_rng(a * x[i] + b + mu[i], sigma_sq);
    muhat[i] = mu[i];
    lhat[i] = a * x[i] + b;
  }
}
"""
#generated quantities {
#  vector[N] yhat;
#
#  for(n in 1:N) {
#    yhat[n] <- normal_rng(a * x[n] + b, sigma);
#  }
#}
#"""

smg = pystan.StanModel(model_code = model_code)

#%%
#%%
#print fit

a = fit.extract()['a']
#b = fit.extract()['b']
yhat = fit.extract()['yhat']
muhat = fit.extract()['muhat']
lhat = fit.extract()['lhat']
#gp = fit.extract()['gp']
#res = fit.extract()['res']

#plt.plot(a, b, '*')
#plt.show()

plt.plot(x, y, '--')
for i in range(10):
    #plt.plot(x, a[-i] * x + b[-i], '.')
    plt.plot(x, yhat[-i, :], '-b')
    plt.plot(x, muhat[-i, :], '-g')
    plt.plot(x, lhat[-i, :], '-r')
    #plt.plot(x, gp[-i, :], '.')
    #plt.plot(x, y - yhat[-i, :], '--')
plt.show()
#%%

import GPy

kernel = GPy.kern.Bias(1, name = 'bias') + GPy.kern.RBF(input_dim=1, variance=1., lengthscale=1.) + GPy.kern.LinearSlopeBasisFuncKernel(1, 0.0, 1.0, name = 'linear_slopes')

m = GPy.models.GPRegression(x.reshape(-1, 1), y.reshape(-1, 1), kernel)

print m
m.plot()

#%%
m.optimize()
m.optimize_restarts(num_restarts = 2)

print m

#print m.kern.bias.posterior_inf()
print m.kern.linear_slopes.posterior_inf()
#%%
m.plot([0.0, 1.0], plot_density = True)

print m

#%%

m.plot()