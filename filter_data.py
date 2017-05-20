#%%

import scipy
import matplotlib.pyplot as plt
import pystan
import numpy
import seaborn
import pandas
import pickle

files = [
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000A_200MPa_noTertiary.mat', 1.0 / 2.0, 1.0, 2000, 200, 0.0, False],
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_150MPa.mat', 1.0 / 2, 1.0, 2000, 150, 0.0, False],
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_170MPa.mat', 0.0, 1.0, 2000, 170, 0.0, False],
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_190MPa.mat', 0.0, -1.0 / 4.0, 2000, 190, 0.0, False],
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_200MPa.mat', 0.0, 1.0, 2000, 200, 0.0, False],
    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200A_245MPa.mat', 1.0 / 2.0, -1.0 / 4.0, 200, 245, 0.0, False],
    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200A_260MPa.mat', 0.0, 1.0, 200, 260, 0.0, False],
    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200A_290MPa.mat', 0.0, 1.0 / 4.0, 200, 290, 0.0, False],
    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200B_275MPa.mat', 1.0 / 2.0, 1.0, 200, 275, 0.0, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65C_320MPa.mat', 1.0 / 2.0, 1.0, 65, 320, 3.0, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65C_350MPa.mat', 0.0, 1.0 / 2.0, 65, 350, 71.0, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_275MPa.mat', 1.0 / 2.0, 1.0, 65, 275, 4.3, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_290MPa.mat', 0.0, 1.0, 65, 290, 103.1, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_320MPa.mat', 0.0, 1.0, 65, 320, 193.4, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_335MPa.mat', 0.0, 1.0, 65, 335, 266.1, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_230MPa.mat', 0.0, 1.0, 65, 230, 447.8, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_245MPa.mat', 0.0, 1.0, 65, 245, 514.6, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_260MPa.mat', 0.0, 1.0, 65, 260, 593.1, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_290MPa.mat', 0.0, 1.0, 65, 290, 641.6, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_320MPa.mat', 0.0, 1.0, 65, 320, 721.3, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_350MPa.mat', 0.0, 1.0, 65, 350, 827.5, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_370MPa_noTertiary.mat', 0.0, 1.0, 65, 370, 895.8, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65HT_320MPa.mat', 1.0 / 2.0, 1.0, 65, 320, 723.0, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65HT_350MPa_noTertiary.mat', 0.0, 1.0, 65, 350, 813.3, True],
    ['/home/bbales2/CreepInfo/CreepInfo/500nm-500A-500B-500C/creep_500A_160MPa_noPrimary.mat', 0.0, 1.0, 500, 160, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/500nm-500A-500B-500C/creep_500A_180MPa.mat', 0.2, 1.0, 500, 180, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/500nm-500A-500B-500C/creep_500A_200MPa.mat', 0.0, 1.0, 500, 200, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/500nm-500A-500B-500C/creep_500A_215MPa_noTertiary.mat', 0.0, 1.0, 500, 215, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/500nm-500A-500B-500C/creep_500B_175MPa_noPrimary.mat', 0.0, 1.0, 500, 175, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/500nm-500A-500B-500C/creep_500B_190MPa.mat', 0.0, 1.0, 500, 190, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/500nm-500A-500B-500C/creep_500B_205MPa_noTertiary.mat', 0.0, 1.0, 500, 205, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/500nm-500A-500B-500C/creep_500C_200MPa_noPrimaryorTertiary.mat', 0.0, 1.0, 500, 200, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/200nm-200A-200B-200D/creep_200D_215MPa.mat', 0.2, 1.0, 200, 215, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/200nm-200A-200B-200D/creep_200D_275MPa_noTertiary.mat', 0.0, 1.0, 200, 275, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/2000nm-2B-2A-2D/creep_2000D_175MPa_noPrimaryorTertiary.mat', 0.0, 1.0, 2000, 175, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/65nm-65C-65D-65H/creep_65H_295MPa_noPrimaryorTertiary.mat', 0.0, 1.0, 65, 295, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/200nm-200A-200B-200D/creep_200D_215MPa.mat', 0.2, 1.0, 200, 215, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/200nm-200A-200B-200D/creep_200D_245MPa.mat', 0.0, 1.0, 200, 245, 0.0, False],
    ['/home/bbales2/CreepInfo/CreepInfo/200nm-200A-200B-200D/creep_200D_275MPa_noTertiary.mat', 0.0, 1.0, 200, 275, 0.0, False]]

files = sorted(files, key = lambda x : x[-3])

df = pandas.DataFrame(files, columns = ['file', 'iminf', 'imaxf', 'thickness', 'stress', 'heat_treatment', 'treated'])

df = df[df['treated'] == False]
df = df.reset_index(drop = True)

slopes = []
#print "{0:10s}, {1:10s}, {2:12s}, {3:10s}, {4:10s}, {5}".format("thickness", "stress", "heat treated", "avg minimum strain rate", "std. deviation of minimum strain rate est.", "noise level (same units at strain rate)")
for idx, row in df.iterrows():
    filename, imin, imax, thickness, stress, _, heat_treated = row.values

    data = scipy.io.loadmat(filename)
    data = scipy.io.loadmat(filename)
    for key in data:
        if hasattr(data[key], "shape"):
            data = data[key]
            break

    imin = int(len(data) * imin)
    imax = int(len(data) * imax)
    data = data[imin : imax]

    if len(data) > 100:
        if (data.shape[0] - 99) / 100 > 0:
            data = data[::(data.shape[0] - 99) / 100]

    numpy.savetxt(filename + '.csv', data, delimiter = ',', header = 'time, stress', comments = '')

df.to_csv('/home/bbales2/CreepInfo/creep.csv', index = False)

#%%

model_code = """
data {
  //int<lower=1> N; // Number of single samples
  int<lower=1> L;
  int<lower=1> T;
  int<lower=1> labels[L];
  vector<lower=0.0>[L] stress;
  vector<lower=0.0>[L] thickness;
  vector<lower=0.0>[T] stress0;
  //int<lower=1> thickness[L];
  //vector<lower=0.0>[N] y;
  vector[L] mus;
  vector<lower=0.0>[L] sigmas;
}

parameters {
  //real mus[L];
  //real<lower=0.0> sigmas[L];

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
  real tmp[L];
  //sigmas ~ cauchy(0.0, 10.0);
  //sigma ~ cauchy(0.0, 10.0);
  a_sigma ~ normal(0.0, 10.01);
  b_sigma ~ normal(0.0, 10.01);
  a ~ normal(a_mu, a_sigma);
  b ~ normal(b_mu, b_sigma);
  sigma ~ normal(0.0, 5.0);
  //b ~ cauchy(0.0, 10.0);
  //c ~ normal(-2.2, 5.0);

  //for(n in 1:N) {
  //  y[n] ~ lognormal(mus[labels[n]], sigmas[labels[n]]);
  //}

  for(l in 1:L) {
    //mus[l] ~ normal(a * log(tmp[l] / min(stress)) + c, sigma);
    mus[l] ~ normal(a[labels[l]] * log(stress[l]) + b[labels[l]] * log(1 / thickness[l]) + c, sigma[labels[l]]);
    //mus[l] ~ normal(a * log(stress[l]) + b * log(1 / thickness[l]) + c, sigma[labels[l]]);
  }
}

generated quantities {
  vector[L] yhat;

  {
    for(l in 1:L) {
      yhat[l] <- lognormal_rng(normal_rng(a[labels[l]] * log(stress[l]) + b[labels[l]] * log(1 / thickness[l]) + c, sigma[labels[l]]), sigmas[l]);
      //yhat[l] <- lognormal_rng(normal_rng(a * log(stress[l]) + b * log(1 / thickness[l]) + c, sigma[labels[l]]), sigmas[l]);
    }
  }
}
"""

sm2 = pystan.StanModel(model_code = model_code)
#%%
