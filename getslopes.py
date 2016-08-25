#%%

import scipy
import matplotlib.pyplot as plt
import pystan
import numpy
import seaborn

#data1 = scipy.io.loadmat()[]
files = [
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000A_200MPa_noTertiary.mat', 'creep2000A_200MPa_notertiary', 1.0 / 2.0, 1.0, 2000, 200, 0.0],
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_150MPa.mat', 'creep_2B_150_good', 1.0 / 2, 1.0, 2000, 150, 0.0],
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_170MPa.mat', 'creep_2B_170_good', 0.0, 1.0, 2000, 170, 0.0],
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_190MPa.mat', 'creep_2B_190_full', 0.0, -1.0 / 4.0, 2000, 190, 0.0],
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_200MPa.mat', 'creep_2B_200_goodyay', 0.0, 1.0, 2000, 200, 0.0],
    ['/home/bbales2/CreepInfo/8000nm-8A/creep_8000A_140MPa.mat', 'creepinfo_8A_140', 1.0 / 2.0, 1.0, 8000, 140, 0.0],
    ['/home/bbales2/CreepInfo/8000nm-8A/creep_8000A_150MPa.mat', 'creepinfo_8A_150', 1.0 / 2.0, 1.0, 8000, 150, 0.0],
    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200A_245MPa.mat', 'creepinfo_200a_245_full', 1.0 / 2.0, -1.0 / 4.0, 200, 245, 0.0],
    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200A_260MPa.mat', 'creepinfo_200a_260', 0.0, 1.0, 200, 260, 0.0],
    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200A_290MPa.mat', 'creepinfo_200a_290_full', 0.0, 1.0 / 4.0, 200, 290, 0.0],
    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200B_275MPa.mat', 'creepinfo_200B_275MPa', 1.0 / 2.0, 1.0, 200, 275, 0.0],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65C_320MPa.mat', 'creep65C_320full', 1.0 / 2.0, 1.0, 65, 320, 3.0],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65C_350MPa.mat', 'creep65C_350full', 0.0, 1.0 / 2.0, 65, 350, 71.0],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_275MPa.mat', 'CreepInfo_275corr', 1.0 / 2.0, 1.0, 65, 275, 4.3],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_290MPa.mat', 'CreepInfo_290corr', 0.0, 1.0, 65, 290, 103.1],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_320MPa.mat', 'CreepInfo_320corr', 0.0, 1.0, 65, 320, 193.4],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_335MPa.mat', 'CreepInfo_335corrected', 0.0, 1.0, 65, 335, 266.1],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_230MPa.mat', 'CreepInfo65B_230full', 0.0, 1.0, 65, 230, 447.8],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_245MPa.mat', 'CreepInfo65B_245full', 0.0, 1.0, 65, 245, 514.6],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_260MPa.mat', 'CreepInfo65B_260full', 0.0, 1.0, 65, 260, 593.1],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_290MPa.mat', 'CreepInfo65B_290full', 0.0, 1.0, 65, 290, 641.6],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_320MPa.mat', 'CreepInfo65B_320full', 0.0, 1.0, 65, 320, 721.3],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_350MPa.mat', 'CreepInfo65B_350full', 0.0, 1.0, 65, 350, 827.5],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_370MPa_noTertiary.mat', 'CreepInfo65B_370beforetertiary', 0.0, 1.0, 65, 370, 895.8],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65HT_320MPa.mat', 'creep65HT_320full', 1.0 / 2.0, 1.0, 65, 320, 723.0],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65HT_350MPa_noTertiary.mat', 'creep65HT_350notertiary', 0.0, 1.0, 65, 350, 813.3]
]
#%%

model_code = """
data {
  int<lower=1> N; // Number of single samples
  vector<lower=0.0>[N] t;
  vector<lower=0.0>[N] y;
}

parameters {
  real<lower=0.0> sigma;
  real a;
  real b;
}

model {
  y ~ normal(a * t + b, sigma);
}

generated quantities {
  vector[N] yhat;

  for(n in 1:N) {
    yhat[n] <- normal_rng(a * t[n] + b, sigma);
  }
}
"""

sm = pystan.StanModel(model_code = model_code)

#%%


#%%
a
seaborn.distplot(a['a'])
#%%
slopes = []
for filename, variable, imin, imax, _, _, _ in files[:]:
    data = scipy.io.loadmat(filename)
    data = data[variable]
    imin = int(len(data) * imin)
    imax = int(len(data) * imax)
    data = data[imin : imax]

    if len(data) > 100:
        data = data[sorted(numpy.random.choice(range(len(data)), 100, replace = False))]

    data[:, 0] /= 1e6
    #data[:, 1] = 10.0

    fit = sm.sampling(data = {
      'N' : len(data),
      't' : data[:, 0],
      'y' : data[:, 1]
    })

    slope_samples = fit.extract()['a'][numpy.random.choice(range(2000, 4000), 100, replace = False)]

    slopes.append(slope_samples)

    plt.plot(data[:, 0], data[:, 1], '-*')

    for i in numpy.random.choice(range(2000, 4000), 5, replace = False):
        plt.plot(data[:, 0], fit.extract()['yhat'][i], '--*')
    plt.show()

    seaborn.distplot(slope_samples)
    plt.show()

    print numpy.std(slope_samples), numpy.mean(slope_samples)

    print len(data)
#%%

model_code = """
data {
  int<lower=1> N; // Number of single samples
  int<lower=1> L;
  int<lower=1> T;
  int<lower=1> labels[N];
  vector<lower=0.0>[L] stress;
  //vector<lower=0.0>[L] thickness;
  int<lower=1> thickness[L];
  vector<lower=0.0>[N] y;
}

parameters {
  real mus[L];
  real<lower=0.0> sigmas[L];

  real<lower=0.0> sigma;
  real a;
  real b;
  real c[L];
  //real d;
}

model {
  for(n in 1:N) {
    y[n] ~ normal(mus[labels[n]], sigmas[labels[n]]);
  }

  //for(l in 1:L) {
  //  mus[l] ~ lognormal(a * log(stress[l] + c[thickness[l]]) + b, sigma);
  //}
}

//generated quantities {
//  vector[N] yhat;
//
//  for(n in 1:N) {
//    yhat[n] <- lognormal_rng(a * log(stress[n] + c / sqrt(thickness[n]) + d) + b, sigma);
//  }
//}
"""

sm2 = pystan.StanModel(model_code = model_code)
#%%
_, _, _, _, thicknesses, _, _ = zip(*files)

thicknessLabels = dict([(v, i + 1) for i, v in enumerate(sorted(list(set(thicknesses))))])
#%%

ys = []
thicknesses = []
stresses = []
labels = []
for e, (slope_samples, (_, _, _, _, thickness, stress, heat_treatment)) in enumerate(zip(slopes, files)):
    if heat_treatment == 447.8:
        break

    thicknesses.append(thicknessLabels[thickness])
    stresses.append(stress)
    ys.extend(slope_samples)
    labels.extend([e + 1] * len(slope_samples))
#%%
sm2.sampling(data = {
    'N' : len(ys),
    'L' : len(stresses),
    'T' : len(thicknesses),
    'y' : ys,
    'thickness' : thicknesses,
    'stress' : stresses,
    'labels' : labels
})
#%%
#data2 = scipy.io.loadmat('jackie/CreepInfo_90MPa_corr_NoBlip.mat')['CreepInfo_090corr_NoBlip'][::20]
#data3 = scipy.io.loadmat('jackie/CreepInfo_105MPa_corr.mat')['CreepInfo_105corr'][::20]

plt.plot(data1[len(data1) / 2:-len(data1) / 4, 0], data1[len(data1) / 2:-len(data1) / 4, 1], '*')
plt.title('From 200nm-200A-200B/creep_200B_245MPa.mat')
plt.show()