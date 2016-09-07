#%%

import scipy
import matplotlib.pyplot as plt
import pystan
import numpy
import seaborn
import pandas
import pickle
#['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_230MPa.mat', 'CreepInfo65B_230full', 0.0, 1.0, 65, 230, 447.8, True, 3],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_245MPa.mat', 'CreepInfo65B_245full', 0.0, 1.0, 65, 245, 514.6, True, 3],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_260MPa.mat', 'CreepInfo65B_260full', 0.0, 1.0, 65, 260, 593.1, True, 3],

#data1 = scipy.io.loadmat()[]
files = [
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65C_320MPa.mat', 'creep65C_320full', 1.0 / 2.0, 1.0, 65, 320, 3.0, False, 1],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65C_350MPa.mat', 'creep65C_350full', 0.0, 1.0 / 2.0, 65, 350, 71.0, False, 1],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_275MPa.mat', 'CreepInfo_275corr', 1.0 / 2.0, 1.0, 65, 275, 4.3, False, 2],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_290MPa.mat', 'CreepInfo_290corr', 0.0, 1.0, 65, 290, 103.1, False, 2],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_320MPa.mat', 'CreepInfo_320corr', 0.0, 1.0, 65, 320, 193.4, False, 2],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_335MPa.mat', 'CreepInfo_335corrected', 0.0, 1.0, 65, 335, 266.1, False, 2],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_290MPa.mat', 'CreepInfo65B_290full', 0.0, 1.0, 65, 290, 641.6, True, 3],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_320MPa.mat', 'CreepInfo65B_320full', 0.0, 1.0, 65, 320, 721.3, True, 3],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_350MPa.mat', 'CreepInfo65B_350full', 0.0, 1.0, 65, 350, 827.5, True, 3],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_370MPa_noTertiary.mat', 'CreepInfo65B_370beforetertiary', 0.0, 1.0, 65, 370, 895.8, True, 3],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65HT_320MPa.mat', 'creep65HT_320full', 1.0 / 2.0, 1.0, 65, 320, 723.0, True, 4],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65HT_350MPa_noTertiary.mat', 'creep65HT_350notertiary', 0.0, 1.0, 65, 350, 813.3, True, 4]
]

files = sorted(files, key = lambda x : x[-3])

df = pandas.DataFrame(files, columns = ['file', 'variable', 'iminf', 'imaxf', 'thickness', 'stress', 'heat_treatment', 'treated', 'labels'])

#df = df[df['treated'] == True]
#df = df.query('stress > 290')#not (thickness == 200 and (stress == 260 or stress == 290))')
df = df.reset_index(drop = True)
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

slopes = []
print "{0:10s}, {1:10s}, {2:12s}, {3}".format("thickness", "stress", "heat treated", "avg minimum strain rate")
for idx, row in df.iterrows():
    filename, variable, imin, imax, thickness, stress, _, heat_treated, label = row.values

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

    slope_samples = fit.extract()['a'][numpy.random.choice(range(2000, 4000), 1000, replace = False)]

    new_slopes = []
    for slope in slope_samples:
        if slope > 0.0:
            new_slopes.append(slope)

    slopes.append(new_slopes)

    #plt.plot(data[:, 0], data[:, 1], '-*')

    #for i in numpy.random.choice(range(2000, 4000), 50, replace = False):
    #    plt.plot(data[:, 0], fit.extract()['yhat'][i], '--*')
    #plt.show()

    #seaborn.distplot(slope_samples)
    #plt.show()

    print "{0:10.0f}, {1:10.0f}, {2:12s}, {3:10.4e}".format(thickness, stress, str(heat_treated), numpy.mean(slope_samples) / 1.0e6)

    #print len(data)numpy.std(slope_samples),

#%%

model_code = """
data {
  int<lower=1> N; // Number of single samples
  int<lower=1> L;
  int<lower=1> labels[N];
  vector<lower=0.0>[N] stress;
  vector<lower=0.0>[N] ht;
  vector<lower=0.0>[L] stress0;
  vector[N] mus;
  vector<lower=0.0>[N] sigmas;
  real c_sigma_sigma;
  real a_sigma_sigma;
  real b_sigma_sigma;
}

parameters {
  real<lower=0.0> sigma[L];

  real a_mu;
  real<lower = 0.0> a_sigma;

  real b_mu;
  real<lower = 0.0> b_sigma;

  real c_mu;
  real<lower = 0.0> c_sigma;

  real a[L];
  real b[L];
  real c[L];
}

model {
  real tmp[L];

  a_mu ~ normal(11.5, 0.1);

  a_sigma ~ normal(0.0, a_sigma_sigma);
  b_sigma ~ normal(0.0, b_sigma_sigma);
  c_sigma ~ normal(0.0, c_sigma_sigma);
  a ~ normal(a_mu, a_sigma);
  b ~ normal(b_mu, b_sigma);
  c ~ normal(c_mu, c_sigma);

  sigma ~ normal(0.0, 5.0);

  for(n in 1:N) {
    mus[n] ~ normal(a[labels[n]] * log(stress[n] / stress0[labels[n]]) +
                    b[labels[n]] * log(ht[n]) +
                    c[labels[n]], sigma[labels[n]]);
  }
}

generated quantities {
  vector[N] yhat;//log

  for(n in 1:N) {
    yhat[n] <- lognormal_rng(normal_rng(a[labels[n]] * log(stress[n] / stress0[labels[n]]) +
                                        b[labels[n]] * log(ht[n]) +
                                        c[labels[n]], sigma[labels[n]]), sigmas[n]);
  }
}
"""

sm2 = pystan.StanModel(model_code = model_code)
#%%
#%%

model_code = """
data {
  int<lower=1> N; // Number of single samples
  int<lower=1> L;
  int<lower=1> labels[N];
  vector<lower=0.0>[N] stress;
  vector<lower=0.0>[N] ht;
  vector<lower=0.0>[L] stress0;
  vector[N] mus;
  vector<lower=0.0>[N] sigmas;
  real c_sigma_sigma;
  real a_sigma_sigma;
  real b_sigma_sigma;
}

parameters {
  real<lower=0.0> sigma;

  real a;
  real b;
  real c;
}

model {
  for(n in 1:N) {
    stress[n] ~ lognormal(a * mus[n] +
                        b * ht[n] +
                          c, sigma);
  }
}

generated quantities {
  vector[N] yhat;//log

  for(n in 1:N) {
    yhat[n] <- lognormal_rng(a * mus[n] +
                             b * ht[n] +
                             c, sigma);
  }
}
"""

sm2 = pystan.StanModel(model_code = model_code)
#%%

thicknessLabels = dict([(v, i + 1) for i, v in enumerate(sorted(list(set(df['thickness']))))])
#%%

mus = []
sigmas = []
for idx, row in df.iterrows():
    slope_samples = slopes[idx]

    mus.append(numpy.log(0.01 / numpy.array(slope_samples)).mean())
    sigmas.append(numpy.log(0.01 / numpy.array(slope_samples)).std())

stress0 = []
for label in sorted(list(set(df['labels']))):
    stress0.append(min(df[df['labels'] == label]['stress']))

#%%

fit2 = sm2.sampling(data = {
    'N' : len(df),
    'L' : 4,
    'labels' : df['labels'],
    'stress' : df['stress'],
    'ht' : df['heat_treatment'],
    'stress0' : stress0,
    'mus' : mus,
    'sigmas' : sigmas,
    'a_sigma_sigma' : 10.01,
    'b_sigma_sigma' : 10.01,
    'c_sigma_sigma' : 0.01
})

print fit2

r = fit2.extract()

#%%
print 'P(b[0] < 0) = {0}'.format(sum(r['b'][-2000:, 0] < 0.0) / 2000.0)
print 'P(b[1] < 0) = {0}'.format(sum(r['b'][-2000:, 1] < 0.0) / 2000.0)
print 'P(b[2] < 0) = {0}'.format(sum(r['b'][-2000:, 2] < 0.0) / 2000.0)
print 'P(b[3] < 0) = {0}'.format(sum(r['b'][-2000:, 3] < 0.0) / 2000.0)

#%%

labels = ['65C', '65D', '65B', '65HT']


for idx, row in df.iterrows():
    samples = r['yhat'][-2000:, idx]
    print numpy.mean(samples), numpy.std(samples)
    print labels[row['labels'] - 1], '{0}Mpa'.format(row['stress'])#, '{0} m/s'.format(numpy.exp(mus[idx]) / 1e6)
    print '----'

#%%
plt.plot(r['yhat'][-2000:].mean(axis = 0), mus, 'bo')
plt.plot(df['stress'], mus, 'r*')
ax = plt.gca()
for idx, row in df.iterrows():
    ax.annotate('{0}, {1}hr'.format(labels[row['labels'] - 1], row['heat_treatment']), xy = (df['stress'][idx], mus[idx]))
plt.gcf().set_size_inches((12, 8))
#%%

thickness_ = []
stresses_ = []
slopes_ = []
generated_ = []
for idx, row in df.sort_values('labels', ascending = True).iterrows():
    #if row['labels'] != 2:
    #    continue

    thickness_.extend([row['labels']] * len(slopes[idx]))
    stresses_.extend([row['stress']] * len(slopes[idx]))
    slopes_.extend(numpy.log(slopes[idx]))
    generated_.extend(['Measured, h = {0}'.format(labels[row['labels'] - 1])] * len(slopes[idx]))

    samples = 500

    thickness_.extend([row['labels']] * samples)
    stresses_.extend([row['stress']] * samples)
    slopes_.extend(numpy.log(r['yhat'][-samples:, idx]))
    generated_.extend(['Generated, h = {0}'.format(labels[row['labels'] - 1])] * samples)

df2 = pandas.DataFrame(data = { 'thickness' : thickness_, 'stresses' : stresses_, 'log_slopes' : slopes_, 'generated' : generated_ })

#for thickness in sorted(set(df2['thickness'])):
#    df3 = df2[df2['thickness'] == thickness]

seaborn.boxplot(x = 'stresses', y = 'log_slopes', hue = 'generated', data = df2, linewidth = 0.5, showfliers = False)
plt.gcf().set_size_inches((20, 14))
plt.show()

#%%

plt.plot(r['a'][-500:, 1], r['b'][-500:, 1], 'o')
plt.show()