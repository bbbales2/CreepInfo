#%%

import scipy
import matplotlib.pyplot as plt
import pystan
import numpy
import seaborn
import pandas
import pickle

f = open('/home/bbales2/CreepInfo/dump')
slopes, df = pickle.load(f)
f.close()

#%%

model_code = """
data {
  //int<lower=1> N; // Number of single samples
  int<lower=1> L;
  vector<lower=0.0>[L] stress;
  vector[L] mus;
  vector<lower=0.0>[L] sigmas;
}

parameters {
  real<lower=0.0> sigma;

  real a;
  real b;
}

model {
  sigma ~ cauchy(0.0, 10.0);

  for(l in 1:L) {
    mus[l] ~ normal(a * log(stress[l] / min(stress)) + b, sigma);
  }
}

//generated quantities {
//  vector[L] yhat;//log
//
//  for(l in 1:L) {
//    yhat[l] <- lognormal_rng(normal_rng(a * log(stress[l] / min(stress)) + b, sigma), sigmas[l]);
//  }
//}
"""

sm3 = pystan.StanModel(model_code = model_code)
#%%

ts = sorted(list(set(df['thickness'])))

thicknessLabels = dict([(v, i + 1) for i, v in enumerate(sorted(list(set(df['thickness']))))])
#%%
P = numpy.random.random((5, 5))

v = numpy.random.random(5)

for i in range(5):
    P[:, i] /= sum(P[:, i])

v /= sum(v)
#%%

for t in ts:
    stresses = []
    mus = []
    sigmas = []
    for idx, row in df[df['thickness'] == t].iterrows():
        slope_samples = slopes[idx]

        stresses.append(row['stress'])
        mus.append(numpy.log(slope_samples).mean())
        sigmas.append(numpy.log(slope_samples).std())

    fit2 = sm3.sampling(data = {
        'L' : len(stresses),
        'stress' : stresses,
        'mus' : mus,
        'sigmas' : sigmas
    })

    print 'thickness = {0}'.format(t)
    print fit2

#%%

r = fit2.extract()
#%%
for slope_samples, generated in zip(slopes, r['yhat'].transpose()):
    seaborn.distplot(slope_samples, norm_hist = True)
    seaborn.distplot(generated[-200:], norm_hist = True)
    print slope_samples.mean(), generated[-200:].mean()
    plt.legend(['data', 'generated'])
    plt.show()
#%%
plt.plot(mus, 'r*')
for i in range(1):
    idx = numpy.random.randint(3800, 4000)#3300
    plt.plot(numpy.log(r['yhat'][idx]), 'b*')
    print r['a'][idx], r['b'][idx], r['c'][idx]
plt.show()
#%%
thickness_ = []
stresses_ = []
slopes_ = []
generated_ = []

for idx, row in df.sort_values('thickness', ascending = False).iterrows():
    if row['treated']:
        continue

    thickness_.extend([row['thickness']] * len(slopes[idx]))
    stresses_.extend([row['stress']] * len(slopes[idx]))
    slopes_.extend(numpy.log(slopes[idx]))
    generated_.extend(['Measured, h = {0}'.format(row['thickness'])] * len(slopes[idx]))

    thickness_.extend([row['thickness']] * len(slopes[idx]))
    stresses_.extend([row['stress']] * len(slopes[idx]))
    slopes_.extend(numpy.log(r['yhat'][-len(slopes[idx]):, idx]))
    generated_.extend(['Generated, h = {0}'.format(row['thickness'])] * len(slopes[idx]))

df2 = pandas.DataFrame(data = { 'thickness' : thickness_, 'stresses' : stresses_, 'log_slopes' : slopes_, 'generated' : generated_ })

#for thickness in sorted(set(df2['thickness'])):
#    df3 = df2[df2['thickness'] == thickness]

seaborn.boxplot(x = 'stresses', y = 'log_slopes', hue = 'generated', data = df2, linewidth = 0.5)
plt.gcf().set_size_inches((20, 14))
plt.show()
#%%
import seaborn
import pandas
import matplotlib.pyplot as plt

df3 = pandas.DataFrame({'a' : r['a'][-200:], 'b' : r['b'][-200:], 'c' : r['c'][-200:], 'sigma' : r['sigma'][-200:]})

seaborn.pairplot(df3)
plt.gcf().set_size_inches((12, 8))
plt.show()
#%%
for name, d in [('a', r['a']), ('b', r['b']), ('c', r['c']), ('sigma', r['sigma'])]:
    seaborn.distplot(d[-200:], kde = False, fit = scipy.stats.norm)
    plt.title("Dist. {0} w/ mean {1:0.4f} and std. {2:0.4f}".format(name, numpy.mean(d[-200:]), numpy.std(d[-200:])))
    plt.gcf().set_size_inches((5, 4))
    plt.show()
#%%
#data2 = scipy.io.loadmat('jackie/CreepInfo_90MPa_corr_NoBlip.mat')['CreepInfo_090corr_NoBlip'][::20]
#data3 = scipy.io.loadmat('jackie/CreepInfo_105MPa_corr.mat')['CreepInfo_105corr'][::20]

plt.plot(data1[len(data1) / 2:-len(data1) / 4, 0], data1[len(data1) / 2:-len(data1) / 4, 1], '*')
plt.title('From 200nm-200A-200B/creep_200B_245MPa.mat')
plt.show()