#%%

import scipy
import matplotlib.pyplot as plt
import pystan
import numpy
import seaborn
import pandas
import pickle


#files = [
#    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000A_200MPa_noTertiary.mat', 'creep2000A_200MPa_notertiary', 1.0 / 2.0, 1.0, 2000, 200, 0.0, False],
#    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_150MPa.mat', 'creep_2B_150_good', 1.0 / 2, 1.0, 2000, 150, 0.0, False],
#    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_170MPa.mat', 'creep_2B_170_good', 0.0, 1.0, 2000, 170, 0.0, False],
#    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_190MPa.mat', 'creep_2B_190_full', 0.0, -1.0 / 4.0, 2000, 190, 0.0, False],
#    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_200MPa.mat', 'creep_2B_200_goodyay', 0.0, 1.0, 2000, 200, 0.0, False],
#    ['/home/bbales2/CreepInfo/8000nm-8A/creep_8000A_140MPa.mat', 'creepinfo_8A_140', 1.0 / 2.0, 1.0, 8000, 140, 0.0, False],
#    ['/home/bbales2/CreepInfo/8000nm-8A/creep_8000A_150MPa.mat', 'creepinfo_8A_150', 1.0 / 2.0, 1.0, 8000, 150, 0.0, False],
#    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200A_245MPa.mat', 'creepinfo_200a_245_full', 1.0 / 2.0, -1.0 / 4.0, 200, 245, 0.0, False],
#    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200A_260MPa.mat', 'creepinfo_200a_260', 0.0, 1.0, 200, 260, 0.0, False],
#    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200A_290MPa.mat', 'creepinfo_200a_290_full', 0.0, 1.0 / 4.0, 200, 290, 0.0, False],
#    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200B_275MPa.mat', 'creepinfo_200B_275MPa', 1.0 / 2.0, 1.0, 200, 275, 0.0, False],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65C_320MPa.mat', 'creep65C_320full', 1.0 / 2.0, 1.0, 65, 320, 3.0, False],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65C_350MPa.mat', 'creep65C_350full', 0.0, 1.0 / 2.0, 65, 350, 71.0, False],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_275MPa.mat', 'CreepInfo_275corr', 1.0 / 2.0, 1.0, 65, 275, 4.3, False],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_290MPa.mat', 'CreepInfo_290corr', 0.0, 1.0, 65, 290, 103.1, False],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_320MPa.mat', 'CreepInfo_320corr', 0.0, 1.0, 65, 320, 193.4, False],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_335MPa.mat', 'CreepInfo_335corrected', 0.0, 1.0, 65, 335, 266.1, False],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_230MPa.mat', 'CreepInfo65B_230full', 0.0, 1.0, 65, 230, 447.8, True],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_245MPa.mat', 'CreepInfo65B_245full', 0.0, 1.0, 65, 245, 514.6, True],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_260MPa.mat', 'CreepInfo65B_260full', 0.0, 1.0, 65, 260, 593.1, True],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_290MPa.mat', 'CreepInfo65B_290full', 0.0, 1.0, 65, 290, 641.6, True],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_320MPa.mat', 'CreepInfo65B_320full', 0.0, 1.0, 65, 320, 721.3, True],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_350MPa.mat', 'CreepInfo65B_350full', 0.0, 1.0, 65, 350, 827.5, True],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_370MPa_noTertiary.mat', 'CreepInfo65B_370beforetertiary', 0.0, 1.0, 65, 370, 895.8, True],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65HT_320MPa.mat', 'creep65HT_320full', 1.0 / 2.0, 1.0, 65, 320, 723.0, True],
#    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65HT_350MPa_noTertiary.mat', 'creep65HT_350notertiary', 0.0, 1.0, 65, 350, 813.3, True]
#]

files = [
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000A_200MPa_noTertiary.mat', 'creep2000A_200MPa_notertiary', 1.0 / 2.0, 1.0, 500, 200, 0.0, False],
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_150MPa.mat', 'creep_2B_150_good', 1.0 / 2, 1.0, 500, 150, 0.0, False],
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_170MPa.mat', 'creep_2B_170_good', 0.0, 1.0, 500, 170, 0.0, False],
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_190MPa.mat', 'creep_2B_190_full', 0.0, -1.0 / 4.0, 500, 190, 0.0, False],
    ['/home/bbales2/CreepInfo/2000nm-2B-2A/creep_2000B_200MPa.mat', 'creep_2B_200_goodyay', 0.0, 1.0, 500, 200, 0.0, False],
    ['/home/bbales2/CreepInfo/8000nm-8A/creep_8000A_140MPa.mat', 'creepinfo_8A_140', 1.0 / 2.0, 1.0, 1000, 140, 0.0, False],
    ['/home/bbales2/CreepInfo/8000nm-8A/creep_8000A_150MPa.mat', 'creepinfo_8A_150', 1.0 / 2.0, 1.0, 1000, 150, 0.0, False],
    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200A_245MPa.mat', 'creepinfo_200a_245_full', 1.0 / 2.0, -1.0 / 4.0, 200, 245, 0.0, False],
    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200A_260MPa.mat', 'creepinfo_200a_260', 0.0, 1.0, 200, 260, 0.0, False],
    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200A_290MPa.mat', 'creepinfo_200a_290_full', 0.0, 1.0 / 4.0, 200, 290, 0.0, False],
    ['/home/bbales2/CreepInfo/200nm-200A-200B/creep_200B_275MPa.mat', 'creepinfo_200B_275MPa', 1.0 / 2.0, 1.0, 200, 275, 0.0, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65C_320MPa.mat', 'creep65C_320full', 1.0 / 2.0, 1.0, 65, 320, 3.0, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65C_350MPa.mat', 'creep65C_350full', 0.0, 1.0 / 2.0, 65, 350, 71.0, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_275MPa.mat', 'CreepInfo_275corr', 1.0 / 2.0, 1.0, 65, 275, 4.3, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_290MPa.mat', 'CreepInfo_290corr', 0.0, 1.0, 65, 290, 103.1, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_320MPa.mat', 'CreepInfo_320corr', 0.0, 1.0, 65, 320, 193.4, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/As Rolled/creep_65D_335MPa.mat', 'CreepInfo_335corrected', 0.0, 1.0, 65, 335, 266.1, False],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_230MPa.mat', 'CreepInfo65B_230full', 0.0, 1.0, 65, 230, 447.8, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_245MPa.mat', 'CreepInfo65B_245full', 0.0, 1.0, 65, 245, 514.6, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_260MPa.mat', 'CreepInfo65B_260full', 0.0, 1.0, 65, 260, 593.1, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_290MPa.mat', 'CreepInfo65B_290full', 0.0, 1.0, 65, 290, 641.6, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_320MPa.mat', 'CreepInfo65B_320full', 0.0, 1.0, 65, 320, 721.3, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_350MPa.mat', 'CreepInfo65B_350full', 0.0, 1.0, 65, 350, 827.5, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65B_370MPa_noTertiary.mat', 'CreepInfo65B_370beforetertiary', 0.0, 1.0, 65, 370, 895.8, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65HT_320MPa.mat', 'creep65HT_320full', 1.0 / 2.0, 1.0, 65, 320, 723.0, True],
    ['/home/bbales2/CreepInfo/65nm-65B-65HT-65C-65D/Heat Treatment/creep_65HT_350MPa_noTertiary.mat', 'creep65HT_350notertiary', 0.0, 1.0, 65, 350, 813.3, True]
]

files = sorted(files, key = lambda x : x[-3])

df = pandas.DataFrame(files, columns = ['file', 'variable', 'iminf', 'imaxf', 'thickness', 'stress', 'heat_treatment', 'treated'])

df = df[df['treated'] == False]
#df = df.query('not (thickness == 200 and (stress == 260 or stress == 290))')
df = df.reset_index(drop = True)
#%%
#%%

model_code = """
data {
  int<lower=1> N; // Number of single samples
  vector<lower=0.0>[N] t;
  vector<lower=0.0>[N] y;
  matrix[N, N] y2[N, N];
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
print "{0:10s}, {1:10s}, {2:12s}, {3:10s}, {4}".format("thickness", "stress", "heat treated", "avg minimum strain rate", "std deviation of minimum strain rate")
for idx, row in df.iterrows():
    filename, variable, imin, imax, thickness, stress, _, heat_treated = row.values

    data = scipy.io.loadmat(filename)
    data = data[variable]

    imin = int(len(data) * imin)
    imax = int(len(data) * imax)
    data = data[imin : imax]

    plt.plot(data[:, 0], data[:, 1])
    plt.show()

    if len(data) > 100:
        data = data[sorted(numpy.random.choice(range(len(data)), 100, replace = False))]

    data[:, 0] /= 1e6
    #data[:, 1] = 10.0

    fit = sm.sampling(data = {
      'N' : len(data),
      't' : data[:, 0],
      'y' : data[:, 1]
    })

    sample_idxs = numpy.random.choice(range(2000, 4000), 1000, replace = False)

    slope_samples = 1e-6 * fit.extract()['a'][sample_idxs]
    error_samples = 1e-6 * fit.extract()['sigma'][sample_idxs]

    slopes.append(slope_samples)

    #plt.plot(data[:, 0], data[:, 1], '-*')

    #for i in numpy.random.choice(range(2000, 4000), 50, replace = False):
    #    plt.plot(data[:, 0], fit.extract()['yhat'][i], '--*')
    #plt.show()

    #seaborn.distplot(slope_samples)
    #plt.show()

    print "{0:10.0f}, {1:10.0f}, {2:12s}, {3:10.4e}, {4:10.4e}".format(thickness, stress, str(heat_treated), numpy.mean(slope_samples) / 1.0e6, numpy.mean(error_samples) / 1.0e6)

    #print len(data)numpy.std(slope_samples),
#%%
f = open('/home/bbales2/CreepInfo/dump', 'w')
pickle.dump((slopes, df), f)
f.close()
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
  //real a_mu;
  //real<lower = 0.0> a_sigma;
  //real b_mu;
  //real<lower = 0.0> b_sigma;
  real a;//[T];
  real b;//[T];
  real c;
}

model {
  real tmp[L];
  //sigmas ~ cauchy(0.0, 10.0);
  //sigma ~ cauchy(0.0, 10.0);
  //a_sigma ~ normal(0.0, 10.01);
  //b_sigma ~ normal(0.0, 10.01);
  //a ~ normal(a_mu, a_sigma);
  //b ~ normal(b_mu, b_sigma);
  sigma ~ normal(0.0, 5.0);
  //b ~ cauchy(0.0, 10.0);
  //c ~ normal(-2.2, 5.0);

  //for(n in 1:N) {
  //  y[n] ~ lognormal(mus[labels[n]], sigmas[labels[n]]);
  //}

  for(l in 1:L) {
    //mus[l] ~ normal(a * log(tmp[l] / min(stress)) + c, sigma);
    //mus[l] ~ normal(a[labels[l]] * log(stress[l]) + b[labels[l]] * log(1 / thickness[l]) + c, sigma[labels[l]]);
    mus[l] ~ normal(a * log(stress[l]) + b * log(1 / thickness[l]) + c, sigma[labels[l]]);
  }// / stress0[labels[l]]
}

generated quantities {
  vector[L] yhat;//log

  {
    for(l in 1:L) {
      //yhat[l] <- lognormal_rng(normal_rng(a[labels[l]] * log(stress[l]) + b[labels[l]] * log(1 / thickness[l]) + c, sigma[labels[l]]), sigmas[l]);
      yhat[l] <- lognormal_rng(normal_rng(a * log(stress[l]) + b * log(1 / thickness[l]) + c, sigma[labels[l]]), sigmas[l]);
    }// / stress0[labels[l]]
  }
}
"""

sm2 = pystan.StanModel(model_code = model_code)
#%%

thicknessLabels = dict([(v, i + 1) for i, v in enumerate(sorted(list(set(df['thickness']))))])
#%%

ys = []
thicknesses = []
stresses = []
labels = []
mus = []
sigmas = []
for idx, row in df.iterrows():
    if row['treated']:
        continue

    slope_samples = slopes[idx]

    thicknesses.append(row['thickness'])#
    stresses.append(row['stress'])
    ys.extend(slope_samples)
    #labels.extend([idx + 1] * len(slope_samples))
    labels.append(thicknessLabels[row['thickness']])
    mus.append(numpy.log(slope_samples).mean())
    sigmas.append(numpy.log(slope_samples).std())

stress0 = []
for thickness in sorted(thicknessLabels.keys()):
    stress0.append(min(df[df['thickness'] == thickness]['stress']))
#%%
fit2 = sm2.sampling(data = {
    #'N' : len(ys),
    'L' : len(stresses),
    'T' : len(thicknessLabels),
    #'y' : ys,
    'stress0' : stress0,
    'thickness' : thicknesses,
    'stress' : stresses,
    'labels' : labels,
    'mus' : mus,
    'sigmas' : sigmas
})

print fit2

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

seaborn.boxplot(x = 'stresses', y = 'log_slopes', hue = 'generated', data = df2, linewidth = 0.5, showfliers = False)
plt.gcf().set_size_inches((20, 14))
plt.show()
#%%
import seaborn
import pandas
import matplotlib.pyplot as plt

stuff = { 'c' : r['c'][-500:] }
for i in range(4):
    stuff['a{0}'.format(i)] = r['a'][-500:, i]
    stuff['b{0}'.format(i)] = r['b'][-500:, i]

df3 = pandas.DataFrame(stuff)#{'a' : r['a'][-200:, 0], 'b' : r['b'][-200:, 0], 'c' : r['c'][-200:], 'sigma' : r['b'][-200:, 2]}

seaborn.pairplot(df3)
plt.gcf().set_size_inches((12, 8))
plt.show()
#%%
import seaborn
import pandas
import matplotlib.pyplot as plt

df3 = pandas.DataFrame({'a' : r['a'][-1000:], 'b' : r['b'][-1000:], 'c' : r['c'][-1000:]})#, 'sigma' : r['b'][-500:, 2]

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