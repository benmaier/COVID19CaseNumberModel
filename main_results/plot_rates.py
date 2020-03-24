import pickle

import numpy as np
import simplejson as json
from tabulate import tabulate

from bfmplot import pl

from bfmplot import brewer_qualitative, simple_cycler, markers
colors = simple_cycler(brewer_qualitative)


class REPL(dict):

    def __init__(self, items):
        self.items = items

    def __getitem__(self,i):
        try:
            return self.items[i]
        except KeyError as e:
            return i

fn = 'all_provinces_after_feb_12.p'
fn = 'different_initial_conds_after_feb_12.p'

with open('fit_parameters/'+fn,'rb') as f:
    fit_parameters = pickle.load(f)
#with open('fit_parameters/confirmed_cases_500.p','rb') as f:
#    fit_parameters = pickle.load(f)
#
#with open('fit_parameters/hubei_china.p','rb') as f:
#    fit_parameters.update(pickle.load(f))

with open('../data/all_confirmed_cases_with_population.json','r') as f:
    data = json.load(f)

tuplelist = [ (p, d)  for p, d in data.items()\
                               if max(d['cases']) >= 20\
                               and len(np.where(np.logical_and(np.array(d['times'])<=12,
                                                               np.array(d['cases'])>0))[0])\
                                   >= 8
                             ]

tuplelist = sorted([ t for t in tuplelist ],key=lambda x: -max(x[1]['cases']))
titlemap = REPL({'mainland_china':'All exc. Hubei'})

tabledata = []
ks = []
k0s = []
for province, pdata in tuplelist:
    p = fit_parameters[province]
    k = p['kappa'].value
    k0 = p['kappa0'].value
    b = p['rho'].value
    a = p['eta'].value
    I0f = p['I0_factor'].value

    N = pdata['population'] / 1e6

    print(p['kappa0'],p['kappa'])
    Q = (k+k0)/(b+k+k0)
    P = k0/(k+k0)
    R0eff = a / (b+k+k0)
    tabledata.append([titlemap[province], N, Q, P, R0eff, k, k0, I0f, ])
    ks.append(k)
    k0s.append(k0)

for i, (_k0, _ks) in enumerate(zip(k0s, ks)):
    pl.plot([_k0],[_ks],marker=markers[i],color=colors[i],markersize=7)

pl.xlim([0,0.14])
pl.ylim([-.05,0.55])
pl.gcf().savefig('model_fit_figures/'+fn+'.png',dpi=300)

pl.show()

headers = [
             'Province', 
             '$N/10^6$','$Q$', 
             '$P$', 
             '$R_{0,\mathrm{eff}}$', 
             '$\kappa$ [$\mathrm{d}^{-1}$]',
             '$\kappa_{0}$ [$\mathrm{d}^{-1}$]',
             '$I_0/X_0$', 
          ]
floatfmt = (
               "",
               ".1f", 
               ".2f", 
               ".2f",
               ".2f",
               ".3f",
               ".3f",
               ".2f",
           ) 
latex = str(tabulate(tabledata,headers=headers,tablefmt="latex_raw",floatfmt=floatfmt))
print(latex)


vals = np.array(tabledata)[:,1:]
vals = np.array(vals,dtype=float)
