import pickle

import numpy as np
import simplejson as json
from tabulate import tabulate


class REPL(dict):

    def __init__(self, items):
        self.items = items

    def __getitem__(self,i):
        try:
            return self.items[i]
        except KeyError as e:
            return i

with open('fit_parameters/shutdown_model_all_provinces_after_feb_12.p','rb') as f:
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
    tabledata.append([titlemap[province], N, Q, R0eff, k0, I0f, ])

headers = [
             'Province', 
             '$N/10^6$',
             '$Q$', 
             '$R_{0,\mathrm{eff}}$', 
             '$\kappa_{0}$ [$\mathrm{d}^{-1}$]',
             '$I_0/X_0$', 
          ]
floatfmt = (
               "",
               ".1f", 
               ".2f", 
               ".2f",
               ".3f",
               ".2f",
           ) 
latex = str(tabulate(tabledata,headers=headers,tablefmt="latex_raw",floatfmt=floatfmt))
print(tabulate(tabledata,headers=headers,floatfmt=floatfmt))


with open('fit_parameters/shutdown_parameters_after_Feb_12.tex','w') as f:
    f.write(latex)

vals = np.array(tabledata)[:,1:]
vals = np.array(vals,dtype=float)
