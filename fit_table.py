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

with open('fit_parameters/confirmed_cases_500.p','rb') as f:
    fit_parameters = pickle.load(f)

with open('fit_parameters/hubei_china.p','rb') as f:
    fit_parameters.update(pickle.load(f))

with open('data/all_confirmed_cases_with_population.json','r') as f:
    data = json.load(f)

tuplelist = [ (p, d)  for p, d in data.items()\
                               if max(d['cases']) >= 20\
                               and len(np.where(np.logical_and(np.array(d['times'])<=12,
                                                               np.array(d['cases'])>0))[0])\
                                   >= 8
                             ]

tuplelist = sorted([ t for t in tuplelist ],key=lambda x: -max(x[1]['cases']))
titlemap = REPL({'mainland_china':'all w/o Hubei'})

tabledata = []
for province, pdata in tuplelist[:10]:
    p = fit_parameters[province]
    k = p['kappa'].value
    k0 = p['kappa0'].value
    b = p['rho'].value
    I0f = p['I0_factor'].value

    N = pdata['population'] / 1e6

    Q = (k+k0)/(b+k+k0)
    P = k0/(k+k0)
    tabledata.append([titlemap[province], N, Q, P, I0f])
print(tabulate(tabledata,headers=['Province', '$N$','$Q$', '$P$', '$I_0/Q_0$'],floatfmt=("",".1f", ".2f", ".2f",".2f",".2f") ))
latex = str(tabulate(tabledata,headers=['Province', '$N/10^6$','$Q$', '$P$', '$I_0/Q_0$'],tablefmt="latex_raw",floatfmt=("",".1f", ".2f", ".2f",".2f",".2f") ))


with open('fit_parameters/parameters.tex','w') as f:
    f.write(latex)

vals = np.array(tabledata)[:,1:]
vals = np.array(vals,dtype=float)
print(vals.mean(axis=0))
print(vals.std(axis=0))
