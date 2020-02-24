import sys
sys.path.insert(0,'..')

import numpy as np
from scipy.integrate import ode
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters
import json
from tqdm import tqdm
from bfmplot import pl
from bfmplot import brewer_qualitative, simple_cycler, markers
from SIRX import SIRXConfirmedModel
import pickle

import bfmplot as bp

model = SIRXConfirmedModel()

colors = simple_cycler(brewer_qualitative)

class REPL(dict):

    def __init__(self, items):
        self.items = items

    def __getitem__(self,i):
        try:
            return self.items[i]
        except KeyError as e:
            return i

with open('../data/all_confirmed_csse_cases_with_population.json','r') as f:
    data = json.load(f)
with open('../data/all_confirmed_cases_with_population.json','r') as f:
    dataold = json.load(f)

tuplelist = [ (p, d)  for p, d in data.items()\
                               if max(d['cases']) >= 20\
                               and len(np.where(np.logical_and(np.array(d['times'])<=12,
                                                               np.array(d['cases'])>0))[0])\
                                   >= 8
                             ]

tuplelist = sorted([ t for t in tuplelist ],key=lambda x: -max(x[1]['cases']))

loaded_fits = len(sys.argv) > 1
if loaded_fits:
    pickle_filename = sys.argv[1]


n_fits = len(tuplelist)
n_col = int(np.ceil(np.sqrt(n_fits)))
n_row = n_col
n_col = n_col-1
fig, ax = pl.subplots(n_row,n_col,figsize=(10,10))
ax = ax.flatten()

titlemap = REPL({'mainland_china':'All exc. Hubei'})

if loaded_fits:
    with open(pickle_filename,'rb') as f:
        fit_parameters = pickle.load(f)
else:
    fit_parameters = {}

letter = "abcdefg"
roman = [ "i", "ii", "iii", "iv", "v", "vi"]


i = -1
for province, pdata in tqdm(tuplelist[:]):
    i += 1

    t = np.array(pdata['times'])
    cases = np.array(pdata['cases'])
    dates = np.array(pdata['dates'],np.datetime64)

    told = np.array(dataold[province]['times'])
    casesold = np.array(dataold[province]['cases'])
    datesold = np.array(dataold[province]['dates'],np.datetime64)

    if province == 'Hubei':
        ndx = dates < np.datetime64("2020-02-13")
        t = t[ndx]
        cases = cases[ndx]
        dates = dates[ndx]

    if max(cases) <= 20:
        continue

    i0 = np.where(cases>0)[0][0]
    print(i0)
    t = t[i0:]
    cases = cases[i0:]

    if len(t) < 6:
        continue

    print(pdata['population'])

    if loaded_fits: 
        params = fit_parameters[province]
    else:
        out = model.fit(t,cases,maxfev=100000,N=pdata['population']
                )
        params = out.params
        fit_parameters[province] = params
    print(params)
    N = params['N']

    pl.sca(ax[i])

    tt = np.logspace(np.log(t[0]), np.log(60), base=np.exp(1))
    result = model.SIRX(tt, cases[0], 
                        params['eta'],
                        params['rho'],
                        params['kappa'],
                        params['kappa0'],
                        N,
                        params['I0_factor'],
                        )
    X = result[2,:]*N
    I = result[1,:]*N
#S = result[0,:]*N

    pl.plot(t, cases,marker=markers[i],c=colors[i],label='data',mfc='None')
    #pl.plot(told, casesold,marker=markers[i],c=colors[i],label='data',mfc='None')
    pl.plot(tt, X,c='k',label='$Q_I$ (detected and quarantined)')
    pl.plot(tt, I,'--',c=colors[2],lw=1.5,label='$I$ (undected infected)')

#pl.plot(tt, S,label='model')
    
    _c = i % n_col
    _r = i // n_col
    if _r == n_row-1:
        pl.xlabel('days since Jan. 20th')        
    if _c == 0:
        pl.ylabel('confirmed')
    #pl.title(titlemap[province])
    ax[i].text(0.03,0.97,
            "{}.{}".format(letter[_r], roman[_c]),
            transform=ax[i].transAxes,
            ha='left',
            va='top',
            fontweight='bold',
            fontsize=12,
            bbox={'facecolor':'w','edgecolor':'w','pad':0}
            )
    ax[i].text(0.03,0.8,
            titlemap[province],
            transform=ax[i].transAxes,
            ha='left',
            va='top',
            bbox={'facecolor':'w','edgecolor':'w','pad':0}
            )
    ax[i].text(0.97,0.5,
            r"$Q=%4.2f$" %((params['kappa'].value+params['kappa0'].value)/(params['rho'].value+params['kappa'].value+params['kappa0'].value)),
            transform=ax[i].transAxes,
            ha='right',
            va='bottom',
            bbox={'facecolor':'w','edgecolor':'w','pad':0}
            )
    ax[i].text(0.97,0.38,
            r"$P=%4.2f$" %((params['kappa0'].value)/(params['kappa'].value+params['kappa0'].value)),
            transform=ax[i].transAxes,
            ha='right',
            va='bottom',
            bbox={'facecolor':'w','edgecolor':'w','pad':0}
            )

    #pl.xscale('log')
    #pl.yscale('log')
    ylim = pl.gca().get_ylim()
    min_ylim = 10**np.floor(np.log(ylim[0])/np.log(10))
    max_ylim = 10**np.ceil(np.log(ylim[1])/np.log(10))
    if min_ylim < 1:
        min_ylim = 1
    #pl.ylim([min_ylim, max_ylim])
    #pl.xlim([1,35])
    #if _r < n_row-1:
    #    [ x.set_visible(False) for x in ax[i].xaxis.get_major_ticks() ]
    pl.xlim([0,35])
    pl.xticks([0,10,20,30])
    bp.strip_axis(pl.gca())

pl.gcf().tight_layout()
pl.gcf().subplots_adjust(wspace=0.34,hspace=0.3)
pl.gcf().savefig("model_fit_figures/linlin_all_confirmed_fit_after_feb_12.png",dpi=300)
pl.gcf().savefig("model_fit_figures/SI_Fig_02.png",dpi=300)

if not loaded_fits:
    with open('fit_parameters/all_provinces_after_feb_12.p','wb') as f:
        pickle.dump(fit_parameters,f)

pl.show()
