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

with open('../data/all_confirmed_cases_with_population.json','r') as f:
    data = json.load(f)
with open('../data/all_confirmed_csse_cases_with_population.json','r') as f:
    datafeb = json.load(f)

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
n_row = 2
n_col = 4
fig, ax = pl.subplots(n_row,n_col,figsize=(8,3))
ax = ax.flatten()

titlemap = REPL({'mainland_china':'All w/o Hubei'})

if loaded_fits:
    with open(pickle_filename,'rb') as f:
        fit_parameters = pickle.load(f)
else:
    fit_parameters = {}

letter = "abcdefg"
roman = [ "i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix"]

max_dates = ['Feb. 2nd.', 'Feb. 2nd.', 'Jan. 31st', 'Feb. 1st', 'Feb. 4th', 'Feb. 3rd', 'Feb. 3rd', 'Feb. 1st']
max_dates_pos = [1500, 1500, 1500, 1500, 1500, 1500, 1500, 900]
max_dates_va = 4*['bottom'] + 4*['top']


i = -1
for province, pdata in tqdm(tuplelist[2:10]):
    i += 1

    t = np.array(pdata['times'])
    cases = np.array(pdata['cases'])
    dates = np.array(pdata['dates'],dtype=np.datetime64)
    if pdata['dates'][0] == "2020-01-22 12:00:00":
        t += 14/24
        print("===================Jiangsu=============")

    if max(cases) <= 20:
        continue

    i0 = np.where(cases>0)[0][0]
    tswitch = t[-1]
    t = t[i0:]
    cases = cases[i0:]
    dates = dates[i0:]

    t2 = np.array(datafeb[province]['times']) + (1-5/24) # adjust for shift of t0 in data
    cases2 = np.array(datafeb[province]['cases'])
    dates2 = np.array(datafeb[province]['dates'],np.datetime64)
    print(province)
    print(t[0], dates[0])
    print(t2[0], dates2[0])
    #print(dates2, cases2)
    i0 = np.where(dates2>=np.datetime64("2020-02-13"))[0][0]
    t2 = t2[i0:]
    cases2 = cases2[i0:]
    dates2 = dates2[i0:]


    print(pdata['population'])

    if loaded_fits: 
        params = fit_parameters[province]
    else:
        out = model.fit(t,cases,maxfev=1000,N=pdata['population']
                )
        params = out.params
        fit_parameters[province] = params
    print(params)
    N = params['N']

    pl.sca(ax[i])

    tt = np.logspace(np.log(t[0]), np.log(30), 1000,base=np.exp(1))
    tt1 = tt[tt<=tswitch] 
    tt2 = tt[tt>tswitch] 
    tt_dates = np.array( (tt-t[0]) *24*3600 ,np.timedelta64) + dates[0]
    print(tt[0], tt_dates[0])
    tt1_dates = tt_dates[tt<=tswitch] 
    tt2_dates = tt_dates[tt>tswitch] 
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
    imax = np.argmax(I)
    print(imax)
    max_date = tt_dates[imax]
    max_tt = tt[imax]
    print("=======", province, "max", max_tt, max_date)
    print("=======", province, "max", tt[100], tt_dates[100])


#S = result[0,:]*N

    pl.plot(t, cases,marker=markers[i+2],c=colors[i+2],label='data',mfc='None')
    pl.plot(t2, cases2,marker='o',ms=5,c='grey',label='data',mfc='None')
    pl.plot(tt, X,'-',c='k')
    #pl.plot(tt1, X[tt<=tswitch],'-',c='k')
    #pl.plot(tt2, X[tt>tswitch],'-',c='k')
    pl.plot(tt, I,'--',c=colors[2],lw=1.5)
    pl.plot([max_tt]*2, [0,max_dates_pos[i]],':',c=colors[0],lw=1.5)
    pl.text(max_tt-1, max_dates_pos[i], max_dates[i],
            transform=ax[i].transData,
            ha='left',
            va=max_dates_va[i],
            color=colors[0],
            fontsize=9,
            bbox={'facecolor':'w','edgecolor':'w','pad':0}
            )
    ax[i].plot([22.5,22.5], [0,cases[-1]],':',c=colors[0],lw=1.5)

#pl.plot(tt, S,label='model')
    
    _c = i % n_col
    _r = i // n_col
    if _r == n_row-1:
        pl.xlabel('days since Jan. 20th')        
    if _c == 0 and _r == 0:
        pl.ylabel('confirmed cases')
        pl.gca().yaxis.set_label_coords(-0.3,-0.2)
    #pl.title(titlemap[province])
    ax[i].text(0.03,0.97,
            "{}".format(roman[i]),
            transform=ax[i].transAxes,
            ha='left',
            va='top',
            fontweight='bold',
            fontsize=10,
            bbox={'facecolor':'w','edgecolor':'w','pad':0}
            )
    ax[i].text(0.03,0.8,
            titlemap[province],
            transform=ax[i].transAxes,
            ha='left',
            va='top',
            bbox={'facecolor':'w','edgecolor':'w','pad':0}
            )
    #ax[i].text(0.97,0.15,
    #        r"$P=%4.2f$" %(params['kappa'].value/(params['rho'].value+params['kappa'].value)),
    #        transform=ax[i].transAxes,
    #        ha='right',
    #        va='bottom',
    #        bbox={'facecolor':'w','edgecolor':'w','pad':0}
    #        )
    #ax[i].text(0.97,0.03,
    #        r"$\xi=%4.2f$" %(params['xi'].value),
    #        transform=ax[i].transAxes,
    #        ha='right',
    #        va='bottom',
    #        bbox={'facecolor':'w','edgecolor':'w','pad':0}
    #        )

    #pl.xscale('log')
    #pl.yscale('log')
    ylim = pl.gca().set_ylim([1,1.5e3])
    ylim = pl.gca().get_ylim()
    min_ylim = 10**np.floor(np.log(ylim[0])/np.log(10))
    max_ylim = 10**np.ceil(np.log(ylim[1])/np.log(10))
    if min_ylim < 1:
        min_ylim = 1
    #pl.ylim([min_ylim, max_ylim])
    pl.xlim([0,30])
    if _r < n_row-1:
        ax[i].set_xticklabels('')
    #    [ x.set_visible(False) for x in ax[i].xaxis.get_major_ticks() ]
    ax[i].set_yticks([0,500,1000,1500])
    ax[i].set_yticklabels(['0','0.5k','1.0k','1.5k'])
    bp.strip_axis(pl.gca())
    #bp.humanify_yticks(ax[i],precision=1)

ax[0].text(-0.4,1.1,
           'C',
            transform=ax[0].transAxes,
            ha='left',
            va='top',
            fontweight='bold',
            fontsize=14,
            bbox={'facecolor':'w','edgecolor':'w','pad':0}
          )

pl.gcf().tight_layout()
pl.gcf().subplots_adjust(wspace=0.3,hspace=0.3)
pl.gcf().savefig("model_fit_figures/model_fit_confirmed_500.png",dpi=300)

if not loaded_fits:
    with open('fit_parameters/confirmed_cases_500.p','wb') as f:
        pickle.dump(fit_parameters,f)

pl.show()
