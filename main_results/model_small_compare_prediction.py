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

from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)

import bfmplot as bp

import matplotlib.pyplot as mpl
mpl.rcParams['font.size'] = 11

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
n_row = n_col-1
#ax = ax.flatten()

titlemap = REPL({'mainland_china':'all except Hubei'})

if loaded_fits:
    with open(pickle_filename,'rb') as f:
        fit_parameters = pickle.load(f)
else:
    fit_parameters = {}

letter = "abcdefg"
roman = [ "i", "ii", "iii", "iv", "v", "vi"]

ylims = [[100,100000],[10,20000]]
max_dates = ['Feb. 7th', 'Feb. 1st.']
change_date = np.datetime64('2020-02-12 18:00')
max_dates_pos = [60000,15000]
max_dates_va = ['bottom']*2


i = -1
fig, ax = pl.subplots(1,2,figsize=(3.5,2))
curves = {}
for province, pdata in tqdm(tuplelist[:2]):
    i += 1

    t = np.array(pdata['times'])
    cases = np.array(pdata['cases'])
    dates = np.array(pdata['dates'],np.datetime64)

    if max(cases) <= 20:
        continue

    i0 = np.where(cases>0)[0][0]
    tswitch = t[-1]
    t = t[i0:]
    cases = cases[i0:]

    t2 = np.array(datafeb[province]['times']) + 1 - 5/24
    cases2 = np.array(datafeb[province]['cases'])
    dates2 = np.array(datafeb[province]['dates'],np.datetime64)
    print(dates2, cases2)
    i0 = np.where(dates2>=np.datetime64("2020-02-13"))[0][0]
    t2 = t2[i0:]
    cases2 = cases2[i0:]
    dates2 = dates2[i0:]



    if loaded_fits: 
        params = fit_parameters[province]
    else:
        out = model.fit(t,cases,maxfev=1000,N=pdata['population']
                )
        params = out.params
        fit_parameters[province] = params
    print(params)
    N = params['N']

    #pl.sca(ax[i])

    tt = np.logspace(np.log(t[0]), np.log(50), 1000,base=np.exp(1))
    tt1 = tt[tt<=tswitch] 
    tt2 = tt[tt>tswitch] 
    tt_dates = np.array( (tt-t[0]) *24*3600 ,np.timedelta64) + dates[0]
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
    S = result[0,:]*N
    Z = result[3,:]*N
    imax = np.argmax(I)
    print(imax)
    max_date = tt_dates[imax]
    max_tt = tt[imax]
    print("=======", province, "max", max_tt, max_date)

    curves[i] = {}
    curves[i]['I'] = {'x': tt_dates, 'y': I}
    curves[i]['X'] = {'x': tt_dates, 'y': X}
    curves[i]['S'] = {'x': tt, 'y': S}
    curves[i]['S+Z'] ={ 'x': tt, 'y':S+Z}


    Xlabel = None
    Ilabel = None
    ax[i].plot(t, cases,marker=markers[i],c=colors[i],ms=3,label=None,mfc='None',)
    ax[i].plot(t2, cases2,marker='o',ms=5,c='grey',label=None,mfc='None',)
    ax[i].plot(tt1, X[tt<=tswitch],'-',c='k',label=Xlabel)
    ax[i].plot(tt2, X[tt>tswitch],'-',c='k',label=None)
    ax[i].plot(tt, I,'--',c=colors[2],lw=1.5,label=Ilabel)
    ax[i].plot([max_tt]*2, [0,max_dates_pos[i]],':',c=colors[0],lw=1.5)
    ax[i].plot([23,23], [0,[60000,15000][i]],':',c=colors[0],lw=1.5)


    #ax[0].plot(tt, Q,c='k',label='$Q_I$ (detected and quarantined)')
    #ax[i].text(max_tt, max_dates_pos[i], max_dates[i],
    #        transform=ax[i].transData,
    #        ha='right',
    #        va=max_dates_va[i],
    #        color=colors[0],
    #        fontsize=9,
    #        bbox={'facecolor':'w','edgecolor':'w','pad':0}
    #        )


    #ax[0].plot(tt, Q,c='k',label='$Q_I$ (detected and quarantined)')
    #ax[0].plot(tt, I,'--',c=colors[1],lw=1.5,label='$I$ (undected infected)')
    #ax[i].plot(tt, S,'-.',c=colors[5],lw=1.5,label=None)
    #ax[i].plot(tt, S+Z,':',c=colors[3],lw=1.5,label=None)

#pl.plot(tt, S,label='model')
    
    _c = i % n_col
    _r = i // n_col
    #pl.title(titlemap[province])
    #ax[i].text(0.03,0.97,
    #        "{}.{}".format(letter[_r], roman[_c]),
    #        transform=ax[i].transAxes,
    #        ha='left',
    #        va='top',
    #        fontweight='bold',
    #        fontsize=12,
    #        bbox={'facecolor':'w','edgecolor':'w','pad':0}
    #        )
    #ax[i].text(0.03,0.8,
    #        titlemap[province],
    #        transform=ax[i].transAxes,
    #        ha='left',
    #        va='top',
    #        bbox={'facecolor':'w','edgecolor':'w','pad':0}
    #        )
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

    ax[i].set_xscale('log')
    ax[i].set_yscale('log')
    #ylim = pl.gca().get_ylim()
    #min_ylim = 10**np.floor(np.log(ylim[0])/np.log(10))
    #max_ylim = 10**np.ceil(np.log(ylim[1])/np.log(10))
    #if min_ylim < 1:
    #    min_ylim = 1
    ax[i].set_ylim(ylims[i])
    #ax[i].set_xlim([1,40])
    #if i == 1:
    #    ax[i].set_ylim([10,20000])

    bp.strip_axis(ax[0])
    bp.strip_axis(ax[1])
    #bp.align_legend_right(leg)
    #ax[i].set_yticks([0,10000,20000,30000,40000])
    if i == 0:
        ax[i].set_yticks([1e2, 1e3, 1e4, 1e5])
    else:
        ax[i].set_yticks([1e1, 1e2, 1e3, 1e4])
    #bp.humanify_yticks(ax[i],precision=0 if i == 0 else 1)
    ax[i].set_xlabel('days since Jan. 20th')

#ax[0].text(-0.18,1.03,'A',fontsize=14,fontweight='bold',transform=ax[0].transAxes,va='top')
#ax[1].text(-0.21,1.03,'B',fontsize=14,fontweight='bold',transform=ax[1].transAxes,va='top')

pl.gcf().tight_layout()

#bp.add_curve_label(ax[0],
#                   curve_x=curves[0]['I']['x'],
#                   curve_y=curves[0]['I']['y'],
#                   label='$I$',
#                   label_pos_rel=0.2,
#                   angle=0,
#                   #color=colors[2],
#                   fontsize=12
#                   )
#
#bp.add_curve_label(ax[0],
#                   curve_x=curves[0]['S']['x'],
#                   curve_y=curves[0]['S']['y'],
#                   label='$S$',
#                   label_pos_rel=0.6,
#                   angle=0,
#                   #color=colors[2],
#                   fontsize=12
#                   )
#
#bp.add_curve_label(ax[0],
#                   curve_x=curves[0]['X']['x'],
#                   curve_y=curves[0]['X']['y'],
#                   label='$X$',
#                   label_pos_rel=0.9,
#                   angle=0,
#                   #color=colors[2],
#                   fontsize=12
#                   )
#
#bp.add_curve_label(ax[1],
#                   curve_x=curves[1]['I']['x'],
#                   curve_y=curves[1]['I']['y'],
#                   label='$I$',
#                   label_pos_rel=0.2,
#                   angle=0,
#                   #color=colors[2],
#                   fontsize=12
#                   )
#
#bp.add_curve_label(ax[1],
#                   curve_x=curves[1]['S']['x'],
#                   curve_y=curves[1]['S']['y'],
#                   label='$S$',
#                   label_pos_rel=0.6,
#                   angle=0,
#                   #color=colors[2],
#                   fontsize=12
#                   )
#
#bp.add_curve_label(ax[1],
#                   curve_x=curves[1]['X']['x'],
#                   curve_y=curves[1]['X']['y'],
#                   label='$X$',
#                   label_pos_rel=0.9,
#                   angle=0,
#                   #color=colors[2],
#                   fontsize=12
#                   )


#pl.gcf().subplots_adjust(wspace=0.3,hspace=0.3)
pl.gcf().savefig("model_fit_figures/small_hubei_and_mainland_china.png",dpi=300,transparent=True)

if not loaded_fits:
    with open('fit_parameters/hubei_china.p','wb') as f:
        pickle.dump(fit_parameters,f)

pl.show()
