import numpy as np
from scipy.integrate import ode
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters
import json
from tqdm import tqdm
from bfmplot import pl
from bfmplot import brewer_qualitative, simple_cycler, markers

import bfmplot as bp
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)

#brewer_qualitative[6] = brewer_qualitative[1]

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

tuplelist = [ (p, d)  for p, d in data.items()\
                               if max(d['cases']) >= 20\
                               and len(np.where(np.logical_and(np.array(d['times'])<=12,
                                                               np.array(d['cases'])>0))[0])\
                                   >= 8
                             ]

colors

tuplelist = sorted([ t for t in tuplelist ],key=lambda x: -max(x[1]['cases']))

n_fits = len(tuplelist)
n_col = int(np.ceil(np.sqrt(n_fits)))
n_row = 1
n_col = 2
fig, ax = pl.subplots(n_row,n_col,figsize=(8,3),sharex=True)
ax = ax.flatten()

titlemap = REPL({'mainland_china':'all except Hubei'})

letter = "abcdefg"
roman = [ "i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix"]

i = -1
for province, pdata in tqdm(tuplelist[:2]):
    i += 1
    print(province)

    t = np.array(pdata['times'])
    cases = np.array(pdata['cases'])
    dates = np.array(pdata['dates'],dtype=np.datetime64)
    print(dates)

    if max(cases) <= 20:
        continue


    i0 = np.where(cases>0)[0][0]
    t = t[i0:]
    cases = cases[i0:]
    dates = dates[i0:]

    if province != 'Hubei':
        maxt = 12
    else:
        maxt = 19
    i1 = np.where(t<=maxt)[0][-1]
    _t = t[:i1+1]
    _cases = cases[:i1+1]
    print(t, cases)

    if len(t) < 8:
        continue

    f = lambda x, mu, B: B*x**mu
    p,_ = curve_fit(f, _t, _cases, [1.9,4.5])

    print(p)

    tt = np.logspace(np.log(t[0]), np.log(t[-1]), base=np.exp(1))
    tt_dates = np.array( (tt-1) *24*3600 ,np.timedelta64) + dates[0]
    print(dates, tt_dates)

    pl.sca(ax[i])

    pl.plot_date(dates, cases,marker=markers[i],c=colors[i],label=titlemap[province],mfc='None')
    pl.plot_date(tt_dates, f(tt,*p),'-',c='k',lw=1,label='$t^\mu$, $\mu=%4.2f$' % p[0])

    _c = i % (n_col)
    _r = i // (n_col)
    #if _r == n_row-1:
        #pl.xlabel('days since Jan. 20')
    pl.ylabel('confirmed cases',)
        #pl.gca().yaxis.set_label_coords(-0.3,-0.2)
    #pl.xlim([1,30])
    #pl.xscale('log')
    #pl.yscale('log')

    ylim = ax[i].set_ylim([10,4e4])
    ylim = ax[i].get_ylim()
    if i == 0:
        ax[i].plot([np.datetime64('2020-02-09 16:00')]*2,[10,40000],':')
    else:
        ax[i].plot([np.datetime64('2020-02-02 16:00')]*2,[10,25000],':')



    ax[i].text(-0.18,1.03,
            "{}".format(letter[i].upper()),
            transform=ax[i].transAxes,
            ha='left',
            va='top',
            fontweight='bold',
            fontsize=14,
            bbox={'facecolor':'w','edgecolor':'w','pad':0}
            )
    if i == 1:
        ax[i].text(0.55,0.6,
                "Feb. 2nd".format(p[0]),
                transform=ax[i].transAxes,
                ha='left',
                va='bottom',
                fontsize=10,
                bbox={'facecolor':'w','edgecolor':'w','pad':0}
                )
    else:
        ax[i].text(0.78,0.2,
                "Feb. 9th".format(p[0]),
                transform=ax[i].transAxes,
                ha='center',
                va='bottom',
                fontsize=10,
                bbox={'facecolor':'w','edgecolor':'w','pad':0}
                )
    #ax[i].text([0.2,0.3][i],0.2,
    #       "$\mu={0:4.2f}$".format(p[0]),
    #       transform=ax[i].transAxes,
    #       ha='left',
    #       va='bottom',
    #       bbox={'facecolor':'w','edgecolor':'w','pad':0}
    #       )
           
    #ax[i].text(0.97,0.03,
    #        titlemap[province],
    #        transform=ax[i].transAxes,
    #        ha='right',
    #        va='bottom',
    #        bbox={'facecolor':'w','edgecolor':'w','pad':0}
    #        )
    #if _r < n_row-1:
    #    [ x.set_visible(False) for x in ax[i].xaxis.get_major_ticks() ]
    #ax[i].set_yticks([10,100,1000,10000])

    rule = rrulewrapper(DAILY, interval=4)
    loc = RRuleLocator(rule)
    #ax.set_yscale('log')
    formatter = DateFormatter('%d. %m.')
    ax[i].xaxis.set_major_locator(loc)    
    ax[i].xaxis.set_major_formatter(formatter)
    ax[i].xaxis.set_tick_params(rotation=30, labelsize=10)
    ax[i].set_yticks([0,10000,20000,30000,40000])
    bp.humanify_yticks(ax[i],precision=0)
    bp.humanify_yticks(ax[i],precision=0)
    bp.strip_axis(pl.gca())

ax[0].legend(loc=(0.495,0.755))
leg = ax[1].legend(loc='upper right')
bp.align_legend_right(leg)

#ax[0].text(-0.4,1.1,
#           'C',
#            transform=ax[0].transAxes,
#            ha='left',
#            va='top',
#            fontweight='bold',
#            fontsize=14,
#            bbox={'facecolor':'w','edgecolor':'w','pad':0}
#        )


    
pl.gcf().tight_layout()
pl.gcf().subplots_adjust(wspace=0.3,hspace=0.3)
pl.gcf().savefig("powerlaw_fit_figures/fit_powerlaw_large_hubei_china.png",dpi=300)



pl.show()
