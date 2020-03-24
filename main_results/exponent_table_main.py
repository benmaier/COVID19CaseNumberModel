import pickle

import sys
sys.path.insert(0,'..')
from SIRX import SIRXConfirmedModel
model = SIRXConfirmedModel()

import numpy as np
import simplejson as json
from tabulate import tabulate
import numpy
from bfmplot import pl
from scipy.special import erf

def mu(alpha, beta, kappa, kappa0):
    R0 = alpha/(beta+kappa+kappa0)
#    print(np.log(R0))
#    return ( numpy.e )**( -1/2 * ( R0 )**( -1/2 ) * ( alpha )**( -1 ) * ( \
#kappa0 )**( -1 ) * ( ( alpha + -1 * ( R0 )**( 1/2 ) * ( beta + ( \
#kappa + kappa0 ) ) ) )**( 2 ) ) * ( 2 * numpy.pi )**( -1/2 ) * ( R0 \
#)**( -1/4 ) * ( alpha )**( 1/2 ) * ( kappa0 )**( -1/2 ) * ( ( \
#erf( ( 2 )**( -1/2 ) * ( R0 )**( -1/4 ) * ( alpha )**( -1/2 ) * \
#( kappa0 )**( -1/2 ) * ( -1 * alpha + ( R0 )**( 1/2 ) * ( beta + ( \
#kappa + kappa0 ) ) ) ) + -1 * erf( 1/2 * ( 2 )**( -1/2 ) * ( R0 \
#)**( -1/4 ) * ( alpha )**( -1/2 ) * ( kappa0 )**( -1/2 ) * ( -2 * \
#alpha + ( 2 * ( R0 )**( 1/2 ) * ( beta + ( kappa + kappa0 ) ) + -1 * \
#alpha * numpy.log( R0 ) ) ) ) ) )**( -1 ) * numpy.log( R0 )

    return ( numpy.e )**( -1/2 * ( R0 )**( -1/2 ) * ( alpha )**( -1 ) * ( kappa0 \
)**( -1 ) * ( ( alpha + -1 * ( R0 )**( 1/2 ) * ( beta + ( kappa + \
kappa0 ) ) ) )**( 2 ) ) * ( 2 * numpy.pi )**( -1/2 ) * ( R0 )**( -1/4 \
) * ( alpha )**( 1/2 ) * ( kappa0 )**( -1/2 ) * ( ( erf( ( 2 \
)**( -1/2 ) * ( R0 )**( -1/4 ) * ( alpha )**( -1/2 ) * ( kappa0 )**( \
-1/2 ) * ( -1 * alpha + ( R0 )**( 1/2 ) * ( beta + ( kappa + kappa0 ) \
) ) ) + -1 * erf( 1/2 * ( 2 )**( -1/2 ) * ( R0 )**( -1/4 ) * ( \
alpha )**( -1/2 ) * ( kappa0 )**( -1/2 ) * ( -2 * alpha + ( 2 * ( R0 \
)**( 1/2 ) * ( beta + ( kappa + kappa0 ) ) + -1 * alpha * numpy.log( \
R0 ) ) ) ) ) )**( -1 ) * numpy.log( R0 )

def mu2(alpha, beta, kappa, kappa0):
    R0 = alpha/(beta+kappa+kappa0)
    return ( numpy.e )**( -1/2 * ( alpha )**( -1 ) * ( kappa0 )**( -1 ) * ( ( -1 \
* alpha + ( beta + ( kappa + ( kappa0 + 1/2 * alpha * numpy.log( R0 ) \
) ) ) ) )**( 2 ) ) * ( 2 * numpy.pi )**( -1/2 ) * ( alpha )**( 1/2 ) \
* ( kappa0 )**( -1/2 ) * ( ( -1 * erf( ( 2 )**( -1/2 ) * ( \
alpha )**( -1/2 ) * ( kappa0 )**( -1/2 ) * ( -1 * alpha + ( beta + ( \
kappa + kappa0 ) ) ) ) + erf( ( 2 )**( -1/2 ) * ( alpha )**( \
-1/2 ) * ( kappa0 )**( -1/2 ) * ( -1 * alpha + ( beta + ( kappa + ( \
kappa0 + 1/2 * alpha * numpy.log( R0 ) ) ) ) ) ) ) )**( -1 ) * \
numpy.log( R0 )


class REPL(dict):

    def __init__(self, items):
        self.items = items

    def __getitem__(self,i):
        try:
            return self.items[i]
        except KeyError as e:
            return i

with open('fit_parameters/all_provinces_before_feb_12.p','rb') as f:
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
fit_values = [2.32,1.92,2.04,2.12,2.16,2.09,2.05,2.75,2.20,1.45]
for i, (province, pdata) in enumerate(tuplelist[:10]):
    p = fit_parameters[province]
    k = p['kappa'].value
    k0 = p['kappa0'].value
    b = p['rho'].value
    a = p['eta'].value
    I0f = p['I0_factor'].value

    N = pdata['population'] / 1e6

    Q = (k+k0)/(b+k+k0)
    P = k0/(k+k0)
    R0eff = a / (b+k+k0)
    tabledata.append([titlemap[province], fit_values[i], mu(a,b,k,k0), np.abs( 1-fit_values[i]/mu(a,b,k,k0))])

    pl.figure()
print(tabulate(tabledata,
               headers=['Province',r'$\mu_{\mathrm{\fit}}$', r'$\mu_{\mathrm{approx}}$', 'rel. err.'],
               floatfmt=("",".2f",".2f",".2f") 
               ))
latex = str(tabulate(tabledata,
            tablefmt="latex_raw",
           headers=['Province',r'$\mu_{\mathrm{\fit}}$', r'$\mu_{\mathrm{approx}}$', 'rel. err.'],
           floatfmt=("",".2f",".2f",".2f") 
            ))

print("affected provinces but hubei and aggregated mainland, mu = ", np.array(fit_values[2:]).mean(),"+/-", np.array(fit_values[2:]).std())


with open('fit_parameters/exponents.tex','w') as f:
    f.write(latex)

