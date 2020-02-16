import numpy as np
from scipy.integrate import ode
from lmfit import minimize, Parameters

class SIRXConfirmedModel:

    def __init__(self):
        pass


    # set equation of motion for SIRX dynaics
    def dxdt(self,t,y,eta,rho,kappa,kappa0):

        S = y[0]
        I = y[1]
        X = y[2]
        H = y[3]
     
        dy = np.zeros(4)
        dy[0] = -eta*S*I - kappa0*S
        dy[1] = +eta*S*I - rho*I - kappa*I - kappa0*I
        dy[2] = +kappa*I + kappa0*I
        dy[3] = +kappa0*S


        return dy

    def SIRX(self,t, y0, eta, rho, kappa,kappa0, N, I0_factor):

        X0 = y0 / N
        I0 = X0 * I0_factor
        S0 = 1-X0-I0
        y0 = np.array([S0, I0, X0, 0.0])
        t0 = t[0]

        t = t[1:]

        r = ode(self.dxdt)

        # Runge-Kutta with step size control
        r.set_integrator('dopri5')

        # set initial values
        r.set_initial_value(y0,t0)

        # set transmission rate and recovery rate
        r.set_f_params(eta,rho,kappa,kappa0)

        result = np.zeros((4,len(t)+1))
        result[:,0] = y0

        # loop through all demanded time points
        for it, t_ in enumerate(t):

            # get result of ODE integration
            y = r.integrate(t_)

            # write result to result vector
            result[:,it+1] = y

        return result

    def residual(self,params, x, data):

        eta = params['eta']
        rho = params['rho']
        kappa = params['kappa']
        kappa0 = params['kappa0']
        I0_factor = params['I0_factor']
        #N = 10**params['log10N']
        N = params['N']

        result = self.SIRX(x, data[0], eta, rho, kappa, kappa0, N, I0_factor)
        X = result[2,:]

        residual = X*N - data

        return residual

    def fit(self,t, data,maxfev=1000,params=None,N=None,Nmax=None):

        if params is None:
            params = Parameters()
            R0 = 6.2
            rho = 1/8
            eta = R0*rho
            params.add('eta',value=eta,vary=False)
            params.add('rho',value=rho, vary=False)
            params.add('kappa',value=rho,min=0)        
            params.add('kappa0',value=rho/2,min=0)
            #params.add('xi',value=1/2,min=0.2,max=0.8,vary=True)
            #params.add('xi',value=1/2,min=0,max=1,vary=True)
            params.add('I0_factor', value=10,min=1)
            varyN = N is None
            if varyN:
                N = 1e7
            if Nmax is None:
                Nmax=115000000
            params.add('N',value=N,min=100000,max=Nmax,vary=varyN)

        out = minimize(self.residual, params, args=(t, data, ),maxfev=maxfev)
        return out

class SIRXShutdownModel:

    def __init__(self):
        pass


    # set equation of motion for SIRX dynaics
    def dxdt(self,t,y,eta,rho,kappa,xi):

        S = y[0]
        I = y[1]
        Q = y[2]
        H = y[3]
     
        dy = np.zeros(4)
        dy[0] = -eta*S*I - xi*S
        dy[1] = +eta*S*I - rho*I - xi*I
        dy[2] = +kappa*I
        dy[3] = +kappa*xi*S


        return dy

    def SIRX(self,t, y0, eta, rho, kappa, xi, N, I0_factor):

        Q0 = y0 / N
        I0 = Q0 * I0_factor
        S0 = 1-Q0-I0
        y0 = np.array([S0, I0, Q0, 0.0])
        t0 = t[0]

        t = t[1:]

        r = ode(self.dxdt)

        # Runge-Kutta with step size control
        r.set_integrator('dopri5')

        # set initial values
        r.set_initial_value(y0,t0)

        # set transmission rate and recovery rate
        r.set_f_params(eta,rho,kappa,xi)

        result = np.zeros((4,len(t)+1))
        result[:,0] = y0

        # loop through all demanded time points
        for it, t_ in enumerate(t):

            # get result of ODE integration
            y = r.integrate(t_)

            # write result to result vector
            result[:,it+1] = y

        return result

    def residual(self,params, x, data):

        eta = params['eta']
        rho = params['rho']
        kappa = params['kappa']
        xi = params['xi']
        I0_factor = params['I0_factor']
        #N = 10**params['log10N']
        N = params['N']

        result = self.SIRX(x, data[0], eta, rho, kappa, xi, N, I0_factor)
        Q = result[2,:]

        residual = Q*N - data

        return residual

    def fit(self,t, data,maxfev=1000,params=None,N=None,Nmax=None):

        if params is None:
            params = Parameters()
            R0 = 6.2
            rho = 1/8
            eta = R0*rho
            params.add('eta',value=eta,vary=False)
            params.add('rho',value=rho, vary=False)
            params.add('kappa',value=rho)        
            #params.add('xi',value=1/2,min=0.2,max=0.8,vary=True)
            params.add('xi',value=rho)
            params.add('I0_factor', value=10,min=1)
            varyN = N is None
            if varyN:
                N = 1e7
            if Nmax is None:
                Nmax=115000000
            params.add('N',value=N,min=10000,max=Nmax,vary=varyN)

        out = minimize(self.residual, params, args=(t, data, ),maxfev=maxfev)
        return out
