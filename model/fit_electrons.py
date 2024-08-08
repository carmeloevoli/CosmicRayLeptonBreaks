import numpy as np
from iminuit import Minuit
from jacobi import propagate
import math

def bpl(E, params):
    I0, E0, alpha, Eb, dalpha, s = params
    y = I0 * np.power(E / E0, -alpha)
    y /= np.power(1. + np.power(E / Eb, s), dalpha / s)
    return y
    
def spl(E, params):
    I0, E0, alpha = params
    y = I0 * np.power(E / E0, -alpha)
    return y
    
def gaussian1d(x, mu, sigma):
    sqrt2pi_sigma = np.sqrt(2. * math.pi) * sigma
    return 1. / sqrt2pi_sigma * np.exp(-0.5 * np.power((x - mu) / sigma, 2.))

def load_data(minEnergy, maxEnergy = 1e20):
    filename = 'data/AMS-02_e-_minus_e+_rigidity.txt'
    slope = 3.3
    E, y, err_stat_lo, err_stat_up, err_sys_lo, err_sys_up = np.loadtxt(filename,usecols=(0,1,2,3,4,5),unpack=True)
    y = np.power(E, slope) * y
    y_err_lo = np.power(E, slope) * err_stat_lo
    y_err_up = np.power(E, slope) * err_stat_up
    items = [i for i in range(len(E)) if (E[i] > minEnergy and E[i] < maxEnergy)]
    return E[items], y[items], y_err_lo[items], y_err_up[items]

def fit_bpl(params):
    def experiment_chi2(params):
        xd, yd, errd_lo, errd_up = load_data(20., 700.)
        chi2 = 0.
        for x_i, y_i, err_lo_i, err_up_i in zip(xd, yd, errd_lo, errd_up):
            m = bpl(x_i, params)
            if m > y_i:
                chi2 += -np.log(gaussian1d(m, y_i, err_up_i))
            else:
                chi2 += -np.log(gaussian1d(m, y_i, err_lo_i))
        return chi2
    
    def chi2_function(I0, E0, alpha, Eb, dalpha, s):
        chi2 = 0.
        chi2 += experiment_chi2([I0, E0, alpha, Eb, dalpha, s])
        return chi2

    I0, E0, alpha, Eb, dalpha, s = params

    m = Minuit(chi2_function, I0=I0, E0=E0, alpha=alpha, Eb=Eb, dalpha=dalpha, s=s)
    m.errordef = Minuit.LIKELIHOOD

    m.fixed['E0'] = True

    m.simplex()
    m.migrad()
    m.hesse()

    print(m)

    return m.values, m.errors, m.fval

def fit_spl(params):
    def experiment_chi2(params):
        xd, yd, errd_lo, errd_up = load_data(20., 700.)
        chi2 = 0.
        for x_i, y_i, err_lo_i, err_up_i in zip(xd, yd, errd_lo, errd_up):
            m = spl(x_i, params)
            if m > y_i:
                chi2 += -np.log(gaussian1d(m, y_i, err_up_i))
            else:
                chi2 += -np.log(gaussian1d(m, y_i, err_lo_i))
        return chi2
    
    def chi2_function(I0, E0, alpha):
        chi2 = 0.
        chi2 += experiment_chi2([I0, E0, alpha])
        return chi2

    I0, E0, alpha = params

    m = Minuit(chi2_function, I0=I0, E0=E0, alpha=alpha)
    m.errordef = Minuit.LIKELIHOOD

    m.fixed['E0'] = True

    m.simplex()
    m.migrad()
    m.hesse()

    print(m)

    return m.values, m.errors, m.fval
    
if __name__== "__main__":
    I0, E0, alpha, Eb, dalpha, s = 420.8, 20., 0.022, 40., 0., 10.0

    values, errors, fval = fit_spl([I0, E0, alpha])

    values, errors, fval = fit_bpl([I0, E0, alpha, Eb, dalpha, s])

#    print(fval)
