import numpy as np
from iminuit import Minuit
#from jacobi import propagate
import math

def bpl(E, params):
    I0, E0, alpha, Eb, dalpha, s = params
    y = (I0 / 1e3) * np.power(E / E0, -alpha)
    y /= np.power(1. + np.power(E / Eb, s), dalpha / s)
    return y
    
def spl(E, params):
    I0, E0, alpha = params
    y = (I0 / 1e3) * np.power(E / E0, -alpha)
    return y
    
def chi2_single(x, mu, sigma):
    return np.power((x - mu) / sigma, 2.)

def load_data(minEnergy, maxEnergy = 1e20):
    filename = 'data/AMS-02_e-_minus_e+_rigidity.txt'
    E, y, err_stat_lo, err_stat_up, err_sys_lo, err_sys_up = np.loadtxt(filename,usecols=(0,1,2,3,4,5),unpack=True)
    items = [i for i in range(len(E)) if (E[i] > minEnergy and E[i] < maxEnergy)]
    err_stat_lo = np.sqrt(err_stat_lo * err_stat_lo + err_sys_lo * err_sys_lo)
    err_stat_up = np.sqrt(err_stat_up * err_stat_up + err_sys_up * err_sys_up)
    return E[items], y[items], err_stat_lo[items], err_stat_up[items]

def fit_bpl(params):
    def chi2_function(I0, E0, alpha, Eb, dalpha, s):
        xd, yd, errd_lo, errd_up = load_data(20., 1000.)
        chi2 = 0.
        for x_i, y_i, err_lo_i, err_up_i in zip(xd, yd, errd_lo, errd_up):
            m = bpl(x_i, [I0, E0, alpha, Eb, dalpha, s])
            if m > y_i:
                chi2 += chi2_single(m, y_i, err_up_i)
            else:
                chi2 += chi2_single(m, y_i, err_lo_i)
        return chi2

    I0, E0, alpha, Eb, dalpha, s = params

    m = Minuit(chi2_function, I0=I0, E0=E0, alpha=alpha, Eb=Eb, dalpha=dalpha, s=s)
    m.errordef = Minuit.LEAST_SQUARES
    
    m.fixed['E0'] = True
    m.limits['alpha'] = (3.1, 3.4)
    m.limits['I0'] = (10., 50.)
    m.limits['Eb'] = (20., 100.)
    m.limits['dalpha'] = (-0.1, 0.)
    m.limits['s'] = (5., 15.)

    m.simplex()
    m.migrad()
    m.hesse()

    #print(m)

    print(m.fval, 38. - m.nfit - 1)

    return m.values, m.errors, m.fval

def fit_spl(params):
    def chi2_function(I0, E0, alpha):
        xd, yd, errd_lo, errd_up = load_data(20., 1000.)
        chi2 = 0.
        for x_i, y_i, err_lo_i, err_up_i in zip(xd, yd, errd_lo, errd_up):
            m = spl(x_i, [I0, E0, alpha])
            if m > y_i:
                chi2 += chi2_single(m, y_i, err_up_i)
            else:
                chi2 += chi2_single(m, y_i, err_lo_i)
        return chi2
    
    I0, E0, alpha = params

    m = Minuit(chi2_function, I0=I0, E0=E0, alpha=alpha)
    m.errordef = Minuit.LEAST_SQUARES

    m.fixed['E0'] = True
    m.limits['alpha'] = (3.1, 3.4)
    m.limits['I0'] = (10., 50.)

    m.simplex()
    m.migrad()
    m.hesse()

    #print(m)

    print(m.fval, 38. - m.nfit - 1)

    return m.values, m.errors, m.fval
    
if __name__== "__main__":
    I0, E0, alpha = 21.790, 20., 3.2788
    
    values, errors, fval = fit_spl([I0, E0, alpha])

    I0, E0, alpha, Eb, dalpha, s = 22.02, 20.0, 3.322, 37.2, -0.076, 11

    values, errors, fval = fit_bpl([I0, E0, alpha, Eb, dalpha, s])
