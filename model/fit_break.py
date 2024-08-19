import numpy as np
from iminuit import Minuit
#from jacobi import propagate
import math
import random as rnd

def bpl(E, params):
    I0, E0, alpha, Eb, dalpha, s = params
    y = (I0 / 1e3) * np.power(E / E0, -alpha)
    y /= np.power(1. + np.power(E / Eb, s), dalpha / s)
    return y

def chi2_single(x, mu, sigma):
    return np.power((x - mu) / sigma, 2.)

def load_data(filename, minEnergy, maxEnergy = 1e20):
    E, y, err_stat_lo, err_stat_up, err_sys_lo, err_sys_up = np.loadtxt(filename,usecols=(0,1,2,3,4,5),unpack=True)
    items = [i for i in range(len(E)) if (E[i] > minEnergy and E[i] < maxEnergy)]
    #err_stat_lo = np.sqrt(err_stat_lo * err_stat_lo + err_sys_lo * err_sys_lo)
    #err_stat_up = np.sqrt(err_stat_up * err_stat_up + err_sys_up * err_sys_up)
    return E[items], y[items], err_stat_lo[items], err_stat_up[items]

def fit_bpl(params, filename):
    def chi2_function(I0, E0, alpha, logEb, dalpha, s):
        xd, yd, errd_lo, errd_up = load_data(filename, 80., 1e4)
        chi2 = 0.
        for x_i, y_i, err_lo_i, err_up_i in zip(xd, yd, errd_lo, errd_up):
            m = bpl(x_i, [I0, E0, alpha, np.power(10., logEb), dalpha, s])
            if m > y_i:
                chi2 += chi2_single(m, y_i, err_up_i)
            else:
                chi2 += chi2_single(m, y_i, err_lo_i)
        return chi2

    I0, E0, alpha, logEb, dalpha, s = params

    m = Minuit(chi2_function, I0=I0, E0=E0, alpha=alpha, logEb=logEb, dalpha=dalpha, s=s)
    m.errordef = Minuit.LEAST_SQUARES
    
    m.fixed['E0'] = True
#    m.fixed['s'] = True

    m.limits['alpha'] = (2.9, 3.3)
#    m.limits['I0'] = (10., 50.)
    m.limits['logEb'] = (1., 4.)
#    m.limits['dalpha'] = (0.01, 0.2)
#    m.limits['s'] = (1., 20.)

    m.simplex()
    m.migrad()
    m.hesse()

    #print(m)

 #   print(m.fval, 38. - m.nfit - 1)

    return m.values, m.errors, m.fval
    
def fit_CALET():
    filename = 'data/CALET_e-e+_energy.txt'
    f = open('CALET_fit.txt', 'w')
    I0, E0, alpha, logEb, dalpha, s = 0.291, 80.0, 3.113, 2.882, 0.827, 3.470
    values, errors, fval = fit_bpl([I0, E0, alpha, logEb, dalpha, s], filename)
    for x,e in zip(values, errors):
        print(f'{x:6.3f} {e:6.3f}')
        f.write(f'{x:6.3f} {e:6.3f}\n')
    f.close()
    print('running fit CALET done!')

def fit_DAMPE():
    filename = 'data/DAMPE_e-e+_energy.txt'
    f = open('DAMPE_fit.txt', 'w')
    I0, E0, alpha, logEb, dalpha, s = 0.323, 80.0, 3.088, 3.038, 1.152, 5.
    values, errors, fval = fit_bpl([I0, E0, alpha, logEb, dalpha, s], filename)
    for x,e in zip(values, errors):
        print(f'{x:6.3f} {e:6.3f}')
        f.write(f'{x:6.3f} {e:6.3f}\n')
    f.close()
    print('running fit DAMPE done!')
    
if __name__== "__main__":
    fit_CALET()
    fit_DAMPE()
