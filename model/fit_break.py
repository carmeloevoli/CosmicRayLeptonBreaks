import numpy as np
from iminuit import Minuit
#from jacobi import propagate
import math
import random as rnd

def bpl(E, params):
    I0, E0, alpha, logEb, dalpha, s = params
    Eb = np.power(10., logEb)
    y = (I0 / 1e3) * np.power(E / E0, -alpha)
    y /= np.power(1. + np.power(E / Eb, s), dalpha / s)
    return y

def cutoff(E, params):
    I0, E0, alpha, logEb = params
    Eb = np.power(10., logEb)
    y = (I0 / 1e3) * np.power(E / E0, -alpha)
    y *= np.exp(- E / Eb)
    return y

def chi2_single(x, mu, sigma):
    return np.power((x - mu) / sigma, 2.)

def load_data(filename, minEnergy, maxEnergy = 1e20):
    E, y, err_stat_lo, err_stat_up, err_sys_lo, err_sys_up = np.loadtxt(filename,usecols=(0,1,2,3,4,5),unpack=True)
    items = [i for i in range(len(E)) if (E[i] > minEnergy and E[i] < maxEnergy)]
    #err_stat_lo = np.sqrt(err_stat_lo * err_stat_lo + err_sys_lo * err_sys_lo)
    #err_stat_up = np.sqrt(err_stat_up * err_stat_up + err_sys_up * err_sys_up)
    return E[items], y[items], err_stat_lo[items], err_stat_up[items]

def fit_bpl(params : list, filename: str, doFixS: bool):
    xd, yd, errd_lo, errd_up = load_data(filename, 1e2, 1e4)
    def chi2_function(I0, E0, alpha, logEb, dalpha, s):
        chi2 = 0.
        for x_i, y_i, err_lo_i, err_up_i in zip(xd, yd, errd_lo, errd_up):
            m = bpl(x_i, [I0, E0, alpha, logEb, dalpha, s])
            if m > y_i:
                chi2 += chi2_single(m, y_i, err_up_i)
            else:
                chi2 += chi2_single(m, y_i, err_lo_i)
        return chi2

    I0, E0, alpha, logEb, dalpha, s = params

    m = Minuit(chi2_function, I0=I0, E0=E0, alpha=alpha, logEb=logEb, dalpha=dalpha, s=s)
    m.errordef = Minuit.LEAST_SQUARES
    
    m.fixed['E0'] = True
    if doFixS:
        m.fixed['s'] = True

    m.limits['alpha'] = (2.9, 3.3)
#    m.limits['I0'] = (10., 50.)
    m.limits['logEb'] = (1., 4.)
#    m.limits['dalpha'] = (0.01, 0.2)
#    m.limits['s'] = (1., 20.)

    m.simplex()
    m.migrad()
    m.hesse()

    dof = len(xd) - m.nfit - 1

    #print(m)

    return m.values, m.errors, m.fval, dof
    
def fit_cutoff(params : list, filename: str):
    xd, yd, errd_lo, errd_up = load_data(filename, 1e2, 1e4)
    def chi2_function(I0, E0, alpha, logEb):
        chi2 = 0.
        for x_i, y_i, err_lo_i, err_up_i in zip(xd, yd, errd_lo, errd_up):
            m = cutoff(x_i, [I0, E0, alpha, logEb])
            if m > y_i:
                chi2 += chi2_single(m, y_i, err_up_i)
            else:
                chi2 += chi2_single(m, y_i, err_lo_i)
        return chi2

    I0, E0, alpha, logEb = params

    m = Minuit(chi2_function, I0=I0, E0=E0, alpha=alpha, logEb=logEb)
    m.errordef = Minuit.LEAST_SQUARES
    
    m.fixed['E0'] = True
    m.limits['alpha'] = (2.9, 3.3)
    m.limits['logEb'] = (1., 4.)

    m.simplex()
    m.migrad()
    m.hesse()

    dof = len(xd) - m.nfit - 1

    #print(m)

    return m.values, m.errors, m.fval, dof

def fit_bpl_experiment(datafilename : str, paramsfilename: str, doFixS: bool):
    I0, E0, alpha, logEb, dalpha, s = 0.291, 80.0, 3.113, 2.882, 0.827, 5.0
    values, errors, fval, dof = fit_bpl([I0, E0, alpha, logEb, dalpha, s], datafilename, doFixS)
    f = open(paramsfilename, 'w')
    for x,e in zip(values, errors):
        print(f'{x:6.3f} {e:6.3f}')
        f.write(f'{x:6.3f} {e:6.3f}\n')
    print(f'{fval:4.1f} / {dof}')
    f.close()
    print(f'dump params on {paramsfilename}')

def fit_cutoff_experiment(datafilename : str, paramsfilename: str):
    I0, E0, alpha, logEb = 0.323, 80.0, 3.089, 3.029 
    values, errors, fval, dof = fit_cutoff([I0, E0, alpha, logEb], datafilename)
    f = open(paramsfilename, 'w')
    for x,e in zip(values, errors):
        print(f'{x:6.3f} {e:6.3f}')
        f.write(f'{x:6.3f} {e:6.3f}\n')
    print(f'{fval:4.1f} / {dof}')
    f.close()
    print(f'dump params on {paramsfilename}')
                
if __name__== "__main__":
    print('\n CALET')
    fit_bpl_experiment('data/CALET_e-e+_energy.txt', 'CALET_bpl_fit.txt', False)
    fit_cutoff_experiment('data/CALET_e-e+_energy.txt', 'CALET_cutoff_fit.txt')

    print('\n DAMPE')
    fit_bpl_experiment('data/DAMPE_e-e+_energy.txt', 'DAMPE_bpl_fit.txt', True)
    fit_cutoff_experiment('data/DAMPE_e-e+_energy.txt', 'DAMPE_cutoff_fit.txt')
