import numpy as np
from iminuit import Minuit

def spl(E, params):
    I0, E0, alpha = params
    y = (I0 * 1e2) * np.power(E / E0, alpha)
    return y
    
def chi2_single(x, mu, sigma):
    return np.power((x - mu) / sigma, 2.)

def load_data(minEnergy, maxEnergy = 1e20):
    filename = 'data/AMS-02_H_over_e-_energy.txt'
    E, y, err_stat_lo, err_stat_up, err_sys_lo, err_sys_up = np.loadtxt(filename,usecols=(0,1,2,3,4,5),unpack=True)
    items = [i for i in range(len(E)) if (E[i] > minEnergy and E[i] < maxEnergy)]
    err_stat_lo = np.sqrt(err_stat_lo * err_stat_lo + err_sys_lo * err_sys_lo)
    err_stat_up = np.sqrt(err_stat_up * err_stat_up + err_sys_up * err_sys_up)
    return E[items], y[items], err_stat_lo[items], err_stat_up[items]

def fit_spl(params):
    xd, yd, errd_lo, errd_up = load_data(40., 500.)

    def chi2_function(I0, E0, alpha):
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
    #m.limits['alpha'] = (3.1, 3.4)
    #m.limits['I0'] = (10., 50.)

    m.simplex()
    m.migrad()
    m.hesse()

    #print(m)
    print('Single power-law:')
    print(f'params : {m.values[0]:5.3f} {m.values[1]:5.0f} {m.values[2]:5.3f}')
    print(f'chi2 / dof : {m.fval:5.0f} / {len(xd) - m.nfit - 1}')

    return m.values, m.errors, m.fval
    
if __name__== "__main__":
    I0, E0, alpha = 3., 1e2, 0.35
    
    values, errors, fval = fit_spl([I0, E0, alpha])
