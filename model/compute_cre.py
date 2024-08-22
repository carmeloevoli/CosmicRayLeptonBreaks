import numpy as np
import math

kpc2cm = 3.0857e21
Myr2sec = 3.155815e13
GeV2erg = 624.151

# constants

electronMassC2 = 0.511e-3 # GeV
cLight = 2.99792458e10 / kpc2cm * Myr2sec # kpc/Myr
sigmaTh = 6.66524e-25 / kpc2cm**2. # kpc2
kBoltzmann = 8.617333262e-14 # GeV / K

# diffusion parameters

D_0_over_H = 0.44e28 / kpc2cm**2. * Myr2sec # kpc/Myr
E_0 = 1. # GeV
D_delta = 0.54 # none
D_E_break = 312. # GeV
D_ddelta = 0.2 # none
D_smoothness = 0.1 # none

# Galaxy

R_Galaxy = 15. # kpc
rate_sources = 1. / 50. * 1e6 # Myr-1

# ISRF

CMB = [2.7, 0.260e-9 * kpc2cm**3.0] # K, GeV / kpc3
IR = [33.07, 25.4e-11* kpc2cm**3.0] # K, GeV / kpc3
opt = [313.32,  5.47e-11 * kpc2cm**3.0] # K, GeV / kpc3
UVI = [313.32,  5.47e-11 * kpc2cm**3.0] # K, GeV / kpc3
UVII = [6150.4, 22.9e-11 * kpc2cm**3.0] # K, GeV / kpc3
UVIII = [23209.0, 11.89e-11 * kpc2cm**3.0] # K, GeV / kpc3

# Quantities

def tau_diffusion(E: float, H: float) -> float:
    D_0 = D_0_over_H * H
    D = D_0 * np.power(E / E_0, D_delta)
    D /= np.power(1. + np.power(E / D_E_break, D_ddelta / D_smoothness), D_smoothness)
    return (H * H) / 2. / D 

def dEdt_IC(E: float, T: float, energy_density: float):
    def Y(x: float) -> float:
        value = 0.
        c = [-3.996e-2, -9.100e-1, -1.197e-1, 3.305e-3, 1.044e-3, -7.013e-5, -9.618e-6]
        if (x < 1.5e-3):
            value = np.power(math.pi, 4.) / 15.
        elif (x < 150.):
            c_log = 0.
            for i in range(7):
                c_log += c[i] * np.power(np.log(x), float(i))
            value = np.exp(c_log)
        else:
            value = 3. / 4. * np.power(math.pi / x, 2.) * (np.log(x) - 1.9805)
        return value

    gamma_e = E / electronMassC2 # none
    factor = 20. * cLight * sigmaTh / np.power(math.pi, 4.) # kpc3/Myr
    S_i = Y(4. * gamma_e * kBoltzmann * T / electronMassC2) # none
    return factor * energy_density * gamma_e**2 * S_i # GeV/Myr

def tau_IC(E: float, phField: list) -> float: 
    return E / dEdt_IC(E, phField[0], phField[1]) # Myr

def dEdt_sync(E: float, B: float) -> float:
    gamma_e = E / electronMassC2
    factor = 3. / 4. * cLight * sigmaTh
    energy_density = B * B / 8. / math.pi * GeV2erg * kpc2cm**3.0
    return factor * energy_density * gamma_e**2 # GeV/Myr

def tau_sync(E: float, B: float) -> float: 
    return E / dEdt_sync(E, B) # Myr

def lambda_2(E: float, B: float) -> float:
    b_0 = dEdt_IC(E_0, CMB[0], CMB[1]) # CMB
    b_0 += dEdt_IC(E_0, IR[0], IR[1]) # IR
    b_0 += dEdt_IC(E_0, opt[0], opt[1]) # optical
    b_0 += dEdt_IC(E_0, UVI[0], UVI[1]) # UVI
    b_0 += dEdt_sync(E_0, B) # synchrotron

    D_0 = D_0_over_H * 5. # kp2/Myr
    l2 = 4. * D_0 * E_0 / b_0 / (1. - D_delta) # kpc2
    l2 *= np.power(E / E_0, D_delta - 1.)
    return l2

def nsources(E: float, B: float) -> float:
    l2 = lambda_2(E, B)
    t_l = 1. / (1. / tau_IC(E, CMB) +  1. / tau_IC(E, IR) + 1. / tau_IC(E, opt) + 1. / tau_IC(E, UVI) + 1. / tau_sync(E, B))
    x = min(1., l2 / R_Galaxy**2.0)
    return t_l * rate_sources * x

def lambda_prop(E: float, Eprime: float, delta: float) -> float:
    assert(Eprime >= E)
    assert(delta < 1.)
    l2 = np.power(E / E_0, delta - 1.) - np.power(Eprime / E_0, delta - 1.)
    return np.sqrt(l2)

def qconst(E: float) -> float:
    return 1.

def qcutoff(E: float) -> float:
    return np.exp(- (E / 2e4)**2.0)

def qbreak(E: float) -> float:
    E_break = 0.90e3
    alpha = 2.46
    dalpha = 0.75
    s = 30.
    I0, E0 = 1., 1e2
    y = (I0 / 1e3) * np.power(E / E0, -alpha)
    y /= np.power(1. + np.power(E / E_break, s), dalpha / s)
    return y

def compute_cre(E: float) -> float:
    import scipy.integrate as integrate
    def integrand(lnx: float) -> float:
        x = np.exp(lnx)
        return x * qbreak(x) / lambda_prop(E, x, delta = 0.34)
    ne = integrate.quad(lambda lnx: integrand(lnx), np.log(E), np.log(E) + 20.)
    return ne / np.power(E / E_0, 2.)

def print_timescale(filename: str):
    energies = np.logspace(1, 5, 2000)
    with open(filename, 'w') as f:
        f.write('# E [Myr] - tau_esc(5kpc) - tau_esc(2kpc) - CMB - IR - opt - UV - ISRF - B1muG - B2muG [Myr]\n')
        for E in energies:
            t_d5 = tau_diffusion(E, 5.)
            t_d2 = tau_diffusion(E, 2.)
            t_CMB = tau_IC(E, CMB)
            t_IR = tau_IC(E, IR)
            t_opt = tau_IC(E, opt)
            t_UV = 1. / (1. / tau_IC(E, UVI) + 1. / tau_IC(E, UVII) + 1. / tau_IC(E, UVIII)) 
            t_ISRF = 1. / (1. / t_IR + 1. / t_opt + 1. / t_UV)
            t_B1 = tau_sync(E, 1e-6)
            t_B3 = tau_sync(E, 3e-6)
            f.write(f'{E:5.3e} {t_d5:5.3e} {t_d2:5.3e} {t_CMB:5.3e} {t_IR:5.3e} {t_opt:5.3e} {t_UV:5.3e} {t_ISRF:5.3e} {t_B1:5.3e} {t_B3:5.3e}\n')  
    f.close()
    print(f'dump on {filename}')

def print_horizon(filename: str):
    energies = np.logspace(1, 5, 2000)
    with open(filename, 'w') as f:
        f.write('# E [GeV] - l2(kpc2)\n')
        for E in energies:
            l2_1 = lambda_2(E, 1e-6)
            l2_3 = lambda_2(E, 3e-6)
            N_1 = nsources(E, 1e-6)
            N_3 = nsources(E, 3e-6)
            f.write(f'{E:5.3e} {l2_1:5.3e} {l2_3:5.3e} {N_1:5.3e} {N_3:5.3e}\n')
    f.close()
    print(f'dump on {filename}')

def print_sourceterms(filename: str):
    energies = np.logspace(1, 4, 300)
    with open(filename, 'w') as f:
        f.write('# E [GeV] - Q\n')
        for E in energies:
            q_1 = qconst(E)
            q_2 = qcutoff(E)
            q_3 = qbreak(E)
            f.write(f'{E:5.3e} {q_1:5.3e} {q_2:5.3e} {q_3:5.3e}\n')
    f.close()
    print(f'dump on {filename}')

def print_cre(filename: str):
    energies = np.logspace(1, 4, 1000)
    with open(filename, 'w') as f:
        f.write('# E [GeV] - Q\n')
        for E in energies:
            n_3 = compute_cre(E)
            f.write(f'{E:5.3e} 0 0 {n_3[0]:5.3e}\n')
    f.close()
    print(f'dump on {filename}')

if __name__== "__main__":
    #print_timescale('TeVPA24_timescales.txt')
    #print_horizon('TeVPA24_horizons.txt')
    #print_sourceterms('TeVPA24_sources.txt')
    print_cre('TeVPA24_cre.txt')
