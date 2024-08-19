import numpy as np

def compute_mean_energy(E_min, E_max, slope):
    x_ratio = E_min / (E_max - E_min)
    x_tilde = E_min
    x_tilde *= np.power(x_ratio / (slope - 1.) * (1. - np.power(E_max / E_min, -slope + 1.)), -1. / slope)
    return x_tilde

def quadrature(a, b):
    return np.sqrt(a * a + b * b)
    
def readfile(filename):
    x_min, x_max, y, yStaLo, yStaUp, ySysLo, ySysUp = np.loadtxt(filename, usecols=(0,1,2,3,4,5,6), unpack=True)
    return x_min, x_max, y, yStaLo, yStaUp, ySysLo, ySysUp

def dump(data, filename):
    print(filename)
    f = open('output/' + filename, 'w')
    E_mean, I_E, IStaLo, IStaUp, ISysLo, ISysUp = data
    for i,j,k,l,m,n in zip(E_mean, I_E, IStaLo, IStaUp, ISysLo, ISysUp):
        f.write(f'{i:6.4e} {j:6.4e} {k:6.4e} {l:6.4e} {m:6.4e} {n:6.4e}\n')
    f.close()

def transform_AMS02():
    size = 73
    R_min, R_max, I_R, eStaLo, eStaUp, eSysLo, eSysUp = readfile('lake/AMS-02_e+_rigidity.txt')
    R_mean = compute_mean_energy(R_min, R_max, 3.0)
    data = [R_mean[0:size], I_R[0:size], eStaLo[0:size], eStaUp[0:size], eSysLo[0:size], eSysUp[0:size]]
    dump(data, 'AMS-02_e+_rigidity.txt')

    R_min, R_max, I_R, eStaLo, eStaUp, eSysLo, eSysUp = readfile('lake/AMS-02_e-_rigidity.txt')
    R_mean = compute_mean_energy(R_min, R_max, 3.0)
    data = [R_mean[0:size], I_R[0:size], eStaLo[0:size], eStaUp[0:size], eSysLo[0:size], eSysUp[0:size]]
    dump(data, 'AMS-02_e-_rigidity.txt')

    R_min, R_max, I_e, eStaLo_e, eStaUp_e, eSysLo_e, eSysUp_e = readfile('lake/AMS-02_e-_rigidity.txt')
    R_min, R_max, I_p, eStaLo_p, eStaUp_p, eSysLo_p, eSysUp_p = readfile('lake/AMS-02_e+_rigidity.txt')
    R_mean = compute_mean_energy(R_min[0:size], R_max[0:size], 3.0)
    y = I_e[0:size] - I_p[0:size]
    estaLo = quadrature(eStaLo_e[0:size], eStaLo_p[0:size])
    estaUp = quadrature(eStaUp_e[0:size], eStaUp_p[0:size])
    eSysLo = eSysLo_e[0:size] # + eSysLo_p[0:size]
    eSysUp = eSysUp_e[0:size] # + eSysUp_p[0:size]
    data = [R_mean, y, eStaLo, eStaUp, eSysLo, eSysUp]
    dump(data, 'AMS-02_e-_minus_e+_rigidity.txt')

    y = I_e[0:size] - (I_p[0:size] + eSysUp_p[0:size])
    data = [R_mean, y, eStaLo, eStaUp, eSysLo, eSysUp]
    dump(data, 'AMS-02_e-_minus_e+_statUp_rigidity.txt')

def transform_AMS02_leptons():
    E_min, E_max, I_E, eStaLo, eStaUp, eSysLo, eSysUp = readfile('lake/AMS-02_e-e+_rigidity.txt')
    E_mean = compute_mean_energy(E_min, E_max, 3.0)
    data = [E_mean, I_E, eStaLo, eStaUp, eSysLo, eSysUp]
    dump(data, 'AMS-02_e-e+_energy.txt')
    
def transform_DAMPE():
    E_min, E_max, I_E, eStaLo, eStaUp, eSysLo, eSysUp = readfile('lake/DAMPE_e-e+_energy.txt')
    E_mean = compute_mean_energy(E_min, E_max, 3.0)
    data = [E_mean, I_E, eStaLo, eStaUp, eSysLo, eSysUp]
    dump(data, 'DAMPE_e-e+_energy.txt')

def transform_CALET():
    E_min, E_max, I_E, eStaLo, eStaUp, eSysLo, eSysUp = readfile('lake/CALET_2311.05916.txt')
    E_mean = compute_mean_energy(E_min, E_max, 3.0)
    data = [E_mean, I_E, eStaLo, eStaUp, eSysLo, eSysUp]
    dump(data, 'CALET_e-e+_energy.txt')

def transform_FERMI():
    E_min, E_max, I_E, eStaLo, eStaUp, eSysLo, eSysUp = readfile('lake/FERMI_e-e+_energy.txt')
    E_mean = compute_mean_energy(E_min, E_max, 3.0)
    data = [E_mean, I_E, eStaLo, eStaUp, eSysLo, eSysUp]
    dump(data, 'FERMI_e-e+_energy.txt')
    
def transform_HESS():
    filename = 'lake/HESS_e+e-_totalEnergy.txt'
    E_mean, I_E, eStaLo, eStaUp = np.loadtxt(filename, usecols=(0,1,2,3), unpack=True)
    data = [E_mean, I_E, eStaLo, eStaUp, 0. * eStaLo, 0. * eStaUp]
    dump(data, 'HESS_e-e+_energy.txt')
    
def transform_VERITAS():
    filename = 'lake/VERITAS_e+e-_totalEnergy.txt'
    E_min, E_max, I_E, eSta = np.loadtxt(filename, usecols=(0,1,2,3), unpack=True)
    E_mean = 1e3 * compute_mean_energy(E_min, E_max, 3.0) # TeV -> GeV
    I_E *= 1e4 #
    eSta *= 1e4 #
    data = [E_mean, I_E, eSta, eSta, 0. * eSta, 0. * eSta]
    dump(data, 'VERITAS_e-e+_energy.txt')

if __name__== "__main__":
    transform_AMS02()
    transform_AMS02_leptons()
    transform_CALET()
    transform_DAMPE()
    transform_FERMI()
    transform_HESS()
    transform_VERITAS()

