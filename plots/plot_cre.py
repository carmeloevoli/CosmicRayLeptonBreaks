import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('gryphon.mplstyle')
import numpy as np

def savefig(plt, plotname):
    print (plotname)
    plt.savefig(plotname)

def bpl(E, params):
    I0, E0, alpha, logEb, dalpha, s = params
    Eb = np.power(10., logEb)
    y = (I0 / 1e3) * np.power(E / E0, -alpha)
    y /= np.power(1. + np.power(E / Eb, s), dalpha / s)
    return y

def plot_timescales():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1e1, 1e4])
        ax.set_ylabel(r'timescale [Myr]')
        ax.set_yscale('log')
        ax.set_ylim([0.1, 1e2])

    fig = plt.figure(figsize=(10.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    filename = '../model/TeVPA24_timescales.txt'
    E, tau_d5, tau_d2 = np.loadtxt(filename, usecols=(0,1,2), unpack=True)
    ax.fill_between(E, tau_d2, tau_d5, color='tab:gray', alpha=0.33, zorder=1, label=r'Escape [$2 < H < 5$ kpc]')

    E, tau_B1, tau_B3 = np.loadtxt(filename, usecols=(0,8,9), unpack=True) 
    ax.fill_between(E, tau_B3, tau_B1, color='tab:orange', alpha=0.3, zorder=2, label=r'$1 < B < 3 \mu$G')

    E, tau_CMB, tau_ISRF = np.loadtxt(filename, usecols=(0,3,7), unpack=True) 
    ax.plot(E, tau_CMB, color='r', zorder=3, label='CMB')
    ax.plot(E, tau_ISRF, color='g', zorder=3, label='ISRF')

    tau_all = 1. / (1. / tau_CMB + 1. / tau_ISRF + 1. / tau_B3)
    ax.plot(E, tau_all, ls='--', color='k', label=r'$\tau_{\rm loss}$')
    
    #ax.plot(E, 22. * np.power(E / 10., -0.88))

    ax.legend(fontsize=19, loc='lower left')
    savefig(plt, 'TeVPA24-timescales.pdf')

def plot_horizon():
    def set_axes_top(ax):
        ax.set_xscale('log')
        ax.set_xlim([1e1, 1e4])
        ax.set_ylabel(r'diffusion length [kpc]')
        ax.set_ylim([1, 9])
        ax.set_xticklabels([])

    def set_axes_bottom(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1e1, 1e4])
        ax.set_ylabel(r'number of sources [log$_{10}$]')
        ax.set_ylim([1, 6])
        
    fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(11., 14.))
    ax1, ax2 = axes
    set_axes_top(ax1)
    set_axes_bottom(ax2)

    filename = '../model/TeVPA24_horizons.txt'
    E, l2_1, l2_3 = np.loadtxt(filename, usecols=(0,1,2), unpack=True)
    
    #l2 = lambda_2(E, 1.)
    ax1.fill_between(E, np.sqrt(l2_1), 10., color='tab:blue', alpha=0.1)
    ax1.plot(E, np.sqrt(l2_1), color='tab:blue', lw=4)
    ax1.fill_between(E, np.sqrt(l2_3), 10., color='tab:blue', alpha=0.2)
    
    ax1.vlines(1e3, 1, 10, ls='--', lw=2, color='tab:gray')
    ax1.hlines(3.3, 1e1, 1e4, ls='--', lw=2, color='tab:gray')
    ax1.text(0.8e2, 2.5, r'$\lambda (\text{TeV}) \simeq 3$ kpc')
    ax1.text(15., 5.2, r'$\delta \simeq 0.54$', fontsize=26)
    ax1.text(15., 4.5, r'$H \simeq 5$ kpc', fontsize=26)

    #ax.text(1.5e3, 3.6, r'$B \simeq 1 \, \mu$G', fontsize=22)
    ax1.text(1.5e3, 2.0, r'$B \simeq 3 \, \mu$G', fontsize=23)

    filename = '../model/TeVPA24_horizons.txt'
    E, N_1, N_3 = np.loadtxt(filename, usecols=(0,3,4), unpack=True)

    ax2.fill_between(E, np.log10(N_1), 10., color='tab:red', alpha=0.1)
    ax2.plot(E, np.log10(N_1), color='tab:red', lw=4)
    ax2.fill_between(E, np.log10(N_3), 10., color='tab:red', alpha=0.2)

    ax2.text(15., 3.5, r'$\mathcal R$ = 2 / century', fontsize=26)

    ax2.vlines(1e3, 1, 10, ls='--', lw=2, color='tab:gray')
    ax2.hlines(2.7, 1e1, 1e4, ls='--', lw=2, color='tab:gray')
    ax2.text(0.6e2, 2.3, 'N(TeV) $\sim \mathcal O(10^3)$')

    savefig(plt, 'TeVPA24-horizon.pdf')

def plot_sourceterm():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1e1, 1e4])
        ax.set_yscale('log')
        ax.set_ylabel(r'q(E)')
        #ax.set_ylim([0.1, 2])

    fig = plt.figure(figsize=(10.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    params, errors = np.loadtxt('../model/CALET_bpl_fit.txt', usecols=(0,1), unpack=True)

    E = np.logspace(1, 4, 1000)
    y = np.power(E, 3.1) * bpl(E, params)

    ax.plot(E, y / max(y), zorder=9, color='tab:blue', label='CALET fit')

    # filename = '../model/TeVPA24_sources.txt'
    # E, Q_1, Q_2, Q_3 = np.loadtxt(filename, usecols=(0,1,2,3), unpack=True)

    # #ax.plot(E, Q_3, color='tab:gray')

    filename = '../model/TeVPA24_cre.txt'
    E, n_1, n_2, n_3 = np.loadtxt(filename, usecols=(0,1,2,3), unpack=True)

    y = np.power(E, 3.1) * n_3
    
    ax.plot(E, y / max(y), color='tab:red')

    ax.legend(fontsize=19, loc='lower left')
    savefig(plt, 'TeVPA24-sourceterms.pdf')

if __name__== "__main__":
    plot_timescales()
    plot_horizon()
    #plot_sourceterm()
