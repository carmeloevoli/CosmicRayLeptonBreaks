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
    
def spl(E, params):
    I0, E0, alpha = params
    y = (I0 / 1e3) * np.power(E / E0, -alpha)
    return y

def cutoff(E, params):
    I0, E0, alpha, logEb = params
    Eb = np.power(10., logEb)
    y = (I0 / 1e3) * np.power(E / E0, -alpha)
    y *= np.exp(- E / Eb)
    return y
    
def plot_data(ax, filename, slope, norm, fmt, color, label, zorder=3):
    E, y, err_stat_lo, err_stat_up, err_sys_lo, err_sys_up = np.loadtxt(filename,usecols=(0,1,2,3,4,5),unpack=True)
    y = norm * np.power(E, slope) * y
    y_err_lo = norm * np.power(E, slope) * err_stat_lo
    y_err_up = norm * np.power(E, slope) * err_stat_up
    ax.errorbar(E, y, yerr=[y_err_lo, y_err_up], fmt=fmt, markeredgecolor=color, color=color, label=label,
                capsize=3.5, markersize=6, elinewidth=1.8, capthick=1.8, zorder=zorder)
    y_err_lo = norm * np.power(E, slope) * err_sys_lo
    y_err_up = norm * np.power(E, slope) * err_sys_up
    size = len(E)
    for i in range(size):
        ax.fill_between([0.96 * E[i], 1.04 * E[i]], y[i] - y_err_lo[i], y[i] + y_err_up[i], color='tab:gray', alpha=0.28, zorder=1)
#    ax.fill_between(E, y - y_err_lo, y + y_err_up, color='tab:purple', alpha=0.25, zorder=1)
#    ax.errorbar(E, y, yerr=[y_err_lo, y_err_up], fmt=fmt, markeredgecolor=color, color=color,
#                capsize=3.5, markersize=6, elinewidth=1.8, capthick=1.8, zorder=zorder)

def plot_nosys_data(ax, filename, slope, norm, fmt, color, label, zorder=3):
    E, y, err_stat_lo, err_stat_up, err_sys_lo, err_sys_up = np.loadtxt(filename,usecols=(0,1,2,3,4,5),unpack=True)
    y = norm * np.power(E, slope) * y
    y_err_lo = norm * np.power(E, slope) * err_stat_lo
    y_err_up = norm * np.power(E, slope) * err_stat_up
    ax.errorbar(E, y, yerr=[y_err_lo, y_err_up], fmt=fmt, markeredgecolor=color, color=color, label=label,
                capsize=3.5, markersize=6, elinewidth=1.8, capthick=1.8, zorder=zorder)

def plot_data_leptons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([6, 1e4])
        ax.set_ylabel(r'E$^{3}$ I [GeV$^{2}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        #ax.set_yscale('log')
        ax.set_ylim([0., 280])

    fig = plt.figure(figsize=(10.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    plot_data(ax, 'data/DAMPE_e-e+_energy.txt', 3.0, 1., 'o', 'r', r'DAMPE 2017', 4)
    plot_data(ax, 'data/CALET_e-e+_energy.txt', 3.0, 1., 'o', 'b', r'CALET 2024', 4)
    
    plot_data(ax, 'data/AMS-02_e-e+_energy.txt', 3.0, 1., 'o', 'tab:green', r'AMS-02 2021', 3)
    
    plot_data(ax, 'data/FERMI_e-e+_energy.txt', 3.0, 1., 'o', 'tab:olive', r'FERMI 2017', 2)
   
    plot_nosys_data(ax, 'data/VERITAS_e-e+_energy.txt', 3.0, 1., 's', 'tab:gray', r'VERITAS (only stat) 2018', 1)
    plot_nosys_data(ax, 'data/HESS_e-e+_energy.txt', 3.0, 1., 's', 'tab:orange', r'HESS (only stat) 2008', 1)
    plot_nosys_data(ax, 'data/HESS-LE_e-e+_energy.txt', 3.0, 1., 's', 'tab:orange', r'HESS-LE (only stat) 2009', 1)

    ax.annotate("", xy=(1e3, 190.), xytext=(1e3, 210.), arrowprops=dict(arrowstyle="->", color='blue', lw=5.))
    #ax.annotate("", xy=(50, 190.), xytext=(50, 210.), arrowprops=dict(arrowstyle="->", color='blue', lw=5.))

    ax.legend(fontsize=17, loc='lower left')
    savefig(plt, 'TeVPA24-data-leptons-HE.pdf')

def plot_DAMPE_leptons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([5e1, 1e4])
        ax.set_ylabel(r'E$^{3.1}$ I [GeV$^{2.1}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_yscale('log')
        ax.set_ylim([20, 350])

    fig = plt.figure(figsize=(10.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    E = np.logspace(1, 4, 1000)
    E3 = np.power(E, 3.1)

    params, errors = np.loadtxt('../model/DAMPE_bpl_fit.txt', usecols=(0,1), unpack=True)
    ax.plot(E, E3 * bpl(E, params), zorder=9, color='tab:purple', label='$\chi^2$/dof = 36 / 23')

    logEb, sLogEb = params[3], errors[3]
    ax.vlines(np.power(10., logEb), 1, 1e3, lw=2, color='tab:gray', zorder=1)
    ax.fill_between([np.power(10., logEb - sLogEb), np.power(10., logEb + sLogEb)], 1, 1e3, alpha=.3, color='tab:gray', zorder=1)

    params, errors = np.loadtxt('../model/DAMPE_cutoff_fit.txt', usecols=(0,1), unpack=True)
    ax.plot(E, E3 * cutoff(E, params), zorder=9, ls=':', color='tab:purple', label='$\chi^2$/dof = 49 / 24')

    ax.legend(fontsize=24, loc='lower left')

    plot_data(ax, 'data/DAMPE_e-e+_energy.txt', 3.1, 1., 'o', 'g', r'DAMPE', 3)

    ax.text(60., 170., 'DAMPE', color='g', fontsize=28)
    savefig(plt, 'TeVPA24-DAMPE-leptons.pdf')

def plot_CALET_leptons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([5e1, 1e4])
        ax.set_ylabel(r'E$^{3.1}$ I [GeV$^{2.1}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_yscale('log')
        ax.set_ylim([20, 350])

    fig = plt.figure(figsize=(10.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    E = np.logspace(1, 4, 1000)
    E3 = np.power(E, 3.1)

    params, errors = np.loadtxt('../model/CALET_bpl_fit.txt', usecols=(0,1), unpack=True)
    ax.plot(E, E3 * bpl(E, params), zorder=9, color='tab:blue', label='$\chi^2$/dof = 14.7 / 16')

    logEb, sLogEb = params[3], errors[3]
    ax.vlines(np.power(10., logEb), 1, 1e3, lw=2, color='tab:gray', zorder=1)
    ax.fill_between([np.power(10., logEb - sLogEb), np.power(10., logEb + sLogEb)], 1, 1e3, alpha=.3, color='tab:gray', zorder=1)

    params, errors = np.loadtxt('../model/CALET_cutoff_fit.txt', usecols=(0,1), unpack=True)
    ax.plot(E, E3 * cutoff(E, params), zorder=9, ls=':', color='tab:blue', label='$\chi^2$/dof = 24.3 / 18')

    ax.legend(fontsize=24, loc='lower left')

    plot_data(ax, 'data/CALET_e-e+_energy.txt', 3.1, 1., 'o', 'r', r'DAMPE', 3)

    ax.text(60., 170., 'CALET', color='r', fontsize=28)
    savefig(plt, 'TeVPA24-CALET-leptons.pdf')

def plot_CALET_theory():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([5e1, 3e3])
        ax.set_ylabel(r'E$^{3.1}$ I [GeV$^{2.1}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        #ax.set_yscale('log')
        ax.set_ylim([70, 270])

    fig = plt.figure(figsize=(10.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    E = np.logspace(1, 4, 1000)
    E3 = np.power(E, 3.1)

    params, errors = np.loadtxt('../model/CALET_bpl_fit.txt', usecols=(0,1), unpack=True)
    ax.plot(E, E3 * bpl(E, params), zorder=9, color='tab:blue', ls='--', label='fit')

    filename = '../model/TeVPA24_cre.txt'
    E, n_1, n_2, n_3 = np.loadtxt(filename, usecols=(0,1,2,3), unpack=True)
    y = np.power(E, 3.1) * n_3    
    ax.plot(E, 2.49e2 * y / max(y), color='tab:blue', zorder=5, label=r'$n_{\rm cre}$')

    plot_data(ax, 'data/CALET_e-e+_energy.txt', 3.1, 1., 'o', 'r', r'CALET', 1)
    ax.legend(fontsize=24, loc='lower left')

    #ax.text(60., 170., 'CALET', color='r', fontsize=28)
    savefig(plt, 'TeVPA24-CALET-theory.pdf')

if __name__== "__main__":
    plot_data_leptons()
    plot_DAMPE_leptons()
    plot_CALET_leptons()
    plot_CALET_theory()
    
