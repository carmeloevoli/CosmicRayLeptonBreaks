import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('gryphon.mplstyle')
import numpy as np

def savefig(plt, plotname):
    print (plotname)
    plt.savefig(plotname)
    
def bpl(E, params):
    I0, E0, alpha, Eb, dalpha, s = params
    y = (I0 / 1e3) * np.power(E / E0, -alpha)
    y /= np.power(1. + np.power(E / Eb, s), dalpha / s)
    return y
    
def spl(E, params):
    I0, E0, alpha = params
    y = (I0 / 1e3) * np.power(E / E0, -alpha)
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

def plot_direct_leptons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([8e1, 8e3])
        ax.set_ylabel(r'E$^{3}$ I [GeV$^{2}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_yscale('log')
        ax.set_ylim([8, 200])

    fig = plt.figure(figsize=(11.0, 8.0))
    ax = fig.add_subplot(111)
    set_axes(ax)

    plot_data(ax, 'data/DAMPE_e-e+_energy.txt', 3.0, 1., 'o', 'r', r'DAMPE', 3)
    plot_data(ax, 'data/CALET_e-e+_energy.txt', 3.0, 1., 'o', 'g', r'CALET', 3)
    plot_data(ax, 'data/AMS-02_e-e+_energy.txt', 3.0, 1., 'o', 'tab:orange', r'AMS-02', 2)
    plot_data(ax, 'data/FERMI_e-e+_energy.txt', 3.0, 1., 'o', 'b', r'FERMI', 2)
   
    ax.legend(fontsize=25, loc='lower left')
    savefig(plt, 'direct-leptons.pdf')

def plot_DAMPE_leptons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([8e1, 8e3])
        ax.set_ylabel(r'E$^{3}$ I [GeV$^{2}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_yscale('log')
        ax.set_ylim([8, 200])

    fig = plt.figure(figsize=(11.0, 8.0))
    ax = fig.add_subplot(111)
    set_axes(ax)

    plot_data(ax, 'data/DAMPE_e-e+_energy.txt', 3.0, 1., 'o', 'r', r'DAMPE', 3)

    E = np.logspace(1, 4, 1000)
    E3 = np.power(E, 3.)

    I0, E0, alpha, Eb, dalpha, s = 0.3231, 80.0, 3.089, 1.07e3, 1.2, 5
    ax.plot(E, E3 * bpl(E, [I0, E0, alpha, Eb, dalpha, s]), zorder=9, color='r', label='$\chi^2$/dof = 29 / 32')

#    ax.vlines(1e3, 0, 1000, ls=':')
#    ax.fill_between([10, 20], 10, 230, color='tab:gray', alpha=0.2)

    ax.legend(fontsize=25, loc='lower left')
    savefig(plt, 'DAMPE-leptons.pdf')

def plot_indirect_leptons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([8e1, 8e3])
        ax.set_ylabel(r'E$^{3}$ I [GeV$^{2}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_yscale('log')
        #ax.set_ylim([8, 200])
        #ax.ticklabel_format(axis='y', style='sci', scilimits=(3,3))

    fig = plt.figure(figsize=(11.0, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    plot_data(ax, 'data/VERITAS_e-e+_energy.txt', 3.0, 1., 'o', 'r', r'VERITAS', 3)
    plot_data(ax, 'data/HESS_e-e+_energy.txt', 3.0, 1., 'o', 'g', r'HESS', 3)

    ax.legend(fontsize=25, loc='lower left')
    savefig(plt, 'indirect-leptons.pdf')
    
if __name__== "__main__":
    plot_direct_leptons()
    plot_DAMPE_leptons()
    plot_indirect_leptons()
