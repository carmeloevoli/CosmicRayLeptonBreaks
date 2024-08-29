import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('gryphon.mplstyle')
import numpy as np

def savefig(plt, plotname):
    print (plotname)
    plt.savefig(plotname)
        
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

def plot_positronfraction():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1, 2e3])
        ax.set_ylabel(r'$e^+$/($e^+$+$e^-$)')
        ax.set_ylim([0.001, 0.28])

    fig = plt.figure(figsize=(10.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    plot_data(ax, 'data/AMS-02_pf_energy.txt', 0., 1., 'o', 'r', r'AMS-02', 4)
    plot_data(ax, 'data/FERMI_pf_energy.txt', 0., 1., 'o', 'b', r'FERMI', 3)   
    plot_data(ax, 'data/PAMELA_pf_energy.txt', 0., 1., 'o', 'g', r'PAMELA', 5)

    E = np.logspace(0, 4, 1000)
    y = 0.115 * np.power(E, -0.5)

    ax.fill_between(E, 0., y, color='tab:gray', alpha=0.2, zorder=1)
    ax.plot(E, y, color='tab:gray')
    ax.text(1.8, 0.1, r'$E^{-\delta}$', color='dimgrey', fontsize=26)

    ax.legend(fontsize=20, loc='upper left')
    savefig(plt, 'TeVPA24-positron-fraction.pdf')

def plot_E2positrons():   
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([10, 1e3])
        #ax.set_ylabel(r'$e^+$/($e^+$+$e^-$)')
        ax.set_ylim([0, 0.5])

    fig = plt.figure(figsize=(10.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    plot_data(ax, 'data/AMS-02_e+_energy.txt', 2., 1., 'o', 'r', r'AMS-02', 4)
    ax.vlines(1e2, 0, 0.5, ls=':')
    ax.hlines(0.2, 10, 1e3, ls=':')
    ax.legend(fontsize=20, loc='upper left')
    savefig(plt, 'TeVPA24-E2positrons.pdf')

if __name__== "__main__":
    plot_positronfraction()
    #plot_E2positrons()
