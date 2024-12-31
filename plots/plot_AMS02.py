import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('gryphon.mplstyle')
import numpy as np
from utils import SPL, BPL, BACKGROUND

def savefig(plt, plotname):
    print (plotname)
    plt.savefig(plotname)
    
def plot_data(ax, filename, slope, norm, fmt, color, label, zorder=3):
    E, y, err_stat_lo, err_stat_up, err_sys_lo, err_sys_up = np.loadtxt(filename,usecols=(0,1,2,3,4,5),unpack=True)
    y = norm * np.power(E, slope) * y
    y_err_lo = norm * np.power(E, slope) * err_stat_lo
    y_err_up = norm * np.power(E, slope) * err_stat_up
    ax.errorbar(E, y, yerr=[y_err_lo, y_err_up], fmt=fmt, markeredgecolor=color, color=color,
                capsize=3.5, markersize=6, elinewidth=1.8, capthick=1.8, zorder=zorder, label=label)
    y_err_lo = norm * np.power(E, slope) * err_sys_lo
    y_err_up = norm * np.power(E, slope) * err_sys_up
    ax.bar(E, y_err_up + y_err_lo, bottom=y-y_err_lo, width=.06 * E, color='y', zorder=0, align='center')
    #ax.fill_between(E, y - y_err_lo, y + y_err_up, color='tab:purple', alpha=0.25, zorder=1)
 
def plot_electrons_and_positrons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1e1, 7e2])
        ax.set_ylabel(r'E$^{3}$ I [GeV$^{2}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_ylim([60, 230])

    fig = plt.figure(figsize=(10.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    plot_data(ax, 'data/AMS-02_e-_energy.txt', 3.0, 1., 'o', 'tab:gray', r'e$^-$', 3)
    #ax.text(35., 170., 'electrons', color='tab:gray', fontsize=21)

    plot_data(ax, 'data/AMS-02_e-_minus_e+_energy.txt', 3.0, 1., 'o', 'tab:red', r'e$^-$ - e$^+$', 5)
    #ax.text(35., 170., 'electrons', color='tab:red', fontsize=21)

    plot_data(ax, 'data/AMS-02_e+_energy.txt', 3.0, 10., 'o', 'tab:blue', r'e$^+$ [10x]', 6)
    #ax.text(115., 70., 'positrons [5x]', color='tab:blue', fontsize=21)

    ax.text(22., 215., 'AMS-02', color='tab:gray', fontsize=25)

    ax.fill_between([0, 20], 0, 300, color='tab:gray', alpha=0.18, zorder=1)
    
    ax.legend(fontsize=24, loc='best')
    savefig(plt, 'AMS02-electrons-and-positrons.pdf')

def plot_electrons_minus_positrons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([2e1, 7e2])
        ax.set_ylabel(r'E$^{3.3}$ I [GeV$^{2.3}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_ylim([400, 540])

    fig = plt.figure(figsize=(10.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    ax.text(100., 420., r'$\sigma^{\rm sta}_{e^- - e^+} = \sqrt{(\sigma^{\rm sta}_{e^-})^2 + (\sigma^{\rm sta}_{e^+})^2}$', color='tab:gray', fontsize=22)
    ax.fill_between([0, 20], 400, 600, color='tab:gray', alpha=0.2, zorder=1)

    E = np.logspace(1, 3, 1000)
    E3 = np.power(E, 3.3)
    
    I0, alpha = 21.801, 3.280
    ax.plot(E, E3 * SPL(E, [I0, alpha]), zorder=9, color='r', ls='--', label=r'$\chi^2$/dof = 96 / 35')

    I0, alpha, Eb, dalpha, s = 22.02, 3.321, 37.2, 0.08, 0.007
    ax.plot(E, E3 * BPL(E, [I0, alpha, Eb, dalpha, s]), zorder=9, color='r', label='$\chi^2$/dof = 25 / 32')

    plot_data(ax, 'data/AMS-02_e-_minus_e+_energy.txt', 3.3, 1., 'o', 'tab:blue', r'e$^-$ - e$^+$', zorder=3)
    ax.text(30., 445., r'$e^-$ - $e^+$', color='tab:blue', fontsize=23)

    ax.legend(fontsize=22, loc='best')
    savefig(plt, 'AMS02-electrons-minus-positrons.pdf')
 
def plot_positrons_background():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([2e1, 7e2])
        ax.set_ylabel(r'E$^{3}$ I [GeV$^{2}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_ylim([11, 23])

    fig = plt.figure(figsize=(10.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    plot_data(ax, 'data/AMS-02_e+_energy.txt', 3.0, 1., 'o', 'tab:blue', r'AMS-02 e$^+$', 1)

    E = np.logspace(1, 3, 1000)
    E3 = np.power(E, 3.0)

    b = BACKGROUND(E, [5.799e-01, 3.892e+00, 5.381e+01, 8.515e-02, 2.497e+00, 6.415e+02])
 
    ax.plot(E, E3 * b, zorder=9, color='r', ls='-', label=r'$\chi^2$/dof = 38 / 31')

    ax.fill_between([0, 20], 10, 30, color='tab:gray', alpha=0.2, zorder=1)
    
    ax.legend(fontsize=22, loc='best')
    savefig(plt, 'AMS02-positrons-background-fit.pdf')

def plot_electrons_background():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([2e1, 7e2])
        ax.set_ylabel(r'E$^{3.3}$ I [GeV$^{2.3}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_ylim([400, 680])

    fig = plt.figure(figsize=(10.5, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    plot_data(ax, 'data/AMS-02_e-_energy.txt', 3.3, 1., 'o', 'tab:blue', r'AMS-02 e$^-$', 1)

    E = np.logspace(1, 3, 1000)
    E3 = np.power(E, 3.3)

    y = SPL(E, [2.180e+01, 3.280e+00])
    b = BACKGROUND(E, [5.969e-01, 4.093e+00, 5.385e+01, 8.649e-02, 2.480e+00, 6.676e+02])

    ax.plot(E, E3 * y, zorder=9, color='g', ls=':', label=r'SPL')
    ax.plot(E, E3 * (b + y), zorder=9, color='g', ls='-', label=r'$\chi^2$/dof = 136 / 67')

    y = BPL(E, [2.202e+01, 3.322e+00,  3.758e+01,  7.794e-02,  9.228e-03])
    b = BACKGROUND(E, [6.027e-01, 3.835e+00, 5.387e+01, 8.356e-02, 2.488e+00, 6.411e+02])
    
    ax.plot(E, E3 * y, zorder=9, color='r', ls=':', label=r'SBPL')
    ax.plot(E, E3 * (b + y), zorder=9, color='r', ls='-', label=r'$\chi^2$/dof = 66 / 64')

    ax.fill_between([0, 20], 100, 900, color='tab:gray', alpha=0.2, zorder=1)
    
    ax.legend(fontsize=22, loc='best')
    savefig(plt, 'AMS02-electrons-background-fit.pdf')

if __name__== "__main__":
    #plot_electrons_and_positrons()
    #plot_electrons_minus_positrons()
    #plot_positrons_background()
    plot_electrons_background()