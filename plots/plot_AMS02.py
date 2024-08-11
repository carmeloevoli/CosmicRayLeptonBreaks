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
    ax.fill_between(E, y - y_err_lo, y + y_err_up, color='tab:purple', alpha=0.25, zorder=1)
#    ax.errorbar(E, y, yerr=[y_err_lo, y_err_up], fmt=fmt, markeredgecolor=color, color=color,
#                capsize=3.5, markersize=6, elinewidth=1.8, capthick=1.8, zorder=zorder)

def plot_leptons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1e1, 1e3])
        ax.set_ylabel(r'E$^{3}$ I [GeV$^{2}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        #ax.set_yscale('log')
        ax.set_ylim([10, 230])
        #ax.ticklabel_format(axis='y', style='sci', scilimits=(3,3))

    fig = plt.figure(figsize=(11.0, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)

    plot_data(ax, 'data/AMS-02_e-_rigidity.txt', 3.0, 1., 'o', 'tab:red', r'e$^-$', 1)
    plot_data(ax, 'data/AMS-02_e+_rigidity.txt', 3.0, 4., 'o', 'tab:blue', r'e$^+$ [4x]', 2)

    ax.fill_between([10, 20], 10, 230, color='tab:gray', alpha=0.2)
    ax.legend(fontsize=25, loc='best')
    savefig(plt, 'AMS02-leptons.pdf')

def plot_electrons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1e1, 1e3])
        ax.set_ylabel(r'E$^{3.3}$ I [GeV$^{2.3}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        #ax.set_yscale('log')
        ax.set_ylim([390, 570])
        #ax.ticklabel_format(axis='y', style='sci', scilimits=(3,3))

    fig = plt.figure(figsize=(11.0, 8.5))
    ax = fig.add_subplot(111)
    set_axes(ax)


    E = np.logspace(1, 3, 1000)
    E3 = np.power(E, 3.3)

    I0, E0, alpha = 21.790, 20., 3.2788
    ax.plot(E, E3 * spl(E, [I0, E0, alpha]), zorder=9, color='r', ls='--', label=r'$\chi^2$/dof = 110 / 35')

    I0, E0, alpha, Eb, dalpha, s = 22.02, 20.0, 3.322, 37.2, -0.076, 11
    ax.plot(E, E3 * bpl(E, [I0, E0, alpha, Eb, dalpha, s]), zorder=9, color='r', label='$\chi^2$/dof = 29 / 32')

    ax.fill_between([10, 20], 100, 900, color='tab:gray', alpha=0.2)

    #plot_data(ax, 'data/AMS-02_e-_rigidity.txt', 3.3, 1., 'o', 'tab:gray', r'e$^-$', 1)
    plot_data(ax, 'data/AMS-02_e-_minus_e+_rigidity.txt', 3.3, 1., 'o', 'tab:blue', r'e$^-$ - e$^+$')

    ax.legend(fontsize=22, loc='best')

    savefig(plt, 'AMS02-electrons.pdf')

if __name__== "__main__":
    #plot_leptons()
    plot_electrons()

