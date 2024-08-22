import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('gryphon.mplstyle')
import numpy as np

def savefig(plt, plotname):
    print (plotname)
    plt.savefig(plotname)
        
def plot_data(ax, filename, slope, norm, fmt, color, label, zorder=3):
    def quadrature(a, b):
        return np.sqrt(a * a + b * b)
    
    E, y, err_stat_lo, err_stat_up, err_sys_lo, err_sys_up = np.loadtxt(filename,usecols=(0,1,2,3,4,5),unpack=True)
    y = norm * np.power(E, slope) * y
    y_err_lo = norm * np.power(E, slope) * quadrature(err_stat_lo, err_sys_lo)
    y_err_up = norm * np.power(E, slope) * quadrature(err_stat_up, err_sys_up)
    ax.errorbar(E, y, yerr=[y_err_lo, y_err_up], fmt=fmt, markeredgecolor=color, color=color,
                capsize=3.5, markersize=6, elinewidth=1.8, capthick=1.8, zorder=zorder)

def plot_protons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1, 3e3])
        ax.set_ylabel(r'E$^{2.7}$ I [GeV$^{1.7}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_yscale('log')
        ax.set_ylim([0.1, 200])

    fig = plt.figure(figsize=(10.5, 8.0))
    ax = fig.add_subplot(111)
    set_axes(ax)

    # Is rigidity!
    plot_data(ax, 'data/AMS-02_H_energy.txt', 2.7, 1e-2, 'o', 'tab:blue', r'p', 1)
    ax.text(0.8e3, 0.5e2, '$p$/100', color='tab:blue', fontsize=29)
    plot_data(ax, 'data/AMS-02_e-_energy.txt', 2.7, 1., 'o', 'r', r'e$^-$', 3)  
    ax.text(0.80e3, 11, r'$e^-$', color='r', fontsize=29)
    plot_data(ax, 'data/AMS-02_e+_energy.txt', 2.7, 1., 'o', 'g', r'e$^-$', 3)
    ax.text(0.80e3, 2.5, r'$e^+$', color='g', fontsize=29)
    plot_data(ax, 'data/AMS-02_ap_energy.txt', 2.7, 1., 'o', 'tab:gray', r'e$^+$', 2)
    ax.text(0.55e3, 0.95, r'$\bar p$', color='tab:gray', fontsize=29)

    E, I = np.loadtxt('FERMI_gammas_inner.txt', usecols=(0,1), unpack=True)
    ax.plot(E, np.power(E, 0.7) * I, color='tab:brown')
    ax.text(30, 0.35, r'$\gamma$ (inner Galaxy)', color='tab:brown', fontsize=22)

    #ax.legend(fontsize=25, loc='best')
    savefig(plt, 'TeVPA24-protons.pdf')

def plot_electrons_over_protons():
    def spl(E, params):
        I0, E0, alpha = params
        y = I0 * np.power(E / E0, alpha)
        return y

    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1e1, 1e3])
        ax.set_ylabel(r'H/$e^-$ $\times 10^{-2}$')
        #ax.set_yscale('log')
        ax.set_ylim([0.5, 9])

    fig = plt.figure(figsize=(10.5, 8.0))
    ax = fig.add_subplot(111)
    set_axes(ax)

    plot_data(ax, 'data/AMS-02_H_over_e-_energy.txt', 0., 0.01, 'o', 'tab:blue', r'p', 1)

    I0, E0, alpha = 3.45, 1e2, 0.360
    E = np.logspace(1, 3, 1000)
    ax.plot(E, spl(E, [I0, E0, alpha]), color='r')

    ax.fill_between([10., 50.], 0, 10, color='tab:gray', alpha=0.2)
    ax.fill_between([500., 1e4], 0, 10, color='tab:gray', alpha=0.2)
    ax.text(60., 4.5, r'$\alpha = 0.36$')
    ax.legend(fontsize=25, loc='best')
    savefig(plt, 'TeVPA24-electrons-over-protons.pdf')

if __name__== "__main__":
    #plot_protons()
    plot_electrons_over_protons()