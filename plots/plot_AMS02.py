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
    ax.errorbar(E, y, yerr=[y_err_lo, y_err_up], fmt=fmt, markeredgecolor=color, color=color,
                capsize=3.5, markersize=6, elinewidth=1.8, capthick=1.8, zorder=zorder)
    y_err_lo = norm * np.power(E, slope) * err_sys_lo
    y_err_up = norm * np.power(E, slope) * err_sys_up
#    ax.fill_between(E, y - y_err_lo, y + y_err_up, color='tab:purple', alpha=0.25, zorder=1)
#    ax.errorbar(E, y, yerr=[y_err_lo, y_err_up], fmt=fmt, markeredgecolor=color, color=color,
#                capsize=3.5, markersize=6, elinewidth=1.8, capthick=1.8, zorder=zorder)

def plot_electrons_and_positrons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1e1, 7e2])
        ax.set_ylabel(r'E$^{3}$ I [GeV$^{2}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_ylim([40, 230])

    fig = plt.figure(figsize=(11.0, 8.0))
    ax = fig.add_subplot(111)
    set_axes(ax)

    plot_data(ax, 'data/AMS-02_e-_rigidity.txt', 3.0, 1., 'o', 'tab:red', r'e$^-$', 1)
    ax.text(35., 170., 'electrons', color='tab:red', fontsize=21)

    plot_data(ax, 'data/AMS-02_e+_rigidity.txt', 3.0, 5., 'o', 'tab:blue', r'e$^+$ [5x]', 2)
    ax.text(115., 70., 'positrons [5x]', color='tab:blue', fontsize=21)

    ax.text(130., 210., 'AMS-02 (stat only)', color='tab:gray', fontsize=25)

    ax.fill_between([0, 20], 10, 230, color='tab:gray', alpha=0.2, zorder=1)
    
    #ax.legend(fontsize=25, loc='best')
    savefig(plt, 'AMS02-electrons-and-positrons.pdf')

def plot_electrons_minus_positrons():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1e1, 7e2])
        ax.set_ylabel(r'E$^{3.3}$ I [GeV$^{2.3}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_ylim([400, 540])

    fig = plt.figure(figsize=(11.0, 8.0))
    ax = fig.add_subplot(111)
    set_axes(ax)

    ax.text(100., 420., r'$\sigma^{\rm sta}_{e^- - e^+} = \sigma^{\rm sta}_{e^-} + \sigma^{\rm sta}_{e^+}$', color='tab:gray', fontsize=22)
    ax.fill_between([0, 20], 400, 600, color='tab:gray', alpha=0.2, zorder=1)

    E = np.logspace(1, 3, 1000)
    E3 = np.power(E, 3.3)
    
    I0, E0, alpha = 21.807, 20, 3.281
    ax.plot(E, E3 * spl(E, [I0, E0, alpha]), zorder=9, color='r', ls='--', label=r'$\chi^2$/dof = 59 / 35')

    I0, E0, alpha, Eb, dalpha, s = 22.022, 20, 3.321, 37.131, -0.075, 11.842
    ax.plot(E, E3 * bpl(E, [I0, E0, alpha, Eb, dalpha, s]), zorder=9, color='r', label='$\chi^2$/dof = 15 / 32')

    plot_data(ax, 'data/AMS-02_e-_minus_e+_rigidity.txt', 3.3, 1., 'o', 'tab:blue', r'e$^-$ - e$^+$', zorder=3)
    ax.text(30., 445., r'$e^-$ - $e^+$', color='tab:blue', fontsize=23)

    ax.legend(fontsize=22, loc='best')
    savefig(plt, 'AMS02-electrons-minus-positrons.pdf')
    
def plot_electrons_minus_positrons_up():
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1e1, 7e2])
        ax.set_ylabel(r'E$^{3.3}$ I [GeV$^{2.3}$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        ax.set_ylim([400, 540])

    fig = plt.figure(figsize=(11.0, 8.0))
    ax = fig.add_subplot(111)
    set_axes(ax)

    ax.text(100., 420., r'$\sigma^{\rm sta}_{e^- - e^+} = \sigma^{\rm sta}_{e^-} + \sigma^{\rm sta}_{e^+}$', color='tab:gray', fontsize=22)
    ax.fill_between([0, 20], 400, 600, color='tab:gray', alpha=0.2, zorder=1)

    E = np.logspace(1, 3, 1000)
    E3 = np.power(E, 3.3)
    
    I0, E0, alpha = 21.789, 20, 3.282
    ax.plot(E, E3 * spl(E, [I0, E0, alpha]), zorder=9, color='r', ls='--', label=r'$\chi^2$/dof = 57 / 35')

    I0, E0, alpha, Eb, dalpha, s = 22.001, 20, 3.321, 36.738, -0.072, 13.691
    ax.plot(E, E3 * bpl(E, [I0, E0, alpha, Eb, dalpha, s]), zorder=9, color='r', label='$\chi^2$/dof = 15 / 32')

    plot_data(ax, 'data/AMS-02_e-_minus_e+_statUp_rigidity.txt', 3.3, 1., 'o', 'tab:purple', r'e$^-$ - e$^+$', zorder=3)
    ax.text(23., 448., r'$e^-$ - ($e^+$ + $\sigma_+^{\rm sys}$)', color='tab:purple', fontsize=23)

    ax.legend(fontsize=22, loc='best')
    savefig(plt, 'AMS02-electrons-minus-positrons-up.pdf')

def plot_sigmas():
    def compute_chi2(dof, sigmas):
        import scipy.stats
        ci = scipy.stats.chi2.cdf(sigmas * sigmas, 1)
        return scipy.stats.chi2.ppf(ci, dof)

    def set_axes(ax):
        ax.set_xlabel(r'$\Delta \chi^2$')
#        ax.set_xscale('log')
        ax.set_xlim([0, 60])
        ax.set_ylabel(r'$\sigma$')
        ax.set_ylim([0, 8])
    
    fig = plt.figure(figsize=(11.0, 8.0))
    ax = fig.add_subplot(111)
    set_axes(ax)

    dof = 3.

    sigmas = np.linspace(0, 9, 100)
    
    y = compute_chi2(dof, sigmas)

    ax.plot(y, sigmas)
    
    items = [i for i in range(len(y)) if y[i] < 59 - 15]

    ax.fill_between(y[items], 0., sigmas[items], alpha=0.3)
    
    ax.text(5., 7., 'dof = 3')

    savefig(plt, 'AMS02-sigmas.pdf')

if __name__== "__main__":
    plot_electrons_and_positrons()
    plot_electrons_minus_positrons()
    plot_electrons_minus_positrons_up()
    plot_sigmas()
