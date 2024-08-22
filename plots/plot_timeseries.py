import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('gryphon.mplstyle')
import numpy as np

def savefig(plt, plotname):
    print (plotname)
    plt.savefig(plotname)
    
def plot_data(ax, filename, fmt, color, label, zorder=3):
    E, y, err_stat_lo, err_stat_up, err_sys_lo, err_sys_up = np.loadtxt(filename,usecols=(0,1,2,3,4,5),unpack=True)
    y_err_lo = (err_sys_lo + err_stat_lo) / y * 100.
    y_err_up = (err_sys_up + err_stat_up) / y * 100.

    ax.fill_between(E, -y_err_lo, y_err_up, color=color, alpha=0.3, label=label, zorder=zorder)

def get_data(filename):
    filename = 'timeseries/' + filename
    E_lo, E_up, y, err_stat_lo, err_stat_up, err_sys_lo, err_sys_up = np.loadtxt(filename,usecols=(2,3,4,5,6,7,8),unpack=True)
    ind = np.argsort(E_lo)
    return E_lo[ind], E_up[ind], y[ind], 0.5 * (err_stat_up[ind] + err_stat_lo[ind])
    
def build_timeseries(size = 3300):
    E_lo, E_hi, y, erry = get_data('data_exp1.dat')
    y_series = np.zeros((10, size))
    for i in range(1, size + 1):
        E_lo, E_hi, y, erry = get_data('data_exp' + str(i) + '.dat')
        assert(len(y) == 10)
        y_series[:,i - 1] = y
    return E_lo, E_hi, y_series

def plot_timeaverage():
    from scipy import stats
    
    def set_axes(ax):
        ax.set_xlabel('E [GeV]')
        ax.set_xscale('log')
        ax.set_xlim([1.6e0, 50.])
        ax.set_ylabel(r'median absolute deviation [\%]')
        ax.set_ylim([-32, 32])

    fig = plt.figure(figsize=(11.0, 8.0))
    ax = fig.add_subplot(111)
    set_axes(ax)
    
    plot_data(ax, 'data/AMS-02_e-_rigidity.txt','o', 'tab:blue', r'$\sigma_{\rm sys} + \sigma_{\rm stat}$', 3)
    
    E_lo, E_hi, y_series = build_timeseries()

    y_median, y_mad = np.array([]), np.array([])
    E = np.sqrt(E_lo * E_hi)
    for i in range(10):
        y_median = np.append(y_median, np.median(y_series[i,:]))
        y_mad = np.append(y_mad, stats.median_abs_deviation(y_series[i,:]))
        
    ax.errorbar(E[1:10], 0. * E[1:10], yerr=y_mad[1:10] / y_median[1:10] * 100., fmt='o', markeredgecolor='tab:red', color='tab:red',
                capsize=4.0, markersize=7, elinewidth=2.0, capthick=2.0, label='2011/05/20 - 2021/11/02', zorder=10)

    ax.vlines(20., -30., 30., color='tab:gray', ls=':')
    ax.legend(fontsize=19, loc='upper left')
    savefig(plt, 'AMS02-timeaverage.pdf')

if __name__== "__main__":
    plot_timeaverage()