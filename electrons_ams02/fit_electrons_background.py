from iminuit import Minuit
from utils import load_data, BACKGROUND, SPL, BPL, chi2_single, compute_p_value, compute_sigmas

def fit_background(params: tuple) -> tuple:
    # Load the data
    filename = 'data/AMS-02_e+_energy.txt'
    energy, flux, err_lo, err_up = load_data(filename, add_stat_u=False)

    # Define the chi-squared function
    def chi2_function(C1, alpha1, Eb1, C2, alpha2, Ec):
        chi2 = 0.0
        for e, f, e_lo, e_up in zip(energy, flux, err_lo, err_up):
            model_flux = BACKGROUND(e, [C1, alpha1, Eb1, C2, alpha2, Ec])
            sigma = e_up if model_flux > f else e_lo
            chi2 += chi2_single(model_flux, f, sigma)
        return chi2

    # Unpack initial parameters
    C1, alpha1, Eb1, C2, alpha2, Ec = params

    # Initialize Minuit
    m = Minuit(chi2_function, C1=C1, alpha1=alpha1, Eb1=Eb1, C2=C2, alpha2=alpha2, Ec=Ec)
    m.errordef = Minuit.LEAST_SQUARES

    # Set parameter limits
    m.limits['C1'] = (0.1, 3.)   
    m.limits['alpha1'] = (1., 5.)   
    m.limits['Eb1'] = (10., 100.)   
    m.limits['C2'] = (0.01, 1.)   
    m.limits['alpha2'] = (1., 5.)   
    m.limits['Ec'] = (100., 2000.)   

    # Perform the fit
    m.simplex()
    m.migrad()
    m.hesse()

    # Degrees of freedom
    dof = len(energy) - m.nfit - 1

    return m.values, m.errors, m.fval, dof

def fit_spl(params: tuple) -> tuple:
    # Load positron and electron data
    pos_energy, pos_flux, pos_err_lo, pos_err_up = load_data('data/AMS-02_e+_energy.txt', add_stat_u=False)
    ele_energy, ele_flux, ele_err_lo, ele_err_up = load_data('data/AMS-02_e-_energy.txt', add_stat_u=False)

    # Define chi-squared function
    def chi2_function(C1, alpha1, Eb1, C2, alpha2, Ec, I0, alpha):
        chi2 = 0.0
        # Contribution from positron data (background only)
        for e, f, e_lo, e_up in zip(pos_energy, pos_flux, pos_err_lo, pos_err_up):
            model_flux = BACKGROUND(e, [C1, alpha1, Eb1, C2, alpha2, Ec])
            sigma = e_up if model_flux > f else e_lo
            chi2 += chi2_single(model_flux, f, sigma)
        # Contribution from electron data (SPL + background)
        for e, f, e_lo, e_up in zip(ele_energy, ele_flux, ele_err_lo, ele_err_up):
            model_flux = SPL(e, [I0, alpha]) + BACKGROUND(e, [C1, alpha1, Eb1, C2, alpha2, Ec])
            sigma = e_up if model_flux > f else e_lo
            chi2 += chi2_single(model_flux, f, sigma)
        return chi2

    # Unpack initial parameters
    C1, alpha1, Eb1, C2, alpha2, Ec, I0, alpha = params

    # Initialize Minuit
    m = Minuit(
        chi2_function, 
        C1=C1, alpha1=alpha1, Eb1=Eb1, C2=C2, alpha2=alpha2, Ec=Ec, I0=I0, alpha=alpha
    )
    m.errordef = Minuit.LEAST_SQUARES

    # Set parameter limits
    m.limits['C1'] = (0.1, 3.)   
    m.limits['alpha1'] = (1., 5.)   
    m.limits['Eb1'] = (10., 100.)   
    m.limits['C2'] = (0.01, 1.)   
    m.limits['alpha2'] = (1., 5.)   
    m.limits['Ec'] = (100., 2000.)   

    # Perform the fit
    m.simplex()
    m.migrad()
    m.hesse()

    # Degrees of freedom
    dof = len(pos_energy) + len(ele_energy) - m.nfit - 1

    return m.values, m.errors, m.fval, dof

def fit_bpl(params: tuple) -> tuple:
    # Load positron and electron data
    pos_energy, pos_flux, pos_err_lo, pos_err_up = load_data('data/AMS-02_e+_energy.txt', add_stat_u=True)
    ele_energy, ele_flux, ele_err_lo, ele_err_up = load_data('data/AMS-02_e-_energy.txt', add_stat_u=True)

    # Define chi-squared function
    def chi2_function(C1, alpha1, Eb1, C2, alpha2, Ec, I0, alpha, Eb, dalpha, s):
        chi2 = 0.0
        # Contribution from positron data (background only)
        for e, f, e_lo, e_up in zip(pos_energy, pos_flux, pos_err_lo, pos_err_up):
            model_flux = BACKGROUND(e, [C1, alpha1, Eb1, C2, alpha2, Ec])
            sigma = e_up if model_flux > f else e_lo
            chi2 += chi2_single(model_flux, f, sigma)
        # Contribution from electron data (SPL + background)
        for e, f, e_lo, e_up in zip(ele_energy, ele_flux, ele_err_lo, ele_err_up):
            model_flux = BPL(e, [I0, alpha, Eb, dalpha, s]) + BACKGROUND(e, [C1, alpha1, Eb1, C2, alpha2, Ec])
            sigma = e_up if model_flux > f else e_lo
            chi2 += chi2_single(model_flux, f, sigma)
        return chi2

    # Unpack initial parameters
    C1, alpha1, Eb1, C2, alpha2, Ec, I0, alpha, Eb, dalpha, s = params

    # Initialize Minuit
    m = Minuit(
        chi2_function, 
        C1=C1, alpha1=alpha1, Eb1=Eb1, C2=C2, alpha2=alpha2, Ec=Ec, I0=I0, alpha=alpha, Eb=Eb, dalpha=dalpha, s=s,
    )
    m.errordef = Minuit.LEAST_SQUARES

    # Set parameter limits
    m.limits['C1'] = (0.1, 3.)   
    m.limits['alpha1'] = (1., 5.)   
    m.limits['Eb1'] = (10., 100.)   
    m.limits['C2'] = (0.01, 1.)   
    m.limits['alpha2'] = (1., 5.)   
    m.limits['Ec'] = (100., 2000.)   

    # Perform the fit
    m.simplex()
    m.migrad()
    m.hesse()

    # Degrees of freedom
    dof = len(pos_energy) + len(ele_energy) - m.nfit - 1

    return m.values, m.errors, m.fval, dof

if __name__== "__main__":
    # Initial parameters for BACKGROUND
    params_initial = [0.58, 3.89, 53.8, 0.085, 2.50, 0.64e3]
    bkg_values, bkg_errors, bkg_chi2, bkg_dof = fit_background(params_initial)
    print("Background Model Fit Results:")
    print(f"C1: {bkg_values[0]:4.1f} {bkg_errors[0]:4.1f}")
    print(f"alpha1: {bkg_values[1]:4.1f} {bkg_errors[1]:4.1f}")
    print(f"Eb1: {bkg_values[2]:4.1f} {bkg_errors[2]:4.1f}")
    print(f"C2: {bkg_values[3]:4.2f} {bkg_errors[3]:4.2f}")
    print(f"alpha2: {bkg_values[4]:4.1f} {bkg_errors[4]:4.1f}")
    print(f"Ec: {bkg_values[5]:4.0f} {bkg_errors[5]:4.0f}")
    print(f'Chi2/dof: {bkg_chi2:5.0f} / {bkg_dof}')
    print(f"p-value: {compute_p_value(bkg_chi2, bkg_dof):.3f}")
    print('')
    # Initial parameters for SPL
    spl_initial = [0.58, 3.89, 53.8, 0.085, 2.50, 0.64e3, 21.81, 3.281]
    spl_values, spl_errors, spl_chi2, spl_dof = fit_spl(spl_initial)
    print('SPL + Background Model Fit Results:')
    print(f'I0: {spl_values[6]:5.2f} {spl_errors[6]:5.2f}')
    print(f'alpha: {spl_values[7]:5.2f} {spl_errors[7]:5.2f}')
    print(f'Chi2 / dof : {spl_chi2:5.0f} / {spl_dof}')
    print('')
    bpl_initial = [0.58, 3.89, 53.8, 0.085, 2.50, 0.64e3, 22.02, 3.321, 45, 0.1, 0.01]
    bpl_values, bpl_errors, bpl_chi2, bpl_dof = fit_bpl(bpl_initial)
    print('BPL + Background Model Fit Results:')
    print(f'Chi2 / dof : {bpl_chi2:5.0f} / {bpl_dof}')
    chi_squared = spl_chi2 - bpl_chi2
    dof = spl_dof - bpl_dof
    print(f'{compute_sigmas(chi_squared, dof):5.1f}')