from iminuit import Minuit
from utils import load_data, SPL, BPL, chi2_single, compute_p_value, compute_sigmas

#FILENAME = 'data/AMS-02_e-_minus_e+_energy.txt'
#FILENAME = 'data/AMS-02_e-_minus_e+_statUp_energy.txt'
FILENAME = 'data/AMS-02_e-_minus_e+_statAdd_energy.txt'

# Function to fit data using a Single Power Law (SPL) model
def fit_spl(initial_params):
    # Load filtered data
    energy, flux, err_lo, err_up = load_data(FILENAME)

    # Define chi-squared function for SPL
    def chi2_function(I0, alpha):
        chi2 = 0.
        for e, f, e_lo, e_up in zip(energy, flux, err_lo, err_up):
            model_flux = SPL(e, [I0, alpha])
            sigma = e_up if model_flux > f else e_lo
            chi2 += chi2_single(model_flux, f, sigma)
        return chi2

    # Initialize Minuit with initial parameters
    I0, alpha = initial_params
    m = Minuit(chi2_function, I0=I0, alpha=alpha)
    m.errordef = Minuit.LEAST_SQUARES

    # Set parameter limits
    m.limits['alpha'] = (3.1, 3.4)
    m.limits['I0'] = (10., 50.)

    # Perform minimization
    m.simplex()
    m.migrad()
    m.hesse()

    dof = len(energy) - m.nfit - 1

    return m.values, m.errors, m.fval, dof

# Function to fit data using a Broken Power Law (BPL) model
def fit_bpl(initial_params):
    # Load filtered data
    energy, flux, err_lo, err_up = load_data(FILENAME)

    # Define chi-squared function for BPL
    def chi2_function(I0, alpha, Eb, dalpha, s):
        chi2 = 0.
        for e, f, e_lo, e_up in zip(energy, flux, err_lo, err_up):
            model_flux = BPL(e, [I0, alpha, Eb, dalpha, s])
            sigma = e_up if model_flux > f else e_lo
            chi2 += chi2_single(model_flux, f, sigma)
        return chi2

    # Initialize Minuit with initial parameters
    I0, alpha, Eb, dalpha, s = initial_params
    m = Minuit(chi2_function, I0=I0, alpha=alpha, Eb=Eb, dalpha=dalpha, s=s)
    m.errordef = Minuit.LEAST_SQUARES

    # Set parameter limits
    m.limits['alpha'] = (3.1, 3.6)
    m.limits['I0'] = (10., 50.)
    m.limits['Eb'] = (10., 100.)
    m.limits['dalpha'] = (0., 0.3)
    m.limits['s'] = (0.001, 0.1)

    # Perform minimization
    m.simplex()
    m.migrad()
    m.hesse()

    dof = len(energy) - m.nfit - 1

    return m.values, m.errors, m.fval, dof 

# Main execution
if __name__ == "__main__":
    # Print model name
    print(FILENAME)
    print('')
    # Initial parameters for SPL
    spl_initial = [21.81, 3.281]
    spl_values, spl_errors, spl_chi2, spl_dof = fit_spl(spl_initial)
    # Print results
    print('Single Power Law:')
    print(f'I0: {spl_values[0]:5.2f} {spl_errors[0]:5.2f}')
    print(f'alpha: {spl_values[1]:5.3f} {spl_errors[1]:5.3f}')
    print(f'Chi2/dof: {spl_chi2:4.0f} / {spl_dof}')
    print(f'p-value : {compute_p_value(spl_chi2, spl_dof):5.2e}')
    print('')
    # Initial parameters for BPL
    bpl_initial = [22.02, 3.321, 45, 0.1, 0.01]
    bpl_values, bpl_errors, bpl_chi2, bpl_dof = fit_bpl(bpl_initial)
    # Print results
    print('Broken Power Law:')
    print(f'I0: {bpl_values[0]:5.2f} {bpl_errors[0]:5.2f}')
    print(f'alpha: {bpl_values[1]:5.3f} {bpl_errors[1]:5.3f}')
    print(f'Eb: {bpl_values[2]:5.1f} {bpl_errors[2]:5.1f} ')
    print(f'dalpha: {bpl_values[3]:5.2f} {bpl_errors[3]:5.2f} ')
    print(f's: {bpl_values[4]:5.3f} {bpl_errors[4]:5.3f} ')
    print(f'Chi2 / dof : {bpl_chi2:5.2f} / {bpl_dof}')
    print(f'p-value : {compute_p_value(bpl_chi2, bpl_dof):5.2f}')
    print('')
    print('Difference between models')
    chi_squared = spl_chi2 - bpl_chi2
    dof = spl_dof - bpl_dof
    print(f'{compute_sigmas(chi_squared, dof):5.1f}')