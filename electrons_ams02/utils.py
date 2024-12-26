import numpy as np
from scipy.stats import chi2, norm

# Constants
ENORM = 20.0  # GeV

# Single Power-Law (SPL) Model
def SPL(E: np.ndarray, params: tuple) -> np.ndarray:
    """
    Compute the single power-law model values.

    Parameters:
        E (np.ndarray): Energy array.
        params (tuple): Parameters (I0, alpha).

    Returns:
        np.ndarray: Computed model values.
    """
    I0, alpha = params
    return (I0 / 1e3) * np.power(E / ENORM, -alpha)

# Broken Power-Law (BPL) Model
def BPL(E: np.ndarray, params: tuple) -> np.ndarray:
    """
    Compute the broken power-law model values.

    Parameters:
        E (np.ndarray): Energy array.
        params (tuple): Parameters (I0, alpha, Eb, dalpha, s).

    Returns:
        np.ndarray: Computed model values.
    """
    I0, alpha, Eb, dalpha, s = params
    y = (I0 / 1e3) * np.power(E / ENORM, -alpha)
    y *= np.power(1.0 + np.power(E / Eb, dalpha / s), s)
    return y

# Background Model
def BACKGROUND(E: np.ndarray, params: tuple) -> np.ndarray:
    """
    Compute the background model values.

    Parameters:
        E (np.ndarray): Energy array.
        params (tuple): Parameters (C1, alpha1, Eb1, C2, alpha2, Ec).

    Returns:
        np.ndarray: Computed background model values.
    """
    C1, alpha1, Eb1, C2, alpha2, Ec = params
    y = (C1 / 1e3) * np.power(E / ENORM, -alpha1)
    y += (C2 / 1e3) * np.power(E / Eb1, -alpha2) * np.exp(-E / Ec)
    return y

# Function to load data within a specified energy range
def load_data(
        filename : str,
        min_energy: float = 20., max_energy: float = 1e20, 
        add_stat_u: bool = False ) -> tuple:
    """
    Load and filter data within a specified energy range.

    Parameters:
        filename (string): filename to read
        min_energy (float): Minimum energy threshold.
        max_energy (float): Maximum energy threshold.
        add_stat_u (bool): Whether to add statistical uncertainties.

    Returns:
        tuple: Filtered energy, flux, lower errors, and upper errors.
    """
    E, y, err_stat_lo, err_stat_up, err_sys_lo, err_sys_up = np.loadtxt(
        filename, usecols=(0, 1, 2, 3, 4, 5), unpack=True
    )
    if add_stat_u:
        err_stat_lo = np.sqrt(err_stat_lo**2 + err_sys_lo**2)
        err_stat_up = np.sqrt(err_stat_up**2 + err_sys_up**2)

    # Boolean indexing for efficiency
    mask = (E > min_energy) & (E < max_energy)
    return E[mask], y[mask], err_stat_lo[mask], err_stat_up[mask]

# Function to calculate chi-squared for a single data point
def chi2_single(x: float, mu: float, sigma: float) -> float:
    """
    Compute chi-squared for a single data point.

    Parameters:
        x (float): Observed value.
        mu (float): Expected value.
        sigma (float): Standard deviation.

    Returns:
        float: Chi-squared value for the data point.
    """
    return ((x - mu) / sigma) ** 2

# Function to compute p-value
def compute_p_value(chi_squared: float, dof: int) -> float:
    """
    Compute the p-value for a given chi-squared value and degrees of freedom.

    Parameters:
        chi_squared (float): The chi-squared value.
        dof (int): Degrees of freedom.

    Returns:
        float: The p-value.
    """
    return chi2.sf(chi_squared, dof)

# Function to compute number of sigmas
def compute_sigmas(chi_squared: float, dof: int) -> float:
    """
    Compute the number of sigmas corresponding to a chi-squared value and degrees of freedom.

    Parameters:
        chi_squared (float): The chi-squared value.
        dof (int): Degrees of freedom.

    Returns:
        float: The number of sigmas.
    """
    p_value = compute_p_value(chi_squared, dof)
    return norm.isf(p_value)  # Inverse survival function