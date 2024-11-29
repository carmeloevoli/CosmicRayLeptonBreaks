from scipy.stats import chi2
from scipy.stats import norm

chi2_PL = 59  
chi2_SBPL = 16 
dof_PL = 35   
dof_SBPL = 32     

# ----------------------------------------------------------------------------------------------------
def compute_p_value():

    delta_chi2 = chi2_PL - chi2_SBPL
    delta_dof = dof_PL - dof_SBPL

    if delta_dof <= 0:
        raise ValueError('Degrees of freedom difference must be positive.')

    return chi2.sf(delta_chi2, delta_dof)

# ----------------------------------------------------------------------------------------------------
def p_value_to_sigma():

    p_value = compute_p_value()

    if p_value <= 0 or p_value > 1:
        raise ValueError('The p-value must be in the range (0, 1].')
        
    return norm.ppf(1 - p_value / 2) 

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    print(compute_p_value())
    print(p_value_to_sigma())

# ----------------------------------------------------------------------------------------------------