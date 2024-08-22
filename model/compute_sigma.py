import scipy.stats

def compute_chi2(dof, sigmas):
    ci = scipy.stats.chi2.cdf(sigmas * sigmas, 1)
    chi_squared = scipy.stats.chi2.ppf(ci, dof)
    
    print (f'{dof} {sigmas} {ci} {chi_squared}')

if __name__== "__main__":
    compute_chi2(2., 4.7) # test DAMPE paper