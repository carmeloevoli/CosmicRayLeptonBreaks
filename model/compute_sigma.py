import scipy.stats

df = 3.
    
s = 3.

ci = scipy.stats.chi2.cdf(s * s, 1)
    
chi_squared = scipy.stats.chi2.ppf(ci, df)
    
print (df, s, ci, chi_squared)
