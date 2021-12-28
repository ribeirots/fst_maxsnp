# File to store python3 functions I use on my popgen analyses (more especifically, for now, all the functions I used on the SNP FST Power Analysis project)
# tribeiro@wisc.edu

# Function to calculate p-value
def pvalue_calc(obs_value, null_distr):
    pval = len([i for i in null_distr if i >= obs_value]) / len(null_distr)
    return pval

## need to add the fst function
