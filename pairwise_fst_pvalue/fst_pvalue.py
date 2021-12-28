# combine simulated fst and calculate pvalue for empirical fst
# tribeiro@wisc.edu
# arg 1: path to output

import os.path, sys
sys.path.insert(1, '/home/tiago/snp_fst/fst_distribution/')
from snpfst_funcs import *
from ribeiro_popgen import *

pop1 = 'RG'
pop2 = 'ZI'

modelX = 'NA_X_EG_FR'
model2R = 'NA_2R'
model3L = 'NA_3L'

pop1_x_code = 1
pop2_x_code = 0
pop1_2r_code = 1
pop2_2r_code = 0
pop1_3l_code = 1
pop2_3l_code = 0

n_sim = 501
print(pop1+pop2+'\n')

emp_chrm = 'X'
pvalue_recrates('../../recwindows/recwindows_'+emp_chrm+'_full.txt', '../../empirical_minAllele3/window_'+pop1+'_'+pop2+'_Chr'+emp_chrm+'_summary.csv', modelX, emp_chrm, pop1_x_code, pop2_x_code, 'emp_pvalue_'+pop1+'_'+pop2+'_'+emp_chrm+'_simX.txt', n_sim, emp_fst_key=0) # turn key to 1 to calculate empirical distribution

for ch in ['2L','2R','3L','3R']:
    emp_chrm = ch
    pvalue_recrates('../../recwindows/recwindows_'+emp_chrm+'_full.txt', '../../empirical_minAllele3/window_'+pop1+'_'+pop2+'_Chr'+emp_chrm+'_summary.csv', model2R, emp_chrm, pop1_2r_code, pop2_2r_code, 'emp_pvalue_'+pop1+'_'+pop2+'_'+emp_chrm+'_sim2R.txt', n_sim, emp_fst_key=0)

for ch in ['2L','2R','3L','3R']:
    emp_chrm = ch
    pvalue_recrates('../../recwindows/recwindows_'+emp_chrm+'_full.txt', '../../empirical_minAllele3/window_'+pop1+'_'+pop2+'_Chr'+emp_chrm+'_summary.csv', model3L, emp_chrm, pop1_3l_code, pop2_3l_code, 'emp_pvalue_'+pop1+'_'+pop2+'_'+emp_chrm+'_sim3L.txt', n_sim, emp_fst_key=0)
