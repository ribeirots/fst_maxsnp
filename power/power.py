# power analyses SNP/Window FST
# tribeiro@wisc.edu

import os, sys, re
import numpy as np
# [1]: input name, file with the filename patterns for adpt and neutral simulations. Without file extension, which should be .txt. 
# [2] path to adpt files, [3] path to neutral files, [4]: name of the output file.

output = open(sys.argv[4],'w')
header = 'ModelID\tSNPFST_Power\tWINFST_Power\tSNPFST_AvgPvalue\tWINFST_AvgPvalue\tAdptRepl\tSimRepl\tLastSimRepl\n'
output.write(header)

repl_number = 50
pvalue_thresh = 0.05
window_size_power = 0 # Change to 1 for window size

path_adpt = sys.argv[2]
path_neut = sys.argv[3]

previous_neutral_pattern = 'no-pattern-yet' 

with open(sys.argv[1]) as filename_patterns: # First column: pattern from files that you want to calculate power for. 2nd column: pattern from files with the neutral distribution
    for pattern in filename_patterns:
        patterns = re.split(',',pattern[:-1])
        if window_size_power == 1:
            wlength = re.split('wlength',patterns[1])
            if 'LN' in wlength[0]:
                base_wlength = 100000
            else:
                base_wlength = 5000
            current_wlength = int(wlength[1])
            pvalue_thresh = 0.05 / (base_wlength /current_wlength)
#        print(patterns[0])
#        print((base_wlength/current_wlength))
#        print(pvalue_thresh)
        if patterns[1] != previous_neutral_pattern:
            previous_neutral_pattern = patterns[1]
            nsnp_fst = [] # neutral snp fst
            nwin_fst = []
            sim_repl = 0
            for i in range(0,repl_number): # loops through possible replicate ID
                if os.path.isfile(path_neut+patterns[1]+'_rep'+str(i)+'.txt'): # check if the file exists, if it does open it and get all the SNP and Window FST
                    with open(path_neut+patterns[1]+'_rep'+str(i)+'.txt') as neutral:
                        next(neutral) # only keep this if the input file with FST measures has a header.
                        for l in neutral:
                            if sim_repl < 10000:
                                nsim = list(map(float, l.split()))
                                nsnp_fst.append(nsim[0]) 
                                nwin_fst.append(nsim[1])
                                sim_repl += 1
        
        psnp = [] # list to hold the pvalues for snp and window fst
        pwin = []
        adpt_repl = 0
        j = -1
        if len(nsnp_fst) > 0:
            for j in range(0,repl_number): # loops through possible replicate ID - this time for adpt
                if os.path.isfile(path_adpt+patterns[0]+'_rep'+str(j)+'.txt'):
                    with open(path_adpt+patterns[0]+'_rep'+str(j)+'.txt') as adpt:
                        next(adpt)
                        for l in adpt:
                            if adpt_repl < 10000:
                                adptsim = list(map(float, l.split()))
                                snp_fst = adptsim[0] 
                                win_fst = adptsim[1]
                                psnp.append(len([i for i in nsnp_fst if i >= snp_fst])/len(nsnp_fst))
                                pwin.append(len([i for i in nwin_fst if i >= win_fst])/len(nwin_fst))
                                adpt_repl += 1
                            else:
                                break
                if adpt_repl >= 10000:
                    break
        if len(psnp) > 0:
            meanPValueSNP = sum(psnp)/len(psnp)
            meanPValueWindow = sum(pwin)/len(pwin)
            power_snp = len([i for i in psnp if i <= pvalue_thresh])/len(psnp)
            power_win = len([i for i in pwin if i <= pvalue_thresh])/len(pwin)
        else:
            power_snp = 'NA'
            power_win = 'NA'
            meanPValueSNP = 'NA'
            meanPValueWindow = 'NA'
        outting = [patterns[0], str(power_snp), str(power_win), str(meanPValueSNP), str(meanPValueWindow), str(adpt_repl), str(sim_repl), str(j)]
        output.write('\t'.join(outting)+'\n')

output.close()

