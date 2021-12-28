# File to store python3 functions for the SNP FST Power Analysis project
# tribeiro@wisc.edu

##############################################################################################
##############################################################################################
############################### Hard coded functions #########################################
######################### (specific to their own project) ####################################
##############################################################################################
##############################################################################################

import os.path, sys
sys.path.insert(1, '/home/tiago/snp_fst/fst_distribution/')
from ribeiro_popgen import *


# Function to create dictionary for a given recombination list of files
def sim_recomb_list(model,recrate_dir, recrate_file,chrm, p1_code, p2_code, replicates):
    snp_dict = {}
    win_dict = {}
    poppair = [p1_code, p2_code]
    key_dict = str(p1_code)+'_'+str(p2_code)
    first_reading = 1
    for i in range(0,replicates):
        filepath = '../../gene_conv/fst/'+model+'/'+recrate_dir+'/fst_recrate'+recrate_file+'_'+str(i)+'_chr'+chrm+'.csv'
        if os.path.isfile(filepath):
            with open(filepath) as file_sims:
                for line in file_sims:
                    if line[0] != 'n' and line[0] != 'M':
                        line = line.split()
                        if 'winFST' in line:
                            print(line)
                        elif int(line[5]) in poppair and int(line[6]) in poppair:
                            if first_reading == 1:
                                first_reading = 0
                                snp_dict[key_dict] = []
                                win_dict[key_dict] = []
                            snp_dict[key_dict].append(float(line[0]))
                            win_dict[key_dict].append(float(line[1]))

    dict_lists = []
    dict_lists.append(snp_dict)
    dict_lists.append(win_dict)
    return dict_lists

# Function to calculate de pvalues for a given Sim x Emp comparison
def pvalue_recrates(win_full,empirical_summ, sim, emp, p1, p2, outputfilename, rep_sim_number, emp_fst_key=1, null_fst_key=1): # win_full: path+file to window_chr basic summary.  empirical_summ: to empirical fst per window. sim = model (eg. NA_3L). emp = empirical data (eg. 2R). rep_sim_number = how many replicate files containing the simulations exist. emp_fst_key: turn to 0 if you dont want to create a file with empirical fst and pvalues. null_fst_key: turn to 0 if you dont want to create a file with null distribution fst and pvalues.

    recombination_rate_dir = ['00','05','1','15','2','3','inf']
    recombination_rate_fil = ['0', '0.5', '1', '1.5', '2', '3', 'inf']

    sim_fst_AoA = []
    for i in range(0,len(recombination_rate_dir)):
        dicts_rec = sim_recomb_list(sim, recombination_rate_dir[i], recombination_rate_fil[i], emp, p1, p2, rep_sim_number)
        sim_fst_AoA.append(dicts_rec)

    
    pop_key = str(p1)+'_'+str(p2)
    
    # block to calculate pvalue for empirical data entry 
    if emp_fst_key == 1:
        recs = []
        with open(win_full) as recr: # '../../recwind../recw....txt'
            for line in recr:
                line = line.split()
                recs.append(line)

        win_summ = []
        rec_index = 0
        # window FST = 2, SNP = 3
        with open(empirical_summ) as win_summary: # '../../emp_minalle/win..ChX.csv'
            for line in win_summary:
                if line[0] != 'S':
                    line = line.split()
                    while int(line[0]) > int(recs[rec_index][0]):
                        rec_index += 1
                    line.append(recs[rec_index][-1])
                    win_summ.append(line)

        recs = []

        # calculating p-value
        output = open(outputfilename,'w')
        header = 'Start\tEnd\tWindowFST\tMaxSNPFST\tMaxSNPPosition\tRecRate\tWindowPvalue\tSNPPvalue\n'
        output.write(header)
        for w in win_summ:
            if w[2] != 'NA':
                if float(w[5]) >= 0.5:
                    if float(w[5]) < 1:
                        w.append(str(pvalue_calc(float(w[2]), sim_fst_AoA[2][1][pop_key]))) # window pvalue
                        w.append(str(pvalue_calc(float(w[3]), sim_fst_AoA[2][0][pop_key]))) # max_snp pvalue
                    elif float(w[5]) < 1.5:
                        w.append(str(pvalue_calc(float(w[2]), sim_fst_AoA[3][1][pop_key]))) # window pvalue
                        w.append(str(pvalue_calc(float(w[3]), sim_fst_AoA[3][0][pop_key]))) # max_snp pvalue
                    elif float(w[5]) < 2:
                        w.append(str(pvalue_calc(float(w[2]), sim_fst_AoA[4][1][pop_key]))) # window pvalue
                        w.append(str(pvalue_calc(float(w[3]), sim_fst_AoA[4][0][pop_key]))) # max_snp pvalue
                    elif float(w[5]) < 3:
                        w.append(str(pvalue_calc(float(w[2]), sim_fst_AoA[5][1][pop_key]))) # window pvalue
                        w.append(str(pvalue_calc(float(w[3]), sim_fst_AoA[5][0][pop_key]))) # max_snp pvalue
                    else:
                        w.append(str(pvalue_calc(float(w[2]), sim_fst_AoA[6][1][pop_key]))) # window pvalue
                        w.append(str(pvalue_calc(float(w[3]), sim_fst_AoA[6][0][pop_key]))) # max_snp pvalue
                    output.write('\t'.join(w)+'\n')
        output.close()

    # block to calculate pvalue for the null, simulated distribution
    if null_fst_key == 1:
        # calculation null output
        output = open('sim_'+outputfilename,'w')
        header = 'WindowFST\tMaxSNPFST\tWindowPvalue\tSNPPvalue\tRecRate\n'
        output.write(header)
        for i in range(2, len(sim_fst_AoA)): # starting from recomb rate 1 (>= 0.5)
            for j in range(0, len(sim_fst_AoA[i][0][pop_key])):
                out_row = []
                out_row.append(str(sim_fst_AoA[i][1][pop_key][j])) # window fst
                out_row.append(str(sim_fst_AoA[i][0][pop_key][j])) # max_snp fst
                out_row.append(str(pvalue_calc(sim_fst_AoA[i][1][pop_key][j], sim_fst_AoA[i][1][pop_key]))) # window pvalue
                out_row.append(str(pvalue_calc(sim_fst_AoA[i][0][pop_key][j], sim_fst_AoA[i][0][pop_key]))) # max_snp pvalue
                out_row.append(recombination_rate_fil[i])
                output.write('\t'.join(out_row)+'\n')
        output.close()








