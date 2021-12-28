# script to calculate simulated FST_MaxSNP and FST_Window from ms output for 2 populations
### accepts neutral or adaptive scenarios, with or without migration, complete or incomplete
### complete sweep scenarios need to be simulated with 2 extra samples.
### incomplete sweep scenarios need input of the min and max allele threshold for the target final frequency
# tribeiro@wisc.edu, tiaaago@gmail.com

# arguments:
# arg 1: input file
# arg 2: *minimum ending allele frequency accepted (in number of alleles) # *or 'migr' if simulating migration
# arg 3: maximum ending allele frequency accepted (in number of alleles)
# arg 4: sample size (can downsample to a value lower than the simulated sample size) - if simulating migration, add a random number to arg 3

# Examples:
# For complete sweeps, at least arg 1 and 2 need to be used: python3 ms_2pop_fst.py ms_output.txt complete
# For incomplete sweeps, at least arg 1, 2 and 3 need to be used: python3 ms_2pop_fst.py ms_output_finalfreq50.txt 23 27
# For migration, at least arg 1 and 2 need to be used: python3 ms_2pop_fst.py ms_output.txt migr
# For sample size, all args need to be used (arg 2 and 3 will be considered if relevant): python3 ms_2pop_fst.py ms_output.txt 1 1 30

import re, sys
from random import sample

# This function calculates Reynolds FST for any two populations, it is called inside the calc_fst_pbs function.
def fst(MinorFreq1, MinorFreq2, SampleSize1, SampleSize2):
#    "The arguments for this functions are: MinorFreqAllele from pop A and B (qA, qB), and the sample size from the respective populations."
    MinorFreq1 = float(MinorFreq1)
    MinorFreq2 = float(MinorFreq2)
    denom = MinorFreq1 + MinorFreq2 - 2 * MinorFreq1 * MinorFreq2
    SharedNum = SampleSize1 * (2 * MinorFreq1 - 2 * MinorFreq1 ** 2) + SampleSize2 * (2 * MinorFreq2 - 2 * MinorFreq2 ** 2)
    NumA = (MinorFreq1 - MinorFreq2) ** 2
    FracNum = (SampleSize1 + SampleSize2) * SharedNum
    FracDen = 4 * SampleSize1 * SampleSize2 * (SampleSize1 + SampleSize2 - 1)
    frac = FracNum / FracDen
    WholeNum = NumA - frac
    DenFracNum = (4 * SampleSize1 * SampleSize2 - SampleSize1 - SampleSize2) * SharedNum
    DenFrac = DenFracNum / FracDen
    WholeDen = NumA + DenFrac

    if WholeDen != 0:
        FST = WholeNum / WholeDen
    else:
        FST = 0

    return [FST, WholeNum, WholeDen]

# This function gets a list of output rows from 1 replicate of a ms_simulation, calculates allele freq for each SNP and calls an FST function
def sim_fst(ind_list,pop1,pop2,frag_size, samplesize=0):
    ind_pop1 = ind_list[:pop1]
    ind_pop2 = ind_list[pop1:]
    if samplesize != 0:  
        ind_pop1 = ind_pop1[:samplesize]
        ind_pop2 = ind_pop2[:samplesize]
                
    fst_list = []
    MaxFSTSNP = -999
    SumWholeNum = 0
    SumWholeDen = 0

    for i in range(0,len(ind_list[0])):
        count1 = 0
        count2 = 0
        for j1 in range(0,len(ind_pop1)):
            count1 = count1 + ind_pop1[j1][i]
        freq1 = count1/pop1
        
        for j2 in range(0,len(ind_pop2)):
            count2 = count2 + ind_pop2[j2][i]
        freq2 = count2/pop2

        if samplesize != 0:
            freq1 = count1/samplesize
            freq2 = count2/samplesize
            fst_output = fst(freq1,freq2,samplesize,samplesize) # FST, WholeNum, WholeDen
        else:
            fst_output = fst(freq1,freq2,pop1,pop2) # FST, WholeNum, WholeDen
   
        SumWholeNum = SumWholeNum + fst_output[1]
        SumWholeDen = SumWholeDen + fst_output[2]
        if fst_output[0] > MaxFSTSNP:
            MaxFSTSNP = fst_output[0]

    if SumWholeDen != 0:            
        winFST = SumWholeNum/SumWholeDen
    else:
        winFST = 0

    return [MaxFSTSNP, winFST]

               

# function to obtain information about the simulation: complete/incomplete, selected loci kept/not kept, n1+n2, frag size
def first_row(full_str, f='fst'):
    first_line = full_str.split()
    compl_flag = 0 # if 0 = incomplete (n1=n2)
    n1_flag = 0 # used to flag if -I present
    pop_size_err = 0 # no -I found on command line, wont analyse this scenario
    smark_flag = 1 # flag if selected loci was kept
    
    for i in range(0,len(first_line)):
        if first_line[i] == '-I':
            sample1 = int(first_line[i+2])
            sample2 = int(first_line[i+3])
            len_fragm = int(first_line[i-1])
            n1_flag = 1

    if n1_flag == 0:
        pop_size_err = 1
    if '-Smark' not in first_line:
        smark_flag = 0
    if f == 'fst':
        if smark_flag == 1:
            if sample1 != sample2:
                sample1 = sample1 - 2
                compl_flag = 1

    return [compl_flag, smark_flag, sample1, sample2, len_fragm, pop_size_err]

# if applicable, this function will find the location of the beneficial site (which should be in the middle of the window, position=0.50000'. 
# if there is more than one position 0.50000 it will try to detect which will inform all of them (to be checked in the newcol() function)
def pos_to_check(full_str):
    positions = full_str.split()
    pos2check = positions.count('0.50000')
    if pos2check == 1:
        position_list = []
        for sl in range(0,len(positions)):
            if positions[sl] == '0.50000':
                position_list.append(sl-1) # 1st position is a str
                break
    elif pos2check > 1:
        position_list = []
        for sl in range(0,len(positions)):
            if positions[sl] == '0.50000':
                position_list.append(sl-1) 
    else:
        position_list.append("error, no position to check")

    return position_list

# decides which sample to remove from complete sweeps
def completesweep_ind2rm(locus, sample1, sample2,p2check,rmid=-9, mcheck='no_migr'):
    ind2rm = []
    if p2check == 1:
        s1_count = sum(locus[:(sample1+2)])
        
        if mcheck == 'no_migr':
            if s1_count == sample1+2:
                ind2rm = sample(range(0,sample1+2),2) # remove any 2
            elif s1_count == sample1+1:
                samp = list(range(0,sample1+2))
                for find_0 in range(0,sample1+2):
                    if locus[find_0] == 0:
                        ind2rm.append(find_0) # remove the one with a 0
                        samp.remove(find_0)
                        ind2rm.append(sample(samp,1)[0]) # remove any 1 additionally
            elif s1_count == sample1:
                for find_0 in range(0,sample1+2):
                    if locus[find_0] == 0:
                        ind2rm.append(find_0) # remove the ones with a 0
            else:
                ind2rm.append(-999)

            return ind2rm
        else:   
            s1_count = sum(locus[:(sample1)])      
            s2_count = sum(locus[(sample1):])
            return [s1_count/(sample1),s2_count/(sample2)] # if migr flag is on, will return the freqs of both pops instead of a list of inds to remove
        
    elif p2check == 2:
        if mcheck == 'no_migr':
            if sum(locus[rmid][:(sample1+2)]) == sample1+2:
                ind2rm = sample(range(0,sample1+2),2) # remove any 2
            elif sum(locus[rmid][:(sample1+2)]) == sample1+1:
                samp = list(range(0,sample1+2))
                for find_0 in range(0,sample1+2):
                    if locus[rmid][find_0] == 0:
                        ind2rm.append(find_0) # remove the one with a 0
                        samp.remove(find_0)
                        ind2rm.append(sample(samp,1)[0]) # remove any 1 additionally
            elif sum(locus[rmid][:(sample1+2)]) == sample1:
                for find_0 in range(0,sample1+2):
                    if locus[rmid][find_0] == 0:
                        ind2rm.append(find_0) # remove the ones with a 0
            else:
                ind2rm = [-999]

            return ind2rm
        else:
            s1_count = sum(locus[rmid][:(sample1)])
            s2_count = sum(locus[rmid][(sample1):])
            return [s1_count/(sample1),s2_count/(sample2)] # if migr flag is on, will return the freqs of both pops instead of a list of inds to remove

# This function isolates the locus targeted by selection to verify if it meets the criteria of simulated scenario.
# For complete sweeps, which has 2 extra samples, at least total samples - 2 need to have the beneficial allele.
# For incomplete sweeps, number of samples with the beneficial allele need to match the min and max allele boundaries for the target final frequency
# if there is more than 1 possible location for the beneficial allele, this function will try to identify which one is correct and discard the replicate if that is not possible.
def newcol(locus,p2check):
    if p2check == 1:
        n_column = []

        # ensures multiple mutations are all recorded as 1
        for allele in locus:
            if allele != '0':
                n_column.append(1)
            else:
                n_column.append(0)
        return n_column
    elif p2check == 2: #len(locus) > 1
        s_flag = []
        multi_column = []
        for base_list in range(0,len(locus)):
            n_column = []
            s_list = 0
            for allele in locus[base_list]:
                if allele != '0':
                    n_column.append(1)
                else:
                    n_column.append(0)
                if allele != '0' and allele != '1':
                    s_list = 1
            multi_column.append(n_column)

            if s_list == 1:
                s_flag.append(base_list)

        if len(s_flag) != 1: # no clear selected loci identified
            multi_column = -999

        return [multi_column,s_flag]
            


###############################################
##################### FST #####################
###############################################                    

if sys.argv[2] == 'migr':
    migr_check = 'migr'
else:
    migr_check = 'no_migr'
file = sys.argv[1]

row_count = 0
read_sim = 0
sim_start = 0
neutral_check = 0

sim_rows = []
selected_locus_list = []
selected_locus_listN = []

output_detail = open('fst_'+sys.argv[1],'w')
output_detail.write('MaxSNPFST\tWindowFST\n')

n1 = 0
n2 = 0
length_fragment_ms = 0
n_found = 0

if len(sys.argv) == 5:
    ss = int(sys.argv[4])
else:
    ss = 0
    
with open(file) as sim:
    for l in sim:
        if n_found == 0: # Only runs here on the first line, find the -I flag and get the sample sizes.
            first_r_info = first_row(l) # compl_flag, smark_flag, sample1, sample2, len_fragm, pop_size_err
            if first_r_info[5] == 1: # 5 = pop_size_err
                print("Population size not found.")
                break
            else:
                n_found = 1
            n1 = first_r_info[2]
            n2 = first_r_info[3]
            length_fragment_ms = first_r_info[4]
                
        if read_sim == 1: # start reading the individuals
            row = [c for c in l[:-1]]
            
            if len(row) > 1: # read one individual at a time and store in a list
                if first_r_info[1] == 1: #Smark is present. Selected loci kept. Store info about the individual, but keeping selected loci separately.
                    if len(selected_pos_list) == 1:
                        selected_pos = selected_pos_list[0]
                        selected_locus_list.append(row[selected_pos])
                        row = row[:selected_pos] + row[selected_pos+1:]
                        row = list(map(int,row))
                    elif len(selected_pos_list) > 1:
                        if len(selected_locus_listN) == 0:
                            for sl_loop in range(0,len(selected_pos_list)):
                                selected_locus_listN.append([])
                        for sl_loop in range(0,len(selected_pos_list)):
                            selected_locus_listN[sl_loop].append(row[selected_pos_list[sl_loop]])

                        row = row[:selected_pos_list[0]] + row[selected_pos_list[-1]+1:]
                        row = list(map(int,row))
                else:
                    row = list(map(int,row))
                sim_rows.append(row)
                
            # empty row -> all individuals from this replicate were read. Run SNP/FST analyses. Save. Reset.
            else:
                if len(sim_rows) == (n1+n2+2) or len(sim_rows) == (n1+n2): # right number of individuals probably correct for this simulation. Run the analysis

                    # Detect the kind of sweep and prepare individuals for FST calculation if all pre-reqs are met.
                    
                    if first_r_info[1] == 1: #smark is present. Selected loci kept. Transform selected loci in 0 and 1
                        discard_check = 0 #selection was simulated, check whether meet criteria to discard or not.
                        if len(selected_pos_list) == 1:
                            
                            new_column = newcol(selected_locus_list,1) # ensures multiple mutations are all recorded as 1



                            if first_r_info[0] == 1  and len(sim_rows) == (n1+n2+2): #Complete sweep simulated and right number of individuals found in this simulation. Run the analysis
                                ind_to_remove = completesweep_ind2rm(new_column,n1,n2,1, mcheck=migr_check)

                            elif migr_check == 'migr':
                                ind_to_remove = completesweep_ind2rm(new_column,n1,n2,1, mcheck=migr_check)

                            elif first_r_info[0] == 0  and len(sim_rows) == (n1+n2): #Incomple sweep and enough simulations within the range of the end allele freq
                                if sum(new_column[:n1]) >= int(sys.argv[2]) and sum(new_column[:n1]) <= int(sys.argv[3]):
                                    ind_to_remove = [-1]

                                   
                                else:
                                    ind_to_remove = [-999]
                                    
                            else:
                                discard_check = 1 # not enough individuals in this simulation
                            
                                
                            if ind_to_remove[0] == -999:
                                discard_check = 1 # Not enough individuals with the selected loci simulated
                            else:     
                                for paste_05000_index in range(0,len(new_column)):
                                    sim_rows[paste_05000_index].append(int(new_column[paste_05000_index]))

                                subsample_sim_rows = []
                                for rmving in range(0,len(sim_rows)):
                                    if rmving not in ind_to_remove:
                                        subsample_sim_rows.append(sim_rows[rmving])

                            
                        else:

                            columns_ind2rm = newcol(selected_locus_listN,2) # new columns and columns with selected loci
                            new_columnN = columns_ind2rm[0]


                            
                            if new_columnN == -999:
                                # No clear selected loci
                                # print('a')
                                discard_check = 1 # more than 1 possible selected loci
                            else:
                                #rm indv
                                flag_id = columns_ind2rm[1] # id of column with selected loci
                                

                                
                                for b_list in range(0,len(new_columnN)):
                                    for ind_index in range(0,len(new_columnN[b_list])):
                                        sim_rows[ind_index].append(new_columnN[b_list][ind_index])



                                if first_r_info[0] == 1  and len(sim_rows) == (n1+n2+2): #Complete sweep simulated and right number of individuals found this simulation. Run the analysis
                                    ind_to_remove = completesweep_ind2rm(new_columnN,n1,n2,2,flag_id[0],mcheck=migr_check)
                                elif first_r_info[0] == 0  and len(sim_rows) == (n1+n2): #Incomple sweep
                                    if sum(new_columnN[flag_id[0]][:n1]) >= int(sys.argv[2]) and sum(new_columnN[flag_id[0]][:n1]) <= int(sys.argv[3]):
                                        ind_to_remove = [-1]

                                       
                                    else:
                                        ind_to_remove = [-999]
                                    
                                else:
                                    discard_check = 1 # not enough individuals in this simulation

                                if ind_to_remove[0] == -999:
                                    discard_check = 1 # Not enough individuals with the selected loci simulated
                                elif migr_check == 'no_migr':
                                    subsample_sim_rows = []
                                    for rmving in range(0,len(sim_rows)):
                                        if rmving not in ind_to_remove:
                                            subsample_sim_rows.append(sim_rows[rmving])
                                    
                    else:
                        discard_check = 0 
                        neutral_check = 1 # neutral simulation. not a sweep with select loci recorded
                else:

                    discard_check = 1 # not enough individuals in this simulation


                    
                if discard_check == 0:
                    if neutral_check == 1: # if it is neutral

                        fst_sim_replicate = sim_fst(sim_rows,n1,n2,length_fragment_ms, ss) # SNP \t Window FST
                                  
                    else: # complete or incomplete sweeps in general

                        if len(subsample_sim_rows) == n1 + n2:
                            fst_sim_replicate = sim_fst(subsample_sim_rows,n1,n2,length_fragment_ms, ss)
                        else:
                            fst_sim_replicate = sim_fst(sim_rows,n1,n2,length_fragment_ms, ss)
                            
                    if migr_check == 'migr':
                        fst_sim_replicate = fst_sim_replicate + ind_to_remove

                        
                    outting_detail = list(map(str,fst_sim_replicate))
                    output_detail.write('\t'.join(outting_detail)+'\n')

                read_sim = 0
                sim_rows = []
                selected_locus_list = []
                selected_locus_listN = []

        elif 'positions' in l: # When 'positions' is in the row, that means the row has the individuals.
            # row with SNP positions and potential 0.50000 -Smark. If no -Smark used, just move on to read indvs
            if first_r_info[1] == 1: # -Smark flag used in ms command
                selected_pos_list = pos_to_check(l) # list of positions to check
            read_sim = 1

output_detail.close()
