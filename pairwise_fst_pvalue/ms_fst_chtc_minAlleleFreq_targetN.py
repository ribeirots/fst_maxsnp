# script to calculate snp and window FST from ms sim, with minAllele freq and a fixed N per pop
# tribeiro@wisc.edu

from random import sample
# args:
# arg 1: input file for multipop (full path)
# arg 2: path to output file for multipop
# arg 3: output file
# arg 4: list of target N for pops, in order. Separated with comma and no spaces. E.g. 158,21,4,25,101
import re, sys

# Model NA_2R: ZI        RG        CO        EF        FR
####### Chr2R: 158, 21, 4, 25, 101
####### Chr2L: 144, 18, 4, 25, 87 *downsample empirical EF 2L from 29 to 25 to compare against NA_2R model.
####### Chr3L: 151, 21, 4, 10, 86 *downsample empirical RG 3L from 24 to 21 to compare against NA_2R model.
####### Chr3R: 129, 19, 4, 13, 63

# Model NA_3L: ZI, RG, CO, EF, EG, FR
####### Chr3L: 151, 24, 4, 10, 4, 86
####### Chr2L: 144, 18, 4, 10, 4, 86 *downsample empirical (EF 2L from 29 to 10) and (FR 2L from 87 to 86) to compare against NA_3L model.
####### Chr2R: 151, 21, 4, 10, 4, 86 *downsample empirical (ZI 2R from 158 to 151) and (EF 2R from 25 to 10) and (FR 2R from 101 to 86) to compare against NA_3L model.
####### Chr3R: 129, 19, 4, 10, 4, 63 *downsample empirical EF 3R from 13 to 10 to compare against NA_3L model.

# Model SA_2R: ZI RG CO EA KF SD
####### Chr2R: 158, 21, 4, 4, 26, 14
####### Chr2L: 144, 18, 4, 4, 26, 12 *downsample empirical (KF 2L from 28 to 26) to compare against SA_2R model.
####### Chr3L: 151, 21, 4, 4, 24, 14 *downsample empirical (RG 3L from 24 to 21) and (SD 3L from 17 to 14) to compare against SA_2R model.
####### Chr3R: 129, 19, 4, 4, 26, 9 *downsample empirical (KF 3R from 30 to 26) to compare against SA_2R model.

# Model SA_3L: ZI RG CO EA KF SD
####### Chr3L: 151, 24, 4, 4, 24, 17
####### Chr2L: 144, 18, 4, 4, 24, 12 *downsample empirical (KF 2L from 28 to 24) to compare against NA_2R model.
####### Chr2R: 151, 21, 4, 4, 24, 14 *downsample empirical (ZI 2R from 158 to 151) and (KF 2R from 26 to 24) to compare against NA_2R model.
####### Chr3R: 129, 19, 4, 4, 24, 9 *downsample empirical (KF 3R from 30 to 24) to compare against NA_2R model.

target_Npops = sys.argv[4]
target_n_multipops = [int(numb) for numb in re.split(',',target_Npops)]

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

#    if FST >= 0.99:
#        FST = 0.99 # Maximum FST
#    elif FST < 0:
#        FST = 0

    return [FST, WholeNum, WholeDen]

# This function gets a list of output rows from 1 replicate of a ms_simulation, calculates allele freq for each SNP and calls an FST function
def sim_fst(ind_list,pop1,pop2,frag_size, n_multipops=-999,z=-999,y=-999, targetN = target_n_multipops): #z and y are pop ids for multipop. They are only use to return/print the pop_id correctly
    min_allele_freq = 3
    if n_multipops == -999:
        ind_pop1 = ind_list[:pop1]
        ind_pop2 = ind_list[pop1:]
        fst_list = []
        pi_1_sum = 0
        pi_2_sum = 0
        MaxFSTSNP = -999
        SumWholeNum = 0
        SumWholeDen = 0
        dxy_sum = 0

        for i in range(0,len(ind_list[0])):
            count1 = 0
            count2 = 0
            for j1 in range(0,len(ind_pop1)):
                count1 = count1 + ind_pop1[j1][i]
            freq1 = count1/pop1
            
            pi_1_sum = pi_1_sum + 2*freq1*(1-freq1)
            
            for j2 in range(0,len(ind_pop2)):
                count2 = count2 + ind_pop2[j2][i]
            freq2 = count2/pop2
            pi_2_sum = pi_2_sum + 2*freq2*(1-freq2) # to mimic analyze_ms_data.pl, mean pairwise diff, should be just 2*q1

            if count1 + count2 >= min_allele_freq:

                dxy_count = (freq1 * (1-freq2)) + (freq2 * (1-freq1))
                dxy_sum = dxy_sum + dxy_count
                fst_output = fst(freq1,freq2,pop1,pop2) # FST, WholeNum, WholeDen
                SumWholeNum = SumWholeNum + fst_output[1]
                SumWholeDen = SumWholeDen + fst_output[2]
                if fst_output[0] > MaxFSTSNP:
                    MaxFSTSNP = fst_output[0]

                
        if MaxFSTSNP == -999:
            winFST  = 'NA'
            MaxFSTSNP = 'NA'
            dxy_return = 'NA'
        elif SumWholeDen != 0:            
            winFST = SumWholeNum/SumWholeDen
            dxy_return = dxy_sum/frag_size
        else:
            winFST = 0
            dxy_return = dxy_sum/frag_size
            
        pi1 = pi_1_sum/frag_size
        pi2 = pi_2_sum/frag_size


        if z == -999:
            return [MaxFSTSNP, winFST, pi1, pi2]
        else:
            return [MaxFSTSNP, winFST, pi1, pi2, dxy_return, z, y] # return for multipops
    
    else: # multipops detecteds
        p1_start = 0
        p1_end = 0
        locus_length = frag_size
        list_of_fsts = []
        for p1_i in range(0,len(n_multipops)-1):
            p1_start = p1_end
            p1_end = p1_end + n_multipops[p1_i]
            p2_end = p1_end
            for p2_i in range(p1_i+1,len(n_multipops)):                
                p2_start = p2_end
                p2_end = p2_start + n_multipops[p2_i]
                indv_pop1 = ind_list[p1_start:p1_end]
                indv_pop2 = ind_list[p2_start:p2_end]

                # Downsample to target N # targetN[p1_i] targetN[p2_i]
                downsamp_indv_pop1 = sample(indv_pop1,targetN[p1_i])
                downsamp_indv_pop2 = sample(indv_pop2,targetN[p2_i])

                new_indv_list = downsamp_indv_pop1 + downsamp_indv_pop2
                list_of_fsts.append(sim_fst(new_indv_list, targetN[p1_i], targetN[p2_i],locus_length,z=p1_i,y=p2_i))
        return list_of_fsts            
              
            
                
                

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
            s2_count = sum(locus[(sample1+2):])
            return [s1_count/(sample1+2),s2_count/(sample2)] # if migr flag is on, will return the freqs of both pops instead of a list of inds to remove
        
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
            s1_count = sum(locus[rmid][:(sample1+2)])
            s2_count = sum(locus[rmid][(sample1+2):])
            return [s1_count/(sample1+2),s2_count/(sample2)] # if migr flag is on, will return the freqs of both pops instead of a list of inds to remove

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

        if len(s_flag) == 1: # return condition for sys.argv[0] == 'freq', with clear loci.
            return multi_column[s_flag[0]]
        else:
            return [-999]
        
###############################################
#### MULTIPLE POPULATIONS #####################
###############################################

if 1:
    with open(sys.argv[1]) as sim:
        first_row_flag = 1
        new_sim_found = 0
        pre_sim_row_count = 0
        sim_start = 0
        repl_ind_list = []
        frag_length = 0
        repl_counts = 0
        sums_calcs = []
        outting_dets = open(sys.argv[2]+'details_'+sys.argv[3],'w')
#        det_head = 'MaxFSTSNP\twinFST\tpi1\tpi2\tdxy\tpop1\tpop2\n'
        det_head = 'MaxFSTSNP\twinFST\tpop1\tpop2\n'
        outting_dets.write(det_head)
        
        for l in sim:
            # only goes into this if in the first row, to get number of pops and their sample sizes
            if first_row_flag == 1:
                first_row_flag = 0
                first_line = l.split()
                for i in range(0,len(first_line)):
                    if first_line[i] == '-I':
                        frag_length = int(first_line[i-1])
                        I_index = i
                        npop = int(first_line[i+1])
                pop_sample_sizes = []
                for j in range(0,npop):
                    pop_sample_sizes.append(int(first_line[I_index+2+j]))
                for np1 in range(0,len(pop_sample_sizes)-1):
                    for np2 in range(np1+1,len(pop_sample_sizes)):
                        sums_calcs.append([0,0,0,0,0,0,0])

                        
            # goes into this part when detects new replicate and skips two rows until the 1st individual
            elif len(l) == 1 or new_sim_found == 1 or 'ms' in l:
                new_sim_found = 1
                repl_ind_list = []

                if 'positions' in l:
                    repl_ind_list = []
                    new_sim_found = 0
                    sim_start = 1
                    

            # goes into this part to collect all individuals from a replicate in a list
            elif sim_start == 1:
                indv = [c for c in l[:-1]]
                if '/' in indv:
                    print(indv)
                indv = list(map(int,indv))
                repl_ind_list.append(indv)
                
            if len(repl_ind_list) == sum(pop_sample_sizes):
                two_by_two_fst = sim_fst(repl_ind_list, 0, 0, frag_length, n_multipops = pop_sample_sizes)

                outting_dets.write('new_replicate\n')
                for poppair in two_by_two_fst:
                    poppair = list(map(str,poppair))
                    outting_dets.write('\t'.join(poppair)+'\n')
        
    outting_dets.close()


