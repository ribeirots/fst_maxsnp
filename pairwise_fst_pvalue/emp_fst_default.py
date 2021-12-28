# Script to calculate pi and fst from counts file using a cutoff for minimum sample size (read count) included in the code.
# use AllCounts files. Downsamples as needed.
# To generate files to be used against def
# Author: Tiago Ribeiro - tiaaagosr@gmail.com - March 2021

import re, random, sys

# arg 1: pop1
# arg 2: pop2
# arg 3: arm

import numpy as np

site_counts = []
min_allele_count = 3 # combined counts from both population - only used for FST and dxy

pops = [sys.argv[1], sys.argv[2]] # 'CO','EA', and 'EG' not used due to low sample size
chrm = [sys.argv[3]]

cutoff_dict = { \
'CO': {'X': 9, '2R': 9, '3L':9}, \
'EA': {'X': 13, '2R': 5, '3L':6}, \
'EF': {'X': 41, '2R': 25, '3L':10, '2L':29, '3R':13}, \
'EG': {'X': 18, '2R': 15, '3L':10}, \
'FRht': {'X': 89, '2R': 101, '3L': 86, '2L':87, '3R':63}, \
'KF': {'X': 15, '2R': 26, '3L': 24, '2L':28, '3R':30}, \
'RG': {'X': 21, '2R': 21, '3L':24, '2L':18, '3R':19}, \
'SD': {'X': 40, '2R': 14, '3L':17, '2L':12, '3R':9}, \
'ZI': {'X': 181, '2R': 158, '3L': 151, '2L':144, '3R':129}}

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

    if FST > 1:
        print('Fst larger than 1.')

    return [str(FST), str(WholeNum), str(WholeDen)]

def down_sample(count_4, target_n): # A C G T, and sample size
    all_bases = 'A' * count_4[0] + 'C' * count_4[1] + 'G' * count_4[2] + 'T' * count_4[3]
    all_bases = random.sample(all_bases,target_n)
    return_count = []
    return_count.append(all_bases.count('A'))
    return_count.append(all_bases.count('C'))
    return_count.append(all_bases.count('G'))
    return_count.append(all_bases.count('T'))

    return return_count

for arm in chrm:
    for p1 in range(0,len(pops)):
        if pops[p1] == 'FRht':
            filename = '../invFreeCounts/AllCounts_'+pops[p1]+'_NoInv_Chr'+arm+'_diploids.txt'
        else:
            filename = '../invFreeCounts/AllCounts_'+pops[p1]+'_NoInv_Chr'+arm+'.txt'
        sites_p1 = [] # In this vector, sites with value 0 = not used in pop1
        p1_pi_sites_used = 0 # counts how many sites were used to calculate pi_pop1
        pi_pop1_sum = 0 # stores the sum of 2pq site by site for pop1
        with open(filename) as counts_file:
            print(filename)
            for line in counts_file:
                line = re.split('\t',line[:-1])
                line = list(map(int,line))
                sample_size_full = sum(line)
                if sample_size_full >= cutoff_dict[pops[p1]][arm]:
                    sites_p1.append(line) # for future use on fst_dxy

                    # pi calculation por pop1
                    p1_pi_sites_used += 1
                    p1_pi_min_max = sorted(range(len(line)), key=lambda k: line[k])
                    p1_pi_line = [0,0,line[p1_pi_min_max[2]],line[p1_pi_min_max[3]]]

                    if sum(p1_pi_line) > cutoff_dict[pops[p1]][arm]:
                        pi_site = down_sample(p1_pi_line,cutoff_dict[pops[p1]][arm])
                    else:
                        pi_site = p1_pi_line[:]
                    p1_pi_minor_freq = float(pi_site[2])/sum(pi_site)
                    p1_pi_major_freq = float(pi_site[3])/sum(pi_site)

                    pi_pop1_sum += 2 * p1_pi_minor_freq * p1_pi_major_freq
                else:
                    sites_p1.append([0])
        pi_p1_final = float(pi_pop1_sum)/p1_pi_sites_used
        pi_p1_output = 'pi\tChr'+arm+'\t'+pops[p1]+'\t'+str(pi_p1_final)+'\n'


        if p1 == len(pops)-1: # pi for last pop was calculated, all pairwise comparisons have been calculated already
            break
        
        for p2 in range(p1+1,len(pops)):
            line_count_p2 = -1
            if pops[p2] == 'FRht':
                filename2 = '../invFreeCounts/AllCounts_'+pops[p2]+'_NoInv_Chr'+arm+'_diploids.txt'
            else:
                filename2 = '../invFreeCounts/AllCounts_'+pops[p2]+'_NoInv_Chr'+arm+'.txt'

            outfile_snp_fst = open('empirical_snp_fst_'+pops[p1]+'_'+pops[p2]+'_NoInv_Chr'+arm+'downsampled.txt','w')
            snp_fst_header = 'Position\tSNP_FST\tWholeNum\tWholeDen\n'
            outfile_snp_fst.write(snp_fst_header)
            
            with open(filename2) as counts_file2:
                print('Pop2: ' + pops[p2])
                for line2 in counts_file2:
                    line_count_p2 += 1
                    line2 = re.split('\t',line2[:-1])
                    line2 = list(map(int,line2))
                    sample_size2_full = sum(line2)
                    if sample_size2_full >= cutoff_dict[pops[p2]][arm] and len(sites_p1[line_count_p2]) == 4:
                        min_max_skip = 0
                        l1 = sites_p1[line_count_p2]
                                                                            
                        combined_l = []                        
                        combined_l.append(l1[0]+line2[0])
                        combined_l.append(l1[1]+line2[1])
                        combined_l.append(l1[2]+line2[2])
                        combined_l.append(l1[3]+line2[3])
                        min_max = sorted(range(len(combined_l)), key=lambda k: combined_l[k])

                        l1[min_max[0]] = 0
                        l1[min_max[1]] = 0
                        line2[min_max[0]] = 0
                        line2[min_max[1]] = 0

                        # downsampling to ensure same sampling across all sites and in the simulation
                        if sum(line2) == cutoff_dict[pops[p2]][arm]:
                            line2 = line2
                        elif sum(line2) > cutoff_dict[pops[p2]][arm]:
                            line2 = down_sample(line2,cutoff_dict[pops[p2]][arm])
                        else:
                            min_max_skip = 1

                        if sum(l1) == cutoff_dict[pops[p1]][arm]:
                            l1 = l1
                        elif sum(l1) > cutoff_dict[pops[p1]][arm]:
                            l1 = down_sample(l1,cutoff_dict[pops[p1]][arm])
                        else:
                            min_max_skip = 1

                        if min_max_skip == 0:
                            minor1 = l1[min_max[2]]
                            major1 = l1[min_max[3]]
                            sample_size1 = minor1 + major1
                            
                            minor2 = line2[min_max[2]]
                            major2 = line2[min_max[3]]
                            sample_size2 = minor2 + major2
                            
                            min_freq1 =  float(minor1)/sample_size1
                            maj_freq1 = 1 - min_freq1
                            

                            min_freq2 =  float(minor2)/sample_size2
                            maj_freq2 = 1 - min_freq2
                            
                            if (minor1 + minor2) >= min_allele_count:
                                snp_fst_output = [str(line_count_p2)] + fst(min_freq1, min_freq2, sample_size1, sample_size2)

                                outfile_snp_fst.write('\t'.join(snp_fst_output)+'\n') # Position, SNP_FST, WholeNum, WholeDen

            outfile_snp_fst.close()                



