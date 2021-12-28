# Script to find overlapping regions. The inputs are 2 files with windows (ideally outlier/candidate regions). The first three columns have to be Chrm, Window_Start, Window_End - with a header.
# arg 1: input focal list
# arg 2: input list to be checked against
# author: Tiago Ribeiro - tribeiro@wisc.edu

import re, sys

# reads 1st input file
regions1 = []
with open(sys.argv[1]) as reg1:
    next(reg1)
    for row in reg1:
        row = re.split('\t',row[:-1])
        row[1] = int(row[1])
        row[2] = int(row[2])
        row2 = row[:3]
        regions1.append(row2)


# reads 2nd input file
regions2 = [] # base list of windows, comparisons will be made against this list.
with open(sys.argv[2]) as reg2:
    next(reg2)
    for row in reg2:
        row = re.split('\t',row[:-1])
        row[1] = int(row[1])
        row[2] = int(row[2])
        row2 = row[:3]
        regions2.append(row2)

# function to find windows in 2 that overlap with windows in 1.
def over_windows(r1, r2_list=regions2): # r1 =  the focus window, will return all windows that overlap with it. r2_list: list of windows to be scanned searching for overlaps with r1.
    r2_return = []
    for r2 in r2_list:
        if (r1[0] == r2[0]) and (r1[2] >= r2[1]) and (r1[1] <= r2[2]): # end1 >= start2 and start1 <= end2 -> there is at least 1 base overlap.
            # check if overlap size is at least 50% the size of the  smaller region
            half_len_r1 = (r1[2] - r1[1] + 1)/2
            half_len_r2 = (r2[2] - r2[1] + 1)/2
            len_overlap = r1[2] - r2[1] + 1
            if (len_overlap >= half_len_r1) or (len_overlap >= half_len_r2):
                r2_return.append(r2)
    return r2_return

# output creation
output_id = re.split('\.',sys.argv[1])[0]
out_list1 = open(output_id+'_overlap_list1.txt','w')
out_list2 = open(output_id+'_overlap_list2.txt','w')
out_2in_list1 = open(output_id+'_overlap2_in_list1.txt','w')

# create dictonary of regions from reg2 that overlap reg1.
reg2inreg1_dict = {}
list2_to_output = []
for reg1 in regions1:
    r1_overs = over_windows(reg1)
    if len(r1_overs) > 0:
        reg1out = list(map(str,reg1))
        out_list1.write('\t'.join(reg1out)+'\n')

        full_out = []
        reg1_full_out = '_'.join(reg1out)+':'
        full_out.append(reg1_full_out)

        for obs in r1_overs:
            obs2 = list(map(str,obs))
            obs2_full_out = '_'.join(obs2)
            full_out.append(obs2_full_out)
            if obs not in list2_to_output:
                list2_to_output.append(obs)
                out_list2.write('\t'.join(obs2)+'\n')
        
        full_out = '\t'.join(full_out)+'\n'
        out_2in_list1.write(full_out)

out_list1.close()
out_list2.close()
out_2in_list1.close()



