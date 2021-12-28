# generate msms simulation commands based on empirical window length and recombination
# output bash file to be ran on CHTC
# tribeiro@wisc.edu March2021

# argument 1 = id for the current bash file being generated.
# argument 2 = how many jobs will be submitted. # of simulations per job is given by the script below.
# argument 3 = simulations per job
# argument 4 = recombination rate code based on empirical files

import random, sys

# Write here the name of the model being simulated: South 2R

n_pops = 6
# Ne values for SA and NA models
# Southern Extension, X: 661372.5, 2R: 826813, 3L: 597868.
# Northern Extension: X: 618139, 2R: 823306.5, 3L: 827211.
Ne = 826813
ploidy = 4 # change to 3 for X. change to 4 for autosomes.
n_simulations = int(sys.argv[3]) # use 1000 on chtc
r_rate = sys.argv[4]
geneconversion = 6.25e-8
geneconv_tractlength = 518
perbase_theta = 0.0083 # NA_3L: 0.0077, NA_X_EG_FR = 0.0069, NA_X_EF = 0.00711, NA_2R = 0.0078, SA_X = 0.0077, SA_2R = 0.0083, SA_3L = 0.0056

recrate = []
if ploidy == 4:
    with open('../../recrate'+r_rate+'_auto_lengths.txt') as recwindows:
        next(recwindows)
        for recraterow in recwindows:
            recraterow = recraterow.split() # Arm, Start, End, Length, c
            rec_l = recraterow[3:]
            
            rec_l[0] = int(rec_l[0]) # len
            rec_l[1] = float(rec_l[1]) # c cM/MB. c = 100 -> 1 Morgan/megabase, 1 crossing-over expected every 1mb
            rec_l.append(rec_l[1]/100000000) # 100000000 -> divide by 100 to convert cM to M. Divide by 1m to convert from Comeron's scale (1MB) to per base.
            rec_l[1] = rec_l[2]*rec_l[0]*Ne*ploidy # multiply by 4N to obtain the population-scaled metric, multiply by rec_l[0] (length) to obtain the locus value
            if rec_l[1] == 0:
                rec_l[2] = ploidy * Ne * geneconversion
            recrate.append(rec_l)
else:
    with open('../../recrate'+r_rate+'_X_lengths.txt') as recwindows:
        next(recwindows)
        for recraterow in recwindows:
            recraterow = recraterow.split() # Arm, Start, End, Length, c
            rec_l = recraterow[3:]
            
            rec_l[0] = int(rec_l[0]) # len
            rec_l[1] = float(rec_l[1]) # c cM/MB. c = 100 -> 1 Morgan/megabase, 1 crossing-over expected every 1mb
            rec_l.append(rec_l[1]/100000000) # 100000000 -> divide by 100 to convert cM to M. Divide by 1m to convert from Comeron's scale (1MB) to per base.
            rec_l[1] = rec_l[2]*rec_l[0]*Ne*ploidy # multiply by 4N to obtain the population-scaled metric, multiply by rec_l[0] (length) to obtain the locus value
            if rec_l[1] == 0:
                rec_l[2] = ploidy * Ne * geneconversion
            recrate.append(rec_l)
            
recrate = random.choices(recrate, k=n_simulations) # sampling with replacement

# Start of the bash script to be submitted
bash_output = open('ms_rec_'+sys.argv[1]+'.sh','w')
first_row = '#!/bin/bash\n'
bash_output.write(first_row)
second_row = '\n'
bash_output.write(second_row)

# untar python3
untar = 'tar -xzf python37.tar.gz\n'

exp1 = 'export PATH=$PWD/python/bin:$PATH\n'
exp2 = 'export PYTHONPATH=$PWD/packages\n'
exp3 = 'export HOME=$PWD\n'

bash_output.write(untar)
bash_output.write(exp1)
bash_output.write(exp2)
bash_output.write(exp3)
bash_output.write('\n')


mkdir_row = 'mkdir ms_output\n'
bash_output.write(mkdir_row)
bash_output.write('\n')

bash_output.write('tar -xzf ms_chtc.tar.gz\n')
bash_output.write('cd msdir\n')

bash_output.write('\n')
        
# loop here to run a few simulations per job, always changing the recombination rate and length
for rec_i in range(0,n_simulations):    
    c_size = recrate[rec_i] # randomly obtained from file (with replacement)

    c_size_ms = str(c_size[1])+' '+str(c_size[0])
    
    #gene convertion
    if c_size[1] == 0:
        gn_f = c_size[2]
    else:
        gn_f = geneconversion/c_size[2]

    theta =  perbase_theta * c_size[0] # theta/bp * length
    
    gn_tc = geneconv_tractlength

    gn_code = str(gn_f)+' '+str(gn_tc)
    
    output_ms = 'results_ms'+str(rec_i)+'.txt'

    cmd = './ms 244 1 -t '+str(theta)+' -r '+c_size_ms+' -c '+gn_code+' -I 6 158 24 9 6 30 17 0 -en 0 1 2.329934943 -en 0 2 0.32654905 -en 0 3 0.462832587 -en 0 4 0.273434259 -en 0 5 0.450255983 -en 0 6 2.179843568 -em 0 1 2 2.86E+01 -em 0 2 1 2.86E+01 -em 0 4 1 0 -em 0 4 1 0 -em 0 1 3 1.08E-02 -em 0 3 1 1.08E-02 -em 0 1 5 4.71E+01 -em 0 5 1 4.71E+01 -em 0 1 6 4.47E+01 -em 0 6 1 4.47E+01 -em 0 2 3 2.00E+01 -em 0 3 2 2.00E+01 -em 0 2 4 1.57E+01 -em 0 4 2 1.57E+01 -em 0 2 5 0 -em 0 5 2 0 -em 0 2 6 0 -em 0 6 2 0 -em 0 3 4 1.16E-02 -em 0 4 3 1.16E-02 -em 0 3 5 0 -em 0 5 3 0 -em 0 3 6 0 -em 0 4 5 0 -em 0 5 4 0 -em 0 4 6 0 -em 0 6 4 0 -em 0 5 6 3.60E+01 -em 0 6 5 3.60E+01 -ej 0.004441905 6 1 -en 0.004441905001 1 0.542234459 -em 0.004441905001 1 2 4.73E+01 -em 0.004441905001 2 1 4.73E+01 -em 0.004441905001 1 3 1.07E+00 -em 0.004441905001 1 3 1.07E+00 -em 0.004441905001 1 4 0 -em 0.004441905001 4 1 0 -em 0.004441905001 1 5 4.69E+01 -em 0.004441905001 5 1 4.69E+01 -em 0.004441905001 2 3 2.00E+01 -em 0.004441905001 3 2 2.00E+01 -em 0.004441905001 2 4 1.57E+01 -em 0.004441905001 4 2 1.57E+01 -em 0.004441905001 2 5 0 -em 0.004441905001 5 2 0 -em 0.004441905001 3 4 1.16E-02 -em 0.004441905001 4 3 1.16E-02 -em 0.004441905001 3 5 0 -em 0.004441905001 5 3 0 -em 0.004441905001 4 5 0 -em 0.004441905001 5 4 0 -ej 0.017134013 4 3 -en 0.017134013001 3 0.255253002 -em 0.017134013001 1 2 4.73E+01 -em 0.017134013001 2 1 4.73E+01 -em 0.017134013001 1 3 5.46E+00 -em 0.017134013001 3 1 5.46E+00 -em 0.017134013001 1 5 4.69E+01 -em 0.017134013001 5 1 4.69E+01 -em 0.017134013001 2 3 1.75E-01 -em 0.017134013001 3 2 1.75E-01 -em 0.017134013001 2 5 0 -em 0.056132792001 5 2 0 -em 0.017134013001 3 5 0 -em 0.017134013001 5 3 0 -ej 0.04626303 3 2 -en 0.04626303001 2 0.356825546 -em 0.04626303001 1 2 9.12E-03 -em 0.04626303001 2 1 9.12E-03 -em 0.04626303001 1 5 4.69E+01 -em 0.04626303001 5 1 4.69E+01 -em 0.04626303001 2 5 0 -em 0.04626303001 5 2 0 -ej 0.047890212 1 2 -en 0.047890212001 2 0.282788248 -em 0.047890212001 2 5 2.08E+00 -em 0.047890212001 5 2 2.08E+00 -ej 0.066005403 5 2 -en 0.066005403001 2 1.011163346 -en 0.089409425 2 1 > ' + output_ms

    bash_output.write(cmd+'\n')
    bash_output.write('\n')

bash_output.write('cd ../\n')    
mvrow = 'mv msdir/results_ms*.txt ms_output\n'    
bash_output.write(mvrow)

cd_output = 'cd ms_output\n'
bash_output.write(cd_output)

bash_output.write('\n')
bash_output.write('for z in $(ls *.txt)\n')
bash_output.write('do\n')
bash_output.write('cd ../\n')
bash_output.write('python3 ms_fst_chtc_minAlleleFreq_targetN.py ms_output/$z ./ fst_output_chr2R.out 158,21,5,5,26,14\n')
bash_output.write('python3 ms_fst_chtc_minAlleleFreq_targetN.py ms_output/$z ./ fst_output_chr2L.out 144,18,4,4,28,12\n')
bash_output.write('python3 ms_fst_chtc_minAlleleFreq_targetN.py ms_output/$z ./ fst_output_chr3R.out 129,19,4,4,30,9\n')
bash_output.write('python3 ms_fst_chtc_minAlleleFreq_targetN.py ms_output/$z ./ fst_output_chr3L.out 151,24,4,4,24,17\n')
bash_output.write('cat details_fst_output_chr2R.out >> fst_recrate'+r_rate+'_'+sys.argv[1]+'_chr2R.csv\n')
bash_output.write('cat details_fst_output_chr2L.out >> fst_recrate'+r_rate+'_'+sys.argv[1]+'_chr2L.csv\n')
bash_output.write('cat details_fst_output_chr3R.out >> fst_recrate'+r_rate+'_'+sys.argv[1]+'_chr3R.csv\n')
bash_output.write('cat details_fst_output_chr3L.out >> fst_recrate'+r_rate+'_'+sys.argv[1]+'_chr3L.csv\n')
bash_output.write('rm details_fst_output_chr2R.out\n')
bash_output.write('rm details_fst_output_chr2L.out\n')
bash_output.write('rm details_fst_output_chr3R.out\n')
bash_output.write('rm details_fst_output_chr3L.out\n')
bash_output.write('\n')
bash_output.write('cd ms_output\n')
bash_output.write('cat $z >> recrate'+r_rate+'_'+sys.argv[1]+'.csv\n')
bash_output.write('rm $z\n')
bash_output.write('done\n')



bash_output.write('mv recrate'+r_rate+'_'+sys.argv[1]+'.csv ../\n')

bash_output.close()

if sys.argv[1] == '0':
    output_submit = open('submit_recomb.sub','w')
    output_submit.write('universe = vanilla\n')
    output_submit.write('\n')
    output_submit.write('log = recomb_cmd_$(Process).log\n')
    output_submit.write('error = recomb_cmd_$(Process).err\n')
    output_submit.write('executable = ms_rec_$(Process).sh\n')
    output_submit.write('arguments = $(Process)\n')
    output_submit.write('output = rec_ms_$(Process).out\n')
    output_submit.write('\n')
    output_submit.write('should_transfer_files = YES\n')
    output_submit.write('when_to_transfer_output = ON_EXIT\n')
    output_submit.write('transfer_input_files = http://proxy.chtc.wisc.edu/SQUID/chtc/python37.tar.gz, ../../ms_fst_chtc_minAlleleFreq_targetN.py, ../../ms_chtc.tar.gz\n')
    output_submit.write('\n')
    output_submit.write('request_cpus = 1\n')
    output_submit.write('request_memory = 3GB\n')
    output_submit.write('request_disk = 2MB\n')
    output_submit.write('\n')
    output_submit.write('queue '+sys.argv[2]+'\n')
    output_submit.close()
    
