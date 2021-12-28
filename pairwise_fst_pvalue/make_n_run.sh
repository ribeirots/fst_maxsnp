#!/bin/bash\n

# The following files are needed on chtc:
#### msms_recomb.py
#### recrate*.txt files with rec rates and window lengths
#### msms3.2rc-b163.jar


for i in {0..500}
do
# $i = id. 
# 2nd argv = total number of jobs
# 3rd argv = simulations per job
# 4th argv = recombination rate code/file. Current codes available are:
#### 0
#### 0.5
#### 1
#### 1.5
#### 2
#### 3
#### inf
python3 msms_recomb.py $i 500 100 2
done

condor_submit submit_recomb.sub
