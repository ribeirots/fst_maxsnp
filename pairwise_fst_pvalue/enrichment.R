# Script to generate heatmap plots for enrichment'
# tribeiro@wisc.edu
rm(list=ls())

library(stringr)
library(scales)

setwd('~/fst_enrichment/pvalue_dist/EF_ZI/')
# Function to bin the pvalues windows by pvalue in 100 bins

bin100 <- function(datablist,count_label, target_stat_id=1, enrich=F){
  bins_0to100 <- data.frame(seq(0.005,0.995,0.01))
  bins_0to100$count <- 0
  bin_index <- 1
  for(b in seq(0,0.9900001,0.0100000001)){
    if(b < 0.99){
      bins_0to100$count[bin_index] <- length(datablist[,target_stat_id][datablist[,target_stat_id]>=b & datablist[,target_stat_id] < (b+0.01)])
    } else {
      bins_0to100$count[bin_index] <- length(datablist[,target_stat_id][datablist[,target_stat_id]>=b & datablist[,target_stat_id] < (b+0.01)])
    }
    bin_index <- bin_index + 1
  }
  names(bins_0to100)[1] <- names(datablist[target_stat_id])
  bins_0to100$scale <- count_label
  
  total_obs <- sum(bins_0to100$count)
  
  bins_0to100$prob <- 0
  bin_index = 1
  
  if(enrich==T){
    for(c in bins_0to100$count){
      bins_0to100$prob[bin_index] <- c/(total_obs/100) # enrichment given 100 bins
      bin_index = bin_index + 1
    }
  } else {
    for(c in bins_0to100$count){
      bins_0to100$prob[bin_index] <- c/total_obs # %
      bin_index = bin_index + 1
    }  
  }
  
  return(bins_0to100)
}

# db_empirical, db_sim with count # e.g. sim_re_bin(emp_snp_pv, sim_snp_pv, "SNP")
sim_re_bin <- function(empdb, simdb, class_label, stat_label= "enrichment"){
  bins_1to100 <- data.frame(seq(0.005,0.995,0.01))
  names(bins_1to100)[1] <- 'pv'
  
  bins_1to100$enrichment <- 0
  
  empcount_total <- sum(empdb$count)
  simcount_total <- sum(simdb$count)
  
  for(b in seq(1,100)){
    expected_emp <- (empcount_total*simdb$count[b])/simcount_total
    if(expected_emp == 0 & empdb$count[b] == 0){
      bins_1to100$enrichment[b] <- 1
    } else { 
      bins_1to100$enrichment[b] <- empdb$count[b]/expected_emp
    }
    
  }
  
  bins_1to100$scale <- class_label
  
  return(bins_1to100)
}

# read files
file_emp = read.csv("emp_pvalue_EF_ZI_X_simX.txt", sep='\t') # Start, End, WindowFST, MaxSNPFST, MaxSNPPosition, RecRate, WindowPvalue, SNPPvalue
file_sim = read.csv("sim_emp_pvalue_EF_ZI_X_simX.txt", sep='\t') # WindowFST, MaxSNPFST, WindowPvalue, SNPPvalue, RecRate


# Read files
emp_snp_dist <- data.frame(fst = file_emp$MaxSNPFST)
sim_snp_dist <- data.frame(fst = file_sim$MaxSNPFST)

emp_win_dist <- data.frame(fst = file_emp$WindowFST)
sim_win_dist <- data.frame(fst = file_sim$WindowFST)


# Enrichment #

# Get p-value data
emp_wind_p_dist <- data.frame(pv = file_emp$WindowPvalue)

emp_snp_p_dist <- data.frame(pv = file_emp$SNPPvalue)
sim_snp_p_dist <- data.frame(pv = file_sim$SNPPvalue)


# Make enrichment bins
emp_win_pv <- bin100(emp_wind_p_dist, 'Window', enrich=T)
sim_snp_pv <- bin100(sim_snp_p_dist, 'SNP', enrich=T)
emp_snp_pv <- bin100(emp_snp_p_dist, 'SNP', enrich=T)

# Calculate SNP enrichment based on Sim Expectation
snp_e <- sim_re_bin(emp_snp_pv, sim_snp_pv, "SNP")
# Re-format window enrichment to match SNP
win_pv <- emp_win_pv
names(win_pv)[4] <- 'enrichment'
win_e <- win_pv[,c(1,4)]
win_e$scale <- 'Window'
comb_pvalue <- rbind(win_e,snp_e)

