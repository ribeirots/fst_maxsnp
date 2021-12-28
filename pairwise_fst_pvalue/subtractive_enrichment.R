#!/usr/bin/env Rscript
library(ggplot2)
library(grr)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

# Modify this block for other population pairs, arms, and simulated data
directory = "EF_ZI"
arm = "2L"
simdata = "2R"
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
  return(bins_0to100)
}

# Checks number of windows in the bins "0 to 0.05" versus "0.05 to 0.1" and determines whether the first bin is still enriched
sim_re_bin <- function(empdb_list){
  empdb <- bin100(data.frame(pv=empdb_list))

  bin0_0.05 <- sum(empdb$count[1:5])
  bin0.05_0.1 <- sum(empdb$count[6:10])
  
  if(bin0_0.05 <= bin0.05_0.1){
    return(FALSE) } else {
      return(TRUE)
    }
}


# make file names
filename = paste('emp_pvalue_',directory,'_',arm,'_sim',simdata,'.txt',sep='') # _rec05plus
print(filename)
# read files
file_emp_geral = read.csv(filename, sep='\t') # Start, End, WindowFST, MaxSNPFST, MaxSNPPosition, RecRate, WindowPvalue, SNPPvalue
file_emp <- file_emp_geral

regions_removed = 0
rm_rgions <- c()
rm_rgs_only <- c()


# FST_MaxSNP
while(sim_re_bin(file_emp$SNPPvalue)){  

top_value_match <- grr::matches(max(file_emp$MaxSNPFST), file_emp$MaxSNPFST, all.y=F)

if(length(top_value_match$y) == 1){
  to_rm_i <- top_value_match$y
  print(file_emp$SNPPvalue[to_rm_i])
} else {
  to_rm_i <- sample(top_value_match$y, 1)
  print(file_emp$SNPPvalue[to_rm_i])
}

j = 0
w = 1
while(j < 5){ # 5 = max number of "non-outlier windows" in between "outlier windows"
  if(to_rm_i-w > 1){
    
    if(file_emp$End[to_rm_i-w] == (file_emp$Start[to_rm_i-w+1]-1) ){
      if(file_emp$SNPPvalue[to_rm_i-w] >= 0.1){ # 0.1 = max pvalue in between "outlier windows"
        j = j + 1
      } else {
        j = 0
      }
      w = w + 1
    } else {
      break
    }
    } else {
      break
    }
}
w = w - 1

# Checks the other side
j = 0
z = 1
while(j < 5){
  if(to_rm_i+z <= length(file_emp$SNPPvalue)){ 
    if(file_emp$Start[to_rm_i+z] == (file_emp$End[to_rm_i+z-1]+1) ){
      if(file_emp$SNPPvalue[to_rm_i+z] >= 0.1){ 
        j = j + 1
      } else {
        j = 0
      }
      z = z + 1
    } else {
      break
    }
  } else {
    break
  }
}

z = z - 1

if(to_rm_i == 1){
  start_window <- to_rm_i
} else {
  start_window <- to_rm_i-w
}

if(to_rm_i==length(file_emp$SNPPvalue)){ 
  end_window <- to_rm_i
} else {
  end_window <- to_rm_i+z
}


start_reg <- file_emp$Start[start_window]
end_reg <- file_emp$End[end_window]

rmving_emp <- file_emp[seq(start_window,end_window),]
rm_rgions <- rbind(rm_rgions, rmving_emp)
rm_rgs_only <- rbind(rm_rgs_only, c(start_reg, end_reg))
file_emp <- file_emp[-seq(start_window,end_window),]
regions_removed = regions_removed + 1
print(regions_removed)

}
snp_rmved <- rm_rgions[,1:2]
write.csv(rm_rgions, paste('SNP_enrichment_windows_removed_',arm,'_sim',simdata,'_pvalue0.1.csv',sep=''), row.names=F)
write.csv(rm_rgs_only, paste('SNP_enrichment_regions_removed_',arm,'_sim',simdata,'_pvalue0.1.csv',sep=''), row.names=F)

file_emp_win <- file_emp_geral
regions_removed = 0
rm_rgions <- c()
rm_rgs_only <- c()

# FST_Window
while(sim_re_bin(file_emp_win$WindowPvalue)){  
  
  top_value_match <- grr::matches(max(file_emp_win$WindowFST), file_emp_win$WindowFST, all.y=F) 
  
  if(length(top_value_match$y) == 1){
    to_rm_i <- top_value_match$y
  } else {
    to_rm_i <- sample(top_value_match$y, 1)
  }
  
  j = 0
  w = 1
  while(j < 5){
    if(to_rm_i-w > 1){
      
      if(file_emp_win$End[to_rm_i-w] == (file_emp_win$Start[to_rm_i-w+1]-1) ){
        if(file_emp_win$WindowPvalue[to_rm_i-w] >= 0.1){
          j = j + 1
        } else {
          j = 0
        }
        w = w + 1
      } else {
        break
      }
      
    } else {
      break
    }
  }
  
  w = w -1
  
  j = 0
  z = 1
  while(j < 5){
    if(to_rm_i+z <= length(file_emp_win$WindowPvalue)){ 
      if(file_emp_win$Start[to_rm_i+z] == (file_emp_win$End[to_rm_i+z-1]+1) ){
        if(file_emp_win$WindowPvalue[to_rm_i+z] >= 0.1){ 
          j = j + 1
        } else {
          j = 0
        }
        z = z + 1
      } else {
        break
      }
    } else {
      break
    }
  }
  
  z = z - 1
  
  if(to_rm_i == 1){
    start_window <- to_rm_i
  } else {
    start_window <- to_rm_i-w
  }
  
  if(to_rm_i==length(file_emp_win$WindowPvalue)){ 
    end_window <- to_rm_i
  } else {
    end_window <- to_rm_i+z
  }
  
  
  start_reg <- file_emp_win$Start[start_window]
  end_reg <- file_emp_win$End[end_window]
  
  rmving_emp <- file_emp_win[seq(start_window,end_window),]
  rm_rgions <- rbind(rm_rgions, rmving_emp)
  rm_rgs_only <- rbind(rm_rgs_only, c(start_reg, end_reg))
  file_emp_win <- file_emp_win[-seq(start_window,end_window),]
  regions_removed = regions_removed + 1
  print(regions_removed)
  
}
write.csv(rm_rgions, paste('window_enrichment_windows_removed_',arm,'_sim',simdata,'_pvalue0.1.csv',sep=''), row.names=F)
write.csv(rm_rgs_only, paste('window_enrichment_regions_removed_',arm,'_sim',simdata,'_pvalue0.1.csv',sep=''), row.names=F)

win_rmved <- rm_rgions[,1:2]
comparison <- intersect(snp_rmved,win_rmved)
length(comparison[,1])

length(setdiff(snp_rmved,win_rmved)[,1]) # unique windows removed on SNP
length(setdiff(win_rmved,snp_rmved)[,1]) # unique windows removed on Win




