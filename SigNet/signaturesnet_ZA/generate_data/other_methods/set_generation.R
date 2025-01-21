# Generate test set:

library(SynSigGen)
library(ICAMS)
library(readxl)
library(progress)

generate_realistic_data <- function(cosmic_version, dataset, large_low, experiment_id){
  name_dataset <- paste(dataset, '_realistic_', large_low, sep = '')
  dir_name <- "PCAWG_normalized/"
  #1. Generate synthetic dataset using PCAWG sigProfiler results:
  #SBS <- read.csv("../../data/real_data/signatures_used_PCAWG_v3.csv")
  SBS <- read_xlsx("../../data/data.xlsx")
  SBS <- SBS[,2:ncol(SBS)]
  num_sigs <- ncol(SBS)
  
  E <- read.csv("../../data/real_data/sigprofiler_normalized_PCAWG.csv", row.names = 1)  
  E <- t(E)
  
  S <- ReadCatalog("../../data/real_data/signatures_used_PCAWG_v3.csv", catalog.type = "density.signature")
 
  P <- GetSynSigParamsFromExposures(E)
 
  if(dataset == "train"){
    if(large_low == "large"){
      number_of_samples = 250000
    }else if(large_low == "low"){
      number_of_samples = 350000
    }
  } else if (dataset == "val"){
    number_of_samples = 5000
 } else if (dataset == "test"){
    number_of_samples = 19000
  }
  synthetic.exposures <- GenerateSyntheticExposures(P, num.samples = number_of_samples)
  synthetic.spectra <- CreateAndWriteCatalog(S, synthetic.exposures, my.dir = paste(dir_name, dataset, "_", large_low, sep=""))     
 
 
  # 2. Adapting datasets to signatures-net:
 
  # LABEL:
  if(dataset == "test"){
    num_muts <- c(25, 50, 100, 250, 500, 1000, 5000, 10000, 50000, 100000)
    num_muts <- rep(num_muts, each=1900)
  }
 
  if (dataset == "train"){
    if(large_low == "low"){
      range_muts <- c(15, 50, 100, 250, 500, 1000, 5000, 10000)
      num_samples <- c(50000,50000,50000,50000,50000,50000,50000)
    }else if( large_low == "large"){
      range_muts <- c(1e3, 5e3, 1e4, 5e4, 1e5, 5e5)
      num_samples <- c(50000,50000,50000,50000,50000)
    }
    
    num_muts <- c()
    for(i in 1:length(num_samples)){
      num_muts <- c(num_muts, sample(range_muts[i]:range_muts[i+1], num_samples[i], replace = TRUE))
    }

  }
 
  if (dataset == "val"){
    if(large_low == "low"){
      range_muts <- c(15, 50, 100, 250, 500, 1000)
    }else if( large_low == "large"){
      range_muts <- c(1e3, 5e3, 1e4, 5e5, 1e5, 5e5)
    }
    num_samples <- c(1000, 1000, 1000, 1000, 1000)

    num_muts <- c()
    pb <- progress_bar$new(total = length(num_samples))
    for(i in 1:length(num_samples)){
      num_muts <- c(num_muts, sample(range_muts[i]:range_muts[i+1], num_samples[i], replace = TRUE))
      pb$tick()
    }
  }
 
  test_label <- read.csv(paste(dir_name, dataset, "_", large_low, "/ground.truth.syn.exposures.csv", sep = ""), row.names = 1)
  test_label <- t(test_label)
  test_label <- test_label / rowSums(test_label)
  test_label <- cbind(test_label, num_muts)
 
  test_label_final <- matrix(0, ncol = num_sigs, nrow = nrow(test_label)) 
  colnames(test_label_final) <- colnames(SBS)
  test_label_final <- as.data.frame(test_label_final)
  for(i in 1:ncol(test_label)){
    test_label_final[[colnames(test_label)[i]]] <- test_label[,i]
  }
  write.table(test_label_final, paste(dir_name, name_dataset, "_label.csv", sep = ""), col.names = FALSE, row.names = FALSE, sep = ",") 

  # INPUT:
  test_label <- read.csv(paste(dir_name, name_dataset, "_label.csv", sep = ""), row.names = NULL, header = FALSE) 
  num_muts <- test_label[,num_sigs+1]
  test_label <- test_label[,1:num_sigs]
  test_input_muts <- c()
  pb <- progress_bar$new(total = nrow(test_label))
  for(i in 1:nrow(test_label)){
    prob <- rep(0, 96)
    for(j in 1:ncol(test_label)){
      prob <- prob + test_label[i,j]*SBS[,j]
    }
    if(num_muts[i]>=100000){
      test_input_muts <- rbind(test_input_muts,t(prob*num_muts[i]))
    } 
    else{
      test_input_muts <- rbind(test_input_muts,t(rmultinom(1, size=num_muts[i], t(prob))))
    }
    pb$tick()
  }
  norm_test_input_muts <- test_input_muts/num_muts
  write.table(norm_test_input_muts, paste(dir_name, name_dataset,  "_input.csv", sep = ""), col.names = FALSE, row.names = FALSE, sep = ",") 
  write.table(test_input_muts, paste(dir_name, name_dataset,  "_input_not_normalized.csv", sep = ""), col.names = FALSE, row.names = FALSE, sep = ",") 
}

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

dataset = args[1]
large_low = args[2]

print(paste("Generating ", dataset, large_low, ":"))
generate_realistic_data(cosmic_version, dataset, large_low, experiment_id)
print("DONE!")