
# source("https://bioconductor.org/biocLite.R")
# biocLite("BiocInstaller")

# install.packages("BiocManager")
# BiocManager::install(c("Biostrings","BSgenome","GenomeInfoDb","GenomicRanges"))

# install.packages("remotes")
# remotes::install_github(repo = "steverozen/SynSigGen", ref = "1.0.6-branch")

# install.packages("devtools")
# devtools::install_github("WuyangFF95/SynSigRun")
library(SynSigRun)
 
cosmic_version = "v3"
if(cosmic_version == "v3"){
    path_to_sigs = "../../../data/exp_generator/test_generator/other_methods/data.csv"
    
}else if(cosmic_version == "v2"){
  path_to_sigs = "../../../data/exp_v2/test/other_methods/data_v2.csv"
  
}else{
  print("Error: not implemented for this version of COSMIC")
  break
}

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
method <- args[1]
dataset <- args[2]
replicate <- args[3]

# method <- "all"
# dataset <- "test"
experiment_id <- "exp_all"

# Reffiting only methods:
if(method == "decompTumor2Sig" || method == "all"){
  print("Running decompTumor2Sig")
  start.time <- Sys.time()
  RundecompTumor2SigAttributeOnly(
    input.catalog = paste("../../../data/", experiment_id, "/other_methods/test_ground.truth.syn.catalog_",k,".csv", sep = ""),
    gt.sigs.file = path_to_sigs,
    paste("../../../data/", experiment_id, "/other_methods/decompTumor2Sig_test", sep=""),
    seedNumber = 1,
    test.only = FALSE,
    overwrite = FALSE
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

if(method == "MutationalPatterns" || method == "all"){
  print("Running MutationalPatterns")
  start.time <- Sys.time() 
  RunMutationalPatternsAttributeOnly(
    input.catalog = paste("../../../data/", experiment_id, "/other_methods/test_ground.truth.syn.catalog_",k,".csv", sep = ""),
    gt.sigs.file = path_to_sigs,
    paste("../../../data/", experiment_id, "/other_methods/MutationalPatterns_test", sep=""),
    seedNumber = 1,
    test.only = FALSE,
    overwrite = FALSE
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

if(method == "mutSignatures" || method == "all"){
  print("Running mutSignatures")
  start.time <- Sys.time()
  RunmutSignaturesAttributeOnly(
    input.catalog = paste("../../../data/", experiment_id, "/other_methods/test_ground.truth.syn.catalog_",k,".csv", sep = ""),
    gt.sigs.file = path_to_sigs,
    paste("../../../data/", experiment_id, "/other_methods/mutSignatures_test", sep=""),
    seedNumber = 1,
    test.only = FALSE,
    overwrite = FALSE
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

if(method == "SignatureEstimationQP" || method == "all"){
  print("Running SignatureEstimationQP")
  start.time <- Sys.time()
  RunSignatureEstimationQPAttributeOnly(
    input.catalog = paste("../../../data/", experiment_id, "/other_methods/test_ground.truth.syn.catalog_",k,".csv", sep = ""),
    gt.sigs.file = path_to_sigs,
    paste("../../../data/", experiment_id, "/other_methods/SignatureEstimationQP_test", sep=""),
    seedNumber = 1,
    test.only = FALSE,
    overwrite = FALSE
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

if(method == "SignatureEstimationSA" || method == "all"){
  print("Running SignatureEstimationSA")
  start.time <- Sys.time()
  RunSignatureEstimationSAAttributeOnly(
    input.catalog = paste("../../../data/", experiment_id, "/other_methods/test_ground.truth.syn.catalog_",k,".csv", sep = ""),
    gt.sigs.file = path_to_sigs,
    paste("../../../data/", experiment_id, "/other_methods/SignatureEstimationSA_test", sep=""),
    seedNumber = 1,
    test.only = FALSE,
    overwrite = FALSE
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

if(method == "YAPSA" || method == "all"){
  print("Running YAPSA")
  start.time <- Sys.time()
  RunYAPSAAttributeOnly(
    input.catalog = paste("../../../data/", experiment_id, "/other_methods/test_ground.truth.syn.catalog_",k,".csv", sep = ""),
    gt.sigs.file = path_to_sigs,
    paste("../../../data/", experiment_id, "/other_methods/YAPSA_test", sep=""),
    seedNumber = 1,
    signature.cutoff = NULL,
    test.only = FALSE,
    overwrite = FALSE
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

if(method == "deconstructSigs" || method == "all"){
  print("Running deconstructSigs")
  start.time <- Sys.time()
  RundeconstructSigsAttributeOnly(
    input.catalog = paste("../../../data/", experiment_id, "/other_methods/test_ground.truth.syn.catalog_",k,".csv", sep = ""),
    gt.sigs.file = path_to_sigs,
    paste("../../../data/", experiment_id, "/other_methods/deconstructSigs_test", sep=""),
    seedNumber = 1,
    test.only = FALSE,
    overwrite = FALSE
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

if(method == "sigfit" || method == "all"){
  print("Running sigfit")
  start.time <- Sys.time()
  RunsigfitAttributeOnly(
    input.catalog = paste("../../../data/", experiment_id, "/other_methods/test_ground.truth.syn.catalog_",k,".csv", sep = ""),
    gt.sigs.file = path_to_sigs,
    paste("../../../data/", experiment_id, "/other_methods/sigfit_test", sep=""),
    seedNumber = 1,
    test.only = FALSE,
    overwrite = FALSE
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

if(method == "mSigAct" || method == "all"){
  print("Running mSigAct")
  start.time <- Sys.time() 
  RunmSigActAttributeOnly(
    input.catalog = paste("../../../data/", experiment_id, "/other_methods/test_ground.truth.syn.catalog_",k,".csv", sep = ""),
    gt.sigs.file = path_to_sigs,
    paste("../../../data/", experiment_id, "/other_methods/mSigAct_test", sep=""),
    seedNumber = 1,
    test.only = FALSE,
    overwrite = FALSE
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}


# Extraction + reffiting methods:

if(method == "signeR" || method == "all"){
  print("Running signeR")
  start.time <- Sys.time()
  RunsigneR(
    input.catalog = paste("../../../data/", experiment_id, "/other_methods/test_ground.truth.syn.catalog_",k,".csv", sep = ""),
    gt.sigs.file = path_to_sigs,
    paste("../../../data/", experiment_id, "/other_methods/signeR_test", sep=""),
    seedNumber = 1,
    test.only = FALSE,
    overwrite = FALSE,
    K.range = c(40,60)
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

if(method == "SignatureAnalyzer" || method == "all"){
  print("Running SignatureAnalyzer")
  start.time <- Sys.time()
  RunSignatureAnalyzerAttribution(
    input.catalog = paste("../../../data/", experiment_id, "/other_methods/test_ground.truth.syn.catalog_",k,".csv", sep = ""),
    gt.sigs.file = path_to_sigs,
    paste("../../../data/", experiment_id, "/other_methods/SignatureAnalyzer_test", sep=""),
    seedNumber = 1,
    test.only = FALSE,
    overwrite = FALSE
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

if(method == "mutSpec" || method == "all"){
  print("Running mutSpec")
  start.time <- Sys.time()
  RunmutSpec(
    input.catalog = paste("../../../data/", experiment_id, "/other_methods/test_ground.truth.syn.catalog_",k,".csv", sep = ""),
    gt.sigs.file = path_to_sigs,
    paste("../../../data/", experiment_id, "/other_methods/mutSpec_test", sep=""),
    seedNumber = 1,
    test.only = FALSE,
    overwrite = FALSE,
    K.range = c(40,60)
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

write.table(difftime(end.time, start.time, units='mins'), paste("../../../data/", experiment_id, "/other_methods/times/", method, dataset,"_", replicate, ".txt",sep=""), sep='\t', col.names = FALSE)