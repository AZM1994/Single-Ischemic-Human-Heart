# Convert output files from each method to SigNet format:

experiment_id <- "exp_all"
list_of_methods <- c("decompTumor2Sig", "deconstructSigs", "mSigAct", "MutationalPatterns", "mutSignatures", "SignatureEstimationQP", "YAPSA" )

for(i in 1:length(list_of_methods)){
  final_guess <- c()
  for(j in 1:10){
    test_guess <- read.csv(paste("../../../data/", experiment_id, "/other_methods/", list_of_methods[i], "_test_", j,"/inferred.exposures.csv", sep=""), row.names = 1)
    final_guess <- rbind(final_guess, t(test_guess))
  } 
  write.table(test_guess, paste("../../../data/", experiment_id, "/other_methods/all_results/", list_of_methods[i], "_guess.csv", sep = ""), col.names = FALSE, row.names = FALSE, sep = ",")
}
