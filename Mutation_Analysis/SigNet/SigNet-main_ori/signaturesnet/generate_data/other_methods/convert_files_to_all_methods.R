
########################################################################################################################################################################
# 3. Adapting input dataset to catalog dataset: (if we want to run all the other methods)

experiment_id <- "exp_all"

catalog <- read.csv(paste("../../../data/real_data/ground.truth.syn.catalog.csv", sep = ""))
test_input <- read.csv(paste("../../../data/", experiment_id, "/test_input.csv", sep = ""), header = FALSE, row.names = NULL, sep = ",") 
test_label <- read.csv(paste("../../../data/", experiment_id, "/test_label.csv", sep = ""), header = FALSE, row.names = NULL, sep = ",") 
num_muts <- unique(test_label$V73)

k <- 0
for(num_mut in num_muts){
  test_input_muts <- test_input[test_label$V73 == num_mut,]
  test_input_muts <- t(test_input_muts)
  test_input_muts <- cbind(catalog[,c(1,2)], test_input_muts)
  colnames(test_input_muts) <- c(colnames(catalog)[c(1,2)], paste("sample_", (2780*k+1):(2780*(k+1)), sep = ""))
  k <- k+1
  if(dir.exists(paste("../../../data/", experiment_id, "/other_methods/", sep=""))==FALSE){
    dir.create(paste("../../../data/", experiment_id, "/other_methods/", sep=""))
  }
  write.table(test_input_muts, paste("../../../data/", experiment_id, "/other_methods/", "test_ground.truth.syn.catalog_", k,".csv", sep = ""), row.names = FALSE, sep = ",") 
}
