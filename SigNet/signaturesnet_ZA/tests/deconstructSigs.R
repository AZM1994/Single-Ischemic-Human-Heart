
library(deconstructSigs)
library(readxl)
library(parallel)

cosmic.sigs <- read.table('data/COSMIC_signatures_v3.xlsx', sep="\t")
cosmic.sigs_3 <- read_excel('data/COSMIC_Mutational_Signatures_v3.1.xlsx')
cosmic.sigs_3 <-  t(cosmic.sigs_3)
subs <- cosmic.sigs_3[1,]
cont <- cosmic.sigs_3[2,]
cont1 <- substr(cont, 1, 1)
cont3 <- substr(cont, 3,3)
new_header <- paste(cont1,"[",subs,"]", cont3,sep ="")
colnames(cosmic.sigs) <- new_header


data <- read.csv("data/test_input_w01.csv", header = FALSE)
labels <- read.csv("data/test_label_w01.csv", header = FALSE)
num_mutations <- labels[,73]
data <- data * num_mutations

new_colnames <- colnames(cosmic.sigs) 
order1 <- c(1,2,3,4,25,26,27,28,49,50,51,52,73,74,75,76)
order2 <- order1 + 4
order3 <- order2 + 4
order5 <- c(13,14,15,16,37,38,39,40,61,62,63,64,85,86,87,88)
order6 <- order5 + 4
order7 <- order6 + 4
order_cols <- c(order1,order2,order3, order5, order6, order7)
new_colnames <- new_colnames[order_cols]
rownames(data) <- 1:nrow(data)
colnames(data) <- new_colnames
sigs <- c()
for(i in 1:1000){
  print(i)
  sigs_res <- whichSignatures(data, sample.id=as.character(i), signatures.ref = cosmic.sigs, signature.cutoff = 0,contexts.needed = TRUE, tri.counts.method ="default")
  sigs <- rbind(sigs,sigs_res$weights[1,])
}

whichSigs <- function(x){
  library(deconstructSigs)
  sigs_res <- whichSignatures(data, sample.id=as.character(x), signatures.ref = cosmic.sigs, signature.cutoff = 0,contexts.needed = TRUE, tri.counts.method ="default")
  return(sigs_res$weights[1,])
} 

no_cores <- 8
clust <- makeCluster(no_cores)
clusterExport(clust, c("data", "cosmic.sigs"))
sigs <- parSapply(clust, 1:nrow(data), whichSigs)
stopCluster(clust)

write.csv(t(sigs), file="data/deconstructSigs_test_w01.csv", row.names = FALSE)#, col.names = FALSE)

