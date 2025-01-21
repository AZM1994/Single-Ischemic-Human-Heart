setwd("/n/data1/bch/genetics/lee/shulin/heart_endo/")
#.libPaths("/n/data1/bch/genetics/lee/shulin/heart_endo/Rpackages/")

library(Seurat)
library(ggplot2)
library(stringr)
library(rtracklayer)

# read in data
count_mat <- ReadMtx(mtx = "data/scATACseq/GSE165837/matrix.mtx.gz", 
                features = "data/scATACseq/GSE165837/features.tsv.gz", 
                cells = "data/scATACseq/GSE165837/barcodes.tsv.gz", 
                feature.column = 1)
metadata <- read.table("data/scATACseq/GSE165837/GSE165837_CARE_ATAC_merged_umap_cluster_labels.tsv", header = T)
data <- CreateSeuratObject(counts = count_mat, assay = "ATAC")
data <- data[,metadata$Sample.Barcode]

# extract endo cell
endo_cell <- data[,metadata$Sample.Barcode[metadata$Cluster == "EC"]]
endo_cell_count <- rowSums(endo_cell@assays$ATAC@counts)

# extract position
feature <- str_split(names(endo_cell_count), pattern = ":", simplify = T)
chr <- feature[,1]
table(chr)
pos <- str_split(feature[,2], pattern = "-", simplify = T)
start <- pos[,1]
end <- pos[,2]

# liftover
chainObject <- import.chain("../ref/liftover_chain/hg38ToHg19.over.chain")
grange_object <- GRanges(seqnames = chr, ranges = IRanges(start = as.integer(start), end = as.integer(end)))

liftover_results <- liftOver(grange_object, chainObject)
bedgraph <- c()
for (i in 1:length(liftover_results)){
  if(length(liftover_results[[i]])){
    feature_count <- 
      c(as.character(liftover_results[[i]][1]@seqnames@values),
        as.character(liftover_results[[i]][1]@ranges@start),
        liftover_results[[i]][1]@ranges@start + liftover_results[[i]][1]@ranges@width,
        endo_cell_count[i])
    bedgraph <- rbind(bedgraph, feature_count)
  }
  if(i %% 1000 == 0){print(i)}
}

write.table(bedgraph, "data/scATACseq/GSE165837/endo_cell_count.bedgraph", quote = F, row.names = F, col.names = F, sep = "\t")
bedgraph <- read.table("data/scATACseq/GSE165837/endo_cell_count.bedgraph", sep = "\t")

bedgraph <- as.data.frame(bedgraph)
colnames(bedgraph) <- c("chr", "start", "end", "count")
bedgraph$chr <- factor(bedgraph$chr, levels = 1:24, labels = levels(liftover_results[[501805]][1]@seqnames@values))
options(scipen = 999)
bedgraph$start <- as.integer(bedgraph$start)
bedgraph$end <- as.integer(bedgraph$end)
bedgraph$end <- bedgraph$end - 1

bedgraph_autosome <- bedgraph[bedgraph$chr %in% paste0("chr", 1:22),]

#overlap_region <- c()
for (i in 2:dim(bedgraph_autosome)[1]){
  if(bedgraph_autosome[i,"start"] < bedgraph_autosome[i-1,"end"]){
    #start <- bedgraph_autosome[i-1,"end"]
    #end <- bedgraph_autosome[i,"start"]
    
    bedgraph_autosome[i,"start"] <- bedgraph_autosome[i-1,"end"] + 1
    #bedgraph_autosome[i-1,"end"] <- end
    #feature_count <- c(bedgraph_autosome[i-1,"end"], bedgraph_autosome[i,"start"], bedgraph_autosome[i,"count"] + bedgraph_autosome[i-1,"count"])
  }
}
write.table(bedgraph_autosome, "data/scATACseq/GSE165837/endo_cell_count_autosome.bedgraph", quote = F, row.names = F, col.names = F, sep = "\t")

