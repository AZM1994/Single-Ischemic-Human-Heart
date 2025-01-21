library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DESeq2)
library(enrichplot)
library(ggplot2)
library(stringr)
library(DOSE)

setwd("/Users/zhemingan/Downloads/annovar/heart_PTA_Cases/GOseq_Zheming/")
wd_dir <- getwd()
result_dir <- paste0(wd_dir, "/enrichGO")
dir.create(result_dir)
# dir.create(paste0(result_dir, "/enrichGO"))

### read gene list with metadata
genomic_context_normal <- read.csv("heart_PTA_Cases_Normal.SNV_gene_list.csv", header = TRUE)
heart_PTA_Cases_Normal_vcf <- read.table("heart_PTA_Cases.all_normal_ssnv.vcf", sep = "\t")
genomic_context_normal$Condition <- "Normal"
genomic_context_normal$Cell_ID <- heart_PTA_Cases_Normal_vcf$V8
genomic_context_normal$Case_ID <- str_extract(genomic_context_normal$Cell_ID, "[^_]+")
genomic_context_normal <- genomic_context_normal[c("Chr","Start","End","Ref","Alt","Cell_ID","Case_ID","Condition","Func.refGene","Gene.refGene")]
colnames(genomic_context_normal) <- c("Chr","Start","End","Ref","Alt","Cell_ID","Case_ID","Condition","Type","Gene_symbol")

genomic_context_disease <- read.csv("heart_PTA_Cases_Disease.SNV_gene_list.csv", header = TRUE)
heart_PTA_Cases_Disease_vcf <- read.table("heart_PTA_Cases.all_disease_ssnv.vcf", sep = "\t")
genomic_context_disease$Condition <- "Disease"
genomic_context_disease$Cell_ID <- heart_PTA_Cases_Disease_vcf$V8
genomic_context_disease$Case_ID <- str_extract(genomic_context_disease$Cell_ID, "[^_]+")
genomic_context_disease <- genomic_context_disease[c("Chr","Start","End","Ref","Alt","Cell_ID","Case_ID","Condition","Func.refGene","Gene.refGene")]
colnames(genomic_context_disease) <- c("Chr","Start","End","Ref","Alt","Cell_ID","Case_ID","Condition","Type","Gene_symbol")
genomic_context_all <- rbind(genomic_context_normal, genomic_context_disease)
genomic_context_all <- genomic_context_all[genomic_context_all$Type == "exonic", ]

Gene_list_normal <- genomic_context_normal$Gene_symbol %>% na.omit()
Gene_list_disease <- genomic_context_disease$Gene_symbol %>% na.omit()
Gene_list_all <- genomic_context_all$Gene_symbol %>% na.omit()

#################################################################################
############################## Gene Set Enrichment ##############################
#################################################################################
# Gene_list_all <- genomic_context_all$Start
# names(Gene_list_all) <- as.character(genomic_context_all$Gene_symbol)
# Gene_list_all <- sort(Gene_list_all, decreasing = TRUE)
# 
# gse_results <- gseGO(geneList = Gene_list_all, OrgDb = "org.Hs.eg.db", ont ="ALL", keyType = "SYMBOL",
#                      nPerm = 100, pvalueCutoff = 0.05, pAdjustMethod = "fdr", verbose = TRUE)
# dotplot(gse_results, showCategory=10)
# dotplot(gse_results, showCategory=10, split=".sign") + facet_grid(.~.sign)
# barplot(gse_results, split = "ONTOLOGY")+facet_grid(ONTOLOGY~., scale = "free")+ggtitle("Barplot for GO_ORA_PC")
# emapplot(gse_results, showCategory = 10)
# gseaplot(gse_results, geneSetID = 1, title = gse_results$Description[1])

#################################################################################
######################## GO over-representation analysis ########################
#################################################################################
GO_results <- enrichGO(gene = Gene_list_all, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                       pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")

dotplot(GO_results, showCategory=15)
barplot(GO_results, split = "ONTOLOGY")+facet_grid(ONTOLOGY~., scale = "free")


pdf(paste0(result_dir, "/GO_results_normal.pdf"), width = 10, height = 8)
fit <- plot(barplot(GO_results, showCategory = 20))
dev.off()

pathways <- GO_results_normal$Description
fdr_values <- GO_results_normal$qvalue

# Plot dotplot of pathways with FDR values
dotplot(pathways, fdr_values, xlab = "-log10(FDR)", pch = 19)
