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

genomic_context_normal_exonic <- genomic_context_normal[genomic_context_normal$Type == "exonic", ]
genomic_context_disease_exonic <- genomic_context_disease[genomic_context_disease$Type == "exonic", ]
genomic_context_all_exonic <- genomic_context_all[genomic_context_all$Type == "exonic", ]

Gene_list_normal <- genomic_context_normal$Gene_symbol %>% na.omit()
Gene_list_disease <- genomic_context_disease$Gene_symbol %>% na.omit()
Gene_list_all <- genomic_context_all$Gene_symbol %>% na.omit()

Gene_list_normal_exonic <- genomic_context_normal_exonic$Gene_symbol %>% na.omit()
Gene_list_disease_exonic <- genomic_context_disease_exonic$Gene_symbol %>% na.omit()
Gene_list_all_exonic <- genomic_context_all_exonic$Gene_symbol %>% na.omit()

#################################################################################
######################## GO over-representation analysis ########################
#################################################################################
GO_results_normal <- enrichGO(gene = Gene_list_normal, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                           pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")
GO_results_disease <- enrichGO(gene = Gene_list_disease, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                              pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")
GO_results_all <- enrichGO(gene = Gene_list_all, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                               pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")

GO_results_normal_exonic <- enrichGO(gene = Gene_list_normal_exonic, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                              pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")
GO_results_disease_exonic <- enrichGO(gene = Gene_list_disease_exonic, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                               pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")
GO_results_all_exonic <- enrichGO(gene = Gene_list_all_exonic, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                           pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")

pdf(paste0(result_dir, "/GO_results_normal.pdf"), width = 12, height = 15)
dotplot(GO_results_normal, showCategory=15, title = "GO_results_normal dotplot") 
barplot(GO_results_normal, showCategory=15, title = "GO_results_normal barplot") 
dotplot(GO_results_normal, split = "ONTOLOGY", title = "GO_results_normal dotplot split by ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale = "free")
barplot(GO_results_normal, split = "ONTOLOGY", title = "GO_results_normal barplot split by ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale = "free") 
pairwise_term_GO_results_normal <- pairwise_termsim(GO_results_normal) 
emapplot(pairwise_term_GO_results_normal, showCategory = 15) + ggtitle("Enrichment map for GO_results_normal")
dev.off()

# pdf(paste0(result_dir, "/GO_results_normal_exonic.pdf"), width = 12, height = 15)
# dotplot(GO_results_normal_exonic, showCategory=15, title = "GO_results_normal_exonic dotplot") 
# barplot(GO_results_normal_exonic, showCategory=15, title = "GO_results_normal_exonic barplot") 
# dotplot(GO_results_normal_exonic, split = "ONTOLOGY", title = "GO_results_normal_exonic dotplot split by ONTOLOGY") +
#   facet_grid(ONTOLOGY~., scale = "free")
# barplot(GO_results_normal_exonic, split = "ONTOLOGY", title = "GO_results_normal_exonic barplot split by ONTOLOGY") + 
#   facet_grid(ONTOLOGY~., scale = "free") 
# pairwise_term_GO_results_normal_exonic <- pairwise_termsim(GO_results_normal_exonic) 
# emapplot(pairwise_term_GO_results_normal_exonic, showCategory = 15) + ggtitle("Enrichment map for GO_results_normal_exonic")
# dev.off()

pdf(paste0(result_dir, "/GO_results_disease.pdf"), width = 12, height = 15)
dotplot(GO_results_disease, showCategory=15, title = "GO_results_disease dotplot") 
barplot(GO_results_disease, showCategory=15, title = "GO_results_disease barplot") 
dotplot(GO_results_disease, split = "ONTOLOGY", title = "GO_results_disease dotplot split by ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale = "free")
barplot(GO_results_disease, split = "ONTOLOGY", title = "GO_results_disease barplot split by ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale = "free")
pairwise_term_GO_results_disease <- pairwise_termsim(GO_results_disease)
emapplot(pairwise_term_GO_results_disease, showCategory = 15) + ggtitle("Enrichment map for GO_results_disease")
dev.off()

# pdf(paste0(result_dir, "/GO_results_disease_exonic.pdf"), width = 12, height = 15)
# dotplot(GO_results_disease_exonic, showCategory=15, title = "GO_results_disease_exonic dotplot") 
# barplot(GO_results_disease_exonic, showCategory=15, title = "GO_results_disease_exonic barplot") 
# dotplot(GO_results_disease_exonic, split = "ONTOLOGY", title = "GO_results_disease_exonic dotplot split by ONTOLOGY") +
#   facet_grid(ONTOLOGY~., scale = "free")
# barplot(GO_results_disease_exonic, split = "ONTOLOGY", title = "GO_results_disease_exonic barplot split by ONTOLOGY") + 
#   facet_grid(ONTOLOGY~., scale = "free")
# pairwise_term_GO_results_disease_exonic <- pairwise_termsim(GO_results_disease_exonic)
# emapplot(pairwise_term_GO_results_disease_exonic, showCategory = 15) + ggtitle("Enrichment map for GO_results_disease_exonic")
# dev.off()

pdf(paste0(result_dir, "/GO_results_all.pdf"), width = 12, height = 15)
dotplot(GO_results_all, showCategory=15, title = "GO_results_all dotplot") 
barplot(GO_results_all, showCategory=15, title = "GO_results_all barplot") 
dotplot(GO_results_all, split = "ONTOLOGY", title = "GO_results_all dotplot split by ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale = "free")
barplot(GO_results_all, split = "ONTOLOGY", title = "GO_results_all barplot split by ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale = "free")
pairwise_term_GO_results_all <- pairwise_termsim(GO_results_all)
emapplot(pairwise_term_GO_results_all, showCategory = 15) + ggtitle("Enrichment map for GO_results_all")
dev.off()

# pdf(paste0(result_dir, "/GO_results_all_exonic.pdf"), width = 12, height = 15)
# dotplot(GO_results_all_exonic, showCategory=15, title = "GO_results_all_exonic dotplot") 
# barplot(GO_results_all_exonic, showCategory=15, title = "GO_results_all_exonic barplot") 
# dotplot(GO_results_all_exonic, split = "ONTOLOGY", title = "GO_results_all_exonic dotplot split by ONTOLOGY") +
#   facet_grid(ONTOLOGY~., scale = "free")
# barplot(GO_results_all_exonic, split = "ONTOLOGY", title = "GO_results_all_exonic barplot split by ONTOLOGY") + 
#   facet_grid(ONTOLOGY~., scale = "free")
# pairwise_term_GO_results_all_exonic <- pairwise_termsim(GO_results_all_exonic)
# emapplot(pairwise_term_GO_results_all_exonic, showCategory = 15) + ggtitle("Enrichment map for GO_results_all_exonic")
# dev.off()
