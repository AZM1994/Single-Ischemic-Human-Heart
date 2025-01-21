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
result_dir <- paste0(wd_dir, "/enrichKEGG")
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

kegg_Gene_list_normal <- mapIds(org.Hs.eg.db, Gene_list_normal, 'ENTREZID', 'SYMBOL')
kegg_Gene_list_disease <- mapIds(org.Hs.eg.db, Gene_list_disease, 'ENTREZID', 'SYMBOL')
kegg_Gene_list_all <- mapIds(org.Hs.eg.db, Gene_list_all, 'ENTREZID', 'SYMBOL')

kegg_Gene_list_normal_exonic <- mapIds(org.Hs.eg.db, Gene_list_normal_exonic, 'ENTREZID', 'SYMBOL')
kegg_Gene_list_disease_exonic <- mapIds(org.Hs.eg.db, Gene_list_disease_exonic, 'ENTREZID', 'SYMBOL')
kegg_Gene_list_all_exonic <- mapIds(org.Hs.eg.db, Gene_list_all_exonic, 'ENTREZID', 'SYMBOL')
#################################################################################
######################## GO over-representation analysis ########################
#################################################################################
enrichKEGG_results_normal <- enrichKEGG(gene = kegg_Gene_list_normal, organism = "hsa", pvalueCutoff = 0.05, 
                                        pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "kegg")
enrichKEGG_results_disease <- enrichKEGG(gene = kegg_Gene_list_disease, organism = "hsa", pvalueCutoff = 0.05, 
                                        pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "kegg")
enrichKEGG_results_all <- enrichKEGG(gene = kegg_Gene_list_all, organism = "hsa", pvalueCutoff = 0.05, 
                                         pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "kegg")

enrichKEGG_results_normal_exonic <- enrichKEGG(gene = kegg_Gene_list_normal_exonic, organism = "hsa", pvalueCutoff = 0.05, 
                                        pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "kegg")
enrichKEGG_results_disease_exonic <- enrichKEGG(gene = kegg_Gene_list_disease_exonic, organism = "hsa", pvalueCutoff = 0.05, 
                                         pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "kegg")
enrichKEGG_results_all_exonic <- enrichKEGG(gene = kegg_Gene_list_all_exonic, organism = "hsa", pvalueCutoff = 0.05, 
                                     pAdjustMethod = "fdr", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, keyType = "kegg")

pdf(paste0(result_dir, "/enrichKEGG_results_normal.pdf"), width = 12, height = 15)
dotplot(enrichKEGG_results_normal, showCategory=15, title = "enrichKEGG_results_normal dotplot") 
barplot(enrichKEGG_results_normal, showCategory=15, title = "enrichKEGG_results_normal barplot") 
pairwise_term_enrichKEGG_results_normal <- pairwise_termsim(enrichKEGG_results_normal) 
emapplot(pairwise_term_GO_results_normal, showCategory = 15) + ggtitle("Enrichment map for enrichKEGG_results_normal")
dev.off()

# pdf(paste0(result_dir, "/enrichKEGG_results_normal_exonic.pdf"), width = 12, height = 15)
# dotplot(enrichKEGG_results_normal_exonic, showCategory=15, title = "enrichKEGG_results_normal_exonic dotplot") 
# barplot(enrichKEGG_results_normal_exonic, showCategory=15, title = "enrichKEGG_results_normal_exonic barplot") 
# pairwise_term_GO_results_normal_exonic <- pairwise_termsim(enrichKEGG_results_normal_exonic) 
# emapplot(pairwise_term_GO_results_normal_exonic, showCategory = 15) + ggtitle("Enrichment map for enrichKEGG_results_normal_exonic")
# dev.off()

pdf(paste0(result_dir, "/enrichKEGG_results_disease.pdf"), width = 12, height = 15)
dotplot(enrichKEGG_results_disease, showCategory=15, title = "enrichKEGG_results_disease dotplot") 
barplot(enrichKEGG_results_disease, showCategory=15, title = "enrichKEGG_results_disease barplot") 
pairwise_term_enrichKEGG_results_disease <- pairwise_termsim(enrichKEGG_results_disease)
emapplot(pairwise_term_enrichKEGG_results_disease, showCategory = 15) + ggtitle("Enrichment map for enrichKEGG_results_disease")
dev.off()

# pdf(paste0(result_dir, "/enrichKEGG_results_disease_exonic.pdf"), width = 12, height = 15)
# dotplot(enrichKEGG_results_disease_exonic, showCategory=15, title = "enrichKEGG_results_disease_exonic dotplot") 
# barplot(enrichKEGG_results_disease_exonic, showCategory=15, title = "enrichKEGG_results_disease_exonic barplot") 
# pairwise_term_enrichKEGG_results_disease_exonic <- pairwise_termsim(enrichKEGG_results_disease_exonic)
# emapplot(pairwise_term_enrichKEGG_results_disease_exonic, showCategory = 15) + ggtitle("Enrichment map for enrichKEGG_results_disease_exonic")
# dev.off()

pdf(paste0(result_dir, "/enrichKEGG_results_all.pdf"), width = 12, height = 15)
dotplot(enrichKEGG_results_all, showCategory=15, title = "enrichKEGG_results_all dotplot") 
barplot(enrichKEGG_results_all, showCategory=15, title = "enrichKEGG_results_all barplot") 
pairwise_term_enrichKEGG_results_all <- pairwise_termsim(enrichKEGG_results_all)
emapplot(pairwise_term_enrichKEGG_results_all, showCategory = 15) + ggtitle("Enrichment map for enrichKEGG_results_all")
dev.off()

# pdf(paste0(result_dir, "/enrichKEGG_results_all_exonic.pdf"), width = 12, height = 15)
# dotplot(enrichKEGG_results_all_exonic, showCategory=15, title = "enrichKEGG_results_all_exonic dotplot") 
# barplot(enrichKEGG_results_all_exonic, showCategory=15, title = "enrichKEGG_results_all_exonic barplot") 
# pairwise_term_enrichKEGG_results_all_exonic <- pairwise_termsim(enrichKEGG_results_all_exonic)
# emapplot(pairwise_term_enrichKEGG_results_all_exonic, showCategory = 15) + ggtitle("Enrichment map for enrichKEGG_results_all_exonic")
# dev.off()
