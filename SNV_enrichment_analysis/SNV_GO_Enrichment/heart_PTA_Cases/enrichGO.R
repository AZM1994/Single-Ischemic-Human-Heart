library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DESeq2)
library(enrichplot)
library(ggplot2)
library(stringr)
library(DOSE)
library(dplyr)
library(tidyr)
library(enrichR)

setwd("/Users/zhemingan/Documents/BCH_research/annovar/heart_PTA_Cases/GOseq_Zheming")
wd_dir <- getwd()
result_dir <- paste0(wd_dir, "/enrichGO")
dir.create(result_dir)

### read gene list with metadata
Hypoxia_PTA_Cases_metadata <- readRDS("SCAN2_df.rds") %>% as.data.frame() |> 
  base::`[`(c("cell_ID", "age", "gender", "Case_ID", "condition", "snv.burden", "snv.rate.per.gb")) %>% 
  rename_with(~ c("Cell_ID", "Condition"), .cols = c(1, 5))
Condition_list <- unique(Hypoxia_PTA_Cases_metadata$Condition)
selected_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", 
                       "ExonicFunc.refGene", "AAChange.refGene", "Cell_ID", "Case_ID", "Condition", "mut_type", "age")

genomic_SCAN2_df <- c()
genomic_SCAN2_del_df <- c()
for (condition_temp in Condition_list) {
  for (mutation_type in c("ssnv", "sindel")) {
    cat("Get genomic context for", condition_temp, mutation_type, "...\n")
    heart_PTA_Cases_vcf_temp <- read.table(paste0("heart_PTA_Cases_annovar/heart_PTA_Cases.all_age.", condition_temp, "_", mutation_type, ".vcf"), sep = "\t") %>% mutate(V8 = sub(";.*", "", V8))
    genomic_context_temp <- read.csv(paste0("heart_PTA_Cases_annovar/heart_PTA_Cases.all_age.", condition_temp, "_", mutation_type, ".csv"), header = TRUE) %>%
      mutate(Cell_ID = heart_PTA_Cases_vcf_temp$V8) %>% mutate(Case_ID = str_extract(Cell_ID, "[^_]+")) %>% 
      mutate(Condition = condition_temp) %>% mutate(mut_type = mutation_type)
    
    genomic_context_temp <- genomic_context_temp %>%
      mutate(Cell_ID = ifelse(Cell_ID == "1864_M-E3-2n", "1864_E3", ifelse(Cell_ID == "4402_1_A1-2n", "4402_A1", ifelse(Cell_ID == "4402_1_A3-2n", "4402_A3",
      ifelse(Cell_ID == "4638_1_A1-2n", "4638_A1", ifelse(Cell_ID == "4638_1_A4-2n", "4638_A4", ifelse(Cell_ID == "4638_1_B2-2n", "4638_B2",
      ifelse(Cell_ID == "5657_M-D1-2n", "5657_D1", ifelse(Cell_ID == "5657_M-H1-2n", "5657_H1", ifelse(Cell_ID == "5828_CM_C2_2n", "5828_C2",
      ifelse(Cell_ID == "5828_CM_G2_2n", "5828_G2", ifelse(Cell_ID == "5919_1_C4-2n", "5919_C4", ifelse(Cell_ID == "5919_1_D2-2n", "5919_D2",
      ifelse(Cell_ID == "5919_1_E3-2n", "5919_E3", ifelse(Cell_ID == "5919_1_F6-2n", "5919_F6", ifelse(Cell_ID == "6032_1_A1-2n", "6032_A1",
      ifelse(Cell_ID == "6032_M-C7-2n", "6032_C7", ifelse(Cell_ID == "6032_M-E7-2n", "6032_E7", ifelse(Cell_ID == "1673_CM_A2_2n", "1673_A2",
      ifelse(Cell_ID == "1673_CM_A3_2n", "1673_A3", ifelse(Cell_ID == "1673_CM_D2_2n", "1673_D2", Cell_ID)))))))))))))))))))))
    
    genomic_SCAN2_df_temp <- merge(genomic_context_temp, Hypoxia_PTA_Cases_metadata[c("Cell_ID", "Case_ID", "Condition", "age")]) |> base::`[`(selected_colnames) %>% 
      filter(!grepl(";", Gene.refGene))
    # filter(age >= 40 & age < 80) %>% filter(!grepl(";", Gene.refGene))
    genomic_SCAN2_df_temp_del <- genomic_SCAN2_df_temp %>% 
      filter(Func.refGene %in% c("splicing", "exonic;splicing") | ExonicFunc.refGene %in% c("frameshift deletion", "frameshift insertion", "nonframeshift deletion", "nonsynonymous SNV", "stopgain"))
    
    genomic_SCAN2_df <- rbind(genomic_SCAN2_df, genomic_SCAN2_df_temp)
    genomic_SCAN2_del_df <- rbind(genomic_SCAN2_del_df, genomic_SCAN2_df_temp_del)
  }
}

all_gene_normal <- genomic_SCAN2_df$Gene.refGene[genomic_SCAN2_df$Condition == "Normal"]
all_gene_disease <- genomic_SCAN2_df$Gene.refGene[genomic_SCAN2_df$Condition == "Disease"]
del_gene_normal <- genomic_SCAN2_del_df$Gene.refGene[genomic_SCAN2_del_df$Condition == "Normal"]
del_gene_disease <- genomic_SCAN2_del_df$Gene.refGene[genomic_SCAN2_del_df$Condition == "Disease"]
write.csv(all_gene_normal, "all_gene_normal.csv")
write.csv(all_gene_disease, "all_gene_disease.csv")
write.csv(del_gene_normal, "del_gene_normal.csv")
write.csv(del_gene_disease, "del_gene_disease.csv")
#################################################################################
######################## GO over-representation analysis ########################
#################################################################################
GO_all_normal <- enrichGO(gene = all_gene_normal, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                          pAdjustMethod = "fdr", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")
GO_all_disease <- enrichGO(gene = all_gene_disease, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                           pAdjustMethod = "fdr", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")

GO_del_normal <- enrichGO(gene = del_gene_normal, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                          pAdjustMethod = "fdr", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")
GO_del_disease <- enrichGO(gene = del_gene_disease, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                           pAdjustMethod = "fdr", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")

pdf(paste0(result_dir, "/GO_results.pdf"), width = 12, height = 15)
barplot(GO_all_normal, showCategory = 20, title = "GO_all_normal barplot")
barplot(GO_all_disease, showCategory = 20, title = "GO_all_disease barplot")
barplot(GO_del_normal, showCategory = 20, title = "GO_del_normal barplot")
# barplot(GO_del_disease, showCategory = 20, title = "GO_del_disease barplot")
dev.off()

#################################################################################
############################### Pathway analysis KEGG ################################
#################################################################################
kegg_all_gene_normal <- mapIds(org.Hs.eg.db, all_gene_normal, 'ENTREZID', 'SYMBOL')
kegg_all_gene_disease <- mapIds(org.Hs.eg.db, all_gene_disease, 'ENTREZID', 'SYMBOL')
kegg_del_gene_normal <- mapIds(org.Hs.eg.db, del_gene_normal, 'ENTREZID', 'SYMBOL')
kegg_del_gene_disease <- mapIds(org.Hs.eg.db, del_gene_disease, 'ENTREZID', 'SYMBOL')

enrichKEGG_all_normal <- enrichKEGG(gene = kegg_all_gene_normal, organism = "hsa", pvalueCutoff = 0.05, 
                                    pAdjustMethod = "fdr", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, keyType = "kegg")
enrichKEGG_all_disease <- enrichKEGG(gene = kegg_all_gene_disease, organism = "hsa", pvalueCutoff = 0.05, 
                                    pAdjustMethod = "fdr", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, keyType = "kegg")

enrichKEGG_del_normal <- enrichKEGG(gene = kegg_del_gene_normal, organism = "hsa", pvalueCutoff = 0.05, 
                                    pAdjustMethod = "fdr", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, keyType = "kegg")
enrichKEGG_del_disease <- enrichKEGG(gene = kegg_del_gene_disease, organism = "hsa", pvalueCutoff = 0.05, 
                                     pAdjustMethod = "fdr", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, keyType = "kegg")

pdf(paste0(result_dir, "/enrichKEGG_results.pdf"), width = 12, height = 15)
barplot(enrichKEGG_all_normal, showCategory = 20, title = "KEGG_all_normal barplot")
barplot(enrichKEGG_all_disease, showCategory = 20, title = "KEGG_all_disease barplot")
# barplot(enrichKEGG_del_normal, showCategory = 20, title = "KEGG_del_normal barplot")
# barplot(enrichKEGG_del_disease, showCategory = 20, title = "KEGG_del_disease barplot")
dev.off()

#################################################################################
############################### ENRICHR analysis ################################
#################################################################################
# databases_list <- c("KEGG_2019_Human", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
databases_list <- c("KEGG_2019_Human", "GO_Biological_Process_2018")
enrichR_all_normal <- enrichr(all_gene_normal, databases = databases_list)
enrichR_all_disease <- enrichr(all_gene_disease, databases = databases_list)
enrichR_del_normal <- enrichr(del_gene_normal, databases = databases_list)
enrichR_del_disease <- enrichr(del_gene_disease, databases = databases_list)

pdf(paste0(result_dir, "/enrichR_GO_results.pdf"), width = 15, height = 8)
plotEnrich(enrichR_all_normal[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "GO_all_normal barplot") 
plotEnrich(enrichR_all_disease[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "GO_all_disease barplot")

plotEnrich(enrichR_del_normal[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "GO_del_normal barplot") 
plotEnrich(enrichR_del_disease[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "GO_del_disease barplot")
dev.off()

pdf(paste0(result_dir, "/enrichR_KEGG_results.pdf"), width = 20, height = 8)
plotEnrich(enrichR_all_normal[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "KEGG_all_normal barplot") + 
  plotEnrich(enrichR_all_disease[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "KEGG_all_disease barplot")

plotEnrich(enrichR_del_normal[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "KEGG_del_normal barplot") + 
  plotEnrich(enrichR_del_disease[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "KEGG_del_disease barplot")
dev.off()