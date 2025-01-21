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

wd_dir <- getwd()
result_dir <- paste0(wd_dir, "/results/LOO_DEG_Cardiomyocyte/enrichGO")
dir.create(result_dir, recursive = T)

##### read DEGs
DEG_df <- read.csv(paste0(wd_dir, "/results/LOO_DEG_Cardiomyocyte/DEG_up_down_df.csv"))
DEG_up <- DEG_df$gene[DEG_df$regulation == "up"]
DEG_down <- DEG_df$gene[DEG_df$regulation == "down"]

#################################################################################
######################## GO over-representation analysis ########################
#################################################################################
GO_DEG_up <- enrichGO(gene = DEG_up, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                          pAdjustMethod = "fdr", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")
GO_DEG_down <- enrichGO(gene = DEG_down, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.05, 
                           pAdjustMethod = "fdr", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, keyType = "SYMBOL")

pdf(paste0(result_dir, "/GO_results.pdf"), width = 12, height = 15)
  barplot(GO_DEG_up, showCategory = 30, title = "GO_DEG_up barplot")
  barplot(GO_DEG_down, showCategory = 30, title = "GO_DEG_down barplot")
dev.off()

#################################################################################
############################### Pathway analysis KEGG ################################
#################################################################################
kegg_DEG_up <- mapIds(org.Hs.eg.db, DEG_up, 'ENTREZID', 'SYMBOL')
kegg_DEG_down <- mapIds(org.Hs.eg.db, DEG_down, 'ENTREZID', 'SYMBOL')

enrichKEGG_all_normal <- enrichKEGG(gene = kegg_DEG_up, organism = "hsa", pvalueCutoff = 0.05, 
                                    pAdjustMethod = "fdr", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, keyType = "kegg")
enrichKEGG_all_disease <- enrichKEGG(gene = kegg_DEG_down, organism = "hsa", pvalueCutoff = 0.05, 
                                    pAdjustMethod = "fdr", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, keyType = "kegg")

pdf(paste0(result_dir, "/enrichKEGG_results.pdf"), width = 12, height = 15)
  barplot(enrichKEGG_all_normal, showCategory = 30, title = "KEGG_all_normal barplot")
  barplot(enrichKEGG_all_disease, showCategory = 30, title = "KEGG_all_disease barplot")
dev.off()

#################################################################################
############################### ENRICHR analysis ################################
#################################################################################
# databases_list <- c("KEGG_2019_Human", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
databases_list <- c("KEGG_2019_Human", "GO_Biological_Process_2018")
enrichR_all_normal <- enrichr(DEG_up, databases = databases_list)
enrichR_all_disease <- enrichr(DEG_down, databases = databases_list)

pdf(paste0(result_dir, "/enrichR_GO_results.pdf"), width = 15, height = 8)
  plotEnrich(enrichR_all_normal[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "GO_DEG_up barplot") 
  plotEnrich(enrichR_all_disease[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "GO_DEG_down barplot")
dev.off()

pdf(paste0(result_dir, "/enrichR_KEGG_results.pdf"), width = 20, height = 8)
plotEnrich(enrichR_all_normal[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "KEGG_all_normal barplot") + 
  plotEnrich(enrichR_all_disease[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "KEGG_all_disease barplot")
dev.off()