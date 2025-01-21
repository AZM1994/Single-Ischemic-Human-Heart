library(ggplot2)
library(stringr)
library(readxl)
library(ggsci)
library(dplyr)
library(tidyr)
library(tibble)
library(reshape2)
library(Seurat)
library(pheatmap)
library(readxl)
library(ggpubr)
# ref_genome="BSgenome.Hsapiens.UCSC.hg19"
# library(ref_genome, character.only = T)
library(MutationalPatterns)
# chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

setwd("/Users/zhemingan/Documents/BCH_research/Gene_Expression_Analysis")
group_num <- 8

##### read in metadata
Hypoxia_PTA_Cases_metadata <- readRDS("./data/SCAN2_df.rds") %>%
  as.data.frame() |> base::`[`(c("cell_ID", "age", "gender", "Case_ID", "condition", "snv.burden", "snv.rate.per.gb")) %>%
  rename_with(~ c("Cell_ID", "Condition"), .cols = c(1, 5))
Hypoxia_PTA_Cases_metadata_collapsed <- Hypoxia_PTA_Cases_metadata %>% distinct(Case_ID, .keep_all = TRUE)
Cell_ID_list <- Hypoxia_PTA_Cases_metadata$Cell_ID
Case_ID_list <- Hypoxia_PTA_Cases_metadata_collapsed$Case_ID
Condition_list <- unique(Hypoxia_PTA_Cases_metadata$Condition)
genomic_context_colnames <- c("Cell_ID", "Case_ID", "Condition", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene")

##### read in scRNA-seq data
all_celltype_RNAseq <- readRDS("./data/Seurat.obj_with_annotation.RDS")
CM_cells <- subset(all_celltype_RNAseq, subset = annotated_clusters == "Cardiomyocytes")

analysis_ID <- "generate_low_high_gene_list"
fig_save_dir <- paste0("./figures/transcription_analysis/", analysis_ID, "/")
dir.create(fig_save_dir)

# condition_temp = Condition_list[1]
for (condition_temp in Condition_list){
  ##### get transcription data
  cat("##### Get transcription data for:", condition_temp, "...\n")
  expr_level_temp <- data.frame(AverageExpression(CM_cells, group.by = "condition", slot = "data")$RNA) %>% 
    setNames(Condition_list) |> base::`[`(condition_temp) %>% 
    mutate(gene = row.names(.)) %>% 
    setNames(c("average_expr_level", "gene")) %>% 
    # filter(average_expr_level != 0) %>% 
    mutate(decile = ntile(average_expr_level, n = group_num)) %>% mutate(decile = as.factor(decile))
  
  ###########################################################################
  ##################### raw SCAN2 call mutation analysis ####################
  ###########################################################################
  cat("##### Raw SCAN2 call mutation analysis:", condition_temp, "...\n")
  cat("Get genomic context for", condition_temp, "...\n")
  Case_ID_order <- Hypoxia_PTA_Cases_metadata_collapsed[Hypoxia_PTA_Cases_metadata_collapsed$Condition == condition_temp, "Case_ID"]
  heart_PTA_Cases_vcf_temp <- read.table(paste0("data/heart_PTA_Cases.all_age.", condition_temp, "_ssnv.vcf"), sep = "\t")
  genomic_context_temp <- read.csv(paste0("data/heart_PTA_Cases.all_age.", condition_temp, ".annotation.csv"), header = TRUE) %>%
    mutate(Cell_ID = heart_PTA_Cases_vcf_temp$V8) %>% mutate(Case_ID = str_extract(Cell_ID, "[^_]+")) %>% 
    mutate(Condition = condition_temp) |> base::`[`(genomic_context_colnames)
  
  genomic_context_temp <- genomic_context_temp %>%
    mutate(Cell_ID = ifelse(Cell_ID == "1864_M-E3-2n", "1864_E3", ifelse(Cell_ID == "4402_1_A1-2n", "4402_A1", ifelse(Cell_ID == "4402_1_A3-2n", "4402_A3",
    ifelse(Cell_ID == "4638_1_A1-2n", "4638_A1", ifelse(Cell_ID == "4638_1_A4-2n", "4638_A4", ifelse(Cell_ID == "4638_1_B2-2n", "4638_B2",
    ifelse(Cell_ID == "5657_M-D1-2n", "5657_D1", ifelse(Cell_ID == "5657_M-H1-2n", "5657_H1", ifelse(Cell_ID == "5828_CM_C2_2n", "5828_C2",
    ifelse(Cell_ID == "5828_CM_G2_2n", "5828_G2", ifelse(Cell_ID == "5919_1_C4-2n", "5919_C4", ifelse(Cell_ID == "5919_1_D2-2n", "5919_D2",
    ifelse(Cell_ID == "5919_1_E3-2n", "5919_E3", ifelse(Cell_ID == "5919_1_F6-2n", "5919_F6", ifelse(Cell_ID == "6032_1_A1-2n", "6032_A1",
    ifelse(Cell_ID == "6032_M-C7-2n", "6032_C7", ifelse(Cell_ID == "6032_M-E7-2n", "6032_E7", ifelse(Cell_ID == "1673_CM_A2_2n", "1673_A2",
    ifelse(Cell_ID == "1673_CM_A3_2n", "1673_A3", ifelse(Cell_ID == "1673_CM_D2_2n", "1673_D2", Cell_ID)))))))))))))))))))))
  
  genomic_SCAN2_df_temp <- merge(genomic_context_temp, Hypoxia_PTA_Cases_metadata[c("Cell_ID", "Case_ID", "age", "gender", "Condition")]) %>% 
    rename_with(~ c("Cell_ID", "Case_ID", "Condition", "chr", "start", "end", "ref", "alt", "region", "gene"), .cols = 1:10) %>% 
    mutate(Condition = as.factor(Condition))

  # genic_mutation_temp <- genomic_SCAN2_df_temp
  genic_mutation_temp <- genomic_SCAN2_df_temp[genomic_SCAN2_df_temp$region %in% c("exonic", "exonic;splicing", "intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3"), ] %>%
    mutate(gene = str_remove(gene, "\\(.*\\)$")) %>% filter(!str_detect(gene, ","))
  # table(genomic_SCAN2_df_temp$Case_ID)
  # table(genic_mutation_temp$Case_ID)
  # mutation_num_temp <- data.frame(table(genic_mutation_temp$gene, genic_mutation_temp$Case_ID)) %>% 
  #   setNames(c("gene", "Case_ID", "mut_number"))
  mutation_num_temp <- data.frame(table(genic_mutation_temp$gene)) %>% 
    setNames(c("gene", "mut_number"))
  
  expr_level_mutation_temp <- inner_join(expr_level_temp, mutation_num_temp, by = "gene")
  
  expr_level_mutation_temp_ordered <- expr_level_mutation_temp %>% 
    arrange(average_expr_level)
  # length(unique(expr_level_mutation_temp$gene))
  # length(unique(expr_level_mutation_temp$gene[expr_level_mutation_temp$decile == 1]))
  # condition_temp <- unique(expr_level_mutation_temp$gene[expr_level_mutation_temp$decile == 1])
  write.csv(expr_level_mutation_temp_ordered[, 1:3], paste0(fig_save_dir, "expr_level_SNV.", condition_temp, ".csv"))
}
