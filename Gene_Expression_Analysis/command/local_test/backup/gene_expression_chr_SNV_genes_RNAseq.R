library(ggplot2)
library(ggpubr)
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
ref_genome="BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = T)
library(MutationalPatterns)
# chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

seqlengths_list <- seqlengths(Hsapiens)[1:22]
setwd("/Users/zhemingan/Documents/BCH_research/Gene_Expression_Analysis")
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])
group_num <- 8
batch_size <- 1

##### read hypoxia gene list
# hypoxia_gene_list <- read.csv("data/signet_sup/supplementary_table_5.csv", header = TRUE) %>% 
#   pivot_longer(cols = everything(), names_to = "gene_set", values_to = "genes")
# Elvidge_hypoxia_gene_list <- hypoxia_gene_list$Elvidge

##### read in metadata
Hypoxia_PTA_Cases_metadata <- readRDS("./data/SCAN2_df.rds") %>% 
  as.data.frame() |> base::`[`(c("cell_ID", "age", "gender", "Case_ID", "condition", "snv.burden", "snv.rate.per.gb")) %>%
  rename_with(~ c("Cell_ID", "Condition"), .cols = c(1, 5))
Hypoxia_PTA_Cases_metadata_age_match <- Hypoxia_PTA_Cases_metadata %>% 
  filter(age >= 40 & age < 80)
Hypoxia_PTA_Cases_metadata_collapsed <- Hypoxia_PTA_Cases_metadata %>% distinct(Case_ID, .keep_all = TRUE)
Cell_ID_list <- Hypoxia_PTA_Cases_metadata$Cell_ID
Case_ID_list <- Hypoxia_PTA_Cases_metadata_collapsed$Case_ID
Condition_list <- unique(Hypoxia_PTA_Cases_metadata$Condition)
genomic_context_colnames <- c("Cell_ID", "Case_ID", "Condition", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene")

##### read in scRNA-seq data
all_celltype_RNAseq <- readRDS("./data/Seurat.obj_with_annotation.RDS")
CM_cells <- subset(all_celltype_RNAseq, subset = annotated_clusters == "Cardiomyocytes")

genomic_SCAN2_df <- c()
for (condition_temp in Condition_list){
  cat("Get genomic context for", condition_temp, "...\n")
  heart_PTA_Cases_vcf_temp <- read.table(paste0("data/heart_PTA_Cases_annovar/heart_PTA_Cases.all_age.", condition_temp, "_ssnv.vcf"), sep = "\t")
  genomic_context_temp <- read.csv(paste0("data/heart_PTA_Cases_annovar/heart_PTA_Cases.all_age.", condition_temp, "_ssnv.csv"), header = TRUE) %>%
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
    filter(age >= 40 & age < 80)
  
  genic_mutation_temp <- genomic_SCAN2_df_temp
  # genic_mutation_temp <- genomic_SCAN2_df_temp[genomic_SCAN2_df_temp$region %in% c("exonic", "exonic;splicing", "intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3"), ] %>%
  #   mutate(gene = str_remove(gene, "\\(.*\\)$")) %>% filter(!str_detect(gene, ","))
  genomic_SCAN2_df <- rbind(genomic_SCAN2_df, genomic_SCAN2_df_temp)
}

# condition_temp = Condition_list[1]
cell_chr_expr_all <- c()
cell_chr_expr_all_chr <- c()
for (chr_index in seq(1:22)){
  all_SNV_gene_df_temp <- genomic_SCAN2_df[genomic_SCAN2_df$chr == chr_index, ] %>% filter(!duplicated(gene))
  # print(length(all_SNV_gene_list))
  for (cell_ID in Cell_ID_list){
    cell_chr_SNV_df <- genomic_SCAN2_df[genomic_SCAN2_df$Cell_ID == cell_ID & genomic_SCAN2_df$chr == chr_index, ]
    if (cell_ID %in% Hypoxia_PTA_Cases_metadata$Cell_ID[Hypoxia_PTA_Cases_metadata$Condition == "Normal"]){
      expr_level_temp <- data.frame(AverageExpression(CM_cells, group.by = "condition", slot = "data")$RNA) %>% 
        setNames(Condition_list) |> base::`[`("Normal") %>% 
        mutate(gene = row.names(.)) %>% 
        setNames(c("average_expr_level", "gene")) %>% 
        mutate(decile = ntile(average_expr_level, n = group_num)) %>% mutate(decile = as.factor(decile)) 
      all_SNV_gene_expr_df_temp <- inner_join(expr_level_temp, all_SNV_gene_df_temp) %>% 
        mutate(Condition = "Normal") %>% mutate(Cell_ID = cell_ID)
      cell_chr_expr <- all_SNV_gene_expr_df_temp %>% 
        mutate(average_expr_level = ifelse(gene %in% cell_chr_SNV_df$gene, average_expr_level, 0))
    }else if(cell_ID %in% Hypoxia_PTA_Cases_metadata$Cell_ID[Hypoxia_PTA_Cases_metadata$Condition == "Disease"]){
      expr_level_temp <- data.frame(AverageExpression(CM_cells, group.by = "condition", slot = "data")$RNA) %>% 
        setNames(Condition_list) |> base::`[`("Disease") %>% 
        mutate(gene = row.names(.)) %>% 
        setNames(c("average_expr_level", "gene")) %>% 
        mutate(decile = ntile(average_expr_level, n = group_num)) %>% mutate(decile = as.factor(decile)) 
      all_SNV_gene_expr_df_temp <- inner_join(expr_level_temp, all_SNV_gene_df_temp) %>% 
        mutate(Condition = "Disease") %>% mutate(Cell_ID = cell_ID) %>% mutate(chr = chr_index)
      cell_chr_expr <- all_SNV_gene_expr_df_temp %>% 
        mutate(average_expr_level = ifelse(gene %in% cell_chr_SNV_df$gene, average_expr_level, 0))
    }
    cell_chr_expr_all <- rbind(cell_chr_expr_all, cell_chr_expr)
  }
  cell_chr_expr_all_chr <- rbind(cell_chr_expr_all_chr, cell_chr_expr_all)
}

expr_mutnum_chr_by_cell <- c()
# genomic_SCAN2_df <- c()
# condition_temp = Condition_list[1]
for (condition_temp in Condition_list){
  ##### get transcription data
  cat("##### Get transcription data for:", condition_temp, "...\n")
  expr_level_temp <- data.frame(AverageExpression(CM_cells, group.by = "condition", slot = "data")$RNA) %>% 
    setNames(Condition_list) |> base::`[`(condition_temp) %>% 
    mutate(gene = row.names(.)) %>% 
    setNames(c("average_expr_level", "gene")) %>% 
    mutate(decile = ntile(average_expr_level, n = group_num)) %>% mutate(decile = as.factor(decile))
  
  genomic_SCAN2_df
  
  expr_chr_by_cell_temp <- inner_join(expr_level_temp, genic_mutation_temp, by = "gene") |> base::`[`(c("gene", "average_expr_level", "chr", "Condition", "Cell_ID", "decile"))
    # group_by(Condition, chr, Cell_ID) %>% 
    # summarise("avg_expr_per_chr" = mean(average_expr_level))
  
  expr_mutnum_chr_by_cell_temp <- merge(mutation_num_by_cell_temp, expr_chr_by_cell_temp) %>% arrange(chr)
  
  expr_mutnum_chr_by_cell <- rbind(expr_mutnum_chr_by_cell, expr_mutnum_chr_by_cell_temp)
  genomic_SCAN2_df <- rbind(genomic_SCAN2_df, genomic_SCAN2_df_temp)
}



### 
fig_save_dir <- paste0("./results/local_test/gene_expression_chr/")
dir.create(fig_save_dir, recursive = TRUE)

### by condition, count number of genes in each expression level
expected_num_genes_chr <- expr_mutnum_chr_by_cell %>% 
  group_by(Condition) %>% 
  summarise(count = n()) %>% 
  mutate(count_exp = count / (group_num * 22))

expr_mutnum_summary_by_decile_chr <- expr_mutnum_chr_by_cell %>% 
  group_by(Condition, chr, decile) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  complete(Condition, chr, decile, fill = list(count = 0)) %>% 
  mutate(expected_num_genes_chr = rep(expected_num_genes_chr$count_exp, each = group_num * 22)) %>% 
  mutate(enrich_ratio = count / expected_num_genes_chr)

pdf(paste0(fig_save_dir, "/gene_count_chr.pdf"), width = 8, height = 5.5)
for (chr_index in 1:22){
  data_plot <- expr_mutnum_summary_by_decile_chr[expr_mutnum_summary_by_decile_chr$chr == chr_index, ]
  p_expr_mutnum_summary_by_decile_chr <- ggplot(data_plot, aes(x = decile, y = enrich_ratio, group = Condition, color = Condition)) + 
    geom_hline(yintercept = 1, color = "black", linewidth = 0.6) + 
    geom_line(position = position_dodge(width = 0.1), size = 1) + 
    geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE) + 
    # geom_errorbar(aes(ymin = mean_ER - sd_ER, ymax = mean_ER + sd_ER), width = 0.2, position = position_dodge(width = 0.1)) +
    # geom_smooth(data = merged_summary_plot, aes(x = decile, y = mean_ER, color = Condition, fill = Condition, group = Condition),
    #             method = "lm", se = TRUE, alpha = 0.2, linewidth = 1, linetype = "dashed") +
    theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position="right") +
    theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
    # ylim(c(0.5, 1.6)) + 
    labs(x = "Gene expression levels", y = "Gene number enrichment ratio \n (obs/exp)", color = "Condition", title = chr_index)
  print(p_expr_mutnum_summary_by_decile_chr)
  # ggsave(paste0(fig_save_dir, "/chr", chr_index, "gene_expr_chr.pdf"), width = 8, height = 5)
}
dev.off()

### by condition by chr, calculate the average expression level (median by expression levels) 
# summary_expr_by_condition_chr <- expr_mutnum_chr_by_cell %>% 
#   group_by(Condition, chr) %>% 
#   summarise(count = n(), 
#             mean_expr = mean(average_expr_level), 
#             median_expr = median(average_expr_level), 
#             CI_lower = t.test(average_expr_level)$conf.int[1],
#             CI_upper = t.test(average_expr_level)$conf.int[2])

summary_expr_by_cell_condition_chr <- expr_mutnum_chr_by_cell %>% 
  group_by(Cell_ID, chr) %>% 
  summarise(count = n(), 
            mean_expr = mean(average_expr_level), 
            median_expr = median(average_expr_level)) %>% 
  ungroup() %>% 
  filter(count > 1) %>% 
  mutate(Condition = ifelse(Cell_ID %in% Hypoxia_PTA_Cases_metadata$Cell_ID[Hypoxia_PTA_Cases_metadata$Condition == "Normal"], "Normal", 
                            ifelse(Cell_ID %in% Hypoxia_PTA_Cases_metadata$Cell_ID[Hypoxia_PTA_Cases_metadata$Condition == "Disease"], "Disease", "none"))) %>% 
  mutate(Condition = factor(Condition, levels = c("Normal", "Disease")))

pdf(paste0(fig_save_dir, "/gene_expr_median_by_cell_condition_chr_line.pdf"), width = 8, height = 5.5)
p_expr_by_cell_condition_chr <- ggplot(summary_expr_by_cell_condition_chr, aes(x = chr, y = median_expr, group = Condition, color = Condition, fill = Condition)) + 
  geom_line(stat = "summary", fun = "mean") + 
  geom_point(stat = "summary", fun = "mean") +
  theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right") +
  theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
  # ylim(c(0.5, 1.6)) + 
  labs(x = "Chromosome", y = "Gene expression (median)", color = "Condition", title = "")
print(p_expr_by_cell_condition_chr)
dev.off()

pdf(paste0(fig_save_dir, "/gene_expr_median_by_cell_condition_chr_bar.pdf"), width = 18, height = 5.5)
p_expr_by_cell_condition_chr <- ggplot(summary_expr_by_cell_condition_chr, aes(x = chr, y = median_expr, group = Condition, color = Condition, fill = Condition)) + 
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8)) + 
  geom_errorbar(stat = "summary", fun.data = mean_cl_normal, position = position_dodge(width = 0.8), width = 0.1) +
  stat_compare_means(aes(group = Condition), method = "wilcox.test", label = "p.format", 
                     label.y = 1.02 * max(summary_expr_by_cell_condition_chr$median_expr)) + 
  stat_compare_means(aes(group = Condition), method = "wilcox.test", label = "p.signif", 
                     label.y = 1.06 * max(summary_expr_by_cell_condition_chr$median_expr)) +
  geom_jitter(position = position_dodge(width = 0.05), size = 1) + 
  theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right") +
  theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
  # ylim(c(0.5, 1.6)) + 
  labs(x = "Chromosome", y = "Gene expression (median)", color = "Condition", title = "")
print(p_expr_by_cell_condition_chr)
dev.off()

# ### check gene expression by condition, chr, and decile
# summary_expr_by_condition_chr_decile <- expr_mutnum_chr_by_cell %>% 
#   group_by(Condition, chr, decile) %>% 
#   summarise(count = n(), 
#             mean_expr = mean(average_expr_level), 
#             median_expr = median(average_expr_level)) %>% 
#   ungroup() %>% 
#   complete(Condition, chr, decile, fill = list(count = 0, mean_expr = 0, median_expr = 0)) 
# 
# pdf(paste0(fig_save_dir, "/gene_expr_chr_decile_mean.pdf"), width = 8, height = 5.5)
# for (chr_index in 1:22){
#   data_plot <- summary_expr_by_condition_chr_decile[summary_expr_by_condition_chr_decile$chr == chr_index, ]
#   p_expr_by_condition_chr_decile <- ggplot(data_plot, aes(x = decile, y = mean_expr, group = Condition, color = Condition)) + 
#     geom_line(position = position_dodge(width = 0.1), size = 1) + 
#     geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE) + 
#     theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
#           panel.background = element_rect(fill = "white"), legend.position="right") +
#     theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
#     labs(x = "Gene expression levels", y = "Gene expression (mean)", color = "Condition", title = chr_index)
#   print(p_expr_by_condition_chr_decile)
# }
# dev.off()
# 
# pdf(paste0(fig_save_dir, "/gene_expr_chr_decile_median.pdf"), width = 8, height = 5.5)
# for (chr_index in 1:22){
#   # chr_index = 12
#   data_plot <- summary_expr_by_condition_chr_decile[summary_expr_by_condition_chr_decile$chr == chr_index, ]
#   p_expr_by_condition_chr_decile <- ggplot(data_plot, aes(x = decile, y = median_expr, group = Condition, color = Condition)) + 
#     geom_line(position = position_dodge(width = 0.1), size = 1) + 
#     geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE) + 
#     theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
#           panel.background = element_rect(fill = "white"), legend.position="right") +
#     theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
#     labs(x = "Gene expression levels", y = "Gene expression (median)", color = "Condition", title = chr_index)
#   print(p_expr_by_condition_chr_decile)
# }
# dev.off()
