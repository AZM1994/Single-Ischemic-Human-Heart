# library(ggplot2)
# library(ggpubr)
# library(stringr)
# library(readxl)
# library(ggsci)
# library(dplyr)
# library(tidyr)
# library(tibble)
# library(reshape2)
# library(Seurat)
# library(pheatmap)
# library(readxl)
# ref_genome="BSgenome.Hsapiens.UCSC.hg19"
# library(ref_genome, character.only = T)
# library(MutationalPatterns)
library(ggplot2)
library(ggpubr)
# library(scales)
library(dplyr)
library(Seurat)

# seqlengths_list <- seqlengths(Hsapiens)[1:22]
setwd("/Users/zhemingan/Documents/BCH_research/Gene_Expression_Analysis")
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])
group_num <- 8

##### read in metadata
# Hypoxia_PTA_Cases_metadata <- readRDS("./data/SCAN2_df.rds") %>% 
#   as.data.frame() |> base::`[`(c("cell_ID", "age", "gender", "Case_ID", "condition", "snv.burden", "snv.rate.per.gb")) %>%
#   rename_with(~ c("Cell_ID", "Condition"), .cols = c(1, 5))
# Hypoxia_PTA_Cases_metadata_age_match <- Hypoxia_PTA_Cases_metadata %>% 
#   filter(age >= 40 & age < 80)
# Hypoxia_PTA_Cases_metadata_collapsed <- Hypoxia_PTA_Cases_metadata %>% distinct(Case_ID, .keep_all = TRUE)
# Cell_ID_list <- Hypoxia_PTA_Cases_metadata$Cell_ID
# Case_ID_list <- Hypoxia_PTA_Cases_metadata_collapsed$Case_ID
# Condition_list <- unique(Hypoxia_PTA_Cases_metadata$Condition)
# genomic_context_colnames <- c("Cell_ID", "Case_ID", "Condition", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene")

##### read in scRNA-ATACseq data
atac_normal <- c()
atac_disease <- c()
total_reads_4638 <- 41195236
total_reads_5828 <- 9955703
total_reads_5919 <- 33130651
total_reads_604 <- 17569221

total_reads_normal <- (total_reads_4638 + total_reads_5828 + total_reads_5919) / 3
total_reads_disease <- total_reads_604
reads_ratio <- total_reads_normal / total_reads_disease

atac_4638_2n <- read.table("./data/scATACseq_preprocessing/ATAC_cov_average_hg19_1kb_bin.4638_2n.bed") %>% setNames(c("chr", "start", "end", "id", "score_4638_2n"))
atac_5828_2n <- read.table("./data/scATACseq_preprocessing/ATAC_cov_average_hg19_1kb_bin.5828_2n.bed") %>% setNames(c("chr", "start", "end", "id", "score_5828_2n")) |> base::`[`(c("score_5828_2n"))
atac_5919_all <- read.table("./data/scATACseq_preprocessing/ATAC_cov_average_hg19_1kb_bin.5919_all.bed") %>% setNames(c("chr", "start", "end", "id", "score_5919_all")) |> base::`[`(c("score_5919_all"))
atac_normal <- cbind(atac_4638_2n, atac_5828_2n, atac_5919_all) %>% 
  mutate(score = rowMeans(select(., score_4638_2n, score_5828_2n, score_5919_all), na.rm = TRUE)) |> base::`[`(c("chr", "start", "end", "id", "score")) %>% 
  mutate(score = score / reads_ratio) %>% 
  mutate(group = ntile(score, n = group_num)) %>% 
  mutate(group = as.factor(group)) %>% 
  mutate(Condition = "Normal")
atac_grange_normal <- GRanges(seqnames = atac_normal$chr, 
                              ranges = IRanges(start = atac_normal$start, 
                                               end = atac_normal$end),
                              group = atac_normal$group) 

atac_604_all <- read.table("./data/scATACseq_preprocessing/ATAC_cov_average_hg19_1kb_bin.604_all.bed") %>% setNames(c("chr", "start", "end", "id", "score_604_all"))
atac_disease <- atac_604_all %>% 
  setNames(c("chr", "start", "end", "id", "score")) %>% 
  mutate(group = ntile(score, n = group_num)) %>% 
  mutate(group = as.factor(group)) %>% 
  mutate(Condition = "Disease")
atac_grange_disease <- GRanges(seqnames = atac_disease$chr, 
                               ranges = IRanges(start = atac_disease$start, 
                                                end = atac_disease$end),
                               group = atac_disease$group) 

gene_chr_access <- rbind(atac_normal, atac_disease)

### 
fig_save_dir <- paste0("./results/local_test/gene_expression_chr/")
dir.create(fig_save_dir, recursive = TRUE)

# hypoxia_gene_list <- read.csv("data/signet_sup/supplementary_table_5.csv", header = TRUE) %>% 
#   pivot_longer(cols = everything(), names_to = "gene_set", values_to = "genes") %>% filter(!duplicated(genes))
# hot_gene_list <- read.csv("data/MAFTools/hot_gene_list_heart_pta_02.csv", header = FALSE) %>% filter(!duplicated(V1))

gene_chr_access <- gene_chr_access %>% 
  # filter(chr %in% seq(1:22)) %>% 
  # filter(gene %in% hypoxia_gene_list$genes) %>% 
  mutate(Condition = factor(Condition, levels = c("Normal","Disease"))) %>% 
  mutate(chr = factor(chr, levels = paste0("chr", seq(1:22))))
### by condition, count number of genes in each expression level
expected_num_genes_chr <- gene_chr_access %>% 
  group_by(Condition) %>% 
  summarise(count = n()) %>% 
  mutate(count_exp = count / (group_num * 22))

expr_mutnum_summary_by_group_chr <- gene_chr_access %>% 
  group_by(Condition, chr, group) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  complete(Condition, chr, group, fill = list(count = 0)) %>% 
  mutate(expected_num_genes_chr = rep(expected_num_genes_chr$count_exp, each = group_num * 22)) %>% 
  mutate(enrich_ratio = count / expected_num_genes_chr)

pdf(paste0(fig_save_dir, "/chr_access_count_chr.pdf"), width = 8, height = 5.5)
for (chr_index in 1:22){
  data_plot <- expr_mutnum_summary_by_group_chr[expr_mutnum_summary_by_group_chr$chr == paste0("chr", chr_index), ]
  p_expr_mutnum_summary_by_group_chr <- ggplot(data_plot, aes(x = group, y = enrich_ratio, group = Condition, color = Condition)) + 
    geom_hline(yintercept = 1, color = "black", linewidth = 0.6) + 
    geom_line(position = position_dodge(width = 0.1), size = 1) + 
    geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE) + 
    # geom_errorbar(aes(ymin = mean_ER - sd_ER, ymax = mean_ER + sd_ER), width = 0.2, position = position_dodge(width = 0.1)) +
    # geom_smooth(data = merged_summary_plot, aes(x = group, y = mean_ER, color = Condition, fill = Condition, group = Condition),
    #             method = "lm", se = TRUE, alpha = 0.2, linewidth = 1, linetype = "dashed") +
    theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position="right") +
    theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
    # ylim(c(0.5, 1.6)) + 
    labs(x = "Chromosome accessibility levels", y = "Chromosome accessibility enrichment ratio \n (obs/exp)", color = "Condition", title = chr_index)
  print(p_expr_mutnum_summary_by_group_chr)
  # ggsave(paste0(fig_save_dir, "/chr", chr_index, "gene_expr_chr.pdf"), width = 8, height = 5)
}
dev.off()

### by condition by chr, calculate the average expression level (median by expression levels) 
summary_expr_by_condition_chr <- gene_chr_access %>% 
  group_by(Condition, chr) %>% 
  summarise(count = n(), 
            mean_expr = mean(score), 
            median_expr = median(score))

pdf(paste0(fig_save_dir, "/Chromosome_accessibility_median_by_condition_chr.pdf"), width = 10, height = 5.5)
p_expr_by_condition_chr <- ggplot(summary_expr_by_condition_chr, aes(x = chr, y = median_expr, group = Condition, color = Condition)) + 
  geom_line(position = position_dodge(width = 0.1), size = 1) + 
  geom_point(position = position_dodge(width = 0.1), size = 2) + 
  theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right") +
  theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
  labs(x = "Chromosome", y = "Chromosome accessibility (median)", color = "Condition", title = "")
print(p_expr_by_condition_chr)
dev.off()

pdf(paste0(fig_save_dir, "/Chromosome_accessibility_mean_by_condition_chr.pdf"), width = 10, height = 5.5)
p_expr_by_condition_chr <- ggplot(summary_expr_by_condition_chr, aes(x = chr, y = mean_expr, group = Condition, color = Condition)) + 
  geom_line(position = position_dodge(width = 0.1), size = 1) + 
  geom_point(position = position_dodge(width = 0.1), size = 2) + 
  theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right") +
  theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
  labs(x = "Chromosome", y = "Chromosome accessibility (mean)", color = "Condition", title = "")
print(p_expr_by_condition_chr)
dev.off()

summary_expr_by_cell_condition_chr <- gene_chr_access %>% 
  group_by(Cell_ID, chr) %>% 
  summarise(count = n(), 
            mean_expr = mean(score), 
            median_expr = median(score)) %>% 
  filter(count > 1) %>% 
  mutate(Condition = ifelse(Cell_ID %in% Hypoxia_PTA_Cases_metadata$Cell_ID[Hypoxia_PTA_Cases_metadata$Condition == "Normal"], "Normal", 
                            ifelse(Cell_ID %in% Hypoxia_PTA_Cases_metadata$Cell_ID[Hypoxia_PTA_Cases_metadata$Condition == "Disease"], "Disease", "none"))) %>% 
  mutate(Condition = factor(Condition, levels = c("Normal", "Disease")))

pdf(paste0(fig_save_dir, "/gene_expr_median_by_cell_condition_chr.pdf"), width = 18, height = 5.5)
p_expr_by_cell_condition_chr <- ggplot(summary_expr_by_cell_condition_chr, aes(x = chr, y = median_expr, group = Condition, 
                                                                               color = Condition, fill = Condition)) + 
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8)) + 
  geom_errorbar(stat = "summary", fun.data = mean_cl_normal, position = position_dodge(width = 0.8), width = 0.1) +
  stat_compare_means(aes(group = Condition), method = "wilcox.test", label = "p.format", 
                     label.y = 1.02 * max(summary_expr_by_cell_condition_chr$median_expr)) + 
  stat_compare_means(aes(group = Condition), method = "wilcox.test", label = "p.signif", 
                     label.y = 1.04 * max(summary_expr_by_cell_condition_chr$median_expr)) +
  geom_jitter(position = position_dodge(width = 0.05), size = 1) + 
  theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right") +
  theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
  # ylim(c(0.5, 1.6)) + 
  labs(x = "Chromosome", y = "Gene expression (median)", color = "Condition", title = "")
print(p_expr_by_cell_condition_chr)
dev.off()

### check gene expression by condition, chr, and group
summary_expr_by_condition_chr_group <- gene_chr_access %>% 
  group_by(Condition, chr, group) %>% 
  summarise(count = n(), 
            mean_expr = mean(score), 
            median_expr = median(score))
  # ungroup() %>% 
  # complete(Condition, chr, group, fill = list(count = 0, mean_expr = 0, median_expr = 0)) 

pdf(paste0(fig_save_dir, "/Chromosome_accessibility_condition_chr_group_mean.pdf"), width = 8, height = 5.5)
for (chr_index in 1:22){
  data_plot <- summary_expr_by_condition_chr_group[summary_expr_by_condition_chr_group$chr == paste0("chr", chr_index), ]
  p_expr_by_condition_chr_group <- ggplot(data_plot, aes(x = group, y = mean_expr, group = Condition, color = Condition)) + 
    geom_line(position = position_dodge(width = 0.1), size = 1) + 
    geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE) + 
    theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position="right") +
    theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
    labs(x = "Chromosome accessibility levels", y = "Chromosome accessibility (mean)", color = "Condition", title = chr_index)
  print(p_expr_by_condition_chr_group)
}
dev.off()

pdf(paste0(fig_save_dir, "/gene_expr_chr_group_median.pdf"), width = 8, height = 5.5)
for (chr_index in 1:22){
  # chr_index = 12
  data_plot <- summary_expr_by_condition_chr_group[summary_expr_by_condition_chr_group$chr == paste0("chr", chr_index), ]
  p_expr_by_condition_chr_group <- ggplot(data_plot, aes(x = group, y = median_expr, group = Condition, color = Condition)) + 
    geom_line(position = position_dodge(width = 0.1), size = 1) + 
    geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE) + 
    theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position="right") +
    theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
    labs(x = "Gene expression levels", y = "Gene expression (median)", color = "Condition", title = chr_index)
  print(p_expr_by_condition_chr_group)
}
dev.off()
