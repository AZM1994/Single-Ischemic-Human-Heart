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

expr_mutnum_chr_by_cell <- c()
genomic_SCAN2_df <- c()
# condition_temp = Condition_list[1]
for (condition_temp in Condition_list){
  ##### get transcription data
  cat("##### Get transcription data for:", condition_temp, "...\n")
  expr_level_temp <- data.frame(AverageExpression(CM_cells, group.by = "condition", slot = "data")$RNA) %>% 
    setNames(Condition_list) |> base::`[`(condition_temp) %>% 
    mutate(gene = row.names(.)) %>% 
    setNames(c("average_expr_level", "gene")) %>% 
    mutate(decile = ntile(average_expr_level, n = group_num)) %>% mutate(decile = as.factor(decile))
  
  ###########################################################################
  ##################### raw SCAN2 call mutation analysis ####################
  ###########################################################################
  cat("##### Raw SCAN2 call mutation analysis:", condition_temp, "...\n")
  cat("Get genomic context for", condition_temp, "...\n")
  # Case_ID_order <- Hypoxia_PTA_Cases_metadata_collapsed[Hypoxia_PTA_Cases_metadata_collapsed$Condition == condition_temp, "Case_ID"]
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
    rename_with(~ c("Cell_ID", "Case_ID", "Condition", "chr", "start", "end", "ref", "alt", "region", "gene"), .cols = 1:10)
    # filter(age >= 40 & age < 80)
    # mutate(Condition = as.factor(Condition))

  # genic_mutation_temp <- genomic_SCAN2_df_temp
  genic_mutation_temp <- genomic_SCAN2_df_temp[genomic_SCAN2_df_temp$region %in% c("exonic", "exonic;splicing", "intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3"), ] %>%
    mutate(gene = str_remove(gene, "\\(.*\\)$")) %>% filter(!str_detect(gene, ","))
  
  mutation_num_by_cell_temp <- data.frame(table(genic_mutation_temp$gene, genic_mutation_temp$chr, genic_mutation_temp$Condition, genic_mutation_temp$Cell_ID)) %>% 
    setNames(c("gene", "chr", "Condition", "Cell_ID", "mut_number")) %>% 
    group_by(Condition, chr, Cell_ID) %>% 
    summarise("mut_per_chr" = sum(mut_number))
  
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
summary_expr_by_condition_chr <- expr_mutnum_chr_by_cell %>% 
  group_by(Condition, chr) %>% 
  summarise(count = n(), 
            mean_expr = mean(average_expr_level), 
            median_expr = median(average_expr_level), 
            CI_lower = t.test(average_expr_level)$conf.int[1],
            CI_upper = t.test(average_expr_level)$conf.int[2])


pdf(paste0(fig_save_dir, "/gene_expr_median_by_condition_chr.pdf"), width = 8, height = 5.5)
p_expr_by_condition_chr <- ggplot(summary_expr_by_condition_chr, aes(x = chr, y = median_expr, group = Condition, color = Condition)) + 
  # geom_hline(yintercept = 1, color = "black", linewidth = 0.6) + 
  geom_line(position = position_dodge(width = 0.1), size = 1) + 
  geom_point(position = position_dodge(width = 0.1), size = 2) + 
  # geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper)) + 
  # stat_cor(size = 6, show.legend = FALSE) + 
  # geom_errorbar(aes(ymin = mean_ER - sd_ER, ymax = mean_ER + sd_ER), width = 0.2, position = position_dodge(width = 0.1)) +
  # geom_smooth(data = merged_summary_plot, aes(x = decile, y = mean_ER, color = Condition, fill = Condition, group = Condition),
  #             method = "lm", se = TRUE, alpha = 0.2, linewidth = 1, linetype = "dashed") +
  theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right") +
  theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
  # ylim(c(0.5, 1.6)) + 
  labs(x = "Chromosome", y = "Gene expression (median)", color = "Condition", title = "")
print(p_expr_by_condition_chr)
dev.off()


wilcox.test(expr_mutnum_chr_by_cell$average_expr_level[expr_mutnum_chr_by_cell$Condition == "Normal" & expr_mutnum_chr_by_cell$chr == 13],
expr_mutnum_chr_by_cell$average_expr_level[expr_mutnum_chr_by_cell$Condition == "Disease" & expr_mutnum_chr_by_cell$chr == 13])

median(expr_mutnum_chr_by_cell$average_expr_level[expr_mutnum_chr_by_cell$Condition == "Normal" & expr_mutnum_chr_by_cell$chr == 13])
median(expr_mutnum_chr_by_cell$average_expr_level[expr_mutnum_chr_by_cell$Condition == "Disease" & expr_mutnum_chr_by_cell$chr == 13])

summary_expr_by_cell_condition_chr <- expr_mutnum_chr_by_cell %>% 
  group_by(Cell_ID, chr) %>% 
  summarise(count = n(), 
            mean_expr = mean(average_expr_level), 
            median_expr = median(average_expr_level)) %>% 
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

### check gene expression by condition, chr, and decile
summary_expr_by_condition_chr_decile <- expr_mutnum_chr_by_cell %>% 
  group_by(Condition, chr, decile) %>% 
  summarise(count = n(), 
            mean_expr = mean(average_expr_level), 
            median_expr = median(average_expr_level)) %>% 
  ungroup() %>% 
  complete(Condition, chr, decile, fill = list(count = 0, mean_expr = 0, median_expr = 0)) 

pdf(paste0(fig_save_dir, "/gene_expr_chr_decile_mean.pdf"), width = 8, height = 5.5)
for (chr_index in 1:22){
  data_plot <- summary_expr_by_condition_chr_decile[summary_expr_by_condition_chr_decile$chr == chr_index, ]
  p_expr_by_condition_chr_decile <- ggplot(data_plot, aes(x = decile, y = mean_expr, group = Condition, color = Condition)) + 
    geom_line(position = position_dodge(width = 0.1), size = 1) + 
    geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE) + 
    theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position="right") +
    theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
    labs(x = "Gene expression levels", y = "Gene expression (mean)", color = "Condition", title = chr_index)
  print(p_expr_by_condition_chr_decile)
}
dev.off()

pdf(paste0(fig_save_dir, "/gene_expr_chr_decile_median.pdf"), width = 8, height = 5.5)
for (chr_index in 1:22){
  # chr_index = 12
  data_plot <- summary_expr_by_condition_chr_decile[summary_expr_by_condition_chr_decile$chr == chr_index, ]
  p_expr_by_condition_chr_decile <- ggplot(data_plot, aes(x = decile, y = median_expr, group = Condition, color = Condition)) + 
    geom_line(position = position_dodge(width = 0.1), size = 1) + 
    geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE) + 
    theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position="right") +
    theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
    labs(x = "Gene expression levels", y = "Gene expression (median)", color = "Condition", title = chr_index)
  print(p_expr_by_condition_chr_decile)
}
dev.off()




expr_mutnum_summary_disease <- expr_mutnum_chr_by_cell[expr_mutnum_chr_by_cell$Condition == "Disease", ] %>% 
  mutate(Cell_ID = droplevels(Cell_ID))
  # filter(decile %in% c(5,6,7,8))

# haha <- genomic_SCAN2_df[genomic_SCAN2_df$gene %in% hypoxia_gene_list$genes, ] %>% 
#   group_by(Condition, chr) %>% 
#   summarise(count = n())
  
# haha <- expr_mutnum_summary_normal %>% 
#   group_by(decile, chr) %>% 
#   summarise(average_expr_decile = mean(average_expr_level), 
#             average_mut_decile = mean(mut_per_chr)) %>% 
#   mutate(average_expr_decile = round(average_expr_decile, 4))

summary_normal_df <- expr_mutnum_summary_normal %>% 
  group_by(Cell_ID, chr) %>% 
  summarise(average_expr_level = mean(average_expr_level), 
            average_mut_per_chr = mean(mut_per_chr)) %>% 
  mutate(percentage_expr = average_expr_level / sum(average_expr_level)) %>% 
  mutate(percentage_mut = average_mut_per_chr / sum(average_mut_per_chr)) %>% 
  mutate(condition = "Normal") %>% 
  ungroup() %>% 
  complete(Cell_ID, chr, fill = list(percentage_expr = 0, percentage_mut = 0, average_expr_level = 0, average_mut_per_chr = 0, condition = "Normal")) %>%
  arrange(match(Cell_ID, Hypoxia_PTA_Cases_metadata$Cell_ID)) %>% 
  # mutate(estimated_snv_burden = rep(Hypoxia_PTA_Cases_metadata$snv.burden[Hypoxia_PTA_Cases_metadata$Condition == "Normal"], each = 22)) %>%
  mutate(estimated_snv_burden = rep(Hypoxia_PTA_Cases_metadata_age_match$snv.burden[Hypoxia_PTA_Cases_metadata_age_match$Condition == "Normal"], each = 22)) %>%
  mutate(chr_seqlength = rep(seqlengths_list, length(unique(expr_mutnum_summary_normal$Cell_ID))) / (1024*1024*1024)) %>%
  mutate(Chr_burden = percentage_mut * estimated_snv_burden) %>%
  mutate(estimated_chr_burden_density = Chr_burden / chr_seqlength) %>% 
  mutate(expr_over_mut = average_expr_level / average_mut_per_chr)

summary_disease_df <- expr_mutnum_summary_disease %>% 
  group_by(Cell_ID, chr) %>% 
  summarise(average_expr_level = mean(average_expr_level), 
            average_mut_per_chr = mean(mut_per_chr)) %>% 
  mutate(percentage_expr = average_expr_level / sum(average_expr_level)) %>% 
  mutate(percentage_mut = average_mut_per_chr / sum(average_mut_per_chr)) %>% 
  mutate(condition = "Disease") %>% 
  ungroup() %>% 
  complete(Cell_ID, chr, fill = list(percentage_expr = 0, percentage_mut = 0, average_expr_level = 0, average_mut_per_chr = 0, condition = "Disease")) %>%
  arrange(match(Cell_ID, Hypoxia_PTA_Cases_metadata$Cell_ID)) %>% 
  # mutate(estimated_snv_burden = rep(Hypoxia_PTA_Cases_metadata$snv.burden[Hypoxia_PTA_Cases_metadata$Condition == "Disease"], each = 22)) %>%
  mutate(estimated_snv_burden = rep(Hypoxia_PTA_Cases_metadata_age_match$snv.burden[Hypoxia_PTA_Cases_metadata_age_match$Condition == "Disease"], each = 22)) %>%
  mutate(chr_seqlength = rep(seqlengths_list, length(unique(expr_mutnum_summary_disease$Cell_ID))) / (1024*1024*1024)) %>% 
  mutate(Chr_burden = percentage_mut * estimated_snv_burden) %>% 
  mutate(estimated_chr_burden_density = Chr_burden / chr_seqlength) %>% 
  mutate(expr_over_mut = average_expr_level / average_mut_per_chr)

# summary_all_chr_df <- rbind(summary_normal_df, summary_disease_df) |> base::`[`(c("Cell_ID", "chr", "estimated_snv_burden", "condition", "percentage")) %>% 
#   mutate(condition = factor(condition, level = c("Normal", "Disease")))
summary_all_chr_df <- rbind(summary_normal_df, summary_disease_df)

expr_p.value_list <- numeric()
for (chr_number in 1:22){
  # cat(chr_number)
  chr_wilcox.test_result <- wilcox.test(summary_normal_df$average_expr_level[summary_normal_df$chr == chr_number], 
                                        summary_disease_df$average_expr_level[summary_disease_df$chr == chr_number], "two.sided")
  expr_p.value_list <- c(expr_p.value_list, round(chr_wilcox.test_result$p.value, digits = 3))
}

expr_percentage_p.value_list <- numeric()
for (chr_number in 1:22){
  # cat(chr_number)
  chr_wilcox.test_result <- wilcox.test(summary_normal_df$percentage_expr[summary_normal_df$chr == chr_number], 
                                        summary_disease_df$percentage_expr[summary_disease_df$chr == chr_number], "two.sided")
  expr_percentage_p.value_list <- c(expr_percentage_p.value_list, round(chr_wilcox.test_result$p.value, digits = 3))
}


averages_df <- summary_all_chr_df %>%
  group_by(condition, chr) %>%
  summarize(avg_percentage_expr = mean(percentage_expr), 
            avg_percentage_mut = mean(percentage_mut),
            avg_expr_over_mut = mean(!is.na(expr_over_mut)), 
            average_expr_level = mean(average_expr_level), 
            average_mut_per_chr = mean(average_mut_per_chr))
            # CI_lower = t.test(avg_percentage_expr)$conf.int[1],
            # CI_upper = t.test(avg_percentage_expr)$conf.int[2])

averages_normal_df <- averages_df[averages_df$condition == "Normal", ]
averages_disease_df <- averages_df[averages_df$condition == "Disease", ]
merged_average_df <- merge(averages_normal_df, averages_disease_df, by = "chr") %>% 
  arrange(match(chr, 1:22)) %>% 
  mutate(expr.p.value = expr_p.value_list) %>% 
  mutate(percentage.p.value = percentage_p.value_list)
  # mutate(burden_dot_size = 10/log10(abs(average_burden_density.x - average_burden_density.y))) %>% 
  # # mutate(percentage_dot_size = (abs(average_percentage.x - average_percentage.y)))
  # mutate(percentage_dot_size = 2 - 10/log10(6*abs(average_percentage.x - average_percentage.y)))

analysis_ID <- "gene_expr_chr"
fig_save_dir <- paste0("./results/local_test/RNAseq_results/", analysis_ID, "/")
dir.create(fig_save_dir, recursive = TRUE)

scatter_plot_chr_expr <- ggplot(merged_average_df, aes(x = average_expr_level.x, y = average_expr_level.y, color = expr.p.value)) + 
  geom_segment(aes(x = 0, xend = 5, y = 0, yend = 5), colour = "grey20", linetype = "dashed", size = 1.5) + 
  geom_point(size = 10) + 
  # geom_errorbar(aes(ymin = CI_lower.y, ymax = CI_upper.y)) + 
  # geom_errorbarh(aes(xmin = CI_lower.x, xmax = CI_upper.x)) + 
  geom_text(aes(label = chr), vjust = 0.5, color = "white", size = 5) + 
  annotate("text", size = 15, x = 1.8, y = 4, label = "*", hjust = 0) + 
  annotate("text", size = 15, x = 0.3, y = 1, label = "*", hjust = 0) + 
  scale_x_continuous("avg. expression level in Control") + 
  scale_y_continuous("avg. expression level in IHD") + 
  theme_classic() + 
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), 
        panel.background = element_rect(fill = "white"), legend.position = c(0.9,0.3)) + 
  scale_color_gradient(low = "red", high = "blue", limits = c(0.01, 1), breaks = c(0.01, 0.05, 1), name = "pvalue", transform = "log10")
print(scatter_plot_chr_expr)
ggsave(paste0(fig_save_dir, "/gene_expr_chr.png"),
       plot = scatter_plot_chr_expr, width = 9, height = 8, dpi = 600)

scatter_plot_chr_expr_percentage <- ggplot(merged_average_df, aes(x = avg_percentage_expr.x, y = avg_percentage_expr.y, color = percentage.p.value)) + 
  geom_segment(aes(x = 0, xend = 0.1, y = 0, yend = 0.1), colour = "grey20", linetype = "dashed", size = 1.5) + 
  geom_point(size = 10) + 
  # geom_errorbar(aes(ymin = CI_lower.y, ymax = CI_upper.y)) + 
  # geom_errorbarh(aes(xmin = CI_lower.x, xmax = CI_upper.x)) + 
  geom_text(aes(label = chr), vjust = 0.5, color = "white", size = 5) + 
  annotate("text", size = 15, x = 0.011, y = 0.025, label = "*", hjust = 0) + 
  annotate("text", size = 15, x = 0.052, y = 0.10, label = "*", hjust = 0) + 
  scale_x_continuous("avg. expression level in Control") + 
  scale_y_continuous("avg. expression level in IHD") + 
  theme_classic() + 
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), 
        panel.background = element_rect(fill = "white"), legend.position = c(0.9,0.3)) + 
  scale_color_gradient(low = "red", high = "blue", limits = c(0.01, 1), breaks = c(0.01, 0.05, 1), name = "pvalue", transform = "log10")
print(scatter_plot_chr_expr_percentage)
ggsave(paste0(fig_save_dir, "/gene_expr_percentage_chr.png"),
       plot = scatter_plot_chr_expr_percentage, width = 9, height = 8, dpi = 600)

