library(ggplot2)
library(ggpubr)
library(scales)
library(dplyr)
library(Seurat)

setwd("/Users/zhemingan/Documents/BCH_research/Gene_Expression_Analysis")
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])
group_num <- 8
Condition_list <- c("Control", "IHD")

##### read in scRNA-seq data
all_celltype_RNAseq <- readRDS("./data/Seurat.obj_with_annotation.RDS")
CM_cells <- subset(all_celltype_RNAseq, subset = annotated_clusters == "Cardiomyocytes")
## read in gencode gene annotation
# gtf_data <- import("data/chromosome_enrichment/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf", format = "gtf")
# gene_data <- gtf_data[gtf_data$type == "gene", ]
# gene_name_chr_info <- data.frame(gene_name = gene_data$gene_name, chromosome = as.character(seqnames(gene_data)))
# write.csv(gene_name_chr_info, "data/chromosome_enrichment/gene_name_chr_info.csv", row.names = F)
gene_name_chr_info <- read.csv("data/chromosome_enrichment/gene_name_chr_info.csv") %>% rename_with(~c("gene", "chr"))

gene_chr_expr <- c()
for (condition_temp in Condition_list){
  cat("##### Get RNAseq data for:", condition_temp, "...\n")
  expr_level_temp <- data.frame(AverageExpression(CM_cells, group.by = "condition", slot = "data")$RNA) %>% 
    setNames(Condition_list) |> base::`[`(condition_temp) %>% mutate(gene = row.names(.)) %>% 
    setNames(c("average_expr_level", "gene")) %>% mutate(decile = ntile(average_expr_level, n = group_num)) %>% mutate(decile = as.factor(decile))

  gene_chr_expr_temp <- merge(expr_level_temp, gene_name_chr_info, by = "gene") %>% mutate(Condition = condition_temp)
  gene_chr_expr <- rbind(gene_chr_expr, gene_chr_expr_temp)
}

### 
fig_save_dir <- paste0("./results/local_test/chromosome_enrichment/all_genes_RNAseq/")
dir.create(fig_save_dir, recursive = TRUE)

gene_chr_expr <- gene_chr_expr %>% filter(chr %in% paste0("chr", seq(1:22))) %>% 
  mutate(Condition = factor(Condition, levels = Condition_list)) %>% 
  mutate(chr = factor(chr, levels = paste0("chr", seq(1:22))))
  # filter(average_expr_level > 0)

### by condition by chr, calculate the average expression level (median by expression levels) 
summary_expr_by_condition_chr <- gene_chr_expr %>% group_by(Condition, chr) %>% 
  summarise(count = n(), mean_expr = mean(average_expr_level), median_expr = median(average_expr_level))

pdf(paste0(fig_save_dir, "/gene_expr_mean_by_condition_chr.pdf"), width = 8, height = 5.5)
p_expr_by_condition_chr <- ggplot(summary_expr_by_condition_chr, aes(x = chr, y = mean_expr, group = Condition, color = Condition)) + 
  geom_line(position = position_dodge(width = 0.1), size = 1) + geom_point(position = position_dodge(width = 0.1), size = 2) + 
  theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right") + theme_classic() + 
  scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + labs(x = "Chromosome", y = "Gene expression (mean)", color = "Condition", title = "")
print(p_expr_by_condition_chr)
dev.off()

pdf(paste0(fig_save_dir, "/gene_expr_median_by_condition_chr.pdf"), width = 8, height = 5.5)
p_expr_by_condition_chr <- ggplot(summary_expr_by_condition_chr, aes(x = chr, y = median_expr, group = Condition, color = Condition)) + 
  geom_line(position = position_dodge(width = 0.1), size = 1) + geom_point(position = position_dodge(width = 0.1), size = 2) + 
  theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right") + theme_classic() + 
  scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + labs(x = "Chromosome", y = "Gene expression (median)", color = "Condition", title = "")
print(p_expr_by_condition_chr)
dev.off()

### violin plot by condition and then by chr
Control_median_df <- summary_expr_by_condition_chr[summary_expr_by_condition_chr$Condition == "Control", ] |> base::`[`(c("chr", "median_expr"))
gene_chr_expr_normalize <- merge(gene_chr_expr,  Control_median_df, by = "chr") %>% 
  mutate(average_expr_level_nor = round(average_expr_level / median_expr, 3)) %>% 
  mutate(average_expr_level_nor = ifelse(average_expr_level_nor == 0, 0.001, average_expr_level_nor))

wilcox_test_results <- gene_chr_expr_normalize %>% group_by(chr) %>%
  summarize(p_value = wilcox.test(average_expr_level_nor ~ Condition)$p.value)
gene_chr_expr_normalize <- gene_chr_expr_normalize %>% 
  left_join(wilcox_test_results, by = "chr")

p_gene_chr_expr_normalize <- ggplot(gene_chr_expr_normalize, aes(x = Condition, y = average_expr_level_nor, fill = Condition)) + 
  geom_violin(trim = TRUE, na.rm = FALSE) + geom_boxplot(width=0.2, color="white", alpha=0.2) + 
  # stat_compare_means(aes(group = Condition), method = "wilcox.test", label = "p.format", label.y = 1.02 * max(gene_chr_expr_normalize$median_expr)) +
  # stat_compare_means(aes(group = Condition), method = "wilcox.test", label = "p.signif", label.y = 1.06 * max(gene_chr_expr_normalize$median_expr)) +
  # geom_text(data = gene_chr_expr_normalize, aes(x = 1.5, y = max(average_expr_level_nor), label = paste0("p = ", round(p_value, 2))), vjust = 2, size = 8,
  #           nudge_y = 1.5) +
  # geom_segment(data = gene_chr_expr_normalize, aes(x = 1, xend = 2, y = max(average_expr_level_nor) * 1.5, yend = max(average_expr_level_nor) * 1.5), size = 1.2) +
  # geom_segment(data = gene_chr_expr_normalize, aes(x = 1.02, xend = 1.02, y = max(average_expr_level_nor) * 1.5, yend = max(average_expr_level_nor) * 1.0), size = 1.2) +
  # geom_segment(data = gene_chr_expr_normalize, aes(x = 1.98, xend = 1.98, y = max(average_expr_level_nor) * 1.5, yend = max(average_expr_level_nor) * 1.0), size = 1.2) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  # scale_y_continuous(breaks = c(0.01, 1, 100, 10000)^(1/6), labels = function(x) parse(text = x^6)) + 
  # coord_cartesian(ylim = c(0, 10000)) + 
  facet_wrap(~ chr, nrow = 2, scales = "free_x") + 
  theme_classic() + 
  theme(axis.text.x = element_blank(), text = element_text(size = 25), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
        panel.background = element_rect(fill = "white"), legend.position="right") +
  scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
  labs(title = "", x = " ", y = "normalized average gene expression")
p_gene_chr_expr_normalize
ggsave(paste0(fig_save_dir, "violin_gene_chr_expr_normalize.pdf"), p_gene_chr_expr_normalize, width = 20, height = 10,  dpi = 600)
