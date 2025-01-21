library(tidyr)
library(ggplot2)
library(ggpubr)
# library(scales)
library(dplyr)
library(Seurat)

setwd("/Users/zhemingan/Documents/BCH_research/Gene_Expression_Analysis")
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])
group_num <- 8
Condition_list <- c("Control", "IHD")

##### read in scRNA-ATACseq data
atac_normal <- c()
atac_disease <- c()
total_reads_4638 <- 41195236
total_reads_5828 <- 9955703
total_reads_5919 <- 33130651
total_reads_604 <- 17569221

# total_reads_normal <- (total_reads_4638 + total_reads_5828 + total_reads_5919) / 3
total_reads_normal <- total_reads_5919
total_reads_disease <- total_reads_604
reads_ratio <- total_reads_disease / total_reads_normal

atac_5919_all <- read.table("./data/scATACseq_preprocessing/ATAC_cov_average_hg19_1kb_bin.5919_all.bed") %>% setNames(c("chr", "start", "end", "id", "score_5919_all"))
atac_normal <- atac_5919_all %>% 
  setNames(c("chr", "start", "end", "id", "score")) %>% 
  mutate(group = ntile(score, n = group_num)) %>% 
  mutate(group = as.factor(group)) %>% 
  mutate(Condition = "Control")

atac_604_all <- read.table("./data/scATACseq_preprocessing/ATAC_cov_average_hg19_1kb_bin.604_all.bed") %>% setNames(c("chr", "start", "end", "id", "score_604_all"))
atac_disease <- atac_604_all %>% 
  setNames(c("chr", "start", "end", "id", "score")) %>% 
  mutate(score = score / reads_ratio) %>% 
  mutate(group = ntile(score, n = group_num)) %>% 
  mutate(group = as.factor(group)) %>% 
  mutate(Condition = "IHD")

gene_chr_access <- rbind(atac_normal, atac_disease)

### 
fig_save_dir <- paste0("./results/local_test/chromosome_enrichment/all_genes_ATACseq/")
dir.create(fig_save_dir, recursive = TRUE)

gene_chr_access <- gene_chr_access %>% 
  # filter(score > 0) %>% 
  group_by(Condition, chr) %>% 
  mutate(value_group = ceiling(row_number() / 100)) %>% 
  group_by(Condition, chr, value_group) %>%
  summarize(score = mean(score, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(Condition = factor(Condition, levels = Condition_list)) %>% 
  mutate(chr = factor(chr, levels = paste0("chr", seq(1:22))))
  # filter(score > 0)

# gene_chr_access <- gene_chr_access %>% 
#   mutate(Condition = factor(Condition, levels = Condition_list)) %>% 
#   mutate(chr = factor(chr, levels = paste0("chr", seq(1:22))))
  # filter(score > 0)

### by condition by chr, calculate the average accessibility level (median by expression levels) 
summary_access_by_condition_chr <- gene_chr_access %>% group_by(Condition, chr) %>% 
  summarise(count = n(), mean_access = mean(score), median_access = median(score))

pdf(paste0(fig_save_dir, "/Chromosome_accessibility_mean_by_condition_chr.pdf"), width = 8, height = 5.5)
p_access_by_condition_chr <- ggplot(summary_access_by_condition_chr, aes(x = chr, y = mean_access, group = Condition, color = Condition)) +
  geom_line(position = position_dodge(width = 0.1), size = 1) + geom_point(position = position_dodge(width = 0.1), size = 2) +
  theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
        panel.background = element_rect(fill = "white"), legend.position="right") + theme_classic() +
  scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + labs(x = "Chromosome", y = "Chromosome accessibility (mean)", color = "Condition", title = "")
print(p_access_by_condition_chr)
dev.off()

pdf(paste0(fig_save_dir, "/Chromosome_accessibility_median_by_condition_chr.pdf"), width = 8, height = 5.5)
p_access_by_condition_chr <- ggplot(summary_access_by_condition_chr, aes(x = chr, y = median_access, group = Condition, color = Condition)) +
  geom_line(position = position_dodge(width = 0.1), size = 1) + geom_point(position = position_dodge(width = 0.1), size = 2) +
  theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
        panel.background = element_rect(fill = "white"), legend.position="right") + theme_classic() +
  scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + labs(x = "Chromosome", y = "Chromosome accessibility (median)", color = "Condition", title = "")
print(p_access_by_condition_chr)
dev.off()

### violin plot by condition and then by chr
Control_median_df <- summary_access_by_condition_chr[summary_access_by_condition_chr$Condition == "Control", ] |> base::`[`(c("chr", "median_access"))
gene_chr_access_normalize <- merge(gene_chr_access,  Control_median_df, by = "chr") %>% mutate(score_nor = score / median_access)

t_test_results <- gene_chr_access_normalize %>% group_by(chr) %>%
  summarize(p_value = wilcox.test(score_nor ~ Condition)$p.value)
gene_chr_access_normalize <- gene_chr_access_normalize %>% left_join(t_test_results, by = "chr") |> base::`[`(c("chr", "score_nor", "Condition", "p_value"))

# b1= gene_chr_access_normalize$score_nor[gene_chr_access_normalize$Condition == "Control" & gene_chr_access_normalize$chr == "chr1"]
# b2= gene_chr_access_normalize$score_nor[gene_chr_access_normalize$Condition == "IHD" & gene_chr_access_normalize$chr == "chr11"]
# mean(b1)
# t.test(b1,b2)

p_gene_chr_access_normalize <- ggplot(gene_chr_access_normalize, aes(x = Condition, y = score_nor, fill = Condition)) + 
  geom_violin(trim = TRUE, na.rm = FALSE) + geom_boxplot(width=0.2, color="white", alpha=0.2) + 
  geom_text(data = gene_chr_access_normalize, aes(x = 1.5, y = max(score_nor), label = paste0("p = ", round(p_value, 2))), vjust = 2, size = 8, 
            nudge_y = 1.5) + 
  geom_segment(data = gene_chr_access_normalize, aes(x = 1, xend = 2, y = max(score_nor) * 1.5, yend = max(score_nor) * 1.5), size = 1.2) + 
  geom_segment(data = gene_chr_access_normalize, aes(x = 1.02, xend = 1.02, y = max(score_nor) * 1.5, yend = max(score_nor) * 1.0), size = 1.2) +
  geom_segment(data = gene_chr_access_normalize, aes(x = 1.98, xend = 1.98, y = max(score_nor) * 1.5, yend = max(score_nor) * 1.0), size = 1.2) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  # scale_y_log10(breaks = c(0.01, 1, 100, 10000), labels = function(x) parse(text = paste0("10^", log10(x)))) +
  scale_y_log10() + 
  # coord_cartesian(ylim = c(0, 10000)) + 
  facet_wrap(~ chr, nrow = 2, scales = "free_x") + 
  theme_classic() + 
  theme(axis.text.x = element_blank(), text = element_text(size = 25), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
        panel.background = element_rect(fill = "white"), legend.position="right") +
  scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
  labs(title = "", x = " ", y = "normalized average gene expression")
p_gene_chr_access_normalize
ggsave(paste0(fig_save_dir, "violin_gene_chr_access_normalize.pdf"), p_gene_chr_access_normalize, width = 20, height = 10)
# ggsave(paste0(fig_save_dir, "violin_gene_chr_access_normalize.pdf"), p_gene_chr_access_normalize, width = 20, height = 10,  dpi = 50)
# ggsave(paste0(fig_save_dir, "violin_gene_chr_access_normalize.png"), p_gene_chr_access_normalize, width = 20, height = 10)


