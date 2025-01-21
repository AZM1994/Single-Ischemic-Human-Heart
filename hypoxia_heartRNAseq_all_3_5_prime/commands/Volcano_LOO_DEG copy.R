library(Seurat)
library(ggplot2)
library(tidyverse) 
library(RColorBrewer) 
library(ggrepel)
library(plotly)
library(Rmisc)
library(patchwork)

wd_dir <- getwd()
result_dir <- paste0(wd_dir, "/results/LOO_DEG_Cardiomyocyte")
dir.create(result_dir)

### read DEG data & exclude sex genes
control.disease.combined_CM_only_CM_only  <- readRDS(paste0(wd_dir, "/results/Seurat.obj_sub_clustering_CM_only.RDS"))
sex_gene_file <- read.csv(paste0(wd_dir, "/data/GenesonChrXY.csv"), header = TRUE)
sex_gene_list <- sex_gene_file$SYMBOL
file.remove(paste(result_dir, '/', 'no_lb_top_up_and_down_reg_gene_FC_pV.csv', sep = ""))
file.remove(paste(result_dir, '/', 'with_lb_top_up_and_down_reg_gene_FC_pV.csv', sep = ""))

DEG_data <- read.csv(paste0(result_dir, "/marker_genes_CM_only_LOO.csv")) %>% 
  select(-X) %>% 
  filter(!gene %in% sex_gene_list) %>% 
  mutate(AbsLogFC = abs(avg_log2FC)) %>%  group_by(gene) %>%
  slice_max(AbsLogFC, with_ties = FALSE) %>%
  ungroup() %>% select(-AbsLogFC)

### set astronomically small pvalues to random values that close to the smallest one
pct_bound <- 0.1
FC_bound <- 0.25
P_bound <- 0.05
mean <- 10 # Mean of the distribution
sd <- 5
DEG_data$p_val_adj[DEG_data$p_val_adj == 0] <- min(DEG_data$p_val_adj[DEG_data$p_val_adj != 0])

DEG_data_plot <- DEG_data
# DEG_data_plot <- DEG_data_plot[DEG_data_plot$pct.1 >= pct_bound & DEG_data_plot$pct.2 >= pct_bound, ]
DEG_data_plot$p_val_adj_plot <- DEG_data_plot$p_val_adj
DEG_data_plot$p_val_adj_log10 <- -log10(DEG_data_plot$p_val_adj_plot)
length_of_replacement <- length(DEG_data_plot$p_val_adj_plot[DEG_data_plot$p_val_adj_plot == min(DEG_data_plot$p_val_adj_plot)])
# random_replacement <- rnorm(length_of_replacement, mean = mean, sd = sd)
random_replacement <- runif(length_of_replacement, min = 0, max = 5)
random_replacement_positive <- random_replacement[random_replacement > 0]
while (length(random_replacement_positive) < length_of_replacement) {
  additional_numbers <- rnorm(length_of_replacement - length(random_replacement_positive), mean = mean, sd = sd)
  random_replacement_positive <- c(random_replacement_positive, additional_numbers[additional_numbers > 0])
}
# Trim the excess numbers if the length exceeds length_of_replacement
random_replacement_positive <- head(random_replacement_positive, length_of_replacement)
random_replacement_positive <- random_replacement_positive

DEG_data_plot$p_val_adj_log10[DEG_data_plot$p_val_adj_plot == min(DEG_data_plot$p_val_adj_plot)] <- (-log10(min(DEG_data_plot$p_val_adj_plot))) - random_replacement

log2FoldChange_plot <- DEG_data_plot$avg_log2FC
pvalue_plot <- DEG_data_plot$p_val_adj
pvalue_temp_plot <- pvalue_plot[!is.na(pvalue_plot)]

DEG_data_plot$diffexpressed_plot <- 1
DEG_data_plot$log2FoldChange_plot <- log2FoldChange_plot
DEG_data_plot$pvalue_plot <- pvalue_plot
DEG_data_plot$diffexpressed_plot[DEG_data_plot$log2FoldChange_plot > FC_bound & DEG_data_plot$pvalue_plot < P_bound] <- "UP"
DEG_data_plot$diffexpressed_plot[DEG_data_plot$log2FoldChange_plot < -FC_bound & DEG_data_plot$pvalue_plot < P_bound] <- "DOWN"
# sig_label_plot <- DEG_data_plot$X[DEG_data_plot$diffexpressed_plot != "NA"]
# DEG_data_plot$diffexpressed_plot[DEG_data_plot$pct.1 < pct_bound | DEG_data_plot$pct.2 < pct_bound] <- "low_pct"
# DEG_data_plot$diffexpressed_plot <- factor(x = DEG_data_plot$diffexpressed_plot, levels = c("UP", "DOWN", "NA", "low_pct"))

### find threshold
DEG_data <- DEG_data[DEG_data$pct.1 >= pct_bound & DEG_data$pct.2 >= pct_bound, ]
log2FoldChange <- DEG_data$avg_log2FC
pvalue <- DEG_data$p_val_adj
pvalue_temp <- pvalue[!is.na(pvalue)]

DEG_data$log2FoldChange <- log2FoldChange
DEG_data$pvalue <- pvalue

DEG_data$diffexpressed <- 1
DEG_data$diffexpressed[DEG_data$log2FoldChange > FC_bound & DEG_data$pvalue < P_bound] <- "UP"
DEG_data$diffexpressed[DEG_data$log2FoldChange < -FC_bound & DEG_data$pvalue < P_bound] <- "DOWN"
sig_label <- DEG_data$gene[DEG_data$diffexpressed != "NA"]

####### select top 10 genes
max_FC <- max(range(log2FoldChange))
min_FC <- min(range(log2FoldChange))
max_p <- max(range(-log10(pvalue_temp)))
min_p <- min(range(-log10(pvalue_temp)))

pvalue_weight <- 1
FC_weight <- 1
up_top_ten_bound <- sort(pvalue_weight * (-log10(DEG_data$p_val_adj[DEG_data$avg_log2FC >= 0]) / max_p) + 
                           FC_weight * (DEG_data$avg_log2FC[DEG_data$avg_log2FC >= 0] / max_FC), decreasing = TRUE)[10]
down_top_ten_bound <- sort(pvalue_weight * (-log10(DEG_data$p_val_adj)[DEG_data$avg_log2FC < 0] / max_p) + 
                             FC_weight * (DEG_data$avg_log2FC[DEG_data$avg_log2FC < 0] / min_FC), decreasing = TRUE)[10]
DEG_data$up_topten <- 1
DEG_data$up_topten[pvalue_weight * (-log10(DEG_data$p_val_adj) / max_p) + 
                     FC_weight * (DEG_data$avg_log2FC / max_FC) >= up_top_ten_bound & DEG_data$avg_log2FC >= 0] <- "up_ten"
DEG_data$down_topten <- 1
DEG_data$down_topten[pvalue_weight * (-log10(DEG_data$p_val_adj) / max_p) + 
                       FC_weight * (DEG_data$avg_log2FC / min_FC) >= down_top_ten_bound & DEG_data$avg_log2FC < 0] <- "down_ten"
top_ten_label_mix <- DEG_data$gene[DEG_data$up_topten == "up_ten" | DEG_data$down_topten == "down_ten"]

########## save top 10 up and down regulated genes (ordered by avg_log2FC)
top_upreg_gene_FC_pV <- DEG_data[c("gene", "avg_log2FC", "p_val_adj")][DEG_data$up_topten == "up_ten", ]
top_upreg_gene_FC_pV <- top_upreg_gene_FC_pV[complete.cases(top_upreg_gene_FC_pV), ]
top_upreg_gene_FC_pV <- cbind(top_upreg_gene_FC_pV, "UP")
colnames(top_upreg_gene_FC_pV) <- c("gene", "avg_log2FC", "p_val_adj", "regulation")
top_upreg_gene_FC_pV <- top_upreg_gene_FC_pV[order(top_upreg_gene_FC_pV$avg_log2FC, decreasing = TRUE), ]
top_upreg_gene_name_list <- top_upreg_gene_FC_pV$gene

top_downreg_gene_FC_pV <- DEG_data[c("gene", "avg_log2FC", "p_val_adj")][DEG_data$down_topten == "down_ten", ]
top_downreg_gene_FC_pV <- top_downreg_gene_FC_pV[complete.cases(top_downreg_gene_FC_pV), ]
top_downreg_gene_FC_pV <- cbind(top_downreg_gene_FC_pV, "DOWN")
colnames(top_downreg_gene_FC_pV) <- c("gene", "avg_log2FC", "p_val_adj", "regulation")
top_downreg_gene_FC_pV <- top_downreg_gene_FC_pV[order(top_downreg_gene_FC_pV$avg_log2FC, decreasing = FALSE), ]
top_downreg_gene_name_list <- top_downreg_gene_FC_pV$gene

top_up_and_down_reg_gene_FC_pV <- rbind(top_upreg_gene_FC_pV, top_downreg_gene_FC_pV)
write.table(top_up_and_down_reg_gene_FC_pV, paste(result_dir, '/', '_top_up_and_down_reg_gene_FC_pV.csv', sep = ""),
            append = TRUE, row.names=FALSE, col.names = FALSE, sep=",")

##### feature plots show how top genes expressed in each sample
# pdf(paste0(result_dir, "/DEGs/", plot_type, "/feature_plot/top10_Upregulation_feature_plot_",
#            "by_sample_condition_", file_index,"_", plot_type, ".pdf"), width = 21, height = 15)
# print(FeaturePlot(control.disease.combined_CM_only, features = top_upreg_gene_name_list[1:2], split.by = "donor_condition",
#                   ncol = 4, max.cutoff = NA, cols = c("grey", "red")))
# dev.off()
# 
# pdf(paste0(result_dir, "/DEGs/", plot_type, "/feature_plot/top10_Upregulation_feature_plot_",
#            "by_condition_", file_index,"_", plot_type, ".pdf"), width = 10, height = 25)
# print(FeaturePlot(control.disease.combined_CM_only, features = top_upreg_gene_name_list[1:2], split.by = "condition",
#                   ncol = 4, max.cutoff = NA, cols = c("grey", "red")))
# dev.off()
# 
# pdf(paste0(result_dir, "/DEGs/", plot_type, "/feature_plot/top10_Downregulation_feature_plot_",
#            "by_sample_condition_", file_index,"_", plot_type, ".pdf"), width = 21, height = 15)
# print(FeaturePlot(control.disease.combined_CM_only, features = top_downreg_gene_name_list[1:2], split.by = "donor_condition",
#                   ncol = 4, max.cutoff = NA, cols = c("grey", "red")))
# dev.off()
# 
# pdf(paste0(result_dir, "/DEGs/", plot_type, "/feature_plot/top10_Downregulation_feature_plot_",
#            "by_condition_", file_index,"_", plot_type, ".pdf"), width = 10, height = 25)
# print(FeaturePlot(control.disease.combined_CM_only, features = top_downreg_gene_name_list[1:2], split.by = "condition",
#                   ncol = 4, max.cutoff = NA, cols = c("grey", "red")))
# dev.off()

##### generate violin plot 
control.disease.combined_CM_only$donor_condition <- factor(x = control.disease.combined_CM_only$donor_condition,
                                                     levels = c("936_Ctrl", "5087_Ctrl", "5919_Ctrl", "1156_IHD", "5111_IHD", "5874_IHD"))
Idents(control.disease.combined_CM_only) <- "donor_condition"
# idents_list <- c(paste(file_index, "_Control", sep = ""), paste(file_index, "_Disease", sep = ""))
# pdf(paste0(result_dir, "/DEGs/", plot_type, "/violin_plot/top10_Upreg_violin_plot_by_condition_", file_index,"_", plot_type, ".pdf"), width = 20, height = 10)
# p_v01 <- VlnPlot(control.disease.combined_CM_only, features = top_upreg_gene_name_list[1:2], group.by = "donor_condition", ncol = 5)
# print(p_v01)
# dev.off()
# pdf(paste0(result_dir, "/DEGs/", plot_type, "/violin_plot/top10_Upreg_violin_plot_by_samples.condition_", file_index,"_", plot_type, ".pdf"), width = 30, height = 20)
# p_v02 <- VlnPlot(control.disease.combined_CM_only, features = top_upreg_gene_name_list[1:10], idents = idents_list, group.by = "samples.condition", ncol = 2)
# print(p_v02)
# dev.off()
# pdf(paste0(result_dir, "/DEGs/", plot_type, "/violin_plot/top10_Downreg_violin_plot_by_condition_", file_index,"_", plot_type, ".pdf"), width = 20, height = 10)
# p_v03 <- VlnPlot(control.disease.combined_CM_only, features = top_downreg_gene_name_list[1:10], idents = idents_list, group.by = "condition", ncol = 5)
# print(p_v03)
# dev.off()
# pdf(paste0(result_dir, "/DEGs/", plot_type, "/violin_plot/top10_Downreg_violin_plot_by_samples.condition_", file_index,"_", plot_type, ".pdf"), width = 30, height = 20)
# p_v04 <- VlnPlot(control.disease.combined_CM_only, features = top_downreg_gene_name_list[1:10], idents = idents_list, group.by = "samples.condition", ncol = 2)
# print(p_v04)
# dev.off()

table(DEG_data$diffexpressed)

is_found <- DEG_data_plot$gene %in% c(top_upreg_gene_name_list, top_downreg_gene_name_list)
top_ten_label_mix_plot <- DEG_data_plot$gene
top_ten_label_mix_plot[!is_found] <- NA
DEG_data_plot$diffexpressed_plot <- factor(x = DEG_data_plot$diffexpressed_plot, levels = c("UP", "DOWN", "NA"))
# plot adding up all layers we have seen so far
pdf(paste0(figure_save_path, 'volcano/Volcano_plot_cluster_',file_index,'_', plot_type,'.pdf', sep = ""), width = 10, height = 8)
p <- ggplot(data = DEG_data_plot, aes(x = avg_log2FC, y = p_val_adj_log10, color = diffexpressed_plot, label = top_ten_label_mix_plot)) +
  geom_point(size = 1.5, na.rm = FALSE) + geom_text_repel(max.overlaps = 30, na.rm = FALSE, color = "black") +
  scale_color_manual(values = c("#00AFBB", "#bb0c00", "grey"), labels = c("Upregulated", "Downregulated", "Not significant")) +
  geom_vline(xintercept = c(-FC_bound, FC_bound), col="grey50", linetype = 'dashed') +
  geom_hline(yintercept = -log10(P_bound), col = "grey50", linetype = 'dashed') +
  labs(color = '', x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  xlim(-max_FC - 0.5, max_FC + 0.5) +
  # xlim(-4, 4) + 
  ylim(0, max_p) + 
  theme_classic() + 
  theme(text = element_text(size=24), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), 
        axis.ticks.y.right = element_blank(), panel.background = element_rect(fill = "white"), legend.position="right")
plot(p)
dev.off()
  # }
# }

## top up and down CM DEGs expression in each sample
# control.disease.combined_CM_only  <- readRDS(paste0(result_dir, "/Seurat.obj_with_annotation.RDS"))
head(control.disease.combined_CM_only)
# top10_up_CM_DGEs <- c("NR4A1", "PDK4", "CRYAB", "CORIN", "RORA", "ARID5B", "ARL15", "ENOX1", "MALAT1", "ATP2A2")
top10_up_CM_DGEs <- c("PDK4", "CRYAB", "CORIN", "ARL15")
DoHeatmap(control.disease.combined_CM_only, features = top10_up_CM_DGEs, group.by = "samples.condition")

top10_down_CM_DGEs <- c("PRKAG2", "GALNT17", "DGKG", "FILIP1L", "TNNI3", "MYBPC3", "SCN5A", "INPP4B", "OBSCN", "PLCG2")
DoHeatmap(control.disease.combined_CM_only, features = top10_down_CM_DGEs, group.by = "samples.condition", size = 3)



