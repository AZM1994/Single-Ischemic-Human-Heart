library(Seurat)
library(ggplot2)
# library(tidyverse) 
library(RColorBrewer) 
library(ggrepel)
library(plotly)
library(Rmisc)
library(patchwork)
library(dplyr)

wd_dir <- getwd()
LOO_DEG_Volcano_dir <- paste0(wd_dir, "/results_IHD/5-LOO_DEG_Volcano")
dir.create(LOO_DEG_Volcano_dir, recursive = TRUE, showWarnings = FALSE)

## read DEG data & exclude sex genes
# control.disease.combined_CM_only  <- readRDS(paste0(wd_dir, "/results_IHD/Seurat.obj_sub_clustering_CM_only.RDS"))
sex_gene_file <- read.csv(paste0(wd_dir, "/data/GenesonChrXY.csv"), header = TRUE)
sex_gene_list <- sex_gene_file$SYMBOL
file.remove(paste(LOO_DEG_Volcano_dir, '/', 'no_lb_top_up_and_down_reg_gene_FC_pV.csv', sep = ""))
file.remove(paste(LOO_DEG_Volcano_dir, '/', 'with_lb_top_up_and_down_reg_gene_FC_pV.csv', sep = ""))

pct_bound <- 0.15
FC_bound <- 0.50
P_bound <- 0.05

DEG_data <- read.csv(paste0(LOO_DEG_Volcano_dir, "/LOO_DEG_CM_only_no_filter.csv")) %>% 
  select(-X) %>% filter(!Gene %in% sex_gene_list) %>% 
  mutate(AbsLogFC = abs(avg_log2FC)) %>% group_by(Gene) %>% 
  slice_max(AbsLogFC, with_ties = FALSE) %>%
  ungroup() %>% select(-AbsLogFC) %>% 
  filter(pct.1 >= pct_bound | pct.2 >= pct_bound)

DEG_df <- read.csv(paste0(LOO_DEG_Volcano_dir, "/LOO_DEG_filtered.csv"))
DEG_up <- DEG_df$Gene[DEG_df$regulation == "up"]
DEG_down <- DEG_df$Gene[DEG_df$regulation == "down"]

## set astronomically small pvalues to random values that close to the smallest one
DEG_data$p_val_adj[DEG_data$p_val_adj == 0] <- min(DEG_data$p_val_adj[DEG_data$p_val_adj != 0])
length_of_replacement <- length(DEG_data$p_val_adj[DEG_data$p_val_adj == min(DEG_data$p_val_adj)])
random_replacement <- runif(length_of_replacement, min = 0, max = 5)
DEG_data_plot <- DEG_data %>% mutate(p_val_adj_log10 = -log10(p_val_adj))
DEG_data_plot$p_val_adj_log10[DEG_data_plot$p_val_adj == min(DEG_data_plot$p_val_adj)] <- (-log10(min(DEG_data_plot$p_val_adj))) - random_replacement
DEG_data_plot <- DEG_data_plot %>% mutate(diffexpressed_plot = ifelse(avg_log2FC > FC_bound & p_val_adj < P_bound, "UP", ifelse(avg_log2FC < -FC_bound & p_val_adj < P_bound, "DOWN", "NA")))

## find threshold
DEG_data <- DEG_data %>% mutate(diffexpressed = ifelse(avg_log2FC > FC_bound & p_val_adj < P_bound, "UP", ifelse(avg_log2FC < -FC_bound & p_val_adj < P_bound, "DOWN", "NA")))
sig_label <- DEG_data$Gene[DEG_data$diffexpressed != "NA"]

####### select top 10 genes
max_FC <- max(range(DEG_data$avg_log2FC))
min_FC <- min(range(DEG_data$avg_log2FC))
max_p <- max(range(-log10(DEG_data$p_val_adj)))
min_p <- min(range(-log10(DEG_data$p_val_adj)))
DEG_data <- DEG_data %>% filter(Gene %in% DEG_df$Gene)

pvalue_weight <- 2
FC_weight <- 1
up_top_ten_bound <- sort(pvalue_weight * (-log10(DEG_data$p_val_adj[DEG_data$avg_log2FC >= 0]) / max_p) + 
                           FC_weight * (DEG_data$avg_log2FC[DEG_data$avg_log2FC >= 0] / max_FC), decreasing = TRUE)[20]
down_top_ten_bound <- sort(pvalue_weight * (-log10(DEG_data$p_val_adj)[DEG_data$avg_log2FC < 0] / max_p) + 
                             FC_weight * (DEG_data$avg_log2FC[DEG_data$avg_log2FC < 0] / min_FC), decreasing = TRUE)[20]

DEG_data$up_topten <- NA
DEG_data$up_topten[pvalue_weight * (-log10(DEG_data$p_val_adj) / max_p) + 
                     FC_weight * (DEG_data$avg_log2FC / max_FC) >= up_top_ten_bound & DEG_data$avg_log2FC >= 0] <- "up_ten"
DEG_data$down_topten <- NA
DEG_data$down_topten[pvalue_weight * (-log10(DEG_data$p_val_adj) / max_p) + 
                       FC_weight * (DEG_data$avg_log2FC / min_FC) >= down_top_ten_bound & DEG_data$avg_log2FC < 0] <- "down_ten"

top_ten_label_mix <- DEG_data$Gene[DEG_data$up_topten == "up_ten" | DEG_data$down_topten == "down_ten"]

########## save top 10 up and down regulated genes (ordered by avg_log2FC)
top_upreg_gene_FC_pV <- DEG_data[c("Gene", "avg_log2FC", "p_val_adj")][DEG_data$up_topten == "up_ten", ]
top_upreg_gene_FC_pV <- top_upreg_gene_FC_pV[complete.cases(top_upreg_gene_FC_pV), ]
top_upreg_gene_FC_pV <- cbind(top_upreg_gene_FC_pV, "UP")
colnames(top_upreg_gene_FC_pV) <- c("Gene", "avg_log2FC", "p_val_adj", "regulation")
top_upreg_gene_FC_pV <- top_upreg_gene_FC_pV[order(top_upreg_gene_FC_pV$avg_log2FC, decreasing = TRUE), ]
top_upreg_gene_name_list <- top_upreg_gene_FC_pV$Gene

top_downreg_gene_FC_pV <- DEG_data[c("Gene", "avg_log2FC", "p_val_adj")][DEG_data$down_topten == "down_ten", ]
top_downreg_gene_FC_pV <- top_downreg_gene_FC_pV[complete.cases(top_downreg_gene_FC_pV), ]
top_downreg_gene_FC_pV <- cbind(top_downreg_gene_FC_pV, "DOWN")
colnames(top_downreg_gene_FC_pV) <- c("Gene", "avg_log2FC", "p_val_adj", "regulation")
top_downreg_gene_FC_pV <- top_downreg_gene_FC_pV[order(top_downreg_gene_FC_pV$avg_log2FC, decreasing = FALSE), ]
top_downreg_gene_name_list <- top_downreg_gene_FC_pV$Gene

top_up_and_down_reg_gene_FC_pV <- rbind(top_upreg_gene_FC_pV, top_downreg_gene_FC_pV)
write.table(top_up_and_down_reg_gene_FC_pV, paste(LOO_DEG_Volcano_dir, '/top_LOO_DEG.csv', sep = ""), append = F, row.names = FALSE, col.names = TRUE, sep = ",")

##### feature plots show how top genes expressed in each sample
# pdf(paste0(LOO_DEG_Volcano_dir, "/DEGs/", plot_type, "/feature_plot/top10_Upregulation_feature_plot_",
#            "by_sample_condition_", file_index,"_", plot_type, ".pdf"), width = 21, height = 15)
# print(FeaturePlot(control.disease.combined_CM_only, features = top_upreg_gene_name_list[1:2], split.by = "donor_condition",
#                   ncol = 4, max.cutoff = NA, cols = c("grey", "red")))
# dev.off()

##### generate violin plot 
# control.disease.combined_CM_only$donor_condition <- factor(x = control.disease.combined_CM_only$donor_condition, levels = c("936_Ctrl", "5087_Ctrl", "5919_Ctrl", "1156_IHD", "5111_IHD", "5874_IHD"))
# Idents(control.disease.combined_CM_only) <- "donor_condition"
# pdf(paste0(LOO_DEG_Volcano_dir, "/DEGs/", plot_type, "/violin_plot/top10_Upreg_violin_plot_by_condition_", file_index,"_", plot_type, ".pdf"), width = 20, height = 10)
  # VlnPlot(control.disease.combined_CM_only, features = top_upreg_gene_name_list[1:2], group.by = "donor_condition", ncol = 5)
# dev.off()

table(DEG_data$diffexpressed)
table(DEG_data_plot$diffexpressed_plot)

# is_found <- DEG_data_plot$Gene %in% c(top_upreg_gene_name_list, top_downreg_gene_name_list)
# is_found <- DEG_data_plot$Gene %in% c("EGR1", "ZFP36", "AC023494.1", "CCBE1", "FOS", "BTG2", "DOK6", "NR4A1", "ZNF608",
#                                       "ANKRD2", "SLC35F1", "SYNDIG1", "KCNAB2", "KIF26B", "TMEM120B", "SCN5A")

is_found <- DEG_data_plot$Gene %in% c("EGR1", "ZFP36", "AC023494.1", "CCBE1", "FOS", "BTG2", "DOK6", "NR4A1", "ZNF608", 
                                      "DST", "MYBPC3", "TNNI3", "OBSCN", "LINC01428", "SCN5A", "MYO18A", "PSD3", "ANKRD2", 
                                      "UBA6-AS1", "ZNF106", "DPF3", "SYNDIG1", "SCN5A", "SRSF6", "DDX5", "MET")
# is_found <- DEG_data_plot$Gene %in% DEG_df$Gene

top_ten_label_mix_plot <- DEG_data_plot$Gene
top_ten_label_mix_plot[!is_found] <- NA
DEG_data_plot$diffexpressed_plot <- factor(x = DEG_data_plot$diffexpressed_plot, levels = c("UP", "DOWN", "NA"))

DEG_color <- c("#ff7f11", "#2eb086")
x_max <- max(abs(min_FC), abs(max_FC))
DEG_vocalno <- ggplot(data = DEG_data_plot, aes(x = avg_log2FC, y = p_val_adj_log10, color = diffexpressed_plot, label = top_ten_label_mix_plot)) + 
  geom_point(size = 0.8, na.rm = FALSE) + scale_color_manual(values = c(DEG_color[1], DEG_color[2], "grey70"), labels = c("Up-regulated", "Down-regulated", "Not significant")) + 
  geom_vline(xintercept = c(-FC_bound, FC_bound), col="grey30", linetype = "22") + geom_hline(yintercept = -log10(P_bound), col = "grey30", linetype = "22") + 
  geom_text_repel(max.overlaps = 30, na.rm = FALSE, color = "black", size = 3) +
  # geom_label_repel(max.overlaps = 30, na.rm = FALSE, fill = alpha("dodgerblue", 0.2), color = "black", size = 3, box.padding = 0.1, label.size = 0) + 
  geom_segment(aes(x = -3, y = -70, xend = -8.5, yend = -70), arrow = arrow(length = unit(0.3, "cm")), color = DEG_color[2], size = 1) + 
  annotate("text", x = -8, y = -50, label = "Higher in Control", color = DEG_color[2], size = 7, hjust = 0) + 
  geom_segment(aes(x = 2.8, y = -70, xend = 7.3, yend = -70), arrow = arrow(length = unit(0.3, "cm")), color = DEG_color[1], size = 1) + 
  annotate("text", x = 3, y = -50, label = "Higher in IHD", color = DEG_color[1], size = 7, hjust = 0) + 
  labs(color = '', x = expression("log"[2]*"(Fold Change)"), y = expression("-log"[10]*"(adjusted p-value)")) + 
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, max_p * 1.05), clip = "off") + theme_linedraw() + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                           text = element_text(size = 20), legend.position = "none")
# DEG_vocalno
ggsave(paste0(LOO_DEG_Volcano_dir, "/DEG_vocalno.pdf"), plot = DEG_vocalno, width = 4.5, height = 4, dpi = 300)
