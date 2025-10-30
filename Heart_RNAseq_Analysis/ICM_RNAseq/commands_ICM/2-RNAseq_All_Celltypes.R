library(Seurat)
library(patchwork)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(scCustomize)
library(UpSetR)
library(harmony)
library(sctransform)
library(DoubletFinder)
library(glmGamPoi)
library(pheatmap)

Ctrl_IHD_color <- c("#2D6EA8", "#DD555B")
color_list <- c("red", "#2ca02c", "dodgerblue", "#ff7f0e", "#9467bd", "#e377c2", "#8c564b", "#bcbd22", "#7f7f7f", "#17becf",
                "#393b79", "#5254a3", "#6b6ecf", "#9c9ede", "#637939", "#8ca252", "#b5cf6b", "#cedb9c", "#8c6d31", "#bd9e39",
                "#e7ba52", "#e7969c", "#d6616b", "#ad494a", "#843c39", "#f5f5f5", "#fdd0a2", "#fb6a4a", "#cb181d", "#a50f15",
                "#3182bd", "#6baed6", "#9ecae1", "#c6dbef", "#e6550d", "#fd8d3c", "#fdae6b", "#fdd0a2", "#31a354", "#74c476",
                "#a1d99b", "#c7e9c0", "#756bb1", "#9e9ac8", "#bcbddc")
color_list_02 <- c("#f22c4b", "#9c9ede", "#fd8d3c", "#73B8E1", "#31a354", "#bcbd22", "#e377c2", "#8c564b")
##### set working and results directories
wd_dir <- getwd()
Unsupervised_result_dir <- paste0(wd_dir, "/results_ICM/2-All_Celltypes/Unsupervised_Analysis")
Supervised_result_dir <- paste0(wd_dir, "/results_ICM/2-All_Celltypes/Supervised_Analysis")
dir.create(Unsupervised_result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(Supervised_result_dir, recursive = TRUE, showWarnings = FALSE)

######################################################################################################
######################################### continue from O2 ###########################################
control.disease.combined  <- readRDS(paste0(wd_dir, "/results_ICM/0-Seurat_Object/combined.Seurat.obj_ready_for_annotation.RDS"))
sample_names_ordered <- c("1452", "1650", "1690", "1716", "1739", "1763", "1785", "1801", "1364", "1579", "1693", "1703", "1733", "1773", "1800")
donor_condition_ordered <- c("1452_Control", "1650_Control", "1690_Control", "1716_Control", "1739_Control", "1763_Control", "1785_Control", "1801_Control", 
                             "1364_ICM", "1579_ICM", "1693_ICM", "1703_ICM", "1733_ICM", "1773_ICM", "1800_ICM")

##### basic summary: Calculate the average number of cells and genes per sample
average_metrics <- control.disease.combined@meta.data %>% group_by(samples) %>% summarise(avg_cells_per_sample = n(), avg_genes_per_cell = mean(nFeature_RNA))
total_avg_cells <- nrow(control.disease.combined@meta.data) / length(unique(control.disease.combined@meta.data$samples))
total_avg_genes <- mean(control.disease.combined@meta.data$nFeature_RNA)
num_nuclei <- nrow(control.disease.combined@meta.data)

##### Unsupervised clustering
RNA_resolution <- "RNA_snn_res.1"
Idents(object = control.disease.combined) <- RNA_resolution

dir.create(paste0(Unsupervised_result_dir, "/Clustering"), showWarnings = FALSE)
pdf(paste0(Unsupervised_result_dir, "/Clustering/Unsupervised_clustering.pdf"), width = 6.5, height = 5)
  DimPlot(control.disease.combined, reduction = 'umap', group.by = RNA_resolution, label = TRUE, label.size = 3, cols = color_list, order = T) + ggtitle(NULL) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

pdf(paste0(Unsupervised_result_dir, "/Clustering/Unsupervised_clustering_cond.pdf"), width = 15, height = 8)
  DimPlot(control.disease.combined, reduction = 'umap', split.by = "condition", label = TRUE, label.size = 4, pt.size = 0.5, cols = color_list, order = T) + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 16))
dev.off()

pdf(paste0(Unsupervised_result_dir, "/Clustering/Unsupervised_clustering_donor.pdf"), width = 50, height = 4)
  DimPlot(control.disease.combined, reduction = 'umap', split.by = "donor_condition", label = TRUE, label.size = 3, cols = color_list, order = T) + 
    facet_grid(. ~  factor(donor_condition, level = donor_condition_ordered), scale = "free_x") + theme_linedraw() + ggtitle(NULL) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

######################################################################################################
######################################## Cell type annotation ########################################
######################################################################################################
##### all markers
CM_markers <- c("TNNT2", "ACTN2", "MYOZ2", "MLIP", "MYH7") # Cardiomyocyte
Fib_markers <- c("DCN", "BICC1", "ABCA6") # Fibroblast
Endo_markers <- c("VWF", "PECAM1", "ANO2", "DOCK9") # Endothelial
Pericyte_markers <- c("RGS5", "PDGFRB", "EGFLAM", "GUCY1A2") # Pericyte
VSMC_markers <- c("MYH11", "ACTA2", "TAGLN")
Neuronal_markers <- c("NRXN1", "NRXN3", "CADM2", "KIRREL3") # Neuronal
Immune_markers <- c("PTPRC", "DOCK2", "IQGAP2", "RUNX1") # Immune
Adipocyte_markers <- c("ADIPOQ", "PLIN1", "PPARG") # Adipocyte

all_markers <- c(CM_markers, Fib_markers, Endo_markers, Pericyte_markers, VSMC_markers, Neuronal_markers, Immune_markers, Adipocyte_markers)
marker_groups <- c(rep("Cardiomyocyte markers", length(CM_markers)), rep("Fibroblast markers", length(Fib_markers)), rep("Endothelial markers", length(Endo_markers)), rep("Pericyte markers", length(Pericyte_markers)), 
                   rep("VSMC markers", length(VSMC_markers)), rep("Neuronal markers", length(Neuronal_markers)), rep("Immune markers", length(Immune_markers)), rep("Adipocyte markers", length(Adipocyte_markers)))
gene_group_df <- data.frame(gene = all_markers, group = marker_groups) %>% 
  mutate(group = factor(group, level = c("Cardiomyocyte markers", "Fibroblast markers", "Endothelial markers", "Pericyte markers", "VSMC markers", "Neuronal markers", "Immune markers", "Adipocyte markers")))

##### Dotplot for annotation
dir.create(paste0(Unsupervised_result_dir, "/Dotplot"), showWarnings = FALSE)
control.disease.combined_Control <- subset(control.disease.combined, subset = condition == "Control")
control.disease.combined_Disease <- subset(control.disease.combined, subset = condition == "ICM")

# rm(dotplot_all)
dotplot_all <- DotPlot(control.disease.combined, features = all_markers, group.by = RNA_resolution)
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all <- dotplot_all + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + scale_color_gradient(low = "white", high = "black") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "black", linetype = "dashed"), text = element_text(face = "bold", size = 12)) + RotatedAxis()
ggsave(paste0(Unsupervised_result_dir, "/Dotplot/Unsupervised_dotplot_all.pdf"), plot = dotplot_all, width = 16, height = 8, dpi = 300)

# rm(dotplot_Control)
dotplot_Control <- DotPlot(control.disease.combined_Control, features = all_markers, group.by = RNA_resolution)
dotplot_Control$data$gene_group <- factor(gene_group_df$group[match(dotplot_Control$data$features.plot, gene_group_df$gene)])
dotplot_Control$data$condition <- "Control"
dotplot_Control <- dotplot_Control + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + scale_color_gradient(low = "white", high = "dodgerblue4") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "dodgerblue4", linetype = "dashed"), text = element_text(face = "bold", size = 12)) + RotatedAxis()
ggsave(paste0(Unsupervised_result_dir, "/Dotplot/Unsupervised_dotplot_Control.pdf"), plot = dotplot_Control, width = 13, height = 8, dpi = 300)

# rm(dotplot_Disease)
dotplot_Disease <- DotPlot(control.disease.combined_Disease, features = all_markers, group.by = RNA_resolution)
dotplot_Disease$data$gene_group <- factor(gene_group_df$group[match(dotplot_Disease$data$features.plot, gene_group_df$gene)])
dotplot_Disease$data$condition <- "ICM"
dotplot_Disease <- dotplot_Disease + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + scale_color_gradient(low = "white", high = "firebrick4") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "firebrick4", linetype = "dashed"), text = element_text(face = "bold", size = 12)) + RotatedAxis() + guides(colour = guide_colourbar())
ggsave(paste0(Unsupervised_result_dir, "/Dotplot/Unsupervised_dotplot_Disease.pdf"), plot = dotplot_Disease, width = 13, height = 8, dpi = 300)

##### top10 markers based on FC in each cluster and annotate cell types
# cluster.markers <- FindMarkers(control.disease.combined, ident.1 = 11, min.pct = 0.50) %>% arrange(desc(avg_log2FC))
# row.names(cluster.markers)[1:30]
# top10_FC_markers_each_cluster <- FindAllMarkers(control.disease.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.50)
# write.csv(top10_FC_markers_each_cluster, paste0(result_dir, "/Marker_files/all_markers_each_cluster_unsupervised.csv"), row.names=FALSE)

##### manually annotation
annotated_clusters_list <- c(
  "Fibroblast",    #0 
  "Cardiomyocyte", #1 
  "Endothelial",   #2 
  "Pericyte",      #3 
  "Cardiomyocyte", #4 
  "Immune",        #5 
  "Fibroblast",    #6 
  "Cardiomyocyte", #7 
  "Immune",        #8 
  "Fibroblast",    #9 
  "Endothelial",   #10 
  "Endothelial",   #11
  "VSMC-Pericyte", #12 
  "Endothelial",   #13 
  "Immune",        #14 
  "Endothelial",   #15 
  "Neuronal",      #16 
  "Adipocyte",     #17 
  "Endothelial",   #18
  "Fibroblast",    #19 
  "Cardiomyocyte", #20 
  "Endothelial",   #21 
  "Cardiomyocyte", #22
  "Cardiomyocyte", #23 
  "Immune",        #24 
  "Cardiomyocyte"  #25 
)

names(annotated_clusters_list) <- paste0("cluster_", 0:(length(annotated_clusters_list) - 1))
annotated_clusters <- sapply(control.disease.combined$RNA_snn_res.1, function(x) annotated_clusters_list[paste0("cluster_", x)])
control.disease.combined@meta.data <- cbind(control.disease.combined@meta.data, annotated_clusters)

##### Reorder the cell types in the metadata
cell_type_order <- c("Cardiomyocyte", "Fibroblast", "Endothelial", "Pericyte", "VSMC-Pericyte", "Neuronal", "Immune", "Adipocyte")
control.disease.combined$annotated_clusters <- factor(control.disease.combined$annotated_clusters, levels = cell_type_order)

dir.create(paste0(Supervised_result_dir, "/Clustering"), showWarnings = FALSE)
pdf(paste0(Supervised_result_dir, "/Clustering/Supervised_clustering.pdf"), width = 8, height = 6)
  DimPlot(control.disease.combined, reduction = 'umap', group.by = "annotated_clusters", label = T, label.size = 5, repel = T, cols = color_list_02, order = T) + ggtitle(NULL) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 20))
dev.off()

pdf(paste0(Supervised_result_dir, "/Clustering/Supervised_clustering_cond.pdf"), width = 13, height = 5.5)
  DimPlot(control.disease.combined, reduction = 'umap', group.by = "annotated_clusters", split.by = "condition", label = T, label.size = 5, repel = T, cols = color_list_02, order = T) + ggtitle(NULL) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24))
dev.off()

pdf(paste0(Supervised_result_dir, "/Clustering/Supervised_clustering_sample_donor.pdf"), width = 50, height = 4)
  DimPlot(control.disease.combined, reduction = 'umap', group.by = "annotated_clusters", split.by = "donor_condition", label = T, label.size = 3, cols = color_list_02, order = T) + 
    facet_grid(. ~  factor(donor_condition, level = donor_condition_ordered), scale = "free_x") + ggtitle(NULL) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

dir.create(paste0(Supervised_result_dir, "/Dotplot"), showWarnings = FALSE)
control.disease.combined_Control <- subset(control.disease.combined, subset = condition == "Control")
control.disease.combined_Disease <- subset(control.disease.combined, subset = condition == "ICM")

rm(dotplot_all)
dotplot_all <- DotPlot(control.disease.combined, features = all_markers, group.by = "annotated_clusters")
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all <- dotplot_all + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + theme_linedraw() + RotatedAxis() + scale_color_gradient(low = "white", high = "black") + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.1), panel.grid.minor = element_blank(), 
        panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(paste0(Supervised_result_dir, "/Dotplot/Supervised_dotplot_all.pdf"), plot = dotplot_all, width = 17, height = 3.5, dpi = 300)

rm(dotplot_Control)
dotplot_Control <- DotPlot(control.disease.combined_Control, features = all_markers, group.by = "annotated_clusters")
dotplot_Control$data$gene_group <- factor(gene_group_df$group[match(dotplot_Control$data$features.plot, gene_group_df$gene)])
dotplot_Control$data$condition <- "Control"
dotplot_Control <- dotplot_Control + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + theme_linedraw() + RotatedAxis() + scale_color_gradient(low = "white", high = "dodgerblue4") + guides(colour = guide_colourbar()) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
ggsave(paste0(Supervised_result_dir, "/Dotplot/Supervised_dotplot_Control.pdf"), plot = dotplot_Control, width = 14, height = 4, dpi = 300)

rm(dotplot_Disease)
dotplot_Disease <- DotPlot(control.disease.combined_Disease, features = all_markers, group.by = "annotated_clusters")
dotplot_Disease$data$gene_group <- factor(gene_group_df$group[match(dotplot_Disease$data$features.plot, gene_group_df$gene)])
dotplot_Disease$data$condition <- "ICM"
dotplot_Disease <- dotplot_Disease + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + theme_linedraw() + RotatedAxis() + scale_color_gradient(low = "white", high = "firebrick4") + guides(colour = guide_colourbar()) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
ggsave(paste0(Supervised_result_dir, "/Dotplot/Supervised_dotplot_Disease.pdf"), plot = dotplot_Disease, width = 14, height = 4, dpi = 300)

######################################################################################################
##### statistic for unannotated and annotated clustering
num_per_donor_cluster_annotated <- table(control.disease.combined$donor_condition, control.disease.combined$annotated_clusters)
num_per_donor_cluster_annotated <- num_per_donor_cluster_annotated[unique(control.disease.combined$donor_condition), ]
num_per_sample_cluster_annotated <- table(control.disease.combined$samples, control.disease.combined$annotated_clusters)
num_per_sample_cluster_annotated <- num_per_sample_cluster_annotated[unique(control.disease.combined$samples), ]

pct_donor_cluster_annotated <- num_per_sample_cluster_annotated %>% 
  as.data.frame.matrix() %>% mutate(Total = rowSums(.)) %>% 
  mutate(across(-Total, ~ round(.x / Total *100 , 5))) %>% select(-Total) %>% 
  slice(match(sample_names_ordered, rownames(.))) %>% 
  mutate(condition = c(rep("Control", 8), rep("ICM", 7)))

pct_long <- pct_donor_cluster_annotated %>% melt() %>% 
  setNames(c("condition", "cluster", "ratio")) %>% mutate(condition = factor(condition, level = c("Control", "ICM"))) %>% mutate(ratio = as.numeric(ratio))

pdf(paste0(Supervised_result_dir, "/Dotplot/Supervised_cluster_prop_comparison.pdf"), width = 20, height = 6)
  barplot_annotated_cluster_ratio <- ggbarplot(pct_long, x = "cluster", y = "ratio", color = "condition", 
                                              label = TRUE, lab.nb.digits = 2, add = c("mean_se", "jitter"), position = position_dodge(0.9), palette = c("dodgerblue4", "firebrick4")) +
    stat_compare_means(aes(group = condition), method = "wilcox", label = "p.format", label.y = 1.02 * max(pct_long$ratio, na.rm = TRUE)) + 
    stat_compare_means(aes(group = condition), method = "wilcox", label = "p.signif", label.y = 1.06 * max(pct_long$ratio, na.rm = TRUE)) + 
    labs(x = "Unannotated Cluster", y = "proportion of cell") + theme_bw(base_size = 24) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 20)) + 
    scale_fill_manual(values = c("dodgerblue4", "firebrick4"))
  print(barplot_annotated_cluster_ratio)
dev.off()

##### Cell Type Composition
melt_num_per_donor_cluster_annotated <- melt(num_per_sample_cluster_annotated) %>% mutate(value = value + 1) %>% 
  mutate(condition = ifelse(Var1 %in% sample_names_ordered[1:8], "Control", ifelse(Var1 %in% sample_names_ordered[9:17], "ICM", value)))
melt_pct_donor_cluster_annotated <- melt(as.matrix(pct_donor_cluster_annotated[, -9])) %>% 
  mutate(condition = ifelse(Var1 %in% sample_names_ordered[1:8], "Control", ifelse(Var1 %in% sample_names_ordered[9:17], "ICM", value)))

red_palette <- c("grey", "#FFCCCC", "#FFB2B2", "#FF9999", "#FF6666", "#FF4D4D", "#FF3333","#FF1A1A", "#FF0000", "#CC0000",
                 "#B30000", "#990000", "#800000", "#660000", "#4D0000", "#330000")

pdf(paste0(Supervised_result_dir, "/Dotplot/Cell_count_per_celltype.pdf"), width = 10, height = 8)
  ggplot(melt_num_per_donor_cluster_annotated, aes(x = Var2, y = as.character(Var1), fill = log10(value))) + 
    geom_tile() + geom_text(aes(label = round(value - 1, 1)), color = "white", size = 6) + 
    scale_fill_gradientn(colors = red_palette, name = "# of cells", breaks = log10(c(1, 10, 100, 1000)), labels = c(1, 10, 100, 1000)) + 
    facet_grid(condition ~ ., scales = "free", space = "free") + theme_linedraw() + labs(title = "", x = "Cell type", y = "Donor") + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          text = element_text(face = "bold", size = 16), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()

pdf(paste0(Supervised_result_dir, "/Dotplot/Cell_pct_per_celltype.pdf"), width = 10, height = 9)
  ggplot(melt_pct_donor_cluster_annotated, aes(x = Var2, y = as.character(Var1), fill = value)) +
    geom_tile() + geom_text(aes(label = round(value, 1)), color = "white", size = 6) +
    # scale_fill_gradientn(colors = red_palette, name = "# of cells", breaks = log10(c(1, 10, 100, 1000)), labels = c(1, 10, 100, 1000)) +
    scale_fill_gradientn(colors = red_palette, name = "% of cells") +
    facet_grid(condition ~ ., scales = "free", space = "free") + theme_linedraw() + labs(title = "", x = "Cell type", y = "Donor") +
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          text = element_text(size = 16), panel.border = element_rect(linewidth = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

##### save object
# saveRDS(control.disease.combined, paste0(wd_dir, "/results_ICM/0-Seurat_Object/Seurat.obj_with_annotation_all_celltypes.RDS"))
# control.disease.combined  <- readRDS(paste0(wd_dir, "/results_ICM/0-Seurat_Object/Seurat.obj_with_annotation_all_celltypes.RDS"))
Idents(object = control.disease.combined) <- "annotated_clusters"
control.disease.combined <- subset(control.disease.combined, idents = c("Cardiomyocyte", "Fibroblast"))

##### check the marker expression in cluster
# gene_list <- rownames(control.disease.combined)
# MMR_pathway_marker_list %in% gene_list
# list_to_check <- Collagen_marker_list
# Idents(object = control.disease.combined) <- ""
# VlnPlot(control.disease.combined, group.by = "annotated_clusters", features = list_to_check, split.by = "condition", ncol = 5)
# VlnPlot(control.disease.combined, features = list_to_check, split.by = "condition", ncol = 5)

gene_sets <- list(
  metagene_cardiomyocyte = c("TNNT2", "ACTN2", "MYOZ2", "MYPN", "MLIP"), 
  metagene_fibroblast = c("DCN", "MMP2", "VIM", "ACTA2"), 
  metagene_hypoxia = c("JUN", "NR4A1", "ZEB1", "ZEB2"), 
  metagene_collagen = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2"), 
  metagene_fibrotic = c("PDGFRA", "PDGFRB", "FAP", "TGFBI", "LOX")
  # metagene_MMR = c("MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2"), 
  # metagene_BER = c("APEX1", "FEN1", "LIG1", "LIG3", "MBD4", "NEIL1", "NEIL12", "OGG1", "POLB", "PNKP"), 
  # metagene_NHEJ = c("XRCC5", "XRCC6", "PRKDC", "LIG4", "DCLRE1C", "NHEJ1", "XRCC4", "PAXX", "MRE11", "TP53BP1", "SHLD1", "SHLD2", "SHLD3"), 
  # metagene_HR = c("RAD51", "BLM", "FAN1", "SLX4", "PALB2", "FANCA", "FANCG", "CREBBP", "UBE2I", "MNT")
)

control.disease.combined <- AddModuleScore(control.disease.combined, features = gene_sets, name = names(gene_sets))

metadata_CM_Fib_02 <- control.disease.combined@meta.data
all_metagene_grouped_avg_02 <- metadata_CM_Fib_02 %>% group_by(condition, samples) %>%
  summarize(cardiomyocyte_metagene = mean(metagene_cardiomyocyte1, na.rm = TRUE), fibroblast_metagene = mean(metagene_fibroblast2, na.rm = TRUE), 
            hypoxia_metagene = mean(metagene_hypoxia3, na.rm = TRUE), collagen_metagene = mean(metagene_collagen4, na.rm = TRUE), 
            fibrotic_metagene = mean(metagene_fibrotic5, na.rm = TRUE))

all_metagene_grouped_avg_melt_02 <- all_metagene_grouped_avg_02 %>% 
  melt(id.vars = c("condition", "samples"),
       measure.vars = c("cardiomyocyte_metagene", "fibroblast_metagene", "hypoxia_metagene", "collagen_metagene", "fibrotic_metagene"),
       variable.name = "metagenes", value.name = "avg_expr")

pdf(paste0(Supervised_result_dir, "/metagene_expr_boxplot.pdf"), width = 18, height = 8)
  ggplot(all_metagene_grouped_avg_melt_02, aes(x = metagenes, y = avg_expr, color = condition)) + 
    geom_boxplot(aes(fill = condition), position = position_dodge(), alpha = 0.5, outlier.shape = NA) + 
    geom_jitter(aes(color = condition), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 2) + 
    stat_compare_means(aes(group = condition), label = "p.format", label.y = 1.02 * max(all_metagene_grouped_avg_melt_02$avg_expr, na.rm = TRUE)) + 
    stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.06 * max(all_metagene_grouped_avg_melt_02$avg_expr, na.rm = TRUE)) + 
    scale_fill_manual(values = Ctrl_IHD_color) + scale_color_manual(values = Ctrl_IHD_color) + # facet_wrap(. ~  annotated_clusters) +
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             text = element_text(face = "bold", size = 18), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + 
    labs(x = "", y = "Average metagene expression")
dev.off()
