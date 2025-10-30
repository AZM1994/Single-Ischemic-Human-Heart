library(Seurat)
library(patchwork)
library(dplyr)
library(tidyr)
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
color_list <- c("#e6550d", "#9ecae1", "#a1d99b", "#d6616b", "#17becf", "#fd8d3c", "#9c9ede", "#1f77b4", "#8c6d31", "#d62728",
                "#bcbd22", "#8ca252", "#e7ba52", "#e377c2", "#9467bd", "#2ca02c", "#393b79", "#bd9e39", "#a50f15")
color_list_grey <- c("#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#d62728",
                   "#7f7f7f", "#8ca252", "#7f7f7f", "#7f7f7f", "#9467bd", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f")
color_list_regroup <- c("#E97679", "#EE5A9B", "yellow", "#fd8d3c", "#FFB6C1", "#f22c4b", "#FFC34D", "#FF9999",
                        "#a1d99b", "limegreen", "#66B032", "olivedrab", 
                        "#9ecae1", "dodgerblue", "#9c9ede", "#800080", "#8A2BE2", "#e377c2", "#4169E1")
CM_Fibro_color <- c("#f22c4b", "#66B032", "#9c9ede")

##### set working and results directories
wd_dir <- getwd()
Unsupervised_result_dir <- paste0(wd_dir, "/results_IHD/3-CMFib/Unsupervised_Analysis")
Supervised_result_dir <- paste0(wd_dir, "/results_IHD/3-CMFib/Supervised_Analysis")
dir.create(Unsupervised_result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(Supervised_result_dir, recursive = TRUE, showWarnings = FALSE)

##### read sample metadata and gene lists
sample_metadata <- read.csv(paste0(wd_dir, "/data/sample_meta.csv"), header = TRUE)
sample_names <- sample_metadata$sample_name

################################################################################
############################## extract CM and Fib ##############################
################################################################################
# control.disease.combined <- readRDS(paste0(wd_dir, "/results_IHD/0-Seurat_Object/Seurat.obj_with_annotation_all_celltypes.RDS"))
# all_cell_type <- unique(control.disease.combined$annotated_clusters)
# control.disease.combined_CM_Fib <- subset(control.disease.combined, subset = annotated_clusters %in% c("Cardiomyocyte", "Fibroblast"))
# unique(control.disease.combined_CM_Fib$annotated_clusters)
# control.disease.combined_CM_Fib_test <- RunPCA(control.disease.combined_CM_Fib)
# ElbowPlot(control.disease.combined_CM_Fib_test, ndims = 50, reduction = "pca")

# DefaultAssay(control.disease.combined_CM_Fib) <- 'RNA'
# control.disease.combined_CM_Fib <- RunPCA(control.disease.combined_CM_Fib, npcs = 15, verbose = TRUE) %>% 
#   # RunHarmony(group.by.vars = c("samples", "technique"), lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
#   # RunTSNE(reduction = "harmony", dims = 1:15, verbose = TRUE) %>%
#   RunUMAP(reduction = "harmony", dims = 1:15, verbose = TRUE) %>% 
#   FindNeighbors(reduction = "harmony", dims = 1:15, verbose = TRUE) %>% 
#   FindClusters(resolution = c(0.8, 1.0, 1.2, 1.4), verbose = TRUE)

# saveRDS(control.disease.combined_CM_Fib, paste0(wd_dir, "/results_IHD/0-Seurat_Object/Seurat.obj_extractCM_Fibro_ready_for_subclustering.RDS"))
control.disease.combined_CM_Fib  <- readRDS(paste0(wd_dir, "/results_IHD/0-Seurat_Object/Seurat.obj_extractCM_Fibro_ready_for_subclustering.RDS"))
selected_sample_names <- unique(control.disease.combined_CM_Fib$samples)
selected_donor_names <- unique(control.disease.combined_CM_Fib$donor)

RNA_resolution <- "RNA_snn_res.1"
Idents(object = control.disease.combined_CM_Fib) <- RNA_resolution
# table_stat <- t(table(control.disease.combined_CM_Fib$RNA_snn_res.1, control.disease.combined_CM_Fib$samples))
# table_stat <- table_stat[unique(control.disease.combined_CM_Fib$samples), ]
level_list <- sort(unique(control.disease.combined_CM_Fib$RNA_snn_res.1))
control.disease.combined_CM_Fib$RNA_snn_res.1 <- factor(control.disease.combined_CM_Fib$RNA_snn_res.1, levels = level_list)

control.disease.combined_CM_Fib@meta.data <- control.disease.combined_CM_Fib@meta.data %>%
  mutate(regroup = recode(RNA_snn_res.1, 
                          "0" = "17", "1" = "0", "2" = "13", "3" = "4", "4" = "14", 
                          "5" = "1", "6" = "2", "7" = "5", "8" = "16", "9" = "11", 
                          "10" = "12", "11" = "8", "12" = "3", "13" = "15", "14" = "10", 
                          "15" = "7", "16" = "18", "17" = "6", "18" = "9"), 
         regroup = factor(regroup, levels = level_list))

dir.create(paste0(Unsupervised_result_dir, "/Clustering"), showWarnings = FALSE)
pdf(paste0(Unsupervised_result_dir, "/Clustering/Unsupervised_clustering_sample.pdf"), width = 50, height = 4)
  DimPlot(control.disease.combined_CM_Fib, reduction = "umap", split.by = "samples", label = TRUE, label.size = 3, cols = color_list, order = T) + 
    facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x") + theme_linedraw() + ggtitle(NULL) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

pdf(paste0(Unsupervised_result_dir, "/Clustering/Unsupervised_clustering_donor.pdf"), width = 28, height = 4)
  DimPlot(control.disease.combined_CM_Fib, reduction = "umap", split.by = "donor", label = TRUE, label.size = 3, cols = color_list, order = T) + 
    facet_grid(. ~  factor(donor, level = selected_donor_names), scale = "free_x") + theme_linedraw() + ggtitle(NULL) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

pdf(paste0(Unsupervised_result_dir, "/Clustering/Unsupervised_clustering_cond.pdf"), width = 10, height = 5)
  DimPlot(control.disease.combined_CM_Fib, reduction = "umap", group.by = RNA_resolution, split.by = "condition", label = TRUE, label.size = 4, cols = color_list, order = T) + 
    scale_color_manual(values = color_list, breaks = as.character(level_list)) + theme_linedraw() + ggtitle(NULL) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 18))
dev.off()

pdf(paste0(Unsupervised_result_dir, "/Clustering/Unsupervised_clustering_cond_8_9_10_11.pdf"), width = 10, height = 5)
  DimPlot(control.disease.combined_CM_Fib, reduction = "umap", group.by = "regroup", split.by = "condition", label = TRUE, label.size = 4, cols = color_list_regroup, order = T) + theme_linedraw() + ggtitle(NULL) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24))
dev.off()

################################################################################
##### dotplot for marker genes
CM_markers <- c("TNNT2", "ACTN2", "MYOZ2", "MYPN", "MLIP") # Cardiomyocyte
Fib_markers <- c("DCN", "APOD", "BICC1", "ABCA6") # Fibroblast
all_markers <- c(CM_markers, Fib_markers)
marker_groups <- c(rep("Cardiomyocyte markers", length(CM_markers)), rep("Fibroblast markers", length(Fib_markers)))
gene_group_df <- data.frame(gene = all_markers, group = marker_groups) %>% mutate(group = factor(group, level = c("Cardiomyocyte markers", "Fibroblast markers")))

dir.create(paste0(Unsupervised_result_dir, "/Dotplot"), showWarnings = FALSE)
control.disease.combined_CM_Fib_Control <- subset(control.disease.combined_CM_Fib, subset = condition == "Control")
control.disease.combined_CM_Fib_Disease <- subset(control.disease.combined_CM_Fib, subset = condition == "IHD")

dotplot_all <- DotPlot(control.disease.combined_CM_Fib, features = all_markers, group.by = "regroup")
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all <- dotplot_all + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + scale_color_gradient(low = "white", high = "black") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "black", linetype = "dashed"), text = element_text(face = "bold", size = 12)) + RotatedAxis()
ggsave(paste0(Unsupervised_result_dir, "/Dotplot/Unsupervised_dotplot_all.pdf"), plot = dotplot_all, width = 10, height = 6, dpi = 600)

rm(dotplot_Control)
dotplot_Control <- DotPlot(control.disease.combined_CM_Fib_Control, features = all_markers, group.by = "regroup")
dotplot_Control$data$gene_group <- factor(gene_group_df$group[match(dotplot_Control$data$features.plot, gene_group_df$gene)])
dotplot_Control$data$condition <- "Control"
dotplot_Control <- dotplot_Control + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + scale_color_gradient(low = "black", high = "dodgerblue1") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), text = element_text(size = 12), 
                           legend.background = element_rect(fill = "white", color = NA), panel.border = element_rect(size = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  geom_tile(data = subset(dotplot_Control$data, id %in% c("8", "9", "10", "11")), fill = "yellow", color = "#FFFF004D", alpha = 0.65) + 
  geom_point(aes(color = avg.exp.scaled, size = pct.exp), show.legend = F) + guides(colour = guide_colourbar())
ggsave(paste0(Unsupervised_result_dir, "/Dotplot/Unsupervised_dotplot_Control.pdf"), plot = dotplot_Control, width = 10, height = 6, dpi = 600)

rm(dotplot_Disease)
dotplot_Disease <- DotPlot(control.disease.combined_CM_Fib_Disease, features = all_markers, group.by = "regroup")
dotplot_Disease$data$gene_group <- factor(gene_group_df$group[match(dotplot_Disease$data$features.plot, gene_group_df$gene)])
dotplot_Disease$data$condition <- "IHD"
dotplot_Disease <- dotplot_Disease + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + scale_color_gradient(low = "black", high = "firebrick1") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), text = element_text(size = 12), 
                           panel.border = element_rect(size = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  geom_tile(data = subset(dotplot_Control$data, id %in% c("8", "9", "10", "11")), fill = "yellow", color = "#FFFF004D", alpha = 0.65) + 
  geom_point(aes(color = avg.exp.scaled, size = pct.exp)) + guides(colour = guide_colourbar())
ggsave(paste0(Unsupervised_result_dir, "/Dotplot/Unsupervised_dotplot_Disease.pdf"), plot = dotplot_Disease, width = 10, height = 6, dpi = 600)

# top10_FC_markers_each_cluster_CM <- FindAllMarkers(control.disease.combined_CM_Fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.50)
# top10_FC_markers_each_cluster_CM <- top10_FC_markers_each_cluster_CM %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)
##### further filtering
min_genes <- 200  # Minimum number of genes expressed
max_genes <- 3000 # Maximum number of genes expressed
min_counts <- 500  # Minimum total counts (optional)
control.disease.combined_CM_Fib[["nGene"]] <- Matrix::colSums(GetAssayData(control.disease.combined_CM_Fib, slot = "counts") > 0)
control.disease.combined_CM_Fib[["nCount"]] <- Matrix::colSums(GetAssayData(control.disease.combined_CM_Fib, slot = "counts"))
control.disease.combined_CM_Fib <- subset(control.disease.combined_CM_Fib, subset = nGene > min_genes & nGene < max_genes & nCount > min_counts)

##### simple statistics
num_per_sample_cluster_unannotated <- table(control.disease.combined_CM_Fib$donor, control.disease.combined_CM_Fib$regroup)
num_per_sample_cluster_unannotated <- num_per_sample_cluster_unannotated[as.character(unique(control.disease.combined_CM_Fib$donor)), ]
num_per_sample_cluster_annotated <- table(control.disease.combined_CM_Fib$donor, control.disease.combined_CM_Fib$annotated_clusters)
num_per_sample_cluster_annotated <- num_per_sample_cluster_annotated[as.character(unique(control.disease.combined_CM_Fib$donor)), ]

####################################################################
##### plot metagenes
##### check the marker expression in cluster
# Fib_marker_list <- c("DCN", "LUM", "COL1A1", "COL3A1", "PDGFRA", "MMP2", "VIM")
# gene_list <- rownames(control.disease.combined_CM_Fib)
# Fib_marker_list %in% gene_list
# list_to_check <- Fib_marker_list
# VlnPlot(control.disease.combined_CM_Fib, features = list_to_check, idents = c(8, 11), split.by = "condition", ncol = 5)

gene_sets <- list(
  Cardiomyocyte_metagene = c("TNNT2", "ACTN2", "MYOZ2", "MYPN", "MLIP"), 
  Fibroblast_metagene = c("DCN", "MMP2", "VIM", "ACTA2"), 
  Hypoxia_metagene = c("JUN", "NR4A1", "ZEB1", "ZEB2"), 
  Collagen_metagene = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2"), 
  Fibrotic_metagene = c("PDGFRA", "PDGFRB", "FAP", "TGFBI", "LOX")
  # MMR_metagene = c("MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2"), 
  # BER_metagene = c("APEX1", "FEN1", "LIG1", "LIG3", "MBD4", "NEIL1", "NEIL12", "OGG1", "POLB", "PNKP"), 
  # NHEJ_metagene = c("XRCC5", "XRCC6", "PRKDC", "LIG4", "DCLRE1C", "NHEJ1", "XRCC4", "PAXX", "MRE11", "TP53BP1", "SHLD1", "SHLD2", "SHLD3"), 
  # HR_metagene = c("RAD51", "BLM", "FAN1", "SLX4", "PALB2", "FANCA", "FANCG", "CREBBP", "UBE2I", "MNT")
)

control.disease.combined_CM_Fib <- AddModuleScore(control.disease.combined_CM_Fib, features = gene_sets, name = names(gene_sets))

added_cols <- tail(colnames(control.disease.combined_CM_Fib@meta.data), length(gene_sets))
names(added_cols) <- names(gene_sets)
colnames(control.disease.combined_CM_Fib@meta.data)[match(added_cols, colnames(control.disease.combined_CM_Fib@meta.data))] <- names(added_cols)
metadata_CM_Fib <- control.disease.combined_CM_Fib@meta.data

all_metagene_grouped_avg <- metadata_CM_Fib %>% group_by(regroup, condition, samples) %>% 
  summarize(Cardiomyocyte = mean(Cardiomyocyte_metagene, na.rm = TRUE), Fibroblast = mean(Fibroblast_metagene, na.rm = TRUE), 
            Hypoxia = mean(Hypoxia_metagene, na.rm = TRUE), Collagen = mean(Collagen_metagene, na.rm = TRUE), 
            Fibrotic = mean(Fibrotic_metagene, na.rm = TRUE))

all_metagene_grouped_avg_melt <- all_metagene_grouped_avg %>% 
  filter(regroup %in% c(8, 11)) %>%
  melt(id.vars = c("regroup", "condition", "samples"), 
       measure.vars = c("Cardiomyocyte", "Fibroblast", "Hypoxia", "Collagen", "Fibrotic"), 
       variable.name = "metagenes", value.name = "avg_expr") %>% 
  mutate(regroup.modified = paste0("Cardiomyocyte subcluster ", regroup)) %>% 
  mutate(regroup.modified = factor(regroup.modified, levels = c("Cardiomyocyte subcluster 8", "Cardiomyocyte subcluster 11")))

pdf(paste0(Unsupervised_result_dir, "/metagene_expr_barplot_CMFib.pdf"), width = 21.5, height = 10)
  ggplot(all_metagene_grouped_avg_melt, aes(x = metagenes, y = avg_expr, color = condition)) + 
    geom_boxplot(size = 0.5, position = position_dodge(), outlier.shape = NA) + 
    geom_jitter(aes(color = condition), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 2) + 
    stat_compare_means(aes(group = condition), label = "p.format", size = 8, label.y = 1.1 * max(all_metagene_grouped_avg_melt$avg_expr, na.rm = TRUE)) +
    stat_compare_means(aes(group = condition), label = "p.signif", size = 10, label.y = 0.95 * max(all_metagene_grouped_avg_melt$avg_expr, na.rm = TRUE), show.legend = FALSE) + 
    scale_fill_manual(values = Ctrl_IHD_color) + scale_color_manual(values = Ctrl_IHD_color) + 
    scale_y_continuous(limits = c(-1, 2), breaks = c(-1, 0, 1, 2)) + facet_wrap(. ~  regroup.modified) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + 
    labs(x = "", y = "Normalized Metagene Expression Score")
dev.off()

pdf(paste0(Unsupervised_result_dir, "/metagene_CMFib_feature_plot.pdf"), width = 13, height = 20)
  feature_plot_color <- c("lightgrey", "dodgerblue3")
  p1 <- FeaturePlot(control.disease.combined_CM_Fib, features = "Cardiomyocyte_metagene", cols = feature_plot_color) + ggtitle("Cardiomyocyte_metagene")
  p2 <- FeaturePlot(control.disease.combined_CM_Fib, features = "Fibroblast_metagene", cols = feature_plot_color) + ggtitle("Fibroblast_metagene")
  p3 <- FeaturePlot(control.disease.combined_CM_Fib, features = "Hypoxia_metagene", cols = feature_plot_color) + ggtitle("Hypoxia_metagene")
  p4 <- FeaturePlot(control.disease.combined_CM_Fib, features = "Collagen_metagene", cols = feature_plot_color) + ggtitle("Collagen_metagene")
  p5 <- FeaturePlot(control.disease.combined_CM_Fib, features = "Fibrotic_metagene", cols = feature_plot_color) + ggtitle("Fibrotic_metagene")
  umap_col <- p1 / p2 / p3 / p4 / p5
  v1 <- VlnPlot(control.disease.combined_CM_Fib, features = "Cardiomyocyte_metagene", group.by = "regroup", cols = color_list_regroup, , pt.size = 0) + ggtitle("Cardiomyocyte_metagene") + theme(legend.position = "none")
  v2 <- VlnPlot(control.disease.combined_CM_Fib, features = "Fibroblast_metagene", group.by = "regroup", cols = color_list_regroup, pt.size = 0) + ggtitle("Fibroblast_metagene") + theme(legend.position = "none")
  v3 <- VlnPlot(control.disease.combined_CM_Fib, features = "Hypoxia_metagene", group.by = "regroup", cols = color_list_regroup, pt.size = 0) + ggtitle("Hypoxia_metagene")
  v4 <- VlnPlot(control.disease.combined_CM_Fib, features = "Collagen_metagene", group.by = "regroup", cols = color_list_regroup, pt.size = 0) + ggtitle("Collagen_metagene") + theme(legend.position = "none")
  v5 <- VlnPlot(control.disease.combined_CM_Fib, features = "Fibrotic_metagene", group.by = "regroup", cols = color_list_regroup, pt.size = 0) + ggtitle("Fibrotic_metagene") + theme(legend.position = "none")
  vln_col  <- v1  / v2  / v3  / v4  / v5
  final_plot <- (umap_col | vln_col) + plot_layout(widths = c(1, 2))
  final_plot
dev.off()

##### find marker genes in 9, 11, 14 for conditions, celltypes
# Idents(object = control.disease.combined_CM_Fib) <- RNA_resolution
# table_condition_RNA_snn_res.1 <- table(control.disease.combined_CM_Fib$condition, control.disease.combined_CM_Fib$RNA_snn_res.1)
# 
# subset_Control <- subset(control.disease.combined_CM_Fib, subset = condition == "Control")
# marker_genes_Control <- FindAllMarkers(subset_Control, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
# cluster9_markers_Control <- marker_genes_Control_CM %>% filter(cluster == 9, p_val_adj < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) %>% write.csv(paste0(Unsupervised_result_dir, "/Markers_Cluster9_Control.csv"), row.names = TRUE)
# cluster11_markers_Control <- marker_genes_Control_CM %>% filter(cluster == 11, p_val_adj < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) %>% write.csv(paste0(Unsupervised_result_dir, "/Markers_Cluster11_Control.csv"), row.names = TRUE)
# subset_IHD <- subset(control.disease.combined_CM_Fib, subset = condition == "IHD")
# marker_genes_IHD <- FindAllMarkers(subset_IHD, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
# cluster9_markers_IHD <- marker_genes_IHD %>% filter(cluster == 9, p_val_adj < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) %>% write.csv(paste0(Unsupervised_result_dir, "/Markers_Cluster9_IHD.csv"), row.names = TRUE)
# cluster11_markers_IHD <- marker_genes_IHD %>% filter(cluster == 11, p_val_adj < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) %>% write.csv(paste0(Unsupervised_result_dir, "/Markers_Cluster11_IHD.csv"), row.names = TRUE)

######################################################################################################
######################################## Cell type annotation ########################################
######################################################################################################
##### manually annotation
annotated_clusters_list_CM <- c("Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", 
                                "Fibrotic-cardiomyocyte", "Fibrotic-cardiomyocyte", "Fibrotic-cardiomyocyte", "Fibrotic-cardiomyocyte", 
                                "Fibroblast", "Fibroblast", "Fibroblast", "Fibroblast", "Fibroblast", "Fibroblast", "Fibroblast")

names(annotated_clusters_list_CM) <- paste0("cluster_", 0:(length(annotated_clusters_list_CM) - 1))
annotated_clusters_CM <- sapply(control.disease.combined_CM_Fib$regroup, function(x) annotated_clusters_list_CM[paste0("cluster_", x)])
control.disease.combined_CM_Fib@meta.data <- cbind(control.disease.combined_CM_Fib@meta.data, annotated_clusters_CM)
cell_type_order <- c("Cardiomyocyte", "Fibrotic-cardiomyocyte", "Fibroblast")

##### Reorder the cell types in the metadata
control.disease.combined_CM_Fib$annotated_clusters_CM <- factor(control.disease.combined_CM_Fib$annotated_clusters_CM, levels = cell_type_order)

dir.create(paste0(Supervised_result_dir, "/Clustering"), showWarnings = FALSE)
pdf(paste0(Supervised_result_dir, "/Clustering/Supervised_clustering.pdf"), width = 8.5, height = 5.5)
  DimPlot(control.disease.combined_CM_Fib, reduction = 'umap', group.by = "annotated_clusters_CM", label = T, label.size = 5, cols = CM_Fibro_color, order = T) + ggtitle(NULL) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 20))
dev.off()

pdf(paste0(Supervised_result_dir, "/Clustering/Supervised_clustering_cond.pdf"), width = 12, height = 5)
  DimPlot(control.disease.combined_CM_Fib, reduction = 'umap', group.by = "annotated_clusters_CM", split.by = "condition", label = T, label.size = 5, cols = CM_Fibro_color, order = T) + ggtitle(NULL) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(size = 0.5), text = element_text(size = 24))
dev.off()

pdf(paste0(Supervised_result_dir, "/Clustering/Supervised_clustering_sample.pdf"), width = 50, height = 4)
  DimPlot(control.disease.combined_CM_Fib, reduction = 'umap', group.by = "annotated_clusters_CM", split.by = "samples", label = T, label.size = 3, cols = CM_Fibro_color, order = T) + facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x") + 
    ggtitle(NULL) + theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

pdf(paste0(Supervised_result_dir, "/Clustering/Supervised_clustering_sample_donor.pdf"), width = 28, height = 4)
  DimPlot(control.disease.combined_CM_Fib, reduction = 'umap', group.by = "annotated_clusters_CM", split.by = "donor", label = T, label.size = 3, cols = CM_Fibro_color, order = T) + facet_grid(. ~  factor(donor, level = selected_donor_names), scale = "free_x") + 
    ggtitle(NULL) + theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

dir.create(paste0(Supervised_result_dir, "/Dotplot"), showWarnings = FALSE)
control.disease.combined_CM_Fib_Control <- subset(control.disease.combined_CM_Fib, subset = condition == "Control")
control.disease.combined_CM_Fib_Disease <- subset(control.disease.combined_CM_Fib, subset = condition == "IHD")

rm(dotplot_all)
dotplot_all <- DotPlot(control.disease.combined_CM_Fib, features = all_markers, group.by = "annotated_clusters_CM")
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all <- dotplot_all + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + theme_linedraw() + RotatedAxis() + scale_color_gradient(low = "white", high = "black") + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
ggsave(paste0(Supervised_result_dir, "/Dotplot/Supervised_dotplot_all.pdf"), plot = dotplot_all, width = 14, height = 4, dpi = 300)

rm(dotplot_Control)
dotplot_Control <- DotPlot(control.disease.combined_CM_Fib_Control, features = all_markers, group.by = "annotated_clusters_CM")
dotplot_Control$data$gene_group <- factor(gene_group_df$group[match(dotplot_Control$data$features.plot, gene_group_df$gene)])
dotplot_Control$data$condition <- "Control"
dotplot_Control <- dotplot_Control + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + theme_linedraw() + RotatedAxis() + scale_color_gradient(low = "white", high = "dodgerblue4") + guides(colour = guide_colourbar()) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
ggsave(paste0(Supervised_result_dir, "/Dotplot/Supervised_dotplot_Control.pdf"), plot = dotplot_Control, width = 14, height = 4, dpi = 300)

rm(dotplot_Disease)
dotplot_Disease <- DotPlot(control.disease.combined_CM_Fib_Disease, features = all_markers, group.by = "annotated_clusters_CM")
dotplot_Disease$data$gene_group <- factor(gene_group_df$group[match(dotplot_Disease$data$features.plot, gene_group_df$gene)])
dotplot_Disease$data$condition <- "IHD"
dotplot_Disease <- dotplot_Disease + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + theme_linedraw() + RotatedAxis() + scale_color_gradient(low = "white", high = "firebrick4") + guides(colour = guide_colourbar()) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
ggsave(paste0(Supervised_result_dir, "/Dotplot/Supervised_dotplot_Disease.pdf"), plot = dotplot_Disease, width = 14, height = 4, dpi = 300)

######################################################################################################
##### simple statistic for unannotated and annotated clustering
num_per_donor_cluster_annotated <- table(control.disease.combined_CM_Fib$donor, control.disease.combined_CM_Fib$annotated_clusters_CM)
num_per_donor_cluster_annotated <- num_per_donor_cluster_annotated[as.character(unique(control.disease.combined_CM_Fib$donor)), ]
num_per_sample_cluster_annotated <- table(control.disease.combined_CM_Fib$samples, control.disease.combined_CM_Fib$annotated_clusters_CM)
num_per_sample_cluster_annotated <- num_per_sample_cluster_annotated[unique(control.disease.combined_CM_Fib$samples), ]

pct_sample_cluster_annotated <- num_per_sample_cluster_annotated %>% 
  as.data.frame.matrix() %>% mutate(Total = rowSums(.)) %>% 
  mutate(across(-Total, ~ round(.x / Total *100 , 5))) %>% select(-Total) %>% 
  slice(match(selected_sample_names, rownames(.))) %>% 
  mutate(condition = c(rep("Control", 7), rep("IHD", 7)), ratio  = Cardiomyocyte / Fibroblast)

pct_donor_cluster_annotated <- num_per_donor_cluster_annotated %>% 
  as.data.frame.matrix() %>% mutate(Total = rowSums(.)) %>% 
  mutate(across(-Total, ~ round(.x / Total *100 , 5))) %>% select(-Total) %>% 
  slice(match(selected_donor_names, rownames(.))) %>% 
  mutate(condition = c(rep("Control", 3), rep("IHD", 5)), ratio  = Cardiomyocyte / Fibroblast)

pct_long <- pct_sample_cluster_annotated %>% tibble::rownames_to_column("donor") %>% 
  melt(id.vars = c("donor", "condition"), measure.vars = c("Cardiomyocyte", "Fibrotic-cardiomyocyte", "Fibroblast"), variable.name = "Cell_Type", value.name = "Percentage")

wilcox_test <- wilcox.test(ratio ~ condition, data = pct_sample_cluster_annotated)
p_value <- wilcox_test$p.value
significance <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns")))
max_y <- max(pct_long$Percentage) + 5

wilcox_test_transition <- wilcox.test(`Fibrotic-cardiomyocyte` ~ condition, data = pct_sample_cluster_annotated)
p_value_transition <- wilcox_test_transition$p.value

pdf(paste0(Supervised_result_dir, "/Dotplot/Supervised_cluster_prop_comparison.pdf"), width = 8, height = 5)
  ggplot(pct_long, aes(x = donor, y = Percentage, fill = Cell_Type)) + geom_bar(stat = "identity") + facet_wrap(~ condition, scales = "free_x") + labs(y = "Percentage", x = "Donor", fill = "Cell Type") + 
    theme_linedraw() + labs(title = "", x = "", y = "Percentage of celltype") + scale_fill_manual(values = CM_Fibro_color) + scale_y_continuous(limits = c(0, 120), breaks = c(0, 25, 50, 75, 100)) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 16), 
          panel.border = element_rect(size = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    ggtitle(paste0(round(p_value, 3), ", ", significance))
dev.off()

##### Cell Type Composition
melt_num_per_donor_cluster_annotated <- melt(num_per_donor_cluster_annotated) %>% mutate(value = value + 1) %>% 
  mutate(condition = ifelse(Var1 %in% selected_donor_names[1:3], "Control", ifelse(Var1 %in% selected_donor_names[4:8], "IHD", value)))
melt_pct_donor_cluster_annotated <- melt(as.matrix(pct_donor_cluster_annotated[, -(4:5)])) %>% 
  mutate(condition = ifelse(Var1 %in% selected_donor_names[1:3], "Control", ifelse(Var1 %in% selected_donor_names[4:8], "IHD", value)))
red_palette <- c("grey", "#FFCCCC", "#FFB2B2", "#FF9999", "#FF6666", "#FF4D4D", "#FF3333","#FF1A1A", "#FF0000", "#CC0000", "#B30000", "#990000", "#800000", "#660000", "#4D0000", "#330000")

pdf(paste0(Supervised_result_dir, "/Dotplot/Cell_count_per_celltype.pdf"), width = 6, height = 8)
  ggplot(melt_num_per_donor_cluster_annotated, aes(x = Var2, y = as.character(Var1), fill = value)) + geom_tile() + geom_text(aes(label = round(value - 1, 1)), color = "white", size = 6) + 
    scale_fill_gradientn(colors = red_palette, name = "# of cells") + facet_grid(condition ~ ., scales = "free", space = "free") + theme_linedraw() + labs(title = "", x = "Cell type", y = "Donor") + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 16), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()

pdf(paste0(Supervised_result_dir, "/Dotplot/Cell_pct_per_celltype.pdf"), width = 6, height = 8)
  ggplot(melt_pct_donor_cluster_annotated, aes(x = Var2, y = as.character(Var1), fill = value)) + geom_tile() + geom_text(aes(label = round(value, 1)), color = "white", size = 6) + 
    scale_fill_gradientn(colors = red_palette, name = "% of cells") + facet_grid(condition ~ ., scales = "free", space = "free") + theme_linedraw() + labs(title = "", x = "Cell type", y = "Donor") + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 16), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()

##### save object
# saveRDS(control.disease.combined_CM_Fib, paste0(wd_dir, "/results_IHD/0-Seurat_Object/Seurat.obj_sub_clustering_with_annotation_final.RDS"))
