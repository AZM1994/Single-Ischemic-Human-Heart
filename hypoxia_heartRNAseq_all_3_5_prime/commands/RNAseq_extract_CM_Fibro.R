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

color_list_02 <- c("#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#d62728",
                   "#7f7f7f", "#8ca252", "#7f7f7f", "#7f7f7f", "#9467bd", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f", "#7f7f7f")

color_list_regroup <- c("#E97679", "#EE5A9B", "yellow", "#fd8d3c", "#FFB6C1", "#f22c4b", "#FFC34D", "#FF9999", "#A65A34",
                        "#a1d99b", "#657A51", "#66B032", 
                        "dodgerblue", "#9ecae1", "#9c9ede", "#800080", "#8A2BE2", "#e377c2", "#4169E1")

CM_Fibro_color <- c("#f22c4b", "#66B032", "#9c9ede")

##### set working and results directories
wd_dir <- getwd()
result_dir <- paste0(wd_dir, "/results")
dir.create(paste0(result_dir, "/Annotation_extract_CM_Fib"), recursive = T)
##### read sample metadata and gene lists
sample_metadata <- read.csv(paste0(wd_dir, "/data/sample_meta.csv"), header = TRUE)
sample_names <- sample_metadata$sample_name

################################################################################
################################################################################
## extract CM and Fib
################################################################################
################################################################################
# control.disease.combined <- readRDS(paste0(result_dir, "/Seurat.obj_with_annotation_all_celltypes.RDS"))
# all_cell_type <- unique(control.disease.combined$annotated_clusters)
# control.disease.combined_CM_Fib <- subset(control.disease.combined, subset = annotated_clusters %in% c("Cardiomyocyte", "Fibroblast"))
# unique(control.disease.combined_CM_Fib$annotated_clusters)
# control.disease.combined_CM_Fib_test <- RunPCA(control.disease.combined_CM_Fib)
# ElbowPlot(control.disease.combined_CM_Fib_test, ndims = 50, reduction = "pca")

# DefaultAssay(control.disease.combined_CM_Fib) <- 'RNA'
# control.disease.combined_CM_Fib <- RunPCA(control.disease.combined_CM_Fib, npcs = 15, verbose = TRUE) %>%
#   # RunHarmony(group.by.vars = c("samples", "technique"), lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>%
#   # RunHarmony(group.by.vars = "samples", lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>%
#   # RunTSNE(reduction = "harmony", dims = 1:15, verbose = TRUE) %>%
#   RunUMAP(reduction = "harmony", dims = 1:15, verbose = TRUE) %>%
#   FindNeighbors(reduction = "harmony", dims = 1:15, verbose = TRUE) %>%
#   FindClusters(resolution = c(0.8, 1.0, 1.2, 1.4, 1.6), verbose = TRUE)

# saveRDS(control.disease.combined_CM_Fib, paste0(result_dir, "/Seurat.obj_extractCM_Fibro_ready_for_subclustering.RDS"))
control.disease.combined_CM_Fib  <- readRDS(paste0(result_dir, "/Seurat.obj_extractCM_Fibro_ready_for_subclustering.RDS"))

# pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_clustering/clustree.pdf"), width = 10, height = 8)
#   clustree(control.disease.combined_CM_Fib)
# dev.off()

selected_sample_names <- unique(control.disease.combined_CM_Fib$samples)
selected_donor_names <- unique(control.disease.combined_CM_Fib$donor_condition)

RNA_resolution <- "RNA_snn_res.1"

Idents(object = control.disease.combined_CM_Fib) <- RNA_resolution
unique(control.disease.combined_CM_Fib$RNA_snn_res.1)
# table_stat <- t(table(control.disease.combined_CM_Fib$RNA_snn_res.1, control.disease.combined_CM_Fib$samples))
# table_stat <- table_stat[unique(control.disease.combined_CM_Fib$samples), ]
# level_list <- c(0:18)
# control.disease.combined_CM_Fib$RNA_snn_res.1 <- factor(control.disease.combined_CM_Fib$RNA_snn_res.1, levels = level_list)

control.disease.combined_CM_Fib$regroup <- paste0(
  ifelse(control.disease.combined_CM_Fib$RNA_snn_res.1 %in% c("1", "3", "5", "6", "7","12", "15", "17", "18"), "CM", 
         ifelse(control.disease.combined_CM_Fib$RNA_snn_res.1 %in% c("9", "11", "14"), "CM_to_Fib", 
                "Fib")), control.disease.combined_CM_Fib$RNA_snn_res.1)

control.disease.combined_CM_Fib$regroup_02 <- ifelse(control.disease.combined_CM_Fib$regroup == "CM1", "1", ifelse(control.disease.combined_CM_Fib$regroup == "CM12", "4", 
                                              ifelse(control.disease.combined_CM_Fib$regroup == "CM15", "8", ifelse(control.disease.combined_CM_Fib$regroup == "CM17", "7", 
                                              ifelse(control.disease.combined_CM_Fib$regroup == "CM18", "0", ifelse(control.disease.combined_CM_Fib$regroup == "CM3", "5", 
                                              ifelse(control.disease.combined_CM_Fib$regroup == "CM5", "2", ifelse(control.disease.combined_CM_Fib$regroup == "CM6", "3", 
                                              ifelse(control.disease.combined_CM_Fib$regroup == "CM7", "6", ifelse(control.disease.combined_CM_Fib$regroup == "CM_to_Fib11", "9", 
                                              ifelse(control.disease.combined_CM_Fib$regroup == "CM_to_Fib14", "10", ifelse(control.disease.combined_CM_Fib$regroup == "CM_to_Fib9", "11", 
                                              ifelse(control.disease.combined_CM_Fib$regroup == "Fib0", "17", ifelse(control.disease.combined_CM_Fib$regroup == "Fib10", "12", 
                                              ifelse(control.disease.combined_CM_Fib$regroup == "Fib13", "15", ifelse(control.disease.combined_CM_Fib$regroup == "Fib16", "18", 
                                              ifelse(control.disease.combined_CM_Fib$regroup == "Fib2", "13", ifelse(control.disease.combined_CM_Fib$regroup == "Fib4", "14", 
                                              ifelse(control.disease.combined_CM_Fib$regroup == "Fib8", "16", "NA")))))))))))))))))))

level_list <- c(0:18)
control.disease.combined_CM_Fib$regroup_02 <- factor(control.disease.combined_CM_Fib$regroup_02, levels = level_list)

dir.create(paste0(result_dir, "/Annotation_extract_CM_Fib/Unannotated_clustering"))
pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Unannotated_clustering/Unannotated_clustering_split_sample.pdf"), width = 50, height = 4)
  DimPlot(control.disease.combined_CM_Fib, reduction = "umap", split.by = "samples", label = TRUE, label.size = 3, cols = color_list, order = T) + 
    facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x") + theme_linedraw() + ggtitle(NULL) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Unannotated_clustering/Unannotated_clustering_split_donor.pdf"), width = 28, height = 4)
  DimPlot(control.disease.combined_CM_Fib, reduction = "umap", split.by = "donor_condition", label = TRUE, label.size = 3, cols = color_list, order = T) + 
    facet_grid(. ~  factor(donor_condition, level = selected_donor_names), scale = "free_x") + theme_linedraw() + ggtitle(NULL) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()


pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Unannotated_clustering/Unannotated_clustering_split_condition.pdf"), width = 10, height = 5)
  DimPlot(control.disease.combined_CM_Fib, reduction = "umap", group.by = RNA_resolution, split.by = "condition", label = TRUE, label.size = 4, cols = color_list, order = T) + 
    scale_color_manual(values = color_list, breaks = as.character(0:18)) +
    theme_linedraw() + ggtitle(NULL) + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 18))
dev.off()

pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Unannotated_clustering/Unannotated_clustering_split_condition_highlight_9_11_14.pdf"), width = 10, height = 5)
  DimPlot(control.disease.combined_CM_Fib, reduction = "umap", group.by = "regroup_02", split.by = "condition", label = TRUE, label.size = 4, cols = color_list_regroup, order = T) + 
    theme_linedraw() + ggtitle(NULL) + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                             panel.border = element_rect(size = 0.5), text = element_text(size = 24))
  # DimPlot(control.disease.combined_CM_Fib, reduction = "umap", group.by = RNA_resolution, split.by = "condition", label = TRUE, label.size = 4, cols = color_list_02, order = T) + 
  #   theme_linedraw() + ggtitle(NULL) + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  #                            panel.border = element_rect(size = 0.5), text = element_text(size = 24))
dev.off()

################################################################################
################################################################################
### dotplot for marker genes
marker_set1 <- c("TNNT2", "ACTN2", "MYOZ2", "MYPN", "MLIP") # Cardiomyocyte
marker_set2 <- c("DCN", "APOD", "BICC1", "ABCA6") # Fibroblast
all_markers <- c(marker_set1, marker_set2)
marker_groups <- c(rep("Cardiomyocyte markers", length(marker_set1)), rep("Fibroblast markers", length(marker_set2)))

gene_group_df <- data.frame(gene = all_markers, group = marker_groups) %>% mutate(group = factor(group, level = c("Cardiomyocyte markers", "Fibroblast markers")))

dir.create(paste0(result_dir, "/Annotation_extract_CM_Fib/Unannotated_dotplot"))
control.disease.combined_CM_Fib_Control <- subset(control.disease.combined_CM_Fib, subset = condition == "Control")
control.disease.combined_CM_Fib_Disease <- subset(control.disease.combined_CM_Fib, subset = condition == "IHD")

dotplot_all <- DotPlot(control.disease.combined_CM_Fib, features = all_markers, group.by = "regroup_02")
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all <- dotplot_all + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + scale_color_gradient(low = "white", high = "black") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "black", linetype = "dashed"), text = element_text(face = "bold", size = 12)) + RotatedAxis()
ggsave(paste0(result_dir, "/Annotation_extract_CM_Fib/Unannotated_dotplot/Unannotated_dotplot_all.pdf"), plot = dotplot_all, width = 10, height = 6, dpi = 600)

rm(dotplot_Control)
dotplot_Control <- DotPlot(control.disease.combined_CM_Fib_Control, features = all_markers, group.by = "regroup_02")
dotplot_Control$data$gene_group <- factor(gene_group_df$group[match(dotplot_Control$data$features.plot, gene_group_df$gene)])
dotplot_Control$data$condition <- "Control"
dotplot_Control <- dotplot_Control + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + scale_color_gradient(low = "white", high = "dodgerblue3") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "black"), panel.grid.major = element_blank(), 
                           # panel.grid.major = element_line(color = "dodgerblue3", linetype = "dashed", size = 0.15),
                           text = element_text(size = 12), legend.background = element_rect(fill = "white", color = NA), 
                           panel.border = element_rect(size = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  geom_tile(data = subset(dotplot_Control$data, id %in% c("9", "10", "11")), fill = "yellow", color = "#FFFF004D", alpha = 0.65) + 
  geom_point(aes(color = avg.exp.scaled, size = pct.exp), show.legend = F) + guides(colour = guide_colourbar()) + 
  guides(size = guide_legend(override.aes = list(color = "white", fill = "white")))
ggsave(paste0(result_dir, "/Annotation_extract_CM_Fib/Unannotated_dotplot/Unannotated_dotplot_Control.pdf"), plot = dotplot_Control, width = 10, height = 6, dpi = 600)

rm(dotplot_Disease)
dotplot_Disease <- DotPlot(control.disease.combined_CM_Fib_Disease, features = all_markers, group.by = "regroup_02")
dotplot_Disease$data$gene_group <- factor(gene_group_df$group[match(dotplot_Disease$data$features.plot, gene_group_df$gene)])
dotplot_Disease$data$condition <- "IHD"
dotplot_Disease <- dotplot_Disease + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + scale_color_gradient(low = "white", high = "firebrick3") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "black"), panel.grid.major = element_blank(), 
                           # panel.grid.major = element_line(color = "firebrick3", linetype = "dashed", size = 0.15),
                           text = element_text(size = 12), panel.border = element_rect(size = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  geom_tile(data = subset(dotplot_Control$data, id %in% c("9", "10", "11")), fill = "yellow", color = "#FFFF004D", alpha = 0.65) + 
  geom_point(aes(color = avg.exp.scaled, size = pct.exp)) + guides(colour = guide_colourbar()) + 
  guides(size = guide_legend(override.aes = list(color = "white", fill = "white")))
ggsave(paste0(result_dir, "/Annotation_extract_CM_Fib/Unannotated_dotplot/Unannotated_dotplot_Disease.pdf"), plot = dotplot_Disease, width = 10, height = 6, dpi = 600)

####################################################################
##### plot metagenes
# top10_FC_markers_each_cluster_CM <- FindAllMarkers(control.disease.combined_CM_Fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.50)
# top10_FC_markers_each_cluster_CM <- top10_FC_markers_each_cluster_CM %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)
##### numebr of cells in each cluster
min_genes <- 200  # Minimum number of genes expressed
max_genes <- 3000 # Maximum number of genes expressed
min_counts <- 500  # Minimum total counts (optional)
control.disease.combined_CM_Fib[["nGene"]] <- Matrix::colSums(GetAssayData(control.disease.combined_CM_Fib, slot = "counts") > 0)
control.disease.combined_CM_Fib[["nCount"]] <- Matrix::colSums(GetAssayData(control.disease.combined_CM_Fib, slot = "counts"))
control.disease.combined_CM_Fib <- subset(control.disease.combined_CM_Fib, subset = nGene > min_genes & nGene < max_genes & nCount > min_counts)

# saveRDS(control.disease.combined_CM_Fib, paste0(result_dir, "/Seurat.obj_with_annotation_extractCM_Fibro_filtered.RDS"))
# control.disease.combined_CM_Fib <- readRDS(paste0(result_dir, "/Seurat.obj_with_annotation_extractCM_Fibro_filtered.RDS"))

num_per_sample_cluster_unannotated <- table(control.disease.combined_CM_Fib$donor_condition, control.disease.combined_CM_Fib$RNA_snn_res.1)
num_per_sample_cluster_unannotated <- num_per_sample_cluster_unannotated[unique(control.disease.combined_CM_Fib$donor_condition), ]

num_per_sample_cluster_annotated <- table(control.disease.combined_CM_Fib$donor_condition, control.disease.combined_CM_Fib$annotated_clusters)
num_per_sample_cluster_annotated <- num_per_sample_cluster_annotated[unique(control.disease.combined_CM_Fib$donor_condition), ]

cardiomyocyte_marker_list <- c("TNNT2", "ACTN2", "MYOZ2", "MYPN", "MLIP")
fibroblast_marker_list <- c("DCN", "MMP2", "LUM")
hypoxia_marker_list <- c("JUN", "NR4A1", "ZEB1", "ZEB2")
collagen_marker_list <- c("COL4A2", "COL6A1", "COL6A2", "COL21A1")
# fibrotic_marker_list <- c("ACTA2", "TGFB1", "NFKB1", "SMAD1", "SMAD5", "SMAD9", "SMAD6", "SMAD3", "SMAD2", "SMAD7", "SMAD4", "TNF", "CCN2", "ECM1", "ECM2", "DCN", "VIM", "PDGFRA")
fibrotic_marker_list <- c("COL3A1", "S100A4", "PDGFRB")
MMR_pathway_marker_list <- c("MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2")
BER_pathway_marker_list <- c("APEX1", "FEN1", "LIG1", "LIG3", "MBD4", "NEIL1", "NEIL12", "OGG1", "POLB", "PNKP")
NHEJ_pathway_marker_list <- c("XRCC5", "XRCC6", "PRKDC", "LIG4", "DCLRE1C", "NHEJ1", "XRCC4", "PAXX", "MRE11", "TP53BP1", "SHLD1", "SHLD2", "SHLD3")
# ddr_genes <- c("ATM", "ATR", "BRCA1", "RAD50", "RAD52", "RAD54", "TP53", "XRCC5", "XRCC6")

##### check the marker expression in cluster
# gene_list <- rownames(control.disease.combined_CM_Fib)
# ddr_genes %in% gene_list
# list_to_check <- ddr_genes
# VlnPlot(control.disease.combined_CM_Fib, features = list_to_check, idents = c(9, 11), split.by = "condition", ncol = 5)

all_marker_genes <- list()
all_marker_genes[["metagene_cardiomyocyte"]] <- cardiomyocyte_marker_list
all_marker_genes[["metagene_fibroblast"]] <- fibroblast_marker_list
all_marker_genes[["metagene_hypoxia"]] <- hypoxia_marker_list
all_marker_genes[["metagene_collagen"]] <- collagen_marker_list
all_marker_genes[["metagene_fibrotic"]] <- fibrotic_marker_list
# all_marker_genes[["metagene_MMR"]] <- MMR_pathway_marker_list
# all_marker_genes[["metagene_BER"]] <- BER_pathway_marker_list
# all_marker_genes[["metagene_NHEJ"]] <- NHEJ_pathway_marker_list

for (name in names(all_marker_genes)){
  print(name)
  average_expr <- colMeans(as.data.frame(GetAssayData(control.disease.combined_CM_Fib[all_marker_genes[[name]], ])))
  control.disease.combined_CM_Fib@meta.data[, name] <- average_expr
}

metadata_CM_Fib <- control.disease.combined_CM_Fib@meta.data

# all_metagene_grouped_avg <- metadata_CM_Fib %>% group_by(RNA_snn_res.1, condition, samples) %>% 
#   summarize(cardiomyocyte_metagene = mean(metagene_cardiomyocyte, na.rm = TRUE), fibroblast_metagene = mean(metagene_fibroblast, na.rm = TRUE), 
#             hypoxia_metagene = mean(metagene_hypoxia, na.rm = TRUE), collagen_metagene = mean(metagene_collagen, na.rm = TRUE), 
#             fibrotic_metagene = mean(metagene_fibrotic, na.rm = TRUE), MMR_metagene = mean(metagene_MMR, na.rm = TRUE), 
#             BER_metagene = mean(metagene_BER, na.rm = TRUE), NHEJ_metagene = mean(metagene_NHEJ, na.rm = TRUE))
# 
# all_metagene_grouped_avg_melt <- all_metagene_grouped_avg %>% filter(RNA_snn_res.1 %in% c(9, 11)) %>%
#   melt(id.vars = c("RNA_snn_res.1", "condition", "samples"),
#        measure.vars = c("cardiomyocyte_metagene", "fibroblast_metagene", "hypoxia_metagene", "collagen_metagene", 
#                         "fibrotic_metagene", "MMR_metagene", "BER_metagene", "NHEJ_metagene"),
#        variable.name = "metagenes", value.name = "avg_expr")

all_metagene_grouped_avg <- metadata_CM_Fib %>% group_by(regroup_02, condition, samples) %>% 
  summarize(Cardiomyocyte = mean(metagene_cardiomyocyte, na.rm = TRUE), Fibroblast = mean(metagene_fibroblast, na.rm = TRUE), 
            Hypoxia = mean(metagene_hypoxia, na.rm = TRUE), Collagen = mean(metagene_collagen, na.rm = TRUE), 
            Fibrotic = mean(metagene_fibrotic, na.rm = TRUE))

all_metagene_grouped_avg_melt <- all_metagene_grouped_avg %>% filter(regroup_02 %in% c(9, 11)) %>%
  melt(id.vars = c("regroup_02", "condition", "samples"),
       measure.vars = c("Cardiomyocyte", "Fibroblast", "Hypoxia", "Collagen", "Fibrotic"),
       variable.name = "metagenes", value.name = "avg_expr") %>% 
  mutate(regroup_02.modified = paste0("Cardiomyocyte subcluster ", regroup_02)) %>% 
  mutate(regroup_02.modified = factor(regroup_02.modified, levels = c("Cardiomyocyte subcluster 9", "Cardiomyocyte subcluster 11")))

pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/metagene_expr_barplot.pdf"), width = 18, height = 8.5)
  ggplot(all_metagene_grouped_avg_melt, aes(x = metagenes, y = avg_expr, color = condition)) + 
    geom_boxplot(aes(fill = condition), size = 0.5, position = position_dodge(), alpha = 0.5, outlier.shape = NA) + 
    geom_jitter(aes(color = condition), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 2) + 
    stat_compare_means(aes(group = condition), label = "p.format", size = 8, label.y = 1.1 * max(all_metagene_grouped_avg_melt$avg_expr, na.rm = TRUE)) +
    stat_compare_means(aes(group = condition), label = "p.signif", size = 10, label.y = 1.02 * max(all_metagene_grouped_avg_melt$avg_expr, na.rm = TRUE), show.legend = FALSE) + 
    scale_fill_manual(values = Ctrl_IHD_color) + scale_color_manual(values = Ctrl_IHD_color) + facet_wrap(. ~  regroup_02.modified) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    labs(x = "", y = "Average metagene expression")
dev.off()

# pdf(paste0(result_dir, "/Annotation/extract_CM/metagene_feature_plot.pdf"), width = 8, height = 6)
# FeaturePlot(control.disease.combined_CM_Fib, features = "metagene_CM", order = T)
# FeaturePlot(control.disease.combined_CM_Fib, features = "metagene_Fibro", order = T)
# 
# control.disease.combined_CM_Fib_Control <- subset(control.disease.combined_CM_Fib, subset = condition == "Control")
# FeaturePlot(control.disease.combined_CM_Fib_Control, features = "metagene_CM", order = T)
# FeaturePlot(control.disease.combined_CM_Fib_Control, features = "metagene_Fibro", order = T)
# 
# control.disease.combined_CM_Fib_Disease <- subset(control.disease.combined_CM_Fib, subset = condition == "IHD")
# FeaturePlot(control.disease.combined_CM_Fib_Disease, features = "metagene_CM", order = T)
# FeaturePlot(control.disease.combined_CM_Fib_Disease, features = "metagene_Fibro", order = T)
# # FeaturePlot(control.disease.combined_CM_Fib, features = "DCN", order = T)
# dev.off()

# annotated_clusters_list_CM <- c("Fibroblast", "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", "Fibroblast", # 0-4
#   "Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", # 5-9
#   "Fibroblast", "Cardiomyocyte", "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", # 10-14
#   "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", "Cardiomyocyte") # 15-18

annotated_clusters_list_CM <- c("Fibroblast", "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", "Fibroblast", # 0-4
                                "Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", "Fibroblast", "Fibrotic-cardiomyocyte", # 5-9
                                "Fibroblast", "Fibrotic-cardiomyocyte", "Cardiomyocyte", "Fibroblast", "Fibrotic-cardiomyocyte", # 10-14
                                "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", "Cardiomyocyte") # 15-18

names(annotated_clusters_list_CM) <- paste0("cluster_", 0:(length(annotated_clusters_list_CM) - 1))
annotated_clusters_CM <- sapply(control.disease.combined_CM_Fib$RNA_snn_res.1, function(x) annotated_clusters_list_CM[paste0("cluster_", x)])
control.disease.combined_CM_Fib@meta.data <- cbind(control.disease.combined_CM_Fib@meta.data, annotated_clusters_CM)
cell_type_order <- c("Cardiomyocyte", "Fibrotic-cardiomyocyte", "Fibroblast")
## Reorder the cell types in the metadata
control.disease.combined_CM_Fib$annotated_clusters_CM <- factor(control.disease.combined_CM_Fib$annotated_clusters_CM, levels = cell_type_order)

## find marker genes in 9, 11, 14 for conditions, celltypes
# Idents(object = control.disease.combined_CM_Fib) <- RNA_resolution
# table_condition_RNA_snn_res.1 <- table(control.disease.combined_CM_Fib$condition, control.disease.combined_CM_Fib$RNA_snn_res.1)

# subset_Control_CM <- subset(control.disease.combined_CM_Fib, subset = condition == "Control" & annotated_clusters_CM == "Cardiomyocyte")
# marker_genes_Control_CM <- FindAllMarkers(subset_Control_CM, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
# cluster9_markers_Control_CM <- marker_genes_Control_CM %>% filter(cluster == 9, p_val_adj < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) %>% 
#   write.csv(paste0(result_dir, "/Annotation_extract_CM_Fib/marker_genes/cluster9_markers_Control_CM.csv"), row.names = TRUE)
# cluster11_markers_Control_CM <- marker_genes_Control_CM %>% filter(cluster == 11, p_val_adj < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) %>% 
#   write.csv(paste0(result_dir, "/Annotation_extract_CM_Fib/marker_genes/cluster11_markers_Control_CM.csv"), row.names = TRUE)
# cluster14_markers_Control_CM <- marker_genes_Control_CM %>% filter(cluster == 14, p_val_adj < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) %>% 
#   write.csv(paste0(result_dir, "/Annotation_extract_CM_Fib/marker_genes/cluster14_markers_Control_CM.csv"), row.names = TRUE)

# subset_IHD_CM <- subset(control.disease.combined_CM_Fib, subset = condition == "IHD" & annotated_clusters_CM == "Cardiomyocyte")
# marker_genes_IHD_CM <- FindAllMarkers(subset_IHD_CM, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
# cluster9_markers_IHD_CM <- marker_genes_IHD_CM %>% filter(cluster == 9, p_val_adj < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) %>% 
#   write.csv(paste0(result_dir, "/Annotation_extract_CM_Fib/marker_genes/cluster9_markers_IHD_CM.csv"), row.names = TRUE)
# cluster11_markers_IHD_CM <- marker_genes_IHD_CM %>% filter(cluster == 11, p_val_adj < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) %>% 
#   write.csv(paste0(result_dir, "/Annotation_extract_CM_Fib/marker_genes/cluster11_markers_IHD_CM.csv"), row.names = TRUE)
# cluster14_markers_IHD_CM <- marker_genes_IHD_CM %>% filter(cluster == 14, p_val_adj < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) %>% 
#   write.csv(paste0(result_dir, "/Annotation_extract_CM_Fib/marker_genes/cluster14_markers_IHD_CM.csv"), row.names = TRUE)

# marker_genes_9_IHD <- FindMarkers(object = control.disease.combined_CM_Fib, ident.1 = 9, subset.ident = annotated_clusters_CM == "Cardiomyocyte" & condition == "IHD", logfc.threshold = 0.25) %>% 
#   # filter(pct.1 > 0.10, pct.2 > 0.10) %>% 
#   arrange(desc(avg_log2FC), p_val_adj) %>% 
#   write.csv(paste0(result_dir, "/Annotation_extract_CM_Fib/marker_genes/marker_genes_9_Control_CM.csv"), row.names = TRUE)
# marker_genes_9_IHD_Fib <- FindMarkers(object = control.disease.combined_CM_Fib, ident.1 = 9, subset.ident = annotated_clusters_CM == "Fibroblast" & condition == "IHD", logfc.threshold = 0.25) %>% 
#   # filter(pct.1 > 0.10, pct.2 > 0.10) %>% 
#   arrange(desc(avg_log2FC), p_val_adj) %>% 
#   write.csv(paste0(result_dir, "/Annotation_extract_CM_Fib/marker_genes/marker_genes_9_Control_Fib.csv"), row.names = TRUE)

# subset_9_Control_CM <- subset(control.disease.combined_CM_Fib, subset = annotated_clusters_CM == "Cardiomyocyte" & condition == "Control" & RNA_snn_res.1 == 9)
# subset_9_IHD_CM <- subset(control.disease.combined_CM_Fib, subset = annotated_clusters_CM == "Cardiomyocyte" & condition == "IHD" & RNA_snn_res.1 == 9)

dir.create(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_clustering"))
pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_clustering/Annotated_clustering.pdf"), width = 8, height = 5.5)
  DimPlot(control.disease.combined_CM_Fib, reduction = 'umap', group.by = "annotated_clusters_CM", label = T, label.size = 5, cols = CM_Fibro_color, order = T) + ggtitle(NULL) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 20))
dev.off()

pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_clustering/Annotated_clustering_split_condition.pdf"), width = 12, height = 5)
  DimPlot(control.disease.combined_CM_Fib, reduction = 'umap', group.by = "annotated_clusters_CM", split.by = "condition", label = T, label.size = 5, cols = CM_Fibro_color, order = T) + ggtitle(NULL) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                             panel.border = element_rect(size = 0.5), text = element_text(size = 24))
dev.off()

pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_clustering/Annotated_clustering_split_sample.pdf"), width = 50, height = 4)
DimPlot(control.disease.combined_CM_Fib, reduction = 'umap', group.by = "annotated_clusters_CM", split.by = "samples", label = T, label.size = 3, cols = CM_Fibro_color, order = T) + 
  facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x") + ggtitle(NULL) + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_clustering/Annotated_clustering_split_sample_donor.pdf"), width = 28, height = 4)
DimPlot(control.disease.combined_CM_Fib, reduction = 'umap', group.by = "annotated_clusters_CM", split.by = "donor_condition", label = T, label.size = 3, cols = CM_Fibro_color, order = T) + 
  facet_grid(. ~  factor(donor_condition, level = selected_donor_names), scale = "free_x") + ggtitle(NULL) + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

dir.create(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_dotplot"))
control.disease.combined_CM_Fib_Control <- subset(control.disease.combined_CM_Fib, subset = condition == "Control")
control.disease.combined_CM_Fib_Disease <- subset(control.disease.combined_CM_Fib, subset = condition == "IHD")

dotplot_all <- DotPlot(control.disease.combined_CM_Fib, features = all_markers, group.by = "annotated_clusters_CM")
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all <- dotplot_all + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + 
  theme_linedraw() + RotatedAxis() + scale_color_gradient(low = "white", high = "black") + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
ggsave(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_dotplot/Annotated_dotplot_all.pdf"), plot = dotplot_all, width = 14, height = 4, dpi = 300)

dotplot_Control <- DotPlot(control.disease.combined_CM_Fib_Control, features = all_markers, group.by = "annotated_clusters_CM")
dotplot_Control$data$gene_group <- factor(gene_group_df$group[match(dotplot_Control$data$features.plot, gene_group_df$gene)])
dotplot_Control$data$condition <- "Control"
dotplot_Control <- dotplot_Control + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + 
  theme_linedraw() + RotatedAxis() + scale_color_gradient(low = "white", high = "dodgerblue4") + guides(colour = guide_colourbar()) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
ggsave(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_dotplot/Annotated_dotplot_Control.pdf"), plot = dotplot_Control, width = 14, height = 4, dpi = 300)

dotplot_Disease <- DotPlot(control.disease.combined_CM_Fib_Disease, features = all_markers, group.by = "annotated_clusters_CM")
dotplot_Disease$data$gene_group <- factor(gene_group_df$group[match(dotplot_Disease$data$features.plot, gene_group_df$gene)])
dotplot_Disease$data$condition <- "IHD"
dotplot_Disease <- dotplot_Disease + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + 
  theme_linedraw() + RotatedAxis() + scale_color_gradient(low = "white", high = "firebrick4") + guides(colour = guide_colourbar()) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
ggsave(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_dotplot/Annotated_dotplot_Disease.pdf"), plot = dotplot_Disease, width = 14, height = 4, dpi = 300)

################################################################################################################################################################
### statistic for unannotated and annotated clustering
################################################################################################################################################################
num_per_donor_cluster_annotated <- table(control.disease.combined_CM_Fib$donor_condition, control.disease.combined_CM_Fib$annotated_clusters_CM)
num_per_donor_cluster_annotated <- num_per_donor_cluster_annotated[unique(control.disease.combined_CM_Fib$donor_condition), ]

num_per_sample_cluster_annotated <- table(control.disease.combined_CM_Fib$samples, control.disease.combined_CM_Fib$annotated_clusters_CM)
num_per_sample_cluster_annotated <- num_per_sample_cluster_annotated[unique(control.disease.combined_CM_Fib$samples), ]

percentage_per_donor_cluster_annotated <- num_per_donor_cluster_annotated %>% 
  as.data.frame.matrix() %>% mutate(Total = rowSums(.)) %>% 
  mutate(across(-Total, ~ round(.x / Total *100 , 5))) %>% select(-Total) %>% 
  mutate(condition = c("Control", "Control", "Control", "IHD", "IHD", "IHD", "IHD", "IHD"))

percentage_per_sample_cluster_annotated <- num_per_sample_cluster_annotated %>% 
  as.data.frame.matrix() %>% mutate(Total = rowSums(.)) %>% 
  mutate(across(-Total, ~ round(.x / Total *100 , 5))) %>% select(-Total) %>% 
  # mutate(condition = c("Control", "Control", "Control", "IHD", "IHD", "IHD", "IHD", "IHD"))
  mutate(condition = c("Control", "Control", "Control", "Control", "Control", "Control", "Control", "IHD", "IHD", "IHD", "IHD", "IHD", "IHD", "IHD")) %>% 
  mutate(ratio = Cardiomyocyte / Fibroblast) %>% mutate(sample = row.names(.))

percentage_per_sample_cluster_annotated_melt <- percentage_per_sample_cluster_annotated[, 1:3] %>% melt() %>% 
  setNames(c("condition", "cluster", "ratio")) %>% mutate(condition = factor(condition, level = c("Control", "IHD"))) %>% mutate(ratio = as.numeric(ratio))

pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_dotplot/Annotated_cluster_prop_comparison.pdf"), width = 8, height = 5)
  barplot_annotated_cluster_ratio <- ggbarplot(percentage_per_sample_cluster_annotated_melt, x = "cluster", y = "ratio", color = "condition", 
                                               label = F, size = 2, lab.nb.digits = 2, add = c("mean_se", "jitter"), position = position_dodge(0.9), palette = c("dodgerblue4", "firebrick4")) +
    stat_compare_means(aes(group = condition), method = "wilcox", label = "p.format", size = 8, label.y = 0.90 * max(percentage_per_sample_cluster_annotated_melt$ratio, na.rm = TRUE)) + 
    stat_compare_means(aes(group = condition), method = "wilcox", label = "p.signif", size = 8, label.y = 0.95 * max(percentage_per_sample_cluster_annotated_melt$ratio, na.rm = TRUE)) + 
    labs(x = "Cell type", y = "percentage of cells") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 20)) + 
    scale_fill_manual(values = Ctrl_IHD_color)
  print(barplot_annotated_cluster_ratio)
dev.off()

df_long <- percentage_per_sample_cluster_annotated %>% select(sample, condition, Cardiomyocyte, Fibroblast) %>%
  pivot_longer(cols = c("Cardiomyocyte", "Fibroblast"), names_to = "Cell_Type", values_to = "Percentage")
wilcox_test <- wilcox.test(ratio ~ condition, data = percentage_per_sample_cluster_annotated)
p_value <- wilcox_test$p.value
significance <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns")))
max_y <- max(df_long$Percentage) + 5

pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_dotplot/Annotated_cluster_prop_comparison_02.pdf"), width = 8, height = 5)
p <- ggplot(df_long, aes(x = sample, y = Percentage, fill = Cell_Type)) + 
  geom_bar(stat = "identity") + facet_wrap(~ condition, scales = "free_x") + labs(y = "Percentage", x = "Sample", fill = "Cell Type") + 
  theme_linedraw() + labs(title = "", x = "", y = "Percentage of celltype") + scale_fill_manual(values = CM_Fibro_color) + 
  scale_y_continuous(limits = c(0, 120), breaks = c(0, 25, 50, 75, 100)) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 16), panel.border = element_rect(size = 0.5), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
print(p)
dev.off()
## Cell Type Composition
melt_num_per_donor_cluster_annotated <- melt(num_per_donor_cluster_annotated) %>% mutate(value = value + 1) %>% 
  mutate(condition = ifelse(Var1 %in% c("936_Ctrl", "5087_Ctrl", "5919_Ctrl"), "Control", 
                            ifelse(Var1 %in% c("604_IHD", "1156_IHD", "5111_IHD", "5364_IHD", "5874_IHD"), "IHD", value)))
melt_percentage_per_donor_cluster_annotated <- melt(as.matrix(percentage_per_donor_cluster_annotated[,-3])) %>% 
  mutate(condition = ifelse(Var1 %in% c("936_Ctrl", "5087_Ctrl", "5919_Ctrl"), "Control", 
                            ifelse(Var1 %in% c("604_IHD", "1156_IHD", "5111_IHD", "5364_IHD", "5874_IHD"), "IHD", value)))

red_palette <- c("grey", "#FFCCCC", "#FFB2B2", "#FF9999", "#FF6666", "#FF4D4D", "#FF3333","#FF1A1A", "#FF0000", "#CC0000",
                 "#B30000", "#990000", "#800000", "#660000", "#4D0000", "#330000")

pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_dotplot/cell_count_in_annotated_cluster.pdf"), width = 10, height = 8)
ggplot(melt_num_per_donor_cluster_annotated, aes(x = Var2, y = Var1, fill = value)) + 
  geom_tile() + geom_text(aes(label = round(value - 1, 1)), color = "white", size = 6) + 
  # scale_fill_gradientn(colors = red_palette, name = "# of cells", breaks = log10(c(100, 1000)), labels = c(100, 1000)) + 
  scale_fill_gradientn(colors = red_palette, name = "# of cells") + 
  facet_grid(condition ~ ., scales = "free", space = "free") + theme_linedraw() + labs(title = "", x = "Cell type", y = "Donor") + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(face = "bold", size = 16), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()

pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/Annotated_dotplot/cell_percentage_in_annotated_cluster.pdf"), width = 10, height = 8)
ggplot(melt_percentage_per_donor_cluster_annotated, aes(x = Var2, y = Var1, fill = value)) + 
  geom_tile() + geom_text(aes(label = round(value, 1)), color = "white", size = 6) + 
  # scale_fill_gradientn(colors = red_palette, name = "# of cells", breaks = log10(c(1, 10, 100, 1000)), labels = c(1, 10, 100, 1000)) + 
  scale_fill_gradientn(colors = red_palette, name = "% of cells") + 
  facet_grid(condition ~ ., scales = "free", space = "free") + theme_linedraw() + labs(title = "", x = "Cell type", y = "Donor") + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(face = "bold", size = 16), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()

## save object
saveRDS(control.disease.combined_CM_Fib, paste0(result_dir, "/Seurat.obj_sub_clustering_with_annotation_final.RDS"))

cardiomyocyte_marker_list <- c("TNNT2", "ACTN2", "MYOZ2", "MYPN", "MLIP")
fibroblast_marker_list <- c("DCN", "MMP2", "LUM")
hypoxia_marker_list <- c("JUN", "NR4A1", "ZEB1", "ZEB2")
collagen_marker_list <- c("COL4A2", "COL6A1", "COL6A2", "COL21A1")
# fibrotic_marker_list <- c("ACTA2", "TGFB1", "NFKB1", "SMAD1", "SMAD5", "SMAD9", "SMAD6", "SMAD3", "SMAD2", "SMAD7", "SMAD4", "TNF", "CCN2", "ECM1", "ECM2", "DCN", "VIM", "PDGFRA")
fibrotic_marker_list <- c("COL3A1", "S100A4", "PDGFRB")
MMR_pathway_marker_list <- c("MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2")
BER_pathway_marker_list <- c("APEX1", "FEN1", "LIG1", "LIG3", "MBD4", "NEIL1", "NEIL12", "OGG1", "POLB", "PNKP")
NHEJ_pathway_marker_list <- c("XRCC5", "XRCC6", "PRKDC", "LIG4", "DCLRE1C", "NHEJ1", "XRCC4", "PAXX", "MRE11", "TP53BP1", "SHLD1", "SHLD2", "SHLD3")
# ddr_genes <- c("ATM", "ATR", "BRCA1", "RAD50", "RAD52", "RAD54", "TP53", "XRCC5", "XRCC6")

##### check the marker expression in cluster
# gene_list <- rownames(control.disease.combined_CM_Fib)
# MMR_pathway_marker_list %in% gene_list
# list_to_check <- MMR_pathway_marker_list
# Idents(object = control.disease.combined_CM_Fib) <- "annotated_clusters_CM"
# VlnPlot(control.disease.combined_CM_Fib, group.by = "annotated_clusters_CM", features = list_to_check, idents = c("Cardiomyocytes"), split.by = "condition", ncol = 5)

all_marker_genes <- list()
all_marker_genes[["metagene_cardiomyocyte"]] <- cardiomyocyte_marker_list
all_marker_genes[["metagene_fibroblast"]] <- fibroblast_marker_list
all_marker_genes[["metagene_hypoxia"]] <- hypoxia_marker_list
all_marker_genes[["metagene_collagen"]] <- collagen_marker_list
all_marker_genes[["metagene_fibrotic"]] <- fibrotic_marker_list
all_marker_genes[["metagene_MMR"]] <- MMR_pathway_marker_list
all_marker_genes[["metagene_BER"]] <- BER_pathway_marker_list
all_marker_genes[["metagene_NHEJ"]] <- NHEJ_pathway_marker_list

for (name in names(all_marker_genes)){
  print(name)
  average_expr <- colMeans(as.data.frame(GetAssayData(control.disease.combined_CM_Fib[all_marker_genes[[name]], ])))
  control.disease.combined_CM_Fib@meta.data[, name] <- average_expr
}

metadata_CM_Fib_02 <- control.disease.combined_CM_Fib@meta.data

all_metagene_grouped_avg_02 <- metadata_CM_Fib_02 %>% group_by(annotated_clusters_CM, condition, samples) %>% 
  summarize(cardiomyocyte_metagene = mean(metagene_cardiomyocyte, na.rm = TRUE), fibroblast_metagene = mean(metagene_fibroblast, na.rm = TRUE), 
            hypoxia_metagene = mean(metagene_hypoxia, na.rm = TRUE), collagen_metagene = mean(metagene_collagen, na.rm = TRUE), 
            fibrotic_metagene = mean(metagene_fibrotic, na.rm = TRUE), MMR_metagene = mean(metagene_MMR, na.rm = TRUE), 
            BER_metagene = mean(metagene_BER, na.rm = TRUE), NHEJ_metagene = mean(metagene_NHEJ, na.rm = TRUE))

all_metagene_grouped_avg_melt_02 <- all_metagene_grouped_avg_02 %>% 
  # filter(RNA_snn_res.1 %in% c(9, 11)) %>%
  melt(id.vars = c("annotated_clusters_CM", "condition", "samples"),
       measure.vars = c("cardiomyocyte_metagene", "fibroblast_metagene", "hypoxia_metagene", "collagen_metagene", 
                        "fibrotic_metagene", "MMR_metagene", "BER_metagene", "NHEJ_metagene"),
       variable.name = "metagenes", value.name = "avg_expr")

pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/metagene_expr_barplot_02.pdf"), width = 18, height = 8)
ggplot(all_metagene_grouped_avg_melt_02, aes(x = metagenes, y = avg_expr, color = condition)) + 
  geom_boxplot(aes(fill = condition), position = position_dodge(), alpha = 0.5, outlier.shape = NA) + 
  geom_jitter(aes(color = condition), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 2) + 
  stat_compare_means(aes(group = condition), label = "p.format", label.y = 1.02 * max(all_metagene_grouped_avg_melt_02$avg_expr, na.rm = TRUE)) + 
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.06 * max(all_metagene_grouped_avg_melt_02$avg_expr, na.rm = TRUE)) + 
  scale_fill_manual(values = Ctrl_IHD_color) + scale_color_manual(values = Ctrl_IHD_color) + facet_wrap(. ~  annotated_clusters_CM) + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           text = element_text(face = "bold", size = 18), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + 
  labs(x = "", y = "Average metagene expression")
dev.off()
