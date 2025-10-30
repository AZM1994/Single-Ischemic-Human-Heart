library(Seurat)
library(patchwork)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(scCustomize)
library(harmony)
library(sctransform)
library(glmGamPoi)
library(pheatmap)
library(monocle3)
library(SeuratWrappers)
library(viridisLite)
library(scales)

Ctrl_IHD_color <- c("#2D6EA8", "#DD555B")
color_list_regroup <- c("#E97679", "#EE5A9B", "yellow", "#fd8d3c", "#FFB6C1", "#f22c4b", "#FFC34D", "#FF9999", "#A65A34",
                        "#a1d99b", "#657A51", "#66B032", 
                        "dodgerblue", "#9ecae1", "#9c9ede", "#800080", "#8A2BE2", "#e377c2", "#4169E1")
CM_Fibro_color <- c("#f22c4b", "#66B032", "#9c9ede")

##### set working and results directories
wd_dir <- getwd()
result_dir <- paste0(wd_dir, "/results_ICM/3-CMFib_Trajectory")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

##### load CM and Fib Seurat object
control.disease.combined_CM_Fib  <- readRDS(paste0(wd_dir, "/results_ICM/0-Seurat_Object/Seurat.obj_sub_clustering_with_annotation_final.RDS"))
selected_sample_names <- unique(control.disease.combined_CM_Fib$samples)
selected_donor_names <- unique(control.disease.combined_CM_Fib$donor_condition)

## tranditional way
cds <- as.cell_data_set(control.disease.combined_CM_Fib)
if (is.null(rowData(cds)$gene_short_name)) {
  rowData(cds)$gene_short_name <- row.names(rowData(cds))
}

reducedDims(cds)$UMAP <- control.disease.combined_CM_Fib@reductions$umap@cell.embeddings
# cds <- cluster_cells(cds)
# cds <- learn_graph(cds, use_partition = F, learn_graph_control = list(minimal_branch_len = 5))

# meta <- as.data.frame(colData(cds))
# root_cells <- rownames(meta %>% filter(annotated_clusters_CM == "Cardiomyocyte") %>% arrange(desc(Cardiomyocyte_metagene1)) %>% head(20))
# cds <- order_cells(cds, root_cells = root_cells)
# 
# p_pseudotime_all <- plot_cells(
#   cds, color_cells_by = "pseudotime", label_groups_by_cluster = F, trajectory_graph_color = "green", label_leaves = F,
#   label_branch_points = F, label_roots = T, graph_label_size = 4, label_cell_groups = FALSE, show_trajectory_graph = TRUE,
#   trajectory_graph_segment_size = 1.25, cell_size = 0.5) + scale_color_viridis_c(option = "magma", limits = c(0, 32), direction = -1)
# ggsave(paste0(result_dir, "/p_pseudotime_classic.pdf"), plot = p_pseudotime_all, width = 10, height = 8, dpi = 300)

##### refined approach
##### Condition-stratified pseudotime analysis of cardiomyocyte-to-fibroblast transition
target_cm <- 5000 # healthy CM to keep
target_fb <- 3000 # mature fibro to keep
k_conn_cm    <- 200 # extra CM closest to FCM for connectivity
k_conn_fb    <- 800 # extra FB closest to FCM for connectivity
tip_n     <- 100 # far-end tips for context

plot_condition_trajectory <- function(cds, condition_label, keep_all) {
  
  cds_cond <- cds[, colData(cds)$condition == condition_label]
  emb <- as.data.frame(reducedDims(cds_cond)$UMAP); colnames(emb) <- c("UMAP_1","UMAP_2")
  md  <- as.data.frame(colData(cds_cond))
  lab <- md$annotated_clusters_CM
  fcm_cells <- rownames(md)[lab == "Fibrotic-cardiomyocyte"]
  cm_cells  <- rownames(md)[lab == "Cardiomyocyte"]
  fb_cells  <- rownames(md)[lab == "Fibroblast"]
  if (keep_all == TRUE){
    target_cm <- length(cm_cells)
    target_fb <- length(fb_cells)
  }
  
  # z-scores for CM and Fib
  z <- function(x) as.numeric(scale(x))
  cm_z  <- setNames(z(md$Cardiomyocyte_metagene), rownames(md))
  fib_z <- setNames(z(md$Fibrotic_metagene),      rownames(md))
  col_z <- if ("Collagen_metagene" %in% names(md)) setNames(z(md$Collagen_metagene), rownames(md)) else setNames(rep(0, nrow(md)), rownames(md))
  
  cm_score <-  1.0*cm_z + (-0.5)*fib_z + (-0.5)*col_z
  fb_score <- (-0.5)*cm_z +  1.0*fib_z +  0.5*col_z
  
  # pick top by score
  cm_keep <- names(sort(cm_score[cm_cells], decreasing = TRUE))[seq_len(min(target_cm, length(cm_cells)))]
  fb_keep <- names(sort(fb_score[fb_cells], decreasing = TRUE))[seq_len(min(target_fb, length(fb_cells)))]
  
  # connectors
  fcm_centroid <- colMeans(emb[fcm_cells, , drop = FALSE])
  d_cm <- rowSums((emb[cm_cells, , drop = FALSE] - fcm_centroid)^2)
  d_fb <- rowSums((emb[fb_cells, , drop = FALSE] - fcm_centroid)^2)
  cm_conn <- cm_cells[order(d_cm)][seq_len(min(k_conn_cm, length(cm_cells)))]
  fb_conn <- fb_cells[order(d_fb)][seq_len(min(k_conn_fb, length(fb_cells)))]
  
  # tips
  cm_tips <- cm_cells[order(d_cm, decreasing = TRUE)][seq_len(min(tip_n, length(cm_cells)))]
  fb_tips <- fb_cells[order(d_fb, decreasing = TRUE)][seq_len(min(tip_n, length(fb_cells)))]
  
  # final subset
  cells_final <- unique(c(fcm_cells, cm_keep, fb_keep, cm_conn, fb_conn, cm_tips, fb_tips))
  cds_sub <- cds_cond[, cells_final]
  
  # re-learn graph
  cds_sub <- cluster_cells(cds_sub)
  cds_sub <- learn_graph(cds_sub, use_partition = FALSE, close_loop = FALSE,
                         learn_graph_control = list(prune_graph = TRUE, ncenter = 500, minimal_branch_len = 6))
  
  # root cell
  top_cm50 <- names(sort(cm_score[cm_cells], decreasing = TRUE))[seq_len(min(20, length(cm_cells)))]
  centroid <- colMeans(emb[top_cm50, , drop = FALSE])
  d2 <- rowSums((emb[top_cm50, , drop = FALSE] - matrix(centroid, nrow = length(top_cm50), ncol = 2, byrow = TRUE))^2)
  root_cell <- top_cm50[which.min(d2)]
  cds_sub <- order_cells(cds_sub, root_cells = root_cell)
  
  magma_shift <- magma(256, direction = -1)
  plot_cells(cds_sub, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE, trajectory_graph_color = "#00FFFF", label_leaves = FALSE, 
             label_branch_points = FALSE, label_roots = TRUE, graph_label_size = 3, label_cell_groups = FALSE, show_trajectory_graph = TRUE, 
             trajectory_graph_segment_size = 1, cell_size = 0.1) + facet_grid(. ~  condition, scale = "free_x") + ggtitle(NULL) + 
    scale_x_continuous(limits = c(-8, 7), breaks = c(-5, 0, 5)) + scale_y_continuous(limits = c(-4, 5), breaks = c(-2.5, 0, 2.5, 5)) + theme_linedraw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), , panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24)) + 
    scale_color_gradientn(colors = magma_shift, limits = c(0, 35), values = rescale(c(0, 15, 35)), oob = scales::squish) + labs(color = "Pseudotime")
}

## Run for both conditions
p_ctrl <- plot_condition_trajectory(cds, "Control", keep_all = TRUE) + guides(color = "none")
p_ihd  <- plot_condition_trajectory(cds, "ICM", keep_all = TRUE)
p_cond <- p_ctrl + p_ihd
ggsave(paste0(result_dir, "/p_pseudotime_cond.pdf"), plot = p_cond, width = 12, height = 5, dpi = 300)