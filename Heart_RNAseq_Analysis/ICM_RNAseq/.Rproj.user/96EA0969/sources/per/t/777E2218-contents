library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(reshape2)
##### set working and results directories
wd_dir <- getwd()
LOO_DEG_CM_dir <- paste0(wd_dir, "/results/4-LOO_DEG_CM")
dir.create(LOO_DEG_CM_dir, recursive = TRUE, showWarnings = FALSE)

##### load annotated rds file
control.disease.combined_CM_Fib <- readRDS(paste0(wd_dir, "/results/0-Seurat_Object/Seurat.obj_sub_clustering_with_annotation_final.RDS"))
GenesonChrXY <- read.csv(paste0(wd_dir, "/data/GenesonChrXY.csv"), header = TRUE)
GenesonChrXY <- GenesonChrXY$SYMBOL

##### select the CM clusters for DEG analysis
# annotated_clusters_list_CM <- c("Fibroblast", "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", "Fibroblast", # 0-4
#                                 "Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", # 5-9
#                                 "Fibroblast", "Cardiomyocyte", "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", # 10-14
#                                 "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", "Cardiomyocyte") # 15-19

Idents(object = control.disease.combined_CM_Fib) <- "RNA_snn_res.1"
control.disease.combined_CM_only <- subset(control.disease.combined_CM_Fib, idents = c(1,3,5,6,7,9,11,12,14,15,17,18))
# num_per_sample_cluster_unannotated <- table(control.disease.combined_CM_only$donor_condition, control.disease.combined_CM_only$RNA_snn_res.1)
# num_per_sample_cluster_unannotated <- num_per_sample_cluster_unannotated[unique(control.disease.combined_CM_only$donor_condition), ]

Idents(object = control.disease.combined_CM_only) <- "donor_condition"
control.disease.combined_CM_only <- subset(control.disease.combined_CM_only, idents = c("936_Ctrl", "5087_Ctrl", "5919_Ctrl", "1156_IHD", "5111_IHD", "5874_IHD"))
Ctrl_samples <- c("936_Ctrl",  "5087_Ctrl", "5919_Ctrl")
IHD_samples <- c("1156_IHD", "5111_IHD", "5874_IHD")
# table(control.disease.combined_CM_only$donor_condition)
# saveRDS(control.disease.combined_CM_only, paste0(wd_dir, "/results/0-Seurat_Object/Seurat.obj_sub_clustering_CM_only.RDS"))

##### LOO DEG analysis: Find DEGs for all combinations of leave one out sample combinations
DEG_Celltypes <- unique(control.disease.combined_CM_only$annotated_clusters_CM)

## select min.pct, logfc.threshold for LOO
min.pct_list <- c(0.1, 0.2, 0.3, 0.4, 0.5)
logfc.threshold_up_list <- c(0.1, 0.2, 0.3, 0.4, 0.5)
logfc.threshold_dn_list <- c(0.1, 0.2, 0.3, 0.4, 0.5)
p_val_adj_thr <- 0.05
param_set <- expand.grid(min.pct_thr = min.pct_list, logfc.threshold_up_thr = logfc.threshold_up_list, logfc.threshold_dn_thr = logfc.threshold_dn_list)
param_set_02 <- expand.grid(min.pct_thr = min.pct_list, logfc.threshold_thr = logfc.threshold_up_list)
min.pct_lb <- min(min.pct_list)
logfc.threshold_lb <- min(logfc.threshold_up_list)

## all LOO sample combinations
all_sample_combinations <- expand.grid(Ctrl = Ctrl_samples, IHD = IHD_samples)
LOO_sample_list <- apply(all_sample_combinations, 1, function(comb) {
  remaining_Ctrl <- Ctrl_samples[Ctrl_samples != comb["Ctrl"]]
  remaining_IHD <- IHD_samples[IHD_samples != comb["IHD"]]
  return(list(remaining_samples = c(remaining_Ctrl, remaining_IHD)))
})

LOO_DEG_list <- list()
for (loop_idx in 1:length(LOO_sample_list)) {
  LOO_samples <- unlist(LOO_sample_list[[loop_idx]])
  cat(loop_idx, "Identifying LOO DEGs for samples:", paste(LOO_samples, collapse = ", "), "\n")
  loop_idx <- loop_idx + 1 
  
  ## Subset seurat object by samples
  Idents(object = control.disease.combined_CM_only) <- "donor_condition"
  control.disease.combined_CM_only_LOO <- subset(control.disease.combined_CM_only, idents = LOO_samples)
        
  ## find DEG for each cell type
  DEG_list_idx <- c()
  for (celltype in DEG_Celltypes) {
    cat("Cell type:", as.character(celltype), ", ")
    LOO_idx_celltype <- control.disease.combined_CM_only_LOO[, control.disease.combined_CM_only_LOO$annotated_clusters_CM %in% celltype]
    Idents(LOO_idx_celltype) <- LOO_idx_celltype$condition
    DEG_list_idx_celltype <- FindMarkers(LOO_idx_celltype, ident.1 = "IHD", ident.2 = "Control", min.pct = min.pct_lb, logfc.threshold = logfc.threshold_lb, test.use = "wilcox")
    # latent.vars = "num_cells_per_donor.condition")
    DEG_list_idx_celltype <- DEG_list_idx_celltype %>% filter(p_val_adj < p_val_adj_thr) %>% mutate(Gene = rownames(.), Celltype = celltype)
    cat("Number of DEGs:", nrow(DEG_list_idx_celltype), "\n")
    
    DEG_list_idx <- rbind(DEG_list_idx, DEG_list_idx_celltype)
  }
  LOO_DEG_list[[paste(LOO_samples, collapse = ".")]] <- DEG_list_idx
}
saveRDS(LOO_DEG_list, file = paste0(LOO_DEG_CM_dir, "/LOO_DEG_CM_only.RDS"))

##### Identify the DEGs for different parameter thresholds
LOO_DEG_list <- readRDS(paste0(LOO_DEG_CM_dir, "/LOO_DEG_CM_only.RDS"))
Target_Celltype <- "Cardiomyocyte"
LOO_DEG_table <- do.call(rbind, LOO_DEG_list)
write.csv(LOO_DEG_table, paste0(LOO_DEG_CM_dir, "/LOO_DEG_CM_only.csv"), row.names = TRUE)

LOO_combo_thr <- 9
LOO_DEG_count_table <- c()
# for (param_idx in 1:nrow(param_set_02)) {
for (param_idx in 1) {
  # min.pct_thr = param_set_02[param_idx, 1]
  # logfc.threshold_thr_up = param_set_02[param_idx, 2]
  # logfc.threshold_thr_down = param_set_02[param_idx, 2]
  min.pct_thr <- 0.10
  logfc.threshold_thr_up <- 0.30
  logfc.threshold_thr_down <- -0.25
  cat("Evaluating min.pct_thr:", min.pct_thr, ", logfc.threshold_thd:", logfc.threshold_thr_up, logfc.threshold_thr_down, "...\n")
  
  ## extract all LOO DEGs in Target Celltype and filter with given thresholds
  LOO_DEGs_target_celltype_filtered <- purrr::map_dfr(LOO_DEG_list, ~ .x %>% filter(Celltype %in% Target_Celltype)) %>% 
    mutate(pct_diff = abs(as.numeric(pct.1) - as.numeric(pct.2))) %>% 
    filter(avg_log2FC > logfc.threshold_thr_up | avg_log2FC < logfc.threshold_thr_down, pct.1 > min.pct_thr | pct.2 > min.pct_thr)
  # ggplot(LOO_DEGs_target_celltype[LOO_DEGs_target_celltype$avg_log2FC >=0,], aes(x = celltype, y = avg_log2FC)) + geom_violin(trim = TRUE, na.rm = FALSE)
  # ggplot(LOO_DEGs_target_celltype, aes(x = Celltype, y = pct.2)) + geom_violin(trim = TRUE, na.rm = FALSE)
  
  for (LOO_combo_thr_idx in LOO_combo_thr) {
    ## Up regulated DEGs, ordered by avg_log2FC
    LOO_DEGs_up_ordered <- LOO_DEGs_target_celltype_filtered %>% 
      filter(avg_log2FC > 0) %>% count(Gene) %>% 
      filter(n >= LOO_combo_thr_idx) %>% 
      inner_join(LOO_DEGs_target_celltype_filtered %>% group_by(Gene) %>% summarise(mean_log2FC = mean(avg_log2FC), .groups = 'drop'), by = "Gene") %>% 
      arrange(mean_log2FC) %>% pull(Gene)
    
    # LOO_DEGs_up_ordered_ext <- LOO_DEGs_target_celltype_filtered %>% filter(Gene %in% LOO_DEGs_up_ordered) %>% group_by(Gene) %>% 
    #   summarise(mean_log2FC = mean(avg_log2FC)) %>% arrange(mean_log2FC)
    
    ## Down regulated DEGs, ordered by avg_log2FC
    LOO_DEGs_dn_ordered <- LOO_DEGs_target_celltype_filtered %>% 
      filter(as.numeric(avg_log2FC) < 0) %>% count(Gene) %>% 
      filter(n >= LOO_combo_thr_idx) %>%
      inner_join(LOO_DEGs_target_celltype_filtered %>% group_by(Gene) %>% summarise(mean_log2FC = mean(avg_log2FC), .groups = 'drop'), by = "Gene") %>% 
      arrange(mean_log2FC) %>% pull(Gene)
    
    # LOO_DEGs_dn_ordered_ext <- LOO_DEGs_target_celltype_filtered %>% filter(Gene %in% LOO_DEGs_dn_ordered) %>% group_by(Gene) %>% 
    #   summarise(mean_log2FC = mean(avg_log2FC)) %>% arrange(mean_log2FC)
    
    LOO_DEG_count_df_idx <- c(min.pct_thr, logfc.threshold_thr_up, logfc.threshold_thr_down, LOO_combo_thr_idx, 
                       length(LOO_DEGs_up_ordered), length(LOO_DEGs_dn_ordered), length(LOO_DEGs_up_ordered) + length(LOO_DEGs_dn_ordered))
    LOO_DEG_count_table <- rbind(LOO_DEG_count_table, LOO_DEG_count_df_idx) %>% 
      `rownames<-`(NULL) %>%
      `colnames<-`(c("pct_thr", "log2fc_thr_up", "log2fc_thr_down", "LOO_thr", "Num_of_DEGs_up", "Num_of_DEGs_dn", "Total_DEGs"))
  }
}

## Heatmap plot after identifying the best thresholds combination
control.disease.combined_CM_only$donor_condition <- factor(x = control.disease.combined_CM_only$donor_condition, levels = c(Ctrl_samples, IHD_samples))
LOO_DEGs_up_filtered <- setdiff(LOO_DEGs_up_ordered, GenesonChrXY)
# valid_genes_up <- intersect(LOO_DEGs_up_filtered, rownames(control.disease.combined_CM_only[["RNA"]]@scale.data))
LOO_DEGs_dn_filtered <- setdiff(LOO_DEGs_dn_ordered, GenesonChrXY)
# valid_genes_dn <- intersect(LOO_DEGs_dn_filtered, rownames(control.disease.combined_CM_only[["RNA"]]@scale.data))
LOO_DEGs_all_filtered <- c(LOO_DEGs_up_filtered, LOO_DEGs_dn_filtered)

##### using pheatmap
LOO_DEG_expr_matrix <- as.matrix(GetAssayData(object = control.disease.combined_CM_only, assay = "RNA", slot = "count")[LOO_DEGs_all_filtered, ])
sample_metadata <- control.disease.combined_CM_only@meta.data
LOO_DEG_scaled <- t(apply(LOO_DEG_expr_matrix, 1, scale))
colnames(LOO_DEG_scaled) <- rownames(sample_metadata)

##### count of number of cells expressed LOO DEGs in each donor_condition
# LOO_DEGs_count_donor_condition <- data.frame(donor_condition = c("936_Ctrl", "5087_Ctrl", "5919_Ctrl", "1156_IHD", "5111_IHD", "5874_IHD"))
# for (gene in LOO_DEGs_all_filtered) {
#   expr_data <- GetAssayData(control.disease.combined_CM_only, assay = "RNA", slot = "count")[gene, ]
#   expr_data <- as.matrix(expr_data)
#   cell_counts <- data.frame(cell = rownames(expr_data), expression = expr_data)
#   cell_counts$donor_condition <- sample_metadata$donor_condition
#   cell_counts <- cell_counts[cell_counts$expression != 0, ]
#   count_summary <- as.data.frame(table(cell_counts$donor_condition)) %>% select(-Var1) %>% setNames(gene)
#   LOO_DEGs_count_donor_condition <- cbind(LOO_DEGs_count_donor_condition, count_summary)
# }

##### draw violin plot for gene expression for selected DEGs
# seleted_DEGs_up <- c("EGR1", "ZFP36", "AC023494.1", "CCBE1", "FOS", "BTG2", "DOK6", "NR4A1", "ZNF608", "PDGFD")
# seleted_DEGs_dn <- c("ANKRD2", "SLC35F1", "SYNDIG1", "KCNAB2", "KIF26B", "TMEM120B", "SCN5A")
# expr_data_up <- as.matrix(GetAssayData(object = control.disease.combined_CM_only, assay = "RNA", slot = "count")[LOO_DEGs_up_filtered, ])
# expr_data_down <- as.matrix(GetAssayData(object = control.disease.combined_CM_only, assay = "RNA", slot = "count")[LOO_DEGs_dn_filtered, ])
# expr_data_up <- as.matrix(GetAssayData(object = control.disease.combined_CM_only, assay = "RNA", slot = "count")[seleted_DEGs_up, ]) %>%
#   t() %>% as.data.frame() %>% mutate(donor_condition = sample_metadata$donor_condition)
# expr_data_down <- as.matrix(GetAssayData(object = control.disease.combined_CM_only, assay = "RNA", slot = "count")[seleted_DEGs_dn, ]) %>%
#   t() %>% as.data.frame() %>% mutate(donor_condition = sample_metadata$donor_condition)
# expr_data_up_long <- melt(expr_data_up, id.vars = "donor_condition", variable.name = "gene", value.name = "expression")
# expr_data_down_long <- melt(expr_data_down, id.vars = "donor_condition", variable.name = "gene", value.name = "expression")

# p_DEG_up_violin <- ggplot(expr_data_up_long, aes(x = donor_condition, y = expression, fill = donor_condition)) + 
#   geom_violin(scale = "width", adjust = 1, trim = T) + geom_boxplot(width=0.2, color="white", alpha=0.2) + facet_wrap(~ gene, scales = "free", ncol = 5) + scale_y_log10() +
#   labs(title = "", x = "Donor Condition", y = "Gene Expression (counts)") + theme_linedraw() + ggtitle(NULL) + RotatedAxis() +
#   theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
# ggsave(paste0(LOO_DEG_CM_dir, "/LOO_DEG_Cardiomyocyte/p_DEG_up_violin.pdf"), plot = p_DEG_up_violin, width = 35, height = 25, dpi = 300)

# p_DEG_down_violin <- ggplot(expr_data_down_long, aes(x = donor_condition, y = expression, fill = donor_condition)) +
#   geom_violin(scale = "width", adjust = 1, trim = T) + geom_boxplot(width=0.2, color="white", alpha=0.2) + facet_wrap(~ gene, scales = "free", ncol = 5) + scale_y_log10() +
#   labs(title = "", x = "Donor Condition", y = "Gene Expression (counts)") + theme_linedraw() + ggtitle(NULL) + RotatedAxis() +
#   theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
# ggsave(paste0(LOO_DEG_CM_dir, "/LOO_DEG_Cardiomyocyte/p_DEG_down_violin.pdf"), plot = p_DEG_down_violin, width = 35, height = 25, dpi = 300)

# FeaturePlot(control.disease.combined_CM_only, features = "PIK3R1", split.by = "donor_condition", ncol = 4, max.cutoff = NA, cols = c("grey", "red"))
# FeaturePlot(control.disease.combined_CM_only, features = "GALNT17", split.by = "donor_condition", ncol = 4, max.cutoff = NA, cols = c("grey", "red"))

##### Shuffle cells within each sample
sample_info <- sapply(strsplit(colnames(LOO_DEG_scaled), "_"), function(x) {paste0(x[1], "_", x[3])})
LOO_DEG_scaled <- do.call(cbind, lapply(unique(sample_info), function(sample) {
  sample_cols <- LOO_DEG_scaled[, sample_info == sample]
  sample_cols_shuffled <- sample_cols[, sample(1:ncol(sample_cols))]
  return(sample_cols_shuffled)}))

annotation_cols <- data.frame(Condition = sample_metadata$condition, Donor = sample_metadata$donor_condition, row.names = rownames(sample_metadata))
annotation_rows <- data.frame(DEG.Group = factor(rep(c("Higher in IHD", "Higher in Control"), c(length(LOO_DEGs_up_filtered), length(LOO_DEGs_dn_filtered)))), row.names = LOO_DEGs_all_filtered)
annotation_colors <- list(Condition = c(Control = "dodgerblue2", IHD = "firebrick2"), 
                          Donor = c("936_Ctrl" = "#ADD8E6", "5087_Ctrl" = "#2D6EA9", "5919_Ctrl" = "#00008B", "1156_IHD" = "#FFCCCC", "5111_IHD" = "#DD555D", "5874_IHD" = "#660000"), 
                          DEG.Group = c("Higher in IHD" = "#66A61E", "Higher in Control" = "#e377c2"))

# gaps <- cumsum(table(sample_info))
# heatmap_colors <- colorRampPalette(c("purple", "black", "yellow"))(100)
# breaks <- seq(-2.5, 2.5, length.out = 101)

purple_colors <- colorRampPalette(c("purple2", "purple4"))(2)
yellow_colors <- colorRampPalette(c("yellow3", "yellow1"))(3)
heatmap_colors <- colorRampPalette(c(purple_colors, "black", yellow_colors))(100)
breaks <- seq(-2, 3, length.out = 101)

p <- pheatmap(mat = LOO_DEG_scaled, color = heatmap_colors, annotation_col = annotation_cols, annotation_row = annotation_rows, annotation_colors = annotation_colors, 
  cluster_rows = F, cluster_cols = F, scale = "none", fontsize_row = 8, fontsize_col = 8, border_color = "white", 
  show_colnames = F, show_rownames = F, legend = TRUE, annotation_legend = TRUE, gaps_col = NULL, gaps_row = NULL, breaks = breaks, useRaster=F)
ggsave(paste0(LOO_DEG_CM_dir, "/LOO_DEG_pheatmap.png"), plot = p, width = 15, height = 6, dpi = 600)

LOO_DEG_up_down_df <- as.data.frame(rbind(cbind(LOO_DEGs_up_filtered, "up"), cbind(LOO_DEGs_dn_filtered, "down"))) %>% setNames(c("Gene", "regulation"))
LOO_DEG_up_down_df_ext <- LOO_DEGs_target_celltype_filtered %>% filter(Gene %in% LOO_DEGs_all_filtered) %>% 
  group_by(Gene) %>% 
  summarise(p_val = mean(p_val), avg_log2FC = mean(avg_log2FC), pct.1 = mean(pct.1), pct.2 = mean(pct.2), 
            p_val_adj = mean(p_val_adj), pct_diff = mean(pct_diff), Celltype = first(Celltype)) %>% 
  left_join(LOO_DEG_up_down_df, by = "Gene") %>% arrange(desc(avg_log2FC), p_val_adj) %>% ungroup()
write.csv(LOO_DEG_up_down_df_ext, paste0(LOO_DEG_CM_dir, "/LOO_DEG_up_down_df.csv"))
