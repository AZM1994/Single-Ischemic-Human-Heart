library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(reshape2)
##### set working and results directories
wd_dir <- getwd()
result_dir <- paste0(wd_dir, "/results")
dir.create(paste0(result_dir, "/LOO_DEG_Cardiomyocyte"), recursive = T)
GenesonChrXY <- read.csv(paste0(wd_dir, "/data/GenesonChrXY.csv"), header = TRUE)
GenesonChrXY <- GenesonChrXY$SYMBOL

##### load annotated rds file
control.disease.combined_CM_Fib <- readRDS(paste0(result_dir, "/Seurat.obj_sub_clustering_with_annotation_final.RDS"))

num_per_sample_cluster_unannotated <- table(control.disease.combined_CM_Fib$donor_condition, control.disease.combined_CM_Fib$RNA_snn_res.1)
num_per_sample_cluster_unannotated <- num_per_sample_cluster_unannotated[unique(control.disease.combined_CM_Fib$donor_condition), ]
# annotated_clusters_list_CM <- c("Fibroblast", "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", "Fibroblast", # 0-4
#                                 "Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", # 5-9
#                                 "Fibroblast", "Cardiomyocyte", "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", # 10-14
#                                 "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", "Cardiomyocyte") # 15-19


Idents(object = control.disease.combined_CM_Fib) <- "RNA_snn_res.1"
# control.disease.combined_CM_only <- control.disease.combined_CM_Fib
control.disease.combined_CM_only <- subset(control.disease.combined_CM_Fib, idents = c(1,3,5,6,7,9,11,12,14,15,17,18))
# control.disease.combined_CM_only <- subset(control.disease.combined_CM_Fib, idents = c(1,3,5,6,7,9,11,12))
# control.disease.combined_CM_only <- subset(control.disease.combined_CM_Fib, idents = c(1,3,5,6,7))

num_per_sample_cluster_unannotated <- table(control.disease.combined_CM_only$donor_condition, control.disease.combined_CM_only$RNA_snn_res.1)
num_per_sample_cluster_unannotated <- num_per_sample_cluster_unannotated[unique(control.disease.combined_CM_only$donor_condition), ]

Idents(object = control.disease.combined_CM_only) <- "donor_condition"
control.disease.combined_CM_only <- subset(control.disease.combined_CM_only, idents = c("936_Ctrl", "5087_Ctrl", "5919_Ctrl", "1156_IHD", "5111_IHD", "5874_IHD"))

Control_samples <- c("936_Ctrl",  "5087_Ctrl", "5919_Ctrl")
# IHD_samples <- c("604_IHD", "1156_IHD", "5111_IHD", "5364_IHD", "5874_IHD")
# Control_samples <- c("5087_Ctrl", "5919_Ctrl")
IHD_samples <- c("1156_IHD", "5111_IHD", "5874_IHD")

table(control.disease.combined_CM_only$donor_condition)

celltypes <- unique(control.disease.combined_CM_only$annotated_clusters_CM)
# saveRDS(control.disease.combined_CM_only, paste0(result_dir, "/Seurat.obj_sub_clustering_CM_only.RDS"))

##### DEG identification using LOO
## Parameters
min.pct_thd <- c(0.1, 0.2, 0.3, 0.4, 0.5)
logfc.threshold_thd <- c(0.1, 0.2, 0.3, 0.4, 0.5)
# min.pct_thd <- 0.05
# logfc.threshold_thd <- 0
p_val_adj_thd <- 0.05
param_set <- expand.grid(min.pct_thd = min.pct_thd, logfc.threshold_thd = logfc.threshold_thd)

##### Find DEGs for all combinations of leave one out sample combinations
pct.n <- min(param_set[, "min.pct_thd"])
threshold_thd.n <- min(param_set[, "logfc.threshold_thd"])
cat("Processing min.pct_thd:", pct.n, "\t logfc.threshold_thd:", threshold_thd.n, "...\n")

marker_genes_LOO <- list()
iter_index <- 0

for (i in 1:length(Control_samples)) {
# for (i in 1) {
    samples.i <- Control_samples[-i]
    # samples.i <- Control_samples
    for (j in 1:length(IHD_samples)) {
    # for (j in 1:length(c(Control_samples, IHD_samples))) {
    # for (j in 1) {
        samples.j <- IHD_samples[-j]
      # samples.j <- IHD_samples
        
        # LOO samples
        iter_index <- iter_index + 1
        samples_LOO <- c(samples.i, samples.j)
        # samples_LOO <- c(Control_samples, IHD_samples)[-j]
        # cat(iter_index, "Calculating LOO DEGs for ", samples_LOO, ")\n")
        cat(iter_index, "Calculating LOO DEGs for Control samples: (", samples.i, ") and IHD samples: (", samples.j, ")\n")
        
        # Subset seurat object by samples
        # control.disease.combined_CM_only_LOO <- control.disease.combined_CM_only[, control.disease.combined_CM_only$donor_condition %in% samples_LOO]
        Idents(object = control.disease.combined_CM_only) <- "donor_condition"
        control.disease.combined_CM_only_LOO <- subset(control.disease.combined_CM_only, idents = samples_LOO)
        
        # Marker genes for each cell type
        marker_genes <- c()
        for (celltype in celltypes) {
            # Identify marker genes for cell type by comparing IHD vs Control
            cat(celltype, ":\n")
            control.disease.combined_CM_only_LOO_celltype <- control.disease.combined_CM_only_LOO[, control.disease.combined_CM_only_LOO$annotated_clusters_CM %in% celltype]
            Idents(control.disease.combined_CM_only_LOO_celltype) <- control.disease.combined_CM_only_LOO_celltype$condition
            marker_gene_celltype <- FindMarkers(control.disease.combined_CM_only_LOO_celltype, ident.1 = "IHD", ident.2 = "Control",
                                                min.pct = pct.n, logfc.threshold = threshold_thd.n,
                                                test.use = "wilcox")
                                                # latent.vars = "num_cells_per_donor.condition")
            # marker_gene_celltype <- FindMarkers(control.disease.combined_CM_only_LOO_celltype, ident.1 = "Control", ident.2 = "IHD", 
            #                                     min.pct = pct.n, logfc.threshold = threshold_thd.n, 
            #                                     test.use = "wilcox")
            # latent.vars = "num_cells_per_donor.condition")
            marker_gene_celltype <- marker_gene_celltype[as.numeric(marker_gene_celltype$p_val_adj)<p_val_adj_thd, ]
            if (nrow(marker_gene_celltype) == 0) {
                cat("Number of marker genes:", nrow(marker_gene_celltype), "\n")
                next
            }
            marker_gene_celltype <- cbind(cbind(marker_gene_celltype, rownames(marker_gene_celltype)), celltype)
            colnames(marker_gene_celltype)[c(6,7)] <- c("gene", "celltype")
            cat("Number of marker genes:", nrow(marker_gene_celltype), "\n")

            # Update marker_genes
            marker_genes <- rbind(marker_genes, marker_gene_celltype)
        }
        marker_genes_LOO[[paste(samples_LOO, collapse = "-")]] <- marker_genes
    }
}
saveRDS(marker_genes_LOO, file = paste0(result_dir, "/LOO_DEG_Cardiomyocyte/marker_genes_CM_only_LOO.RDS"))

##### Identify the DEGs for different parameter threshold
celltype_of_interest <- "Cardiomyocyte"
marker_genes_LOO <- readRDS(paste0(result_dir, "/LOO_DEG_Cardiomyocyte/marker_genes_CM_only_LOO.RDS"))
marker_genes_LOO_table <- do.call(rbind, marker_genes_LOO)
# write.csv(marker_genes_LOO_table, paste0(result_dir, "/LOO_DEG_Cardiomyocyte/marker_genes_CM_only_LOO.csv"), row.names = TRUE)

## Process and identify DEGs in all combinations for the cell type of interest
LOO_found_thd <- 9
Num_of_DEGs <- c()
pct.n <- 0.10
threshold_thd.up <- 0.3
threshold_thd.down <- -0.25
# threshold_thd.up <- 0.5
# threshold_thd.down <- -0.5
# for (n in 1:nrow(param_set)) {
for (n in 1) {
    # pct.n <- param_set[n, "min.pct_thd"]
    # threshold_thd.n <- param_set[n, "logfc.threshold_thd"]
    # pct.n <- 0.1
    # threshold_thd.n <- 0.3
    cat("Processing min.pct_thd:", pct.n, "\t logfc.threshold_thd:", threshold_thd.up, threshold_thd.down, "...\n")

    DEGs_celltype <- c()
    for (i in 1:length(marker_genes_LOO)) {
        # Marker genes for cell type in LOO i
        DEGs_celltype.i <- marker_genes_LOO[[i]]
        DEGs_celltype.i <- DEGs_celltype.i[DEGs_celltype.i$celltype %in% celltype_of_interest, ]
        # Update marker_genes_celltype
        DEGs_celltype <- rbind(DEGs_celltype, DEGs_celltype.i)
    }
    DEGs_celltype$pct_diff <- abs(as.numeric(DEGs_celltype$pct.1) - as.numeric(DEGs_celltype$pct.2))
    
    # Filtering
    # DEGs_celltype_filtered <- DEGs_celltype[abs(as.numeric(DEGs_celltype$avg_log2FC))>threshold_thd.n, ]
    DEGs_celltype_filtered <- DEGs_celltype[as.numeric(DEGs_celltype$avg_log2FC) > threshold_thd.up | as.numeric(DEGs_celltype$avg_log2FC) < threshold_thd.down, ]
    DEGs_celltype_filtered <- DEGs_celltype_filtered[as.numeric(DEGs_celltype_filtered$pct.1)>pct.n | as.numeric(DEGs_celltype_filtered$pct.2)>pct.n, ]
    # DEGs_celltype_filtered <- DEGs_celltype_filtered[as.numeric(DEGs_celltype_filtered$pct.1)>pct.n & as.numeric(DEGs_celltype_filtered$pct.2)>pct.n, ]
    # ggplot(DEGs_celltype[DEGs_celltype$avg_log2FC >=0,], aes(x = celltype, y = avg_log2FC)) + geom_violin(trim = TRUE, na.rm = FALSE)
    
    # ggplot(DEGs_celltype, aes(x = celltype, y = pct.2)) + geom_violin(trim = TRUE, na.rm = FALSE)
    
    for (thd in LOO_found_thd) {
        # Up regulated DEGs
        DEGs_celltype_up <- DEGs_celltype_filtered[as.numeric(DEGs_celltype_filtered$avg_log2FC)>0, ]
        DEGs_up_freq <- table(DEGs_celltype_up$gene)
        DEGs_up_celltype_filtered <- names(DEGs_up_freq)[which(DEGs_up_freq >= thd)]
        if (length(DEGs_up_celltype_filtered)>1) {
            mean_log2FC <- tapply(DEGs_celltype_up$avg_log2FC, DEGs_celltype_up$gene, mean)
            DEGs_up_celltype_filtered_avgfc <- mean_log2FC[DEGs_up_celltype_filtered]
            DEGs_up_celltype_filtered <- DEGs_up_celltype_filtered[order(DEGs_up_celltype_filtered_avgfc, decreasing = F)]
        }
        df_up <- DEGs_celltype_up[DEGs_celltype_up$gene %in% DEGs_up_celltype_filtered, ]

        # Down regulated DEGs
        DEGs_celltype_dn <- DEGs_celltype_filtered[as.numeric(DEGs_celltype_filtered$avg_log2FC)<0, ]
        DEGs_dn_freq <- table(DEGs_celltype_dn$gene)
        DEGs_dn_celltype_filtered <- names(DEGs_dn_freq)[which(DEGs_dn_freq >= thd)]
        if (length(DEGs_dn_celltype_filtered)>1) {
            mean_log2FC <- tapply(DEGs_celltype_dn$avg_log2FC, DEGs_celltype_dn$gene, mean)
            DEGs_dn_celltype_filtered_avgfc <- mean_log2FC[DEGs_dn_celltype_filtered]
            DEGs_dn_celltype_filtered <- DEGs_dn_celltype_filtered[order(DEGs_dn_celltype_filtered_avgfc, decreasing = F)]
        }
        df_down <- DEGs_celltype_dn[DEGs_celltype_dn$gene %in% DEGs_dn_celltype_filtered, ]

        # cat("LOO_found_thd:", thd, "\n")
        # cat("Number of Up-regulated DEGs per cell type after filtering:", length(DEGs_up_celltype_filtered), "\n")
        # cat("Number of Down-regulated DEGs per cell type after filtering:", length(DEGs_dn_celltype_filtered), "\n")
        # cat("Total number of DEGs per cell type after filtering:", length(DEGs_up_celltype_filtered) + length(DEGs_dn_celltype_filtered), "\n")
        Num_of_DEGs.i <- c(pct.n, threshold_thd.up, threshold_thd.down, thd, length(DEGs_up_celltype_filtered), length(DEGs_dn_celltype_filtered), length(DEGs_up_celltype_filtered) + length(DEGs_dn_celltype_filtered))
        Num_of_DEGs <- rbind(Num_of_DEGs, Num_of_DEGs.i)
    }
}
rownames(Num_of_DEGs) <- NULL
colnames(Num_of_DEGs) <- c("pct_thd", "log2fc_thd_up", "log2fc_thd_down", "LOO_thd", "Num_of_DEGs_up", "Num_of_DEGs_dn", "Num_of_DEGs")
Num_of_DEGs

control.disease.combined_CM_only$donor_condition <- factor(x = control.disease.combined_CM_only$donor_condition, levels = c(Control_samples, IHD_samples))

DEGs_up_celltype_filtered <- setdiff(DEGs_up_celltype_filtered, GenesonChrXY)
valid_genes_up <- intersect(DEGs_up_celltype_filtered, rownames(control.disease.combined_CM_only[["RNA"]]@scale.data))
# DEGs_up_celltype_filtered_avgfc
DEGs_dn_celltype_filtered <- setdiff(DEGs_dn_celltype_filtered, GenesonChrXY)
valid_genes_dn <- intersect(DEGs_dn_celltype_filtered, rownames(control.disease.combined_CM_only[["RNA"]]@scale.data))
# DEGs_dn_celltype_filtered_avgfc
DEG_up_down <- c(DEGs_up_celltype_filtered, DEGs_dn_celltype_filtered)

##### using pheatmap
deg_expression_matrix <- as.matrix(GetAssayData(object = control.disease.combined_CM_only, assay = "RNA", slot = "count")[DEG_up_down, ])
sample_metadata <- control.disease.combined_CM_only@meta.data
scaled_data <- t(apply(deg_expression_matrix, 1, scale))
colnames(scaled_data) <- rownames(sample_metadata)

##### count of number of cells expressed certain gene in each donor_condition
# results <- data.frame(donor_condition = c("936_Ctrl", "5087_Ctrl", "5919_Ctrl", "1156_IHD", "5111_IHD", "5874_IHD"))
# for (gene in DEG_up_down) {
#   expr_data <- GetAssayData(control.disease.combined_CM_only, assay = "RNA", slot = "count")[gene, ]
#   expr_data <- as.matrix(expr_data)
#   cell_counts <- data.frame(cell = rownames(expr_data), expression = expr_data)
#   cell_counts$donor_condition <- sample_metadata$donor_condition
#   cell_counts <- cell_counts[cell_counts$expression != 0, ]
#   count_summary <- as.data.frame(table(cell_counts$donor_condition)) %>% select(-Var1) %>% setNames(gene)
#   results <- cbind(results, count_summary)
# }

##### draw violin plot for gene expression for all DEGs
top_up_DEGs_Sangita <- c("EGR1", "ZFP36", "AC023494.1", "CCBE1", "FOS", "BTG2", "DOK6", "NR4A1", "ZNF608", "PDGFD")
top_down_DEGs_Sangita <- c("ANKRD2", "SLC35F1", "SYNDIG1", "KCNAB2", "KIF26B", "TMEM120B", "SCN5A")
              
# expr_data_up <- as.matrix(GetAssayData(object = control.disease.combined_CM_only, assay = "RNA", slot = "count")[DEGs_up_celltype_filtered, ])
# expr_data_down <- as.matrix(GetAssayData(object = control.disease.combined_CM_only, assay = "RNA", slot = "count")[DEGs_dn_celltype_filtered, ])

expr_data_up <- as.matrix(GetAssayData(object = control.disease.combined_CM_only, assay = "RNA", slot = "count")[top_up_DEGs_Sangita, ])
expr_data_down <- as.matrix(GetAssayData(object = control.disease.combined_CM_only, assay = "RNA", slot = "count")[top_down_DEGs_Sangita, ])

expr_data_up <- as.data.frame(t(expr_data_up)) # 1:30, 31:60, 61:83
expr_data_up$donor_condition <- sample_metadata$donor_condition
expr_data_down <- as.data.frame(t(expr_data_down)) # 1:30, 31:60, 61:80
expr_data_down$donor_condition <- sample_metadata$donor_condition

expr_data_up_long <- melt(expr_data_up, id.vars = "donor_condition", variable.name = "gene", value.name = "expression")
expr_data_down_long <- melt(expr_data_down, id.vars = "donor_condition", variable.name = "gene", value.name = "expression")

# p_DEG_up_violin <- ggplot(expr_data_up_long, aes(x = donor_condition, y = expression, fill = donor_condition)) +
#   geom_violin(scale = "width", adjust = 1, trim = T) + geom_boxplot(width=0.2, color="white", alpha=0.2) + 
#   facet_wrap(~ gene, scales = "free", ncol = 5) + scale_y_log10() + 
#   labs(title = "", x = "Donor Condition", y = "Gene Expression (counts)") + theme_linedraw() + ggtitle(NULL) + RotatedAxis() +
#   theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
# ggsave(paste0(result_dir, "/LOO_DEG_Cardiomyocyte/p_DEG_up_violin.pdf"), plot = p_DEG_up_violin, width = 35, height = 25, dpi = 300)
# ggsave(paste0(result_dir, "/LOO_DEG_Cardiomyocyte/p_DEG_up_violin_S.pdf"), plot = p_DEG_up_violin, width = 15, height = 6, dpi = 300)

# p_DEG_down_violin <- ggplot(expr_data_down_long, aes(x = donor_condition, y = expression, fill = donor_condition)) +
#   geom_violin(scale = "width", adjust = 1, trim = T) +
#   facet_wrap(~ gene, scales = "free", ncol = 4) + scale_y_log10() + 
#   labs(title = "", x = "Donor Condition", y = "Gene Expression (counts)") + theme_linedraw() + ggtitle(NULL) + RotatedAxis() +
#   theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
# ggsave(paste0(result_dir, "/LOO_DEG_Cardiomyocyte/p_DEG_down_violin.pdf"), plot = p_DEG_down_violin, width = 35, height = 25, dpi = 300)
# ggsave(paste0(result_dir, "/LOO_DEG_Cardiomyocyte/p_DEG_down_violin_S.pdf"), plot = p_DEG_down_violin, width = 12, height = 6, dpi = 300)

# FeaturePlot(control.disease.combined_CM_only, features = "PIK3R1", split.by = "donor_condition",
#                   ncol = 4, max.cutoff = NA, cols = c("grey", "red"))
# FeaturePlot(control.disease.combined_CM_only, features = "GALNT17", split.by = "donor_condition",
#             ncol = 4, max.cutoff = NA, cols = c("grey", "red"))


##### Shuffle cells within each sample
sample_info <- sapply(strsplit(colnames(scaled_data), "_"), function(x) {paste0(x[1], "_", x[3])})
scaled_data <- do.call(cbind, lapply(unique(sample_info), function(sample) {
  sample_cols <- scaled_data[, sample_info == sample]
  sample_cols_shuffled <- sample_cols[, sample(1:ncol(sample_cols))]
  return(sample_cols_shuffled)}))

annotation_cols <- data.frame(condition = sample_metadata$condition, donor.condition = sample_metadata$donor_condition, row.names = rownames(sample_metadata))
annotation_rows <- data.frame(DEG.Group = factor(rep(c("Up.regulated", "Down.regulated"), c(length(DEGs_up_celltype_filtered), length(DEGs_dn_celltype_filtered)))),
                              row.names = DEG_up_down)

annotation_colors <- list(condition = c(Control = "dodgerblue2", IHD = "firebrick2"), 
                          donor.condition = c("936_Ctrl" = "#ADD8E6", "5087_Ctrl" = "#2D6EA9", "5919_Ctrl" = "#00008B", "1156_IHD" = "#FFCCCC", "5111_IHD" = "#DD555D", "5874_IHD" = "#660000"), 
                          DEG.Group = c("Up.regulated" = "#66A61E", "Down.regulated" = "#e377c2"))

# gaps <- cumsum(table(sample_info))
# heatmap_colors <- colorRampPalette(c("purple", "black", "yellow"))(100)
# breaks <- seq(-2.5, 2.5, length.out = 101)

# purple_colors <- colorRampPalette(c("purple1", "purple4"))(1)
# yellow_colors <- colorRampPalette(c("yellow2", "yellow"))(4)
# heatmap_colors <- colorRampPalette(c(purple_colors, "black", yellow_colors))(100)
# breaks <- seq(-1.5, 6, length.out = 101)

# purple_colors <- colorRampPalette(c("purple1", "purple4"))(1)
# yellow_colors <- colorRampPalette(c("yellow2", "yellow"))(2)
# heatmap_colors <- colorRampPalette(c(purple_colors, "black", yellow_colors))(100)
# breaks <- seq(-1.5, 3, length.out = 101)

# best
purple_colors <- colorRampPalette(c("purple2", "purple4"))(2)
yellow_colors <- colorRampPalette(c("yellow3", "yellow1"))(3)
heatmap_colors <- colorRampPalette(c(purple_colors, "black", yellow_colors))(100)
breaks <- seq(-2, 3, length.out = 101)

Num_of_DEGs
p <- pheatmap(mat = scaled_data, color = heatmap_colors, annotation_col = annotation_cols, annotation_row = annotation_rows, annotation_colors = annotation_colors, 
  cluster_rows = F, cluster_cols = F, scale = "none", fontsize_row = 8, fontsize_col = 8, border_color = "white", 
  show_colnames = F, show_rownames = F, legend = TRUE, annotation_legend = TRUE, gaps_col = NULL, gaps_row = NULL, breaks = breaks, useRaster=F)
ggsave(paste0(result_dir, "/LOO_DEG_Cardiomyocyte/DEG_pheatmap.png"), plot = p, width = 15, height = 6, dpi = 600)

DEG_up_down_df <- as.data.frame(rbind(cbind(DEGs_up_celltype_filtered, "up"), cbind(DEGs_dn_celltype_filtered, "down"))) %>% 
  setNames(c("gene", "regulation"))
write_csv(DEG_up_down_df, paste0(result_dir, "/LOO_DEG_Cardiomyocyte/DEG_up_down_df.csv"))
