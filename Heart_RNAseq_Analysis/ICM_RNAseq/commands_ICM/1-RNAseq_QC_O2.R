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
color_list <- c("#EA7F82", "#EF6BA5", "#FE9E4F", "#FFC9DA", "#F4445E", "#FFD27B", "#FFA8A8", "#B3DCB2", 
                "#7CC657", "#97B86D", "#B6D9ED", "#87C0F4", "#B1B1E1", "#9F45A1", "#963EE3", "#EB88C9", "#6990E3")

##### set working and results directories
wd_dir <- getwd()
result_dir <- paste0(wd_dir, "/results_ICM")
dir.create(paste0(result_dir, "/0-Seurat_Object"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(result_dir, "/1-QC_Control"), recursive = TRUE, showWarnings = FALSE)

##### read counts matrices
ICM_counts <- ReadMtx(
  mtx = paste0(wd_dir, "/data/ICM_Expression_Matrix_raw_counts_V1.mtx"),
  features = paste0(wd_dir, "/data/ICM_Expression_Matrix_genes_V1.tsv"),
  cells = paste0(wd_dir, "/data/ICM_Expression_Matrix_barcodes_V1.tsv")
)

##### Create merged Seurat object
control.disease.combined <- CreateSeuratObject(counts = ICM_counts, project = "ICM")
##### read sample metadata
sample_names_ordered <- c("1452", "1650", "1690", "1716", "1739", "1763", "1785", "1801", "1364", "1579", "1693", "1703", "1733", "1773", "1800")
donor_condition_ordered <- c("1452_Control", "1650_Control", "1690_Control", "1716_Control", "1739_Control", "1763_Control", "1785_Control", "1801_Control", 
                             "1364_ICM", "1579_ICM", "1693_ICM", "1703_ICM", "1733_ICM", "1773_ICM", "1800_ICM")
sample_metadata <- read.delim(paste0(wd_dir, "/data/ICM_MetaData_V1.txt"), header = TRUE, sep = "\t") %>% slice(-1) %>% 
  mutate(condition = recode(disease__ontology_label, "cardiomyopathy" = "ICM", "normal" = "Control")) %>% 
  mutate(nUMI = as.numeric(n_umi), nGene = as.numeric(n_genes), samples = biosample_id, donor = biosample_id, donor_condition = paste0(donor, "_", condition), 
         samples = factor(samples, levels = sample_names_ordered), donor_condition = factor(donor_condition, levels = donor_condition_ordered))
control.disease.combined <- AddMetaData(control.disease.combined, sample_metadata)

# add more info to the merged object
control.disease.combined[["percent_mito"]] <- PercentageFeatureSet(control.disease.combined, pattern = "^MT-")
control.disease.combined[["percent_ribo"]] <- PercentageFeatureSet(control.disease.combined, pattern = "^RPL|^RPS")

######################################################################################################
################### Quality control and filtering for control.disease.combined #######################
######################################################################################################
control.disease.combined<- Add_Cell_Complexity(object = control.disease.combined, overwrite = TRUE)
all.equal(colnames(control.disease.combined), row.names(control.disease.combined@meta.data))

remove_point_outline <- function(p) {
  for (i in seq_along(p$layers)) {
    if ("GeomPoint" %in% class(p$layers[[i]]$geom)) {
      p$layers[[i]]$aes_params$stroke <- 0
    }
  }
  p
}

##### create plots of QC metrics before filtering
pdf(paste0(result_dir, "/1-QC_Control/QC_pre_filtering.pdf"), width = 20, height = 3.5)
  ## Unique Molecular Identifiers, absolute number of observed transcripts
  p_preQC1 <- QC_Plots_Genes(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "Genes Per Cell") + scale_fill_manual(values = color_list)
  p_preQC2 <- QC_Plots_UMIs(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "UMIs Per Cell") + scale_fill_manual(values = color_list)
  p_preQC3 <- QC_Plots_Mito(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "Mito Gene % Per Cell") + scale_fill_manual(values = color_list)
  p_preQC4 <- QC_Plots_Feature(seurat_object = control.disease.combined, group.by = "donor_condition", feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell") + scale_fill_manual(values = color_list)
  wrap_plots(remove_point_outline(p_preQC1), remove_point_outline(p_preQC2), remove_point_outline(p_preQC3), remove_point_outline(p_preQC4), ncol = 4)
dev.off()

################ Cell counts: The cell counts are determined by the number of unique cellular barcodes detected.
pdf(paste0(result_dir, "/1-QC_Control/nCells_per_donor_condition.pdf"), width = 12, height = 6)
  p_nCell_1 <- sample_metadata %>% ggplot(aes(x = donor_condition, fill = donor_condition)) + geom_bar() + scale_fill_manual(values = color_list) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Sample ID") + ggtitle("nCells per sample")
  p_nCell_2 <- sample_metadata %>% ggplot(aes(x = condition, fill = condition)) + geom_bar() + scale_fill_manual(values = Ctrl_IHD_color) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Condition") + ggtitle("nCells per condition")
  wrap_plots(p_nCell_1, p_nCell_2, ncol = 2)
dev.off()

##### UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500. If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply.
# Visualize the number UMIs/transcripts per cell
pdf(paste0(result_dir, "/1-QC_Control/nUMI_per_donor_condition.pdf"), width = 12, height = 6)
  p_nUMI_1 <- sample_metadata %>% ggplot(aes(x = nUMI, color = donor_condition, fill = donor_condition)) + geom_density(alpha = 0.2) + scale_x_log10() + 
    scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + 
    ylab("UMI Count") + geom_vline(xintercept = 500) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme(legend.title = element_blank()) + ggtitle("nUMI per sample")
  p_nUMI_2 <- sample_metadata %>% ggplot(aes(x = nUMI, color = condition, fill= condition)) + geom_density(alpha = 0.2) + scale_x_log10() + 
    scale_fill_manual(values = Ctrl_IHD_color) + scale_color_manual(values = Ctrl_IHD_color) + theme_classic() + 
    ylab("UMI Count") + geom_vline(xintercept = 500) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme(legend.title = element_blank()) + ggtitle("nUMI per condition")
  wrap_plots(p_nUMI_1, p_nUMI_2, ncol = 2)
dev.off()

##### Genes detected per cell
pdf(paste0(result_dir, "/1-QC_Control/nGene_per_donor_condition.pdf"), width = 12, height = 6)
  ## Visualize the distribution of genes detected per cell via histogram
  p_nGene_1 <- sample_metadata %>% ggplot(aes(x=nGene, color = donor_condition, fill = donor_condition)) + geom_density(alpha = 0.2) + scale_x_log10() + 
    scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + 
    geom_vline(xintercept = 300) + ggtitle("nGenes per cell for each sample") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(legend.title = element_blank())
  # Visualize the distribution of genes detected per cell via boxplot
  p_nGene_2 <- sample_metadata %>% ggplot(aes(x=condition, y=log10(nGene), fill=condition)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = Ctrl_IHD_color) + 
    theme(legend.title=element_blank()) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.title = element_text(hjust=0.5, face="bold")) + ggtitle("nGenes per cell for each condition")
  wrap_plots(p_nGene_1, p_nGene_2, ncol = 2)
dev.off()

##### Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
sample_metadata$log10GenesPerUMI <- control.disease.combined$log10GenesPerUMI
pdf(paste0(result_dir, "/1-QC_Control/nGenes_per_UMI_sample_donor_condition.pdf"), width = 12, height = 6)
  p_Genes_per_UMI_1 <- sample_metadata %>% ggplot(aes(x = log10GenesPerUMI, color = donor_condition, fill = donor_condition)) + geom_density(alpha = 0.2) +
    scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + theme(legend.title=element_blank()) + geom_vline(xintercept = 0.8)
  p_Genes_per_UMI_2 <- sample_metadata %>% ggplot(aes(x = log10GenesPerUMI, color = condition, fill = condition)) + geom_density(alpha = 0.2) +
    scale_fill_manual(values = Ctrl_IHD_color) + scale_color_manual(values = Ctrl_IHD_color) + theme_classic() + theme(legend.title=element_blank()) + geom_vline(xintercept = 0.8)
  wrap_plots(p_Genes_per_UMI_1, p_Genes_per_UMI_2, ncol = 2)
dev.off()

##### Determine cells that will be lost by each QC filtering step
nFeature_RNA_five_q <- quantile(control.disease.combined$nFeature_RNA, 0.05)
nFeature_RNA_ninty_five_q <- quantile(control.disease.combined$nFeature_RNA, 0.95)
nCount_RNA_five_q <- quantile(control.disease.combined$nCount_RNA, 0.05)
nCount_RNA_ninty_five_q <- quantile(control.disease.combined$nCount_RNA, 0.95)

cells_filtered_pct_mito <- row.names(control.disease.combined@meta.data %>% filter(percent_mito > 10))
cells_filtered_pct_ribo <- row.names(control.disease.combined@meta.data %>% filter(percent_ribo > 5))
cells_filtered_nFeatures <- row.names(control.disease.combined@meta.data %>% filter(nFeature_RNA < nFeature_RNA_five_q | nFeature_RNA > nFeature_RNA_ninty_five_q))
cells_filtered_nCounts <- row.names(control.disease.combined@meta.data %>% filter(nCount_RNA < nCount_RNA_five_q | nCount_RNA > nCount_RNA_ninty_five_q))
cells_filtered_log10GenesPerUMI <- row.names(control.disease.combined@meta.data %>% filter(log10GenesPerUMI < 0.80))
filtering_list <- list("pct_mito" = cells_filtered_pct_mito, "pct_ribo" = cells_filtered_pct_ribo, "nFeatures" = cells_filtered_nFeatures, 
                       "nCounts" = cells_filtered_nCounts, "log10GenesPerUMI" = cells_filtered_log10GenesPerUMI)
print("Number of cells filtered per criteria: ")
print(lapply(filtering_list, length))
print(length(unique(do.call(c, filtering_list))))
print(length(unique(do.call(c, filtering_list)))/ncol(control.disease.combined))

pdf(paste0(result_dir, "/1-QC_Control/cells_filtered_upset_plot.pdf"))
  upset(fromList(filtering_list))
dev.off()

##### filter cells 
control.disease.combined <- subset(control.disease.combined, subset = percent_mito <= 10)
control.disease.combined <- subset(control.disease.combined, subset = percent_ribo <= 5)
control.disease.combined <- subset(control.disease.combined, subset = nFeature_RNA >= nFeature_RNA_five_q & nFeature_RNA <= nFeature_RNA_ninty_five_q)
control.disease.combined <- subset(control.disease.combined, subset = nCount_RNA >= nCount_RNA_five_q & nCount_RNA <= nCount_RNA_ninty_five_q)
control.disease.combined <- subset(control.disease.combined, subset = log10GenesPerUMI >= 0.80)

##### Create plots of QC metrics after filtering
pdf(paste0(result_dir, "/1-QC_Control/QC_post_filtering.pdf"), width = 20, height = 3.5)
  ## Unique Molecular Identifiers, absolute number of observed transcripts
  p_postQC1 <- QC_Plots_Genes(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "Genes Per Cell") + scale_fill_manual(values = color_list)
  p_postQC2 <- QC_Plots_UMIs(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "UMIs Per Cell") + scale_fill_manual(values = color_list)
  p_postQC3 <- QC_Plots_Mito(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "Mito Gene % Per Cell") + scale_fill_manual(values = color_list)
  p_postQC4 <- QC_Plots_Feature(seurat_object = control.disease.combined, group.by = "donor_condition", feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell") + scale_fill_manual(values = color_list)
  wrap_plots(remove_point_outline(p_postQC1), remove_point_outline(p_postQC2), remove_point_outline(p_postQC3), remove_point_outline(p_postQC4), ncol = 4)
dev.off()

saveRDS(control.disease.combined, paste0(result_dir, "/0-Seurat_Object/merged.Seurat.obj_filtered.RDS"))

######################################################################################################
################## Normalization, scaling, dimensionality reduction, clustering ######################
######################################################################################################
# control.disease.combined <- readRDS(paste0(result_dir, "/0-Seurat_Object/merged.Seurat.obj_filtered.RDS"))
n_cells <- ncol(control.disease.combined)

DefaultAssay(control.disease.combined) <- 'RNA'
control.disease.combined <- NormalizeData(control.disease.combined)
control.disease.combined <- FindVariableFeatures(control.disease.combined, selection.method = "vst", nfeatures = 3000)
control.disease.combined <- ScaleData(control.disease.combined, vars.to.regress = c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA"))

saveRDS(control.disease.combined, paste0(result_dir, "/0-Seurat_Object/merged.Seurat.obj_filtered_Scaled_01.RDS"))
# control.disease.combined <- readRDS(paste0(result_dir, "/0-Seurat_Object/merged.Seurat.obj_filtered_Scaled_01.RDS"))
# check_PC <- RunPCA(control.disease.combined, verbose = TRUE)
# ElbowPlot(check_PC, ndims = 50, reduction = "pca")

control.disease.combined <- RunPCA(control.disease.combined, npcs = 30, verbose = TRUE) %>% 
  RunHarmony(group.by.vars = "samples", lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  # RunHarmony(group.by.vars = c("samples", "technique"), lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  FindClusters(resolution = c(0.5, 0.6, 0.8, 1.0, 1.2), verbose = TRUE)

pdf(paste0(result_dir, "/1-QC_Control/Dim_Plot_before_doublet.pdf"), width = 7, height = 5)
  ElbowPlot(control.disease.combined, ndims = 30, reduction = "pca")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'donor_condition') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'condition') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'donor_condition') + ggtitle("Harmony")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'condition') + ggtitle("Harmony")
  
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE) + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE) + ggtitle("Harmony")
dev.off()

saveRDS(control.disease.combined, paste0(result_dir, "/0-Seurat_Object/combined.Seurat.obj_before_doublet.RDS"))
# control.disease.combined <- readRDS(paste0(result_dir, "/0-Seurat_Object/combined.Seurat.obj_before_doublet.RDS"))

######################################################################################################
########################################## Doublet finder ############################################
######################################################################################################
control.disease.combined$doublet_prediction <- ifelse(control.disease.combined$doublet_score > 0.3, "Doublet", "Singlet")
# table(control.disease.combined@meta.data$doublet_prediction)
control.disease.combined <- subset(control.disease.combined, subset = doublet_prediction == "Singlet")

control.disease.combined <- NormalizeData(control.disease.combined)
control.disease.combined <- FindVariableFeatures(control.disease.combined, selection.method = "vst", nfeatures = 3000)
control.disease.combined <- ScaleData(control.disease.combined, vars.to.regress = c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA"))

saveRDS(control.disease.combined, paste0(result_dir, "/0-Seurat_Object/merged.Seurat.obj_filtered_Scaled_02.RDS"))
# control.disease.combined <- readRDS(paste0(result_dir, "/0-Seurat_Object/merged.Seurat.obj_filtered_Scaled_02.RDS"))
# check_PC_02 <- RunPCA(control.disease.combined, verbose = TRUE)
# ElbowPlot(check_PC_02, ndims = 50, reduction = "pca")

control.disease.combined <- RunPCA(control.disease.combined, npcs = 30, verbose = TRUE) %>% 
  RunHarmony(group.by.vars = "samples", lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  FindClusters(resolution = c(0.5, 0.6, 0.8, 1.0, 1.2), verbose = TRUE)

control.disease.combined <- RunPCA(control.disease.combined, npcs = 30, verbose = TRUE) %>% 
  RunHarmony(group.by.vars = "samples", lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony")

control.disease.combined <- RunUMAP(control.disease.combined, reduction = "harmony", dims = 1:30, 
                                    n.neighbors = 30, min.dist = 0.5, verbose = TRUE)

control.disease.combined <- FindNeighbors(control.disease.combined, reduction = "harmony", dims = 1:30, verbose = TRUE)
control.disease.combined <- FindClusters(control.disease.combined, resolution = c(0.5, 0.6, 0.8, 1.0, 1.2), verbose = TRUE)

pdf(paste0(result_dir, "/1-QC_Control/Dim_Plot_after_doublet.pdf"), width = 7, height = 5)
  ElbowPlot(control.disease.combined, ndims = 30, reduction = "pca")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'donor_condition') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'donor_condition') + ggtitle("Harmony")
  
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE) + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE) + ggtitle("Harmony")
dev.off()

# saveRDS(control.disease.combined, paste0(result_dir, "/0-Seurat_Object/combined.Seurat.obj_ready_for_annotation.RDS"))
