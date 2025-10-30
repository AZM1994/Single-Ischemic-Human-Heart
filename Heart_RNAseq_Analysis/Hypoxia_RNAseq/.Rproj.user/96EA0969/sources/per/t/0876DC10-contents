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
result_dir <- paste0(wd_dir, "/results_IHD")
dir.create(paste0(result_dir, "/0-Seurat_Object"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(result_dir, "/1-QC_Control"), recursive = TRUE, showWarnings = FALSE)

##### read sample metadata
sample_metadata <- read.csv(paste0(wd_dir, "/data/sample_meta.csv"), header = TRUE)
sample_names <- sample_metadata$sample_name

##### read in all 10X files as matrices 
cellranger_data.dir <- sapply(sample_names, function(x) paste0(wd_dir, "/data/cellranger_results/",x, "/"))
list.cellranger.result <- list()
for (dir in cellranger_data.dir){
  print(basename(dir))
  list.cellranger.result[[basename(dir)]] <-  Read10X(data.dir = paste0(dir, "filtered_feature_bc_matrix/"))
}

##### create individual Seurat objects from matrices
options(Seurat.object.assay.version = "v3")
list.of.Seurat.objs <- lapply(seq_along(list.cellranger.result), function(x) 
  CreateSeuratObject(counts = list.cellranger.result[[x]], min.cells = 3, min.features = 200, project = names(list.cellranger.result)[x])
)

names(list.of.Seurat.objs) <- names(list.cellranger.result)
for (obj_index in seq_along(list.of.Seurat.objs)){
  obj <- list.of.Seurat.objs[[obj_index]]
  name <- names(list.of.Seurat.objs)[obj_index]
  obj <- RenameCells(obj, add.cell.id = name)
  list.of.Seurat.objs[[obj_index]] <- obj                                    
}

##### Create merged Seurat object
control.disease.combined <- merge(list.of.Seurat.objs[[1]], y = list.of.Seurat.objs[2:17], project = "heart_merged_analysis")
# add more info to the merged object
control.disease.combined[["percent_mito"]] <- PercentageFeatureSet(control.disease.combined, pattern = "^MT-")
control.disease.combined[["percent_ribo"]] <- PercentageFeatureSet(control.disease.combined, pattern = "^RPL|^RPS")
control.disease.combined$samples <- sapply(control.disease.combined$orig.ident, function(x) sample_names[which(sample_names == x)])
control.disease.combined$condition <- sapply(control.disease.combined$orig.ident, function(x) sample_metadata$conditions[which(sample_names == x)])
control.disease.combined$donor <- sapply(control.disease.combined$orig.ident, function(x) sample_metadata$donor[which(sample_names == x)])
control.disease.combined$technique <- sapply(control.disease.combined$orig.ident, function(x) sample_metadata$technique[which(sample_names == x)])
control.disease.combined$donor_condition <- sapply(control.disease.combined$orig.ident, function(x) sample_metadata$donor_condition[which(sample_names == x)])
# head(control.disease.combined)

donor_condition_ordered <- unique(control.disease.combined$donor_condition)
control.disease.combined$samples <- factor(control.disease.combined$samples, levels = sample_names)
control.disease.combined$donor_condition <- factor(control.disease.combined$donor_condition, levels = donor_condition_ordered)

metadata <- control.disease.combined@meta.data
metadata <- metadata %>% rename(nUMI = nCount_RNA, nGene = nFeature_RNA)

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
  p_preQC11 <- QC_Plots_Genes(seurat_object = control.disease.combined, group.by = "samples", plot_title = "Genes Per Cell") + scale_fill_manual(values = color_list)
  p_preQC21 <- QC_Plots_UMIs(seurat_object = control.disease.combined, group.by = "samples", plot_title = "UMIs Per Cell") + scale_fill_manual(values = color_list)
  p_preQC31 <- QC_Plots_Mito(seurat_object = control.disease.combined, group.by = "samples", plot_title = "Mito Gene % Per Cell") + scale_fill_manual(values = color_list)
  p_preQC41 <- QC_Plots_Feature(seurat_object = control.disease.combined, group.by = "samples", feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell") + scale_fill_manual(values = color_list)
  wrap_plots(remove_point_outline(p_preQC11), remove_point_outline(p_preQC21), remove_point_outline(p_preQC31), remove_point_outline(p_preQC41), ncol = 4)
  
  p_preQC12 <- QC_Plots_Genes(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "Genes Per Cell") + scale_fill_manual(values = color_list)
  p_preQC22 <- QC_Plots_UMIs(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "UMIs Per Cell") + scale_fill_manual(values = color_list)
  p_preQC32 <- QC_Plots_Mito(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "Mito Gene % Per Cell") + scale_fill_manual(values = color_list)
  p_preQC42 <- QC_Plots_Feature(seurat_object = control.disease.combined, group.by = "donor_condition", feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell") + scale_fill_manual(values = color_list)
  wrap_plots(remove_point_outline(p_preQC12), remove_point_outline(p_preQC22), remove_point_outline(p_preQC32), remove_point_outline(p_preQC42), ncol = 4)
dev.off()

################ Cell counts: The cell counts are determined by the number of unique cellular barcodes detected.
pdf(paste0(result_dir, "/1-QC_Control/nCells_per_sample_donor_condition.pdf"), width = 20, height = 6)
  p_nCell_1 <- metadata %>% ggplot(aes(x = samples, fill = samples)) + geom_bar() + scale_fill_manual(values = color_list) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Sample ID") + ggtitle("nCells per sample")
  p_nCell_2 <- metadata %>% ggplot(aes(x = donor_condition, fill = donor_condition)) + geom_bar() + scale_fill_manual(values = color_list) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Donor") + ggtitle("nCells per donor")
  p_nCell_3 <- metadata %>% ggplot(aes(x = condition, fill = condition)) + geom_bar() + scale_fill_manual(values = Ctrl_IHD_color) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Condition") + ggtitle("nCells per condition")
  wrap_plots(p_nCell_1, p_nCell_2, p_nCell_3, ncol = 3)
dev.off()

##### UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500. If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply.
# Visualize the number UMIs/transcripts per cell
pdf(paste0(result_dir, "/1-QC_Control/nUMI_per_sample_donor_condition.pdf"), width = 20, height = 6)
  p_nUMI_1 <- metadata %>% ggplot(aes(x = nUMI, color = samples, fill = samples)) + geom_density(alpha = 0.2) + scale_x_log10() + 
    scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + 
    ylab("UMI Count") + geom_vline(xintercept = 500) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme(legend.title = element_blank()) + ggtitle("nUMI per sample")
  p_nUMI_2 <- metadata %>% ggplot(aes(x = nUMI, color = donor_condition, fill = donor_condition)) + geom_density(alpha = 0.2) + scale_x_log10() + 
    scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + 
    ylab("UMI Count") + geom_vline(xintercept = 500) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme(legend.title = element_blank()) + ggtitle("nUMI per donor")
  p_nUMI_3 <- metadata %>% ggplot(aes(x = nUMI, color = condition, fill= condition)) + geom_density(alpha = 0.2) + scale_x_log10() + 
    scale_fill_manual(values = Ctrl_IHD_color) + scale_color_manual(values = Ctrl_IHD_color) + theme_classic() + 
    ylab("UMI Count") + geom_vline(xintercept = 500) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme(legend.title = element_blank()) + ggtitle("nUMI per condition")
  wrap_plots(p_nUMI_1, p_nUMI_2, p_nUMI_3, ncol = 3)
dev.off()

##### Genes detected per cell
pdf(paste0(result_dir, "/1-QC_Control/nGene_per_sample_donor_condition.pdf"), width = 20, height = 6)
  ## Visualize the distribution of genes detected per cell via histogram
  p_nGene_1 <- metadata %>% ggplot(aes(x = nGene, color = samples, fill = samples)) + geom_density(alpha = 0.2) + scale_x_log10() + 
    scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + 
    geom_vline(xintercept = 300) + ggtitle("nGenes per cell for each sample") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(legend.title = element_blank())
  p_nGene_2 <- metadata %>% ggplot(aes(x = nGene, color = donor_condition, fill = donor_condition)) + geom_density(alpha = 0.2) + scale_x_log10() + 
    scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + 
    geom_vline(xintercept = 300) + ggtitle("nGenes per cell for each donor") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(legend.title = element_blank())
  # Visualize the distribution of genes detected per cell via boxplot
  p_nGene_3 <- metadata %>% ggplot(aes(x = condition, y=log10(nGene), fill = condition)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = Ctrl_IHD_color) + 
    theme(legend.title=element_blank()) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.title = element_text(hjust=0.5, face="bold")) + ggtitle("nGenes per cell for each condition")
  wrap_plots(p_nGene_1, p_nGene_2, p_nGene_3, ncol = 3)
dev.off()

##### Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata$log10GenesPerUMI <- control.disease.combined$log10GenesPerUMI
pdf(paste0(result_dir, "/1-QC_Control/nGenes_per_UMI_sample_donor_condition.pdf"), width = 20, height = 6)
  p_Genes_per_UMI_1 <- metadata %>% ggplot(aes(x = log10GenesPerUMI, color = samples, fill = samples)) + geom_density(alpha = 0.2) +
    scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + theme(legend.title=element_blank()) + geom_vline(xintercept = 0.8)
  p_Genes_per_UMI_2 <- metadata %>% ggplot(aes(x = log10GenesPerUMI, color = donor_condition, fill = donor_condition)) + geom_density(alpha = 0.2) +
    scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + theme(legend.title=element_blank()) + geom_vline(xintercept = 0.8)
  p_Genes_per_UMI_3 <- metadata %>% ggplot(aes(x = log10GenesPerUMI, color = condition, fill = condition)) + geom_density(alpha = 0.2) +
    scale_fill_manual(values = Ctrl_IHD_color) + scale_color_manual(values = Ctrl_IHD_color) + theme_classic() + theme(legend.title=element_blank()) + geom_vline(xintercept = 0.8)
  wrap_plots(p_Genes_per_UMI_1, p_Genes_per_UMI_2, p_Genes_per_UMI_3, ncol = 3)
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
control.disease.combined<- subset(control.disease.combined, subset = percent_mito <= 10)
control.disease.combined <- subset(control.disease.combined, subset = percent_ribo <= 5)
control.disease.combined <- subset(control.disease.combined, subset = nFeature_RNA >= nFeature_RNA_five_q & nFeature_RNA <= nFeature_RNA_ninty_five_q)
control.disease.combined <- subset(control.disease.combined, subset = nCount_RNA >= nCount_RNA_five_q & nCount_RNA <= nCount_RNA_ninty_five_q)
control.disease.combined <- subset(control.disease.combined, subset = log10GenesPerUMI >= 0.80)
##### remove 840_Ctrl, 1039_Ctrl, 5828_Ctrl high mito genes
control.disease.combined <- subset(control.disease.combined, subset = samples %in% setdiff(sample_names , c("840_3_Ctrl", "1039_3_Ctrl", "5828_3_Ctrl")))

##### Create plots of QC metrics after filtering
pdf(paste0(result_dir, "/1-QC_Control/QC_post_filtering.pdf"), width = 20, height = 3.5)
  ## Unique Molecular Identifiers, absolute number of observed transcripts
  p_postQC11 <- QC_Plots_Genes(seurat_object = control.disease.combined, group.by = "samples", plot_title = "Genes Per Cell") + scale_fill_manual(values = color_list)
  p_postQC21 <- QC_Plots_UMIs(seurat_object = control.disease.combined, group.by = "samples", plot_title = "UMIs Per Cell") + scale_fill_manual(values = color_list)
  p_postQC31 <- QC_Plots_Mito(seurat_object = control.disease.combined, group.by = "samples", plot_title = "Mito Gene % Per Cell") + scale_fill_manual(values = color_list)
  p_postQC41 <- QC_Plots_Feature(seurat_object = control.disease.combined, group.by = "samples", feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell") + scale_fill_manual(values = color_list)
  wrap_plots(remove_point_outline(p_postQC11), remove_point_outline(p_postQC21), remove_point_outline(p_postQC31), remove_point_outline(p_postQC41), ncol = 4)
  
  p_postQC12 <- QC_Plots_Genes(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "Genes Per Cell") + scale_fill_manual(values = color_list)
  p_postQC22 <- QC_Plots_UMIs(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "UMIs Per Cell") + scale_fill_manual(values = color_list)
  p_postQC32 <- QC_Plots_Mito(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "Mito Gene % Per Cell") + scale_fill_manual(values = color_list)
  p_postQC42 <- QC_Plots_Feature(seurat_object = control.disease.combined, group.by = "donor_condition", feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell") + scale_fill_manual(values = color_list)
  wrap_plots(remove_point_outline(p_postQC12), remove_point_outline(p_postQC22), remove_point_outline(p_postQC32), remove_point_outline(p_postQC42), ncol = 4)
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
  # RunHarmony(group.by.vars = "samples", lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  RunHarmony(group.by.vars = c("samples", "technique"), lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  RunTSNE(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = TRUE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = TRUE) %>%
  FindClusters(resolution = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6), verbose = TRUE)

pdf(paste0(result_dir, "/1-QC_Control/Dim_Plot_before_doublet.pdf"), width = 7, height = 5)
  ElbowPlot(control.disease.combined, ndims = 30, reduction = "pca")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'samples') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'condition') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'technique') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'samples') + ggtitle("Harmony")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'condition') + ggtitle("Harmony")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'technique') + ggtitle("Harmony")
  
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE) + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE) + ggtitle("Harmony")
dev.off()

saveRDS(control.disease.combined, paste0(result_dir, "/0-Seurat_Object/combined.Seurat.obj_before_doublet.RDS"))
# control.disease.combined <- readRDS(paste0(result_dir, "/0-Seurat_Object/combined.Seurat.obj_before_doublet.RDS"))

######################################################################################################
########################################## Doublet finder ############################################
######################################################################################################
sweep.res.list_combined <- paramSweep_v3(control.disease.combined, PCs = 1:30, sct = FALSE)
sweep.stats_combined <- summarizeSweep(sweep.res.list_combined, GT = FALSE)
bcmvn_combined <- find.pK(sweep.stats_combined)

##### Plot the bcmvn output to find the peak
pdf(paste0(result_dir, "/1-QC_Control/Optimal_pK_doubletFinder.pdf"), width = 12, height = 5)
  ggplot(bcmvn_combined, aes(x = pK, y = BCmetric)) + geom_bar(stat = "identity") + theme_minimal() + labs(title = "Optimal pK Value", x = "pK", y = "BCmvn")
dev.off()

homotypic.prop <- modelHomotypic(control.disease.combined@meta.data$seurat_clusters) 
nExp_poi <- round(0.075*nrow(control.disease.combined@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

control.disease.combined <- doubletFinder_v3(control.disease.combined, PCs = 1:30, pK = 0.040, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
head(control.disease.combined)
saveRDS(control.disease.combined, paste0(result_dir, "/0-Seurat_Object/combined.Seurat.obj_doubletFinder_v3.RDS"))
# control.disease.combined  <- readRDS(paste0(result_dir, "/0-Seurat_Object/combined.Seurat.obj_doubletFinder_v3.RDS"))

pdf(paste0(result_dir, "/1-QC_Control/DF_classifications.pdf"), width = 7, height = 5)
  DimPlot(control.disease.combined, reduction = "umap", group.by = "DF.classifications_0.25_0.06_5066")
dev.off()

# table(control.disease.combined@meta.data$DF.classifications_0.25_0.06_5066)

control.disease.combined <- subset(control.disease.combined, subset = DF.classifications_0.25_0.06_5066 == "Singlet")
saveRDS(control.disease.combined, paste0(result_dir, "/0-Seurat_Object/combined.Seurat.obj_remove_doublet.RDS"))
# control.disease.combined  <- readRDS(paste0(result_dir, "/0-Seurat_Object/combined.Seurat.obj_remove_doublet.RDS"))

control.disease.combined <- NormalizeData(control.disease.combined)
control.disease.combined <- FindVariableFeatures(control.disease.combined, selection.method = "vst", nfeatures = 3000)
control.disease.combined <- ScaleData(control.disease.combined, vars.to.regress = c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA"))

# saveRDS(control.disease.combined, paste0(result_dir, "/0-Seurat_Object/merged.Seurat.obj_filtered_Scaled_02.RDS"))
# control.disease.combined <- readRDS(paste0(result_dir, "/0-Seurat_Object/merged.Seurat.obj_filtered_Scaled_02.RDS"))
# check_PC_02 <- RunPCA(control.disease.combined, verbose = TRUE)
# ElbowPlot(check_PC_02, ndims = 50, reduction = "pca")

control.disease.combined <- RunPCA(control.disease.combined, npcs = 30, verbose = TRUE) %>% 
  # RunHarmony(group.by.vars = "samples", lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  RunHarmony(group.by.vars = c("samples", "technique"), lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  # RunTSNE(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  FindClusters(resolution = c(0.6, 0.8, 1.0, 1.2, 1.4), verbose = TRUE)

pdf(paste0(result_dir, "/1-QC_Control/Dim_Plot_after_doublet.pdf"), width = 7, height = 5)
  ElbowPlot(control.disease.combined, ndims = 30, reduction = "pca")
  # DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'samples') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'condition') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'technique') + ggtitle("UMAP")
  # DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'samples') + ggtitle("Harmony")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'condition') + ggtitle("Harmony")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'technique') + ggtitle("Harmony")
  
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE) + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE) + ggtitle("Harmony")
dev.off()

saveRDS(control.disease.combined, paste0(result_dir, "/0-Seurat_Object/combined.Seurat.obj_ready_for_annotation.RDS"))
