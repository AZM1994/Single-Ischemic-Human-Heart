library(Seurat)
library(patchwork)
library(dplyr)
library(data.table)
library(ggplot2)
# library(ggpubr)
library(scCustomize)
library(UpSetR)
library(harmony)
library(sctransform)
library(DoubletFinder)
library(glmGamPoi)
library(pheatmap)

Ctrl_IHD_color <- c("#2D6EA8", "#DD555B")
color_list <- c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                "#393b79", "#5254a3", "#6b6ecf", "#9c9ede", "#637939", "#8ca252", "#b5cf6b", "#cedb9c", "#8c6d31", "#bd9e39",
                "#e7ba52", "#e7969c", "#d6616b", "#ad494a", "#843c39", "#f5f5f5", "#fdd0a2", "#fb6a4a", "#cb181d", "#a50f15",
                "#3182bd", "#6baed6", "#9ecae1", "#c6dbef", "#e6550d", "#fd8d3c", "#fdae6b", "#fdd0a2", "#31a354", "#74c476",
                "#a1d99b", "#c7e9c0", "#756bb1", "#9e9ac8", "#bcbddc")

##### set working and results directories
wd_dir <- getwd()
# result_dir <- paste0(wd_dir, "/results")
result_dir <- paste0(wd_dir, "/results_5000_PC30")
dir.create(paste0(result_dir, "/QC_Control"), recursive = TRUE)
dir.create(paste0(result_dir, "/Annotation"), recursive = TRUE)

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
head(control.disease.combined)

donor_condition_ordered <- unique(control.disease.combined$donor_condition)
control.disease.combined$samples <- factor(control.disease.combined$samples, levels = sample_names)
control.disease.combined$donor_condition <- factor(control.disease.combined$donor_condition, levels = donor_condition_ordered)

metadata <- control.disease.combined@meta.data
metadata <- metadata %>% rename(nUMI = nCount_RNA, nGene = nFeature_RNA)

######################################################################################################
######################################################################################################
################### Quality control and filtering for control.disease.combined #######################
######################################################################################################
######################################################################################################
control.disease.combined<- Add_Cell_Complexity(object = control.disease.combined, overwrite = TRUE)
all.equal(colnames(control.disease.combined), row.names(control.disease.combined@meta.data))

##### create plots of QC metrics before filtering
pdf(paste0(result_dir, "/QC_Control/QC_pre_filtering.pdf"), width = 20, height = 3.5)
## Unique Molecular Identifiers, absolute number of observed transcripts
p_preQC11 <- QC_Plots_Genes(seurat_object = control.disease.combined, group.by = "samples", plot_title = "Genes Per Cell") + scale_fill_manual(values = color_list)
p_preQC21 <- QC_Plots_UMIs(seurat_object = control.disease.combined, group.by = "samples", plot_title = "UMIs Per Cell") + scale_fill_manual(values = color_list)
p_preQC31 <- QC_Plots_Mito(seurat_object = control.disease.combined, group.by = "samples", plot_title = "Mito Gene % Per Cell") + scale_fill_manual(values = color_list)
p_preQC41 <- QC_Plots_Feature(seurat_object = control.disease.combined, group.by = "samples", feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell") + scale_fill_manual(values = color_list)
wrap_plots(p_preQC11, p_preQC21, p_preQC31, p_preQC41, ncol = 4)

p_preQC12 <- QC_Plots_Genes(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "Genes Per Cell") + scale_fill_manual(values = color_list)
p_preQC22 <- QC_Plots_UMIs(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "UMIs Per Cell") + scale_fill_manual(values = color_list)
p_preQC32 <- QC_Plots_Mito(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "Mito Gene % Per Cell") + scale_fill_manual(values = color_list)
p_preQC42 <- QC_Plots_Feature(seurat_object = control.disease.combined, group.by = "donor_condition", feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell") + scale_fill_manual(values = color_list)
wrap_plots(p_preQC12, p_preQC22, p_preQC32, p_preQC42, ncol = 4)
dev.off()

################ Cell counts: The cell counts are determined by the number of unique cellular barcodes detected.
pdf(paste0(result_dir, "/QC_Control/nCells_per_sample_donor_condition.pdf"), width = 20, height = 6)
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
pdf(paste0(result_dir, "/QC_Control/nUMI_per_sample_donor_condition.pdf"), width = 20, height = 6)
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
pdf(paste0(result_dir, "/QC_Control/nGene_per_sample_donor_condition.pdf"), width = 20, height = 6)
## Visualize the distribution of genes detected per cell via histogram
p_nGene_1 <- metadata %>% ggplot(aes(x=nGene, color = samples, fill = samples)) + geom_density(alpha = 0.2) + scale_x_log10() + 
  scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + 
  geom_vline(xintercept = 300) + ggtitle("nGenes per cell for each sample") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(legend.title = element_blank())
p_nGene_2 <- metadata %>% ggplot(aes(x=nGene, color = donor_condition, fill = donor_condition)) + geom_density(alpha = 0.2) + scale_x_log10() + 
  scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + 
  geom_vline(xintercept = 300) + ggtitle("nGenes per cell for each donor") + theme(plot.title = element_text(hjust = 0.5, face="bold")) + theme(legend.title = element_blank())
# Visualize the distribution of genes detected per cell via boxplot
p_nGene_3 <- metadata %>% ggplot(aes(x=condition, y=log10(nGene), fill=condition)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = Ctrl_IHD_color) + 
  theme(legend.title=element_blank()) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.title = element_text(hjust=0.5, face="bold")) + ggtitle("nGenes per cell for each condition")
wrap_plots(p_nGene_1, p_nGene_2, p_nGene_3, ncol = 3)
dev.off()

##### Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata$log10GenesPerUMI <- control.disease.combined$log10GenesPerUMI
pdf(paste0(result_dir, "/QC_Control/Genes_per_UMI_sample_donor_condition.pdf"), width = 20, height = 6)
p_Genes_per_UMI_1 <- metadata %>% ggplot(aes(x = log10GenesPerUMI, color = samples, fill = samples)) + geom_density(alpha = 0.2) +
  scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + theme(legend.title=element_blank()) + geom_vline(xintercept = 0.8)
p_Genes_per_UMI_2 <- metadata %>% ggplot(aes(x = log10GenesPerUMI, color = donor_condition, fill = donor_condition)) + geom_density(alpha = 0.2) +
  scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + theme(legend.title=element_blank()) + geom_vline(xintercept = 0.8)
p_Genes_per_UMI_3 <- metadata %>% ggplot(aes(x = log10GenesPerUMI, color = condition, fill = condition)) + geom_density(alpha = 0.2) +
  scale_fill_manual(values = color_list) + scale_color_manual(values = color_list) + theme_classic() + theme(legend.title=element_blank()) + geom_vline(xintercept = 0.8)
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

pdf(paste0(result_dir, "/QC_Control/cells_filtered_upset_plot.pdf"))
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
pdf(paste0(result_dir, "/QC_Control/QC_post_filtering.pdf"), width = 20, height = 3.5)
  ## Unique Molecular Identifiers, absolute number of observed transcripts
  p_postQC11 <- QC_Plots_Genes(seurat_object = control.disease.combined, group.by = "samples", plot_title = "Genes Per Cell") + scale_fill_manual(values = color_list)
  p_postQC21 <- QC_Plots_UMIs(seurat_object = control.disease.combined, group.by = "samples", plot_title = "UMIs Per Cell") + scale_fill_manual(values = color_list)
  p_postQC31 <- QC_Plots_Mito(seurat_object = control.disease.combined, group.by = "samples", plot_title = "Mito Gene % Per Cell") + scale_fill_manual(values = color_list)
  p_postQC41 <- QC_Plots_Feature(seurat_object = control.disease.combined, group.by = "samples", feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell") + scale_fill_manual(values = color_list)
  wrap_plots(p_postQC11, p_postQC21, p_postQC31, p_postQC41, ncol = 4)
  
  p_postQC12 <- QC_Plots_Genes(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "Genes Per Cell") + scale_fill_manual(values = color_list)
  p_postQC22 <- QC_Plots_UMIs(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "UMIs Per Cell") + scale_fill_manual(values = color_list)
  p_postQC32 <- QC_Plots_Mito(seurat_object = control.disease.combined, group.by = "donor_condition", plot_title = "Mito Gene % Per Cell") + scale_fill_manual(values = color_list)
  p_postQC42 <- QC_Plots_Feature(seurat_object = control.disease.combined, group.by = "donor_condition", feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell") + scale_fill_manual(values = color_list)
  wrap_plots(p_postQC12, p_postQC22, p_postQC32, p_postQC42, ncol = 4)
dev.off()

saveRDS(control.disease.combined, paste0(result_dir, "/merged.Seurat.obj_filtered.RDS"))

######################################################################################################
######################################################################################################
################## Normalization, scaling, dimensionality reduction, clustering ######################
######################################################################################################
######################################################################################################
# control.disease.combined <- readRDS(paste0(result_dir, "/merged.Seurat.obj_filtered.RDS"))
n_cells <- ncol(control.disease.combined)

DefaultAssay(control.disease.combined) <- 'RNA'
control.disease.combined <- NormalizeData(control.disease.combined)
control.disease.combined <- FindVariableFeatures(control.disease.combined, selection.method = "vst", nfeatures = 5000)
control.disease.combined <- ScaleData(control.disease.combined, vars.to.regress = c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA"))

saveRDS(control.disease.combined, paste0(result_dir, "/merged.Seurat.obj_filtered_Scaled_01.RDS"))
# control.disease.combined <- readRDS(paste0(result_dir, "/merged.Seurat.obj_filtered_Scaled_01.RDS"))
# check_PC <- RunPCA(control.disease.combined, verbose = TRUE)
# ElbowPlot(check_PC, ndims = 50, reduction = "pca")

control.disease.combined <- RunPCA(control.disease.combined, npcs = 30, verbose = TRUE) %>% 
  # RunHarmony(group.by.vars = "samples", lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  RunHarmony(group.by.vars = c("samples", "technique"), lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  RunTSNE(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = TRUE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = TRUE) %>%
  FindClusters(resolution = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6), verbose = TRUE)

pdf(paste0(result_dir, "/Annotation/Dim_Plot_before_doublet.pdf"), width = 7, height = 5)
  ElbowPlot(control.disease.combined, ndims = 30, reduction = "pca")
  DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE, group.by = 'samples') + ggtitle("PCA")
  DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE, group.by = 'condition') + ggtitle("PCA")
  DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE, group.by = 'technique') + ggtitle("PCA")
  DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE, group.by = 'samples') + ggtitle("TSNE")
  DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE, group.by = 'condition') + ggtitle("TSNE")
  DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE, group.by = 'technique') + ggtitle("TSNE")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'samples') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'condition') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'technique') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'samples') + ggtitle("Harmony")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'condition') + ggtitle("Harmony")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'technique') + ggtitle("Harmony")
  
  DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE) + ggtitle("PCA")
  DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE) + ggtitle("TSNE")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE) + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE) + ggtitle("Harmony")
dev.off()

saveRDS(control.disease.combined, paste0(result_dir, "/combined.Seurat.obj_before_doublet.RDS"))
control.disease.combined <- readRDS(paste0(result_dir, "/combined.Seurat.obj_before_doublet.RDS"))

######################################################################################################
######################################################################################################
########################################## Doublet finder ############################################
######################################################################################################
######################################################################################################
sweep.res.list_combined <- paramSweep_v3(control.disease.combined, PCs = 1:30, sct = FALSE)
sweep.stats_combined <- summarizeSweep(sweep.res.list_combined, GT = FALSE)
bcmvn_combined <- find.pK(sweep.stats_combined)

##### Plot the bcmvn output to find the peak
pdf(paste0(result_dir, "/QC_Control/Optimal_pK_doubletFinder.pdf"), width = 12, height = 5)
ggplot(bcmvn_combined, aes(x = pK, y = BCmetric)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Optimal pK Value", x = "pK", y = "BCmvn")
dev.off()

homotypic.prop <- modelHomotypic(control.disease.combined@meta.data$seurat_clusters) 
nExp_poi <- round(0.075*nrow(control.disease.combined@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# control.disease.combined <- doubletFinder_v3(control.disease.combined, PCs = 1:30, pK = 0.040, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
control.disease.combined <- doubletFinder_v3(control.disease.combined, PCs = 1:30, pK = 0.060, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
head(control.disease.combined)
saveRDS(control.disease.combined, paste0(result_dir, "/combined.Seurat.obj_doubletFinder_v3.RDS"))
# control.disease.combined  <- readRDS(paste0(result_dir, "/combined.Seurat.obj_doubletFinder_v3.RDS"))

pdf(paste0(result_dir, "/QC_Control/DF_classifications.pdf"), width = 7, height = 5)
DimPlot(control.disease.combined, reduction = "umap", group.by = "DF.classifications_0.25_0.06_5066")
dev.off()

table(control.disease.combined@meta.data$DF.classifications_0.25_0.06_5066)
# Doublet Singlet 
# 5078   66630
# 5145   66563
# 5119   66589 
# 5066   66642

control.disease.combined <- subset(control.disease.combined, subset = DF.classifications_0.25_0.06_5066 == "Singlet")
# saveRDS(control.disease.combined, paste0(result_dir, "/combined.Seurat.obj_remove_doublet.RDS"))
# control.disease.combined  <- readRDS(paste0(result_dir, "/combined.Seurat.obj_remove_doublet.RDS"))

control.disease.combined <- NormalizeData(control.disease.combined)
control.disease.combined <- FindVariableFeatures(control.disease.combined, selection.method = "vst", nfeatures = 5000)
control.disease.combined <- ScaleData(control.disease.combined, vars.to.regress = c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA"))

saveRDS(control.disease.combined, paste0(result_dir, "/merged.Seurat.obj_filtered_Scaled_02.RDS"))
# control.disease.combined <- readRDS(paste0(result_dir, "/merged.Seurat.obj_filtered_Scaled_02.RDS"))
# check_PC_02 <- RunPCA(control.disease.combined, verbose = TRUE)
# ElbowPlot(check_PC_02, ndims = 50, reduction = "pca")

control.disease.combined <- RunPCA(control.disease.combined, npcs = 30, verbose = TRUE) %>% 
  # RunHarmony(group.by.vars = "samples", lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  RunHarmony(group.by.vars = c("samples", "technique"), lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  RunTSNE(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  FindClusters(resolution = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6), verbose = TRUE)

pdf(paste0(result_dir, "/Annotation/Dim_Plot_after_doublet.pdf"), width = 7, height = 5)
  ElbowPlot(control.disease.combined, ndims = 30, reduction = "pca")
  DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE, group.by = 'samples') + ggtitle("PCA")
  DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE, group.by = 'condition') + ggtitle("PCA")
  DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE, group.by = 'technique') + ggtitle("PCA")
  DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE, group.by = 'samples') + ggtitle("TSNE")
  DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE, group.by = 'condition') + ggtitle("TSNE")
  DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE, group.by = 'technique') + ggtitle("TSNE")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'samples') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'condition') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'technique') + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'samples') + ggtitle("Harmony")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'condition') + ggtitle("Harmony")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'technique') + ggtitle("Harmony")
  
  DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE) + ggtitle("PCA")
  DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE) + ggtitle("TSNE")
  DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE) + ggtitle("UMAP")
  DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE) + ggtitle("Harmony")
dev.off()

saveRDS(control.disease.combined, paste0(result_dir, "/combined.Seurat.obj_ready_for_annotation.RDS"))

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
### continue from O2
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

##### clustering results summary
control.disease.combined  <- readRDS(paste0(result_dir, "/combined.Seurat.obj_ready_for_annotation.RDS"))
sample_names <- unique(control.disease.combined$samples)
RNA_resolution <- "RNA_snn_res.0.8"
Idents(object = control.disease.combined) <- RNA_resolution

pdf(paste0(result_dir, "/Annotation/Unannotated_clustering.pdf"), width = 6, height = 5)
DimPlot(control.disease.combined, reduction = 'umap', group.by = RNA_resolution, label = TRUE, label.size = 3, cols = color_list, order = T)
DimPlot(control.disease.combined, reduction = 'tsne', group.by = RNA_resolution, label = TRUE, label.size = 3, cols = color_list, order = T)
dev.off()

pdf(paste0(result_dir, "/Annotation/Unannotated_clustering_split_sample.pdf"), width = 24, height = 5)
DimPlot(control.disease.combined, reduction = 'umap', split.by = "samples", label = TRUE, label.size = 3, cols = color_list, order = T) + 
  facet_grid(. ~  factor(samples, level = sample_names), scale = "free_x")
DimPlot(control.disease.combined, reduction = 'tsne', split.by = "samples", label = TRUE, label.size = 3, cols = color_list, order = T) + 
  facet_grid(. ~  factor(samples, level = sample_names), scale = "free_x")
dev.off()

pdf(paste0(result_dir, "/Annotation/Unannotated_clustering_split_condition.pdf"), width = 10, height = 5)
DimPlot(control.disease.combined, reduction = 'umap', split.by = "condition", label = TRUE, label.size = 3, cols = color_list, order = T)
DimPlot(control.disease.combined, reduction = 'tsne', split.by = "condition", label = TRUE, label.size = 3, cols = color_list, order = T)
dev.off()

################################################################################################################################################################
### cell type annotation
################################################################################################################################################################
# feature_gene_list <- c("TNNT2", "ACTN2", "TTN", "MYOZ2", "RYR2", "MYPN", "MLIP", #CM 
#                        "VWF", "FLT1", "PECAM1", "FABP5", "ID1", "CDH5", "CD9", # Endo 
#                        "DCN", "APOD", "GSN", "BICC1", "ABCA6", "CDH19", # fib 
#                        "RGS5", "EPS8", "PDGFRB", "EGFLAM", "GUCY1A2", # Pericytes
#                        "NRXN1", "XKR4", "NRXN3", "CADM2", "KIRREL3", "ERBB3", "PCSK2",  # neurons
#                        "F13A1", "MRC1", "RBM47", "MS4A6A", "ATP8B4", # Myeloid
#                        "PTPRC", "SKAP1", "AOAH", "CD247", "IKZF1",  # NK-Cells
#                        "SLC24A3", "IL18R1", "KIT", "CPA3") # mast

marker_set1 <- c("TNNT2", "ACTN2", "MYOZ2", "RYR2", "MYPN", "MLIP")
marker_set2 <- c("VWF", "FLT1", "PECAM1", "FABP5", "ID1", "CDH5", "CD9")
marker_set3 <- c("DCN", "APOD", "GSN", "BICC1", "ABCA6", "CDH19")
marker_set4 <- c("RGS5", "EPS8", "PDGFRB", "EGFLAM", "GUCY1A2")
marker_set5 <- c("NRXN1", "XKR4", "NRXN3", "CADM2", "KIRREL3", "ERBB3", "PCSK2")
marker_set6 <- c("F13A1", "MRC1", "RBM47", "MS4A6A", "ATP8B4")
marker_set7 <- c("PTPRC", "SKAP1", "AOAH", "CD247", "IKZF1")
marker_set8 <- c("SLC24A3", "IL18R1", "KIT", "CPA3")
marker_set9 <- c("PPARG", "ADIPOQ", "MEST", "GPAM", "FASN")

all_markers <- c(marker_set1, marker_set2, marker_set3, marker_set4, marker_set5, marker_set6, marker_set7, marker_set8, marker_set9)
marker_groups <- c(rep("CM markers", length(marker_set1)),        rep("Endothelial markers", length(marker_set2)), rep("Fibroblasts markers", length(marker_set3)),
                   rep("Pericytes markers", length(marker_set4)), rep("Neurons markers", length(marker_set5)),     rep("Myeloid markers", length(marker_set6)),
                   rep("NK-Cell markers", length(marker_set7)),   rep("Mast markers", length(marker_set8)),        rep("Adipocyte markers", length(marker_set9)))

gene_group_df <- data.frame(gene = all_markers, group = marker_groups) %>% 
  mutate(group = factor(group, level = c("CM markers", "Endothelial markers", "Fibroblasts markers", 
                                         "Pericytes markers", "Neurons markers", "Myeloid markers", 
                                         "NK-Cell markers", "Mast markers", "Adipocyte markers")))

pdf(paste0(result_dir, "/Annotation/Unannotated_dotplot.pdf"), width = 20, height = 10)
dotplot_all <- DotPlot(control.disease.combined, features = all_markers, group.by = RNA_resolution) + RotatedAxis() + scale_color_gradient(low = "lightblue", high = "darkblue")
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all <- dotplot_all + facet_grid(~ gene_group, scales = "free_x", space = "free_x")
dotplot_all

control.disease.combined_Control <- subset(control.disease.combined, subset = condition == "Control")
dotplot_Control <- DotPlot(control.disease.combined_Control, features = all_markers, group.by = RNA_resolution) + RotatedAxis() + scale_color_gradient(low = "lightblue", high = "darkblue")
dotplot_Control$data$gene_group <- factor(gene_group_df$group[match(dotplot_Control$data$features.plot, gene_group_df$gene)])
dotplot_Control <- dotplot_Control + facet_grid(~ gene_group, scales = "free_x", space = "free_x")
dotplot_Control

control.disease.combined_Disease <- subset(control.disease.combined, subset = condition == "Disease")
dotplot_Disease <- DotPlot(control.disease.combined_Disease, features = all_markers, group.by = RNA_resolution) + RotatedAxis() + scale_color_gradient(low = "lightblue", high = "darkblue")
dotplot_Disease$data$gene_group <- factor(gene_group_df$group[match(dotplot_Disease$data$features.plot, gene_group_df$gene)])
dotplot_Disease <- dotplot_Disease + facet_grid(~ gene_group, scales = "free_x", space = "free_x")
dotplot_Disease
dev.off()

### top10 markers based on FC in each cluster and annotate cell types
# cluster0.markers <- FindMarkers(control.disease.combined, ident.1 = 0, min.pct = 0.25)
# print(slice_max(cluster0.markers, n = 15, order_by = avg_log2FC))
# top10_FC_markers_each_cluster <- FindAllMarkers(control.disease.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.50)
# write.csv(top10_FC_markers_each_cluster, paste0(result_dir, "/Marker_files/all_markers_each_cluster_unsupervised.csv"), row.names=FALSE)
# 
# top10_FC_markers_each_cluster %>%
#   group_by(cluster) %>%
#   slice_max(n = 30, order_by = avg_log2FC) %>%
#   print(n=2500) %>%
#   write.csv(paste0(result_dir, "/Marker_files/top10_FC_markers_each_cluster_unsupervised.csv"), row.names=FALSE)

############ manually annotation
annotated_clusters_list <- c(
  "Cardiomyocytes", #0 ca
  "Endothelial",    #1 en
  "Fibroblasts",    #2 fi
  "Pericytes",      #3 pe 
  "Fibroblasts",    #4 fi
  "Endothelial",    #5 en
  "Cardiomyocytes", #6 ca
  "Endothelial",    #7 en
  "Endothelial",    #8 en
  "Pericytes",      #9 pe
  "Cardiomyocytes", #10 ca
  "Neurons",        #11 neuron
  "Cardiomyocytes", #12 ca and en
  "Myeloid",        #13 Myeloid
  "NK-Cells",       #14 NK-T-Cells
  "Cardiomyocytes", #15 ca
  "Cardiomyocytes", #16 ca
  "Pericytes",      #17 pe
  "Cardiomyocytes", #18 ca
  "Endothelial",    #19 en
  "Cardiomyocytes", #20 ca
  "Mast",           #21 MAST
  "Cardiomyocytes", #22 ca
  "Endothelial",    #23 en
  "Adipocyte"            #24
)

names(annotated_clusters_list) <- paste0("cluster_", 0:(length(annotated_clusters_list) - 1))
annotated_clusters <- sapply(control.disease.combined$RNA_snn_res.0.8, function(x) annotated_clusters_list[paste0("cluster_", x)])
control.disease.combined@meta.data <- cbind(control.disease.combined@meta.data, annotated_clusters)
cell_type_order <- c("Cardiomyocytes", "Fibroblasts", "Endothelial", "Pericytes", "Neurons", "Myeloid", "NK-Cells", "Mast", "Adipocyte")
## Reorder the cell types in the metadata
control.disease.combined$annotated_clusters <- factor(control.disease.combined$annotated_clusters, levels = cell_type_order)

pdf(paste0(result_dir, "/Annotation/Annotated_clusters.pdf"), width = 8, height = 6)
DimPlot(control.disease.combined, reduction = 'umap', group.by = "annotated_clusters", label = T, label.size = 3, cols = color_list, order = T) 
DimPlot(control.disease.combined, reduction = 'tsne', group.by = "annotated_clusters", label = T, label.size = 3, cols = color_list, order = T)
dev.off()

pdf(paste0(result_dir, "/Annotation/Annotated_clustering_split_condition.pdf"), width = 12, height = 6)
DimPlot(control.disease.combined, reduction = 'umap', group.by = "annotated_clusters", split.by = "condition", label = T, label.size = 3, cols = color_list, order = T) 
DimPlot(control.disease.combined, reduction = 'tsne', group.by = "annotated_clusters", split.by = "condition", label = T, label.size = 3, cols = color_list, order = T) 
dev.off()

pdf(paste0(result_dir, "/Annotation/Annotated_clustering_split_sample.pdf"), width = 30, height = 5)
DimPlot(control.disease.combined, reduction = 'umap', group.by = "annotated_clusters", split.by = "samples", 
        label = T, label.size = 3, cols = color_list, order = T) + facet_grid(. ~  factor(samples, level = sample_names), scale = "free_x")
DimPlot(control.disease.combined, reduction = 'tsne', group.by = "annotated_clusters", split.by = "samples",
        label = T, label.size = 3, cols = color_list, order = T) + facet_grid(. ~  factor(samples, level = sample_names), scale = "free_x")
dev.off()

gene_group_df <- data.frame(gene = all_markers, group = marker_groups) %>% 
  mutate(group = factor(group, level = c("CM markers", "Fibroblasts markers", "Endothelial markers", 
                                         "Pericytes markers", "Neurons markers", "Myeloid markers", 
                                         "NK-Cell markers", "Mast markers", "Adipocyte markers")))

pdf(paste0(result_dir, "/Annotation/Annotated_dotplot.pdf"), width = 18, height = 8)
dotplot_all <- DotPlot(control.disease.combined, features = all_markers, group.by = "annotated_clusters") + RotatedAxis() + scale_color_gradient(low = "lightblue", high = "darkblue")
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all <- dotplot_all + facet_grid(~ gene_group, scales = "free_x", space = "free_x")
dotplot_all

control.disease.combined_Control <- subset(control.disease.combined, subset = condition == "Control")
dotplot_Control <- DotPlot(control.disease.combined_Control, features = all_markers, group.by = "annotated_clusters") + RotatedAxis() + scale_color_gradient(low = "lightblue", high = "darkblue")
dotplot_Control$data$gene_group <- factor(gene_group_df$group[match(dotplot_Control$data$features.plot, gene_group_df$gene)])
dotplot_Control <- dotplot_Control + facet_grid(~ gene_group, scales = "free_x", space = "free_x")
dotplot_Control

control.disease.combined_Disease <- subset(control.disease.combined, subset = condition == "Disease")
dotplot_Disease <- DotPlot(control.disease.combined_Disease, features = all_markers, group.by = "annotated_clusters") + RotatedAxis() + scale_color_gradient(low = "lightblue", high = "darkblue")
dotplot_Disease$data$gene_group <- factor(gene_group_df$group[match(dotplot_Disease$data$features.plot, gene_group_df$gene)])
dotplot_Disease <- dotplot_Disease + facet_grid(~ gene_group, scales = "free_x", space = "free_x")
dotplot_Disease
dev.off()

################################################################################################################################################################
### statistic for unannotated and annotated clustering
################################################################################################################################################################
num_per_sample_cluster_unannotated <- table(control.disease.combined$samples, control.disease.combined$RNA_snn_res.0.8)
num_per_sample_cluster_unannotated <- num_per_sample_cluster_unannotated[unique(control.disease.combined$samples), ]

percentage_per_sample_cluster_unannotated <- num_per_sample_cluster_unannotated %>% 
  as.data.frame.matrix() %>% mutate(Total = rowSums(.)) %>% 
  mutate(across(-Total, ~ round(.x / Total, 5))) %>% select(-Total) %>% 
  mutate(condition = c("Control", "Control", "Control", "Control", "Disease", "Disease", "Disease", "Disease", "Disease"))

percentage_per_sample_cluster_unannotated_melt <- percentage_per_sample_cluster_unannotated %>% melt() %>% 
  setNames(c("condition", "cluster", "ratio")) %>% mutate(condition = factor(condition, level = c("Control", "Disease"))) %>% mutate(ratio = as.numeric(ratio))

pdf(paste0(result_dir, "/Annotation/Unannotated_cluster_ratio_comparison.pdf"), width = 26, height = 6)
barplot_strand_ratio_age_match <- ggbarplot(percentage_per_sample_cluster_unannotated_melt, x = "cluster", y = "ratio", color = "condition", 
                                            label = TRUE, lab.nb.digits = 2, add = c("mean_se", "jitter"), position = position_dodge(0.9), palette = c("dodgerblue3", "#FC4E07")) +
  stat_compare_means(aes(group = condition), method = "wilcox", label = "p.format", label.y = 1.02 * max(percentage_per_sample_cluster_unannotated_melt$ratio, na.rm = TRUE)) + 
  stat_compare_means(aes(group = condition), method = "wilcox", label = "p.signif", label.y = 1.06 * max(percentage_per_sample_cluster_unannotated_melt$ratio, na.rm = TRUE)) + 
  labs(x = "Unannotated Cluster", y = "ratio of cell #") + theme_classic(base_size = 24) + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values=c("dodgerblue3", "#FC4E07"))
print(barplot_strand_ratio_age_match)
dev.off()


num_per_sample_cluster_annotated <- table(control.disease.combined$samples, control.disease.combined$annotated_clusters)
num_per_sample_cluster_annotated <- num_per_sample_cluster_annotated[unique(control.disease.combined$samples), ]

percentage_per_sample_cluster_annotated <- num_per_sample_cluster_annotated %>% 
  as.data.frame.matrix() %>% mutate(Total = rowSums(.)) %>% 
  mutate(across(-Total, ~ round(.x / Total *100 , 5))) %>% select(-Total) %>% 
  mutate(condition = c("Control", "Control", "Control", "Control", "Disease", "Disease", "Disease", "Disease", "Disease"))

percentage_per_sample_cluster_annotated_melt <- percentage_per_sample_cluster_annotated %>% melt() %>% 
  setNames(c("condition", "cluster", "ratio")) %>% mutate(condition = factor(condition, level = c("Control", "Disease"))) %>% mutate(ratio = as.numeric(ratio))

pdf(paste0(result_dir, "/Annotation/Annotated_cluster_ratio_comparison.pdf"), width = 20, height = 6)
barplot_strand_ratio_age_match <- ggbarplot(percentage_per_sample_cluster_annotated_melt, x = "cluster", y = "ratio", color = "condition", 
                                            label = TRUE, lab.nb.digits = 2, add = c("mean_se", "jitter"), position = position_dodge(0.9), palette = c("dodgerblue3", "#FC4E07")) +
  stat_compare_means(aes(group = condition), method = "wilcox", label = "p.format", label.y = 1.02 * max(percentage_per_sample_cluster_annotated_melt$ratio, na.rm = TRUE)) + 
  stat_compare_means(aes(group = condition), method = "wilcox", label = "p.signif", label.y = 1.06 * max(percentage_per_sample_cluster_annotated_melt$ratio, na.rm = TRUE)) + 
  labs(x = "Unannotated Cluster", y = "ratio of cell #") + theme_classic(base_size = 24) + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values=c("dodgerblue3", "#FC4E07"))
print(barplot_strand_ratio_age_match)
dev.off()

## Cell Type Composition
donor_conditions <- data.frame(Condition = percentage_per_sample_cluster_annotated$condition)
rownames(donor_conditions) <- rownames(num_per_sample_cluster_annotated)

pdf(paste0(result_dir, "/Annotation/cell_count_in_annotated_cluster.pdf"), width = 10, height = 6)
pheatmap(num_per_sample_cluster_annotated, annotation_row = donor_conditions, 
         display_numbers = TRUE, number_format = "%.0f", number_color = "black", fontsize_number = 15,
         cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

pdf(paste0(result_dir, "/Annotation/cell_percentage_in_annotated_cluster.pdf"), width = 10, height = 6)
pheatmap(percentage_per_sample_cluster_annotated[,-10], annotation_row = donor_conditions, 
         display_numbers = TRUE, number_format = "%.1f", number_color = "black", fontsize_number = 15,
         cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

## save object
control.disease.combined$celltype.condition <- paste(control.disease.combined$annotated_clusters, control.disease.combined$condition, sep = "_")
control.disease.combined$samples.condition <- paste0(control.disease.combined$samples, "-", control.disease.combined$condition)
saveRDS(control.disease.combined, paste0(result_dir, "/Seurat.obj_with_annotation.RDS"))

