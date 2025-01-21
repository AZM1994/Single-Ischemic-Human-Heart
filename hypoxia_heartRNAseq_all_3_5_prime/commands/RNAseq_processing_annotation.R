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
color_list <- c("red", "#2ca02c", "dodgerblue", "#ff7f0e", "#9467bd", "#e377c2", "#8c564b", "#7f7f7f", "#bcbd22", "#17becf",
                "#393b79", "#5254a3", "#6b6ecf", "#9c9ede", "#637939", "#8ca252", "#b5cf6b", "#cedb9c", "#8c6d31", "#bd9e39",
                "#e7ba52", "#e7969c", "#d6616b", "#ad494a", "#843c39", "#f5f5f5", "#fdd0a2", "#fb6a4a", "#cb181d", "#a50f15",
                "#3182bd", "#6baed6", "#9ecae1", "#c6dbef", "#e6550d", "#fd8d3c", "#fdae6b", "#fdd0a2", "#31a354", "#74c476",
                "#a1d99b", "#c7e9c0", "#756bb1", "#9e9ac8", "#bcbddc")
color_list_02 <- c("#f22c4b", "#9c9ede", "#fd8d3c", "#73B8E1", "#bcbd22", "#e377c2", "#8c564b", "#7f7f7f")
##### set working and results directories
wd_dir <- getwd()
result_dir <- paste0(wd_dir, "/results")
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
control.disease.combined <- FindVariableFeatures(control.disease.combined, selection.method = "vst", nfeatures = 3000)
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

control.disease.combined <- doubletFinder_v3(control.disease.combined, PCs = 1:30, pK = 0.040, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
# control.disease.combined <- doubletFinder_v3(control.disease.combined, PCs = 1:30, pK = 0.060, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
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
control.disease.combined <- FindVariableFeatures(control.disease.combined, selection.method = "vst", nfeatures = 3000)
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

######################################################################################################
######################################################################################################
######################################### continue from O2 ###########################################
######################################################################################################
######################################################################################################
##### clustering results summary
control.disease.combined  <- readRDS(paste0(result_dir, "/combined.Seurat.obj_ready_for_annotation.RDS"))

# Calculate the average number of cells and genes per sample
average_metrics <- control.disease.combined@meta.data %>%
  group_by(samples) %>% summarise(avg_cells_per_sample = n(), avg_genes_per_cell = mean(nFeature_RNA))

total_avg_cells <- nrow(control.disease.combined@meta.data) / length(unique(control.disease.combined@meta.data$samples))
total_avg_genes <- mean(control.disease.combined@meta.data$nFeature_RNA)
num_nuclei <- nrow(control.disease.combined@meta.data)

selected_sample_names <- unique(control.disease.combined$samples)
selected_donor_names <- unique(control.disease.combined$donor_condition)
RNA_resolution <- "RNA_snn_res.1"
Idents(object = control.disease.combined) <- RNA_resolution

dir.create(paste0(result_dir, "/Annotation/Unannotated_clustering"))
pdf(paste0(result_dir, "/Annotation/Unannotated_clustering/Unannotated_clustering.pdf"), width = 6.5, height = 5)
  DimPlot(control.disease.combined, reduction = 'umap', group.by = RNA_resolution, label = TRUE, label.size = 3, cols = color_list, order = T) + ggtitle(NULL) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

pdf(paste0(result_dir, "/Annotation/Unannotated_clustering/Unannotated_clustering_split_condition.pdf"), width = 15, height = 8)
  DimPlot(control.disease.combined, reduction = 'umap', split.by = "condition", label = TRUE, label.size = 4, pt.size = 0.5, cols = color_list, order = T) + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 16))
dev.off()

pdf(paste0(result_dir, "/Annotation/Unannotated_clustering/Unannotated_clustering_split_sample.pdf"), width = 50, height = 4)
  DimPlot(control.disease.combined, reduction = 'umap', split.by = "samples", label = TRUE, label.size = 3, cols = color_list, order = T) + 
    facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x") + theme_linedraw() + ggtitle(NULL) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

pdf(paste0(result_dir, "/Annotation/Unannotated_clustering/Unannotated_clustering_split_donor.pdf"), width = 28, height = 4)
  DimPlot(control.disease.combined, reduction = 'umap', split.by = "donor_condition", label = TRUE, label.size = 3, cols = color_list, order = T) + 
    facet_grid(. ~  factor(donor_condition, level = selected_donor_names), scale = "free_x") + theme_linedraw() + ggtitle(NULL) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

################################################################################################################################################################
### cell type annotation
################################################################################################################################################################
marker_set1 <- c("TNNT2", "ACTN2", "MYOZ2", "MYPN", "MLIP") # Cardiomyocyte
marker_set2 <- c("DCN", "APOD", "BICC1", "ABCA6") # Fibroblast
marker_set3 <- c("VWF", "PECAM1", "FABP5", "ID1") # Endothelial
marker_set4 <- c("RGS5", "PDGFRB", "EGFLAM", "GUCY1A2") # Pericyte
marker_set5 <- c("NRXN1", "NRXN3", "CADM2", "KIRREL3") # Neuronal
marker_set6 <- c("PTPRC", "DOCK2", "IQGAP2", "RUNX1") # Immune
marker_set7 <- c("IL18R1", "KIT", "CPA3") # Mast

all_markers <- c(marker_set1, marker_set2, marker_set3, marker_set4, marker_set5, marker_set6, marker_set7)
marker_groups <- c(rep("Cardiomyocyte markers", length(marker_set1)), rep("Fibroblast markers", length(marker_set2)), rep("Endothelial markers", length(marker_set3)), rep("Pericyte markers", length(marker_set4)), 
                   rep("Neuronal markers", length(marker_set5)), rep("Immune markers", length(marker_set6)), rep("Mast markers", length(marker_set7)))

gene_group_df <- data.frame(gene = all_markers, group = marker_groups) %>% 
  mutate(group = factor(group, level = c("Cardiomyocyte markers", "Fibroblast markers", "Endothelial markers", "Pericyte markers", "Neuronal markers", "Immune markers", "Mast markers")))

dir.create(paste0(result_dir, "/Annotation/Unannotated_dotplot"))
control.disease.combined_Control <- subset(control.disease.combined, subset = condition == "Control")
control.disease.combined_Disease <- subset(control.disease.combined, subset = condition == "IHD")

dotplot_all <- DotPlot(control.disease.combined, features = all_markers, group.by = RNA_resolution)
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all <- dotplot_all + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + scale_color_gradient(low = "white", high = "black") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "black", linetype = "dashed"), text = element_text(face = "bold", size = 12)) + RotatedAxis()
ggsave(paste0(result_dir, "/Annotation/Unannotated_dotplot/Unannotated_dotplot_all.pdf"), plot = dotplot_all, width = 13, height = 8, dpi = 300)

dotplot_Control <- DotPlot(control.disease.combined_Control, features = all_markers, group.by = RNA_resolution)
dotplot_Control$data$gene_group <- factor(gene_group_df$group[match(dotplot_Control$data$features.plot, gene_group_df$gene)])
dotplot_Control$data$condition <- "Control"
dotplot_Control <- dotplot_Control + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + scale_color_gradient(low = "white", high = "dodgerblue4") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "dodgerblue4", linetype = "dashed"), text = element_text(face = "bold", size = 12)) + RotatedAxis()
ggsave(paste0(result_dir, "/Annotation/Unannotated_dotplot/Unannotated_dotplot_Control.pdf"), plot = dotplot_Control, width = 13, height = 8, dpi = 300)
  
dotplot_Disease <- DotPlot(control.disease.combined_Disease, features = all_markers, group.by = RNA_resolution)
dotplot_Disease$data$gene_group <- factor(gene_group_df$group[match(dotplot_Disease$data$features.plot, gene_group_df$gene)])
dotplot_Disease$data$condition <- "IHD"
dotplot_Disease <- dotplot_Disease + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + scale_color_gradient(low = "white", high = "firebrick4") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "firebrick4", linetype = "dashed"), text = element_text(face = "bold", size = 12)) + RotatedAxis() + guides(colour = guide_colourbar())
ggsave(paste0(result_dir, "/Annotation/Unannotated_dotplot/Unannotated_dotplot_Disease.pdf"), plot = dotplot_Disease, width = 13, height = 8, dpi = 300)

### top10 markers based on FC in each cluster and annotate cell types
# cluster28.markers <- FindMarkers(control.disease.combined, ident.1 = 28, min.pct = 0.50) %>% arrange(desc(avg_log2FC))
# row.names(cluster28.markers)[1:50]
# top10_FC_markers_each_cluster <- FindAllMarkers(control.disease.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.50)
# write.csv(top10_FC_markers_each_cluster, paste0(result_dir, "/Marker_files/all_markers_each_cluster_unsupervised.csv"), row.names=FALSE)

##### manually annotation
# annotated_clusters_list <- c("Endothelial", "Cardiomyocyte", "Fibroblast", "Fibroblast", "Pericyte", "Endothelial", "Cardiomyocyte", "Pericyte", "Endothelial", "Cardiomyocyte", 
#                              "Endothelial", "Pericyte", "Fibroblast", "Neuronal", "Cardiomyocyte", "Cardiomyocyte", "Pericyte", "Endothelial", "Immune", "Pericyte", 
#                              "Cardiomyocyte", "Immune", "Cardiomyocyte", "Pericyte", "Cardiomyocyte", "Immune-Mast", "Cardiomyocyte", "Neuronal", "Unidentified", "Endothelial")

annotated_clusters_list <- c(
  "Cardiomyocyte", #0 
  "Endothelial",   #1 
  "Fibroblast",    #2 
  "Fibroblast",    #3 
  "Pericyte",      #4 
  "Endothelial",   #5 
  "Cardiomyocyte", #6 
  "Endothelial",   #7 
  "Endothelial",   #8 
  "Pericyte",      #9 
  "Cardiomyocyte", #10 
  "Pericyte",      #11 
  "Fibroblast",    #12 
  "Endothelial",   #13 
  "Pericyte",      #14 
  "Cardiomyocyte", #15 
  "Cardiomyocyte", #16 
  "Immune",        #17 
  "Pericyte",      #18 
  "Endothelial",   #19 
  "Cardiomyocyte", #20 
  "Cardiomyocyte", #21 
  "Neuronal",      #22 
  "Cardiomyocyte", #23 
  "Immune-Mast",   #24 
  "Cardiomyocyte", #25 
  "Fibroblast",    #26 
  "Unidentified",  #27 
  "Neuronal",      #28 
  "Endothelial"    #29 
)
names(annotated_clusters_list) <- paste0("cluster_", 0:(length(annotated_clusters_list) - 1))
annotated_clusters <- sapply(control.disease.combined$RNA_snn_res.1, function(x) annotated_clusters_list[paste0("cluster_", x)])
control.disease.combined@meta.data <- cbind(control.disease.combined@meta.data, annotated_clusters)
##### Reorder the cell types in the metadata
cell_type_order <- c("Cardiomyocyte", "Fibroblast", "Endothelial", "Pericyte", "Neuronal", "Immune", "Immune-Mast", "Unidentified")
control.disease.combined$annotated_clusters <- factor(control.disease.combined$annotated_clusters, levels = cell_type_order)

dir.create(paste0(result_dir, "/Annotation/Annotated_clustering"))
pdf(paste0(result_dir, "/Annotation/Annotated_clustering/Annotated_clustering.pdf"), width = 8, height = 6)
  DimPlot(control.disease.combined, reduction = 'umap', group.by = "annotated_clusters", label = T, label.size = 5, repel = T, cols = color_list_02, order = T) + ggtitle(NULL) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 20))
dev.off()

pdf(paste0(result_dir, "/Annotation/Annotated_clustering/Annotated_clustering_split_condition.pdf"), width = 13, height = 5.5)
  DimPlot(control.disease.combined, reduction = 'umap', group.by = "annotated_clusters", split.by = "condition", label = T, label.size = 5, repel = T, cols = color_list_02, order = T) + ggtitle(NULL) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                             panel.border = element_rect(size = 0.5), text = element_text(size = 24))
dev.off()

pdf(paste0(result_dir, "/Annotation/Annotated_clustering/Annotated_clustering_split_sample.pdf"), width = 50, height = 4)
  DimPlot(control.disease.combined, reduction = 'umap', group.by = "annotated_clusters", split.by = "samples", label = T, label.size = 3, cols = color_list_02, order = T) + 
    facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x") + ggtitle(NULL) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

pdf(paste0(result_dir, "/Annotation/Annotated_clustering/Annotated_clustering_split_sample_donor.pdf"), width = 28, height = 4)
DimPlot(control.disease.combined, reduction = 'umap', group.by = "annotated_clusters", split.by = "donor_condition", label = T, label.size = 3, cols = color_list_02, order = T) + 
  facet_grid(. ~  factor(donor_condition, level = selected_donor_names), scale = "free_x") + ggtitle(NULL) + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
dev.off()

dir.create(paste0(result_dir, "/Annotation/Annotated_dotplot"))
control.disease.combined_Control <- subset(control.disease.combined, subset = condition == "Control")
control.disease.combined_Disease <- subset(control.disease.combined, subset = condition == "IHD")

rm(dotplot_all)
dotplot_all <- DotPlot(control.disease.combined, features = all_markers, group.by = "annotated_clusters")
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all <- dotplot_all + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + 
  theme_linedraw() + RotatedAxis() + scale_color_gradient(low = "white", high = "black") + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(paste0(result_dir, "/Annotation/Annotated_dotplot/Annotated_dotplot_all.pdf"), plot = dotplot_all, width = 17, height = 3.5, dpi = 300)

dotplot_Control <- DotPlot(control.disease.combined_Control, features = all_markers, group.by = "annotated_clusters")
dotplot_Control$data$gene_group <- factor(gene_group_df$group[match(dotplot_Control$data$features.plot, gene_group_df$gene)])
dotplot_Control$data$condition <- "Control"
dotplot_Control <- dotplot_Control + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + 
  theme_linedraw() + RotatedAxis() + scale_color_gradient(low = "white", high = "dodgerblue4") + guides(colour = guide_colourbar()) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
ggsave(paste0(result_dir, "/Annotation/Annotated_dotplot/Annotated_dotplot_Control.pdf"), plot = dotplot_Control, width = 14, height = 4, dpi = 300)

dotplot_Disease <- DotPlot(control.disease.combined_Disease, features = all_markers, group.by = "annotated_clusters")
dotplot_Disease$data$gene_group <- factor(gene_group_df$group[match(dotplot_Disease$data$features.plot, gene_group_df$gene)])
dotplot_Disease$data$condition <- "IHD"
dotplot_Disease <- dotplot_Disease + facet_grid(~ gene_group, scales = "free_x", space = "free_x") + 
  theme_linedraw() + RotatedAxis() + scale_color_gradient(low = "white", high = "firebrick4") + guides(colour = guide_colourbar()) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
ggsave(paste0(result_dir, "/Annotation/Annotated_dotplot/Annotated_dotplot_Disease.pdf"), plot = dotplot_Disease, width = 14, height = 4, dpi = 300)

################################################################################################################################################################
### statistic for unannotated and annotated clustering
################################################################################################################################################################
num_per_donor_cluster_annotated <- table(control.disease.combined$donor_condition, control.disease.combined$annotated_clusters)
num_per_donor_cluster_annotated <- num_per_donor_cluster_annotated[unique(control.disease.combined$donor_condition), ]

num_per_sample_cluster_annotated <- table(control.disease.combined$samples, control.disease.combined$annotated_clusters)
num_per_sample_cluster_annotated <- num_per_sample_cluster_annotated[unique(control.disease.combined$samples), ]

percentage_per_donor_cluster_annotated <- num_per_donor_cluster_annotated %>% 
  as.data.frame.matrix() %>% mutate(Total = rowSums(.)) %>% 
  mutate(across(-Total, ~ round(.x / Total *100 , 5))) %>% select(-Total) %>% 
  mutate(condition = c("Control", "Control", "Control", "IHD", "IHD", "IHD", "IHD", "IHD"))
  
percentage_per_sample_cluster_annotated <- num_per_sample_cluster_annotated %>% 
  as.data.frame.matrix() %>% mutate(Total = rowSums(.)) %>% 
  mutate(across(-Total, ~ round(.x / Total *100 , 5))) %>% select(-Total) %>% 
  # mutate(condition = c("Control", "Control", "Control", "IHD", "IHD", "IHD", "IHD", "IHD"))
  mutate(condition = c("Control", "Control", "Control", "Control", "Control", "Control", "Control", 
                       "IHD", "IHD", "IHD", "IHD", "IHD", "IHD", "IHD"))

percentage_per_sample_cluster_annotated_melt <- percentage_per_sample_cluster_annotated %>% melt() %>% 
  setNames(c("condition", "cluster", "ratio")) %>% mutate(condition = factor(condition, level = c("Control", "IHD"))) %>% mutate(ratio = as.numeric(ratio))

pdf(paste0(result_dir, "/Annotation/Annotated_dotplot/Annotated_cluster_prop_comparison.pdf"), width = 20, height = 6)
  barplot_annotated_cluster_ratio <- ggbarplot(percentage_per_sample_cluster_annotated_melt, x = "cluster", y = "ratio", color = "condition", 
                                              label = TRUE, lab.nb.digits = 2, add = c("mean_se", "jitter"), position = position_dodge(0.9), palette = c("dodgerblue4", "firebrick4")) +
    stat_compare_means(aes(group = condition), method = "wilcox", label = "p.format", label.y = 1.02 * max(percentage_per_sample_cluster_annotated_melt$ratio, na.rm = TRUE)) + 
    stat_compare_means(aes(group = condition), method = "wilcox", label = "p.signif", label.y = 1.06 * max(percentage_per_sample_cluster_annotated_melt$ratio, na.rm = TRUE)) + 
    labs(x = "Unannotated Cluster", y = "proportion of cell") + theme_bw(base_size = 24) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 20)) + 
    scale_fill_manual(values = c("dodgerblue4", "firebrick4"))
  print(barplot_annotated_cluster_ratio)
dev.off()

## Cell Type Composition
melt_num_per_donor_cluster_annotated <- melt(num_per_donor_cluster_annotated) %>% mutate(value = value + 1) %>% 
  mutate(condition = ifelse(Var1 %in% c("936_Ctrl", "5087_Ctrl", "5919_Ctrl"), "Control", 
                            ifelse(Var1 %in% c("604_IHD", "1156_IHD", "5111_IHD", "5364_IHD", "5874_IHD"), "IHD", value)))
melt_percentage_per_donor_cluster_annotated <- melt(as.matrix(percentage_per_donor_cluster_annotated[,-9])) %>% 
  mutate(condition = ifelse(Var1 %in% c("936_Ctrl", "5087_Ctrl", "5919_Ctrl"), "Control", 
                            ifelse(Var1 %in% c("604_IHD", "1156_IHD", "5111_IHD", "5364_IHD", "5874_IHD"), "IHD", value)))

red_palette <- c("grey", "#FFCCCC", "#FFB2B2", "#FF9999", "#FF6666", "#FF4D4D", "#FF3333","#FF1A1A", "#FF0000", "#CC0000",
                 "#B30000", "#990000", "#800000", "#660000", "#4D0000", "#330000")

pdf(paste0(result_dir, "/Annotation/Annotated_dotplot/cell_count_in_annotated_cluster.pdf"), width = 10, height = 8)
ggplot(melt_num_per_donor_cluster_annotated, aes(x = Var2, y = Var1, fill = log10(value))) + 
  geom_tile() + geom_text(aes(label = round(value - 1, 1)), color = "white", size = 6) + 
  scale_fill_gradientn(colors = red_palette, name = "# of cells", breaks = log10(c(1, 10, 100, 1000)), labels = c(1, 10, 100, 1000)) + 
  facet_grid(condition ~ ., scales = "free", space = "free") + theme_linedraw() + labs(title = "", x = "Cell type", y = "Donor") + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(face = "bold", size = 16), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()

pdf(paste0(result_dir, "/Annotation/Annotated_dotplot/cell_percentage_in_annotated_cluster.pdf"), width = 10, height = 8)
ggplot(melt_percentage_per_donor_cluster_annotated, aes(x = Var2, y = Var1, fill = value)) + 
  geom_tile() + geom_text(aes(label = round(value, 1)), color = "white", size = 6) + 
  # scale_fill_gradientn(colors = red_palette, name = "# of cells", breaks = log10(c(1, 10, 100, 1000)), labels = c(1, 10, 100, 1000)) + 
  scale_fill_gradientn(colors = red_palette, name = "% of cells") + 
  facet_grid(condition ~ ., scales = "free", space = "free") + theme_linedraw() + labs(title = "", x = "Cell type", y = "Donor") + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 16), panel.border = element_rect(size = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

##### save object
saveRDS(control.disease.combined, paste0(result_dir, "/Seurat.obj_with_annotation_all_celltypes.RDS"))


control.disease.combined  <- readRDS(paste0(result_dir, "/Seurat.obj_with_annotation_all_celltypes.RDS"))
Idents(object = control.disease.combined) <- "annotated_clusters"
control.disease.combined <- subset(control.disease.combined, idents = c("Cardiomyocyte", "Fibroblast"))

cardiomyocyte_marker_list <- c("TNNT2", "ACTN2", "MYOZ2", "MYPN", "MLIP")
# fibroblast_marker_list <- c("DCN", "MMP2", "LUM")
fibroblast_marker_list <- c("DCN", "APOD", "BICC1", "ABCA6")
# fibroblast_marker_list <- c("DCN", "MMP2")
# hypoxia_marker_list <- c("JUN", "NR4A1", "ZEB1", "ZEB2")
hypoxia_marker_list <- c("JUN", "NR4A1")
# collagen_marker_list <- c("COL4A2", "COL6A1", "COL6A2", "COL21A1")
collagen_marker_list <- c("COL21A1", "COL6A2")
# fibrotic_marker_list <- c("ACTA2", "TGFB1", "NFKB1", "SMAD1", "SMAD5", "SMAD9", "SMAD6", "SMAD3", "SMAD2", "SMAD7", "SMAD4", "TNF", "CCN2", "ECM1", "ECM2", "DCN", "VIM", "PDGFRA")
# fibrotic_marker_list <- c("COL3A1", "S100A4", "PDGFRB")
fibrotic_marker_list <- c("SMAD3")
# MMR_pathway_marker_list <- c("MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2")
MMR_pathway_marker_list <- c("MLH1")
BER_pathway_marker_list <- c("APEX1", "FEN1", "LIG1", "LIG3", "MBD4", "NEIL1", "NEIL12", "OGG1", "POLB", "PNKP")
# BER_pathway_marker_list <- c("APEX1")
# NHEJ_pathway_marker_list <- c("XRCC5", "XRCC6", "PRKDC", "LIG4", "DCLRE1C", "NHEJ1", "XRCC4", "PAXX", "MRE11", "TP53BP1", "SHLD1", "SHLD2", "SHLD3")
NHEJ_pathway_marker_list <- c("SHLD1")
# ddr_genes <- c("ATM", "ATR", "BRCA1", "RAD50", "RAD52", "RAD54", "TP53", "XRCC5", "XRCC6")

##### check the marker expression in cluster
gene_list <- rownames(control.disease.combined)
MMR_pathway_marker_list %in% gene_list
list_to_check <- collagen_marker_list
Idents(object = control.disease.combined) <- ""
# VlnPlot(control.disease.combined, group.by = "annotated_clusters", features = list_to_check, split.by = "condition", ncol = 5)
VlnPlot(control.disease.combined, features = list_to_check, split.by = "condition", ncol = 5)

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
  average_expr <- colMeans(as.data.frame(GetAssayData(control.disease.combined[all_marker_genes[[name]], ])))
  control.disease.combined@meta.data[, name] <- average_expr
}

metadata_CM_Fib_02 <- control.disease.combined@meta.data
# head(control.disease.combined)
# all_metagene_grouped_avg_02 <- metadata_CM_Fib_02 %>% group_by(annotated_clusters, condition, samples) %>%
all_metagene_grouped_avg_02 <- metadata_CM_Fib_02 %>% group_by(condition, samples) %>%
  summarize(cardiomyocyte_metagene = mean(metagene_cardiomyocyte, na.rm = TRUE), fibroblast_metagene = mean(metagene_fibroblast, na.rm = TRUE), 
            hypoxia_metagene = mean(metagene_hypoxia, na.rm = TRUE), collagen_metagene = mean(metagene_collagen, na.rm = TRUE), 
            fibrotic_metagene = mean(metagene_fibrotic, na.rm = TRUE), MMR_metagene = mean(metagene_MMR, na.rm = TRUE), 
            BER_metagene = mean(metagene_BER, na.rm = TRUE), NHEJ_metagene = mean(metagene_NHEJ, na.rm = TRUE))

all_metagene_grouped_avg_melt_02 <- all_metagene_grouped_avg_02 %>% 
  # filter(RNA_snn_res.1 %in% c(9, 11)) %>%
  # melt(id.vars = c("annotated_clusters", "condition", "samples"),
  melt(id.vars = c("condition", "samples"),
       measure.vars = c("cardiomyocyte_metagene", "fibroblast_metagene", "hypoxia_metagene", "collagen_metagene", 
                        "fibrotic_metagene", "MMR_metagene", "BER_metagene", "NHEJ_metagene"),
       variable.name = "metagenes", value.name = "avg_expr")

pdf(paste0(result_dir, "/Annotation_extract_CM_Fib/metagene_expr_barplot_all_celltype.pdf"), width = 18, height = 8)
ggplot(all_metagene_grouped_avg_melt_02, aes(x = metagenes, y = avg_expr, color = condition)) + 
  geom_boxplot(aes(fill = condition), position = position_dodge(), alpha = 0.5, outlier.shape = NA) + 
  geom_jitter(aes(color = condition), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 2) + 
  stat_compare_means(aes(group = condition), label = "p.format", label.y = 1.02 * max(all_metagene_grouped_avg_melt_02$avg_expr, na.rm = TRUE)) + 
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.06 * max(all_metagene_grouped_avg_melt_02$avg_expr, na.rm = TRUE)) + 
  scale_fill_manual(values = Ctrl_IHD_color) + scale_color_manual(values = Ctrl_IHD_color) + 
  # facet_wrap(. ~  annotated_clusters) +
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           text = element_text(face = "bold", size = 18), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + 
  labs(x = "", y = "Average metagene expression")
dev.off()

