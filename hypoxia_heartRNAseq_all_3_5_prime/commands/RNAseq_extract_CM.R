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

# color_list <- c("#D18E68","#F22020","#FCFF5D","#7DFC00","#0EC434","#228C68","#8AD8E8","#235B54","#29BDAB","#3998F5",
#                 "#37294F","#277DA7","#3750DB","#991919","#FFCBA5","#E68F66","#C56133","#96341C","#632819",
#                 "#FFC413","#F47A22","#2F2AA0","#B732CC","#772B9D","#F07CAB","#D30B94","#EDEFF3","#C3A5B4","#946AA2","#5D4C86")
# color_list <- c("#F22020","#3998F5","#FCFF5D","#7DFC00","#8AD8E8","#235B54","#29BDAB","#D18E68",
#                 "#37294F","#277DA7","#3750DB","#991919","#FFCBA5","#E68F66","#C56133","#96341C","#632819",
#                 "#FFC413","#F47A22","#2F2AA0","#B732CC","#772B9D","#F07CAB","#D30B94","#C3A5B4","#946AA2","#5D4C86","#EDEFF3")

color_list <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                "#393b79", "#5254a3", "#6b6ecf", "#9c9ede", "#637939", "#8ca252", "#b5cf6b", "#cedb9c", "#8c6d31", "#bd9e39",
                "#e7ba52", "#e7969c", "#d6616b", "#ad494a", "#843c39", "#f5f5f5", "#fdd0a2", "#fb6a4a", "#cb181d", "#a50f15",
                "#3182bd", "#6baed6", "#9ecae1", "#c6dbef", "#e6550d", "#fd8d3c", "#fdae6b", "#fdd0a2", "#31a354", "#74c476",
                "#a1d99b", "#c7e9c0", "#756bb1", "#9e9ac8", "#bcbddc")

##### set working and results directories
setwd("/Users/zhemingan/Documents/BCH_research/heart_RNAseq/LOODEG")
wd_dir <- getwd()
result_dir <- paste0(wd_dir, "/five_Control_five_Disease_Results_new")
dir.create(result_dir)
# dir.create(paste0(result_dir, "/Annotation"))
dir.create(paste0(result_dir, "/Annotation/extract_CM"), recursive = T)
dir.create(paste0(result_dir, "/Marker_files"))
dir.create(paste0(result_dir, "/DEGs"))
dir.create(paste0(result_dir, "/DEGs/no_lb"))
dir.create(paste0(result_dir, "/DEGs/with_lb"))
##### read sample metadata
sample_metadata <- read.csv(paste0(wd_dir, "/sample_meta.csv"), header = TRUE)
sample_names <- sample_metadata$sample_name
sample_condition <- sample_metadata$Conditions
both_sample_condition <- paste0(sample_metadata$Conditions, "_", sample_metadata$sample_name)
cellranger_data.dir <- sapply(sample_names, function(x) paste0(wd_dir, "/cellranger_results/",x, "/"))

##### read in all 10X files as matrices 
list.cellranger.result <- list()
for (dir in cellranger_data.dir){
  print(basename(dir))
  list.cellranger.result[[basename(dir)]] <-  Read10X(data.dir = paste0(dir, "filtered_feature_bc_matrix/"))
}

##### create individual Seurat objects from matrices
options(Seurat.object.assay.version = "v3")
list.of.Seurat.objs <- lapply(seq_along(list.cellranger.result), function(x) 
  CreateSeuratObject(counts = list.cellranger.result[[x]], min.cells = 3, 
                     min.features = 200, project = names(list.cellranger.result)[x])) 
names(list.of.Seurat.objs) <- names(list.cellranger.result)
for (obj_index in seq_along(list.of.Seurat.objs)){
  obj <- list.of.Seurat.objs[[obj_index]]
  name <- names(list.of.Seurat.objs)[obj_index]
  obj <- RenameCells(obj, add.cell.id = name)
  list.of.Seurat.objs[[obj_index]] <- obj                                    
}

##### create merged Seurat object
control.disease.combined <- merge(list.of.Seurat.objs[[1]], y = list.of.Seurat.objs[2:10], project = "heart_merged_analysis")
head(control.disease.combined)
selected_sample_names <- unique(control.disease.combined$orig.ident)
# add more info to the merged object
control.disease.combined[["percent_mito"]] <- PercentageFeatureSet(control.disease.combined, pattern = "^MT-")
control.disease.combined[["percent_ribo"]] <- PercentageFeatureSet(control.disease.combined, pattern = "^RPL|^RPS")
## merge 5919_2n and 5919_CM and assign new sample name column
control.disease.combined$samples <- sapply(control.disease.combined$orig.ident, function(x) sample_names[which(sample_names == x)])
control.disease.combined$samples <- ifelse(control.disease.combined$samples %in% c("5919_2n", "5919_CM"), "5919", control.disease.combined$samples)
control.disease.combined$condition <- sapply(control.disease.combined$orig.ident, function(x) sample_condition[which(sample_names == x)])
control.disease.combined$both_sample_condition <- paste(control.disease.combined$samples, control.disease.combined$condition, sep = "_")
# control.disease.combined$both_sample_condition <- sapply(control.disease.combined$samples, function(x) both_sample_condition[which(sample_names == x)])
metadata <- control.disease.combined@meta.data
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)
head(control.disease.combined@meta.data, 5)

##### QUALITY CONTROL AND FILTERING for control.disease.combined
dir.create(paste0(result_dir, "/QC_Control"))
control.disease.combined<- Add_Cell_Complexity(object = control.disease.combined)
all.equal(colnames(control.disease.combined), row.names(control.disease.combined@meta.data))

##### create plots of QC metrics before filtering
pdf(paste0(result_dir, "/QC_Control/QC_pre_filtering.pdf"), width = 14, height = 3.5)
Idents(control.disease.combined) <- "orig.ident"
# Idents(control.disease.combined) <- "samples"
p1 <- QC_Plots_Genes(seurat_object = control.disease.combined, plot_title = "Genes Per Cell")
# Unique Molecular Identifiers, absolute number of observed transcripts
p2 <- QC_Plots_UMIs(seurat_object = control.disease.combined, plot_title = "UMIs Per Cell") 
p3 <- QC_Plots_Mito(seurat_object = control.disease.combined, plot_title = "Mito Gene % Per Cell")
p4 <- QC_Plots_Feature(seurat_object = control.disease.combined, feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell")
wrap_plots(p1, p2, p3, p4, ncol = 4)
dev.off()

################ Cell counts: The cell counts are determined by the number of unique cellular barcodes detected.
pdf(paste0(result_dir, "/QC_Control/nCells_per_sample_condition.pdf"), width = 15, height = 7)
p1 <- metadata %>% 
  ggplot(aes(x=factor(orig.ident, level = selected_sample_names), fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  labs(x = "Sample ID") +
  ggtitle("nCells per sample")
p2 <- metadata %>% 
  ggplot(aes(x=condition, fill=condition)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCells per condition")
wrap_plots(p1, p2, ncol = 2)
dev.off()

################ UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500. 
# If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply.
# Visualize the number UMIs/transcripts per cell
pdf(paste0(result_dir, "/QC_Control/nUMI_per_sample_condition.pdf"), width = 15, height = 7)
p3 <- metadata %>% 
  ggplot(aes(color=factor(orig.ident, level = selected_sample_names), x=nUMI, 
             fill= factor(orig.ident, level = selected_sample_names))) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.title=element_blank()) +
  ggtitle("nUMI per sample")
p4 <- metadata %>% 
  ggplot(aes(color=condition, x=nUMI, fill= condition)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.title=element_blank()) +
  ggtitle("nUMI per condition")
#annotate("text", x=500, y=-0.05, label="500", angle=0)
wrap_plots(p3, p4, ncol = 2)
dev.off()

##### Genes detected per cell
pdf(paste0(result_dir, "/QC_Control/nGene_per_sample_condition.pdf"), width = 15, height = 7)
# Visualize the distribution of genes detected per cell via histogram
p5 <- metadata %>% 
  ggplot(aes(color=factor(orig.ident, level = selected_sample_names), x=nGene, 
             fill= factor(orig.ident, level = selected_sample_names))) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) +
  ggtitle("nGenes per cell for each sample") + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) + 
  theme(legend.title=element_blank())
# Visualize the distribution of genes detected per cell via boxplot
p6 <- metadata %>% 
  ggplot(aes(x=condition, y=log10(nGene), fill=condition)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(legend.title=element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nGenes per cell for each condition")
wrap_plots(p5, p6, ncol = 2)
dev.off()

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata$log10GenesPerUMI <- control.disease.combined$log10GenesPerUMI
pdf(paste0(result_dir, "/QC_Control/Genes_per_UMI.pdf"), width = 15, height = 7)
p7 <- metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = factor(orig.ident, level = selected_sample_names), 
             fill=factor(orig.ident, level = selected_sample_names))) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  theme(legend.title=element_blank()) +
  geom_vline(xintercept = 0.8)
print(p7)
p8 <- metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = condition, fill=condition)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  theme(legend.title=element_blank()) +
  geom_vline(xintercept = 0.8)
wrap_plots(p7, p8, ncol = 2)
dev.off()

# metadata %>% 
#   ggplot(aes(color=samples, x=percent.mt, fill=samples)) + 
#   geom_density(alpha = 0.2) + 
#   scale_x_log10() + 
#   theme_classic() +
#   geom_vline(xintercept = 0.2)

#determine cells that will be lost by each QC filtering step
nFeature_RNA_five_q <- quantile(control.disease.combined$nFeature_RNA, 0.05)
nFeature_RNA_ninty_five_q <- quantile(control.disease.combined$nFeature_RNA, 0.95)
nCount_RNA_five_q <- quantile(control.disease.combined$nCount_RNA, 0.05)
nCount_RNA_ninty_five_q <- quantile(control.disease.combined$nCount_RNA, 0.95)

cells_filtered_pct_mito <- row.names(control.disease.combined@meta.data %>% filter(percent_mito > 5))
cells_filtered_pct_ribo <- row.names(control.disease.combined@meta.data %>% filter(percent_ribo > 5))
cells_filtered_nFeatures <- row.names(control.disease.combined@meta.data %>% filter(nFeature_RNA < nFeature_RNA_five_q | nFeature_RNA > nFeature_RNA_ninty_five_q))
cells_filtered_nCounts <- row.names(control.disease.combined@meta.data %>% filter(nCount_RNA < nCount_RNA_five_q | nCount_RNA > nCount_RNA_ninty_five_q))
cells_filtered_log10GenesPerUMI <- row.names(control.disease.combined@meta.data %>% filter(log10GenesPerUMI < 0.80))
filtering_list <- list("pct_mito" = cells_filtered_pct_mito, "pct_ribo" = cells_filtered_pct_ribo, 
                       "nFeatures" = cells_filtered_nFeatures, "nCounts" = cells_filtered_nCounts,
                       "log10GenesPerUMI" = cells_filtered_log10GenesPerUMI)
print("NUMBER OF CELLS FILTERED PER CRITERIA: ")
print(lapply(filtering_list, length))
print(length(unique(do.call(c, filtering_list))))
print(length(unique(do.call(c, filtering_list)))/ncol(control.disease.combined))

pdf(paste0(result_dir, "/QC_Control/cells_filtered_upset_plot.pdf"))
upset(fromList(filtering_list))
dev.off()

#filter cells 
control.disease.combined<- subset(control.disease.combined, subset = percent_mito <= 5)
control.disease.combined <- subset(control.disease.combined, subset = percent_ribo <= 5)
control.disease.combined <- subset(control.disease.combined, subset = nFeature_RNA >= nFeature_RNA_five_q & nFeature_RNA <= nFeature_RNA_ninty_five_q)
control.disease.combined <- subset(control.disease.combined, subset = nCount_RNA >= nCount_RNA_five_q & nCount_RNA <= nCount_RNA_ninty_five_q)
control.disease.combined <- subset(control.disease.combined, subset = log10GenesPerUMI >= 0.80)

#create plots of qc metrics after filtering
pdf(paste0(result_dir, "/QC_Control/QC_post_filtering.pdf"), width = 14, height = 3.5)
Idents(control.disease.combined) <- "orig.ident"
p11 <- QC_Plots_Genes(seurat_object = control.disease.combined, plot_title = "Genes Per Cell")
p12 <- QC_Plots_UMIs(seurat_object = control.disease.combined, plot_title = "UMIs Per Cell")
p13 <- QC_Plots_Mito(seurat_object = control.disease.combined, plot_title = "Mito Gene % Per Cell")
p14 <- QC_Plots_Feature(seurat_object = control.disease.combined, feature = "percent_ribo", plot_title = "Ribo Gene % Per Cell")
wrap_plots(p11, p12, p13, p14, ncol = 4)
dev.off()
saveRDS(control.disease.combined, paste0(result_dir, "/merged.Seurat.obj_filtered.RDS"))

##########################################################################
#NORMALIZATION, scaling with scTransform,
#DIMENSIONALITY REDUCTION, AND CLUSTER DATA
control.disease.combined <- readRDS(paste0(result_dir, "/merged.Seurat.obj_filtered.RDS"))

DefaultAssay(control.disease.combined) <- 'RNA'
control.disease.combined <- NormalizeData(control.disease.combined)
control.disease.combined <- FindVariableFeatures(control.disease.combined, selection.method = "vst", nfeatures = 3000)
control.disease.combined <- ScaleData(control.disease.combined, vars.to.regress = c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA"))

# check_PC <- RunPCA(control.disease.combined, verbose = TRUE)
# ElbowPlot(check_PC, ndims = 50, reduction = "pca")

control.disease.combined <- RunPCA(control.disease.combined, npcs = 30, verbose = TRUE) %>% 
  RunHarmony(group.by.vars = "samples", lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  RunTSNE(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = TRUE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = TRUE) %>%
  FindClusters(resolution = c(0.8, 1.0, 1.2), verbose = TRUE)
# FindClusters(resolution = c(0.4, 0.6, 1.0, 1.4, 1.6), verbose = TRUE)

pdf(paste0(result_dir, "/Annotation/Dim_Plot_before_doublet.pdf"), width = 7, height = 5)
ElbowPlot(control.disease.combined, ndims = 30, reduction = "pca")
DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE, group.by = 'samples') + ggtitle("PCA")
DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE, group.by = 'condition') + ggtitle("PCA")
DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE, group.by = 'samples') + ggtitle("TSNE")
DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE, group.by = 'condition') + ggtitle("TSNE")
DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'samples') + ggtitle("UMAP")
DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'condition') + ggtitle("UMAP")
DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'samples') + ggtitle("UMAP")
DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'condition') + ggtitle("UMAP")

DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE) + ggtitle("PCA")
DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE) + ggtitle("TSNE")
DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE) + ggtitle("UMAP")
DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE) + ggtitle("UMAP")
dev.off()

saveRDS(control.disease.combined, paste0(result_dir, "/combined.Seurat.obj_for_annotation.RDS"))
control.disease.combined <- readRDS(paste0(result_dir, "/combined.Seurat.obj_for_annotation.RDS"))

################# Doublet finder
sweep.res.list_combined <- paramSweep_v3(control.disease.combined, PCs = 1:30, sct = FALSE)
sweep.stats_combined <- summarizeSweep(sweep.res.list_combined, GT = FALSE)
bcmvn_combined <- find.pK(sweep.stats_combined)

# Plot the bcmvn output to find the peak
pdf(paste0(result_dir, "/QC_Control/Optimal_pK_doubletFinder.pdf"), width = 12, height = 5)
ggplot(bcmvn_combined, aes(x = pK, y = BCmetric)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Optimal pK Value", x = "pK", y = "BCmvn")
dev.off()

homotypic.prop <- modelHomotypic(control.disease.combined@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.060*nrow(control.disease.combined@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

control.disease.combined <- doubletFinder_v3(control.disease.combined, PCs = 1:30, pK = 0.060, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(control.disease.combined)

pdf(paste0(result_dir, "/QC_Control/DF_classifications.pdf"), width = 7, height = 5)
DimPlot(control.disease.combined, reduction = "umap", group.by = "DF.classifications_0.25_0.06_2823")
dev.off()

table(control.disease.combined@meta.data$DF.classifications_0.25_0.06_2823)
# Doublet Singlet 
# 2823   44232

control.disease.combined <- subset(control.disease.combined, subset = DF.classifications_0.25_0.06_2823 == "Singlet")
control.disease.combined <- NormalizeData(control.disease.combined)
control.disease.combined <- FindVariableFeatures(control.disease.combined, selection.method = "vst", nfeatures = 3000)
control.disease.combined <- ScaleData(control.disease.combined, vars.to.regress = c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA"))
control.disease.combined <- RunPCA(control.disease.combined, npcs = 30, verbose = TRUE) %>% 
  RunHarmony(group.by.vars = "samples", lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  RunTSNE(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = TRUE) %>% 
  FindClusters(resolution = c(0.8, 1.0, 1.2), verbose = TRUE)
# FindClusters(resolution = c(0.4, 0.6, 1.0, 1.4, 1.6), verbose = TRUE)

pdf(paste0(result_dir, "/Annotation/Dim_Plot_after_doublet.pdf"), width = 7, height = 5)
ElbowPlot(control.disease.combined, ndims = 30, reduction = "pca")
DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE, group.by = 'samples') + ggtitle("PCA")
DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE, group.by = 'condition') + ggtitle("PCA")
DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE, group.by = 'samples') + ggtitle("TSNE")
DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE, group.by = 'condition') + ggtitle("TSNE")
DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'samples') + ggtitle("UMAP")
DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE, group.by = 'condition') + ggtitle("UMAP")
DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'samples') + ggtitle("UMAP")
DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE, group.by = 'condition') + ggtitle("UMAP")

DimPlot(object = control.disease.combined, reduction = 'pca', label = TRUE) + ggtitle("PCA")
DimPlot(object = control.disease.combined, reduction = 'tsne', label = TRUE) + ggtitle("TSNE")
DimPlot(object = control.disease.combined, reduction = 'umap', label = TRUE) + ggtitle("UMAP")
DimPlot(object = control.disease.combined, reduction = 'harmony', label = TRUE) + ggtitle("harmony")
dev.off()

saveRDS(control.disease.combined, paste0(result_dir, "/combined.Seurat.obj_remove_doublet.RDS"))

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
### continue from O2
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

##### clustering results summary
control.disease.combined  <- readRDS(paste0(result_dir, "/combined.Seurat.obj_remove_doublet.RDS"))
selected_sample_names <- unique(control.disease.combined$samples)
RNA_resolution <- "RNA_snn_res.1.2"
Idents(object = control.disease.combined) <- RNA_resolution

pdf(paste0(result_dir, "/Annotation/Unannotated_clustering.pdf"), width = 6, height = 5)
DimPlot(control.disease.combined, reduction = 'umap', group.by = RNA_resolution, label = TRUE, label.size = 3, cols = color_list, order = T)
DimPlot(control.disease.combined, reduction = 'tsne', group.by = RNA_resolution, label = TRUE, label.size = 3, cols = color_list, order = T)
dev.off()

pdf(paste0(result_dir, "/Annotation/Unannotated_clustering_split_sample.pdf"), width = 24, height = 5)
DimPlot(control.disease.combined, reduction = 'umap', split.by = "samples", label = TRUE, label.size = 3, cols = color_list, order = T) + 
  facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x")
DimPlot(control.disease.combined, reduction = 'tsne', split.by = "samples", label = TRUE, label.size = 3, cols = color_list, order = T) + 
  facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x")
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
  "Fibroblasts",    #0 fi
  "Cardiomyocytes", #1 ca
  "Endothelial",    #2 en
  "Fibroblasts",    #3 fi
  "Cardiomyocytes", #4 ca
  "Endothelial",    #5 en
  "Pericytes",      #6 pe
  "Pericytes",      #7 pe
  "Endothelial",    #8 en
  "Endothelial",    #9 en
  "Endothelial",    #10 en
  "Pericytes",      #11 pe
  "Cardiomyocytes", #12 ca
  "Neurons",        #13 neuro
  "Myeloid",        #14 Myeloid
  "NK-Cells",       #15 NK-Cells
  "Cardiomyocytes", #16 ca
  "Cardiomyocytes", #17 ca
  "Pericytes",      #18 pe
  "Cardiomyocytes", #19 ca
  "Fibroblasts",    #20 fi
  "Cardiomyocytes", #21 ca
  "Endothelial",    #22 en
  "Cardiomyocytes", #23 ca
  "Mast",           #24 mast
  "Cardiomyocytes", #25 ca
  "Pericytes",      #26 pe
  "Cardiomyocytes", #27 ca
  "Endothelial",    #28 en
  "Adipocyte"       #29 Adipocyte
)

names(annotated_clusters_list) <- paste0("cluster_", 0:(length(annotated_clusters_list) - 1))
annotated_clusters <- sapply(control.disease.combined$RNA_snn_res.1.2, function(x) annotated_clusters_list[paste0("cluster_", x)])
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
        label = T, label.size = 3, cols = color_list, order = T) + facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x")
DimPlot(control.disease.combined, reduction = 'tsne', group.by = "annotated_clusters", split.by = "samples",
        label = T, label.size = 3, cols = color_list, order = T) + facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x")
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
temp1 <- table(control.disease.combined$RNA_snn_res.1.2)
temp_sample <- table(control.disease.combined$samples)
temp_sample <- temp_sample[unique(control.disease.combined$samples)]
temp2 <- t(table(control.disease.combined$RNA_snn_res.1.2, control.disease.combined$condition))
temp3 <- t(table(control.disease.combined$RNA_snn_res.1.2, control.disease.combined$samples))
temp3 <- temp3[unique(control.disease.combined$samples), ]

ratio_by_cluster_by_sample <- as.data.frame.matrix(temp3)
ratio_by_cluster_by_sample <- ratio_by_cluster_by_sample / temp1
ratio_by_cluster_by_sample <- as.data.frame.matrix(ratio_by_cluster_by_sample)
ratio_by_cluster_by_sample$condition <- c("Control", "Control", "Control", "Control", "Disease", "Disease", "Disease", "Disease", "Disease")
ratio_by_cluster_by_sample_melt <- ratio_by_cluster_by_sample %>% 
  melt() %>% 
  setNames(c("condition", "cluster", "ratio")) %>% 
  mutate(condition = factor(condition, level = c("Control", "Disease")))
# table(control.disease.combined_CM$RNA_snn_res.2.6, control.disease.combined_CM$samples)

pdf(paste0(result_dir, "/Annotation/Unannotated_cluster_ratio_comparison.pdf"), width = 26, height = 6)
barplot_strand_ratio_age_match <- ggbarplot(ratio_by_cluster_by_sample_melt, x = "cluster", y = "ratio", color = "condition", 
                                            label = TRUE, lab.nb.digits = 2, add = c("mean_se", "jitter"), position = position_dodge(0.9), 
                                            palette = c("dodgerblue3", "#FC4E07")) +
  stat_compare_means(aes(group = condition), label = "p.format", label.y = 1.02 * max(ratio_by_cluster_by_sample_melt$ratio, na.rm = TRUE)) + 
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.06 * max(ratio_by_cluster_by_sample_melt$ratio, na.rm = TRUE)) + 
  labs(x = "Unannotated Cluster", y = "ratio of cell #") + theme_classic(base_size = 24) + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values=c("dodgerblue3", "#FC4E07"))
print(barplot_strand_ratio_age_match)
dev.off()


temp12 <- table(control.disease.combined$annotated_clusters)
temp22 <- t(table(control.disease.combined$annotated_clusters, control.disease.combined$condition))
temp32 <- t(table(control.disease.combined$annotated_clusters, control.disease.combined$samples))
temp32 <- temp32[unique(control.disease.combined$samples), ]



#save object
control.disease.combined$celltype.condition <- paste(control.disease.combined$annotated_clusters, control.disease.combined$condition, sep = "_")
control.disease.combined$samples.condition <- paste0(control.disease.combined$samples, "-", control.disease.combined$condition)
saveRDS(control.disease.combined, paste0(result_dir, "/Seurat.obj_with_annotation.RDS"))

################################################################################
################################################################################
## extract CM and Fib
################################################################################
################################################################################
control.disease.combined <- readRDS(paste0(result_dir, "/Seurat.obj_with_annotation.RDS"))
all_cell_type <- unique(control.disease.combined$annotated_clusters)
control.disease.combined_CM <- subset(control.disease.combined, subset = annotated_clusters %in% c("Cardiomyocytes", "Fibroblasts"))
# control.disease.combined_CM_test <- RunPCA(control.disease.combined_CM)
# ElbowPlot(control.disease.combined_CM_test, ndims = 30, reduction = "pca")

DefaultAssay(control.disease.combined_CM) <- 'RNA'
control.disease.combined_CM <- RunPCA(control.disease.combined_CM, npcs = 15, verbose = TRUE) %>% 
  RunHarmony(group.by.vars = "samples", lambda = 0.1, plot_convergence = TRUE, reduction.save = "harmony") %>% 
  RunTSNE(reduction = "harmony", dims = 1:15, verbose = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:15, verbose = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:15, verbose = TRUE) %>% 
  FindClusters(resolution = c(0.7, 0.8, 0.9, 1.0, 1.5, 2.4), verbose = TRUE)
  # FindClusters(resolution = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6), verbose = TRUE)

saveRDS(control.disease.combined_CM, paste0(result_dir, "/Seurat.obj_with_annotation_extractCM_Fibro.RDS"))
control.disease.combined_CM  <- readRDS(paste0(result_dir, "/Seurat.obj_with_annotation_extractCM_Fibro.RDS"))

selected_sample_names <- unique(control.disease.combined_CM$samples)
RNA_resolution <- "RNA_snn_res.0.8"
Idents(object = control.disease.combined_CM) <- RNA_resolution
unique(control.disease.combined_CM$RNA_snn_res.0.8)
table_stat <- t(table(control.disease.combined_CM$RNA_snn_res.0.8, control.disease.combined_CM$samples))
table_stat <- table_stat[unique(control.disease.combined_CM$samples), ]
level_list <- c(0:16)
control.disease.combined_CM$RNA_snn_res.0.8 <- factor(control.disease.combined_CM$RNA_snn_res.0.8, levels = level_list)

pdf(paste0(result_dir, "/Annotation/extract_CM/TSNE_UMAP_split_condition.pdf"), width = 7, height = 5)
DimPlot(control.disease.combined_CM, reduction = "tsne", group.by = "condition", label = TRUE, label.size = 4, cols = color_list) +
  theme(text = element_text(size = 6), axis.text = element_text(size = 6)) + ggtitle("TSNE")
DimPlot(control.disease.combined_CM, reduction = "tsne", group.by = RNA_resolution, label = TRUE, label.size = 4, cols = color_list) +
  theme(text = element_text(size = 6), axis.text = element_text(size = 6)) + ggtitle("TSNE")

DimPlot(control.disease.combined_CM, reduction = "umap", group.by = "condition", label.size = 4, cols = color_list) +
  theme(text = element_text(size = 6), axis.text = element_text(size = 6)) + ggtitle("UMAP")
DimPlot(control.disease.combined_CM, reduction = "umap", group.by = RNA_resolution, label = TRUE, label.size = 4, cols = color_list) +
  theme(text = element_text(size = 6), axis.text = element_text(size = 6)) + ggtitle("UMAP")
  # FeaturePlot(control.disease.combined_CM, features = c("DCN"), order = T)
dev.off()

temp13 <- table(control.disease.combined_CM$RNA_snn_res.0.8)
temp23 <- t(table(control.disease.combined_CM$RNA_snn_res.0.8, control.disease.combined_CM$condition))
temp33 <- t(table(control.disease.combined_CM$RNA_snn_res.0.8, control.disease.combined_CM$samples))
temp33 <- temp33[unique(control.disease.combined_CM$samples), ]

pdf(paste0(result_dir, "/Annotation/extract_CM/Unannotated_clustering_split_sample.pdf"), width = 24, height = 5)
DimPlot(control.disease.combined_CM, reduction = "umap", group.by = RNA_resolution, label = TRUE, split.by = "samples", cols = color_list) + 
  facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x")
DimPlot(control.disease.combined_CM, reduction = "tsne", group.by = RNA_resolution, label = TRUE, split.by = "samples", cols = color_list) + 
  facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x")
dev.off()

pdf(paste0(result_dir, "/Annotation/extract_CM/Unannotated_clustering_split_condition.pdf"), width = 8, height = 4)
DimPlot(control.disease.combined_CM, reduction = "umap", group.by = RNA_resolution, label = TRUE, split.by = "condition", label.size = 4, cols = color_list)
DimPlot(control.disease.combined_CM, reduction = "tsne", group.by = RNA_resolution, label = TRUE, split.by = "condition", label.size = 4, cols = color_list)
dev.off()

################################################################################
################################################################################
### dotplot for marker genes
marker_set1 <- c("TNNT2", "ACTN2", "MYOZ2", "RYR2", "MYPN", "MLIP")
marker_set3 <- c("DCN", "APOD", "GSN", "BICC1", "ABCA6", "CDH19")
all_markers <- c(marker_set1, marker_set3)
marker_groups <- c(rep("CM markers", length(marker_set1)), rep("Fibroblasts markers", length(marker_set3)))

gene_group_df <- data.frame(gene = all_markers, group = marker_groups) %>% mutate(group = factor(group, level = c("CM markers", "Fibroblasts markers")))

pdf(paste0(result_dir, "/Annotation/extract_CM/Unannotated_dotplot.pdf"), width = 15, height = 10)
# DotPlot(control.disease.combined, features = feature_gene_list) + RotatedAxis()
dotplot_all <- DotPlot(control.disease.combined_CM, features = all_markers, group.by = RNA_resolution) + RotatedAxis() + scale_color_gradient(low = "lightblue", high = "darkblue")
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all <- dotplot_all + facet_grid(~ gene_group, scales = "free_x", space = "free_x")
dotplot_all

control.disease.combined_CM_Control <- subset(control.disease.combined_CM, subset = condition == "Control")
dotplot_Control <- DotPlot(control.disease.combined_CM_Control, features = all_markers, group.by = c(RNA_resolution)) + RotatedAxis() + scale_color_gradient(low = "lightblue", high = "darkblue")
dotplot_Control$data$gene_group <- factor(gene_group_df$group[match(dotplot_Control$data$features.plot, gene_group_df$gene)])
dotplot_Control <- dotplot_Control + facet_grid(~ gene_group, scales = "free_x", space = "free_x")
dotplot_Control

control.disease.combined_CM_Disease <- subset(control.disease.combined_CM, subset = condition == "Disease")
dotplot_Disease <- DotPlot(control.disease.combined_CM_Disease, features = all_markers, group.by = c(RNA_resolution)) + RotatedAxis() + scale_color_gradient(low = "lightblue", high = "darkblue")
dotplot_Disease$data$gene_group <- factor(gene_group_df$group[match(dotplot_Disease$data$features.plot, gene_group_df$gene)])
dotplot_Disease <- dotplot_Disease + facet_grid(~ gene_group, scales = "free_x", space = "free_x")
dotplot_Disease
dev.off()

pdf(paste0(result_dir, "/Annotation/extract_CM/Unannotated_clustering_feature_plot_onlyCM.pdf"), width = 10, height = 6)
FeaturePlot(control.disease.combined_CM, features = marker_set1, order = T)
FeaturePlot(control.disease.combined_CM, features = marker_set3, order = T)
FeaturePlot(control.disease.combined_CM, reduction = "tsne", features = marker_set1, order = T)
FeaturePlot(control.disease.combined_CM, reduction = "tsne", features = marker_set3, order = T)
dev.off()

####################################################################
pdf(paste0(result_dir, "/Annotation/extract_CM/Unannotated_clustering_feature_plot_onlyCM_Control.pdf"), width = 10, height = 8)
FeaturePlot(control.disease.combined_CM_Control, features = marker_set1, order = T)
FeaturePlot(control.disease.combined_CM_Control, features = marker_set3, order = T)
FeaturePlot(control.disease.combined_CM_Control, reduction = "tsne", features = marker_set1, order = T)
FeaturePlot(control.disease.combined_CM_Control, reduction = "tsne", features = marker_set3, order = T)
dev.off()

pdf(paste0(result_dir, "/Annotation/extract_CM/Unannotated_clustering_feature_plot_onlyCM_Disease.pdf"), width = 10, height = 8)
FeaturePlot(control.disease.combined_CM_Disease, features = marker_set1, order = T)
FeaturePlot(control.disease.combined_CM_Disease, features = marker_set3, order = T)
FeaturePlot(control.disease.combined_CM_Disease, reduction = "tsne", features = marker_set1, order = T)
FeaturePlot(control.disease.combined_CM_Disease, reduction = "tsne", features = marker_set3, order = T)
dev.off()
####################################################################
####################################################################

temp14 <- table(control.disease.combined_CM$RNA_snn_res.0.8)
temp_sample <- table(control.disease.combined_CM$samples)
temp_sample <- temp_sample[unique(control.disease.combined_CM$samples)]
temp24 <- t(table(control.disease.combined_CM$RNA_snn_res.0.8, control.disease.combined_CM$condition))
temp34 <- t(table(control.disease.combined_CM$RNA_snn_res.0.8, control.disease.combined_CM$samples))
temp34 <- temp34[unique(control.disease.combined_CM$samples), ]

ratio_by_cluster_by_sample <- as.data.frame.matrix(temp34)
ratio_by_cluster_by_sample <- ratio_by_cluster_by_sample / temp_sample
ratio_by_cluster_by_sample <- as.data.frame.matrix(ratio_by_cluster_by_sample)
ratio_by_cluster_by_sample$condition <- c("Control", "Control", "Control", "Control", "Disease", "Disease", "Disease", "Disease", "Disease")
ratio_by_cluster_by_sample_melt <- ratio_by_cluster_by_sample %>% 
  melt() %>% 
  setNames(c("condition", "cluster", "ratio")) %>% 
  mutate(condition = factor(condition, level = c("Control", "Disease")))
# table(control.disease.combined_CM$RNA_snn_res.0.8, control.disease.combined_CM$samples)

pdf(paste0(result_dir, "/Annotation/extract_CM/Cluster_ratio_comparison.pdf"), width = 26, height = 6)
barplot_strand_ratio_age_match <- ggbarplot(ratio_by_cluster_by_sample_melt, x = "cluster", y = "ratio", color = "condition", 
                                            label = TRUE, lab.nb.digits = 2, add = c("mean_se", "jitter"), position = position_dodge(0.9), 
                                            palette = c("dodgerblue3", "#FC4E07")) +
  stat_compare_means(aes(group = condition), label = "p.format", label.y = 1.02 * max(ratio_by_cluster_by_sample_melt$ratio, na.rm = TRUE)) + 
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.06 * max(ratio_by_cluster_by_sample_melt$ratio, na.rm = TRUE)) + 
  labs(x = "Unannotated Cluster", y = "ratio of cell #") + theme_classic(base_size = 24) + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values=c("dodgerblue3", "#FC4E07"))
print(barplot_strand_ratio_age_match)
dev.off()

##### plot metagenes
# top10_FC_markers_each_cluster_CM <- FindAllMarkers(control.disease.combined_CM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.50)
# top10_FC_markers_each_cluster_CM <- top10_FC_markers_each_cluster_CM %>%
#   group_by(cluster) %>%
#   slice_max(n = 5, order_by = avg_log2FC)
all_marker_genes <- list()
all_marker_genes[["metagene_CM"]] <- marker_set1
all_marker_genes[["metagene_Fibro"]] <- marker_set3

for (name in names(all_marker_genes))
{
  print(name)
  average_expr <- colMeans(as.data.frame(GetAssayData(control.disease.combined_CM[all_marker_genes[[name]], ])))
  control.disease.combined_CM@meta.data[, name] <- average_expr
}

pdf(paste0(result_dir, "/Annotation/extract_CM/metagene_feature_plot.pdf"), width = 8, height = 6)
FeaturePlot(control.disease.combined_CM, features = "metagene_CM", order = T)
FeaturePlot(control.disease.combined_CM, features = "metagene_Fibro", order = T)

control.disease.combined_CM_Control <- subset(control.disease.combined_CM, subset = condition == "Control")
FeaturePlot(control.disease.combined_CM_Control, features = "metagene_CM", order = T)
FeaturePlot(control.disease.combined_CM_Control, features = "metagene_Fibro", order = T)

control.disease.combined_CM_Disease <- subset(control.disease.combined_CM, subset = condition == "Disease")
FeaturePlot(control.disease.combined_CM_Disease, features = "metagene_CM", order = T)
FeaturePlot(control.disease.combined_CM_Disease, features = "metagene_Fibro", order = T)
# FeaturePlot(control.disease.combined_CM, features = "DCN", order = T)
dev.off()

annotated_clusters_list_CM <- c(
  "Cardiomyocytes", "Cardiomyocytes", "Fibroblasts", "Fibroblasts", "Fibroblasts",
  "Fibroblasts", "Cardiomyocytes", "Fibroblasts", "Cardiomyocytes", "Cardiomyocytes",
  "Cardiomyocytes", "Fibroblasts", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes",
  "Fibroblasts", "Cardiomyocytes")

names(annotated_clusters_list_CM) <- paste0("cluster_", 0:(length(annotated_clusters_list_CM) - 1))
annotated_clusters_CM <- sapply(control.disease.combined_CM$RNA_snn_res.0.8, function(x) annotated_clusters_list_CM[paste0("cluster_", x)])
control.disease.combined_CM@meta.data <- cbind(control.disease.combined_CM@meta.data, annotated_clusters_CM)
cell_type_order <- c("Cardiomyocytes", "Fibroblasts")
## Reorder the cell types in the metadata
control.disease.combined_CM$annotated_clusters_CM <- factor(control.disease.combined_CM$annotated_clusters_CM, levels = cell_type_order)

pdf(paste0(result_dir, "/Annotation/extract_CM/Annotated_clusters.pdf"), width = 5.5, height = 4)
DimPlot(control.disease.combined_CM, reduction = "umap", label = T, repel = T, group.by = "annotated_clusters_CM", cols = color_list)
DimPlot(control.disease.combined_CM, reduction = "tsne", label = T, repel = T, group.by = "annotated_clusters_CM", cols = color_list)
dev.off()

Idents(object = control.disease.combined_CM) <- "annotated_clusters_CM"
pdf(paste0(result_dir, "/Annotation/extract_CM/Annotated_clustering_split_condition.pdf"), width = 8, height = 4)
DimPlot(control.disease.combined_CM, label = TRUE, split.by = "condition", cols = color_list)
dev.off()

pdf(paste0(result_dir, "/Annotation/extract_CM/Annotated_clustering_split_sample.pdf"), width = 20, height = 3)
DimPlot(control.disease.combined_CM, label = TRUE, split.by = "samples", cols = color_list) + 
  facet_grid(. ~  factor(samples, level = selected_sample_names), scale = "free_x")
# DimPlot(control.disease.combined, label = TRUE, split.by = "condition")
dev.off()

pdf(paste0(result_dir, "/Annotation/extract_CM/Annotated_dotplot.pdf"), width = 8, height = 3)
# DotPlot(control.disease.combined, features = feature_gene_list) + RotatedAxis()
dotplot_all <- DotPlot(control.disease.combined_CM, features = c(marker_set1, marker_set3), group.by = "annotated_clusters_CM") + RotatedAxis() + scale_color_gradient(low = "lightblue", high = "darkblue")
dotplot_all$data$gene_group <- factor(gene_group_df$group[match(dotplot_all$data$features.plot, gene_group_df$gene)])
dotplot_all <- dotplot_all + facet_grid(~ gene_group, scales = "free_x", space = "free_x")
dotplot_all

control.disease.combined_CM_Control <- subset(control.disease.combined_CM, subset = condition == "Control")
dotplot_Control <- DotPlot(control.disease.combined_CM_Control, features = c(marker_set1, marker_set3), group.by = "annotated_clusters_CM") + RotatedAxis() + scale_color_gradient(low = "lightblue", high = "darkblue")
dotplot_Control$data$gene_group <- factor(gene_group_df$group[match(dotplot_Control$data$features.plot, gene_group_df$gene)])
dotplot_Control <- dotplot_Control + facet_grid(~ gene_group, scales = "free_x", space = "free_x")
dotplot_Control

control.disease.combined_CM_Disease <- subset(control.disease.combined_CM, subset = condition == "Disease")
dotplot_Disease <- DotPlot(control.disease.combined_CM_Disease, features = c(marker_set1, marker_set3), group.by = "annotated_clusters_CM") + RotatedAxis() + scale_color_gradient(low = "lightblue", high = "darkblue")
dotplot_Disease$data$gene_group <- factor(gene_group_df$group[match(dotplot_Disease$data$features.plot, gene_group_df$gene)])
dotplot_Disease <- dotplot_Disease + facet_grid(~ gene_group, scales = "free_x", space = "free_x")
dotplot_Disease
dev.off()

temp15 <- table(control.disease.combined_CM$annotated_clusters_CM)
temp25 <- t(table(control.disease.combined_CM$annotated_clusters_CM, control.disease.combined_CM$condition))
temp35 <- t(table(control.disease.combined_CM$annotated_clusters_CM, control.disease.combined_CM$samples))
temp35 <- temp35[unique(control.disease.combined_CM$samples), ]

#save object
control.disease.combined_CM$celltype.condition <- paste(control.disease.combined_CM$annotated_clusters_CM, control.disease.combined_CM$condition, sep = "_")
control.disease.combined_CM$samples.condition <- paste0(control.disease.combined_CM$samples, "-", control.disease.combined_CM$condition)
saveRDS(control.disease.combined_CM, paste0(result_dir, "/Seurat.obj_with_annotation_CM_Fibro_final.RDS"))
control.disease.combined_CM  <- readRDS(paste0(result_dir, "/Seurat.obj_with_annotation_CM_Fibro_final.RDS"))

################################################################################
### only 16 and 17
Idents(object = control.disease.combined_CM) <- RNA_resolution
# control.disease.combined_CM_1617 <- subset(control.disease.combined_CM, idents = c(16, 17))
# control.disease.combined_CM_all_CM <- subset(control.disease.combined_CM, idents = c(0:29, 31, 33:39))
# control.disease.combined_CM_weak_transfer <- subset(control.disease.combined_CM, idents = c(32))
# control.disease.combined_CM_strong_transfer <- subset(control.disease.combined_CM, idents = c(30,41))
# control.disease.combined_CM_all_fib <- subset(control.disease.combined_CM, idents = c(40,42))

haha <- t(table(control.disease.combined_CM$RNA_snn_res.0.8, control.disease.combined_CM$samples))
haha <- haha[unique(control.disease.combined_CM$samples), ]

control.disease.combined_CM$transition_state <- NA
control.disease.combined_CM$transition_state[control.disease.combined_CM$RNA_snn_res.0.8 %in% c(0,1,8,9,14,16)] <- "CM"
control.disease.combined_CM$transition_state[control.disease.combined_CM$RNA_snn_res.0.8 %in% c(15)] <- "late_stage_tranfer"
control.disease.combined_CM$transition_state[control.disease.combined_CM$RNA_snn_res.0.8 %in% c(10,12)] <- "reverse_tranfer"
control.disease.combined_CM$transition_state[control.disease.combined_CM$RNA_snn_res.0.8 %in% c(6,13)] <- "early_stage_transfer"
control.disease.combined_CM$transition_state[control.disease.combined_CM$RNA_snn_res.0.8 %in% c(2:5,7,11)] <- "Fib"

# control.disease.combined_CM$transition_state <- NA
# control.disease.combined_CM$transition_state[control.disease.combined_CM$RNA_snn_res.0.8 %in% c(1,2,5:9,14:15,19,22:23,26,32:35)] <- "CM"
# control.disease.combined_CM$transition_state[control.disease.combined_CM$RNA_snn_res.0.8 %in% c(31)] <- "late_stage_tranfer"
# control.disease.combined_CM$transition_state[control.disease.combined_CM$RNA_snn_res.0.8 %in% c(27,28,29)] <- "reverse_tranfer"
# control.disease.combined_CM$transition_state[control.disease.combined_CM$RNA_snn_res.0.8 %in% c(16)] <- "early_stage_transfer"
# control.disease.combined_CM$transition_state[control.disease.combined_CM$RNA_snn_res.0.8 %in% c(0,3,4,10:13,17:18,20:21,24:25,30)] <- "Fib"

pdf(paste0(result_dir, "/Annotation/extract_CM/Unannotated_clustering_transition_state.pdf"), width = 20, height = 6)
DimPlot(control.disease.combined_CM, reduction = "umap", group.by = "transition_state", label = TRUE, label.size = 2, cols = color_list, order = T) + 
  FeaturePlot(control.disease.combined_CM, features = c("metagene_CM", "metagene_Fibro"), order = T)
control.disease.combined_CM_Control <- subset(control.disease.combined_CM, subset = condition == "Control")
DimPlot(control.disease.combined_CM_Control, reduction = "umap", group.by = "transition_state", label = TRUE, label.size = 2, cols = color_list, order = T) + 
  FeaturePlot(control.disease.combined_CM_Control, features = c("metagene_CM", "metagene_Fibro"), order = T)
control.disease.combined_CM_Disease <- subset(control.disease.combined_CM, subset = condition == "Disease")
DimPlot(control.disease.combined_CM_Disease, reduction = "umap", group.by = "transition_state", label = TRUE, label.size = 2, cols = color_list, order = T) + 
  FeaturePlot(control.disease.combined_CM_Disease, features = c("metagene_CM", "metagene_Fibro"), order = T)

DimPlot(control.disease.combined_CM, reduction = "tsne", group.by = "transition_state", label = TRUE, label.size = 2, cols = color_list, order = T) + 
  FeaturePlot(control.disease.combined_CM, reduction = "tsne", features = c("metagene_CM", "metagene_Fibro"), order = T)
control.disease.combined_CM_Control <- subset(control.disease.combined_CM, subset = condition == "Control")
DimPlot(control.disease.combined_CM_Control, reduction = "tsne", group.by = "transition_state", label = TRUE, label.size = 2, cols = color_list, order = T) + 
  FeaturePlot(control.disease.combined_CM_Control, reduction = "tsne", features = c("metagene_CM", "metagene_Fibro"), order = T)
control.disease.combined_CM_Disease <- subset(control.disease.combined_CM, subset = condition == "Disease")
DimPlot(control.disease.combined_CM_Disease, reduction = "tsne", group.by = "transition_state", label = TRUE, label.size = 2, cols = color_list, order = T) + 
  FeaturePlot(control.disease.combined_CM_Disease, reduction = "tsne", features = c("metagene_CM", "metagene_Fibro"), order = T)
dev.off()

##### meatagene expression for each cluster
average_metagene_CM_expr_per_sample <- aggregate(as.matrix(control.disease.combined_CM$metagene_CM), by = list(control.disease.combined_CM$samples), FUN = mean)
average_metagene_CM_expr_per_sample <- average_metagene_CM_expr_per_sample[order(match(average_metagene_CM_expr_per_sample$Group.1, unique(control.disease.combined_CM$samples))), ]
t.test(average_metagene_CM_expr_per_sample$V1[1:4], average_metagene_CM_expr_per_sample$V1[5:9])

average_metagene_Fibro_expr_per_sample <- aggregate(as.matrix(control.disease.combined_CM$metagene_Fibro), by = list(control.disease.combined_CM$samples), FUN = mean)
average_metagene_Fibro_expr_per_sample <- average_metagene_Fibro_expr_per_sample[order(match(average_metagene_Fibro_expr_per_sample$Group.1, unique(control.disease.combined_CM$samples))), ]
t.test(average_metagene_Fibro_expr_per_sample$V1[1:4], average_metagene_Fibro_expr_per_sample$V1[5:9])

average_metagene_CM_expr_per_sample_control <- aggregate(as.matrix(control.disease.combined_CM_Control$metagene_CM), by = list(control.disease.combined_CM_Control$RNA_snn_res.0.8), FUN = mean)
average_metagene_CM_expr_per_sample_disease <- aggregate(as.matrix(control.disease.combined_CM_Disease$metagene_CM), by = list(control.disease.combined_CM_Disease$RNA_snn_res.0.8), FUN = mean)
average_metagene_Fibro_expr_per_sample_control <- aggregate(as.matrix(control.disease.combined_CM_Control$metagene_Fibro), by = list(control.disease.combined_CM_Control$RNA_snn_res.0.8), FUN = mean)
average_metagene_Fibro_expr_per_sample_disease <- aggregate(as.matrix(control.disease.combined_CM_Disease$metagene_Fibro), by = list(control.disease.combined_CM_Disease$RNA_snn_res.0.8), FUN = mean)


grouping <- paste(control.disease.combined_CM_Control$samples, control.disease.combined_CM_Control$RNA_snn_res.0.8, sep = "_")
average_expr_per_group <- aggregate(as.matrix(control.disease.combined_CM_Control$metagene_CM), by = list(grouping), FUN = mean)
colnames(average_expr_per_group)[1] <- "sample_cluster"
long_expr_data1 <- melt(average_expr_per_group, id.vars = "sample_cluster", variable.name = "gene", value.name = "average_expression")
long_expr_data1 <- tidyr::separate(long_expr_data1, sample_cluster, into = c("sample", "cluster"), sep = "_")

grouping <- paste(control.disease.combined_CM_Disease$samples, control.disease.combined_CM_Disease$RNA_snn_res.0.8, sep = "_")
average_expr_per_group <- aggregate(as.matrix(control.disease.combined_CM_Disease$metagene_CM), by = list(grouping), FUN = mean)
colnames(average_expr_per_group)[1] <- "sample_cluster"
long_expr_data2 <- melt(average_expr_per_group, id.vars = "sample_cluster", variable.name = "gene", value.name = "average_expression")
long_expr_data2 <- tidyr::separate(long_expr_data2, sample_cluster, into = c("sample", "cluster"), sep = "_")

wilcox.test(long_expr_data1[long_expr_data1$cluster == 6, 4], long_expr_data2[long_expr_data2$cluster == 6, 4])
# Plot the average expression for each gene across samples and clusters
# ggplot(long_expr_data, aes(x = sample, y = average_expression, fill = gene)) +
#   geom_point() + 
#   facet_wrap(~ cluster, scales = "free_x") + 
#   labs(x = "Sample", y = "Average Gene Expression", 
#        title = "Average Gene Expression per Sample and Cluster") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the average expression for each gene across samples and clusters
# ggplot(long_expr_data, aes(x = sample, y = average_expression, fill = gene)) +
#   geom_point() + 
#   facet_wrap(~ cluster, scales = "free_x") + 
#   labs(x = "Sample", y = "Average Gene Expression", 
#        title = "Average Gene Expression per Sample and Cluster") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf(paste0(result_dir, "/Annotation/extract_CM/metagene_control_vs_disease.pdf"), width = 8, height = 5)
ggplot(average_metagene_CM_expr_per_sample_control, aes(x = Group.1, y = V1)) +
  geom_point(size = 3, color = "blue") + geom_point(data = average_metagene_CM_expr_per_sample_disease, aes(x = Group.1, y = V1), color = "red", size = 3) + 
  ggtitle("metagene_CM_control_vs_disease")
ggplot(average_metagene_Fibro_expr_per_sample_control, aes(x = Group.1, y = V1)) +
  geom_point(size = 3, color = "blue") + geom_point(data = average_metagene_Fibro_expr_per_sample_disease, aes(x = Group.1, y = V1), color = "red", size = 3) + 
  ggtitle("metagene_Fibro_control_vs_disease")
dev.off()
# pdf(paste0(result_dir, "/Annotation/extract_CM/Unannotated_clustering_transition.pdf"), width = 6, height = 5)
# DimPlot(control.disease.combined_CM_all_CM, reduction = "umap", group.by = RNA_resolution, label = TRUE, label.size = 4, cols = color_list, order = T) + xlim(-8, 8) + ylim(-8, 8)
# DimPlot(control.disease.combined_CM_weak_transfer, reduction = "umap", group.by = RNA_resolution, label = TRUE, label.size = 4, cols = color_list, order = T) + xlim(-8, 8) + ylim(-8, 8)
# DimPlot(control.disease.combined_CM_strong_transfer, reduction = "umap", group.by = RNA_resolution, label = TRUE, label.size = 4, cols = color_list, order = T) + xlim(-8, 8) + ylim(-8, 8)
# DimPlot(control.disease.combined_CM_all_fib, reduction = "umap", group.by = RNA_resolution, label = TRUE, label.size = 4, cols = color_list, order = T) + xlim(-8, 8) + ylim(-8, 8)
# 
# DimPlot(control.disease.combined_CM_all_CM, reduction = "tsne", group.by = RNA_resolution, label = TRUE, label.size = 4, cols = color_list, order = T)
# DimPlot(control.disease.combined_CM_weak_transfer, reduction = "tsne", group.by = RNA_resolution, label = TRUE, label.size = 4, cols = color_list, order = T)
# DimPlot(control.disease.combined_CM_strong_transfer, reduction = "tsne", group.by = RNA_resolution, label = TRUE, label.size = 4, cols = color_list, order = T)
# DimPlot(control.disease.combined_CM_all_fib, reduction = "tsne", group.by = RNA_resolution, label = TRUE, label.size = 4, cols = color_list, order = T)
# dev.off()


# pdf(paste0(result_dir, "/Annotation/extract_CM/Unannotated_clustering_only1617.pdf"), width = 6, height = 5)
# DimPlot(control.disease.combined_CM_1617, reduction = "umap", group.by = RNA_resolution, label = TRUE, label.size = 4, cols = color_list, order = T) +
#   theme(text = element_text(size = 6), axis.text = element_text(size = 6))
# DimPlot(control.disease.combined_CM_1617, reduction = "tsne", group.by = RNA_resolution, label = TRUE, label.size = 4, cols = color_list) +
#   theme(text = element_text(size = 6), axis.text = element_text(size = 6))
# dev.off()
# pdf(paste0(result_dir, "/Annotation/extract_CM/Unannotated_clustering_split_condition_only1617.pdf"), width = 8, height = 4)
# DimPlot(control.disease.combined_CM_1617, reduction = "umap", group.by = RNA_resolution, label = TRUE, split.by = "condition", label.size = 4, cols = color_list)
# DimPlot(control.disease.combined_CM_1617, reduction = "tsne", group.by = RNA_resolution, label = TRUE, split.by = "condition", label.size = 4, cols = color_list)
# dev.off()
# 
# pdf(paste0(result_dir, "/Annotation/extract_CM/Unannotated_clustering_feature_plot_only1617.pdf"), width = 10, height = 6)
# FeaturePlot(control.disease.combined_CM_1617, features = marker_set1, order = T)
# FeaturePlot(control.disease.combined_CM_1617, features = marker_set3, order = T)
# FeaturePlot(control.disease.combined_CM_1617, reduction = "tsne", features = marker_set1, order = T)
# FeaturePlot(control.disease.combined_CM_1617, reduction = "tsne", features = marker_set3, order = T)
# dev.off()
# 
# 
# control.disease.combined_CM_16 <- subset(control.disease.combined_CM, idents = c(16))
# control.disease.combined_CM_17 <- subset(control.disease.combined_CM, idents = c(17))
# Idents(control.disease.combined_CM_16) <- control.disease.combined_CM_16$condition
# Idents(control.disease.combined_CM_17) <- control.disease.combined_CM_17$condition
# marker_gene_celltype_CM_16 <- FindMarkers(control.disease.combined_CM_16, only.pos = TRUE, ident.1 = "Disease", ident.2 = "Control", min.pct = 0.1, logfc.threshold = 0.5)
# marker_gene_celltype_CM_17 <- FindMarkers(control.disease.combined_CM_17, only.pos = TRUE, ident.1 = "Disease", ident.2 = "Control", min.pct = 0.1, logfc.threshold = 0.5)
# control.disease.combined_CM_16$samples.condition <- factor(x = control.disease.combined_CM_16$samples.condition, 
#                                                            levels = c("936-Control", "5087-Control", "5828-Control", "5919-Control", 
#                                                                       "1156-Disease", "604-Disease", "5111-Disease", "5364-Disease", "5874-Disease"))
# control.disease.combined_CM_17$samples.condition <- factor(x = control.disease.combined_CM_17$samples.condition, 
#                                                            levels = c("936-Control", "5087-Control", "5828-Control", "5919-Control", 
#                                                                       "1156-Disease", "604-Disease", "5111-Disease", "5364-Disease", "5874-Disease"))
# 
# p_16 <- DoHeatmap(control.disease.combined_CM_16, features = rownames(marker_gene_celltype_CM_16)[1:9], group.by = "samples.condition", slot = "count") + labs(title = "")
# p_17 <- DoHeatmap(control.disease.combined_CM_17, features = rownames(marker_gene_celltype_CM_17)[1:3], group.by = "samples.condition", slot = "count") + labs(title = "")
