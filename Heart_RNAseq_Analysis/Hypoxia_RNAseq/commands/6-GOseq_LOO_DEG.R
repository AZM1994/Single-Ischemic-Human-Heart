library(goseq)
library(ggplot2)
library(reshape2)
library(stringr)
library(biomaRt)
library(DOSE)
library(dplyr)
# library(tidyr)
library(ggrepel)
library(patchwork)
library(Seurat)

wd_dir <- getwd()
LOO_DEG_GOseq_dir <- paste0(wd_dir, "/results_IHD/6-LOO_DEG_GOseq")
dir.create(LOO_DEG_GOseq_dir, recursive = TRUE, showWarnings = FALSE)

numDEInCat_threshold = 2
numInCat_threshold = 1000
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])

##### read DEGs
DEG_df <- read.csv(paste0(LOO_DEG_GOseq_dir, "/LOO_DEG_filtered.csv"))
DEG_up <- DEG_df$Gene[DEG_df$regulation == "up"]
DEG_down <- DEG_df$Gene[DEG_df$regulation == "down"]
### read RNAseq CM data
Seurat.obj_sub_clustering_CM_only <- readRDS(paste0(wd_dir, "/results_IHD/0-Seurat_Object/Seurat.obj_sub_clustering_CM_only.RDS"))
expr_data <- as.matrix(GetAssayData(object = Seurat.obj_sub_clustering_CM_only, assay = "RNA", slot = "count"))
expressed_genes <- row.names(expr_data)
gene_expression_percentage <- Matrix::rowSums(expr_data > 0) / ncol(expr_data) * 100
# expressed_genes <- names(gene_expression_percentage[gene_expression_percentage > 0])

################################################################################
################################################################################
##### GOseq and filtering for UP
DEG_up_list <- as.integer(expressed_genes %in% DEG_up)
names(DEG_up_list) <- expressed_genes
# pwf_up <- nullp(DEG_up, "hg19", bias.data = gene_length_deduped[, gene_length_type])
pwf_up <- nullp(DEG_up_list, "hg19", id = "geneSymbol", bias.data = NULL)
DEG_up_GO <- goseq(pwf_up, "hg19", "geneSymbol") %>% 
  mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% 
  mutate(FDR = p.adjust(over_represented_pvalue, method = "fdr")) %>% 
  mutate(Regulation = "UP")

## find gene names in each GO term: Control
go_term_list_up <- DEG_up_GO$category  # GO term IDs from goseq
go_genes_list_up <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list_up, keytype = "GOALL", columns = c("SYMBOL"))
genes_in_GO_up <- go_genes_list_up[go_genes_list_up$SYMBOL %in% DEG_up, ]
genes_in_GO_up_collapsed <- aggregate(SYMBOL ~ GOALL, data = genes_in_GO_up, function(x) paste(unique(x), collapse = ", "))
DEG_up_GO_with_genes <- merge(DEG_up_GO, genes_in_GO_up_collapsed, by.x = "category", by.y = "GOALL") 
DEG_up_GO_filtered_with_genes <- DEG_up_GO_with_genes %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
  filter(over_represented_pvalue < 0.05)

write.csv(DEG_up_GO_filtered_with_genes, paste0(LOO_DEG_GOseq_dir, "/DEG_up_GO_filtered_with_genes.csv"))

################################################################################
##### GOseq and filtering for DOWN
DEG_down_list <- as.integer(expressed_genes %in% DEG_down)
names(DEG_down_list) <- expressed_genes
# pwf_down <- nullp(DEG_down_list, bias.data = gene_length_deduped[, gene_length_type])
pwf_down <- nullp(DEG_down_list, "hg19", id = "geneSymbol", bias.data = NULL)
DEG_down_GO <- goseq(pwf_down, "hg19", "geneSymbol") %>% 
  mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% 
  mutate(FDR = p.adjust(over_represented_pvalue, method = "fdr")) %>% 
  mutate(Regulation = "DOWN")

## find gene names in each GO term: IHD
go_term_list_down <- DEG_down_GO$category  # GO term IDs from goseq
go_genes_list_down <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list_down, keytype = "GOALL", columns = c("SYMBOL"))
genes_in_GO_down <- go_genes_list_down[go_genes_list_down$SYMBOL %in% DEG_down, ]
genes_in_GO_down_collapsed <- aggregate(SYMBOL ~ GOALL, data = genes_in_GO_down, function(x) paste(unique(x), collapse = ", "))
DEG_down_GO_with_genes <- merge(DEG_down_GO, genes_in_GO_down_collapsed, by.x = "category", by.y = "GOALL")
DEG_down_GO_filtered_with_genes <- DEG_down_GO_with_genes %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>% 
  filter(over_represented_pvalue < 0.05)

write.csv(DEG_down_GO_filtered_with_genes, paste0(LOO_DEG_GOseq_dir, "/DEG_down_GO_filtered_with_genes.csv"))

################################################################################
### generate volcano plots for BP
nature_colors <- c("#a2d2e7", "#67a8cd", "#ffc17f", "#cf9f88", "#8c564b", "#6fb3a8", "#b3e19b","#50aa4b","#ff9d9f","#f36569","#3581b7","#cdb6da",
                   "#704ba3", "#9a7fbd", "#e377c2",  "#dba9a8", "#e43030", "#e99b78", "#ff8831", "#bcbd22")
warm_colors <- c("#FFD700", "#FF6F61", "#BFA2DB", "#98DD62", "#7EC8E3")
cool_colors <- c("#FF9445", "#E64545", "#B151A2", "#6A974F", "#3F7CA7")

warm_colors1 <- c("#FFD700", "#FF6F61", "#BFA2DB", "#98DD62", "#7EC8E3")
cool_colors1 <- c("#FF9445", "#E64545", "#B151A2", "#6A974F", "#3F7CA7")

## DEG UP
Oxidative_Stress_and_Inflammatory_Response_up <- c(3, 10, 13, 33, 35, 47, 54, 55, 56, 58, 61, 74, 75, 110, 116, 117, 118, 119, 120, 124, 126, 127, 135, 137, 
                                                   144, 145, 153, 154, 162, 188, 189, 192, 196, 224, 235, 236, 237, 238, 259, 273, 309)
Angiogenesis_and_Vascular_Remodeling_up <- c(9, 15, 16, 17, 18, 57, 59, 60, 64, 65, 66, 68, 69, 71, 72, 142, 180, 202, 204, 215, 225, 227, 228, 239, 241, 
                                             244, 245, 246, 248, 249, 253, 254, 255, 295)
Autophagy_and_Proteostasis_up <- c(1, 4, 5, 32, 36, 37, 38, 63, 67, 79, 80, 81, 82, 97, 99, 108, 109, 155, 167, 190, 217, 218, 219, 257, 289, 292, 294, 297, 300, 
                                   301, 302, 313)
Metabolic_Reprogramming_up <- c(25, 26, 28, 29, 30, 31, 46, 49, 51, 52, 53, 77, 78, 84, 85, 103, 111, 114, 115, 121, 122, 125, 131, 132, 133, 136, 139, 
                                140, 141, 151, 152, 156, 171, 172, 174, 177, 178, 184, 185, 191, 199, 240, 252, 256, 270, 271, 284, 286, 305, 310, 311, 312)
Signaling_and_Structural_Remodeling_up <- c(2, 6, 7, 11, 12, 19, 20, 21, 22, 23, 24, 27, 34, 39, 40, 41, 42, 43, 44, 45, 48, 50, 62, 70, 73, 76, 83, 86, 87, 
                                            88, 89, 90, 91, 92, 98, 100, 101, 102, 104, 105, 106, 107, 128, 129, 130, 134, 138, 143, 146, 147, 148, 149, 150, 
                                            157, 158, 159, 160, 161, 163, 164, 165, 166, 168, 169, 170, 173, 175, 176, 179, 181, 182, 183, 186, 187, 193, 194, 
                                            195, 197, 198, 200, 201, 205, 206, 209, 211, 213, 216, 220, 221, 222, 223, 226, 229, 230, 231, 232, 233, 234, 260, 
                                            261, 262, 263, 264, 265, 266, 267, 268, 269, 272, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 285, 287, 288, 
                                            290, 291, 293, 296, 298, 299, 303, 304, 306, 307, 308)

GO_BP_up <- DEG_up_GO_filtered_with_genes %>% 
  mutate(GO_Group = case_when(row.names(.) %in% Oxidative_Stress_and_Inflammatory_Response_up ~ "Oxidative Stress & Inflammatory Response", 
                              row.names(.) %in% Angiogenesis_and_Vascular_Remodeling_up ~ "Angiogenesis & Vascular Remodeling", 
                              row.names(.) %in% Autophagy_and_Proteostasis_up ~ "Autophagy & Proteostasis", 
                              row.names(.) %in% Metabolic_Reprogramming_up ~ "Metabolic Reprogramming", 
                              row.names(.) %in% Signaling_and_Structural_Remodeling_up ~ "Signaling & Structural Remodeling", TRUE ~ NA_character_), 
         regulation = "IHD Up-regulated") %>% 
  # filter(row.names(.) %in% c(1:300)) %>% dplyr::select(term, GO_Group)
  filter(!is.na(GO_Group), hitsPerc >= 2, ontology == "BP") %>% arrange(GO_Group)

## DEG DOWN
Muscle_Contraction_and_Sarcomere_Organization_down <- c(19, 20, 21, 22, 23, 24, 25, 26, 27, 37, 40, 41, 52, 53, 54, 55, 63, 74, 77, 78, 99, 100, 105, 109, 114, 
                                                        115, 116, 117, 119, 122, 132, 139, 142, 148, 157, 183, 185, 187, 188, 189, 190, 191, 193, 194, 195, 204, 213, 214)
Cytoskeleton_and_Cell_Adhesion_Remodeling_down <- c(6, 29, 32, 33, 42, 43, 44, 56, 57, 58, 64, 65, 66, 68, 81, 101, 102, 103, 104, 106, 107, 108, 111, 112, 121, 
                                                    124, 125, 126, 127, 143, 145, 146, 147, 149, 150, 151, 158, 166, 174, 181, 182, 184, 202, 203, 215, 216, 218, 221, 228, 231)
Cardiac_and_Muscle_Development_Morphogenesis_down <- c(1, 2, 3, 4, 7, 8, 17, 18, 28, 59, 60, 61, 62, 75, 79, 83, 91, 92, 93, 94, 95, 96, 110, 123, 128, 129, 130, 
                                                       133, 134, 152, 153, 154, 155, 156, 159, 162, 163, 170, 171, 172, 175, 176, 177, 186, 198, 199, 200, 210, 
                                                       211, 212, 219, 226, 227, 229, 237, 238, 239)
Metabolic_Transport_and_Proteostatic_Processes_down <- c(10, 34, 35, 36, 45, 46, 47, 48, 49, 50, 51, 69, 70, 71, 76, 80, 82, 84, 85, 86, 87, 88, 89, 90, 97, 120, 
                                                         135, 136, 144, 192, 201, 207, 208, 224, 232, 233, 234, 235, 236)
Immune_Stress_Regulation_and_Apoptosis_down <- c(9, 11, 12, 13, 14, 15, 16, 113, 137, 138, 140, 160, 161, 164, 165, 167, 168, 173, 178, 179, 180, 196, 197, 205, 
                                                 206, 209, 220, 222, 223, 225, 230)

GO_BP_down <- DEG_down_GO_filtered_with_genes %>% 
  mutate(GO_Group = case_when(row.names(.) %in% Muscle_Contraction_and_Sarcomere_Organization_down ~ "Muscle Contraction & Sarcomere Organization", 
                              row.names(.) %in% Cytoskeleton_and_Cell_Adhesion_Remodeling_down ~ "Cytoskeleton & Cell Adhesion Remodeling", 
                              row.names(.) %in% Cardiac_and_Muscle_Development_Morphogenesis_down ~ "Cardiac & Muscle Development/Morphogenesis", 
                              row.names(.) %in% Metabolic_Transport_and_Proteostatic_Processes_down ~ "Metabolic, Transport & Proteostatic Processes", 
                              row.names(.) %in% Immune_Stress_Regulation_and_Apoptosis_down ~ "Immune Stress Regulation & Apoptosis", TRUE ~ NA_character_), 
         regulation = "IHD Down-regulated") %>% 
  # filter(row.names(.) %in% c(1:300)) %>% dplyr::select(term, GO_Group)
  filter(!is.na(GO_Group), hitsPerc >= 2, ontology == "BP") %>% arrange(GO_Group)

GO_BP_up_labels <- c("response to reactive oxygen species", "interleukin-1 production", # Oxidative Stress & Inflammatory Response
                     "positive regulation of vascular endothelial growth factor production", # Angiogenesis & Vascular Remodeling
                     "autophagy of mitochondrion", # Autophagy & Proteostasis
                     "mRNA alternative polyadenylation", "regulation of fatty acid oxidation", # Metabolic Reprogramming
                     "regulation of cardiac conduction", "postsynaptic cytoskeleton organization") # Signaling & Structural Remodeling

GO_BP_down_labels <- c("myofibril assembly", "striated muscle cell differentiation", # Muscle Contraction & Sarcomere Organization
                     "actin cytoskeleton organization", "actomyosin structure organization", # Cytoskeleton & Cell Adhesion Remodeling
                     "muscle cell development", "heart morphogenesis", # Cardiac & Muscle Development/Morphogenesis
                     "carnitine shuttle", "gamma-aminobutyric acid biosynthetic process", # Metabolic, Transport & Proteostatic Processes
                     "susceptibility to natural killer cell mediated cytotoxicity", "hepatocyte growth factor receptor signaling pathway") # Immune Stress Regulation & Apoptosis

write.csv(GO_BP_up, paste0(LOO_DEG_GOseq_dir, "/IHD_UP_GO.csv"))
write.csv(GO_BP_down, paste0(LOO_DEG_GOseq_dir, "/IHD_DOWN_GO.csv"))

nature_colors <- c("#a2d2e7", "#67a8cd", "#ffc17f", "#cf9f88", "#8c564b", "#6fb3a8", "#b3e19b","#50aa4b","#ff9d9f","#f36569","#3581b7","#cdb6da",
                   "#704ba3", "#9a7fbd", "#e377c2",  "#dba9a8", "#e43030", "#e99b78", "#ff8831", "#bcbd22")
warm_colors <- c("#FFD700", "#FF6F61", "#cf9f88", "#98DD62", "#7EC8E3")
cool_colors <- c("#E64545", "#ffc17f", "#B151A2", "#bcbd22", "#3F7CA7")

p_GO_IHD_up <- ggplot(GO_BP_up, aes(x = hitsPerc, y = -log10(over_represented_pvalue), color = GO_Group)) +
  geom_point(size = 3.5) + geom_vline(xintercept = 2, color = "gray40", linetype = "22") + geom_hline(yintercept = -log10(0.05), color = "gray40", linetype = "22") +
  geom_text_repel(data = subset(GO_BP_up, term %in% GO_BP_up_labels), aes(label = term),
                  nudge_x = c(1.5, 1.2, 0.5, 0.2, 1.8, -0.25, 0.5, 0.5),
                  nudge_y = c(0.2, 0.2, 0.3, 0.0, 0.1, 0.25, 0.2, 0.2),
                  size = 5, color = "black", show.legend = F) +
  scale_color_manual(values = warm_colors) + labs(x = "Hits Percentage (%)", y = "-log10(p-value)", color = "GO Group", title = "") + facet_wrap(~ regulation, scales = "fixed") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(1.3, 5.5), breaks = c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5)) +
  scale_x_continuous(trans = "log2", limits = c(1.97, 105), breaks = c(2, 4, 8, 16, 32, 64)) +
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5),
                           text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
ggsave(paste0(LOO_DEG_GOseq_dir, "/GO_IHD_up.pdf"), plot = p_GO_IHD_up, width = 14, height = 7.5, dpi = 600)

p_GO_IHD_down <- ggplot(GO_BP_down, aes(x = hitsPerc, y = -log10(over_represented_pvalue), color = GO_Group)) +
  geom_point(size = 3.5) + geom_vline(xintercept = 2, color = "gray40", linetype = "22") + geom_hline(yintercept = -log10(0.05), color = "gray40", linetype = "22") +
  geom_text_repel(data = subset(GO_BP_down, term %in% GO_BP_down_labels), aes(label = term), 
                  nudge_x = c(1.2, 1.2, 1.4, 1.0, -1.5, 0.8, -0.8, -0.25, 0.5, 0.5), 
                  nudge_y = c(0.1, 0.1, -0.1, -0.3, 0.8, -0.3, 0.0, 0.25, 0.2, 0.2), 
                  size = 5, color = "black", show.legend = F) +
  scale_color_manual(values = cool_colors) + labs(x = "Hits Percentage (%)", y = "-log10(p-value)", color = "GO Group", title = "") + facet_wrap(~ regulation, scales = "fixed") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(1.3, 5.5), breaks = c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5) ) +
  scale_x_continuous(trans = "log2", limits = c(1.97, 105), breaks = c(2, 4, 8, 16, 32, 64)) +
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5),
                           text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
ggsave(paste0(LOO_DEG_GOseq_dir, "/GO_IHD_down.pdf"), plot = p_GO_IHD_down, width = 15, height = 8, dpi = 600)
