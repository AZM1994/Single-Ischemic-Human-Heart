library(goseq)
library(ggplot2)
library(reshape2)
library(stringr)
library(biomaRt)
library(DOSE)
library(dplyr)
library(tidyr)
library(ggrepel)
library(patchwork)
library(Seurat)

wd_dir <- getwd()
result_dir <- paste0(wd_dir, "/results/LOO_DEG_Cardiomyocyte/GOseq")
dir.create(result_dir, recursive = T)
numDEInCat_threshold = 2
numInCat_threshold = 1000
# gene_length_type = "Gene_length"
# gene_length_type = "Exon_length"
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])

##### read gene length
# gene_length <- read.delim(paste0(wd_dir, "/data/hg19_refGene.length.tsv"), header = F)
# colnames(gene_length) <- c("Gene_symbol", "Transcript", "Gene_length", "Exon_length")
# gene_length_deduped <- gene_length[!duplicated(gene_length$Gene_symbol),]
##### read DEGs
DEG_df <- read.csv(paste0(wd_dir, "/results/LOO_DEG_Cardiomyocyte/DEG_up_down_df.csv"))
DEG_up <- DEG_df$gene[DEG_df$regulation == "up"]
DEG_down <- DEG_df$gene[DEG_df$regulation == "down"]
### read RNAseq CM data
Seurat.obj_sub_clustering_CM_only <- readRDS(paste0(wd_dir, "/results/Seurat.obj_sub_clustering_CM_only.RDS"))
expr_data <- as.matrix(GetAssayData(object = Seurat.obj_sub_clustering_CM_only, assay = "RNA", slot = "count"))
expressed_genes <- row.names(expr_data)
gene_expression_percentage <- Matrix::rowSums(expr_data > 0) / ncol(expr_data) * 100
expressed_genes <- names(gene_expression_percentage[gene_expression_percentage > 0])

################################################################################
################################################################################
##### GOseq and filtering for UP
DEG_up_list <- as.integer(expressed_genes %in% DEG_up)
names(DEG_up_list) <- expressed_genes
# pwf_up <- nullp(DEG_up, "hg19", bias.data = gene_length_deduped[, gene_length_type])
pwf_up <- nullp(DEG_up_list, "hg19", id = "geneSymbol", bias.data = NULL)
DEG_up_GO <- goseq(pwf_up, "hg19", "geneSymbol") %>% mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% mutate(Regulation = "UP")
# mut_gene_normal_GO_filtered <- DEG_up_GO %>%
#   filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
#   filter(over_represented_pvalue < 0.05) %>%
#   mutate(Condition = "Control")

## find gene names in each GO term: Control
go_term_list_up <- DEG_up_GO$category  # GO term IDs from goseq
go_genes_list_up <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list_up, keytype = "GOALL", columns = c("SYMBOL"))
genes_in_GO_up <- go_genes_list_up[go_genes_list_up$SYMBOL %in% DEG_up, ]
genes_in_GO_up_collapsed <- aggregate(SYMBOL ~ GOALL, data = genes_in_GO_up, function(x) paste(unique(x), collapse = ", "))
DEG_up_GO_with_genes <- merge(DEG_up_GO, genes_in_GO_up_collapsed, by.x = "category", by.y = "GOALL") 
DEG_up_GO_filtered_with_genes <- DEG_up_GO_with_genes %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
  filter(over_represented_pvalue < 0.05)

# write.csv(DEG_up_GO_with_genes, paste0(result_dir, "/DEG_up_GO_with_genes.csv"))
write.csv(DEG_up_GO_filtered_with_genes, paste0(result_dir, "/DEG_up_GO_filtered_with_genes.csv"))

################################################################################
################################################################################
##### GOseq and filtering for DOWN
DEG_down_list <- as.integer(expressed_genes %in% DEG_down)
names(DEG_down_list) <- expressed_genes
# pwf_down <- nullp(DEG_down_list, bias.data = gene_length_deduped[, gene_length_type])
pwf_down <- nullp(DEG_down_list, "hg19", id = "geneSymbol", bias.data = NULL)
DEG_down_GO <- goseq(pwf_down, "hg19", "geneSymbol") %>% mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% mutate(Regulation = "DOWN")
# mut_gene_disease_GO_filtered <- DEG_down_GO %>%
#   filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
#   filter(over_represented_pvalue < 0.05) %>%
#   mutate(Condition = "IHD")

## find gene names in each GO term: IHD
go_term_list_down <- DEG_down_GO$category  # GO term IDs from goseq
go_genes_list_down <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list_down, keytype = "GOALL", columns = c("SYMBOL"))
genes_in_GO_down <- go_genes_list_down[go_genes_list_down$SYMBOL %in% DEG_down, ]
# genes_in_GO_down_collapsed <- genes_in_GO_down %>% group_by(GOALL) %>% 
#   filter(!duplicated(SYMBOL)) %>% 
#   summarise(genes = paste(SYMBOL, collapse = ", "))
genes_in_GO_down_collapsed <- aggregate(SYMBOL ~ GOALL, data = genes_in_GO_down, function(x) paste(unique(x), collapse = ", "))
DEG_down_GO_with_genes <- merge(DEG_down_GO, genes_in_GO_down_collapsed, by.x = "category", by.y = "GOALL")
DEG_down_GO_filtered_with_genes <- DEG_down_GO_with_genes %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>% 
  filter(over_represented_pvalue < 0.05)

# write.csv(DEG_down_GO_with_genes, paste0(result_dir, "/DEG_down_GO_with_genes.csv"))
write.csv(DEG_down_GO_filtered_with_genes, paste0(result_dir, "/DEG_down_GO_filtered_with_genes.csv"))

################################################################################
################################################################################
### generate volcano plots for BP
# nature_colors <- c("#a6761d", "#228C68", "#1f77b4", "#ff7f0e", "#e377c2", "#F22020")
# nature_colors <- c("#4d79a6", "#a1cde6", "#e05759", "#ff9d99", "#f1cd63", "#79706d","#489893","#8cd17d","#00afba","#e6b800","#fb4d07")
color_list <- c("#d62728", "#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#e377c2", "#8c564b", "#7f7f7f", "#bcbd22", "#17becf",
                "#393b79", "#5254a3", "#6b6ecf", "#9c9ede", "#637939", "#8ca252", "#b5cf6b", "#cedb9c", "#8c6d31", "#bd9e39",
                "#e7ba52", "#e7969c", "#d6616b", "#ad494a", "#843c39", "#f5f5f5", "#fdd0a2", "#fb6a4a", "#cb181d", "#a50f15",
                "#3182bd", "#6baed6", "#9ecae1", "#c6dbef", "#e6550d", "#fd8d3c", "#fdae6b", "#fdd0a2", "#31a354", "#74c476",
                "#a1d99b", "#c7e9c0", "#756bb1", "#9e9ac8", "#bcbddc")

color_list <- c("#e6550d", "#9ecae1", "#a1d99b", "#d6616b", "#17becf", "#fd8d3c", "#9c9ede", "#1f77b4", "#8c6d31", "#d62728",
                "#bcbd22", "#8ca252", "#e7ba52", "#e377c2", "#9467bd", "#2ca02c", "#393b79", "#bd9e39", "#a50f15", "black")

nature_colors <- c("#a2d2e7", "#67a8cd", "#ffc17f", "#cf9f88", "#8c564b", "#6fb3a8", "#b3e19b","#50aa4b","#ff9d9f","#f36569","#3581b7","#cdb6da",
                   "#704ba3", "#9a7fbd", "#e377c2",  "#dba9a8", "#e43030", "#e99b78", "#ff8831", "#bcbd22")
## DEG UP
Ribosomal_activity_Up <- c(236, 235, 3, 2, 8, 325, 149, 213, 214, 223, 198, 169, 55, 146, 91, 114, 102, 113)
# Splicesosome_Up <- c(142, 182, 6, 5)
# Translation_Up <- c(169, 55, 146, 91, 114, 102, 113)
Hypoxia_Up <- c(263, 85, 262, 15, 261, 254, 294, 4, 60, 172, 308, 181, 260, 160, 259, 148, 147, 247)
miRNA_regulation_Up <- c(319, 343, 318, 246, 142, 182, 6, 5)
Cytokine_regulation_Up <- c(21, 20, 19)
Immune_response_Up <- c(39, 43, 38, 40, 31, 30, 37, 29, 153, 204, 344, 33, 16, 217, 336, 82, 216, 168, 166, 215, 138, 251, 206, 137, 41, 34, 279, 278, 165, 32, 35, 36, 42, 167, 205, 192, 218, 225, 219, 321, 65, 324, 125, 184, 170, 227, 316, 328, 315, 237, 327, 195, 322, 332, 185)
Metabolic_pathways_Up <- c(75, 77, 106, 78, 107, 275, 76, 57, 233, 53, 276, 255, 252, 296, 171, 56, 54, 277, 145, 339, 155, 64, 62, 11, 104, 207, 98, 97, 157, 228, 193, 90, 89, 88, 197)
Differentiation_and_migration_Up <- c(122, 58, 121, 115, 116, 248, 248, 144, 123, 80, 238, 118, 119, 128, 129, 120, 342, 250, 286, 220)
Protein_Catabolism_and_localization_Up <- c(134, 334, 281, 337, 280, 229, 331, 330, 63, 81, 201)
# Growth_factors_Up <- c(155, 64, 62, 11, 104, 207, 98, 97, 157, 228, 193, 90, 89, 88, 197)
Membrane_transport_Up <- c(302, 291, 293, 307, 301, 326, 249, 208, 297, 202)
Muscle_organization_Up <- c(320, 191, 126, 232, 156, 96, 141, 68, 83, 224, 222, 92, 69, 241, 240, 299, 239, 309, 257, 341, 164, 139, 140, 211, 295, 304, 303, 25, 93, 126, 300, 230, 282, 26, 94)
# T_cell_differentiation_Up <- c(167, 205, 192, 218, 225, 219, 321, 65, 324, 125, 184, 170, 227, 316, 328, 315, 237, 327, 195, 322, 332, 185)
# Cytoskeletal_function_Up <- c(341, 164, 139, 140, 211, 295, 304, 303, 25, 93, 126, 300, 230, 282, 26, 94)
# Signaling_Up <- c(1, 95, 79, 256, 253, 173, 61, 158, 183, 306, 131, 23, 162)
# Cell_death_pathways_Up <- c(323, 333, 74, 132, 314)
# Muscle_development_Up <- c(68, 83, 224, 222, 92, 69, 241, 240, 299, 239, 309, 257)
# Wnt_Signalling_Up <- c(103, 305, 124, 234)
# Circardian_rhythm_Up <- c(194, 210, 71, 72, 209, 135)
Signaling_and_pathways_Up <- c(323, 333, 74, 132, 314, 103, 305, 124, 234, 1, 95, 79, 256, 253, 173, 61, 158, 183, 306, 131, 23, 162)

GO_BP_up <- DEG_up_GO_filtered_with_genes %>% 
  mutate(GO_Group = ifelse(row.names(.) %in% Ribosomal_activity_Up, "Ribosomal Activity", ifelse(row.names(.) %in% Hypoxia_Up, "Hypoxia",
                    ifelse(row.names(.) %in% miRNA_regulation_Up, "miRNA Regulation", ifelse(row.names(.) %in% Cytokine_regulation_Up, "Cytokine Regulation", 
                    ifelse(row.names(.) %in% Immune_response_Up, "Immune Response", ifelse(row.names(.) %in% Metabolic_pathways_Up, "Metabolic Pathways", 
                    ifelse(row.names(.) %in% Differentiation_and_migration_Up, "Differentiation & Migration", ifelse(row.names(.) %in% Protein_Catabolism_and_localization_Up, "Protein Catabolism & Localization", 
                    ifelse(row.names(.) %in% Membrane_transport_Up, "Membrane Transport", ifelse(row.names(.) %in% Muscle_organization_Up, "Muscle Organization", 
                    ifelse(row.names(.) %in% Signaling_and_pathways_Up, "Signaling Pathways", NA)))))))))))) %>% filter(!is.na(GO_Group)) %>% mutate(regulation = "IHD Up-regulated")

## DEG DOWN
Heart_development_Down <- c(100, 99, 70, 69, 16, 46, 14, 152, 68, 40, 15, 22, 18, 163, 42, 162, 17, 43, 41, 190, 159, 21, 20, 19, 101, 158, 158, 160, 77, 65, 61, 62, 166, 47, 171, 48, 167, 45, 33, 126, 123, 169, 185, 180, 182, 178, 173, 170, 183, 184, 181, 172, 177, 179, 44, 192, 168, 157, 156, 129, 153, 153, 23, 131, 164, 34, 155, 116, 105, 110, 215, 111)
# Development_and_morphogenesis_Down <- c(61, 62, 166, 47, 171, 48, 167, 45, 33, 126, 123, 169, 185, 180, 182, 178, 173, 170, 183, 184, 181, 172, 177, 179, 44, 192, 168, 157, 156, 129, 153, 153, 23, 131, 164, 34, 155, 116, 105, 110, 215, 111)
RNA_biogenesis_Down <- c(95, 35, 84, 127, 25, 124, 4, 3, 5, 2)
Cellular_transport_Down <- c(219, 218, 67, 36, 37, 217, 64, 222, 133, 211, 130, 28, 31, 30, 72, 79, 73, 212, 71, 80, 24, 76, 39, 148, 29, 75, 74, 38, 209, 112, 82, 221, 220, 216, 147, 81, 78)
Metabolic_pathways_Down <- c(58, 57, 114, 115, 90, 144, 188, 187, 60, 59, 189, 56)
Differentiation_Down <- c(143, 164, 141, 142)
Cell_cycle_Down <- c(138, 223, 213, 224, 214, 208)
Signaling_Down <- c(154, 210, 121, 120, 204, 175, 202, 119, 196, 93, 91, 198, 197, 176, 193, 92, 86, 85, 203, 89, 51, 55, 87, 54, 106, 32, 102, 52, 200)
Ischemia_Down <- c(13)
Heart_organization_Down <- c(139, 201, 199, 132, 117, 194, 125, 195, 113, 128, 107, 97, 96, 108, 98, 10)

GO_BP_down <- DEG_down_GO_filtered_with_genes %>% 
  mutate(GO_Group = ifelse(row.names(.) %in% Heart_development_Down, "Heart Development", 
                                  ifelse(row.names(.) %in% RNA_biogenesis_Down, "RNA Biogenesis",ifelse(row.names(.) %in% Cellular_transport_Down, "Cellular Transport",
                                                ifelse(row.names(.) %in% Metabolic_pathways_Down, "Metabolic Pathways", ifelse(row.names(.) %in% Differentiation_Down, "Differentiation", 
                                                              ifelse(row.names(.) %in% Differentiation_Down, "Differentiation", ifelse(row.names(.) %in% Cell_cycle_Down, "Cell Cycle", 
                                                                            ifelse(row.names(.) %in% Signaling_Down, "Signaling", ifelse(row.names(.) %in% Ischemia_Down, "Ischemia", 
                                                                                          ifelse(row.names(.) %in% Heart_organization_Down, "Heart Organization", NA))))))))))) %>% filter(!is.na(GO_Group)) %>% mutate(regulation = "IHD Down-regulated")

GO_BP_all <- rbind(GO_BP_up, GO_BP_down) 
GO_BP_all <- GO_BP_all %>% mutate(term = ifelse(term == "positive regulation of protein catabolic process in the vacuole", "positive regulation of protein catabolic processes", 
                                                ifelse(term == "regulation of T cell activation via T cell receptor contact with antigen bound to MHC molecule on antigen presenting cell", "regulation of T-cell activation",
                                                    ifelse(term == "regulation of extrinsic apoptotic signaling pathway via death domain receptors", "regulation of apoptotic signaling pathway via DDR", 
                                                           ifelse(term == "gamma-aminobutyric acid biosynthetic process", "aminobutyric acid biosynthetic process", 
                                                                  ifelse(term == "regulation of toll-like receptor 9 signaling pathway", "regulation of TLR 9 signaling pathway", 
                                                                         ifelse(term == "regulation of cardiac muscle contraction by calcium ion signaling", " regulation of cardiac muscle contraction", term)))))))

GO_BP_all_labels <- c("regulation of transepithelial transport", "positive regulation of protein catabolic processes", 
                      "regulation of T-cell activation", "skeletal muscle cell differentiation", 
                      "regulation of apoptotic signaling pathway via DDR", "cellular response to reactive oxygen species", 
                      "amino acid:sodium symporter activity", "fatty acid transmembrane transport", 
                      "cardiac myofibril", "aminobutyric acid biosynthetic process", "contractile fiber",  
                      "regulation of TLR 9 signaling pathway", "skeletal muscle fiber development", " regulation of cardiac muscle contraction")


write.csv(GO_BP_all, paste0(result_dir, "/GO_BP_all.csv"))

pdf(paste0(result_dir, "/DEG_GO_volcano_Up.pdf"), width = 15, height = 10)
  ggplot(GO_BP_up, aes(x = hitsPerc, y = -log10(over_represented_pvalue), color = GO_Group)) + 
    geom_point(size = 5) + geom_hline(yintercept = -log10(0.05), lty = "dashed", lwd = 0.5, alpha = 0.8) + geom_vline(xintercept = 1, lty = "dashed", lwd = 0.5, alpha = 0.8) + 
    # geom_text_repel(data = subset(GO_BP_all, term %in% GO_BP_all_labels), aes(label = term), size = 8, show.legend = F) + 
    geom_text_repel(data = subset(GO_BP_all, term %in% GO_BP_all_labels & Regulation == "UP"), aes(label = term), nudge_x = c(35, 21, -20, 32, 0, 10), nudge_y = c(0.1, 0.2, 0.2, 0, 0.2, -0.2), size = 8, show.legend = F) +
    scale_color_manual(values = nature_colors[1:11]) + labs(x = "Hits Percentage (%)", y = "-log10(p-value)", color = "GO Group", title = "") + 
    facet_wrap(~ regulation, scales = "fixed") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + xlim(0, 68) + ylim(1.3, 6.3) + 
    theme_linedraw() + ggtitle(NULL) + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                             panel.border = element_rect(size = 0.5), text = element_text(size = 24), strip.text = element_text(size = 30))
dev.off()

pdf(paste0(result_dir, "/DEG_GO_volcano_Down.pdf"), width = 13, height = 10)
  ggplot(GO_BP_down, aes(x = hitsPerc, y = -log10(over_represented_pvalue), color = GO_Group)) + 
    geom_point(size = 5) + geom_hline(yintercept = -log10(0.05), lty = "dashed", lwd = 0.5, alpha = 0.8) + geom_vline(xintercept = 0.80, lty = "dashed", lwd = 0.5, alpha = 0.8) +
    geom_text_repel(data = subset(GO_BP_all, term %in% GO_BP_all_labels  & Regulation == "DOWN"), aes(label = term), 
                    nudge_x = c(25,5,30,23,12,15,0,25), nudge_y = c(-0.1,0.2,-0.2,-0.2,0,0.61,0.2,0), size = 8, show.legend = F) + 
    scale_color_manual(values = nature_colors[1:10]) + labs(x = "Hits Percentage (%)", y = "-log10(p-value)", color = "GO Group", title = "") +
    facet_wrap(~ regulation, scales = "fixed") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + xlim(0, 68) + ylim(1.3, 6.3) + 
    theme_linedraw() + ggtitle(NULL) + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                             panel.border = element_rect(size = 0.5), text = element_text(size = 24), strip.text = element_text(size = 30))
dev.off()

