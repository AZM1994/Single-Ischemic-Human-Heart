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

setwd("/Users/zhemingan/Documents/BCH_research/SNV_enrichment_analysis/SNV_GO_Enrichment/heart_PTA_Cases")
figure_save_dir <- "GOseq_results/deleterious_mutation_noFDR"
dir.create(figure_save_dir, recursive = T)
numDEInCat_threshold = 2
numInCat_threshold = 1000
# gene_length_type = "Gene_length"
gene_length_type = "Exon_length"
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])

### read gene list with metadata
Hypoxia_PTA_Cases_metadata <- readRDS("SCAN2_df.rds") %>% as.data.frame() |> base::`[`(c("cell_ID", "age", "gender", "Case_ID", "condition", "snv.burden", "snv.rate.per.gb")) %>% 
  rename_with(~ c("Cell_ID", "Condition"), .cols = c(1, 5))
Condition_list <- unique(Hypoxia_PTA_Cases_metadata$Condition)
selected_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "Cell_ID", "Case_ID", "Condition", "mut_type", "age")
rename_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Type", "Gene_symbol", "ExonicFunc.refGene", "Cell_ID", "Case_ID", "Condition", "mut_type", "age")
genic_region <- c("exonic", "exonic;splicing", "intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3")
deleterious_mutation <- c("splicing", "exonic;splicing", "frameshift deletion", "frameshift insertion", "nonframeshift deletion", "nonsynonymous SNV", "stopgain")

### read RNAseq CM data
Seurat.obj_sub_clustering_CM_only <- readRDS("Seurat.obj_sub_clustering_CM_only.RDS")
expr_data <- as.matrix(GetAssayData(object = Seurat.obj_sub_clustering_CM_only, assay = "RNA", slot = "count"))
expressed_genes <- row.names(expr_data)

gene_expression_percentage <- Matrix::rowSums(expr_data > 0) / ncol(expr_data) * 100
genes_expressed_15_percent <- names(gene_expression_percentage[gene_expression_percentage > 15])

genomic_SCAN2_df <- c()
for (condition_temp in Condition_list) {
  for (mutation_type in c("ssnv", "sindel")) {
  # for (mutation_type in c("ssnv")) {
    cat("Get genomic context for", condition_temp, mutation_type, "...\n")
    heart_PTA_Cases_vcf_temp <- read.table(paste0("heart_PTA_Cases_annovar/heart_PTA_Cases.all_age.", condition_temp, "_", mutation_type, ".vcf"), sep = "\t") %>% mutate(V8 = sub(";.*", "", V8))
    genomic_context_temp <- read.csv(paste0("heart_PTA_Cases_annovar/heart_PTA_Cases.all_age.", condition_temp, "_", mutation_type, ".csv"), header = TRUE) %>%
      mutate(Cell_ID = heart_PTA_Cases_vcf_temp$V8) %>% mutate(Case_ID = str_extract(Cell_ID, "[^_]+")) %>% 
      mutate(Condition = condition_temp) %>% mutate(mut_type = mutation_type)
    
    genomic_context_temp <- genomic_context_temp %>%
      mutate(Cell_ID = ifelse(Cell_ID == "1864_M-E3-2n", "1864_E3", ifelse(Cell_ID == "4402_1_A1-2n", "4402_A1", ifelse(Cell_ID == "4402_1_A3-2n", "4402_A3",
      ifelse(Cell_ID == "4638_1_A1-2n", "4638_A1", ifelse(Cell_ID == "4638_1_A4-2n", "4638_A4", ifelse(Cell_ID == "4638_1_B2-2n", "4638_B2",
      ifelse(Cell_ID == "5657_M-D1-2n", "5657_D1", ifelse(Cell_ID == "5657_M-H1-2n", "5657_H1", ifelse(Cell_ID == "5828_CM_C2_2n", "5828_C2",
      ifelse(Cell_ID == "5828_CM_G2_2n", "5828_G2", ifelse(Cell_ID == "5919_1_C4-2n", "5919_C4", ifelse(Cell_ID == "5919_1_D2-2n", "5919_D2",
      ifelse(Cell_ID == "5919_1_E3-2n", "5919_E3", ifelse(Cell_ID == "5919_1_F6-2n", "5919_F6", ifelse(Cell_ID == "6032_1_A1-2n", "6032_A1",
      ifelse(Cell_ID == "6032_M-C7-2n", "6032_C7", ifelse(Cell_ID == "6032_M-E7-2n", "6032_E7", ifelse(Cell_ID == "1673_CM_A2_2n", "1673_A2",
      ifelse(Cell_ID == "1673_CM_A3_2n", "1673_A3", ifelse(Cell_ID == "1673_CM_D2_2n", "1673_D2", Cell_ID)))))))))))))))))))))
    
    genomic_SCAN2_df_temp <- merge(genomic_context_temp, Hypoxia_PTA_Cases_metadata[c("Cell_ID", "Case_ID", "Condition", "age")]) |> 
      base::`[`(selected_colnames) %>% rename_with(~rename_colnames) %>% 
      # filter(Type %in% genic_region) %>%
      # filter(age >= 40 & age < 80) %>% 
      filter(Type %in% deleterious_mutation[1:2] | ExonicFunc.refGene %in% deleterious_mutation[3:7]) %>% 
      mutate(Gene_symbol = str_remove(Gene_symbol, "\\(.*\\)$")) %>% 
      filter(!grepl(";", Gene_symbol))
    
    genomic_SCAN2_df <- rbind(genomic_SCAN2_df, genomic_SCAN2_df_temp)
  }
}

genomic_SCAN2_df <- genomic_SCAN2_df %>% filter(Gene_symbol %in% genes_expressed_15_percent)

genomic_context_normal <- genomic_SCAN2_df[genomic_SCAN2_df$Condition == "Normal", ] %>% filter(!duplicated(Gene_symbol))
genomic_context_disease <- genomic_SCAN2_df[genomic_SCAN2_df$Condition == "Disease", ] %>% filter(!duplicated(Gene_symbol))

##### read gene length
gene_length <- read.delim("hg19_refGene.length.tsv", header = F)
colnames(gene_length) <- c("Gene_symbol", "Transcript", "Gene_length", "Exon_length")
gene_length_deduped <- gene_length[!duplicated(gene_length$Gene_symbol),]

################################################################################
################################################################################
##### GOseq and filtering for Control
mut_gene_normal <- gene_length_deduped$Gene_symbol %in% genomic_context_normal$Gene_symbol
names(mut_gene_normal) <- gene_length_deduped$Gene_symbol
pwf_normal <- nullp(mut_gene_normal, "hg19", bias.data = gene_length_deduped[, gene_length_type])
mut_gene_normal_GO <- goseq(pwf_normal, "hg19", "geneSymbol") %>% mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% mutate(Condition = "Control")
# mut_gene_normal_GO_filtered <- mut_gene_normal_GO %>%
#   filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
#   filter(over_represented_pvalue < 0.05) %>%
#   mutate(Condition = "Control")

## find gene names in each GO term: Control
go_term_list_normal <- mut_gene_normal_GO$category  # GO term IDs from goseq
go_genes_list_normal <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list_normal, keytype = "GOALL", columns = c("SYMBOL"))
genes_in_GO_normal <- go_genes_list_normal[go_genes_list_normal$SYMBOL %in% genomic_context_normal$Gene_symbol, ]
genes_in_GO_normal_collapsed <- genes_in_GO_normal %>% group_by(GOALL) %>% 
  filter(!duplicated(SYMBOL)) %>% 
  summarise(genes = paste(SYMBOL, collapse = ", "))
mut_gene_normal_GO_with_genes <- merge(mut_gene_normal_GO, genes_in_GO_normal_collapsed, by.x = "category", by.y = "GOALL") 
mut_gene_normal_GO_filtered_with_genes <- mut_gene_normal_GO_with_genes %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
  filter(over_represented_pvalue < 0.05)

# write.csv(mut_gene_normal_GO_with_genes, paste0(figure_save_dir, "/mut_gene_normal_GO_with_genes.csv"))
write.csv(mut_gene_normal_GO_filtered_with_genes, paste0(figure_save_dir, "/mut_gene_normal_GO_filtered_with_genes.csv"))

################################################################################
################################################################################
##### GOseq and filtering for IHD
mut_gene_disease <- gene_length_deduped$Gene_symbol %in% genomic_context_disease$Gene_symbol
names(mut_gene_disease) <- gene_length_deduped$Gene_symbol
pwf_disease <- nullp(mut_gene_disease, bias.data = gene_length_deduped[, gene_length_type])
mut_gene_disease_GO <- goseq(pwf_disease, "hg19", "geneSymbol") %>% mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% mutate(Condition = "IHD")
# mut_gene_disease_GO_filtered <- mut_gene_disease_GO %>%
#   filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
#   filter(over_represented_pvalue < 0.05) %>%
#   mutate(Condition = "IHD")

## find gene names in each GO term: IHD
go_term_list_disease <- mut_gene_disease_GO$category  # GO term IDs from goseq
go_genes_list_disease <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list_disease, keytype = "GOALL", columns = c("SYMBOL"))
genes_in_GO_disease <- go_genes_list_disease[go_genes_list_disease$SYMBOL %in% genomic_context_disease$Gene_symbol, ]
genes_in_GO_disease_collapsed <- genes_in_GO_disease %>% group_by(GOALL) %>% 
  filter(!duplicated(SYMBOL)) %>% 
  summarise(genes = paste(SYMBOL, collapse = ", "))
mut_gene_disease_GO_with_genes <- merge(mut_gene_disease_GO, genes_in_GO_disease_collapsed, by.x = "category", by.y = "GOALL")
mut_gene_disease_GO_filtered_with_genes <- mut_gene_disease_GO_with_genes %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>% 
  filter(over_represented_pvalue < 0.05)

# write.csv(mut_gene_disease_GO_with_genes, paste0(figure_save_dir, "/mut_gene_disease_GO_with_genes.csv"))
write.csv(mut_gene_disease_GO_filtered_with_genes, paste0(figure_save_dir, "/mut_gene_disease_GO_filtered_with_genes.csv"))

################################################################################
################################################################################
### generate volcano plots for BP
# nature_colors <- c("#a6761d", "#228C68", "#1f77b4", "#ff7f0e", "#e377c2", "#F22020")
# nature_colors <- c("#4d79a6", "#a1cde6", "#e05759", "#ff9d99", "#f1cd63", "#79706d","#489893","#8cd17d","#00afba","#e6b800","#fb4d07")
nature_colors <- c("#a2d2e7", "#67a8cd", "#ffc17f", "#cf9f88", "#6fb3a8", "#b3e19b","#50aa4b","#ff9d9f","#f36569","#3581b7","#cdb6da",
                   "#704ba3", "#9a7fbd", "#dba9a8", "#e43030", "#e99b78", "#ff8831")
## in Control
Heart_Muscle_development_Ctrl <- c(19, 54, 17, 52, 48, 51, 42, 44, 2, 22, 16, 4, 3, 53, 18)
Contractility_Ctrl <- c(50, 13, 49, 12, 9)
Cytoskeletal_organization_Ctrl <- c(28, 59, 58, 35, 8, 40, 25, 38, 27, 31, 45, 26, 7, 29, 1, 36, 33, 60, 46, 47, 37)
mRNA_binding_Ctrl <- c(6, 5)
Cell_junction_Ctrl <- c(63, 10, 55)
Morphogenesis_Ctrl <- c(43, 20)
Tranportation_Ctrl <- c(23, 64, 62)

# GO_BP_Ctrl <- mut_gene_normal_GO_filtered_with_genes[mut_gene_normal_GO_filtered_with_genes$ontology == "BP", ] %>% 
GO_BP_Ctrl <- mut_gene_normal_GO_filtered_with_genes %>% 
  mutate(GO_Group = ifelse(row.names(.) %in% Heart_Muscle_development_Ctrl, "Heart Muscle Development", 
                           ifelse(row.names(.) %in% Contractility_Ctrl, "Contractility", 
                                  ifelse(row.names(.) %in% Cytoskeletal_organization_Ctrl, "Cytoskeletal Organization", 
                                         ifelse(row.names(.) %in% mRNA_binding_Ctrl, "mRNA Binding", 
                                                ifelse(row.names(.) %in% Cell_junction_Ctrl, "Cell Junction", 
                                                       ifelse(row.names(.) %in% Morphogenesis_Ctrl, "Morphogenesis", 
                                                              ifelse(row.names(.) %in% Tranportation_Ctrl, "Tranportation", NA)))))))) %>% filter(!is.na(GO_Group))

## in IHD
Cellular_respose_to_stress_IHD <- c(53, 46, 54, 18, 17)
Plasma_membrane_IHD <- c(55, 58, 43, 59)
RNA_splicing_IHD <- c(35, 6, 38, 41, 24, 7, 37, 40, 34, 57, 3, 4, 2, 21, 56, 16, 9, 5)
Metabolic_process_IHD <- c(26, 10, 25, 15, 32)
Protein_binding_IHD <- c(12, 29, 52, 22, 49, 33, 50, 51)
Golgi_vesicle_transport_IHD <- c(8, 47, 1, 11, 39)
Microtubule_IHD <- c(13, 14)
Ribonuclei_protein_IHD <- c(30, 28)
Histone_binding_IHD <- c(45, 27, 31)
Cellular_respose_to_stress_IHD <- c(53, 46, 54, 18, 17)

# GO_BP_IHD <- mut_gene_disease_GO_filtered_with_genes[mut_gene_disease_GO_filtered_with_genes$ontology == "BP", ] %>% 
GO_BP_IHD <- mut_gene_disease_GO_filtered_with_genes %>% 
  mutate(GO_Group = ifelse(row.names(.) %in% Cellular_respose_to_stress_IHD, "Cellular Respose to Stress", ifelse(row.names(.) %in% Plasma_membrane_IHD, "Plasma Membrane", 
                                  ifelse(row.names(.) %in% RNA_splicing_IHD, "RNA Splicing", ifelse(row.names(.) %in% Metabolic_process_IHD, "Metabolic Process", 
                                                ifelse(row.names(.) %in% Protein_binding_IHD, "Protein Binding", ifelse(row.names(.) %in% Golgi_vesicle_transport_IHD, "Golgi Vesicle Transport", 
                                                              ifelse(row.names(.) %in% Microtubule_IHD, "Microtubule", ifelse(row.names(.) %in% Ribonuclei_protein_IHD, "Ribonuclei Protein", 
                                                                            ifelse(row.names(.) %in% Histone_binding_IHD, "Histone Binding", 
                                                                                   ifelse(row.names(.) %in% Cellular_respose_to_stress_IHD, "Cellular Respose to Stress", NA))))))))))) %>% filter(!is.na(GO_Group))

GO_BP_all <- rbind(GO_BP_Ctrl, GO_BP_IHD)
GO_BP_all_labels <- c("muscle organ development", "cardiac muscle contraction", "alpha-actinin binding", 
                      "cytoplasmic stress granule", "positive regulation of mRNA processing", "acetylation-dependent protein binding", "Golgi vesicle transport")

write.csv(GO_BP_all, paste0(figure_save_dir, "/GO_BP_all.csv"))

p_del_GO_volcano_control <- ggplot(GO_BP_Ctrl, aes(x = hitsPerc, y = -log10(over_represented_pvalue), color = GO_Group)) + 
  geom_point(size = 6) + geom_hline(yintercept = -log10(0.05), lty = "dashed", lwd = 0.5, alpha = 0.8) + geom_vline(xintercept = 0.33, lty = "dashed", lwd = 0.5, alpha = 0.8) + 
  geom_text_repel(data = subset(GO_BP_Ctrl, term %in% GO_BP_all_labels), aes(label = term), nudge_x = c(4.5,1,4.8), nudge_y = c(-0.2, 0.2, 0.1), size = 12, show.legend = F) + 
  scale_color_manual(values = nature_colors[1:7]) + labs(x = "Hits Percentage (%)", y = "-log10(p-value)", color = "GO Group", title = "") + 
  facet_wrap(~ Condition, scales = "fixed") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + xlim(0, 13) + ylim(1.3, 5.5) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5))
ggsave(paste0(figure_save_dir, "/del_GO_volcano_control.pdf"), plot = p_del_GO_volcano_control, width = 15, height = 10, dpi = 600)

p_del_GO_volcano_IHD <- ggplot(GO_BP_IHD, aes(x = hitsPerc, y = -log10(over_represented_pvalue), color = GO_Group)) + 
  geom_point(size = 6) + geom_hline(yintercept = -log10(0.05), lty = "dashed", lwd = 0.5, alpha = 0.8) + geom_vline(xintercept = 0.33, lty = "dashed", lwd = 0.5, alpha = 0.8) + 
  geom_text_repel(data = subset(GO_BP_IHD, term %in% GO_BP_all_labels), aes(label = term), nudge_x = c(5.5,5.5,1.5,1), nudge_y = c(0.12, -0.1, 0.35, -0.2), size = 12, show.legend = F) + 
  scale_color_manual(values = nature_colors[8:16]) + labs(x = "Hits Percentage (%)", y = "-log10(p-value)", color = "GO Group", title = "") + 
  facet_wrap(~ Condition, scales = "fixed") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + xlim(0, 13) + ylim(1.3, 5.5) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5))
ggsave(paste0(figure_save_dir, "/del_GO_volcano_IHD.pdf"), plot = p_del_GO_volcano_IHD, width = 15, height = 10, dpi = 600)
