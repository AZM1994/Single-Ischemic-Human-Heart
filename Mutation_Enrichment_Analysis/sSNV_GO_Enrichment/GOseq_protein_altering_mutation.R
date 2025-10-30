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

setwd("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration/Mutation_Enrichment_Analysis/")
result_dir <- "sSNV_GO_Enrichment/GOseq/protein_altering"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
numDEInCat_threshold = 2
numInCat_threshold = 1000
gene_length_type = "Exon_length"

##### read gene list with metadata
SCAN2_df <- readRDS("data/1-sSNV_SCAN2_df_filtered.rds") %>% as.data.frame() %>% dplyr::select(Cell_ID, Age, Gender, Case_ID, Condition, snv.burden, snv.rate.per.gb)
selected_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene")
genic_region <- c("exonic", "exonic;splicing", "intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3")
# protein_altering_mutation <- c("splicing", "exonic;splicing", "frameshift deletion", "frameshift insertion", "nonframeshift deletion", "nonsynonymous SNV", "stopgain", "stoploss")
protein_altering_mutation <- c("splicing", "exonic;splicing", "nonsynonymous SNV", "stopgain", "stoploss")

genomic_SCAN2_df <- c()
genomic_SCAN2_alt_df <- c()
for (cond_temp in unique(SCAN2_df$Condition)) {
  for (mutation_type in c("ssnv")) {
    cat("Get genomic context for", cond_temp, mutation_type, "...\n")
    heart_PTA_all_temp_vcf <- read.table(paste0("data/annotation_results/all_age_", mutation_type, "/", cond_temp, "/heart_PTA_all.all_age.", cond_temp, "_", mutation_type, ".vcf"), sep = "\t")
    genomic_context_temp <- read.csv(paste0("data/annotation_results/all_age_", mutation_type, "/", cond_temp, "/heart_PTA_all.all_age.", cond_temp, "_", mutation_type, ".hg19_multianno.csv"), header = TRUE) %>% 
      mutate(has_second_gene = str_detect(Gene.refGene, ";"), 
             first_gene = str_extract(Gene.refGene, "^[^;]+"),
             second_gene = if_else(has_second_gene, str_extract(Gene.refGene, "(?<=;)[^;]+"), NA_character_),
             first_dist = as.numeric(str_extract(GeneDetail.refGene, "(?<=dist=)[0-9]+")),
             second_dist = if_else(has_second_gene, as.numeric(str_extract(GeneDetail.refGene, "(?<=;dist=)[0-9]+")), NA_real_),
             Gene.refGene = case_when(!has_second_gene ~ Gene.refGene, first_dist < second_dist ~ first_gene, TRUE ~ second_gene)) %>% 
      dplyr::select(all_of(selected_colnames)) %>% 
      mutate(Cell_ID = heart_PTA_all_temp_vcf$V8, Case_ID = str_extract(Cell_ID, "[^_]+"), Condition = cond_temp, mut_type = mutation_type)
    
    genomic_SCAN2_df_temp <- genomic_context_temp %>% merge(SCAN2_df %>% dplyr::select(Cell_ID, Condition, Gender, Age), by = c("Cell_ID", "Condition"))
    genomic_SCAN2_df_temp_alt <- genomic_SCAN2_df_temp %>% 
      # filter(Func.refGene %in% c("splicing", "exonic;splicing") | ExonicFunc.refGene %in% c("frameshift deletion", "frameshift insertion", "nonframeshift deletion", "nonsynonymous SNV", "stopgain", "stoploss"))
      filter(Func.refGene %in% c("splicing", "exonic;splicing") | ExonicFunc.refGene %in% c("nonsynonymous SNV", "stopgain", "stoploss"))
    
    genomic_SCAN2_df <- rbind(genomic_SCAN2_df, genomic_SCAN2_df_temp)
    genomic_SCAN2_alt_df <- rbind(genomic_SCAN2_alt_df, genomic_SCAN2_df_temp_alt)
  }
}

genomic_context_ctrl <- genomic_SCAN2_alt_df[genomic_SCAN2_alt_df$Condition == "Control", ] %>% filter(!duplicated(Gene.refGene))
genomic_context_dis <- genomic_SCAN2_alt_df[genomic_SCAN2_alt_df$Condition == "IHD", ] %>% filter(!duplicated(Gene.refGene))
write.csv(genomic_context_ctrl$Gene.refGene, paste0(result_dir, "/Control_protein_altering_mutation_gene.csv"))
write.csv(genomic_context_dis$Gene.refGene, paste0(result_dir, "/IHD_protein_altering_mutation_gene.csv"))
write.csv(genomic_context_ctrl, paste0(result_dir, "/Control_protein_altering_sSNV_annotation.csv"))
write.csv(genomic_context_dis, paste0(result_dir, "/IHD_protein_altering_sSNV_annotation.csv"))

##### read gene length
gene_length <- read.delim("data/hg19_refGene.length.tsv", header = F)
colnames(gene_length) <- c("Gene.refGene", "Transcript", "Gene_length", "Exon_length")
gene_length_deduped <- gene_length[!duplicated(gene_length$Gene.refGene), ]

################################################################################
##### GOseq and filtering for Control
mut_gene_ctrl <- gene_length_deduped$Gene.refGene %in% genomic_context_ctrl$Gene.refGene
names(mut_gene_ctrl) <- gene_length_deduped$Gene.refGene
pwf_ctrl <- nullp(mut_gene_ctrl, "hg19", bias.data = gene_length_deduped[, gene_length_type])
mut_gene_ctrl_GO <- goseq(pwf_ctrl, "hg19", "geneSymbol") %>% 
  mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% 
  mutate(FDR = p.adjust(over_represented_pvalue, method = "fdr")) %>% 
  mutate(Condition = "Control")

##### find gene names in each GO term: Control
go_term_list_ctrl <- mut_gene_ctrl_GO$category  # GO term IDs from goseq
go_genes_list_ctrl <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list_ctrl, keytype = "GOALL", columns = c("SYMBOL"))
genes_in_GO_ctrl <- go_genes_list_ctrl[go_genes_list_ctrl$SYMBOL %in% genomic_context_ctrl$Gene.refGene, ]
genes_in_GO_ctrl_collapsed <- genes_in_GO_ctrl %>% group_by(GOALL) %>% 
  filter(!duplicated(SYMBOL)) %>% 
  summarise(genes = paste(SYMBOL, collapse = ", "))
mut_gene_ctrl_GO_with_genes <- merge(mut_gene_ctrl_GO, genes_in_GO_ctrl_collapsed, by.x = "category", by.y = "GOALL") 
mut_gene_ctrl_GO_filtered_with_genes <- mut_gene_ctrl_GO_with_genes %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
  filter(over_represented_pvalue < 0.05)
  # filter(FDR < 0.05)
# write.csv(mut_gene_ctrl_GO_filtered_with_genes, paste0(result_dir, "/Control_GO.protein_altering_mutation0.csv"))

################################################################################
##### GOseq and filtering for IHD
mut_gene_dis <- gene_length_deduped$Gene.refGene %in% genomic_context_dis$Gene.refGene
names(mut_gene_dis) <- gene_length_deduped$Gene.refGene
pwf_dis <- nullp(mut_gene_dis, bias.data = gene_length_deduped[, gene_length_type])
mut_gene_dis_GO <- goseq(pwf_dis, "hg19", "geneSymbol") %>% 
  mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% 
  mutate(FDR = p.adjust(over_represented_pvalue, method = "fdr")) %>% 
  mutate(Condition = "IHD")

##### find gene names in each GO term: IHD
go_term_list_dis <- mut_gene_dis_GO$category  # GO term IDs from goseq
go_genes_list_dis <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list_dis, keytype = "GOALL", columns = c("SYMBOL"))
genes_in_GO_dis <- go_genes_list_dis[go_genes_list_dis$SYMBOL %in% genomic_context_dis$Gene.refGene, ]
genes_in_GO_dis_collapsed <- genes_in_GO_dis %>% group_by(GOALL) %>% 
  filter(!duplicated(SYMBOL)) %>% 
  summarise(genes = paste(SYMBOL, collapse = ", "))
mut_gene_dis_GO_with_genes <- merge(mut_gene_dis_GO, genes_in_GO_dis_collapsed, by.x = "category", by.y = "GOALL")
mut_gene_dis_GO_filtered_with_genes <- mut_gene_dis_GO_with_genes %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>% 
  filter(over_represented_pvalue < 0.05)
  # filter(FDR < 0.05)
# write.csv(mut_gene_dis_GO_filtered_with_genes, paste0(result_dir, "/IHD_GO.protein_altering_mutation0.csv"))

################################################################################
################################################################################
### generate volcano plots for BP
warm_colors <- c("#FFD700", "#FF6F61", "#BFA2DB", "#98DD62", "#7EC8E3")
cool_colors <- c("#FF9445", "#E64545", "#B151A2", "#6A974F", "#3F7CA7")

## in Control
Heart_Muscle_Development_and_Structure_Ctrl <- c(8, 9, 14, 15, 63, 65, 71, 76, 77, 84, 100)
Cellular_Stress_and_Stimulus_Response_Ctrl <- c(5, 17, 18, 21, 24, 39, 43, 44, 53, 54, 57, 79, 80, 81, 82, 83, 87, 88, 99, 103, 106)
Signal_Processing_and_Homeostasis_Ctrl <- c(1, 2, 3, 4, 6, 7, 10, 12, 13, 16, 19, 20, 22, 25, 26, 27, 28, 31, 32, 33,
                                            34, 36, 37, 38, 45, 46, 48, 49, 52, 56, 58, 60, 61, 62, 64, 66, 67, 68,
                                            69, 70, 72, 73, 74, 85, 86, 89, 93, 94, 96, 97, 98, 101, 102, 105)
Immune_Response_and_Modulation_Ctrl <- c(29, 30, 40, 51)
Transportation_and_Trafficking_Ctrl <- c(11, 23, 35, 50, 55, 90, 91, 92, 95)

mut_gene_ctrl_GO_filtered_with_genes <- mut_gene_ctrl_GO_filtered_with_genes %>% filter(hitsPerc >= 2, ontology == "BP")
GO_BP_Ctrl <- mut_gene_ctrl_GO_filtered_with_genes %>% 
  mutate(GO_Group = case_when(row.names(.) %in% Heart_Muscle_Development_and_Structure_Ctrl ~ "Heart Muscle Development", 
                              row.names(.) %in% Cellular_Stress_and_Stimulus_Response_Ctrl ~ "Cellular Stress & Stimulus Response", 
                              row.names(.) %in% Signal_Processing_and_Homeostasis_Ctrl ~ "Signal Processing and Sensory Transduction", 
                              row.names(.) %in% Immune_Response_and_Modulation_Ctrl ~ "Immune Response & Modulation", 
                              row.names(.) %in% Transportation_and_Trafficking_Ctrl ~ "Transportation", 
                              TRUE ~ NA_character_)) %>% filter(!is.na(GO_Group), hitsPerc >= 2, ontology == "BP") %>% arrange(GO_Group)

## in IHD
Metabolic_and_RNA_Process_IHD <- c(1, 2, 17, 18, 19, 20, 37, 38, 39, 40, 67, 68, 75, 84, 147, 148, 150, 151, 169, 170)
Cell_Signaling_and_Gene_Expression_IHD <- c(4, 5, 42, 43, 48, 49, 53, 60, 61, 62,
                                            63, 85, 86, 87, 88, 89, 90, 93, 94, 111, 112,
                                            113, 114, 116, 117, 118, 119, 120, 121, 129,
                                            130, 131, 132, 133, 134, 135, 138, 139, 140,
                                            141, 153, 158, 159, 160, 161, 163, 171)
Cellular_Stress_and_Apoptosis_IHD <- c(28, 29, 54, 58, 59, 65, 66, 72, 76, 122, 123, 136, 150, 151, 152)
Immune_and_Inflammatory_Modulation_IHD <- c(6, 7, 8, 9, 10, 30, 32, 33, 34, 50, 51, 52,
                                            69, 70, 91, 92, 124, 125, 126, 168)
Genome_Stability_and_Chromosomal_Regulation <- c(3, 12, 13, 14, 15, 16, 21, 22, 23, 31, 35,
                                                 36, 41, 44, 45, 55, 56, 57, 64, 71, 77,
                                                 78, 79, 81, 82, 95, 96, 97, 98, 104, 105, 106, 164)

mut_gene_dis_GO_filtered_with_genes <- mut_gene_dis_GO_filtered_with_genes %>% filter(hitsPerc >= 2, ontology == "BP")
GO_BP_IHD <- mut_gene_dis_GO_filtered_with_genes %>% 
  mutate(GO_Group = case_when(row.names(.) %in% Metabolic_and_RNA_Process_IHD ~ "Metabolic Process & RNA Metabolism", 
                              row.names(.) %in% Cell_Signaling_and_Gene_Expression_IHD ~ "Cell Signaling & Gene Expression", 
                              row.names(.) %in% Cellular_Stress_and_Apoptosis_IHD ~ "Cellular Stress & Apoptosis", 
                              row.names(.) %in% Immune_and_Inflammatory_Modulation_IHD ~ "Immune & Inflammatory Modulation", 
                              row.names(.) %in% Genome_Stability_and_Chromosomal_Regulation ~ "Genome Stability & Chromosomal Regulation", 
                              TRUE ~ NA_character_)) %>% filter(!is.na(GO_Group), hitsPerc >= 2, ontology == "BP") %>% arrange(GO_Group)

GO_BP_all_labels <- c("regulation of skeletal muscle tissue development", "cellular response to prostaglandin E stimulus", 
                      "ionotropic glutamate receptor signaling pathway", "negative regulation of type I interferon production", 
                      "vesicle cytoskeletal trafficking", 
                      "fatty-acyl-CoA metabolic process", "coronary artery morphogenesis", 
                      "intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator", 
                      "positive regulation of T cell chemotaxis", "telomere maintenance via recombination")

write.csv(GO_BP_Ctrl, paste0(result_dir, "/Control_GO.protein_altering_mutation.csv"))
write.csv(GO_BP_IHD, paste0(result_dir, "/IHD_GO.protein_altering_mutation.csv"))

p_alt_GO_volcano_control <- ggplot(GO_BP_Ctrl, aes(x = hitsPerc, y = -log10(over_represented_pvalue), color = GO_Group)) + 
  geom_point(size = 5) + geom_vline(xintercept = 2, color = "gray40", linetype = "22") + geom_hline(yintercept = -log10(0.05), color = "gray40", linetype = "22") + 
  geom_text_repel(data = subset(GO_BP_Ctrl, term %in% GO_BP_all_labels), aes(label = term), 
                  nudge_x = c(0.1, 0.1, 0.1, 0.1, 0.1), nudge_y = c(0.1, 0.3, 0.1, 0.1, 0.3), size = 5, color = "black", show.legend = F) +
  scale_color_manual(values = warm_colors) + labs(x = "Hits Percentage (%)", y = "-log10(p-value)", color = "GO Group", title = "") + facet_wrap(~ Condition, scales = "fixed") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(1.3, 4.0), breaks = c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0) ) + 
  scale_x_continuous(trans = "log2", limits = c(1.97, 40), breaks = c(2, 4, 8, 16, 32)) + 
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), 
                           text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
ggsave(paste0(result_dir, "/alt_GO_volcano_control.pdf"), plot = p_alt_GO_volcano_control, width = 15, height = 8, dpi = 600)

p_alt_GO_volcano_IHD <- ggplot(GO_BP_IHD, aes(x = hitsPerc, y = -log10(over_represented_pvalue), color = GO_Group)) +
  geom_point(size = 5) + geom_vline(xintercept = 2, color = "gray40", linetype = "22") + geom_hline(yintercept = -log10(0.05), color = "gray40", linetype = "22") + 
  geom_text_repel(data = subset(GO_BP_IHD, term %in% GO_BP_all_labels), aes(label = term), 
                  nudge_x = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nudge_y = c(0.1, 0.3, 0.1, 0.1, 0.3, 0.3), size = 5, color = "black", show.legend = F) + 
  scale_color_manual(values = cool_colors) + labs(x = "Hits Percentage (%)", y = "-log10(p-value)", color = "GO Group", title = "") + facet_wrap(~ Condition, scales = "fixed") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(1.3, 4.0), breaks = c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0) ) + 
  scale_x_continuous(trans = "log2", limits = c(1.97, 40), breaks = c(2, 4, 8, 16, 32)) + 
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), 
                           text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
ggsave(paste0(result_dir, "/alt_GO_volcano_IHD.pdf"), plot = p_alt_GO_volcano_IHD, width = 15, height = 8, dpi = 600)
