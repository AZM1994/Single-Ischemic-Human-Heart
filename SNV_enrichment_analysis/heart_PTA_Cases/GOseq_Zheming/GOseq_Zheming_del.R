library(goseq)
library(ggplot2)
library(reshape2)
library(stringr)
library(biomaRt)
library(DOSE)
library(dplyr)
library(tidyr)

setwd("/Users/zhemingan/Documents/BCH_research/annovar/heart_PTA_Cases/GOseq_Zheming")
figure_save_dir <- "GO_results/deleterious_mutation_noFDR"
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
mut_gene_normal_GO <- goseq(pwf_normal, "hg19", "geneSymbol") %>% mutate(hitsPerc = numDEInCat * 100 / numInCat)
mut_gene_normal_GO_filtered <- mut_gene_normal_GO %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>% 
  filter(over_represented_pvalue < 0.05) %>% 
  mutate(Condition = "Control")

## find gene names in each GO term: Control
go_term_list_normal <- mut_gene_normal_GO_filtered$category  # GO term IDs from goseq
go_genes_list_normal <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list_normal, keytype = "GOALL", columns = c("SYMBOL"))
genes_in_GO_normal <- go_genes_list_normal[go_genes_list_normal$SYMBOL %in% genomic_context_normal$Gene_symbol, ]
genes_in_GO_normal_collapsed <- genes_in_GO_normal %>% group_by(GOALL) %>% 
  filter(!duplicated(SYMBOL)) %>% 
  summarise(genes = paste(SYMBOL, collapse = ", "))
mut_gene_normal_GO_filtered_with_genes <- merge(mut_gene_normal_GO_filtered, genes_in_GO_normal_collapsed, by.x = "category", by.y = "GOALL") 
write.csv(mut_gene_normal_GO_filtered_with_genes, paste0(figure_save_dir, "/mut_gene_normal_GO_filtered_with_genes.csv"))

################################################################################
################################################################################
##### GOseq and filtering for IHD
mut_gene_disease <- gene_length_deduped$Gene_symbol %in% genomic_context_disease$Gene_symbol
names(mut_gene_disease) <- gene_length_deduped$Gene_symbol
pwf_disease <- nullp(mut_gene_disease, bias.data = gene_length_deduped[, gene_length_type])
mut_gene_disease_GO <- goseq(pwf_disease, "hg19", "geneSymbol") %>% mutate(hitsPerc = numDEInCat * 100 / numInCat)
mut_gene_disease_GO_filtered <- mut_gene_disease_GO %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>% 
  filter(over_represented_pvalue < 0.05) %>% 
  mutate(Condition = "IHD")

## find gene names in each GO term: IHD
go_term_list_disease <- mut_gene_disease_GO_filtered$category  # GO term IDs from goseq
go_genes_list_disease <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list_disease, keytype = "GOALL", columns = c("SYMBOL"))
genes_in_GO_disease <- go_genes_list_disease[go_genes_list_disease$SYMBOL %in% genomic_context_disease$Gene_symbol, ]
genes_in_GO_disease_collapsed <- genes_in_GO_disease %>% group_by(GOALL) %>% 
  filter(!duplicated(SYMBOL)) %>% 
  summarise(genes = paste(SYMBOL, collapse = ", "))
mut_gene_disease_GO_filtered_with_genes <- merge(mut_gene_disease_GO_filtered, genes_in_GO_disease_collapsed, by.x = "category", by.y = "GOALL")
write.csv(mut_gene_disease_GO_filtered_with_genes, paste0(figure_save_dir, "/mut_gene_disease_GO_filtered_with_genes.csv"))

################################################################################
################################################################################
### generate plots for Control and IHD comparison
GO_final <- rbind(mut_gene_normal_GO_filtered_with_genes, mut_gene_disease_GO_filtered_with_genes)
GO_final_melt <- melt(GO_final[, c(2,4,6,7,8,9)], id = c("term", "ontology", "Condition", "hitsPerc", "numDEInCat"), value.name = "over_represented_pvalue")
GO_final_melt_BP <- GO_final_melt[GO_final_melt$ontology == c("BP"), ] %>% 
  mutate(Condition = factor(Condition, level = c("Control", "IHD"))) %>% 
  complete(term, ontology, Condition, fill = list(over_represented_pvalue = NA))
GO_final_melt_MF <- GO_final_melt[GO_final_melt$ontology == c("MF"), ] %>% 
  mutate(Condition = factor(Condition, level = c("Control", "IHD"))) %>% 
  complete(term, ontology, Condition, fill = list(over_represented_pvalue = NA))
GO_final_melt_CC <- GO_final_melt[GO_final_melt$ontology == c("CC"), ] %>% 
  mutate(Condition = factor(Condition, level = c("Control", "IHD"))) %>% 
  complete(term, ontology, Condition, fill = list(over_represented_pvalue = NA))

pdf(paste0(figure_save_dir, "/del_GO_barplot.pdf"), width = 20, height = 25)
p11 <- ggplot(GO_final_melt_BP, aes(x = reorder(term, -log10(over_represented_pvalue)), y = -log10(over_represented_pvalue), fill = Condition)) +
  geom_col(position = position_dodge2(reverse = TRUE), width = 0.8) + facet_wrap(ontology ~ ., scales="free_x") + scale_fill_manual(values = color_set) + 
  theme_classic() + coord_flip() + labs(title = "", x = "GO Term", y = "-log10(p-value)") + theme(text = element_text(size=20))
p11
p21 <- ggplot(GO_final_melt_MF, aes(x = reorder(term, -log10(over_represented_pvalue)), y = -log10(over_represented_pvalue), fill = Condition)) +
  geom_col(position = position_dodge2(reverse = TRUE), width = 0.8) + facet_wrap(ontology ~ ., scales="free_x") + scale_fill_manual(values = color_set) + 
  theme_classic() + coord_flip() + labs(title = "", x = "GO Term", y = "-log10(p-value)") + theme(text = element_text(size=20))
p21
p31 <- ggplot(GO_final_melt_CC, aes(x = reorder(term, -log10(over_represented_pvalue)), y = -log10(over_represented_pvalue), fill = Condition)) +
  geom_col(position = position_dodge2(reverse = TRUE), width = 0.8) + facet_wrap(ontology ~ ., scales="free_x") + scale_fill_manual(values = color_set) + 
  theme_classic() + coord_flip() + labs(title = "", x = "GO Term", y = "-log10(p-value)") + theme(text = element_text(size=20))
p31
dev.off()

pdf(paste0(figure_save_dir, "/del_GO_dotplot.pdf"), width = 35, height = 25)
p12 <- ggplot(GO_final_melt_BP, aes(x = reorder(term, -log10(over_represented_pvalue)), y = hitsPerc, size = numDEInCat, color = over_represented_pvalue)) + 
  geom_point(size = 8) + facet_wrap(~ Condition, scales = "free") + coord_flip() + scale_color_gradient(low = "red", high = "blue") + 
  labs(title = "BP", x = "GO Term", y = "Hits (%)", size = "Gene count", color = "p-value") + theme(text = element_text(size=20))
p12

p22 <- ggplot(GO_final_melt_MF, aes(x = reorder(term, -log10(over_represented_pvalue)), y = hitsPerc, size = numDEInCat, color = over_represented_pvalue)) + 
  geom_point(size = 8) + facet_wrap(~ Condition, scales = "free") + coord_flip() + scale_color_gradient(low = "red", high = "blue") + 
  labs(title = "MF", x = "GO Term", y = "Hits (%)", size = "Gene count", color = "p-value") + theme(text = element_text(size=20))
p22
p32 <- ggplot(GO_final_melt_CC, aes(x = reorder(term, -log10(over_represented_pvalue)), y = hitsPerc, size = numDEInCat, color = over_represented_pvalue)) + 
  geom_point(size = 8) + facet_wrap(~ Condition, scales = "free") + coord_flip() + scale_color_gradient(low = "red", high = "blue") + 
  labs(title = "CC", x = "GO Term", y = "Hits (%)", size = "Gene count", color = "p-value") + theme(text = element_text(size=20))
p32
dev.off()

pdf(paste0(figure_save_dir, "/del_GO_volcano.pdf"), width = 10, height = 10)
GO_final_melt_BP <- GO_final_melt_BP %>% filter(!is.na(numDEInCat)) %>% 
  mutate(color_scale = numDEInCat / max(numDEInCat[!is.na(numDEInCat)]) + (-log10(over_represented_pvalue) / max(-log10(over_represented_pvalue))))
p13 <- ggplot(GO_final_melt_BP, aes(x = numDEInCat, y = -log10(over_represented_pvalue), color = color_scale)) +
  geom_point(size = 6) + scale_color_gradient(low = "blue", high = "red") + 
  theme_classic() + labs(x = "Gene Count", y = "-log10(p-value)", color = "-log10(p-value)", title = "") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
p13

GO_final_melt_MF <- GO_final_melt_MF %>% filter(!is.na(numDEInCat)) %>% 
  mutate(color_scale = numDEInCat / max(numDEInCat[!is.na(numDEInCat)]) + (-log10(over_represented_pvalue) / max(-log10(over_represented_pvalue))))
p23 <- ggplot(GO_final_melt_MF, aes(x = numDEInCat, y = -log10(over_represented_pvalue), color = -log10(over_represented_pvalue))) +
  geom_point(size = 6) + scale_color_gradient(low = "blue", high = "red") + 
  theme_classic() + labs(x = "Gene Count", y = "-log10(p-value)", color = "-log10(p-value)", title = "") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
p23

GO_final_melt_CC <- GO_final_melt_CC %>% filter(!is.na(numDEInCat)) %>% 
  mutate(color_scale = numDEInCat / max(numDEInCat[!is.na(numDEInCat)]) + (-log10(over_represented_pvalue) / max(-log10(over_represented_pvalue))))
p33 <- ggplot(GO_final_melt_CC, aes(x = numDEInCat, y = -log10(over_represented_pvalue), color = -log10(over_represented_pvalue))) +
  geom_point(size = 6) + scale_color_gradient(low = "blue", high = "red") + 
  theme_classic() + labs(x = "Gene Count", y = "-log10(p-value)", color = "-log10(p-value)", title = "") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
p33
dev.off()