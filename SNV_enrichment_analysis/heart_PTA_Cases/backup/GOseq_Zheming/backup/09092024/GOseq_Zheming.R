library(goseq)
library(ggplot2)
library(reshape2)
library(stringr)
library(biomaRt)
library(DOSE)

setwd("/Users/zhemingan/Documents/BCH_research/annovar/heart_PTA_Cases/GOseq_Zheming")
numDEInCat_threshold = 2
numInCat_threshold = 1000
gene_length_type = "Gene_length"
# gene_length_type = "Exon_length"

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
      # filter(Type %in% c("exonic")) %>%
      # filter(Type %in% genic_region) %>%
      # filter(age >= 40 & age < 80) %>% 
      # filter(Func.refGene %in% deleterious_mutation[1:2] | ExonicFunc.refGene %in% deleterious_mutation[3:7])
      mutate(Gene_symbol = str_remove(Gene_symbol, "\\(.*\\)$"))
      # filter(!grepl(";", Gene_symbol))
    
    genomic_SCAN2_df <- rbind(genomic_SCAN2_df, genomic_SCAN2_df_temp)
  }
}

genomic_context_normal <- genomic_SCAN2_df[genomic_SCAN2_df$Condition == "Normal", ] %>% filter(!duplicated(Gene_symbol))
genomic_context_disease <- genomic_SCAN2_df[genomic_SCAN2_df$Condition == "Disease", ] %>% filter(!duplicated(Gene_symbol))

##### read gene length
gene_length <- read.delim("hg19_refGene.length.tsv", header = F)
colnames(gene_length) <- c("Gene_symbol", "Transcript", "Gene_length", "Exon_length")
gene_length_deduped <- gene_length[!duplicated(gene_length$Gene_symbol),]

# a <- getlength(genomic_context_normal$Gene_symbol, "hg19", "geneSymbol")
# hg19.geneSymbol.LENGTH <- hg19.geneSymbol.LENGTH %>% filter(!duplicated(Gene))
##### GOseq and filtering for Normal
# all_gene_normal <- hg19.geneSymbol.LENGTH$Gene %in% genomic_context_normal$Gene_symbol
# names(all_gene_normal) <- hg19.geneSymbol.LENGTH$Gene
# pwf <- nullp(all_gene_normal, "hg19", bias.data = hg19.geneSymbol.LENGTH$Length)
all_gene_normal <- gene_length_deduped$Gene_symbol %in% genomic_context_normal$Gene_symbol
names(all_gene_normal) <- gene_length_deduped$Gene_symbol
pwf <- nullp(all_gene_normal, "hg19", bias.data = gene_length_deduped[, gene_length_type])
# pwf <- nullp(all_gene_normal, "hg19", "geneSymbol")
all_gene_normal_GO <- goseq(pwf, "hg19", "geneSymbol") %>% mutate(hitsPerc = numDEInCat * 100 / numInCat)
all_gene_normal_GO_filtered <- all_gene_normal_GO %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
  mutate(Condition = "Control") %>% 
  mutate(over_represented_pvalue_corrected = p.adjust(over_represented_pvalue, method = "fdr")) %>% 
  filter(over_represented_pvalue_corrected < 0.05)

### find gene names in each GO term
go_term_list <- all_gene_normal_GO_filtered$category  # GO term IDs from goseq
go_genes_list <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list, keytype = "GO", columns = c("SYMBOL"))
go_genes_all_gene_normal <- go_genes_list[go_genes_list$SYMBOL %in% row.names(as.data.frame(all_gene_normal)), ]
go_genes_collapsed <- go_genes_all_gene_normal %>% group_by(GO) %>% summarise(genes = paste(SYMBOL, collapse = ", "))
all_gene_normal_GO_filtered_with_genes <- merge(all_gene_normal_GO_filtered, go_genes_collapsed, by.x = "category", by.y = "GO") 
# all_gene_normal_GO_filtered_with_genes <- merge(all_gene_normal_GO_filtered, go_genes_collapsed, by.x = "category", by.y = "GO") %>% 
#   group_by(ontology) %>% slice_max(order_by = -over_represented_pvalue, n = 30)

write.csv(all_gene_normal_GO_filtered_with_genes, "GO_results/all_gene_normal_GO_filtered_with_genes.csv")

pdf("GO_results/all_gene_normal_GO.pdf", width = 20, height = 45)
ggplot(all_gene_normal_GO_filtered_with_genes, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = -log10(over_represented_pvalue_corrected), fill = over_represented_pvalue_corrected)) + 
  geom_bar(stat = "identity") + facet_wrap(~ ontology, scales = "free_y", ncol = 1) + 
  coord_flip() + scale_fill_gradient(low = "red", high = "blue") + 
  labs(x = "GO Term", y = "-log10(pvalue_adjusted)", fill = "pvalue_adjusted") + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1)) 

ggplot(all_gene_normal_GO_filtered_with_genes, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = hitsPerc, size = numDEInCat, color = over_represented_pvalue_corrected)) + 
  geom_point() + 
  facet_wrap(~ ontology, scales = "free_y", ncol = 1) + 
  coord_flip() + scale_color_gradient(low = "red", high = "blue") + 
  labs(title = "", x = "GO Term", y = "Hits (%)", size = "Gene count", color = "adjusted pvalue") + 
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5))

p31 <- ggplot(all_gene_normal_GO_filtered_with_genes, aes(x = numDEInCat, y = -log10(over_represented_pvalue_corrected), color = over_represented_pvalue_corrected)) +
  geom_point(size = 6) +  # Adjust size for better visibility
  scale_color_gradient(low = "red", high = "blue") +  # Color gradient for p-values
  # theme_classic() +
  labs(x = "Gene Count", y = "-log10(pvalue_adjusted)", color = "pvalue_adjusted", title = "") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("GO_results/all_gene_normal_GO_02.pdf", plot = p31, width = 12, height = 10, dpi = 300)
dev.off()

##### GOseq and filtering for IHD
# all_gene_disease <- hg19.geneSymbol.LENGTH$Gene %in% genomic_context_disease$Gene_symbol
# names(all_gene_disease) <- hg19.geneSymbol.LENGTH$Gene
# pwf <- nullp(all_gene_disease, "hg19", bias.data = hg19.geneSymbol.LENGTH$Length)
all_gene_disease <- gene_length_deduped$Gene_symbol %in% genomic_context_disease$Gene_symbol
names(all_gene_disease) <- gene_length_deduped$Gene_symbol
pwf <- nullp(all_gene_disease, bias.data = gene_length_deduped[, gene_length_type])
all_gene_disease_GO <- goseq(pwf, "hg19", "geneSymbol") %>% mutate(hitsPerc = numDEInCat * 100 / numInCat)
all_gene_disease_GO_filtered <- all_gene_disease_GO %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
  mutate(Condition = "IHD") %>% 
  mutate(over_represented_pvalue_corrected = p.adjust(over_represented_pvalue, method = "fdr")) %>% 
  filter(over_represented_pvalue_corrected < 0.05)

### find gene names in each GO term
go_term_list <- all_gene_disease_GO_filtered$category  # GO term IDs from goseq
go_genes_list <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list, keytype = "GO", columns = c("SYMBOL"))
go_genes_all_gene_disease <- go_genes_list[go_genes_list$SYMBOL %in% row.names(as.data.frame(all_gene_disease)), ]
go_genes_collapsed <- go_genes_all_gene_disease %>% group_by(GO) %>% summarise(genes = paste(SYMBOL, collapse = ", "))
all_gene_disease_GO_filtered_with_genes <- merge(all_gene_disease_GO_filtered, go_genes_collapsed, by.x = "category", by.y = "GO")
# all_gene_disease_GO_filtered_with_genes <- merge(all_gene_disease_GO_filtered, go_genes_collapsed, by.x = "category", by.y = "GO") %>% 
#   group_by(ontology) %>% slice_max(order_by = -over_represented_pvalue, n = 30)

write.csv(all_gene_disease_GO_filtered_with_genes, "GO_results/all_gene_disease_GO_filtered_with_genes.csv")

pdf("GO_results/all_gene_disease_GO.pdf", width = 20, height = 35)
ggplot(all_gene_disease_GO_filtered_with_genes, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = -log10(over_represented_pvalue_corrected), fill = over_represented_pvalue_corrected)) + 
  geom_bar(stat = "identity") + facet_wrap(~ ontology, scales = "free_y", ncol = 1) + 
  coord_flip() + scale_fill_gradient(low = "red", high = "blue") + 
  labs(x = "GO Term", y = "-log10(pvalue_adjusted)", fill = "pvalue_adjusted") + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1)) 

ggplot(all_gene_disease_GO_filtered_with_genes, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = hitsPerc, size = numDEInCat, color = over_represented_pvalue_corrected)) + 
  geom_point() + 
  facet_wrap(~ ontology, scales = "free_y", ncol = 1) + 
  coord_flip() + scale_color_gradient(low = "red", high = "blue") + 
  labs(title = "", x = "GO Term", y = "Hits (%)", size = "Gene count", color = "adjusted pvalue") + 
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5))

p3 <- ggplot(all_gene_disease_GO_filtered_with_genes, aes(x = numDEInCat, y = -log10(over_represented_pvalue_corrected), color = over_represented_pvalue_corrected)) +
  geom_point(size = 6) +  # Adjust size for better visibility
  scale_color_gradient(low = "red", high = "blue") +  # Color gradient for p-values
  # theme_classic() +
  labs(x = "Gene Count", y = "-log10(pvalue_adjusted)", color = "pvalue_adjusted", title = "") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("GO_results/all_gene_disease_GO_02.pdf", plot = p3, width = 12, height = 10, dpi = 300)
dev.off()

all_GO_filtered_with_genes <- rbind(all_gene_normal_GO_filtered_with_genes, all_gene_disease_GO_filtered_with_genes)
ggplot(all_GO_filtered_with_genes, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = -log10(over_represented_pvalue_corrected), fill = over_represented_pvalue_corrected)) + 
  geom_bar(stat = "identity") + facet_wrap(~ ontology, scales = "free_y", ncol = 1) + 
  coord_flip() + scale_fill_gradient(low = "red", high = "blue") + 
  facet_grid(Condition ~ .,space="free",scales="free")
  labs(x = "GO Term", y = "-log10(pvalue_adjusted)", fill = "pvalue_adjusted") + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1)) 


all_GO_filtered_with_genes_melt <- melt(all_GO_filtered_with_genes[, c(4,6,7,8,9,10)], id = c("term", "ontology", "Condition", "hitsPerc", "numDEInCat"), value.name = "over_represented_pvalue_corrected")
all_GO_filtered_with_genes_melt_01 <- all_GO_filtered_with_genes_melt[all_GO_filtered_with_genes_melt$ontology == c("BP"),] %>% 
  mutate(Condition = factor(Condition, level = c("Control", "IHD")))
all_GO_filtered_with_genes_melt_02 <- all_GO_filtered_with_genes_melt[all_GO_filtered_with_genes_melt$ontology == c("MF"),] %>% 
  mutate(Condition = factor(Condition, level = c("Control", "IHD")))
all_GO_filtered_with_genes_melt_03 <- all_GO_filtered_with_genes_melt[all_GO_filtered_with_genes_melt$ontology == c("CC"),] %>% 
  mutate(Condition = factor(Condition, level = c("Control", "IHD")))
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])

pdf("GO_results/all_gene_GO_barplot.pdf", width = 20, height = 15)
p11 <- ggplot(all_GO_filtered_with_genes_melt_01, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = -log10(over_represented_pvalue_corrected), fill = Condition)) +
  geom_col(position = position_dodge(preserve = "single")) + facet_wrap(ontology ~ ., scales="free_x") + scale_fill_manual(values = color_set) + 
  theme_classic() + coord_flip() + labs(title = "", x = "Group", y = "-log10(FDR-adjusted P-value)") + theme(text = element_text(size=20))
p11
p21 <- ggplot(all_GO_filtered_with_genes_melt_02, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = -log10(over_represented_pvalue_corrected), fill = Condition)) +
  geom_col(position = position_dodge(preserve = "single")) + facet_wrap(ontology ~ ., scales="free_x") + scale_fill_manual(values = color_set) + 
  theme_classic() + coord_flip() + labs(title = "", x = "Group", y = "-log10(FDR-adjusted P-value)") + theme(text = element_text(size=20))
p21
p31 <- ggplot(all_GO_filtered_with_genes_melt_03, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = -log10(over_represented_pvalue_corrected), fill = Condition)) +
  geom_col(position = position_dodge(preserve = "single")) + facet_wrap(ontology ~ ., scales="free_x") + scale_fill_manual(values = color_set) + 
  theme_classic() + coord_flip() + labs(title = "", x = "Group", y = "-log10(FDR-adjusted P-value)") + theme(text = element_text(size=20))
p31
dev.off()

##### dotplot and volcano plot for Control
all_GO_filtered_with_genes_melt_01 <- all_GO_filtered_with_genes_melt[all_GO_filtered_with_genes_melt$ontology == c("BP") & all_GO_filtered_with_genes_melt$Condition == "Control",]
all_GO_filtered_with_genes_melt_02 <- all_GO_filtered_with_genes_melt[all_GO_filtered_with_genes_melt$ontology == c("MF") & all_GO_filtered_with_genes_melt$Condition == "Control",]
all_GO_filtered_with_genes_melt_03 <- all_GO_filtered_with_genes_melt[all_GO_filtered_with_genes_melt$ontology == c("CC") & all_GO_filtered_with_genes_melt$Condition == "Control",]

pdf("GO_results/all_gene_GO_dotplot_Control.pdf", width = 20, height = 15)
p12 <- ggplot(all_GO_filtered_with_genes_melt_01, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = hitsPerc, size = numDEInCat, color = over_represented_pvalue_corrected)) + 
  geom_point() + facet_wrap(~ ontology, scales = "free_y", ncol = 1) + coord_flip() + scale_color_gradient(low = "red", high = "blue") + 
  labs(title = "", x = "GO Term", y = "Hits (%)", size = "Gene count", color = "FDR-adjusted P-value") + theme(text = element_text(size=20))
p12

p22 <- ggplot(all_GO_filtered_with_genes_melt_02, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = hitsPerc, size = numDEInCat, color = over_represented_pvalue_corrected)) + 
  geom_point() + facet_wrap(~ ontology, scales = "free_y", ncol = 1) + coord_flip() + scale_color_gradient(low = "red", high = "blue") + 
  labs(title = "", x = "GO Term", y = "Hits (%)", size = "Gene count", color = "FDR-adjusted P-value") + theme(text = element_text(size=20))
p22
p32 <- ggplot(all_GO_filtered_with_genes_melt_03, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = hitsPerc, size = numDEInCat, color = over_represented_pvalue_corrected)) + 
  geom_point() + facet_wrap(~ ontology, scales = "free_y", ncol = 1) + coord_flip() + scale_color_gradient(low = "red", high = "blue") + 
  labs(title = "", x = "GO Term", y = "Hits (%)", size = "Gene count", color = "FDR-adjusted P-value") + theme(text = element_text(size=20))
p32
dev.off()

pdf("GO_results/all_gene_GO_volcano_Control.pdf", width = 14, height = 10)
p13 <- ggplot(all_GO_filtered_with_genes_melt_01, aes(x = numDEInCat, y = -log10(over_represented_pvalue_corrected), color = -log10(over_represented_pvalue_corrected))) +
  geom_point(size = 6) + scale_color_gradient(low = "blue", high = "red") + 
  theme_classic() + labs(x = "Gene Count", y = "-log10(FDR-adjusted P-value)", color = "-log10(FDR-adjusted P-value)", title = "") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1))
p13

p23 <- ggplot(all_GO_filtered_with_genes_melt_02, aes(x = numDEInCat, y = -log10(over_represented_pvalue_corrected), color = -log10(over_represented_pvalue_corrected))) +
  geom_point(size = 6) + scale_color_gradient(low = "blue", high = "red") + 
  theme_classic() + labs(x = "Gene Count", y = "-log10(FDR-adjusted P-value)", color = "-log10(FDR-adjusted P-value)", title = "") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1))
p23
p33 <- ggplot(all_GO_filtered_with_genes_melt_03, aes(x = numDEInCat, y = -log10(over_represented_pvalue_corrected), color = -log10(over_represented_pvalue_corrected))) +
  geom_point(size = 6) + scale_color_gradient(low = "blue", high = "red") + 
  theme_classic() + labs(x = "Gene Count", y = "-log10(FDR-adjusted P-value)", color = "-log10(FDR-adjusted P-value)", title = "") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1))
p33
dev.off()


##### dotplot and volcano plot for IHD
all_GO_filtered_with_genes_melt_01 <- all_GO_filtered_with_genes_melt[all_GO_filtered_with_genes_melt$ontology == c("BP") & all_GO_filtered_with_genes_melt$Condition == "IHD",]
all_GO_filtered_with_genes_melt_02 <- all_GO_filtered_with_genes_melt[all_GO_filtered_with_genes_melt$ontology == c("MF") & all_GO_filtered_with_genes_melt$Condition == "IHD",]
all_GO_filtered_with_genes_melt_03 <- all_GO_filtered_with_genes_melt[all_GO_filtered_with_genes_melt$ontology == c("CC") & all_GO_filtered_with_genes_melt$Condition == "IHD",]

pdf("GO_results/all_gene_GO_dotplot_IHD.pdf", width = 20, height = 15)
p12 <- ggplot(all_GO_filtered_with_genes_melt_01, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = hitsPerc, size = numDEInCat, color = over_represented_pvalue_corrected)) + 
  geom_point() + facet_wrap(~ ontology, scales = "free_y", ncol = 1) + coord_flip() + scale_color_gradient(low = "red", high = "blue") + 
  labs(title = "", x = "GO Term", y = "Hits (%)", size = "Gene count", color = "FDR-adjusted P-value") + theme(text = element_text(size=20))
p12

p22 <- ggplot(all_GO_filtered_with_genes_melt_02, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = hitsPerc, size = numDEInCat, color = over_represented_pvalue_corrected)) + 
  geom_point() + facet_wrap(~ ontology, scales = "free_y", ncol = 1) + coord_flip() + scale_color_gradient(low = "red", high = "blue") + 
  labs(title = "", x = "GO Term", y = "Hits (%)", size = "Gene count", color = "FDR-adjusted P-value") + theme(text = element_text(size=20))
p22
p32 <- ggplot(all_GO_filtered_with_genes_melt_03, aes(x = reorder(term, -log10(over_represented_pvalue_corrected)), y = hitsPerc, size = numDEInCat, color = over_represented_pvalue_corrected)) + 
  geom_point() + facet_wrap(~ ontology, scales = "free_y", ncol = 1) + coord_flip() + scale_color_gradient(low = "red", high = "blue") + 
  labs(title = "", x = "GO Term", y = "Hits (%)", size = "Gene count", color = "FDR-adjusted P-value") + theme(text = element_text(size=20))
p32
dev.off()

pdf("GO_results/all_gene_GO_volcano_IHD.pdf", width = 14, height = 10)
p13 <- ggplot(all_GO_filtered_with_genes_melt_01, aes(x = numDEInCat, y = -log10(over_represented_pvalue_corrected), color = -log10(over_represented_pvalue_corrected))) +
  geom_point(size = 6) + scale_color_gradient(low = "blue", high = "red") + 
  theme_classic() + labs(x = "Gene Count", y = "-log10(FDR-adjusted P-value)", color = "-log10(FDR-adjusted P-value)", title = "") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1))
p13

p23 <- ggplot(all_GO_filtered_with_genes_melt_02, aes(x = numDEInCat, y = -log10(over_represented_pvalue_corrected), color = -log10(over_represented_pvalue_corrected))) +
  geom_point(size = 6) + scale_color_gradient(low = "blue", high = "red") + 
  theme_classic() + labs(x = "Gene Count", y = "-log10(FDR-adjusted P-value)", color = "-log10(FDR-adjusted P-value)", title = "") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1))
p23
p33 <- ggplot(all_GO_filtered_with_genes_melt_03, aes(x = numDEInCat, y = -log10(over_represented_pvalue_corrected), color = -log10(over_represented_pvalue_corrected))) +
  geom_point(size = 6) + scale_color_gradient(low = "blue", high = "red") + 
  theme_classic() + labs(x = "Gene Count", y = "-log10(FDR-adjusted P-value)", color = "-log10(FDR-adjusted P-value)", title = "") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1))
p33
dev.off()