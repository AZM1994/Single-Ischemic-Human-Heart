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

wd_dir <- getwd()
result_dir <- paste0(wd_dir, "/results/LOO_DEG_Cardiomyocyte/GOseq")
dir.create(result_dir, recursive = T)
numDEInCat_threshold = 2
numInCat_threshold = 1000
# gene_length_type = "Gene_length"
gene_length_type = "Exon_length"
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])

##### read gene length
gene_length <- read.delim(paste0(wd_dir, "/data/hg19_refGene.length.tsv"), header = F)
colnames(gene_length) <- c("Gene_symbol", "Transcript", "Gene_length", "Exon_length")
gene_length_deduped <- gene_length[!duplicated(gene_length$Gene_symbol),]
##### read DEGs
DEG_df <- read.csv(paste0(wd_dir, "/results/LOO_DEG_Cardiomyocyte/DEG_up_down_df.csv"))
DEG_up <- DEG_df$gene[DEG_df$regulation == "up"]
DEG_down <- DEG_df$gene[DEG_df$regulation == "down"]

################################################################################
################################################################################
##### GOseq and filtering for UP
mut_gene_up <- gene_length_deduped$Gene_symbol %in% DEG_up
names(mut_gene_up) <- gene_length_deduped$Gene_symbol
pwf_up <- nullp(mut_gene_up, "hg19", bias.data = gene_length_deduped[, gene_length_type])
mut_gene_up_GO <- goseq(pwf_up, "hg19", "geneSymbol") %>% mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% mutate(Regulation = "UP")
# mut_gene_normal_GO_filtered <- mut_gene_up_GO %>%
#   filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
#   filter(over_represented_pvalue < 0.05) %>%
#   mutate(Condition = "Control")

## find gene names in each GO term: Control
go_term_list_up <- mut_gene_up_GO$category  # GO term IDs from goseq
go_genes_list_up <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list_up, keytype = "GOALL", columns = c("SYMBOL"))
genes_in_GO_up <- go_genes_list_up[go_genes_list_up$SYMBOL %in% DEG_up, ]
genes_in_GO_up_collapsed <- aggregate(SYMBOL ~ GOALL, data = genes_in_GO_up, function(x) paste(unique(x), collapse = ", "))
mut_gene_up_GO_with_genes <- merge(mut_gene_up_GO, genes_in_GO_up_collapsed, by.x = "category", by.y = "GOALL") 
mut_gene_up_GO_filtered_with_genes <- mut_gene_up_GO_with_genes %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
  filter(over_represented_pvalue < 0.05)

write.csv(mut_gene_up_GO_with_genes, paste0(result_dir, "/DEG_up_GO_with_genes.csv"))
write.csv(mut_gene_up_GO_filtered_with_genes, paste0(result_dir, "/DEG_up_GO_filtered_with_genes.csv"))

################################################################################
################################################################################
##### GOseq and filtering for DOWN
mut_gene_down <- gene_length_deduped$Gene_symbol %in% DEG_down
names(mut_gene_down) <- gene_length_deduped$Gene_symbol
pwf_down <- nullp(mut_gene_down, bias.data = gene_length_deduped[, gene_length_type])
mut_gene_down_GO <- goseq(pwf_down, "hg19", "geneSymbol") %>% mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% mutate(Regulation = "DOWN")
# mut_gene_disease_GO_filtered <- mut_gene_down_GO %>%
#   filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>%
#   filter(over_represented_pvalue < 0.05) %>%
#   mutate(Condition = "IHD")

## find gene names in each GO term: IHD
go_term_list_down <- mut_gene_down_GO$category  # GO term IDs from goseq
go_genes_list_down <- AnnotationDbi::select(org.Hs.eg.db, keys = go_term_list_down, keytype = "GOALL", columns = c("SYMBOL"))
genes_in_GO_down <- go_genes_list_down[go_genes_list_down$SYMBOL %in% DEG_down, ]
# genes_in_GO_down_collapsed <- genes_in_GO_down %>% group_by(GOALL) %>% 
#   filter(!duplicated(SYMBOL)) %>% 
#   summarise(genes = paste(SYMBOL, collapse = ", "))
genes_in_GO_down_collapsed <- aggregate(SYMBOL ~ GOALL, data = genes_in_GO_down, function(x) paste(unique(x), collapse = ", "))
mut_gene_down_GO_with_genes <- merge(mut_gene_down_GO, genes_in_GO_down_collapsed, by.x = "category", by.y = "GOALL")
mut_gene_down_GO_filtered_with_genes <- mut_gene_down_GO_with_genes %>% 
  filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>% 
  filter(over_represented_pvalue < 0.05)

write.csv(mut_gene_down_GO_with_genes, paste0(result_dir, "/DEG_down_GO_with_genes.csv"))
write.csv(mut_gene_down_GO_filtered_with_genes, paste0(result_dir, "/DEG_down_GO_filtered_with_genes.csv"))

################################################################################
################################################################################
### generate volcano plots for BP
nature_colors <- c("#a6761d", "#228C68", "#1f77b4", "#ff7f0e", "#e377c2", "#F22020")
# nature_colors <- c("#4d79a6", "#a1cde6", "#e05759", "#ff9d99", "#f1cd63", "#79706d","#489893","#8cd17d","#00afba","#e6b800","#fb4d07")
# nature_colors <- c("#a2d2e7", "#67a8cd", "#ffc17f", "#cf9f88", "#6fb3a8", "#b3e19b","#50aa4b","#ff9d9f","#f36569","#3581b7","#cdb6da",
#                    "#704ba3", "#9a7fbd", "#dba9a8", "#e43030", "#e99b78", "#ff8831")
## DEG UP
# cell_cycle_up <- c("G1/S transition of mitotic cell cycle")
# cytoskeletal_up <- c("actin binding", "muscle contraction", "sarcomere organization")
# morphogenesis_development_up <- c("muscle structure development")
# stress_response_up <- c("stress fiber", "activation of immune response", "response to ischemia")
# transporter_activity_up <- c("amino acid transport", "vascular transport", "regulation of metal ion transport", "monoatomic cation transport", "metal ion transport")
# membrane_organization_up <- c("phosphatidylinositol-3,4,5-trisphosphate binding", "nuclear envelope", "mitochondrial envelope", "mitochondrial outer membrane", "focal adhesion")
# signaling_pathways_up <- c("calcium-mediated signaling", "second-messenger-mediated signaling", "calcineurin-NFAT signaling cascade", "toll-like receptor 9 signaling pathway")

GO_BP_up <- mut_gene_up_GO_filtered_with_genes[mut_gene_up_GO_filtered_with_genes$ontology == "BP", ]
# GO_BP_up <- mut_gene_up_GO_filtered_with_genes[mut_gene_up_GO_filtered_with_genes$ontology == "BP", ] %>% 
GO_BP_up <- mut_gene_up_GO_with_genes[mut_gene_up_GO_with_genes$ontology == "BP", ] %>%
  mutate(BP_group = ifelse(term %in% cell_cycle_up, "Cell Cycle",
                           ifelse(term %in% cytoskeletal_up, "Cytoskeletal",
                                  ifelse(term %in% morphogenesis_development_up, "Morphogenesis & Development",
                                         ifelse(term %in% stress_response_up, "Stress Response",
                                                ifelse(term %in% transporter_activity_up, "Transporter Activity",
                                                       ifelse(term %in% membrane_organization_up, "Membrane Organization",
                                                              ifelse(term %in% signaling_pathways_up, "Signaling Pathways", NA)))))))) %>%
  filter(!is.na(BP_group))

## DEG DOWN
cell_cycle_down <- c("G1/S transition of mitotic cell cycle", "regulation of G1/S transition of mitotic cell cycle", "negative regulation of G1/S transition of mitotic cell cycle")
cytoskeletal_down <- c("actin binding", "muscle contraction", "sarcomere organization", "regulation of striated muscle contraction")
morphogenesis_development_down <- c("muscle structure development")
RNA_splicing_down <- c("spliceosomal complex assembly", "alternative mRNA splicing, via spliceosome", "regulation of alternative mRNA splicing, via spliceosome", "RNA helicase activity", "rRNA processing")
stress_response_down <- c("stress fiber", "activation of immune response", "response to ischemia")
transporter_activity_down <- c("amino acid transport", "vascular transport", "regulation of metal ion transport", "monoatomic cation transport", "metal ion transport")
membrane_organization_down <- c("phosphatidylinositol-3,4,5-trisphosphate binding", "nuclear envelope", "mitochondrial envelope", "mitochondrial outer membrane", "focal adhesion")
signaling_pathways_down <- c("calcium-mediated signaling", "second-messenger-mediated signaling", "calcineurin-NFAT signaling cascade", "toll-like receptor 9 signaling pathway")

GO_BP_down <- mut_gene_down_GO_filtered_with_genes[mut_gene_down_GO_filtered_with_genes$ontology == "BP", ]
# # GO_BP_down <- mut_gene_down_GO_filtered_with_genes[mut_gene_down_GO_filtered_with_genes$ontology == "BP", ] %>% 
GO_BP_down <- mut_gene_down_GO_with_genes[mut_gene_down_GO_with_genes$ontology == "BP", ] %>% 
  mutate(BP_group = ifelse(term %in% cell_cycle_up, "Cell Cycle",
                           ifelse(term %in% cytoskeletal_up, "Cytoskeletal",
                                  ifelse(term %in% morphogenesis_development_up, "Morphogenesis & Development",
                                         ifelse(term %in% stress_response_up, "Stress Response",
                                                ifelse(term %in% transporter_activity_up, "Transporter Activity",
                                                       ifelse(term %in% membrane_organization_up, "Membrane Organization", 
                                                              ifelse(term %in% signaling_pathways_up, "Signaling Pathways", NA)))))))) %>% filter(!is.na(BP_group))


GO_BP_all <- rbind(GO_BP_up, GO_BP_down)
# GO_BP_all <- rbind(GO_BP_up, GO_BP_down)
# GO_BP_all_filtered <- GO_BP_all %>% 
#   filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>% filter(over_represented_pvalue < 0.05)
write.csv(GO_BP_all, paste0(result_dir, "/GO_BP_all.csv"))
# GO_BP_all_display <- GO_BP_all[]
pdf(paste0(result_dir, "/GO_BP_volcano.pdf"), width = 12, height = 6)
  p_GO_BP_all <- ggplot(GO_BP_all, aes(x = hitsPerc, y = -log10(over_represented_pvalue), color = Regulation)) + 
    geom_point(size = 2) + geom_hline(yintercept = -log10(0.05), lty = 4, lwd = 0.6, alpha = 0.8) + geom_vline(xintercept = 5, lty = 4, lwd = 0.6, alpha = 0.8) + 
    geom_text_repel(data = subset(GO_BP_all, hitsPerc > 0 | -log10(over_represented_pvalue) > 3), aes(label = term), nudge_y = 0.035, size = 5) +
    scale_color_manual(values = c("#2D6EA8", "#DD555B")) + labs(x = "Hits Percentage (%)", y = "-log10(p-value)", color = "BP Category", title = "") + 
    facet_wrap(~ factor(Regulation, levels = c("DOWN", "UP")), scales = "fixed") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12))
  p_GO_BP_all
dev.off()

pdf(paste0(result_dir, "/GO_BP_dotplot.pdf"), width = 50, height = 50)
p12 <- ggplot(GO_BP_all, aes(x = reorder(term, -log10(over_represented_pvalue)), y = hitsPerc, size = hitsPerc, color = over_represented_pvalue)) +
  geom_point(size = 8) + facet_wrap(~ Regulation, scales = "free") + coord_flip() + scale_color_gradient(low = "red", high = "blue") +
  # scale_size_continuous(name = "Size Legend") +
  labs(title = "BP", x = "GO Term", y = "Hits (%)", size = "Gene count", color = "p-value") +
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), text = element_text(face = "bold", size = 12))
p12
dev.off()


################################################################################
################################################################################

# GO_final <- rbind(mut_gene_up_GO_with_genes, mut_gene_down_GO_with_genes)
# GO_final_melt <- melt(GO_final[, c(2,4,5,6,7,8,9)], id = c("term", "ontology", "Regulation", "hitsPerc", "numDEInCat", "numInCat"), value.name = "over_represented_pvalue")
# write.csv(GO_final_melt_MF, paste0(result_dir, "/GO_final_melt_MF.csv"))
# GO_final_melt_MF <- GO_final_melt[GO_final_melt$ontology == c("MF"), ] %>% mutate(Regulation = factor(Regulation, level = c("DOWN", "UP"))) %>% 
#   # filter(numDEInCat >= numDEInCat_threshold & numInCat <= numInCat_threshold) %>% 
#   # filter(over_represented_pvalue < 0.05) %>% 
#   complete(term, ontology, Regulation, fill = list(over_represented_pvalue = 1)) %>%
#   # mutate(MF_group = ifelse(term %in% Kinase_receptor_IHD, "Kinase Receptor", ifelse(term %in% Calcium_channel_activity_IHD, "Calcium Channel Activity",
#   #                   ifelse(term %in% Enzymatic_activity_IHD, "Enzymatic Activity", ifelse(term %in% Protein_binding_IHD, "Protein Binding", 
#   #                   ifelse(term %in% Rcceptor_activitie_Ctrl, "Receptor Activitie", 
#   #                   ifelse(term %in% Cell_Adhesion_mediator_activity_Ctrl, "Cell Adhesion Mediator Activity", "not selected"))))))) %>%
#   # filter(MF_group != c("not selected")) %>% 
#   # mutate(enrich_group = ifelse(MF_group %in% c("Kinase Receptor", "Calcium Channel Activity", "Enzymatic Activity", "Protein Binding"), "IHD Enriched",
#   #                                              ifelse(MF_group %in% c("Receptor Activitie", "Cell Adhesion Mediator Activity"), "Control Enriched", NA))) %>%
#   mutate(enrich_group = factor(enrich_group, levels = c("DOWN Enriched", "UP Enriched")))

# pdf(paste0(result_dir, "/del_GO_MF_barplot.pdf"), width = 15, height = 20)
# p21 <- ggplot(GO_final_melt_MF, aes(x = reorder(term, -log10(over_represented_pvalue)), y = -log10(over_represented_pvalue), fill = Regulation)) +
#   geom_col(position = position_dodge2(reverse = TRUE), width = 0.8) +  geom_hline(yintercept=-log10(0.05),linetype="dashed") + 
#   # facet_grid(enrich_group + MF_group ~ ., space = "free", scales = "free") + 
#   # facet_grid(enrich_group ~ ., space = "free", scales = "free") + 
#   scale_fill_manual(values = color_set) + 
#   coord_flip() + labs(title = "", x = "GO Term", y = "-log10(p-value)") + theme(text = element_text(size=20), legend.position="bottom")
# p21
# dev.off()

pdf(paste0(result_dir, "/GO_dotplot_down.pdf"), width = 10, height = 7)
  mut_gene_down_GO_filtered_with_genes %>% 
    top_n(30, wt=-over_represented_pvalue) %>% 
    mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
    ggplot(aes(x=hitsPerc, y=term, colour=over_represented_pvalue, size=numDEInCat)) + 
    facet_grid(ontology ~ ., space = "free", scales = "free") + 
    geom_point() + expand_limits(x=0) + labs(x="Hits (%)", y="GO term", colour="p value", size="Count") + 
    theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey30", linetype = "dashed"), 
                             panel.grid.minor = element_line(color = "grey30", linetype = "dashed"), text = element_text(face = "bold", size = 12))
dev.off()

pdf(paste0(result_dir, "/GO_dotplot_up.pdf"), width = 15, height = 8)
mut_gene_up_GO_filtered_with_genes %>% filter(term %in% c("cellular response to oxygen levels", "cellular response to hypoxia", "cellular response to hydrogen peroxide", 
                                              "response to hypoxia", "response to reactive oxygen species", "cellular response to reactive oxygen species", 
                                              "response to inorganic substance", "response to gamma radiation", "cytoplasmic stress granule", 
                                              "response to decreased oxygen levels", "cellular response to decreased oxygen levels", "response to hydrogen peroxide") | over_represented_pvalue < 7.587488e-04) %>% 
  # top_n(20, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term, colour=over_represented_pvalue, size=numDEInCat)) + 
  facet_grid(ontology ~ ., space = "free", scales = "free") + 
  geom_point() + expand_limits(x=0) + labs(x="Hits (%)", y="GO term", colour="p value", size="Count") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.minor = element_line(color = "grey", linetype = "dashed"), text = element_text(face = "bold", size = 12))
dev.off()