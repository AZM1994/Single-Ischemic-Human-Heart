library(dplyr)
library(readr)
# library(purrr)
library(ggplot2)
# library(ggpubr)
# library(tibble)
library(tidyr)
# library(pheatmap)
# library(gridExtra)
# library(grid)
library(ComplexHeatmap)
# library(circlize)
# library(dplyr)

setwd("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration/Mutation_Enrichment_Analysis/extract_germline_GATK")
ctrl_color_palette <- colorRampPalette(c("skyblue1","dodgerblue4"))
dis_color_palette <- colorRampPalette(c("pink1","firebrick3"))
ctrl_dis_color <- c(ctrl_color_palette(9)[7], dis_color_palette(4)[3])

protein_altering_mutation <- c("splicing", "exonic;splicing", "frameshift deletion", "frameshift insertion", "nonframeshift deletion", "nonsynonymous SNV", "stopgain")
dna_repair_genes <- c("OGG1","MUTYH","APEX1","APE1","PARP1","XRCC1","XPA","XPC","ERCC1","XPF", "ERCC2","XPD","RAD23B","CSA","CSB",
                      "MLH1","MSH2","MSH3","MSH6","PMS2", "EXO1","PALB2","ATM","RAD51","RAD52","KU70","XRCC6","LIG4","XRCC4")
oxidative_stress_genes <- c("GPX","CYP1A1","HIF1A","HIF2A","HIF3A","ARNT","EPO","SLC2A1","LDHA","CA9","BNIP3","PDK1","EGLN2","EGLN1","EGLN3",
                     "VHL","HIF1AN","ANGPT2","TGFB1","NOS2","SLC2A3","ADM","CXCR4","HK2","PFKFB3","SERPINE1","IGFBP3","LOX","OPHOS")
APOBEC_genes <- c("APOBEC1", "APOBEC2", "APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D", "APOBEC3E", "APOBEC3F", "APOBEC3G", "APOBEC3H", "APOBEC4", "AICDA", "AID")

files <- list.files(path = "5-Germline_by_donor", pattern = "_germline_filtered.csv$", full.names = TRUE)
selected_cols <- cols_only(Chr = col_character(), Start = col_double(), End = col_double(), Ref = col_character(), Alt = col_character(), 
                           Func.refGene = col_character(), Gene.refGene = col_character(), ExonicFunc.refGene = col_character(), AAChange.refGene = col_character(), 
                           AF_popmax = col_character(), AF = col_character(), AF_status = col_character(), donor = col_character(), Condition = col_character())

all_germline <- bind_rows(lapply(files, function(f) {
  donor_id <- sub("_germline_filtered.csv", "", f)
  df <- read_csv(f, col_types = selected_cols, show_col_types = FALSE)})) %>% 
  mutate(donor = sub("^([0-9]+).*", "\\1", donor)) %>% 
  filter(AF_status %in% c("Novel", "Extremely rare (<0.01%)", "Ultra-rare (<0.1%)")) %>% 
  filter(Func.refGene %in% c("splicing", "exonic;splicing") | ExonicFunc.refGene %in% c("stopgain", "stoploss", "nonsynonymous SNV")) %>% 
  filter(Gene.refGene %in% c(dna_repair_genes, oxidative_stress_genes)) %>% 
  mutate(Pathway = case_when(Gene.refGene %in% dna_repair_genes ~ "DNA Repair", 
                             Gene.refGene %in% oxidative_stress_genes  ~ "Oxidative Stress", 
                             # Gene.refGene %in% APOBEC_genes  ~ "APOBEC",
                             TRUE ~ "Other"))

# write.csv(all_germline, "6-results/rare_germline_protein_altering.csv", row.names = FALSE)

################################################################################
## burden analysis
burden_df <- all_germline %>% 
  group_by(donor, Condition, Pathway, AF_status) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  mutate(AF_status = factor(AF_status, levels = c("Ultra-rare (<0.1%)", "Extremely rare (<0.01%)", "Novel")))

p_burden_stacked <- ggplot(burden_df, aes(x = as.character(donor), y = n, fill = AF_status)) +
  geom_bar(stat = "identity") + facet_grid( ~ Condition, scales = "free_x", space = "free_x") +
  labs(y = "Number of germline variants", x = "Donor") +
  scale_fill_manual(values = c("Novel" = "#EE6363", "Extremely rare (<0.01%)" = "#FFB347", "Ultra-rare (<0.1%)" = "#BA55D3")) + theme_linedraw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 24), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("6-Results/burden_stacked.pdf", plot = p_burden_stacked, width = 14, height = 8, dpi = 600)

################################################################################
##### heatmap
## Gene-level stats for ordering
## Number of donors per condition
donor_list <- data.frame(
  donor = c("4402", "1864", "6032", "4638", "1863", "1039", "1940", "5919", "5828", "5657", "1363", "1673", "604", "1743", "1113"), 
  Condition = c("Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "IHD", "IHD", "IHD", "IHD", "IHD"),
  stringsAsFactors = FALSE) %>% 
  mutate(donor_cond = paste0(donor, "_", Condition), 
         donor_cond = factor(donor_cond, levels = c("4402_Control", "1864_Control", "6032_Control", "4638_Control", "1863_Control", "1039_Control", "1940_Control", "5919_Control", 
                                                    "5828_Control", "5657_Control", "1363_IHD", "1673_IHD", "604_IHD", "1743_IHD", "1113_IHD"))) %>% 
  arrange(Condition, donor_cond)

n_ctrl <- donor_list %>% filter(Condition == "Control") %>% nrow()
n_ihd  <- donor_list %>% filter(Condition == "IHD") %>% nrow()

## Gene stats
gene_stats <- all_germline %>% 
  group_by(Gene.refGene, Pathway, Condition) %>% 
  summarise(n_donors = n_distinct(donor), .groups = "drop") %>% 
  pivot_wider(names_from = Condition, values_from = n_donors, values_fill = 0) %>% 
  mutate(prop_control = Control / n_ctrl, prop_ihd = IHD / n_ihd, diff_prop = prop_ihd - prop_control)

## Order genes by pathway then diff
gene_order <- gene_stats %>% arrange(Pathway, desc(diff_prop)) %>% pull(Gene.refGene)

## Full gene by donor grid
presence_long <- expand.grid(Gene.refGene = unique(all_germline$Gene.refGene), donor_cond = donor_list$donor_cond, stringsAsFactors = FALSE) %>% 
  left_join(all_germline %>% mutate(donor_cond = paste0(donor, "_", Condition)) %>% 
              distinct(Gene.refGene, donor_cond) %>% mutate(flag = 1), by = c("Gene.refGene", "donor_cond")) %>% 
  mutate(flag = ifelse(is.na(flag), 0, flag))

## Matrix for heatmap
mat_df <- presence_long %>%
  pivot_wider(names_from = donor_cond, values_from = flag, values_fill = 0) %>% 
  filter(Gene.refGene %in% gene_order) %>% 
  mutate(Gene.refGene = factor(Gene.refGene, levels = gene_order)) %>% 
  arrange(Gene.refGene)

presence_mat <- as.matrix(mat_df[, setdiff(colnames(mat_df), "Gene.refGene")])
rownames(presence_mat) <- as.character(mat_df$Gene.refGene)

## Donor order
donor_order_use <- donor_list$donor_cond[donor_list$donor_cond %in% colnames(presence_mat)]
presence_mat <- presence_mat[, donor_order_use, drop = FALSE]

## Pathway coding
gene_pathway_tbl <- all_germline %>% distinct(Gene.refGene, Pathway)
gene_path <- gene_pathway_tbl$Pathway[match(rownames(presence_mat), gene_pathway_tbl$Gene.refGene)]
code_vec <- ifelse(gene_path == "Oxidative Stress", 1, -1)
mat_pathway <- presence_mat * matrix(code_vec, nrow = nrow(presence_mat), ncol = ncol(presence_mat), byrow = FALSE)
mat_pathway <- mat_pathway[gene_order, donor_order_use]

## Donor-level hits (include zeros)
donor_hits <- colSums(mat_pathway != 0)
donor_condition <- ifelse(grepl("_Control$", colnames(mat_pathway)), "Control", "IHD")

## Gene-level diff and pathway
gene_diff <- gene_stats$diff_prop[match(rownames(mat_pathway), gene_stats$Gene.refGene)]
gene_pathway <- gene_stats$Pathway[match(rownames(mat_pathway), gene_stats$Gene.refGene)]

## Top barplot for donor hits
ha_top <- HeatmapAnnotation(
  Hits = anno_barplot(donor_hits, gp = gpar(fill = ifelse(donor_condition == "Control", ctrl_dis_color[1], ctrl_dis_color[2])), border = FALSE, height = unit(3, "cm")), 
  Condition = donor_condition, col = list(Condition = c("Control" = ctrl_dis_color[1], "IHD" = ctrl_dis_color[2])))

## Right barplot for gene-level diff
ha_right <- rowAnnotation(
  # Diff = anno_barplot(gene_diff, gp = gpar(fill = ifelse(gene_pathway == "DNA_repair", "orange", "#9ACD32")), border = FALSE, width = unit(3, "cm")), 
  Pathway = gene_pathway, col = list(Pathway = c("DNA Repair" = "orange", "Oxidative Stress" = "#9ACD32")))

## Heatmap
ht <- Heatmap(mat_pathway, name = "Hit", col = c("-1" = "orange", "0" = "white", "1" = "#9ACD32"), cluster_rows = FALSE, cluster_columns = FALSE, 
              top_annotation = ha_top, right_annotation = ha_right, show_column_names = TRUE, show_row_names = TRUE, row_names_side = "left", 
              heatmap_legend_param = list(at = c(-1, 0, 1), labels = c("DNA Repair hit", "None", "Oxidative Stress hit"), title = "Category"), 
              border = TRUE, rect_gp = gpar(col = "grey80"), row_split = gene_pathway, gap = unit(1, "mm"))

pdf("6-results/donor_gene_heatmap.pdf", width = 12.5, height = 8)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

## per-gene Fisher's exact and FDR
gene_stats_test <- gene_stats %>% rowwise() %>% 
  mutate(fisher_p = fisher.test(matrix(c(IHD, n_ihd - IHD, Control, n_ctrl - Control), nrow = 2))$p.value) %>% 
  ungroup() %>% mutate(fdr = p.adjust(fisher_p, method = "BH"))

sig_genes <- gene_stats_test %>% filter(fdr < 0.1) %>% arrange(fdr) %>% dplyr::select(Gene.refGene, Pathway, IHD, Control, fisher_p, fdr)
