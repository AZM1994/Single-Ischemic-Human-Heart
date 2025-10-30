library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

setwd("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration/Mutation_Enrichment_Analysis/")
result_dir <- "sSNV_in_GWAS_loci"
disease_type <- "IHD" # IHD, CAD
plot_type <- "all" # all, functional

### read metadata and genomic context data
SCAN2_df <- readRDS("data/1-sSNV_SCAN2_df_filtered.rds") %>% as.data.frame() %>% dplyr::select(Cell_ID, Age, Gender, Case_ID, Condition, snv.burden, snv.rate.per.gb)
selected_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene")

genomic_SCAN2_df <- c()
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
    genomic_SCAN2_df <- rbind(genomic_SCAN2_df, genomic_SCAN2_df_temp)
  }
}

genomic_SCAN2_df <- genomic_SCAN2_df %>% filter(Func.refGene != "intergenic")
reported_loci <- read.csv(paste0("./data/reported_", disease_type, "_loci.csv"), header = FALSE) %>% filter(!duplicated(V1))
all_cells <- genomic_SCAN2_df %>% dplyr::select(Cell_ID, Condition) %>% distinct()
df_long <- genomic_SCAN2_df %>% filter(Gene.refGene %in% reported_loci$V1) %>% right_join(all_cells, by = c("Cell_ID", "Condition"))

### Prepare heatmap data, separate by Condition
heatmap_data <- df_long %>% 
  dplyr::select(Condition, Cell_ID, Gene.refGene, Func.refGene) %>% 
  distinct() %>% 
  mutate(Func.refGene = factor(Func.refGene)) %>% 
  { 
    if (plot_type == "functional") {
      filter(., Func.refGene %in% c("exonic", "ncRNA_exonic", "UTR5", "UTR3"))
    } else {
      .
    }
  } %>% 
  right_join(all_cells, by = c("Condition", "Cell_ID"))

### Order genes by the number of unique Cell_IDs within each condition
gene_order <- heatmap_data %>% 
  group_by(Condition, Gene.refGene) %>% 
  summarise(Num_Cell_IDs = n_distinct(Cell_ID)) %>% 
  spread(key = Condition, value = Num_Cell_IDs, fill = 0) %>% 
  arrange(IHD, Control) %>% 
  mutate(Gene.refGene = factor(Gene.refGene, levels = unique(Gene.refGene))) %>% 
  filter(!is.na(Gene.refGene))
write.csv(gene_order$Gene.refGene, paste0(result_dir, "/2-table/", disease_type, "_loci_found.csv"))

heatmap_data <- heatmap_data %>% 
  mutate(Gene.refGene = factor(Gene.refGene, levels = levels(gene_order$Gene.refGene)))

### Order Cell ID by the number of sSNVs within each condition
cell_order <- heatmap_data %>% 
  group_by(Condition, Cell_ID) %>% 
  summarise(Num_SNVs = sum(!is.na(Gene.refGene)), .groups = "drop") %>% 
  arrange(Condition, desc(Num_SNVs))

heatmap_data <- heatmap_data %>% 
  left_join(cell_order, by = c("Condition", "Cell_ID")) %>% 
  mutate(Cell_ID = factor(Cell_ID, levels = unique(cell_order$Cell_ID)))

### Horizontal bar data (number of Cell_IDs per gene, split by Condition)
horizontal_bar_data <- heatmap_data %>% 
  group_by(Condition, Gene.refGene) %>% 
  summarise(Num_Cell_IDs = n_distinct(Cell_ID)) %>% 
  arrange(Condition, desc(Num_Cell_IDs)) %>% 
  group_by(Condition) %>% 
  mutate(Percentage = case_when(Condition == "Control" ~ Num_Cell_IDs / 47 * 100, 
                                Condition == "IHD" ~ Num_Cell_IDs / 33 * 100), 
         Percentage = paste0(round(Percentage, 1), "%")) %>% 
  ungroup()

### Vertical bar data (number of SNVs per Cell_ID, split by Condition)
vertical_bar_data <- heatmap_data %>% 
  group_by(Condition, Cell_ID) %>% 
  summarise(Num_SNVs = sum(!is.na(Gene.refGene)), .groups = "drop") %>% 
  arrange(Condition, desc(Num_SNVs)) %>% 
  mutate(Condition = factor(Condition, levels = c("Control", "IHD")))
wilcox_res <- wilcox.test(Num_SNVs ~ Condition, data = vertical_bar_data, alternative = c("two.sided"))
wilcox_label <- paste0("Wilcoxon test Control v.s. IHD, P = ", signif(wilcox_res$p.value, 3))

### Create heatmap with facets by Condition
mut_colors <- c("exonic" = "#ff9d9f", "ncRNA_exonic" = "#40bfc1", "UTR3" = "#704ba3", "UTR5" = "orange", 
                "downstream" = "#cb77b2", "upstream" = "#8B4513", "intronic" = "#66B032", "ncRNA_intronic" = "#FFE55F")
Gene.refGene_break <- if (plot_type == "all" & disease_type == "CAD") {
    NULL
  } else {
    unique(gene_order$Gene.refGene)
  }

heatmap_plot <- ggplot(heatmap_data, aes(x = Cell_ID, y = Gene.refGene, fill = factor(Func.refGene, levels = names(mut_colors)))) + 
  geom_tile(color = "white") + scale_fill_manual(values = mut_colors) + 
  scale_y_discrete(limits = unique(gene_order$Gene.refGene), breaks = Gene.refGene_break) + theme_linedraw() + 
  theme(panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.05), panel.grid.minor = element_blank(), 
        text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0), axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = c(0.98, 0.02), legend.justification = c("right", "bottom")) + 
  facet_grid( ~ Condition, scales = "free_x", space = "free_x") + labs(x = "Cell ID", y = paste0("Reported ", disease_type, "-associated loci"), fill = "Func.refGene") + guides(fill = guide_legend(title = "Mutation Category")) 

### Create horizontal bar plot (faceted by Condition)
horizontal_bar_plot <- ggplot(horizontal_bar_data, aes(x = Num_Cell_IDs, y = Gene.refGene, fill = Condition)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = c("Control" = "#2D6EA8", "IHD" = "#DD555B")) + 
  scale_y_discrete(limits = unique(gene_order$Gene.refGene), breaks = Gene.refGene_break) + theme_linedraw() + 
  theme(panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.05), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), 
        text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0), legend.position = "none") + 
  facet_grid( ~ Condition, scales = "free_x", space = "free_x") + labs(x = "Number of Cells", y = "") + scale_x_continuous(breaks = seq(0, max(horizontal_bar_data$Num_Cell_IDs), by = 2))

### Create vertical bar plot (faceted by Condition)
vertical_bar_plot <- ggplot(vertical_bar_data, aes(x = Cell_ID, y = Num_SNVs, fill = Condition)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = c("Control" = "#2D6EA8", "IHD" = "#DD555B")) + theme_linedraw() + 
  theme(panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.05), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), 
        text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0), axis.text.x = element_blank(), 
        legend.position = c(0.98, 0.98), legend.justification = c("right", "top")) + 
  facet_grid( ~ Condition, scales = "free_x", space = "free_x") + labs(x = "", y = "Number of sSNVs per cell", subtitle = wilcox_label)

### Combine plots with patchwork
combined_plot <- vertical_bar_plot + plot_spacer() + heatmap_plot + horizontal_bar_plot + 
  plot_layout(ncol = 2, heights = c(1, 4), widths = c(8, 1))
  # plot_layout(ncol = 2, heights = c(1, 6), widths = c(8, 1))

ggsave(paste0(result_dir, "/sSNV_in_", disease_type, "_loci_", plot_type,".pdf"), plot = combined_plot, width = 38, height = 24, dpi = 600) # CAD, all
# ggsave(paste0(result_dir, "/sSNV_in_", disease_type, "_loci_", plot_type,".pdf"), plot = combined_plot, width = 38, height = 24, dpi = 600) # CAD, functional
# ggsave(paste0(result_dir, "/sSNV_in_", disease_type, "_loci_", plot_type,".pdf"), plot = combined_plot, width = 38, height = 18, dpi = 600) # IHD, all