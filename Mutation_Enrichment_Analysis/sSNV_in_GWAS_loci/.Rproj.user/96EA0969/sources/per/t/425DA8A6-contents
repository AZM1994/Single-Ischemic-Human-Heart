library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# setwd("/Users/zhemingan/Documents/BCH_research/SNV_enrichment_analysis/heart_PTA_Cases")
##### read in normal and disease data
Control.SNV_gene_list <- read.csv("./data/annovar_annotation/normal/heart_PTA_Cases.all_age.normal.hg19_multianno.csv", header = TRUE)
heart_PTA_Cases_Control_vcf <- read.table("./data/annovar_annotation/normal/heart_PTA_Cases.all_age.normal_ssnv.vcf", sep = "\t")
Control.SNV_gene_list$Cell_ID <- heart_PTA_Cases_Control_vcf$V8
Control.SNV_gene_list$Case_ID <- str_extract(Control.SNV_gene_list$Cell_ID, "[^_]+")
Control.SNV_gene_list$Condition <- "Control"
Control.SNV_gene_list <- Control.SNV_gene_list[c("Chr", "Start", "End", "Ref", "Alt", "Cell_ID", "Case_ID", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "Condition")] %>% 
  filter(Func.refGene != "intergenic")
Control.SNV_gene_list_exonic <- Control.SNV_gene_list[Control.SNV_gene_list$Func.refGene %in% c("exonic", "exonic;splicing", "splicing"), ]
Control.SNV_gene_list_nonsynonymous <- Control.SNV_gene_list[Control.SNV_gene_list$ExonicFunc.refGene %in% c("nonsynonymous SNV", "stopgain"), ]

IHD.SNV_gene_list <- read.csv("./data/annovar_annotation/disease/heart_PTA_Cases.all_age.disease.hg19_multianno.csv", header = TRUE)
heart_PTA_Cases_IHD_vcf <- read.table("./data/annovar_annotation/disease/heart_PTA_Cases.all_age.disease_ssnv.vcf", sep = "\t")
IHD.SNV_gene_list$Cell_ID <- heart_PTA_Cases_IHD_vcf$V8
IHD.SNV_gene_list$Case_ID <- str_extract(IHD.SNV_gene_list$Cell_ID, "[^_]+")
IHD.SNV_gene_list$Condition <- "IHD"
IHD.SNV_gene_list <- IHD.SNV_gene_list[c("Chr", "Start", "End", "Ref", "Alt", "Cell_ID", "Case_ID", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "Condition")] %>% 
  filter(Func.refGene != "intergenic")
IHD.SNV_gene_list_exonic <- IHD.SNV_gene_list[IHD.SNV_gene_list$Func.refGene %in% c("exonic", "exonic;splicing", "splicing"), ]
IHD.SNV_gene_list_nonsynonymous <- IHD.SNV_gene_list[IHD.SNV_gene_list$ExonicFunc.refGene %in% c("nonsynonymous SNV", "stopgain"), ]

combined.SNV_gene_list <- rbind(Control.SNV_gene_list, IHD.SNV_gene_list)
combined.SNV_gene_list_exonic <- rbind(Control.SNV_gene_list_exonic, IHD.SNV_gene_list_exonic)
combined.SNV_gene_list_nonsynonymous <- rbind(Control.SNV_gene_list_nonsynonymous, IHD.SNV_gene_list_nonsynonymous)

reported_IHD_loci <- read.csv("./data/hot_gene_list_heart_pta_02.csv", header = FALSE)
reported_IHD_loci <- reported_IHD_loci[!duplicated(reported_IHD_loci$V1), ]

### Split Gene.refGene by ";" and filter for reported_IHD_loci
df_long <- combined.SNV_gene_list %>% 
  separate_rows(Gene.refGene, sep = ";") %>% 
  filter(Gene.refGene %in% reported_IHD_loci)

### Prepare heatmap data, separate by Condition
heatmap_data <- df_long %>% 
  select(Condition, Cell_ID, Gene.refGene, Func.refGene) %>% 
  distinct() %>% 
  mutate(Func.refGene = factor(Func.refGene))

### Order genes by the number of unique Cell_IDs within each condition
gene_order <- heatmap_data %>% 
  group_by(Condition, Gene.refGene) %>% 
  summarise(Num_Cell_IDs = n_distinct(Cell_ID)) %>% 
  spread(key = Condition, value = Num_Cell_IDs, fill = 0) %>% 
  arrange(IHD, Control) %>% 
  mutate(Gene.refGene = factor(Gene.refGene, levels = unique(Gene.refGene)))

heatmap_data <- heatmap_data %>% 
  mutate(Gene.refGene = factor(Gene.refGene, levels = levels(gene_order$Gene.refGene)))

### Order Cell ID by the number of sSNVs within each condition
cell_order <- heatmap_data %>% 
  group_by(Condition, Cell_ID) %>% 
  summarise(Num_SNVs = n()) %>% 
  arrange(Condition, desc(Num_SNVs)) %>% 
  mutate(Cell_ID = factor(Cell_ID, levels = unique(Cell_ID)))

heatmap_data <- heatmap_data %>%
  left_join(cell_order, by = c("Condition", "Cell_ID")) %>%
  mutate(Cell_ID = factor(Cell_ID, levels = unique(cell_order$Cell_ID)))

### Horizontal bar data (number of Cell_IDs per gene, split by Condition)
horizontal_bar_data <- heatmap_data %>%
  group_by(Condition, Gene.refGene) %>%
  summarise(Num_Cell_IDs = n_distinct(Cell_ID)) %>%
  arrange(Condition, desc(Num_Cell_IDs))

### Vertical bar data (number of SNVs per Cell_ID, split by Condition)
vertical_bar_data <- heatmap_data %>%
  group_by(Condition, Cell_ID) %>%
  summarise(Num_SNVs = n()) %>%
  arrange(Condition, desc(Num_SNVs))

### Create heatmap with facets by Condition
heatmap_plot <- ggplot(heatmap_data, aes(x = Cell_ID, y = Gene.refGene, fill = Func.refGene)) + 
  geom_tile(color = "white") + scale_fill_manual(values = c("#66B032", "#FFC34D")) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.05), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0), axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(~Condition, scales = "free_x", space = "free_x") + labs(x = "Cell ID", y = "Reported IHD-associated loci", fill = "Func.refGene")

### Create horizontal bar plot (faceted by Condition)
horizontal_bar_plot <- ggplot(horizontal_bar_data, aes(x = Num_Cell_IDs, y = Gene.refGene, fill = Condition)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = c("Control" = "#2D6EA8", "IHD" = "#DD555B")) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.05), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0)) + 
  facet_grid(~Condition, scales = "free_x", space = "free_x") + labs(x = "Number of Cells", y = "")

### Create vertical bar plot (faceted by Condition)
vertical_bar_plot <- ggplot(vertical_bar_data, aes(x = Cell_ID, y = Num_SNVs, fill = Condition)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = c("Control" = "#2D6EA8", "IHD" = "#DD555B")) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.05), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0), axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(~Condition, scales = "free_x", space = "free_x") + labs(x = "", y = "Number of sSNVs per cell")

### Combine plots with patchwork
combined_plot <- vertical_bar_plot + plot_spacer() + heatmap_plot + horizontal_bar_plot + 
  plot_layout(ncol = 2, heights = c(1, 4), widths = c(4, 1))

ggsave("./results/cell_specific_SNV_count_in_reported_IHD_loci.pdf", plot = combined_plot, width = 35, height = 20, dpi = 600)
