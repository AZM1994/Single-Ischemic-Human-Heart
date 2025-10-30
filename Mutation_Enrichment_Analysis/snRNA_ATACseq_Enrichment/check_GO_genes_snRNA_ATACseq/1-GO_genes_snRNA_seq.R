library(ggplot2)
library(dplyr)
library(tidyr)
library(Seurat)

setwd("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration/Mutation_Enrichment_Analysis")
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])

##### read in metadata
Hypoxia_PTA_all_SCAN2 <- readRDS("data/1-sSNV_SCAN2_df_filtered.rds") %>% 
  dplyr::select(c("Cell_ID", "Age", "Gender", "Case_ID", "Condition", "snv.burden", "snv.rate.per.gb"))
Hypoxia_PTA_all_SCAN2_collapsed <- Hypoxia_PTA_all_SCAN2 %>% 
  distinct(Case_ID, .keep_all = TRUE)
Condition_list <- unique(Hypoxia_PTA_all_SCAN2$Condition)

GO_BP_Ctrl <- read.csv("sSNV_GO_Enrichment/GOseq/deleterious/Control_GO.deleterious_mutation.csv", header = TRUE)
GO_BP_IHD <- read.csv("sSNV_GO_Enrichment/GOseq/deleterious/IHD_GO.deleterious_mutation.csv", header = TRUE)
Ctrl_GO_genes <- unique(unlist(strsplit(GO_BP_Ctrl$genes, ",\\s*")))
IHD_GO_genes <- unique(unlist(strsplit(GO_BP_IHD$genes, ",\\s*")))

##### read in CM scRNA-seq data
Seurat_CM <- readRDS("data/Seurat.obj_sub_clustering_CM_only.RDS")
group_num <- 8

expr_level_all <- data.frame(AverageExpression(Seurat_CM, group.by = "condition", slot = "data")$RNA) %>% 
  setNames(Condition_list) %>% 
  mutate(Gene.refGene = rownames(.)) %>% 
  pivot_longer(cols = -Gene.refGene, names_to = "Condition", values_to = "average_expr_level") %>% 
  filter(average_expr_level > 0) %>%
  group_by(Condition) %>% 
  mutate(decile = ntile(average_expr_level, n = group_num)) %>% 
  ungroup()

# Create GO gene list by condition
GO_genes_df <- bind_rows(
  data.frame(Gene = unique(unlist(strsplit(GO_BP_Ctrl$genes, ",\\s*"))), Condition = "Control"), 
  data.frame(Gene = unique(unlist(strsplit(GO_BP_IHD$genes, ",\\s*"))), Condition = "IHD")
)

# Merge GO gene info with expression deciles
go_expr_all <- expr_level_all %>% 
  inner_join(GO_genes_df, by = c("Gene.refGene" = "Gene", "Condition"))

# Get decile counts and expected uniform distribution
decile_counts <- go_expr_all %>%
  group_by(Condition, decile) %>%
  summarise(n = n(), .groups = "drop")

expected_counts <- go_expr_all %>%
  group_by(Condition) %>%
  summarise(expected = n() / 8)

decile_counts <- decile_counts %>%
  left_join(expected_counts, by = "Condition")

# Chi-square test per condition
pvals <- decile_counts %>%
  group_by(Condition) %>%
  summarise(chisq_p = chisq.test(n, p = rep(1/length(n), length(n)))$p.value)

p_barplot <- ggplot(decile_counts, aes(x = factor(decile), y = n, fill = Condition)) + 
  geom_bar(stat = "identity") +
  geom_hline(aes(yintercept = expected), linetype = "dashed", color = "black") +
  facet_wrap(~ Condition) +
  geom_text(data = pvals, aes(x = 4.5, y = max(decile_counts$n) * 0.95, label = paste0("p = ", signif(chisq_p, 3))), inherit.aes = FALSE, size = 6) + 
  scale_fill_manual(values = c("Control" = color_set[1], "IHD" = color_set[2]), guide = "legend") + 
  labs(x = "Gene expression level", y = "Number of GO Genes") + theme_linedraw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), 
        text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
ggsave("snRNA_ATACseq_Enrichment/check_GO_genes_snRNA_ATACseq/GO_genes_snRNAseq.pdf", plot = p_barplot, width = 10, height = 5, dpi = 600)
