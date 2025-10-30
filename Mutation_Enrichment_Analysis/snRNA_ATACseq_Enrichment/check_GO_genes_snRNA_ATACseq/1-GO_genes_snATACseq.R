library(ggplot2)
library(dplyr)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

setwd("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration/Mutation_Enrichment_Analysis/snRNA_ATACseq_Enrichment")
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])

GO_BP_Ctrl <- read.csv("check_GO_genes_snRNA_ATACseq/Control_GO.deleterious_mutation.csv", header = TRUE)
GO_BP_IHD <- read.csv("check_GO_genes_snRNA_ATACseq/IHD_GO.deleterious_mutation.csv", header = TRUE)
Ctrl_GO_genes <- unique(unlist(strsplit(GO_BP_Ctrl$genes, ",\\s*")))
IHD_GO_genes <- unique(unlist(strsplit(GO_BP_IHD$genes, ",\\s*")))

##### read in scRNA-ATACseq data
group_num <- 8
atac_ctrl <- c()
atac_dis <- c()
ATAC_header <- c("chr", "start", "end", "id")
atac_4638_2n <- read.table("./data/scATACseq_preprocessing/ATAC_cov_average_hg19_1kb_bin.4638_2n.bed") %>% setNames(c(ATAC_header, "score_4638_2n"))
atac_5828_2n <- read.table("./data/scATACseq_preprocessing/ATAC_cov_average_hg19_1kb_bin.5828_2n.bed") %>% setNames(c(ATAC_header, "score_5828_2n")) %>% dplyr::select(c("score_5828_2n"))
atac_5919_all <- read.table("./data/scATACseq_preprocessing/ATAC_cov_average_hg19_1kb_bin.5919_all.bed") %>% setNames(c(ATAC_header, "score_5919_all")) %>% dplyr::select(c("score_5919_all"))
atac_ctrl <- cbind(atac_4638_2n, atac_5828_2n, atac_5919_all) %>% 
  mutate(score = rowMeans(dplyr::select(., score_4638_2n, score_5828_2n, score_5919_all), na.rm = TRUE)) %>% 
  dplyr::select(c(ATAC_header, "score")) %>% 
  mutate(group = ntile(score, n = group_num), group = as.factor(group), Condition = "Control")
atac_grange_ctrl <- GRanges(seqnames = atac_ctrl$chr, ranges = IRanges(start = atac_ctrl$start, end = atac_ctrl$end), group = atac_ctrl$group)

atac_604_all <- read.table("./data/scATACseq_preprocessing/ATAC_cov_average_hg19_1kb_bin.604_all.bed") %>% setNames(c("chr", "start", "end", "id", "score_604_all"))
atac_dis <- atac_604_all %>% 
  setNames(c("chr", "start", "end", "id", "score")) %>% 
  mutate(group = ntile(score, n = group_num), group = as.factor(group), Condition = "IHD")
atac_grange_dis <- GRanges(seqnames = atac_dis$chr, ranges = IRanges(start = atac_dis$start, end = atac_dis$end), group = atac_dis$group)

##### annotate the gene symbols
peakAnno_ctrl <- annotatePeak(atac_grange_ctrl, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, tssRegion = c(-2000, 2000), annoDb = "org.Hs.eg.db")
peakAnno_dis <- annotatePeak(atac_grange_dis, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, tssRegion = c(-2000, 2000), annoDb = "org.Hs.eg.db")

anno_ctrl_df <- as.data.frame(peakAnno_ctrl) %>% 
  left_join(atac_ctrl, by = c("seqnames" = "chr", "start" = "start", "end"   = "end", "group" = "group")) %>% 
  dplyr::select(score, group, SYMBOL)

anno_dis_df  <- as.data.frame(peakAnno_dis) %>% 
  left_join(atac_dis, by = c("seqnames" = "chr", "start" = "start", "end"   = "end", "group" = "group")) %>% 
  dplyr::select(score, group, SYMBOL)

atac_gene_ctrl <- anno_ctrl_df %>%
  group_by(SYMBOL) %>%
  summarise(accessibility = mean(score, na.rm = TRUE), decile = as.integer(median(as.numeric(names(which.max(table(group)))))), .groups = "drop") %>%
  mutate(Condition = "Control")

atac_gene_dis <- anno_dis_df %>%
  group_by(SYMBOL) %>%
  summarise(accessibility = mean(score, na.rm = TRUE), decile = as.integer(median(as.numeric(names(which.max(table(group)))))), .groups = "drop") %>%
  mutate(Condition = "IHD")

atac_gene_all <- bind_rows(atac_gene_ctrl, atac_gene_dis)

GO_genes_df <- bind_rows(
  data.frame(Gene = unique(unlist(strsplit(GO_BP_Ctrl$genes, ",\\s*"))), Condition = "Control"), 
  data.frame(Gene = unique(unlist(strsplit(GO_BP_IHD$genes, ",\\s*"))), Condition = "IHD")
)

# Merge GO gene info with expression deciles
go_atac_all <- atac_gene_all %>% 
  inner_join(GO_genes_df, by = c("SYMBOL" = "Gene", "Condition"))

# Get decile counts and expected uniform distribution
decile_counts <- go_atac_all %>% group_by(Condition, decile) %>% summarise(n = n(), .groups = "drop")
expected_counts <- go_atac_all %>% group_by(Condition) %>% summarise(expected = n() / 8)
decile_counts <- decile_counts %>% left_join(expected_counts, by = "Condition")

# Chi-square test per condition
pvals <- decile_counts %>% group_by(Condition) %>%
  summarise(chisq_p = chisq.test(n, p = rep(1/length(n), length(n)))$p.value)

p_barplot <- ggplot(decile_counts, aes(x = factor(decile), y = n, fill = Condition)) + 
  geom_bar(stat = "identity") + geom_hline(aes(yintercept = expected), linetype = "dashed", color = "black") + facet_wrap(~ Condition) +
  geom_text(data = pvals, aes(x = 4.5, y = max(decile_counts$n) * 0.95, label = paste0("p = ", signif(chisq_p, 3))), inherit.aes = FALSE, size = 6) + 
  scale_fill_manual(values = c("Control" = color_set[1], "IHD" = color_set[2]), guide = "legend") + 
  labs(x = "Chromatin accessibility level", y = "Number of GO Genes") + theme_linedraw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), 
        text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
ggsave("check_GO_genes_snRNA_ATACseq/GO_genes_snATACseq.pdf", plot = p_barplot, width = 10, height = 5, dpi = 600)
