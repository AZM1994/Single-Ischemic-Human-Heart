library("dndscv")
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)

setwd("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration/Mutation_Enrichment_Analysis")
result_dir <- "dNdS"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
group_labels <- c("Age-matched Control", "IHD")
Ctrl_IHD_color <- c("#2D6EA8", "#DD555B")
genomic_factor <- 5.845
ctrl_name <- "Control"
dis_name <- "IHD"
age_line_middle <- 30
age_line_senior <- 70

##### genomic_context for Control, IHD, and germline
SCAN2_df <- readRDS("data/1-sSNV_SCAN2_df_filtered.rds") %>% rbind(list("germline_bulk", 2, 0, 0, 0, 0, 0, 0, 0, 0, "FALSE", 999, "bulk", "bulk", "bulk", "bulk", 0, 0, 0, 0, FALSE, FALSE))
SCAN2_df_filtered <- SCAN2_df %>% dplyr::select(Cell_ID, snv.rate.per.gb, Condition, Age) %>%
  mutate(Group = case_when(Condition == ctrl_name & Age >= age_line_middle & Age < age_line_senior ~ group_labels[1],
                           Condition == dis_name ~ group_labels[2], Cell_ID == "germline_bulk" ~ "Germline", TRUE ~ NA_character_)) %>%
  mutate(Group = factor(Group, levels = c(group_labels, "Germline"))) %>% drop_na(Group)

heart_PTA_all_Control_vcf <- read.table("data/annotation_results/all_age_sSNV/Control/heart_PTA_all.all_age.Control_ssnv.vcf", sep = "\t")
heart_PTA_all_IHD_vcf <- read.table("data/annotation_results/all_age_sSNV/IHD/heart_PTA_all.all_age.IHD_ssnv.vcf", sep = "\t")

genomic_context_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Cell_ID", "Func.refGene", "Gene.refGene")

genomic_context_ctrl <- read.csv("data/annotation_results/all_age_sSNV/Control/heart_PTA_all.all_age.Control_ssnv.hg19_multianno.csv", header = TRUE) %>% 
  mutate(Gene.refGene.split = strsplit(as.character(Gene.refGene), ";"), 
         dist_vals = str_extract_all(GeneDetail.refGene, "(?<=dist=)[0-9]+"), 
         dist_numeric = lapply(dist_vals, as.numeric)) %>% 
  rowwise() %>% 
  mutate(closest_idx = ifelse(length(dist_numeric) > 0, which.min(dist_numeric), NA_integer_), 
         Gene.refGene = ifelse(!is.na(closest_idx), Gene.refGene.split[[closest_idx]], Gene.refGene), 
         ClosestDist = ifelse(!is.na(closest_idx), dist_numeric[closest_idx], NA), 
         Func.refGene = case_when(Func.refGene == "upstream;downstream" & closest_idx == 1 ~ "upstream", 
                                  Func.refGene == "upstream;downstream" & closest_idx == 2 ~ "downstream", TRUE ~ Func.refGene), 
         Func.refGene = gsub("exonic;splicing", "exonic", Func.refGene), Func.refGene = gsub("UTR5;UTR3", "UTR5", Func.refGene)) %>% 
  ungroup() %>% 
  mutate(Cell_ID = heart_PTA_all_Control_vcf$V8) %>% dplyr::select(all_of(genomic_context_colnames))

genomic_context_dis <- read.csv("data/annotation_results/all_age_sSNV/IHD/heart_PTA_all.all_age.IHD_ssnv.hg19_multianno.csv", header = TRUE) %>% 
  mutate(Gene.refGene.split = strsplit(as.character(Gene.refGene), ";"), 
         dist_vals = str_extract_all(GeneDetail.refGene, "(?<=dist=)[0-9]+"), 
         dist_numeric = lapply(dist_vals, as.numeric)) %>% 
  rowwise() %>% 
  mutate(closest_idx = ifelse(length(dist_numeric) > 0, which.min(dist_numeric), NA_integer_), 
         Gene.refGene = ifelse(!is.na(closest_idx), Gene.refGene.split[[closest_idx]], Gene.refGene), 
         ClosestDist = ifelse(!is.na(closest_idx), dist_numeric[closest_idx], NA), 
         Func.refGene = case_when(Func.refGene == "upstream;downstream" & closest_idx == 1 ~ "upstream", 
                                  Func.refGene == "upstream;downstream" & closest_idx == 2 ~ "downstream", TRUE ~ Func.refGene), 
         Func.refGene = gsub("exonic;splicing", "exonic", Func.refGene), Func.refGene = gsub("UTR5;UTR3", "UTR5", Func.refGene)) %>% 
  ungroup() %>% 
  mutate(Cell_ID = heart_PTA_all_IHD_vcf$V8) %>% dplyr::select(all_of(genomic_context_colnames))

genomic_context_germline <- read.delim("data/heart_germline.genomic_context.random200000.tsv", header = F) %>% 
  filter(V8 != "upstream;downstream") %>% 
  rename_with(~ c("Chr","Start","End","Ref","Alt","Cell_ID","Context","Func.refGene","Gene.refGene"), everything()) |> 
  dplyr::select(all_of(genomic_context_colnames))

genomic_context <- rbind(genomic_context_ctrl, genomic_context_dis, genomic_context_germline)
genomic_SCAN2_df <- genomic_context |> inner_join(SCAN2_df %>% dplyr::select(Cell_ID, Case_ID, Age, Gender, Condition), by = "Cell_ID") %>% 
  mutate(Group = case_when(Condition == ctrl_name & Age >= age_line_middle & Age < age_line_senior ~ group_labels[1],
                           Condition == dis_name ~ group_labels[2], Cell_ID == "germline_bulk" ~ "Germline", TRUE ~ NA_character_)) %>%
  mutate(Group = factor(Group, levels = c(group_labels, "Germline"))) %>% drop_na(Group)

# sorted_counts <- sort(table(genomic_SCAN2_df$Gene.refGene), decreasing = TRUE)[1:20]

SNV_ctrl <- genomic_SCAN2_df[genomic_SCAN2_df$Condition == "Control", ] %>% 
  dplyr::select(c("Cell_ID", "Chr", "Start", "Ref", "Alt")) %>% 
  setNames(c("sampleID", "chr", "pos", "ref", "mut"))
SNV_dis <- genomic_SCAN2_df[genomic_SCAN2_df$Condition == "IHD", ] %>% 
  dplyr::select(c("Cell_ID", "Chr", "Start", "Ref", "Alt")) %>% 
  setNames(c("sampleID", "chr", "pos", "ref", "mut"))

clean_and_deduplicate_snv <- function(snv_df) {
  # Remove exact duplicate mutations within a sample
  snv_df <- snv_df[!duplicated(snv_df[, c("sampleID", "chr", "pos", "ref", "mut")]), ]
  
  # Remove adjacent (contiguous) mutations within the same sample
  snv_df <- snv_df[order(snv_df$sampleID, snv_df$chr, snv_df$pos), ]
  snv_df$prev_pos <- c(NA, snv_df$pos[-nrow(snv_df)])
  snv_df$prev_sample <- c(NA, snv_df$sampleID[-nrow(snv_df)])
  snv_df$prev_chr <- c(NA, snv_df$chr[-nrow(snv_df)])
  non_adjacent <- with(snv_df, 
                       !(sampleID == prev_sample & chr == prev_chr & abs(pos - prev_pos) == 1))
  snv_df <- snv_df[non_adjacent | is.na(non_adjacent), ]
  snv_df$prev_pos <- NULL
  snv_df$prev_sample <- NULL
  snv_df$prev_chr <- NULL
  
  # Remove shared mutations across samples
  # snv_df <- snv_df[!duplicated(snv_df[, c("chr", "pos", "ref", "mut")]), ]
  return(snv_df)
}

SNV_ctrl_clean <- clean_and_deduplicate_snv(SNV_ctrl)
SNV_dis_clean <- clean_and_deduplicate_snv(SNV_dis)

# Now rerun dndscv
dndsout_ctrl <- dndscv(SNV_ctrl_clean, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
dndsout_dis  <- dndscv(SNV_dis_clean, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)

# hot_gene_list <- read.csv("data/hot_gene_list_heart_pta_02.csv", header = FALSE)
# hot_gene_list <- hot_gene_list[!duplicated(hot_gene_list$V1), ]
# dndsout_ctrl <- dndscv(SNV_ctrl, gene_list = hot_gene_list, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
# dndsout_dis <- dndscv(SNV_dis, gene_list = hot_gene_list, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)

############ Extract significant positively selected genes (q value)
sel_genes_ctrl <- dndsout_ctrl$sel_cv %>% 
  mutate(Condition = "Control", Type = "Missense", Selection = ifelse(wmis_cv > 1, "Positive", "Negative")) %>% dplyr::rename(gene = 1) 
pos_ctrl <- sel_genes_ctrl %>% 
  filter((qmis_cv < 0.05 & wmis_cv > 1) | (qtrunc_cv < 0.05 & (wnon_cv > 1 | wspl_cv > 1)) | (qallsubs_cv < 0.05 & (wmis_cv > 1 | wnon_cv > 1 | wspl_cv > 1)))
neg_ctrl <- sel_genes_ctrl %>% 
  filter((qmis_cv < 0.05 & wmis_cv < 1) | (qtrunc_cv < 0.05 & wnon_cv < 1 & wspl_cv < 1) | (qallsubs_cv < 0.05 & wmis_cv < 1 & wnon_cv < 1 & wspl_cv < 1))
sel_genes_combined_ctrl <- bind_rows(pos_ctrl, neg_ctrl)

sel_genes_dis <- dndsout_dis$sel_cv %>% 
  mutate(Condition = "IHD", Type = "Missense", Selection = ifelse(wmis_cv > 1, "Positive", "Negative")) %>% dplyr::rename(gene = 1) 
pos_dis <- sel_genes_dis %>%
  filter((qmis_cv < 0.05 & wmis_cv > 1) | (qtrunc_cv < 0.05 & (wnon_cv > 1 | wspl_cv > 1)) | (qallsubs_cv < 0.05 & (wmis_cv > 1 | wnon_cv > 1 | wspl_cv > 1)))
neg_dis <- sel_genes_dis %>% 
  filter((qmis_cv < 0.05 & wmis_cv < 1) | (qtrunc_cv < 0.05 & wnon_cv < 1 & wspl_cv < 1) | (qallsubs_cv < 0.05 & wmis_cv < 1 & wnon_cv < 1 & wspl_cv < 1))
sel_genes_combined_dis <- bind_rows(pos_dis, neg_dis)

sig_selected_genes_df <- bind_rows(sel_genes_combined_ctrl, sel_genes_combined_dis)
sig_selected_genes_df_02 <- bind_rows(sel_genes_ctrl, sel_genes_dis) %>% 
  mutate(Condition = factor(Condition, levels = c("Control", "IHD")), 
         Selection = factor(Selection, levels = c("Positive", "Negative")), 
         Significance = if_else(qmis_cv < 0.05, "Significant", "Not Significant"))
sig_gene_with_CI_ctrl <- geneci(dndsout_ctrl, gene_list = sel_genes_combined_ctrl$gene)
sig_gene_with_CI_dis <- geneci(dndsout_dis, gene_list = sel_genes_combined_dis$gene)
sig_gene_with_CI_all <- bind_rows(sig_gene_with_CI_ctrl, sig_gene_with_CI_dis) %>% 
  merge(sig_selected_genes_df) %>% 
  mutate(Condition = factor(Condition, levels = c("Control", "IHD")), Selection = factor(Selection, levels = c("Positive", "Negative")))
write.csv(sig_gene_with_CI_all, paste0(result_dir, "/significant_selected_genes.csv"))
write.csv(sig_selected_genes_df_02, paste0(result_dir, "/dNdScv_analysis.csv"))

# genes under significant positive selection for missense mutations
p_sig <- ggplot(sig_selected_genes_df_02, aes(x = gene, y = -log10(qmis_cv), fill = Condition, color = Condition, shape = Selection)) + 
  geom_point(size = 2) + geom_hline(yintercept = -log10(0.05), color = "gray40", linetype = "22") + geom_vline(xintercept = 1, color = "#A1B84C") + 
  geom_label_repel(data = sig_selected_genes_df, aes(label = gene), nudge_x = c(-1.0, -1.17, -1.43, 0.75), nudge_y = c(0.08, 0.08, 0.1, -0.1), fill = NA, color = "black", size = 3, label.size = 0.5) + 
  scale_fill_manual(values = Ctrl_IHD_color) + scale_color_manual(values = Ctrl_IHD_color) + labs(x = "Gene", y = "-log10(FDR-adjusted P-value)") + theme_linedraw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), 
        text = element_text(size = 18), axis.title.x = element_text(hjust = 0.5, vjust = 0), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(paste0(result_dir, "/dNdScv_significant_genes.pdf"), plot = p_sig, width = 5, height = 5.4)
