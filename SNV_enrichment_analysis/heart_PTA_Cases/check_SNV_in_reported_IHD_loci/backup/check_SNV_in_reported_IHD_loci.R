library(ggplot2)
library(reshape2)
library(stringr)
library(dplyr)
library(stats)
library(data.table)
library(ggpubr)
library(stringr)
library(melt)
library(ggforce)

# setwd("/Users/zhemingan/Documents/BCH_research/SNV_enrichment_analysis/heart_PTA_Cases")
##### read in normal and disease data
Control.SNV_gene_list <- read.csv("./data/annovar_annotation/normal/heart_PTA_Cases.all_age.normal.hg19_multianno.csv", header = TRUE)
heart_PTA_Cases_Control_vcf <- read.table("./data/annovar_annotation/normal/heart_PTA_Cases.all_age.normal_ssnv.vcf", sep = "\t")
Control.SNV_gene_list$Cell_ID <- heart_PTA_Cases_Control_vcf$V8
Control.SNV_gene_list$Case_ID <- str_extract(Control.SNV_gene_list$Cell_ID, "[^_]+")
Control.SNV_gene_list$Condition <- "Control"
Control.SNV_gene_list <- Control.SNV_gene_list[c("Chr", "Start", "End", "Ref", "Alt", "Cell_ID", "Case_ID", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "Condition")]
Control.SNV_gene_list_exonic <- Control.SNV_gene_list[Control.SNV_gene_list$Func.refGene %in% c("exonic", "exonic;splicing", "splicing"), ]
Control.SNV_gene_list_nonsynonymous <- Control.SNV_gene_list[Control.SNV_gene_list$ExonicFunc.refGene %in% c("nonsynonymous SNV", "stopgain"), ]

IHD.SNV_gene_list <- read.csv("./data/annovar_annotation/disease/heart_PTA_Cases.all_age.disease.hg19_multianno.csv", header = TRUE)
heart_PTA_Cases_IHD_vcf <- read.table("./data/annovar_annotation/disease/heart_PTA_Cases.all_age.disease_ssnv.vcf", sep = "\t")
IHD.SNV_gene_list$Cell_ID <- heart_PTA_Cases_IHD_vcf$V8
IHD.SNV_gene_list$Case_ID <- str_extract(IHD.SNV_gene_list$Cell_ID, "[^_]+")
IHD.SNV_gene_list$Condition <- "IHD"
IHD.SNV_gene_list <- IHD.SNV_gene_list[c("Chr", "Start", "End", "Ref", "Alt", "Cell_ID", "Case_ID", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "Condition")]
IHD.SNV_gene_list_exonic <- IHD.SNV_gene_list[IHD.SNV_gene_list$Func.refGene %in% c("exonic", "exonic;splicing", "splicing"), ]
IHD.SNV_gene_list_nonsynonymous <- IHD.SNV_gene_list[IHD.SNV_gene_list$ExonicFunc.refGene %in% c("nonsynonymous SNV", "stopgain"), ]

combined.SNV_gene_list <- rbind(Control.SNV_gene_list, IHD.SNV_gene_list)
combined.SNV_gene_list_exonic <- rbind(Control.SNV_gene_list_exonic, IHD.SNV_gene_list_exonic)
combined.SNV_gene_list_nonsynonymous <- rbind(Control.SNV_gene_list_nonsynonymous, IHD.SNV_gene_list_nonsynonymous)

reported_IHD_loci <- read.csv("./data/hot_gene_list_heart_pta_02.csv", header = FALSE)
reported_IHD_loci <- reported_IHD_loci[!duplicated(reported_IHD_loci$V1), ]

#####################################################################
##### split rows which have multiple genes and save to a matrix #####
combined.count_summary <- matrix(nrow = 1, ncol = length(reported_IHD_loci) + 1)
# combined.count_summary <- data.frame(t(c(reported_IHD_loci, "condition")))
# colnames(combined.count_summary) <- c(reported_IHD_loci, "condition")
# for (condition_index in c("Normal")){
# for (condition_index in c("Disease")){
for (condition_index in c("Control", "IHD")){
  input_df <- combined.SNV_gene_list[combined.SNV_gene_list$Condition == condition_index, ]
  # input_df <- combined.SNV_gene_list_exonic[combined.SNV_gene_list_exonic$condition == condition_index, ]
  # input_df <- combined.SNV_gene_list_nonsynonymous[combined.SNV_gene_list_nonsynonymous$condition == condition_index, ]
  
  split_strings <- strsplit(input_df$Gene.refGene, ";")
  max_split <- max(sapply(split_strings, length))
  split_columns <- do.call(rbind, lapply(split_strings, function(x) {
    length(x) <- max_split
    replace(x, is.na(x), " ")
  }))
  colnames(split_columns) <- paste("split", 1:max_split, sep = "")
  split_columns_df <- as.data.frame(split_columns)
  split_columns_df <- cbind(input_df$Cell_ID, split_columns_df)
  colnames(split_columns_df)[1] <- "Cell_ID"
  
  #####################################################################
  ##### count the number of SNVs in IHD_loci list (cell specific) #####
  count_summary <- matrix(nrow = length(reported_IHD_loci), ncol = 1)
  for (donor_index in unique(split_columns_df$Cell_ID)){
    print(donor_index)
    split_columns_temp <- split_columns_df[split_columns_df$Cell_ID == donor_index, 2 : (max_split + 1)]
    count <- table(factor(as.matrix(split_columns_temp), levels = reported_IHD_loci))
    count_summary <- cbind(count_summary, count)
  }
  count_summary <- count_summary[, -1]
  colnames(count_summary) <- unique(split_columns_df$Cell_ID)
  # sum(count_summary)
  count_summary <- data.frame(t(count_summary))
  # melt_count_summary <- melt(count_summary)
  count_summary$Condition <- condition_index
  combined.count_summary <- rbind(combined.count_summary, as.matrix(count_summary))
}
combined.count_summary <- as.data.frame(combined.count_summary[-1, ])

rownames(combined.count_summary) <- c("1039", "1864_E3", "4402_A1", "4402_A3", "4638_A1", "4638_A4", "4638_B2", "5657_D1",  "5657_H1", 
                                      "5828_C2", "5828_G2", "5919_C4", "5919_D2", "5919_E3", "5919_F6", "6032_A1", "6032_C7", "6032_E7",
                                      "604_B2", "604_B3", "604_B6", "1113_D1", "1113_E1", "1113_F1", "1363_A4", "1363_D4", "1363_H4", 
                                      "1743_A3", "1743_C3", "1743_F3", "1673_A2", "1673_A3", "1673_D2")
combined.count_summary$Cell_ID <- rownames(combined.count_summary)
combined.count_summary$Case_ID <- c("1039", "1864", "4402", "4402", "4638", "4638", "4638", "5657", "5657", 
                                    "5828", "5828", "5919", "5919", "5919", "5919", "6032", "6032", "6032",
                                    "604", "604", "604", "1113", "1113", "1113", "1363", "1363", "1363",
                                    "1743", "1743", "1743", "1673", "1673", "1673")
combined.count_summary_save <- combined.count_summary %>% 
  select((ncol(.)-2):ncol(.), everything()) %>% select(where(~ any(. != 0))) 

write.csv(combined.count_summary_save, "./results/SNV_count_in_reported_IHD_loci.csv", row.names = F)

melt_combined.count_summary <- melt(combined.count_summary[,-83], id.vars = c("Cell_ID", "Condition")) %>% 
  setNames(c("Cell_ID", "Condition", "IHD_loci", "SNV_count")) %>% filter(SNV_count > 0) %>% mutate(SNV_count = as.numeric(SNV_count))
label_context <- melt_combined.count_summary$SNV_count
p_by_cell <- ggplot(melt_combined.count_summary, aes(Cell_ID, IHD_loci, fill= SNV_count)) + 
  scale_fill_gradient(low = "dodgerblue", high = "firebrick3", na.value = "gray80", trans = 'log2', name = "sSNV count") +
  geom_tile(colour = "grey50") + geom_text(aes(label = label_context), size = 3, color = "white") + theme_linedraw() + labs(x = "Cell ID", y = "Reported IHD-associated loci") + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  facet_grid(. ~  factor(Condition, level = c("Control", "IHD")), scales = "free_x", space = "free_x")
ggsave("./results/cell_specific_SNV_count_in_reported_IHD_loci.pdf", plot = p_by_cell, width = 14, height = 10, dpi = 600)

melt_combined.count_summary <- melt(combined.count_summary[,-82], id.vars = c("Case_ID", "Condition")) %>% 
  setNames(c("Case_ID", "Condition", "IHD_loci", "SNV_count")) %>% filter(SNV_count > 0) %>% mutate(SNV_count = as.numeric(SNV_count))
sum_by_Case_ID <- melt_combined.count_summary %>% group_by(Case_ID, IHD_loci, Condition) %>% summarise(SNV_count = sum(SNV_count))
label_context <- sum_by_Case_ID$SNV_count
p_by_case <- ggplot(sum_by_Case_ID, aes(Case_ID, IHD_loci, fill= SNV_count)) + 
  scale_fill_gradient(low = "#6EA8FE", high = "#FF6E6E", na.value = "gray80", trans = 'log2', name = "sSNV count") +
  geom_tile(colour = "grey50") + geom_text(aes(label = label_context), size = 3, color = "white") + theme_linedraw() + labs(x = "Case ID", y = "Reported IHD-associated loci") + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + 
  facet_grid(. ~  factor(Condition, level = c("Control", "IHD")), scales = "free_x", space = "free_x")
ggsave("./results/donor_specific_SNV_count_in_reported_IHD_loci.pdf", plot = p_by_case, width = 9, height = 10, dpi = 600)
