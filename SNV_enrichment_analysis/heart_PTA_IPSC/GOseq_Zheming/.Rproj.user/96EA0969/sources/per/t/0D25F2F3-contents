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

setwd("/Users/zhemingan/Downloads/annovar/heart_PTA_Cases")
wd_dir <- getwd()
##### read in normal and disease data
Disease.SNV_gene_list <- read.csv(paste0(wd_dir, "/disease_results/gene_list/", "heart_PTA_Cases_Disease.SNV_gene_list.csv"), header = TRUE)
heart_PTA_Cases_Disease_vcf <- read.table("heart_PTA_Cases.all_disease_ssnv.vcf", sep = "\t")
Disease.SNV_gene_list$Cell_ID <- heart_PTA_Cases_Disease_vcf$V8
Disease.SNV_gene_list$Case_ID <- str_extract(Disease.SNV_gene_list$Cell_ID, "[^_]+")
Disease.SNV_gene_list$condition <- "Disease"
Disease.SNV_gene_list <- Disease.SNV_gene_list[c("Chr", "Start", "End", "Ref", "Alt", "Cell_ID", "Case_ID", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "condition")]
Disease.SNV_gene_list_exonic <- Disease.SNV_gene_list[Disease.SNV_gene_list$Func.refGene %in% c("exonic", "exonic;splicing", "splicing"), ]
Disease.SNV_gene_list_nonsynonymous <- Disease.SNV_gene_list[Disease.SNV_gene_list$ExonicFunc.refGene %in% c("nonsynonymous SNV", "stopgain"), ]

Normal.SNV_gene_list <- read.csv(paste0(wd_dir, "/normal_results/gene_list/", "heart_PTA_Cases_Normal.SNV_gene_list.csv"), header = TRUE)
heart_PTA_Cases_Normal_vcf <- read.table("heart_PTA_Cases.all_normal_ssnv.vcf", sep = "\t")
Normal.SNV_gene_list$Cell_ID <- heart_PTA_Cases_Normal_vcf$V8
Normal.SNV_gene_list$Case_ID <- str_extract(Normal.SNV_gene_list$Cell_ID, "[^_]+")
Normal.SNV_gene_list$condition <- "Normal"
Normal.SNV_gene_list <- Normal.SNV_gene_list[c("Chr", "Start", "End", "Ref", "Alt", "Cell_ID", "Case_ID", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "condition")]
Normal.SNV_gene_list_exonic <- Normal.SNV_gene_list[Normal.SNV_gene_list$Func.refGene %in% c("exonic", "exonic;splicing", "splicing"), ]
Normal.SNV_gene_list_nonsynonymous <- Normal.SNV_gene_list[Normal.SNV_gene_list$ExonicFunc.refGene %in% c("nonsynonymous SNV", "stopgain"), ]

combined.SNV_gene_list <- rbind(Normal.SNV_gene_list, Disease.SNV_gene_list)
combined.SNV_gene_list_exonic <- rbind(Normal.SNV_gene_list_exonic, Disease.SNV_gene_list_exonic)
combined.SNV_gene_list_nonsynonymous <- rbind(Normal.SNV_gene_list_nonsynonymous, Disease.SNV_gene_list_nonsynonymous)

hot_gene_list <- read.csv(paste0(wd_dir, "/check_SNV_in_hotgenes/hot_gene_list_heart_pta_02.csv"), header = FALSE)
hot_gene_list <- hot_gene_list[!duplicated(hot_gene_list$V1), ]

#####################################################################
##### split rows which have multiple genes and save to a matrix #####
combined.count_summary <- matrix(nrow = 1, ncol = length(hot_gene_list) + 1)
# combined.count_summary <- data.frame(t(c(hot_gene_list, "condition")))
# colnames(combined.count_summary) <- c(hot_gene_list, "condition")
# for (condition_index in c("Normal")){
# for (condition_index in c("Disease")){
for (condition_index in c("Normal", "Disease")){
  # input_df <- combined.SNV_gene_list[combined.SNV_gene_list$condition == condition_index, ]
  # input_df <- combined.SNV_gene_list_exonic[combined.SNV_gene_list_exonic$condition == condition_index, ]
  input_df <- combined.SNV_gene_list_nonsynonymous[combined.SNV_gene_list_nonsynonymous$condition == condition_index, ]
  
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
  ##### count the number of SNVs in hotgene list (cell specific) #####
  count_summary <- matrix(nrow = length(hot_gene_list), ncol = 1)
  for (donor_index in unique(split_columns_df$Cell_ID)){
    print(donor_index)
    split_columns_temp <- split_columns_df[split_columns_df$Cell_ID == donor_index, 2 : (max_split + 1)]
    count <- table(factor(as.matrix(split_columns_temp), levels = hot_gene_list))
    count_summary <- cbind(count_summary, count)
  }
  count_summary <- count_summary[, -1]
  colnames(count_summary) <- unique(split_columns_df$Cell_ID)
  # sum(count_summary)
  count_summary <- data.frame(t(count_summary))
  # melt_count_summary <- melt(count_summary)
  count_summary$condition <- condition_index
  combined.count_summary <- rbind(combined.count_summary, as.matrix(count_summary))
}
combined.count_summary <- as.data.frame(combined.count_summary[-1, ])
rownames(combined.count_summary) <- c("1039", "5828_C2", "5828_G2", "5919_C4", "5919_D2", "5919_E3", "5919_F6",
                                      "604_B2", "604_B3", "604_B6", "1113_D1", "1113_E1", "1113_F1", "1363_A4", "1363_D4", "1363_H4", 
                                      "1743_A3", "1743_C3", "1743_F3", "1673_A2", "1673_A3", "1673_D2")
combined.count_summary$Cell_ID <- rownames(combined.count_summary)
combined.count_summary$Case_ID <- c("1039", "5828", "5828", "5919", "5919", "5919", "5919",
                                    "604", "604", "604", "1113", "1113", "1113", "1363", "1363", "1363",
                                    "1743", "1743", "1743", "1673", "1673", "1673")
write.csv(combined.count_summary, paste0(wd_dir, "/check_SNV_in_hotgenes/SNV_count_in_hotgene_list.csv"))

pdf(paste0(wd_dir, "/check_SNV_in_hotgenes/cell_specific_SNV_count_in_hotgene_list.pdf"), width = 15, height = 20)
melt_combined.count_summary <- melt(combined.count_summary[,-83], id.vars = c("Cell_ID", "condition"))
colnames(melt_combined.count_summary) <- c("Cell_ID", "condition", "hotgene", "SNV_count")
melt_combined.count_summary$SNV_count <- gsub(" 0", "0", melt_combined.count_summary$SNV_count)
melt_combined.count_summary$SNV_count <- as.numeric(melt_combined.count_summary$SNV_count)
melt_combined.count_summary$SNV_count[melt_combined.count_summary$SNV_count == 0] <- NA
# unique(melt_combined.count_summary$SNV_count)
label_context <- melt_combined.count_summary$SNV_count
# label_context[label_context == 0] <- NA
# my_limits <- c(0, max(melt_combined.count_summary$SNV_count))
# my_breaks <- seq(0, max(melt_combined.count_summary$SNV_count), by = 2)
p_combined.count_summary_cell <- ggplot(melt_combined.count_summary, aes(Cell_ID, hotgene, fill= SNV_count)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = "gray80", 
                      # limits = my_limits, breaks = my_breaks, 
                      trans = 'log2', name = "SNV count") +
  geom_tile(colour = "grey50") + 
  geom_text(aes(label = label_context), size = 3, color = "white") +
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.5)) + 
  facet_grid(. ~  factor(condition, level = c("Normal", "Disease")), scale = "free_x")
print(p_combined.count_summary_cell)
dev.off()

pdf(paste0(wd_dir, "/check_SNV_in_hotgenes/donor_specific_SNV_count_in_hotgene_list.pdf"), width = 15, height = 20)
melt_combined.count_summary <- melt(combined.count_summary[,-82], id.vars = c("Case_ID", "condition"))
colnames(melt_combined.count_summary) <- c("Case_ID", "condition", "hotgene", "SNV_count")
melt_combined.count_summary$SNV_count <- gsub(" 0", "0", melt_combined.count_summary$SNV_count)
melt_combined.count_summary$SNV_count <- as.numeric(melt_combined.count_summary$SNV_count)
sum_by_Case_ID <- melt_combined.count_summary %>%
  group_by(Case_ID, hotgene, condition) %>%
  summarise(SNV_count = sum(SNV_count))

sum_by_Case_ID$SNV_count[sum_by_Case_ID$SNV_count == 0] <- NA
# unique(melt_combined.count_summary$SNV_count)
label_context <- sum_by_Case_ID$SNV_count
# label_context[label_context == 0] <- NA
# my_limits <- c(0, max(melt_combined.count_summary$SNV_count))
# my_breaks <- seq(0, max(melt_combined.count_summary$SNV_count), by = 2)
p_combined.count_summary_case <- ggplot(sum_by_Case_ID, aes(Case_ID, hotgene, fill= SNV_count)) + 
  scale_fill_gradient(low = "blue", high = "red", na.value = "gray80", 
                      # limits = my_limits, breaks = my_breaks, 
                      trans = 'log2', name = "SNV count") +
  geom_tile(colour = "grey50") + 
  geom_text(aes(label = label_context), size = 3, color = "white") +
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.5)) + 
  facet_grid(. ~  factor(condition, level = c("Normal", "Disease")), scale = "free_x")
print(p_combined.count_summary_case)

dev.off()

