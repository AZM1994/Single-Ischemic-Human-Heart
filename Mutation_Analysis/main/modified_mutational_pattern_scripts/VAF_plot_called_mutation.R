library(MutationalPatterns)
library(BSgenome)
library(readr)
library(ggplot2)

# head(available.genomes())
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
# ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
# library(ref_genome, character.only = TRUE)
setwd("/Users/zhemingan/My Drive/BCH/Huang Lab/Heart_hypoxia/scan2/PTA_IPSC/Extract_Scan2_Results/IPSC_Scan2_Analysis")

####### read all snv and sindel vcfs and metadata
metadata_df <- read.csv("meta_data.csv", header = TRUE)
snv_vcf_files <- list.files(path = "all_vcfs", pattern = "ssnv_list", full.names = TRUE)
sindel_vcf_files <- list.files(path = "all_vcfs", pattern = "sindel_list", full.names = TRUE)
all_vcf_names <- basename(snv_vcf_files)
sample_names <- unlist(strsplit(all_vcf_names, "\\."))[3*(1:length(all_vcf_names))-1]

raw_sample_name <- c('A1','A2','NormoxiaiPSCA3','NormoxiaiPSCA6','NormoxiaiPSCA7',
                     'NormoxiaiPSCA8','NormoxiaiPSCA9','NormoxiaiPSCA11','Normoxia2n_IPSC',
                     'HypoxiaiPSCB3','B4','HypoxiaiPSCB7','HypoxiaiPSCB8','HypoxiaiPSCB9',
                     'HypoxiaiPSCB10','HypoxiaiPSCB11','B12','Hypoxia2n_IPSC')
raw_sample_name_01 <- c('Hypoxia2n_IPSC')

snv_VCF_df <- NULL
for (f1 in snv_vcf_files) {
  snv_data_temp <- read.table(f1)
  snv_VCF_df <- rbind(snv_VCF_df, snv_data_temp)
}
snv_VCF_df$mut_type <- "ssnv"

sindel_VCF_df <- NULL
for (f2 in sindel_vcf_files) {
  sindel_data_temp <- read.table(f2)
  sindel_VCF_df <- rbind(sindel_VCF_df, sindel_data_temp)
}
sindel_VCF_df$mut_type <- "sindel"

all_VCF_df <- rbind(snv_VCF_df, sindel_VCF_df)

pdf("mutation_rate/vaf_plot_called_mutation.pdf", width = 15, height = 10)
for (raw_name_temp in raw_sample_name){
  df_temp <- all_VCF_df[all_VCF_df$V8 == raw_name_temp,]
  df_temp$idu <- 1:nrow(df_temp)
  df_temp$vaf <- df_temp$V7 / (df_temp$V6 + df_temp$V7)
  df_temp$depth <- df_temp$V6 + df_temp$V7
  number_of_ssnv <-  nrow(df_temp[df_temp$mut_type == 'ssnv', ])
  number_of_sindel <-  nrow(df_temp[df_temp$mut_type == 'sindel', ])
  total_mutation <- number_of_ssnv + number_of_sindel
  # vaf_called <- df_temp$V7 / (df_temp$V6 + df_temp$V7)
  p_i <- ggplot(df_temp, aes(x=V1, y = vaf, color = mut_type, shape = mut_type)) + 
    theme(text = element_text(size=20)) +
    scale_color_manual(values = c("ssnv" = "red", "sindel" = "blue")) + 
    scale_x_continuous(breaks = 1:22) +
    geom_point(size = 3) + 
    geom_hline(yintercept = mean(df_temp$vaf), color="blue") +
    annotate("text", x = 0, y = mean(df_temp$vaf) +0.01, label = round(mean(df_temp$vaf), digits = 2)) +
    labs(x = "Chromosome", y = "VAF") + 
    ggtitle(paste(raw_name_temp, ', ssnv=', number_of_ssnv, ', sindel=', number_of_sindel, ', total=', total_mutation)) + 
  # print(p_i)
  # plot(1:length(vaf_called), vaf_called, main=raw_name_temp,
  #      xlab="Called_mutation_index", ylab="VAF", pch=19)
}
plot_layout(ncol = 2)
dev.off()

pdf("mutation_rate/depth_plot_called_mutation.pdf", width = 15, height = 10)
for (raw_name_temp in raw_sample_name){
  df_temp <- all_VCF_df[all_VCF_df$V8 == raw_name_temp,]
  df_temp$idu <- 1:nrow(df_temp)
  df_temp$vaf <- df_temp$V7 / (df_temp$V6 + df_temp$V7)
  df_temp$depth <- df_temp$V6 + df_temp$V7
  # vaf_called <- df_temp$V7 / (df_temp$V6 + df_temp$V7)
  p_i <- ggplot(df_temp, aes(x=V1, y = depth, color = mut_type, shape = mut_type)) + 
    theme(text = element_text(size=20)) +
    scale_color_manual(values = c("ssnv" = "red", "sindel" = "blue")) + 
    scale_x_continuous(breaks = 1:22) +
    geom_point(size = 3) + 
    geom_hline(yintercept = mean(df_temp$depth), color="blue") +
    annotate("text", x = 0, y = mean(df_temp$depth) + 1, label = round(mean(df_temp$depth), digits = 2)) +
    ggtitle(raw_name_temp) +
    labs(x = "Chromosome", y = "Depth") + ggtitle(raw_name_temp)
  print(p_i)
  # plot(1:length(vaf_called), vaf_called, main=raw_name_temp,
  #      xlab="Called_mutation_index", ylab="VAF", pch=19)
}
dev.off()

pdf("mutation_rate/hist_vaf_plot_called_mutation.pdf", width = 15, height = 10)
for (raw_name_temp in raw_sample_name){
  df_temp <- all_VCF_df[all_VCF_df$V8 == raw_name_temp,]
  df_temp$idu <- 1:nrow(df_temp)
  df_temp$vaf <- df_temp$V7 / (df_temp$V6 + df_temp$V7)
  df_temp$depth <- df_temp$V6 + df_temp$V7
  # vaf_called <- df_temp$V7 / (df_temp$V6 + df_temp$V7)
  hist(df_temp$vaf, main = raw_name_temp, xlab="VAF")
  # plot(1:length(vaf_called), vaf_called, main=raw_name_temp,
  #      xlab="Called_mutation_index", ylab="VAF", pch=19)
}
dev.off()

pdf("mutation_rate/hist_depth_plot_called_mutation.pdf", width = 15, height = 10)
for (raw_name_temp in raw_sample_name){
  df_temp <- all_VCF_df[all_VCF_df$V8 == raw_name_temp,]
  df_temp$idu <- 1:nrow(df_temp)
  df_temp$vaf <- df_temp$V7 / (df_temp$V6 + df_temp$V7)
  df_temp$depth <- df_temp$V6 + df_temp$V7
  # vaf_called <- df_temp$V7 / (df_temp$V6 + df_temp$V7)
  hist(log10(df_temp$depth), main = raw_name_temp, xlab="log10(depth)")
  # plot(1:length(vaf_called), vaf_called, main=raw_name_temp,
  #      xlab="Called_mutation_index", ylab="VAF", pch=19)
}
dev.off()

