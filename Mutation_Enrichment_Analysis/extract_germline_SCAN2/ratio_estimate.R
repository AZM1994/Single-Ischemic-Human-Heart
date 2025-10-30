library(dplyr)
library(readr)
library(stringr)
setwd("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration/Mutation_Enrichment_Analysis/extract_germline_SCAN2")

vcf_dir <- "resampled_training_variants/vcf_annovar"
files <- list.files(vcf_dir, pattern = "_resampled_variants.hg19_multianno.csv$", full.names = TRUE)
genomic_context_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Cell_ID", "Func.refGene", "Gene.refGene")

## for exo_intro_ratio
exo_intro_ratio_df <- data.frame(Cell_ID = character(), exo_intro_germline_ratio = numeric(), stringsAsFactors = FALSE)
for (f in files) {
  Cell_ID <- str_replace(basename(f), "_resampled_variants.hg19_multianno.csv", "")
  message("Processing: ", Cell_ID)
  
  df <- read_csv(f, show_col_types = FALSE) %>% 
    filter(Func.refGene %in% c("exonic", "exonic;splicing", "intronic")) %>% 
    mutate(Func.refGene = gsub("exonic;splicing", "exonic", Func.refGene))
  
  counts <- table(df$Func.refGene)
  exonic_count <- counts["exonic"]
  intronic_count <- counts["intronic"]
  exo_intro_ratio <- ifelse(intronic_count > 0, exonic_count / intronic_count, NA)
  exo_intro_ratio_df <- rbind(exo_intro_ratio_df, data.frame(Cell_ID = Cell_ID, exo_intro_germline_ratio = exo_intro_ratio))
}
write.csv(exo_intro_ratio_df, "resampled_training_variants/exo_intro_ratio.csv", row.names = F)

## for dNdS ratio
dnds_df <- data.frame(Cell_ID = character(), germline_dNdS = numeric(), stringsAsFactors = FALSE)
for (f in files) {
  Cell_ID <- str_replace(basename(f), "_resampled_variants.hg19_multianno.csv", "")
  message("Processing: ", Cell_ID)
  
  df <- read_csv(f, show_col_types = FALSE) %>% 
    filter(Func.refGene %in% c("exonic", "exonic;splicing")) %>% 
    mutate(Func.refGene = gsub("exonic;splicing", "exonic", Func.refGene))
  
  dN <- sum(df$ExonicFunc.refGene %in% c("nonsynonymous SNV"), na.rm = TRUE)
  dS <- sum(df$ExonicFunc.refGene %in% c("synonymous SNV"), na.rm = TRUE)
  
  dNdS <- ifelse(dS > 0, dN / dS, NA)
  
  dnds_df <- rbind(dnds_df, data.frame(Cell_ID = Cell_ID, germline_dNdS = dNdS))
}
write.csv(dnds_df, "resampled_training_variants/dnds.csv", row.names = F)
