##### load metadata, sample information
##### vcf files for SCAN2 SNV and indel calls
control_name <- "Control"
disease_name <- "IHD"

metadata_df <- read.csv(paste0(project_dir, "/data/meta_data.csv"), header = TRUE) %>% 
  group_by(factor(Condition, levels = c(control_name, disease_name))) %>% arrange(Age, .by_group = TRUE)
metadata_df_age_match <- metadata_df[metadata_df$Age >= 40 & metadata_df$Age < 80,]
Cell_ID_list <- metadata_df$Cell_ID
Cell_ID_list_age_match <- metadata_df_age_match$Cell_ID

##### for SNVs
snv_vcf_file_list <- list.files(path = paste0(project_dir, "/data/vcfs/vcfs_by_cell"), pattern = "ssnv_list", full.names = TRUE)
base_vcf_names <- basename(snv_vcf_file_list)
sample_names <- unlist(strsplit(base_vcf_names, "\\."))[3 * (1:length(base_vcf_names)) - 1]
grl <- read_vcfs_as_granges(snv_vcf_file_list, sample_names, ref_genome)
seqlengths_list <- seqlengths(Hsapiens)[1:22]
seqlengths(grl) <- seqlengths_list
chromosomes <- seqnames(get(ref_genome))[1:22]