##### load metadata, sample information
##### vcf files for SCAN2 SNV and indel calls
control_name <- "Control"
disease_name <- "IHD"

metadata_df <- read.csv(paste0(project_dir, "/data/meta_data.csv"), header = TRUE) %>% group_by(factor(Condition, levels = c(control_name, disease_name))) %>% arrange(Age, .by_group = TRUE)
metadata_df_age_match <- metadata_df[metadata_df$Age >= 40 & metadata_df$Age < 80,]
Cell_ID_list <- metadata_df$Cell_ID
Cell_ID_list_age_match <- metadata_df_age_match$Cell_ID

##### for SNVs
snv_vcf_file_list <- list.files(path = paste0(project_dir, "/data/vcfs/vcfs_by_cell/ssnv_vcfs"), pattern = "ssnv", full.names = TRUE)
snv_base_vcf_names <- basename(snv_vcf_file_list)
snv_sample_names <- unlist(strsplit(snv_base_vcf_names, "\\."))[3 * (1:length(snv_base_vcf_names)) - 1]
snv_grl <- read_vcfs_as_granges(snv_vcf_file_list, snv_sample_names, ref_genome)
seqlengths_list <- seqlengths(Hsapiens)[1:22]
seqlengths(snv_grl) <- seqlengths_list
chromosomes <- seqnames(get(ref_genome))[1:22]

## for indels
# indel_vcf_files_list <- list.files(path = paste0(project_dir, "/data/vcfs/vcfs_by_cell/sindel_vcfs"), pattern = "sindel", full.names = TRUE)
# indel_base_vcf_names <- basename(indel_vcf_files_list)
# indel_sample_names <- unlist(strsplit(indel_base_vcf_names, "\\."))[3 * (1:length(indel_base_vcf_names)) - 1]
# indel_grl <- read_vcfs_as_granges(indel_vcf_files_list, indel_sample_names, ref_genome, type = "indel")

# indel_grl <- get_indel_context(indel_grl, ref_genome)
# # seqlengths(indel_grl) <- seqlengths_list
# indel_counts <- count_indel_contexts(indel_grl)
# indel_counts <- indel_counts[ , Cell_ID_list]