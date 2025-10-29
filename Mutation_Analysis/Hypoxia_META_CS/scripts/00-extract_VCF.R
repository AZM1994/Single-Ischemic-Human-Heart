##### load all_final_calls
Condition_list <- c("Control", "IHD")
Ploidy_list <- c("2n", "4n")
Sensitivity_FNR_df <- data.frame(Cell_ID = character(), Sensitivity_FNR = numeric())
for (ploidy in Ploidy_list){
  for (condition_temp in Condition_list){
    vcf_output_dir <- paste0(project_dir, "/data/vcfs/", ploidy, "/ssnv_vcfs/", condition_temp)
    dir.create(paste0(vcf_output_dir, "/DS"), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(vcf_output_dir, "/SS"), recursive = TRUE, showWarnings = FALSE)
    final_calls_path_list_cond <- list.files(path = paste0(project_dir, "/data/all_final_calls/", ploidy, "/", condition_temp), pattern = "txt", full.names = TRUE)
    
    for (final_calls_path_temp in final_calls_path_list_cond){
      final_calls_txt_temp <- read.csv(final_calls_path_temp, header = FALSE, skip = 13, sep = c("\t", " "), na.strings = "")
      Sensitivity_FNR <- ifelse(
        any(final_calls_txt_temp$V1 == "Sensitivity_FNR"),
        final_calls_txt_temp$V4[final_calls_txt_temp$V1 == "Sensitivity_FNR"],
        NA
        )
      
      filename_base <- basename(final_calls_path_temp)
      prefix <- sub("_.*", "", filename_base)
      suffix <- sub(".*[^0-9]([0-9]+)\\.?mypileup.*", "\\1", filename_base)
      Cell_ID <- paste0(prefix, "_", ploidy, "_", suffix)
      
      if (is.na(Sensitivity_FNR)) {
        print(Cell_ID)
      }
      
      Sensitivity_FNR_df <- rbind(Sensitivity_FNR_df, cbind(Cell_ID, Sensitivity_FNR))
      
      DSB_calls <- final_calls_txt_temp %>% filter(V1 == "NV") %>% select(2:5) %>% 
        mutate(Cell_ID = Cell_ID, dot_col1 = ".", dot_col2 = ".", dot_col3 = ".") %>%
        select(V2, V3, dot_col1, V4, V5, dot_col2, dot_col3, Cell_ID)  
      SSB_calls <- final_calls_txt_temp %>% filter(V1 == "DV") %>% select(2:5) %>% 
        mutate(Cell_ID = Cell_ID, dot_col1 = ".", dot_col2 = ".", dot_col3 = ".") %>%
        select(V2, V3, dot_col1, V4, V5, dot_col2, dot_col3, Cell_ID)  
      
      vcf_header <- paste0("##fileformat=VCFv4.0\n#", paste(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), collapse = "\t"), "\n")
      
      vcf_file_path_DSB <- paste0(vcf_output_dir, "/DS/DS_ssnv.", Cell_ID, ".vcf")
      cat(vcf_header, file = vcf_file_path_DSB)
      write.table(DSB_calls, file = vcf_file_path_DSB, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
      
      vcf_file_path_SSB <- paste0(vcf_output_dir,"/SS/SS_ssnv.", Cell_ID, ".vcf")
      cat(vcf_header, file = vcf_file_path_SSB)
      write.table(SSB_calls, file = vcf_file_path_SSB, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
}
write.csv(Sensitivity_FNR_df, paste0(project_dir, "/data/Sensitivity_FNR.csv"), row.names = F)
