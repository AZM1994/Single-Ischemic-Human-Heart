working_dir <- setwd("/Users/zhemingan/My Drive/BCH/Huang Lab/Heart_hypoxia/Mutation_Analysis/Hypoxia_PTA_Cases/QC_and_confounding_factors/Coverage_and_Depth/")
coverage_files <- list.files(path = paste0(working_dir, "/results"), pattern = "coverage", full.names = TRUE)
all_coverage_names <- basename(coverage_files)
# all_sample_names <- unlist(strsplit(all_coverage_names, "\\_\\."))[3*(1:length(all_coverage_names))-1]

coverage_depth_summary <- data.frame()
for (i in coverage_files){
  coverage_df <- read.table(i)
  coverage <- sum(coverage_df$V3[1:22] * coverage_df$V6[1:22] / 100) / sum(coverage_df$V3[1:22])
  depth <- mean(coverage_df$V7[1:22])
  row_temp <- (c(basename(i), round(coverage, digits = 3), round(depth, digits = 3)))
  coverage_depth_summary <- rbind(coverage_depth_summary, row_temp)
  # print(c(basename(i), round(coverage, digits = 3), round(depth, digits = 3)))
}
colnames(coverage_depth_summary) <- c('Cell_ID', 'Coverage', 'Depth')

coverage_depth_summary_df <- as.data.frame(coverage_depth_summary) %>% 
  mutate(Cell_ID = sub("^[^_]*_(.*)\\..*", "\\1", Cell_ID)) %>% 
  mutate(Cell_ID = gsub("_Cases", "", Cell_ID)) %>% 
  mutate(Cell_ID = ifelse(Cell_ID == "1864_M-E3-2n", "1864_E3", ifelse(Cell_ID == "4402_1_A1-2n", "4402_A1", ifelse(Cell_ID == "4402_1_A3-2n", "4402_A3", 
  ifelse(Cell_ID == "4638_1_A1-2n", "4638_A1", ifelse(Cell_ID == "4638_1_A4-2n", "4638_A4", ifelse(Cell_ID == "4638_1_B2-2n", "4638_B2", 
  ifelse(Cell_ID == "5657_M-D1-2n", "5657_D1", ifelse(Cell_ID == "5657_M-H1-2n", "5657_H1", ifelse(Cell_ID == "5828_CM_C2_2n", "5828_C2", 
  ifelse(Cell_ID == "5828_CM_G2_2n", "5828_G2", ifelse(Cell_ID == "5919_1_C4-2n", "5919_C4", ifelse(Cell_ID == "5919_1_D2-2n", "5919_D2", 
  ifelse(Cell_ID == "5919_1_E3-2n", "5919_E3", ifelse(Cell_ID == "5919_1_F6-2n", "5919_F6", ifelse(Cell_ID == "6032_1_A1-2n", "6032_A1", 
  ifelse(Cell_ID == "6032_M-C7-2n", "6032_C7", ifelse(Cell_ID == "6032_M-E7-2n", "6032_E7", ifelse(Cell_ID == "1673_CM_A2_2n", "1673_A2", 
  ifelse(Cell_ID == "1673_CM_A3_2n", "1673_A3", ifelse(Cell_ID == "1673_CM_D2_2n", "1673_D2", Cell_ID)))))))))))))))))))))

write.table(coverage_depth_summary_df, paste0(working_dir, "/Coverage_Depth_summary.csv"), sep = ",", quote = F, row.name = F)
