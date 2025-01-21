library(stringr)
library(readxl)
library(dplyr)
library(data.table)

working_dir <- setwd("/Users/zhemingan/My Drive/BCH/Huang Lab/Heart_hypoxia/Mutation_Analysis/Hypoxia_PTA_Cases/QC_and_confounding_factors/MAPD/")

##### extract MAPD
MAPD_file_list <- list.files(path = paste0(working_dir, "/MAPD_results"), 
                             pattern = ".hg19.50k.k50.varbin.mapd.txt", full.names = TRUE, recursive = TRUE)
MAPD_df <- rbindlist(sapply(MAPD_file_list, fread, simplify = FALSE), idcol = 'Cell_ID') %>% 
  mutate(Cell_ID = sub(".*/([^/]+)/[^/]+", "\\1", Cell_ID)) %>% 
  mutate(Cell_ID = ifelse(Cell_ID == "1864_M-E3-2n", "1864_E3", ifelse(Cell_ID == "4402_1_A1-2n", "4402_A1", ifelse(Cell_ID == "4402_1_A3-2n", "4402_A3", 
  ifelse(Cell_ID == "4638_1_A1-2n", "4638_A1", ifelse(Cell_ID == "4638_1_A4-2n", "4638_A4", ifelse(Cell_ID == "4638_1_B2-2n", "4638_B2", 
  ifelse(Cell_ID == "5657_M-D1-2n", "5657_D1", ifelse(Cell_ID == "5657_M-H1-2n", "5657_H1", ifelse(Cell_ID == "5828_CM_C2_2n", "5828_C2", 
  ifelse(Cell_ID == "5828_CM_G2_2n", "5828_G2", ifelse(Cell_ID == "5919_1_C4-2n", "5919_C4", ifelse(Cell_ID == "5919_1_D2-2n", "5919_D2", 
  ifelse(Cell_ID == "5919_1_E3-2n", "5919_E3", ifelse(Cell_ID == "5919_1_F6-2n", "5919_F6", ifelse(Cell_ID == "6032_1_A1-2n", "6032_A1", 
  ifelse(Cell_ID == "6032_M-C7-2n", "6032_C7", ifelse(Cell_ID == "6032_M-E7-2n", "6032_E7", ifelse(Cell_ID == "1673_CM_A2_2n", "1673_A2", 
  ifelse(Cell_ID == "1673_CM_A3_2n", "1673_A3", ifelse(Cell_ID == "1673_CM_D2_2n", "1673_D2", Cell_ID))))))))))))))))))))) %>% 
  as.data.frame() |> base::`[`(c("Cell_ID", "MAPD")) %>% filter(row_number() %% 25 == 0)

##### extract CoV
CoV_file_list <- list.files(path = paste0(working_dir, "/MAPD_results"), 
                             pattern = "hg19.50k.k50.varbin.normBinStat.txt", full.names = TRUE, recursive = TRUE)
CoV_df <- rbindlist(sapply(CoV_file_list, fread, simplify = FALSE), idcol = 'Cell_ID') %>% 
  mutate(Cell_ID = sub(".*/([^/]+)/[^/]+", "\\1", Cell_ID)) %>% 
  mutate(Cell_ID = ifelse(Cell_ID == "1864_M-E3-2n", "1864_E3", ifelse(Cell_ID == "4402_1_A1-2n", "4402_A1", ifelse(Cell_ID == "4402_1_A3-2n", "4402_A3", 
  ifelse(Cell_ID == "4638_1_A1-2n", "4638_A1", ifelse(Cell_ID == "4638_1_A4-2n", "4638_A4", ifelse(Cell_ID == "4638_1_B2-2n", "4638_B2", 
  ifelse(Cell_ID == "5657_M-D1-2n", "5657_D1", ifelse(Cell_ID == "5657_M-H1-2n", "5657_H1", ifelse(Cell_ID == "5828_CM_C2_2n", "5828_C2", 
  ifelse(Cell_ID == "5828_CM_G2_2n", "5828_G2", ifelse(Cell_ID == "5919_1_C4-2n", "5919_C4", ifelse(Cell_ID == "5919_1_D2-2n", "5919_D2", 
  ifelse(Cell_ID == "5919_1_E3-2n", "5919_E3", ifelse(Cell_ID == "5919_1_F6-2n", "5919_F6", ifelse(Cell_ID == "6032_1_A1-2n", "6032_A1", 
  ifelse(Cell_ID == "6032_M-C7-2n", "6032_C7", ifelse(Cell_ID == "6032_M-E7-2n", "6032_E7", ifelse(Cell_ID == "1673_CM_A2_2n", "1673_A2", 
  ifelse(Cell_ID == "1673_CM_A3_2n", "1673_A3", ifelse(Cell_ID == "1673_CM_D2_2n", "1673_D2", Cell_ID))))))))))))))))))))) %>% 
  as.data.frame() |> base::`[`(c("Cell_ID", "CoV")) %>% filter(row_number() %% 25 == 0)

MAPD_CoV_summary <- merge(MAPD_df, CoV_df)
write.table(MAPD_CoV_summary, paste0(working_dir, "/MAPD_CoV_summary.csv"), sep = ",", quote = F, row.name = F)
