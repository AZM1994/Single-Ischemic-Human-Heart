################################################################################
##################### load SCAN2 results and add metadata ######################
################################################################################

##### load SCAN2 results
SCAN2_file_list <- list.files(path = paste0(project_dir, "/data/scan2_results"), pattern = ".csv", full.names = TRUE)
SCAN2_df <- rbindlist(sapply(SCAN2_file_list, fread, simplify = FALSE), idcol = 'cell_ID')
SCAN2_df <- SCAN2_df[SCAN2_df$V1 == 2]
SCAN2_df$cell_ID <- sub(".*\\\\([^.]*).*", "\\1", gsub(".*/(.+).csv*", "\\1", SCAN2_df$cell_ID))

##### create directory for figures saving
sSNV_figure_dir <- paste0(project_dir, "/sSNV_figures")
dir.create(sSNV_figure_dir)
main_figure_dir <- paste0(sSNV_figure_dir, "/main_figures")
dir.create(main_figure_dir)
suppl_figure_dir <- paste0(sSNV_figure_dir, "/suppl_figures")
dir.create(suppl_figure_dir)
other_figure_dir <- paste0(sSNV_figure_dir, "/other_figures")
dir.create(other_figure_dir)

##### compare original and revised scan2 results
source(paste0(project_dir, "/scripts/2-original_vs_revised_burden.R"))

##### read MAPD and CoV
MAPD_CoV_df <- read.csv(paste0(project_dir, "/data/QC_and_confounding_factors/MAPD/MAPD_CoV_summary.csv")) %>% arrange(match(Cell_ID, Cell_ID_list))

##### read Coverage and Depth
Coverage_Depth_df <- read.csv(paste0(project_dir, "/data/QC_and_confounding_factors/Coverage_and_Depth/Coverage_Depth_summary.csv")) %>% arrange(match(Cell_ID, Cell_ID_list))

##### add metadata
SCAN2_df <- revised_df %>% 
  arrange(match(cell_ID, Cell_ID_list)) %>% 
  mutate(age = metadata_df$Age) %>% 
  mutate(gender = metadata_df$Gender) %>%
  mutate(Case_ID = metadata_df$Case_ID) %>% 
  mutate(condition = metadata_df$Condition) %>% 
  mutate(MAPD = MAPD_CoV_df$MAPD) %>% 
  mutate(CoV = MAPD_CoV_df$CoV) %>% 
  mutate(Coverage = Coverage_Depth_df$Coverage) %>% 
  mutate(Depth = Coverage_Depth_df$Depth) 
  
control_range <- seq(min(which(SCAN2_df$condition == control_name)), max(which(SCAN2_df$condition == control_name)))
disease_range <- seq(min(which(SCAN2_df$condition == disease_name)), max(which(SCAN2_df$condition == disease_name)))

SCAN2_df_age_match <- SCAN2_df[SCAN2_df$cell_ID %in% Cell_ID_list_age_match]

control_range_age_match <- seq(min(which(SCAN2_df_age_match$condition == control_name)), 
                               max(which(SCAN2_df_age_match$condition == control_name)))
disease_range_age_match <- seq(min(which(SCAN2_df_age_match$condition == disease_name)), 
                               max(which(SCAN2_df_age_match$condition == disease_name)))
