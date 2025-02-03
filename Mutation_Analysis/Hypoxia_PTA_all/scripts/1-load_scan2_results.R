################################################################################
##################### load SCAN2 results and add metadata ######################
################################################################################

##### load SCAN2 results
SCAN2_file_list <- list.files(path = paste0(project_dir, "/data/scan2_results/csv"), pattern = ".csv", full.names = TRUE)
SCAN2_df <- SCAN2_file_list %>% sapply(fread, simplify = FALSE) %>% rbindlist(idcol = "Cell_ID") %>% .[V1 == 2] %>% .[, Cell_ID := sub(".*\\\\([^.]*).*", "\\1", gsub(".*/(.+).csv*", "\\1", Cell_ID))]

##### create directory for figures saving
sSNV_figure_dir <- paste0(project_dir, "/sSNV_figures")
dir.create(sSNV_figure_dir)
main_figure_dir <- paste0(sSNV_figure_dir, "/main_figures")
dir.create(main_figure_dir)
suppl_figure_dir <- paste0(sSNV_figure_dir, "/suppl_figures")
dir.create(suppl_figure_dir)
other_figure_dir <- paste0(sSNV_figure_dir, "/other_figures")
dir.create(other_figure_dir)
table_dir <- paste0(sSNV_figure_dir, "/tables")
dir.create(table_dir)

##### compare original and revised scan2 results
source(paste0(project_dir, "/scripts/2-original_vs_revised_burden.R"))

##### read MAPD, CoV, Coverage, and Depth
MAPD_CoV_df <- read.csv(paste0(project_dir, "/data/QC_and_confounding_factors/MAPD_CoV_summary.csv")) %>% arrange(match(Cell_ID, Cell_ID_list))
Coverage_Depth_df <- read.csv(paste0(project_dir, "/data/QC_and_confounding_factors/Coverage_and_Depth/Coverage_Depth_summary.csv")) %>% arrange(match(Cell_ID, Cell_ID_list))

all_QC_metrics <- merge(MAPD_CoV_df, Coverage_Depth_df) %>% 
  mutate(MAPD_norm = 1 - (MAPD / max(MAPD)), CoV_norm = 1 - (CoV / max(CoV)), Coverage_norm = Coverage / max(Coverage), Depth_norm = Depth / max(Depth),
  Summary_Score = 0.2 * MAPD_norm + 0.2 * CoV_norm + 0.3 * Coverage_norm + 0.3 * Depth_norm) %>% filter(Summary_Score > 0.5) %>% arrange(match(Cell_ID, Cell_ID_list))

revised_df <- revised_df %>% filter(Cell_ID %in% all_QC_metrics$Cell_ID)
metadata_df <- metadata_df %>% filter(Cell_ID %in% all_QC_metrics$Cell_ID)
metadata_df_age_match <- metadata_df[metadata_df$Age >= 40 & metadata_df$Age < 80,]
Cell_ID_list <- metadata_df$Cell_ID
Cell_ID_list_age_match <- metadata_df_age_match$Cell_ID

##### add metadata
SCAN2_df <- revised_df %>% arrange(match(Cell_ID, Cell_ID_list)) %>% 
  mutate(Age = metadata_df$Age, Gender = metadata_df$Gender, Case_ID = metadata_df$Case_ID, Condition = metadata_df$Condition, 
         MAPD = all_QC_metrics$MAPD, CoV = all_QC_metrics$CoV, Coverage = all_QC_metrics$Coverage, Depth = all_QC_metrics$Depth)
control_range <- seq(min(which(SCAN2_df$Condition == control_name)), max(which(SCAN2_df$Condition == control_name)))
disease_range <- seq(min(which(SCAN2_df$Condition == disease_name)), max(which(SCAN2_df$Condition == disease_name)))

SCAN2_df_age_match <- SCAN2_df[SCAN2_df$Cell_ID %in% Cell_ID_list_age_match]
control_range_age_match <- seq(min(which(SCAN2_df_age_match$Condition == control_name)), max(which(SCAN2_df_age_match$Condition == control_name)))
disease_range_age_match <- seq(min(which(SCAN2_df_age_match$Condition == disease_name)), max(which(SCAN2_df_age_match$Condition == disease_name)))

##### Customized color palettes
ctrl_color_palette <- colorRampPalette(c("skyblue1","dodgerblue4"))
dis_color_palette <- colorRampPalette(c("pink1","firebrick3"))
ctrl_dis_color <- c(ctrl_color_palette(9)[7], dis_color_palette(4)[3])
SCAN2_df$Color <- c(SCAN2_df %>% filter(Condition == "Control") %>% count(Case_ID) %>% {rep(ctrl_color_palette(nrow(.)), .$n)}, 
                    SCAN2_df %>% filter(Condition == "IHD") %>% count(Case_ID) %>% {rep(dis_color_palette(nrow(.)), .$n)})

##### add metadata to mutation call vcfs
# genomic_SCAN2_df <- c()
# for (condition_temp in c(control_name, disease_name)) {
#   # for (mutation_type in c("ssnv", "sindel")) {
#     for (mutation_type in c("ssnv")) {
#     cat("Get genomic context for", condition_temp, mutation_type, "...\n")
#     heart_PTA_Cases_vcf_temp <- data.frame(read.table(paste0(project_dir, "/data/annotation_results/all_age/", condition_temp, "/heart_PTA_all.all_age.", condition_temp, ".", mutation_type, ".vcf"), sep = "\t"))
#     genomic_context_temp <- read.csv(paste0(project_dir, "/data/annotation_results/all_age/", condition_temp, "/heart_PTA_all.all_age.", condition_temp, ".", mutation_type, ".csv"), header = TRUE) %>%
#       mutate(Cell_ID = heart_PTA_Cases_vcf_temp$V8, mutation_type, .before = 1)
# 
#     SCAN2_df_temp <- SCAN2_df[, c("Cell_ID", "Case_ID", "Condition", "Age", "Gender")]
#     genomic_SCAN2_df_temp <- SCAN2_df_temp %>% merge(genomic_context_temp)
#     # genomic_SCAN2_df_temp <- genomic_SCAN2_df_temp[, c("Cell_ID", "Case_ID", "Condition", "Age", "Gender", "mutation_type", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene")]
#     genomic_SCAN2_df <- rbind(genomic_SCAN2_df, genomic_SCAN2_df_temp)
#     write.csv(genomic_SCAN2_df_temp, paste0(project_dir, "/data/annotation_results/all_age/", "genomic_SCAN2_", condition_temp, "_", mutation_type, ".csv"), row.names = F)
#   }
# }