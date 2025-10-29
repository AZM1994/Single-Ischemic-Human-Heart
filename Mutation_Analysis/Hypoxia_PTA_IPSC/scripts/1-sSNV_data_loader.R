################################################################################
##################### load SCAN2 results and add metadata ######################
################################################################################

##### create directory for figures saving
sSNV_figure_dir <- paste0(project_dir, "/sSNV_figures")
dir.create(sSNV_figure_dir, showWarnings = FALSE)
main_figure_dir <- paste0(sSNV_figure_dir, "/1-main_figures")
dir.create(main_figure_dir, showWarnings = FALSE)
suppl_figure_dir <- paste0(sSNV_figure_dir, "/2-suppl_figures")
dir.create(suppl_figure_dir, showWarnings = FALSE)
other_figure_dir <- paste0(sSNV_figure_dir, "/3-other_figures")
dir.create(other_figure_dir, showWarnings = FALSE)
table_dir <- paste0(sSNV_figure_dir, "/4-table_rds")
dir.create(table_dir, showWarnings = FALSE)

ctrl_name <- "Normoxia"
dis_name <- "Hypoxia"

##### load SCAN2 vcf files for SNV calls
ssnv_vcf_file_list <- list.files(path = paste0(project_dir, "/data/vcfs/vcfs_by_cell"), pattern = "ssnv_list", full.names = TRUE)
base_vcf_names <- basename(ssnv_vcf_file_list)
sample_names <- unlist(strsplit(base_vcf_names, "\\."))[3*(1:length(base_vcf_names))-1]
ssnv_grl <- read_vcfs_as_granges(ssnv_vcf_file_list, sample_names, ref_genome)
seqlengths_list <- seqlengths(Hsapiens)[1:22]
seqlengths(ssnv_grl) <- seqlengths_list
chromosomes <- seqnames(get(ref_genome))[1:22]

##### load SCAN2 burden results and metadata
metadata_df <- read.csv(paste0(project_dir, "/data/meta_data.csv"), header = TRUE)
Cell_ID_list <- metadata_df$Cell_ID

##### read MAPD, CoV, Coverage, and Depth and filter by combined QC and sensitivity
MAPD_CoV_df <- read.csv(paste0(project_dir, "/data/QC_and_confounding_factors/MAPD_CoV_summary.csv")) %>% 
  arrange(match(Cell_ID, Cell_ID_list))
Coverage_Depth_df <- read.csv(paste0(project_dir, "/data/QC_and_confounding_factors/Coverage_and_Depth/Coverage_Depth_summary.csv")) %>% 
  arrange(match(Cell_ID, Cell_ID_list))
all_QC_metrics <- merge(MAPD_CoV_df, Coverage_Depth_df) %>% 
  mutate(MAPD_norm = 1 - (MAPD / max(MAPD)), CoV_norm = 1 - (CoV / max(CoV)), 
         Coverage_norm = Coverage / max(Coverage), Depth_norm = Depth / max(Depth), 
         Summary_Score = 0.2 * MAPD_norm + 0.2 * CoV_norm + 0.3 * Coverage_norm + 0.3 * Depth_norm) %>% 
  arrange(match(Cell_ID, Cell_ID_list))
all_QC_metrics_save <- all_QC_metrics %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(all_QC_metrics_save, paste0(table_dir, "/1-all_QC_metrics.csv"), row.names = F)
all_QC_metrics <- all_QC_metrics %>% filter(Summary_Score > 0.4)

##### compare original and revised scan2 results
scan2_results_file_list <- list.files(path = paste0(project_dir, "/data/scan2_results"), pattern = ".csv", full.names = TRUE)
sSNV_SCAN2_df <- rbindlist(sapply(scan2_results_file_list, fread, simplify = FALSE), idcol = 'Cell_ID')
sSNV_SCAN2_df <- sSNV_SCAN2_df[sSNV_SCAN2_df$V1 == 2]
sSNV_SCAN2_df$Cell_ID <- sub(".*\\\\([^.]*).*", "\\1", gsub(".*/(.+).csv*", "\\1", sSNV_SCAN2_df$Cell_ID))

source(paste0(project_dir, "/scripts/2-original_vs_revised_burden.R"))
sSNV_SCAN2_df <- revised_df

## add metadata
sSNV_SCAN2_df <- slice(sSNV_SCAN2_df, match(Cell_ID_list, sSNV_SCAN2_df$Cell_ID)) %>% 
  mutate(Condition = metadata_df$Condition) %>% 
  filter(Cell_ID %in% all_QC_metrics$Cell_ID, snv.somatic.sens > 0.05) %>% 
  mutate(MAPD = all_QC_metrics$MAPD, CoV = all_QC_metrics$CoV, 
         Coverage = all_QC_metrics$Coverage, Depth = all_QC_metrics$Depth)

##### customized parameters for nice plot
ctrl_color_palette <- colorRampPalette(c("skyblue1","dodgerblue4"))
dis_color_palette <- colorRampPalette(c("pink1","firebrick3"))
ctrl_dis_color <- c(ctrl_color_palette(9)[7], dis_color_palette(4)[3])
sSNV_SCAN2_df$Color <- rep(ctrl_dis_color, each = 9)
sSNV_SCAN2_df$Outline_Color <- ifelse(sSNV_SCAN2_df$Condition == ctrl_name, "#0D3F70", "firebrick4")

ctrl_range <- seq(min(which(sSNV_SCAN2_df$Condition == ctrl_name)), max(which(sSNV_SCAN2_df$Condition == ctrl_name)))
dis_range <- seq(min(which(sSNV_SCAN2_df$Condition == dis_name)), max(which(sSNV_SCAN2_df$Condition == dis_name)))

# saveRDS(sSNV_SCAN2_df, file = paste0(project_dir, "/data/PTA_IPSC_burden.rds"))
# saveRDS(metadata_df, file = paste0(project_dir, "/data/PTA_IPSC_metadata.rds"))
##### add metadata to mutation call vcfs
# genomic_PTA_IPSC_burden_df <- c()
# for (condition_temp in c(ctrl_name, dis_name)) {
#     for (mutation_type in c("ssnv")) {
#       cat("Get genomic context for", condition_temp, mutation_type, "...\n")
#       heart_PTA_IPSC_vcf_temp <- data.frame(read.table(paste0(project_dir, "/data/annotation_results/all_", mutation_type, "/", condition_temp, "/heart_PTA_IPSC.", condition_temp, "_", mutation_type, ".vcf"), sep = "\t")) %>% 
#         mutate(V8 = case_when(V8 == "NormoxiaiPSCA1" ~ "A1", V8 == "NormoxiaiPSCA2" ~ "A2", V8 == "NormoxiaiPSCA3" ~ "A3", V8 == "NormoxiaiPSCA6" ~ "A6", 
#                               V8 == "NormoxiaiPSCA7" ~ "A7", V8 == "NormoxiaiPSCA8" ~ "A8", V8 == "NormoxiaiPSCA9" ~ "A9", 
#                               V8 == "NormoxiaiPSCA11" ~ "A11", V8 == "Normoxia2n_IPSC" ~ "A2n", 
#                               V8 == "HypoxiaiPSCB3" ~ "B3", V8 == "HypoxiaiPSCB4" ~ "B4", V8 == "HypoxiaiPSCB7" ~ "B7", V8 == "HypoxiaiPSCB8" ~ "B8", 
#                               V8 == "HypoxiaiPSCB9" ~ "B9", V8 == "HypoxiaiPSCB10" ~ "B10", V8 == "HypoxiaiPSCB11" ~ "B11", 
#                               V8 == "B12" ~ "B12", V8 == "Hypoxia2n_IPSC" ~ "B2n", TRUE ~ V8))
#       genomic_context_temp <- read.csv(paste0(project_dir, "/data/annotation_results/all_", mutation_type, "/", condition_temp, "/heart_PTA_IPSC.", condition_temp, ".hg19_multianno.csv"), header = TRUE) %>%
#         mutate(Cell_ID = heart_PTA_IPSC_vcf_temp$V8, mutation_type, .before = 1)
# 
#       genomic_PTA_IPSC_burden_df_temp <- metadata_df %>% 
#         select(c("Cell_ID", "Condition")) %>% merge(genomic_context_temp) %>% arrange(factor(Cell_ID, levels = metadata_df$Cell_ID))
#       genomic_PTA_IPSC_burden_df <- rbind(genomic_PTA_IPSC_burden_df, genomic_PTA_IPSC_burden_df_temp)
#       write.csv(genomic_PTA_IPSC_burden_df_temp, paste0(project_dir, "/data/annotation_results/all_", mutation_type, "/", "genomic_context_", condition_temp, "_", mutation_type, ".csv"), row.names = F)
#     }
# }
