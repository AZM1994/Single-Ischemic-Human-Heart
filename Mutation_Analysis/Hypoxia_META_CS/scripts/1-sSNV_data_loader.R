################################################################################
##################### load SCAN2 results and add metadata ######################
factor_hg37 <- 5.845001

##### create directory for figures saving
ploidy_temp <- "2n"
# strand_temp <- "DS"
sSNV_figure_dir <- paste0(project_dir, "/sSNV_figures")
dir.create(sSNV_figure_dir, recursive = T, showWarnings = FALSE)
main_figure_dir <- paste0(sSNV_figure_dir, "/1-main_figures")
dir.create(main_figure_dir, recursive = T, showWarnings = FALSE)
suppl_figure_dir <- paste0(sSNV_figure_dir, "/2-suppl_figures")
dir.create(suppl_figure_dir, recursive = T, showWarnings = FALSE)
other_figure_dir <- paste0(sSNV_figure_dir, "/3-other_figures")
dir.create(other_figure_dir, recursive = T, showWarnings = FALSE)
table_dir <- paste0(sSNV_figure_dir, "/4-tables")
dir.create(table_dir, recursive = T, showWarnings = FALSE)

ctrl_name <- "Control"
dis_name <- "IHD"

##### load metadata
metadata_df <- read.csv(paste0(project_dir, "/data/meta_data.csv"), header = TRUE)
Cell_ID_list <- metadata_df$Cell_ID

##### load vcf files for SNV calls
ssnv_vcf_file_list_DS <- list.files(path = paste0(project_dir, "/data/vcfs/ssnv_vcfs/all_DS"), pattern = "vcf", full.names = TRUE)
ssnv_base_vcf_names_DS <- basename(ssnv_vcf_file_list_DS)
ssnv_sample_names_DS <- unlist(strsplit(ssnv_base_vcf_names_DS, "\\."))[3 * (1:length(ssnv_base_vcf_names_DS)) - 1]
ssnv_grl_DS <- read_vcfs_as_granges(ssnv_vcf_file_list_DS, ssnv_sample_names_DS, ref_genome)
seqlengths_list <- seqlengths(Hsapiens)[1:22]
seqlengths_list <- seqlengths_list[levels(seqnames(unlist(ssnv_grl_DS)))]
seqlengths(ssnv_grl_DS) <- seqlengths_list
chromosomes <- seqnames(get(ref_genome))[1:22]

##### read all QC metrics (MAPD, CoV, Coverage, Depth, Duplex coverage)
all_QC_metrics <- read.csv(paste0(project_dir, "/data/QC/all_META_CS_QC.csv")) %>% arrange(match(Cell_ID, Cell_ID_list)) %>% 
  mutate(Ploidy = sub("^[^_]+_([^_]+)_.*$", "\\1", Cell_ID)) %>% 
  filter(Ploidy == ploidy_temp)

##### burden calculation and filter by Summary_Score
col_name_order <- c("Cell_ID", "Case_ID", "Condition", "Age", "Gender", "Ploidy", "raw_mut_Count", "Duplex_Coverage", "MAPD", "CoV", "Coverage", "Depth")

META_CS_burden_df <- data.frame(Cell_ID = names(ssnv_grl_DS), raw_mut_Count = lengths(ssnv_grl_DS)) %>% 
  inner_join(metadata_df, by = c("Cell_ID")) %>% 
  inner_join(all_QC_metrics, by = c("Cell_ID", "Ploidy")) %>% 
  slice(match(metadata_df$Cell_ID, Cell_ID)) %>% 
  dplyr::select(all_of(col_name_order)) %>% 
  filter(Duplex_Coverage != 0) %>%
  mutate(MAPD_norm = 1 - (MAPD / max(MAPD)), CoV_norm = 1 - (CoV / max(CoV)), 
         Coverage_norm = Coverage / max(Coverage), Depth_norm = Depth / max(Depth), Duplex_Coverage_norm = Duplex_Coverage / max(Duplex_Coverage), 
         Summary_Score = 0.25 * MAPD_norm + 0.25 * CoV_norm + 0.20 * Coverage_norm + 0.20 * Depth_norm + 0.10 * Duplex_Coverage_norm) %>% 
  filter(Summary_Score > 0.50)

a <- data.frame(Cell_ID = names(ssnv_grl_DS), raw_mut_Count = lengths(ssnv_grl_DS)) %>% 
  inner_join(metadata_df, by = c("Cell_ID")) %>% 
  inner_join(all_QC_metrics, by = c("Cell_ID", "Ploidy")) %>% 
  slice(match(metadata_df$Cell_ID, Cell_ID)) %>% 
  dplyr::select(all_of(col_name_order))

META_CS_burden_df_save <- a %>% left_join(META_CS_burden_df, by = col_name_order) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(META_CS_burden_df_save, paste0(table_dir, "/1-all_QC_metrics_META_CS.csv"), row.names = F)

##### filter the outliers in Control cells based on estimated burden
# burden_df_ctrl <- META_CS_burden_df[, c("raw_mut_Count", "Cell_ID", "Case_ID", "Age", "Condition", "Ploidy")] %>% filter(Condition == "Control")
# burden_age_model_ctrl <- lmer(raw_mut_Count ~ Age + (1|Case_ID), burden_df_ctrl, REML = FALSE)
# outliers_cook <- burden_df_ctrl %>% mutate(cooksd = cooks.distance(burden_age_model_ctrl)) %>% filter(cooksd > quantile(cooksd, 0.95))
# geom_line_data <- tibble(x_fit = range(burden_df_ctrl$Age), y_fit = x_fit * fixef(burden_age_model_ctrl)[2] + fixef(burden_age_model_ctrl)[1])
# p_SNV_burden_filter <- ggplot(burden_df_ctrl, aes(x = Age, y = raw_mut_Count)) +
#   geom_point(size = 3) + geom_line(data = geom_line_data, aes(x = x_fit, y = y_fit), color = "dodgerblue3", linetype = "dashed") +
#   geom_point(data = outliers_cook, aes(x = Age, y = raw_mut_Count), color = "red", size = 3) + theme_linedraw() +
#   theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(),
#         panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
# ggsave(paste0(other_figure_dir, "/2-sSNV_burden_filter.pdf"), plot = p_SNV_burden_filter, width = 8, height = 6, dpi = 600)

# outliers_df <- META_CS_burden_df %>% filter(Cell_ID %in% outliers_cook$Cell_ID)
# META_CS_burden_df <- META_CS_burden_df %>% filter(!Cell_ID %in% outliers_cook$Cell_ID)

##### Customized color palettes
ctrl_range <- seq(min(which(META_CS_burden_df$Condition == ctrl_name)), max(which(META_CS_burden_df$Condition == ctrl_name)))
dis_range <- seq(min(which(META_CS_burden_df$Condition == dis_name)), max(which(META_CS_burden_df$Condition == dis_name)))
META_CS_burden_df_AMG <- META_CS_burden_df %>% filter(Age >= 30 & Age <= 70)
ctrl_range_AMG <- seq(min(which(META_CS_burden_df_AMG$Condition == ctrl_name)), max(which(META_CS_burden_df_AMG$Condition == ctrl_name)))
dis_range_AMG <- seq(min(which(META_CS_burden_df_AMG$Condition == dis_name)), max(which(META_CS_burden_df_AMG$Condition == dis_name)))

metadata_df <- metadata_df %>% filter(Cell_ID %in% META_CS_burden_df$Cell_ID)
metadata_df_AMG <- metadata_df %>% filter(Cell_ID %in% META_CS_burden_df_AMG$Cell_ID)
Cell_ID_list <- metadata_df$Cell_ID
Cell_ID_list_AMG <- metadata_df_AMG$Cell_ID
Case_ID_list <- as.character(unique(metadata_df$Case_ID))
Case_ID_list_AMG <- as.character(unique(metadata_df_AMG$Case_ID))

ctrl_color_palette <- colorRampPalette(c("skyblue1","dodgerblue4"))
dis1_color_palette <- colorRampPalette(c("pink1","firebrick3"))
ctrl_dis_color <- c(ctrl_color_palette(9)[7], dis1_color_palette(4)[3])
META_CS_burden_df$Color <- c(META_CS_burden_df %>% filter(Condition == "Control") %>% count(Case_ID) %>% {rep(ctrl_color_palette(nrow(.)), .$n)}, 
                             META_CS_burden_df %>% filter(Condition == "IHD") %>% count(Case_ID) %>% {rep(dis1_color_palette(nrow(.)), .$n)})
META_CS_burden_df$Outline_Color <- ifelse(META_CS_burden_df$Condition == "Control", "#0D3F70", "firebrick4")

# saveRDS(META_CS_burden_df, file = paste0(table_dir, "/1-META_CS_burden_df_filtered.rds"))
# saveRDS(metadata_df, file = paste0(table_dir, "/1-metadata_df_filtered.rds"))
##### add metadata to mutation call vcfs
# genomic_META_CS_burden_df <- c()
# for (condition_temp in c(ctrl_name, dis_name)) {
#     for (mutation_type in c("ssnv")) {
#       cat("Get genomic context for", condition_temp, mutation_type, "...\n")
#       heart_META_CS_vcf_temp <- data.frame(read.table(paste0(project_dir, "/data/annotation_results/all_age_", mutation_type, "/", condition_temp, "/heart_META_CS.all_age.", condition_temp, "_", mutation_type, ".vcf"), sep = "\t"))
#       genomic_context_temp <- read.csv(paste0(project_dir, "/data/annotation_results/all_age_", mutation_type, "/", condition_temp, "/heart_META_CS.all_age.", condition_temp, "_", mutation_type, ".hg19_multianno.csv"), header = TRUE) %>%
#         mutate(Cell_ID = heart_META_CS_vcf_temp$V8, mutation_type, .before = 1)
# 
#       genomic_META_CS_burden_df_temp <- read.csv(paste0(project_dir, "/data/meta_data.csv"), header = TRUE) %>% 
#         select(c("Cell_ID", "Case_ID", "Condition", "Age", "Gender")) %>% 
#         merge(genomic_context_temp)
#       genomic_META_CS_burden_df <- rbind(genomic_META_CS_burden_df, genomic_META_CS_burden_df_temp)
#       write.csv(genomic_META_CS_burden_df_temp, paste0(project_dir, "/data/annotation_results/all_age_", mutation_type, "/", "genomic_context_", condition_temp, "_", mutation_type, ".csv"), row.names = F)
#     }
# }
