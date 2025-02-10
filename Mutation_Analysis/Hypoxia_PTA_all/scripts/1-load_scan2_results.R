################################################################################
##################### load SCAN2 results and add metadata ######################
################################################################################
##### create directory for figures saving
sSNV_figure_dir <- paste0(project_dir, "/sSNV_figures")
dir.create(sSNV_figure_dir, showWarnings = FALSE)
main_figure_dir <- paste0(sSNV_figure_dir, "/main_figures")
dir.create(main_figure_dir, showWarnings = FALSE)
suppl_figure_dir <- paste0(sSNV_figure_dir, "/suppl_figures")
dir.create(suppl_figure_dir, showWarnings = FALSE)
other_figure_dir <- paste0(sSNV_figure_dir, "/other_figures")
dir.create(other_figure_dir, showWarnings = FALSE)
table_dir <- paste0(sSNV_figure_dir, "/tables")
dir.create(table_dir, showWarnings = FALSE)

ctrl_name <- "Control"
dis_name <- "IHD"

################################################################################
##### load SCAN2 vcf files for SNV and indel calls
##### for SNVs
snv_vcf_file_list <- list.files(path = paste0(project_dir, "/data/vcfs/vcfs_by_cell/ssnv_vcfs"), pattern = "ssnv", full.names = TRUE)
snv_base_vcf_names <- basename(snv_vcf_file_list)
snv_sample_names <- unlist(strsplit(snv_base_vcf_names, "\\."))[3 * (1:length(snv_base_vcf_names)) - 1]
snv_grl <- read_vcfs_as_granges(snv_vcf_file_list, snv_sample_names, ref_genome)
seqlengths_list <- seqlengths(Hsapiens)[1:22]
seqlengths(snv_grl) <- seqlengths_list
chromosomes <- seqnames(get(ref_genome))[1:22]

################################################################################
## for indels
indel_vcf_files_list <- list.files(path = paste0(project_dir, "/data/vcfs/vcfs_by_cell/sindel_vcfs"), pattern = "sindel", full.names = TRUE)
indel_base_vcf_names <- basename(indel_vcf_files_list)
indel_sample_names <- unlist(strsplit(indel_base_vcf_names, "\\."))[3 * (1:length(indel_base_vcf_names)) - 1]
indel_grl <- read_vcfs_as_granges(indel_vcf_files_list, indel_sample_names, ref_genome, type = "indel")

indel_grl <- get_indel_context(indel_grl, ref_genome)
# seqlengths(indel_grl) <- seqlengths_list
indel_counts <- count_indel_contexts(indel_grl)
list_01 <- setdiff(Cell_ID_list, c("1863_C1", "1940_B1"))
indel_counts <- indel_counts[ , list_01]

################################################################################
##### load SCAN2 burden results and metadata
metadata_df <- read.csv(paste0(project_dir, "/data/meta_data.csv"), header = TRUE) %>% group_by(factor(Condition, levels = c(ctrl_name, dis_name))) %>% arrange(Age, .by_group = TRUE)
metadata_df_AMG <- metadata_df[metadata_df$Age >= 40 & metadata_df$Age < 80,]
Cell_ID_list <- metadata_df$Cell_ID
Cell_ID_list_AMG <- metadata_df_AMG$Cell_ID

##### read MAPD, CoV, Coverage, and Depth and filter by combined QC and sensitivity
MAPD_CoV_df <- read.csv(paste0(project_dir, "/data/QC_and_confounding_factors/MAPD_CoV_summary.csv")) %>% arrange(match(Cell_ID, Cell_ID_list))
Coverage_Depth_df <- read.csv(paste0(project_dir, "/data/QC_and_confounding_factors/Coverage_and_Depth/Coverage_Depth_summary.csv")) %>% arrange(match(Cell_ID, Cell_ID_list))
all_QC_metrics <- merge(MAPD_CoV_df, Coverage_Depth_df) %>% 
  mutate(MAPD_norm = 1 - (MAPD / max(MAPD)), CoV_norm = 1 - (CoV / max(CoV)), Coverage_norm = Coverage / max(Coverage), Depth_norm = Depth / max(Depth), 
         Summary_Score = 0.2 * MAPD_norm + 0.2 * CoV_norm + 0.3 * Coverage_norm + 0.3 * Depth_norm) %>% 
  filter(Summary_Score > 0.5) %>% arrange(match(Cell_ID, Cell_ID_list))

##### compare original and revised scan2 results
source(paste0(project_dir, "/scripts/2-original_vs_revised_burden.R"))
input_dir <- paste0(project_dir, "/data/revised_scan2_results")
SCAN2_df <- original_revised_SCAN2_call_func(input_dir) %>% arrange(match(Cell_ID, Cell_ID_list)) %>% 
  mutate(Age = metadata_df$Age, Gender = metadata_df$Gender, Case_ID = metadata_df$Case_ID, Condition = metadata_df$Condition) %>% 
  filter(Cell_ID %in% all_QC_metrics$Cell_ID, snv.somatic.sens > 0.05) %>% 
  mutate(MAPD = all_QC_metrics$MAPD, CoV = all_QC_metrics$CoV, Coverage = all_QC_metrics$Coverage, Depth = all_QC_metrics$Depth)

##### filter the outliers in Control cells based on estimated burden
burden_df_ctrl <- SCAN2_df[, c("snv.rate.per.gb", "Cell_ID", "Case_ID", "Age", "Condition")] %>% filter(Condition == "Control")
burden_age_model_ctrl <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_ctrl, REML = FALSE)
outliers_cook <- burden_df_ctrl %>% mutate(cooksd = cooks.distance(burden_age_model_ctrl)) %>% filter(cooksd > 4 * mean(cooksd))
geom_line_data <- tibble(x_fit = range(burden_df_ctrl$Age), y_fit = x_fit * fixef(burden_age_model_ctrl)[2] + fixef(burden_age_model_ctrl)[1])
p_SNV_burden_filter <- ggplot(burden_df_ctrl, aes(x = Age, y = snv.rate.per.gb)) + 
  geom_point(size = 3) + geom_line(data = geom_line_data, aes(x = x_fit, y = y_fit), color = "dodgerblue3", linetype = "dashed") +
  geom_point(data = outliers_cook, aes(x = Age, y = snv.rate.per.gb), color = "red", size = 3) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
ggsave(paste0(other_figure_dir, "/2-SNV_burden_filter.pdf"), plot = p_SNV_burden_filter, width = 8, height = 6, dpi = 600)
SCAN2_df <- SCAN2_df %>% filter(!Cell_ID %in% outliers_cook$Cell_ID)

##### Customized color palettes
ctrl_range <- seq(min(which(SCAN2_df$Condition == ctrl_name)), max(which(SCAN2_df$Condition == ctrl_name)))
dis_range <- seq(min(which(SCAN2_df$Condition == dis_name)), max(which(SCAN2_df$Condition == dis_name)))
SCAN2_df_AMG <- SCAN2_df %>% filter(Cell_ID %in% Cell_ID_list_AMG)
ctrl_range_AMG <- seq(min(which(SCAN2_df_AMG$Condition == ctrl_name)), max(which(SCAN2_df_AMG$Condition == ctrl_name)))
dis_range_AMG <- seq(min(which(SCAN2_df_AMG$Condition == dis_name)), max(which(SCAN2_df_AMG$Condition == dis_name)))

metadata_df <- metadata_df %>% filter(Cell_ID %in% SCAN2_df$Cell_ID)
metadata_df_AMG <- metadata_df_AMG %>% filter(Cell_ID %in% SCAN2_df$Cell_ID)
Cell_ID_list <- metadata_df$Cell_ID
Cell_ID_list_AMG <- metadata_df_AMG$Cell_ID

ctrl_color_palette <- colorRampPalette(c("skyblue1","dodgerblue4"))
dis_color_palette <- colorRampPalette(c("pink1","firebrick3"))
ctrl_dis_color <- c(ctrl_color_palette(9)[7], dis_color_palette(4)[3])
SCAN2_df$Color <- c(SCAN2_df %>% filter(Condition == "Control") %>% count(Case_ID) %>% {rep(ctrl_color_palette(nrow(.)), .$n)}, 
                    SCAN2_df %>% filter(Condition == "IHD") %>% count(Case_ID) %>% {rep(dis_color_palette(nrow(.)), .$n)})

##### add metadata to mutation call vcfs
# genomic_SCAN2_df <- c()
# for (condition_temp in c(ctrl_name, dis_name)) {
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