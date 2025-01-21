#############################################################################
################ compare the original and revised SCAN2 call ################
#############################################################################

#############################################################################
##### read revised snv and indel burden tables
input_dir <- paste0(project_dir, "/data/revised_scan2_results")

refined_burden_table_func <- function(input_path, version_type, mutation_type, post_fix){
  burden_file_list <- list.files(path = paste0(input_path, "/", version_type, "/", mutation_type), pattern = ".csv", full.names = TRUE)
  integrated_df <- rbindlist(sapply(burden_file_list, fread, simplify = FALSE), idcol = 'cell_ID')
  integrated_df <- integrated_df[integrated_df$V1 == 2]
  
  refined_integrated_df <- integrated_df %>% 
    mutate(cell_ID = sub(".*\\\\([^.]*).*", "\\1", gsub(".*/(.+).csv*", "\\1", cell_ID))) %>% mutate(cell_ID = gsub(paste0(mutation_type, '_gatk_table_'), "", cell_ID)) %>% 
    mutate(cell_ID = gsub(paste0('_', post_fix), "", cell_ID)) %>% mutate(cell_ID = gsub('1_', "", cell_ID)) %>% mutate(cell_ID = gsub('-2n', "", cell_ID)) %>% 
    mutate(cell_ID = gsub('M-', "", cell_ID)) %>% mutate(cell_ID = gsub('_2n', "", cell_ID)) %>% mutate(cell_ID = gsub('CM_', "", cell_ID))
  
  return(refined_integrated_df)
}

revised_snv_df <- refined_burden_table_func(input_dir, "revised", "snv", "revised")
revised_indel_df <- refined_burden_table_func(input_dir, "revised", "indel", "revised")
original_snv_df <- refined_burden_table_func(input_dir, "original", "snv", "scan2")
original_indel_df <- refined_burden_table_func(input_dir, "original", "indel", "scan2")

##### integrate snv and indel burden table
original_df <- cbind(original_snv_df, original_indel_df[,-(1:2)])
colnames(original_df) <- colnames(SCAN2_df)
revised_df <- cbind(revised_snv_df, revised_indel_df[,-(1:2)])
colnames(revised_df) <- colnames(SCAN2_df)
revised_df <- revised_df %>% mutate(age = metadata_df$Age) %>% mutate(gender = metadata_df$Gender) %>% mutate(Case_ID = metadata_df$Case_ID)

snv_matrix <- as.data.frame(cbind(original_df$snv.rate.per.gb, revised_df$snv.rate.per.gb)) %>% setNames(c("original", "revised"))
indel_matrix <- as.data.frame(cbind(original_df$indel.rate.per.gb, revised_df$indel.rate.per.gb)) %>% setNames(c("original", "revised")) %>% 
  mutate(original = ifelse(original <= 0, NA, original)) %>% mutate(revised = ifelse(revised <= 0, NA, revised))

geom_line_data_snv <- as.data.frame(matrix(c(10, 3500, 10, 3500), nrow = 2, ncol = 2))
##### scatter plot compare the original and revised SCAN2 results
p_original_vs_revised_burden_snv <- ggplot(snv_matrix, aes(x = original, y = revised)) + 
  geom_point(size = 1, na.rm = TRUE) + geom_line(aes(x = V1, y = V2), size = 0.5, colour = "red", data = geom_line_data_snv) + 
  scale_x_continuous("original SCAN2 sSNV/GB", transform = 'log10') + scale_y_continuous("revised SCAN2 sSNV/GB", transform = 'log10') + 
  theme_linedraw() + coord_fixed(ratio = 1) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), 
        panel.grid.minor = element_blank(), panel.border = element_rect(size = 0.5), text = element_text(size = 12))
ggsave(paste0(suppl_figure_dir, "/2-original_vs_revised_burden_snv.pdf"), plot = p_original_vs_revised_burden_snv, width = 3, height = 3, units = "in", dpi = 600)

revised_df_subset <- revised_df[c(20,21,28,29,30),]
revised_df_subset$data_source <- c("recovered missing burden", "modified SCAN2 burden", "modified SCAN2 burden", "recovered missing burden", "modified SCAN2 burden")
p_recovered_missing_burden_snv <- ggplot(revised_df_subset, aes(x = factor(cell_ID, level = Cell_ID_list), y = snv.rate.per.gb, shape = data_source, color = data_source)) + 
  geom_point(size = 3, na.rm = TRUE, show.legend = TRUE) + labs(x = "Cell ID", y = "sSNV rate per GB") + 
  scale_color_manual(values = c("recovered missing burden" = "green", "modified SCAN2 burden" = "black"), guide = "legend") + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), 
        panel.grid.minor = element_blank(), panel.border = element_rect(size = 0.5), text = element_text(size = 12), 
        legend.position=c(0.3, 0.8), legend.background = element_rect(fill = "grey90"))
ggsave(paste0(suppl_figure_dir, "/2-recovered_missing_burden_snv.pdf"), plot = p_recovered_missing_burden_snv, width = 5, height = 3, dpi = 600)

# geom_line_data_indel <- as.data.frame(matrix(c(3, 160, 3, 160), nrow = 2, ncol = 2))
# p_original_vs_revised_burden_indel <- ggplot(indel_matrix, aes(x = original, y = revised)) +
#   geom_point(size = 3, na.rm = TRUE) + 
#   geom_line(aes(x = V1, y = V2), colour = "#FC4E07", data = geom_line_data_indel) + 
#   scale_x_continuous("original sindel rate per GB", transform = 'log10') + scale_y_continuous("revised sindel rate per GB", transform = 'log10') +
#   theme_classic() + theme(text = element_text(size=24), axis.text.x = element_text(hjust = 0.5), panel.background = element_rect(fill = "white"))
# ggsave(paste0(suppl_figure_dir, "/2-original_vs_revised_burden_indel.pdf"), plot = p_original_vs_revised_burden_indel, width = 12, height = 8, dpi = 600)
