################################################################################
################# compare the original and revised SCAN2 call ##################

original_revised_SCAN2_call_func <- function (input_dir, metadata_df) {
  #############################################################################
  ##### read revised snv and indel burden tables
  refined_burden_table_func <- function(input_path, version_type, mutation_type, post_fix){
    burden_file_list <- list.files(path = paste0(input_path, "/", version_type, "/", mutation_type), pattern = ".csv", full.names = TRUE)
    refined_integrated_df <- burden_file_list %>% sapply(fread, simplify = FALSE) %>% rbindlist(idcol = "Cell_ID") %>% .[V1 == 2] %>% 
      mutate(Cell_ID = sub(".*\\\\([^.]*).*", "\\1", gsub(".*/(.+).csv*", "\\1", Cell_ID)), Cell_ID = gsub(paste0(mutation_type, '_gatk_table_'), "", Cell_ID), 
             Cell_ID = gsub(paste0('_', post_fix), "", Cell_ID)) %>% rename_with(~ paste(mutation_type, ., sep = "."), -1:-2)
    return(refined_integrated_df)
  }

  revised_snv_df <- refined_burden_table_func(input_dir, "revised", "snv", "revised") %>% arrange(match(Cell_ID, metadata_df$Cell_ID))
  revised_indel_df <- refined_burden_table_func(input_dir, "revised", "indel", "revised") %>% arrange(match(Cell_ID, metadata_df$Cell_ID))
  original_snv_df <- refined_burden_table_func(input_dir, "original", "snv", "scan2") %>% arrange(match(Cell_ID, metadata_df$Cell_ID))
  original_indel_df <- refined_burden_table_func(input_dir, "original", "indel", "scan2") %>% arrange(match(Cell_ID, metadata_df$Cell_ID))
  
  ##### integrate original and revised snv burden table
  snv_df_01 <- cbind(original_snv_df, metadata_df$Case_ID, "original")
  snv_df_02 <- cbind(revised_snv_df, metadata_df$Case_ID, "revised")
  snv_plot_df <- full_join(snv_df_01 %>% select(Cell_ID, V2, snv.rate.per.gb) %>% rename(original = snv.rate.per.gb), snv_df_02 %>% select(Cell_ID, V2, snv.rate.per.gb) %>% rename(revised = snv.rate.per.gb), by = "Cell_ID") %>% 
    mutate(data_source = case_when(!is.na(original) & !is.na(revised) ~ "modified SCAN2 burden", is.na(original) & !is.na(revised) ~ "recovered missing burden", TRUE ~ NA_character_)) %>% drop_na(data_source)
  
  geom_line_data_snv <- as.data.frame(matrix(c(10, 7500, 10, 7500), nrow = 2, ncol = 2))
  ##### scatter plot compare the original and revised SCAN2 results
  p_original_vs_revised_burden_snv <- ggplot(snv_plot_df, aes(x = original, y = revised)) + 
    geom_line(aes(x = V1, y = V2), colour = "grey50", data = geom_line_data_snv) +
    geom_point(pch = 21, color = "#8B4500", fill = "#EE9A00", size = 1) + 
    scale_x_continuous("original SCAN2 sSNV/GB", transform = 'log10') + scale_y_continuous("revised SCAN2 sSNV/GB", transform = 'log10') + 
    theme_linedraw() + coord_fixed(ratio = 1) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), 
          panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12), plot.margin = margin(t = 5, r = 15, b = 5, l = 5))
  ggsave(paste0(suppl_figure_dir, "/2-original_vs_revised_burden_ssnv.pdf"), plot = p_original_vs_revised_burden_snv, width = 3, height = 3, dpi = 600)
  
  snv_plot_df <- snv_plot_df %>% filter(V2.x %in% unique(V2.x[data_source == "recovered missing burden"]))
  p_recovered_missing_burden_snv <- ggplot(snv_plot_df, aes(x = Cell_ID, y = revised, shape = data_source, color = data_source)) + 
    geom_point(size = 3, na.rm = TRUE, show.legend = TRUE) + labs(x = "Cell ID", y = "sSNV rate per GB") + 
    facet_grid(. ~  V2.x, scales = "free_x", space = "free_x") + scale_y_continuous(trans = "log2", breaks = c(50, 200, 500, 2000)) + 
    scale_color_manual(values = c("recovered missing burden" = "#6b8e23", "modified SCAN2 burden" = "black"), guide = "legend") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), 
          panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(paste0(suppl_figure_dir, "/2-recovered_missing_burden_ssnv.pdf"), plot = p_recovered_missing_burden_snv, width = 8, height = 4, dpi = 600)
  
  ##### integrate original and revised indel burden table
  indel_df_01 <- cbind(original_indel_df, metadata_df$Case_ID, "original")
  indel_df_02 <- cbind(revised_indel_df, metadata_df$Case_ID, "revised")
  indel_plot_df <- full_join(indel_df_01 %>% select(Cell_ID, V2, indel.rate.per.gb) %>% rename(original = indel.rate.per.gb), indel_df_02 %>% select(Cell_ID, V2, indel.rate.per.gb) %>% rename(revised = indel.rate.per.gb), by = "Cell_ID") %>% 
    mutate(data_source = case_when(!is.na(original) & !is.na(revised) ~ "modified SCAN2 burden", is.na(original) & !is.na(revised) ~ "recovered missing burden", TRUE ~ NA_character_)) %>% drop_na(data_source)
  
  geom_line_data_indel <- as.data.frame(matrix(c(0, 820, 0, 820), nrow = 2, ncol = 2))
  ##### scatter plot compare the original and revised SCAN2 results
  p_original_vs_revised_burden_indel <- ggplot(indel_plot_df, aes(x = original, y = revised)) + 
    geom_line(aes(x = V1, y = V2), colour = "grey50", data = geom_line_data_indel) + 
    geom_point(pch = 21, color = "#8B4500", fill = "#EE9A00", size = 1) + 
    scale_x_continuous(limits = c(0, 820), "original SCAN2 sindel/GB") + scale_y_continuous(limits = c(0, 820), "revised SCAN2 sindel/GB") + 
    theme_linedraw() + coord_fixed(ratio = 1) + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), 
          panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12), plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
  ggsave(paste0(suppl_figure_dir, "/2-original_vs_revised_burden_sindel.pdf"), plot = p_original_vs_revised_burden_indel, width = 3, height = 3, dpi = 600)
  
  indel_plot_df <- indel_plot_df %>% filter(V2.x %in% unique(V2.x[data_source == "recovered missing burden"]))
  p_recovered_missing_burden_indel <- ggplot(indel_plot_df, aes(x = Cell_ID, y = revised, shape = data_source, color = data_source)) + 
    geom_point(size = 3, na.rm = TRUE, show.legend = TRUE) + labs(x = "Cell ID", y = "sSNV rate per GB") + facet_grid(. ~  V2.x, scales = "free_x", space = "free_x") + 
    scale_color_manual(values = c("recovered missing burden" = "#6b8e23", "modified SCAN2 burden" = "black"), guide = "legend") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), 
          panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(paste0(suppl_figure_dir, "/2-recovered_missing_burden_sindel.pdf"), plot = p_recovered_missing_burden_indel, width = 8, height = 3, dpi = 600)
  
  return(list(revised_snv_df = revised_snv_df, revised_indel_df = revised_indel_df, original_snv_df = original_snv_df, original_indel_df = original_indel_df))
}