##### Mutation spectrum
factor_hg37 <- 5.845001
type_occur_raw <- mut_type_occurrences(ssnv_grl, ref_genome) %>% filter(row.names(.) %in% Cell_ID_list) %>% arrange(match(row.names(.), Cell_ID_list))

##### Table S1
# type_occur_S1 <- type_occur_raw %>%
#   filter(row.names(.) %in% sSNV_SCAN2_df$Cell_ID) %>%
#   mutate(ture_call_sum = rowSums(select(., 1:6)),
#          Cell_ID = row.names(.))
# sSNV_SCAN2_df_S1 <- sSNV_SCAN2_df[, c("Cell_ID", "Case_ID", "Condition", "Age", "Gender", "snv.rate.per.gb", "snv.burden", "Coverage", "Depth", "MAPD", "CoV")] %>%
#   merge(type_occur_S1) %>% arrange(match(Cell_ID, Cell_ID_list))
# write.csv(sSNV_SCAN2_df_S1, paste0(table_dir, "/3-Table_S1.csv"), row.names = F)

##### calculate estimated burden (per GB, age matched)
type_occur_raw <- type_occur_raw %>% filter(row.names(.) %in% Cell_ID_list_AMG) %>% 
  mutate(across(everything(), ~ . / factor_hg37), Condition = metadata_df_AMG$Condition)
est_ssnv_burden <- sSNV_SCAN2_df_AMG$snv.rate.per.gb
type_occur_est <- type_occur_raw %>% 
  mutate(true_call_sum = rowSums(across(1:6)), snv.rate.per.gb = est_ssnv_burden)

sentivity_ratio <- type_occur_est %>% mutate(sentivity_ratio = true_call_sum / snv.rate.per.gb) %>% 
  dplyr::select(c("Condition", "true_call_sum", "snv.rate.per.gb", "sentivity_ratio"))

est_total_ssnv_ctrl <- sum(est_ssnv_burden[ctrl_range_AMG])
est_total_ssnv_dis <- sum(est_ssnv_burden[dis_range_AMG])
called_total_ssnv_ctrl <- sum(type_occur_est$true_call_sum[ctrl_range_AMG])
called_total_ssnv_dis <- sum(type_occur_est$true_call_sum[dis_range_AMG])
summary_df <- data.frame(est_total_ssnv_ctrl, est_total_ssnv_dis, called_total_ssnv_ctrl, called_total_ssnv_dis)
ctrl_ratio <- est_total_ssnv_ctrl / called_total_ssnv_ctrl
dis_ratio <- est_total_ssnv_dis / called_total_ssnv_dis

type_occur_est[ctrl_range_AMG, 1:8] <- type_occur_est[ctrl_range_AMG, 1:8] * ctrl_ratio
type_occur_est[dis_range_AMG, 1:8] <- type_occur_est[dis_range_AMG, 1:8] * dis_ratio

write.csv(type_occur_est, file = paste0(table_dir, "/3-type_occur_est.csv"))

##### calculate the net change between the conditions and plot
plot_type_list <- c("raw", 'est')
for (plot_type in plot_type_list) {
  if (plot_type == "raw") {
    ## by condition and calculate the net change
    type_occur_plot <- type_occur_raw
    type_avg_per_cell_diff <- colMeans(type_occur_plot[dis_range_AMG, 1:8]) - colMeans(type_occur_plot[ctrl_range_AMG, 1:8])
    # type_occur_plot[nrow(type_occur_plot) + 1, 1:8] <- type_avg_per_cell_diff
    # type_occur_plot[[9]] <- as.character(type_occur_plot[[9]])
    # type_occur_plot[nrow(type_occur_plot), 9] <- "Net change (IHD - Control)"
    
    p_type_occ_raw <- plot_spectrum_customize(type_occur_plot, by = type_occur_plot$Condition, CT = TRUE, indv_points = TRUE, error_bars = "none") +
      facet_wrap(factor(by, unique(type_occur_plot$Condition)) ~ .) + labs(y = "Relative Contribution") + theme_linedraw() +
      theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
            panel.border = element_rect(linewidth = 0.5), text = element_text(size = 15), axis.ticks.x = element_blank())
    ggsave(paste0(other_figure_dir, "/3-sSNV_type_occur_", plot_type, ".pdf"), plot = p_type_occ_raw, width = 10, height = 5.5, dpi = 600)
  }
  else if (plot_type == "est") {
    type_occur_plot <- type_occur_est[, 1:9]
    type_avg_per_cell_diff <- colMeans(type_occur_plot[dis_range_AMG, 1:8]) - colMeans(type_occur_plot[ctrl_range_AMG, 1:8])
    type_occur_plot[nrow(type_occur_plot) + 1, 1:8] <- type_avg_per_cell_diff
    type_occur_plot[[9]] <- as.character(type_occur_plot[[9]])
    type_occur_plot[nrow(type_occur_plot), 9] <- "Net change (IHD - Control)"
    
    p_type_occur_est <- plot_spectrum_absolute(type_occur_plot, by = type_occur_plot$Condition, CT = TRUE, indv_points = TRUE, error_bars = "none") +
      facet_wrap(factor(by, unique(type_occur_plot$Condition)) ~ total_mutations) + labs(y = "sSNV rate per GB") + theme_linedraw() + 
      scale_y_continuous(trans = "sqrt", breaks = c(0, 20, 50, 200, 500, 1000, 2000)) + 
      theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(),
            panel.border = element_rect(linewidth = 0.5), text = element_text(size = 15), axis.ticks.x = element_blank())
    ggsave(paste0(main_figure_dir, "/3-sSNV_type_occur_", plot_type, ".pdf"), plot = p_type_occur_est, width = 10, height = 5.5, dpi = 600)
  } else {
    print('unknown plot type')
  }
  ## test between the two conditions
  melt_type_occur_test <- type_occur_plot %>% dplyr::select(1, 2, 7, 8, 4, 5, 6, 9) %>% 
    slice(1:(n() - 1)) %>% reshape2::melt() %>% 
    mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% 
    setNames(c("Condition", "type", "no_mutation"))
  
  barplot_type_occur <- ggbarplot(melt_type_occur_test, x = "type", y = "no_mutation", color = "Condition", label = TRUE, lab.size = 6, 
                                  lab.nb.digits = 0, add = c("mean_se", "jitter"), position = position_dodge(0.9), palette = ctrl_dis_color) + 
    stat_compare_means(aes(group = Condition, label = sprintf("p = %1.2e", as.numeric(..p.format..))), label.y = 1.04 * max(melt_type_occur_test$no_mutation), size = 6) +
    stat_compare_means(aes(group = Condition), label = "p.signif", label.y = 1.10 * max(melt_type_occur_test$no_mutation), size = 6) +
    labs(x = "Point mutation type", y = "sSNV rate per GB") + theme_linedraw() + facet_grid(. ~ type, scales = "free_x", space = "free_x") + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 15), axis.ticks.x = element_blank(), axis.text.x = element_blank())
  ggsave(paste0(other_figure_dir, "/3-sSNV_type_occur_test_", plot_type, ".pdf"), plot = barplot_type_occur, width = 15, height = 5.5, dpi = 600)
}
