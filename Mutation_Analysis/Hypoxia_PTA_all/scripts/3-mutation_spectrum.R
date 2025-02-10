##### Mutation spectrum
factor_hg37 <- 5.845001
type_occurrences <- mut_type_occurrences(snv_grl, ref_genome) %>% filter(row.names(.) %in% Cell_ID_list) %>% arrange(match(row.names(.), Cell_ID_list))

##### Table S1
# type_occurrences_S1 <- type_occurrences %>% filter(row.names(.) %in% SCAN2_df$Cell_ID) %>% mutate(ture_call_sum = rowSums(select(., 1:6)), Cell_ID = row.names(.))
# SCAN2_df_S1 <- SCAN2_df[, c("Cell_ID", "Case_ID", "Condition", "Age", "Gender", "snv.rate.per.gb", "snv.burden", "Coverage", "Depth", "MAPD", "CoV")] %>% merge(type_occurrences_S1) %>% arrange(match(Cell_ID, Cell_ID_list))
# write.csv(SCAN2_df_S1, paste0(table_dir, "/Table_S1.csv"), row.names = F)

##### calculate estimated burden (per GB, age matched)
type_occurrences <- type_occurrences %>% filter(row.names(.) %in% Cell_ID_list_AMG) %>% mutate(across(everything(), ~ . / factor_hg37), Condition = metadata_df_AMG$Condition)
est_ssnv_burden <- SCAN2_df_AMG$snv.rate.per.gb
type_occurrences_est <- type_occurrences %>% mutate(true_call_sum = rowSums(across(1:6)), snv.rate.per.gb = est_ssnv_burden)
sentivity_ratio <- type_occurrences_est %>% mutate(sentivity_ratio = true_call_sum / snv.rate.per.gb) |> base::`[`(c("Condition", "true_call_sum", "snv.rate.per.gb", "sentivity_ratio"))

est_total_ssnv_ctrl <- sum(est_ssnv_burden[ctrl_range_AMG])
est_total_ssnv_dis <- sum(est_ssnv_burden[dis_range_AMG])
called_total_ssnv_ctrl <- sum(type_occurrences_est$true_call_sum[ctrl_range_AMG])
called_total_ssnv_dis <- sum(type_occurrences_est$true_call_sum[dis_range_AMG])
summary_df <- data.frame(est_total_ssnv_ctrl, est_total_ssnv_dis, called_total_ssnv_ctrl, called_total_ssnv_dis)
ctrl_ratio <- est_total_ssnv_ctrl / called_total_ssnv_ctrl
dis_ratio <- est_total_ssnv_dis / called_total_ssnv_dis

type_occurrences_est[ctrl_range_AMG, 1:8] <- type_occurrences_est[ctrl_range_AMG, 1:8] * ctrl_ratio
type_occurrences_est[dis_range_AMG, 1:8] <- type_occurrences_est[dis_range_AMG, 1:8] * dis_ratio

##### calculate diff matrix and plot
# plot_type_list <- c('true_call', 'estimated')
plot_type_list <- c('estimated')
for (plot_type in plot_type_list) {
  if (plot_type == 'true_call') {
    type_occurrences <- type_occurrences
  }
  else if (plot_type == 'estimated') {
    type_occurrences <- type_occurrences_est[, 1:9]
  } else {
    print('unknown plot type')
  }

  type_avg_per_cell_diff <- colMeans(type_occurrences[dis_range_AMG, 1:8]) - colMeans(type_occurrences[ctrl_range_AMG, 1:8])
  type_occurrences[nrow(type_occurrences_est) + 1, 1:8] <- type_avg_per_cell_diff
  type_occurrences[nrow(type_occurrences_est) + 1, 9] <- "Net change (IHD - Control)"

  p_type_occ_absolute <- plot_spectrum_absolute(type_occurrences, by = type_occurrences$Condition, CT = TRUE, indv_points = TRUE, error_bars = "none") + 
    facet_wrap(factor(by, unique(type_occurrences$Condition)) ~ total_mutations) + labs(y = "sSNV rate per GB") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 15), axis.ticks.x = element_blank())
  ggsave(paste0(main_figure_dir, "/3-sSNV_type_occurrences_by_condition.pdf"), plot = p_type_occ_absolute, width = 10, height = 5.5, dpi = 600)
  
  melt_type_occ_t_test <- type_occurrences %>% dplyr::select(1, 2, 7, 8, 4, 5, 6, 9) %>% slice(1:(n() - 1)) %>% 
    reshape2::melt() %>% mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% setNames(c("Condition", "type", "no_mutation"))
  
  barplot_type_occ_t_test <- ggbarplot(melt_type_occ_t_test, x = "type", y = "no_mutation", color = "Condition", label = TRUE, lab.size = 6, 
                                       lab.nb.digits = 0, add = c("mean_se", "jitter"), position = position_dodge(0.9), palette = ctrl_dis_color) + 
    stat_compare_means(aes(group = Condition, label = sprintf("p = %1.2e", as.numeric(..p.format..))), label.y = 1.04 * max(melt_type_occ_t_test$no_mutation), size = 6) +
    stat_compare_means(aes(group = Condition), label = "p.signif", label.y = 1.10 * max(melt_type_occ_t_test$no_mutation), size = 6) +
    labs(x = "Point mutation type", y = "sSNV rate per GB") + theme_linedraw() + facet_grid(. ~ type, scales = "free_x", space = "free_x") + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 15), axis.ticks.x = element_blank(), axis.text.x = element_blank())
  ggsave(paste0(main_figure_dir, "/3-sSNV_type_occurrences_t_test.pdf"), plot = barplot_type_occ_t_test, width = 15, height = 5.5, dpi = 600)
}