##### Mutation spectrum
factor_hg37 <- 5.845001
type_occurrences <- mut_type_occurrences(grl, ref_genome) %>% 
  arrange(match(row.names(.), Cell_ID_list))
type_occurrences <- type_occurrences[Cell_ID_list_age_match, ]
type_occurrences <- type_occurrences / factor_hg37
type_occurrences$condition <- metadata_df_age_match$Condition

##### calculate estimated burden
estimated_ssnv_burden <- SCAN2_df_age_match$snv.rate.per.gb
# estimated_ssnv_burden <- SCAN2_df_age_match$snv.burden
type_occurrences_estimated <- type_occurrences
type_occurrences_estimated$true_call_sum <- rowSums(type_occurrences_estimated[1:6])
type_occurrences_estimated$snv.rate.per.gb <- estimated_ssnv_burden
sentivity_ratio <- type_occurrences_estimated %>% 
  mutate(sentivity_ratio = true_call_sum / snv.rate.per.gb) |> base::`[`(c("condition", "true_call_sum", "snv.rate.per.gb", "sentivity_ratio"))

est_total_ssnv_control <- sum(estimated_ssnv_burden[control_range_age_match])
est_total_ssnv_disease <- sum(estimated_ssnv_burden[disease_range_age_match])
called_total_ssnv_control <- sum(type_occurrences_estimated$true_call_sum[control_range_age_match])
called_total_ssnv_disease <- sum(type_occurrences_estimated$true_call_sum[disease_range_age_match])
summary_df <- data.frame(est_total_ssnv_control,est_total_ssnv_disease,called_total_ssnv_control,called_total_ssnv_disease)
# summary_df/factor_hg37
control_ratio <- est_total_ssnv_control / called_total_ssnv_control
disease_ratio <- est_total_ssnv_disease / called_total_ssnv_disease

type_occurrences_estimated[control_range_age_match,1:8] <- type_occurrences_estimated[control_range_age_match,1:8] * control_ratio
type_occurrences_estimated[disease_range_age_match,1:8] <- type_occurrences_estimated[disease_range_age_match,1:8] * disease_ratio

##### calculate diff matrix and plot
# plot_type_list <- c('true_call', 'estimated')
plot_type_list <- c('estimated')
for (plot_type in plot_type_list) {
  if (plot_type == 'true_call') {
    type_occurrences <- type_occurrences
    print("true")
  }
  else if (plot_type == 'estimated') {
    type_occurrences <- type_occurrences_estimated[,1:9]
    print("est")
  } else {
    print('unknown plot type')
  }

  type_avg_per_cell_diff <- colMeans(type_occurrences[disease_range_age_match, 1:8]) - colMeans(type_occurrences[control_range_age_match, 1:8])
  type_occurrences[nrow(type_occurrences_estimated) + 1, 1:8] <- type_avg_per_cell_diff
  type_occurrences[nrow(type_occurrences_estimated) + 1, 9] <- "Net change (IHD - Control)"

  yticks_03 <- c(0, 200, 400, 600, 990, 1000, 1790, 1800) 
  p_type_occurrences_absolute <- plot_spectrum_absolute(type_occurrences, by = type_occurrences$condition, CT = TRUE, indv_points = TRUE) + 
    facet_wrap(factor(by, unique(type_occurrences$condition)) ~ total_mutations) + scale_y_break(c(600, 990, 1001, 1790), scale=0.2) + 
    scale_y_continuous(limits = c(0, 1800), breaks = yticks_03, labels = yticks_03) + labs(y = "sSNV rate per GB") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
          panel.border = element_rect(size = 0.5), text = element_text(size = 15), axis.ticks.x = element_blank())
  ggsave(paste0(main_figure_dir, "/3-sSNV_type_occurrences_by_condition.pdf"), plot = p_type_occurrences_absolute, width = 10, height = 5.5, dpi = 600)
  # type_occurrences_t_test_temp <- as.data.frame(type_occurrences[,c(1,2,7,8,4,5,6,9)])
  # melt_type_occurrences_t_test_temp <- reshape2::melt(type_occurrences_t_test_temp)
  # colnames(melt_type_occurrences_t_test_temp) <- c("condition", "type", "no_mutation")
  # melt_type_occurrences_t_test_temp$conditions <- factor(melt_type_occurrences_t_test_temp$condition, level = c(control_name, disease_name))
  # barplot_melt_type_occurrences_t_test_temp <- ggbarplot(melt_type_occurrences_t_test_temp, x = "type", y = "no_mutation", color = "conditions",
  #                                              label = TRUE, lab.size = 6, lab.nb.digits = 0, add = c("mean_se", "jitter"), position = position_dodge(0.9), 
  #                                              palette = dis_ctrl_color) +
  #   stat_compare_means(aes(group = condition, label = sprintf("p = %1.2f", as.numeric(..p.format..))), label.y = 1.04 * max(melt_type_occurrences_t_test_temp$no_mutation), size = 6) + 
  #   stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.08 * max(melt_type_occurrences_t_test_temp$no_mutation), size = 6) +
  #   labs(x = "Point mutation type", y = "sSNV rate per GB") + 
  #   scale_y_break(c(600, 990, 1001, 1790), scale=0.2) +
  #   scale_y_continuous(limits=c(0, 1800), breaks=yticks_03, labels=yticks_03) +
  #   theme_classic(base_size = 24) +
  #   theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
  #         axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white"), 
  #         axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank())
  # # facet_wrap(factor(mut_type, c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) ~ .)
  # print(barplot_melt_type_occurrences_t_test_temp)
}