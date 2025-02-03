##### indel spectra 
indel_counts_age <- indel_counts[,Cell_ID_list_age_match]
indel_counts_conditional <- indel_counts_age
Control <- rowSums(indel_counts_age[,control_range_age_match]) / length(control_range_age_match)
indel_counts_conditional <- cbind(indel_counts_conditional, Control)
Disease <- rowSums(indel_counts_age[,disease_range_age_match]) / length(disease_range_age_match)
indel_counts_conditional <- cbind(indel_counts_conditional, Disease)
Net_change <- Disease - Control
# Net_change[Net_change<0] <- NA
indel_counts_conditional <- cbind(indel_counts_conditional, Net_change)
colnames(indel_counts_conditional)[(ncol(indel_counts_conditional)-2):ncol(indel_counts_conditional)] <- c(control_name, disease_name, 'Net Change')
######### plot indel spectra
indel_counts_deletion <- indel_counts_age[grepl('deletion', row.names(indel_counts_age)),]
indel_counts_insertion <- indel_counts_age[grepl('insertion', row.names(indel_counts_age)),]
total_deletion_normoxia <- sum(indel_counts_deletion[,control_range_age_match])
total_deletion_hypoxia <- sum(indel_counts_deletion[,disease_range_age_match])
total_insertion_normoxia <- sum(indel_counts_insertion[,control_range_age_match])
total_insertion_hypoxia <- sum(indel_counts_insertion[,disease_range_age_match])
C_deletion_normoxia <- sum(indel_counts_conditional[1:6,23])
C_deletion_hypoxia <- sum(indel_counts_conditional[1:6,24])
T_deletion_normoxia <- sum(indel_counts_conditional[7:12,23])
T_deletion_hypoxia <- sum(indel_counts_conditional[7:12,24])

C_deletion_percent_normoxia <- C_deletion_normoxia / (total_deletion_normoxia + total_insertion_normoxia)
T_deletion_percent_normoxia <- T_deletion_normoxia / (total_deletion_normoxia + total_insertion_normoxia)
C_deletion_percent_hypoxia <- C_deletion_hypoxia / (total_deletion_hypoxia + total_insertion_hypoxia)
T_deletion_percent_hypoxia <- T_deletion_hypoxia / (total_deletion_hypoxia + total_insertion_hypoxia)
C_deletion_percent_total <- (C_deletion_normoxia + C_deletion_hypoxia) / (total_deletion_normoxia + total_insertion_normoxia)
T_deletion_percent_total <- (T_deletion_normoxia + T_deletion_hypoxia) / (total_deletion_hypoxia + total_insertion_hypoxia)

deletion_percent_normoxia <- total_deletion_normoxia / (total_deletion_normoxia + total_insertion_normoxia)
deletion_percent_hypoxia <- total_deletion_hypoxia / (total_deletion_hypoxia + total_insertion_hypoxia)
deletion_total <- (total_deletion_normoxia + total_deletion_hypoxia) / (total_deletion_normoxia + total_insertion_normoxia + total_deletion_hypoxia + total_insertion_hypoxia)

pdf(paste0(sindel_plot_dir, "/8-indel_mutation_pattern.pdf"), width = 12, height = 8)
print(plot_indel_contexts(indel_counts_conditional[, (ncol(indel_counts_conditional)-2):ncol(indel_counts_conditional)], condensed = TRUE))
      print(plot_main_indel_contexts(indel_counts_conditional[, (ncol(indel_counts_conditional)-2):ncol(indel_counts_conditional)]))
dev.off()

########## calculate estimated burden
estimated_sindel_burden <- SCAN2_df_age_match$indel.burden

est_total_sindel_control <- sum(estimated_sindel_burden[control_range_age_match])
est_total_sindel_disease <- sum(estimated_sindel_burden[disease_range_age_match])
called_total_sindel_control <- sum(colSums(indel_counts_age)[control_range_age_match])
called_total_sindel_disease <- sum(colSums(indel_counts_age)[disease_range_age_match])
summary_df <- data.frame(est_total_sindel_control,est_total_sindel_disease,called_total_sindel_control,called_total_sindel_disease)
control_ratio <- est_total_sindel_control / called_total_sindel_control
disease_ratio <- est_total_sindel_disease / called_total_sindel_disease

indel_counts_est <- indel_counts_age
indel_counts_est[,control_range_age_match] <- indel_counts_est[,control_range_age_match] * control_ratio
indel_counts_est[,disease_range_age_match] <- indel_counts_est[,disease_range_age_match] * disease_ratio

Control <- rowSums(indel_counts_est[,control_range_age_match]) / length(control_range_age_match)
indel_counts_est <- cbind(indel_counts_est, Control)
Disease <- rowSums(indel_counts_est[,disease_range_age_match]) / length(disease_range_age_match)
indel_counts_est <- cbind(indel_counts_est, Disease)
Net_change <- Disease - Control
# Net_change[Net_change<0] <- NA
indel_counts_est <- cbind(indel_counts_est, Net_change)
colnames(indel_counts_est)[(ncol(indel_counts_est)-2):ncol(indel_counts_est)] <- c(control_name, disease_name, 'Net Change')

########### calculate diff matrix for estimated burden
indel_counts_diff <- indel_counts_est - indel_counts_conditional

pdf(paste0(sindel_plot_dir, "/8-indel_mutation_pattern_estmiated.pdf"), width = 12, height = 8)
print(plot_indel_contexts(round(indel_counts_est[, (ncol(indel_counts_est)-2):ncol(indel_counts_est)], digits = 0), condensed = TRUE))
print(plot_main_indel_contexts(round(indel_counts_est[, (ncol(indel_counts_est)-2):ncol(indel_counts_est)], digits = 0)))
dev.off()

##### relative and absolute mutation for each type between two conditions
type_list <- c("C_deletion", "T_deletion", "C_insertion", "T_insertion", "2bp_deletion", "3bp_deletion", "4bp_deletion", 
               "5+bp_deletion", "2bp_insertion", "3bp_insertion", "4bp_insertion", "5+bp_insertion", "2bp_deletion_with_microhomology",
               "3bp_deletion_with_microhomology", "4bp_deletion_with_microhomology", "5+bp_deletion_with_microhomology")
# round_est <- round(indel_counts_est, digits = 2)
round_est <- indel_counts_est
total_normoxia <- sum(round_est[,23])
total_hypoxia <- sum(round_est[,24])

C_deletion <- colSums(round_est[1:6,])
T_deletion <- colSums(round_est[7:12,])
C_insertion <- colSums(round_est[13:18,])
T_insertion <- colSums(round_est[19:24,])
two_bp_deletion <- colSums(round_est[25:30,])
three_bp_deletion <- colSums(round_est[31:36,])
four_bp_deletion <- colSums(round_est[37:42,])
five_plus_bp_deletion <- colSums(round_est[43:48,])
two_bp_insertion <- colSums(round_est[49:54,])
three_bp_insertion <- colSums(round_est[55:60,])
four_bp_insertion <- colSums(round_est[61:66,])
five_plus_bp_insertion <- colSums(round_est[67:72,])
two_bp_deletion_with_microhomology <- round_est[73,]
three_bp_deletion_with_microhomology <- colSums(round_est[74:75,])
four_bp_deletion_with_microhomology <- colSums(round_est[76:78,])
five_plus_bp_deletion_with_microhomology <- colSums(round_est[79:83,])

condensed_indel_type <- rbind(C_deletion, T_deletion, C_insertion, T_insertion, two_bp_deletion, three_bp_deletion, 
      four_bp_deletion, five_plus_bp_deletion, two_bp_insertion, three_bp_insertion, 
      four_bp_insertion, five_plus_bp_insertion, two_bp_deletion_with_microhomology, 
      three_bp_deletion_with_microhomology, four_bp_deletion_with_microhomology, 
      five_plus_bp_deletion_with_microhomology)

rownames(condensed_indel_type) <- type_list

absolute_normoxia_each <- round(condensed_indel_type[,control_range_age_match] / length(control_range_age_match), digits = 2)
absolute_hypoxia_each <- round(condensed_indel_type[,disease_range_age_match] / length(disease_range_age_match), digits = 2)
proportion_normoxia_each <- round(condensed_indel_type[,control_range_age_match] / (length(control_range_age_match) * total_normoxia), digits = 4)
proportion_hypoxia_each <- round(condensed_indel_type[,disease_range_age_match] / (length(disease_range_age_match) * total_hypoxia), digits = 4)

# absolute_normoxia_each <- round(condensed_indel_type[,control_range_age_match] / length(control_range_age_match), digits = 2)
# absolute_hypoxia_each <- round(condensed_indel_type[,disease_range_age_match] / length(disease_range_age_match), digits = 2)
# proportion_normoxia_each <- round(condensed_indel_type[,control_range_age_match] / (length(control_range_age_match) * total_normoxia), digits = 4)
# proportion_hypoxia_each <- round(condensed_indel_type[disease_range_age_match] / (length(disease_range_age_match) * total_hypoxia), digits = 4)

### statistical test for absolute
melt_condensed_indel_type_absolute <- as.data.frame(t(condensed_indel_type[,c(disease_range_age_match, control_range_age_match)]))
melt_condensed_indel_type_absolute$condition <- metadata_df_age_match$Condition
melt_condensed_indel_type_absolute <- melt(melt_condensed_indel_type_absolute)
colnames(melt_condensed_indel_type_absolute) <- c("condition", "type", "no_mutation")
melt_condensed_indel_type_absolute$conditions <- factor(melt_condensed_indel_type_absolute$condition, level = c(control_name, disease_name))
# melt_matrix_no_mutation$mut_type <- c(rep("C>A",44), rep("C>G",44), rep("C>T",44), rep("T>A",44), rep("T>C",44), rep("T>G",44))
barplot_melt_condensed_indel_type_absolute <- ggbarplot(melt_condensed_indel_type_absolute, x = "type", y = "no_mutation", color = "conditions",
                                                       label = TRUE, lab.size = 6, lab.nb.digits = 0, add = c("mean_se", "jitter"), position = position_dodge(0.9), 
                                                       palette = c("dodgerblue2", "firebrick2")) +
  stat_compare_means(aes(group = condition, label = sprintf("p = %1.3f", as.numeric(..p.format..))), label.y = 1.02 * max(melt_condensed_indel_type_absolute$no_mutation), size = 6) + 
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.07 * max(melt_condensed_indel_type_absolute$no_mutation), size = 6) +
  labs(x = "Point mutation type", y = "number of sSNV") + 
  # scale_y_break(c(3500, 5800, 6000, 9450), scale=0.2) +
  # scale_y_continuous(limits=c(0,12000), breaks=yticks_03, labels=yticks_03) +
  theme_classic(base_size = 24) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white"), 
        axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank())
# facet_wrap(factor(mut_type, c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) ~ .)
print(barplot_melt_condensed_indel_type_absolute)
ggsave(paste0(sindel_plot_dir, "/high_quality/barplot_melt_type_occurrences_t_test_temp.pdf"), plot = barplot_melt_condensed_indel_type_absolute, width = 23, height = 10)

### statistical test for relative
melt_condensed_indel_type_relative <- as.data.frame(t(cbind(proportion_hypoxia_each, proportion_normoxia_each)))
melt_condensed_indel_type_relative$condition <- metadata_df_age_match$Condition
melt_condensed_indel_type_relative <- melt(melt_condensed_indel_type_relative)
colnames(melt_condensed_indel_type_relative) <- c("condition", "type", "no_mutation")
melt_condensed_indel_type_relative$conditions <- factor(melt_condensed_indel_type_relative$condition, level = c(control_name, disease_name))
# melt_matrix_no_mutation$mut_type <- c(rep("C>A",44), rep("C>G",44), rep("C>T",44), rep("T>A",44), rep("T>C",44), rep("T>G",44))
barplot_melt_condensed_indel_type_relative <- ggbarplot(melt_condensed_indel_type_relative, x = "type", y = "no_mutation", color = "conditions",
                                                        label = TRUE, lab.size = 6, lab.nb.digits = 0, add = c("mean_se", "jitter"), position = position_dodge(0.9), 
                                                        palette = c("dodgerblue2", "firebrick2")) +
  stat_compare_means(aes(group = condition, label = sprintf("p = %1.3f", as.numeric(..p.format..))), label.y = 1.02 * max(melt_condensed_indel_type_relative$no_mutation), size = 6) + 
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.07 * max(melt_condensed_indel_type_relative$no_mutation), size = 6) +
  labs(x = "Point mutation type", y = "number of sSNV") + 
  # scale_y_break(c(3500, 5800, 6000, 9450), scale=0.2) +
  # scale_y_continuous(limits=c(0,12000), breaks=yticks_03, labels=yticks_03) +
  theme_classic(base_size = 24) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white"), 
        axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank())
# facet_wrap(factor(mut_type, c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) ~ .)
print(barplot_melt_condensed_indel_type_relative)
ggsave(paste0(sindel_plot_dir, "/high_quality/barplot_melt_condensed_indel_type_relative.pdf"), plot = barplot_melt_condensed_indel_type_relative, width = 23, height = 10)

# pdf(paste0(sindel_plot_dir, "/indel_type_propotion_each_condition.pdf"), width = 15, height = 8)
# plot_contribution_heatmap_cus(all_proportion, cluster_sigs = FALSE, cluster_samples = FALSE) + 
#   scale_fill_gradient(low = "white", high = "blue") + geom_tile(colour = "black") + 
#   geom_text(aes(label = round(Reshape(all_proportion,nrow(all_proportion)*2,1), 2)), size = 6) +
#   theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.5)) + 
#   theme(legend.title=element_blank())
# dev.off()

C_deletion <- colSums(indel_counts[1:6,])
T_deletion <- colSums(indel_counts[7:12,])
C_insertion <- colSums(indel_counts[13:18,])
T_insertion <- colSums(indel_counts[19:24,])
two_bp_deletion <- colSums(indel_counts[25:30,])
three_bp_deletion <- colSums(indel_counts[31:36,])
four_bp_deletion <- colSums(indel_counts[37:42,])
five_plus_bp_deletion <- colSums(indel_counts[43:48,])
two_bp_insertion <- colSums(indel_counts[49:54,])
three_bp_insertion <- colSums(indel_counts[55:60,])
four_bp_insertion <- colSums(indel_counts[61:66,])
five_plus_bp_insertion <- colSums(indel_counts[67:72,])
two_bp_deletion_with_microhomology <- indel_counts[73,]
three_bp_deletion_with_microhomology <- colSums(indel_counts[74:75,])
four_bp_deletion_with_microhomology <- colSums(indel_counts[76:78,])
five_plus_bp_deletion_with_microhomology <- colSums(indel_counts[79:83,])
condensed_indel_type <- rbind(C_deletion, T_deletion, C_insertion, T_insertion, two_bp_deletion, three_bp_deletion, 
                              four_bp_deletion, five_plus_bp_deletion, two_bp_insertion, three_bp_insertion, 
                              four_bp_insertion, five_plus_bp_insertion, two_bp_deletion_with_microhomology, 
                              three_bp_deletion_with_microhomology, four_bp_deletion_with_microhomology, 
                              five_plus_bp_deletion_with_microhomology)

rownames(condensed_indel_type) <- type_list

C_T_deletion_together <- colSums(indel_counts[1:12,])
C_T_insertion_together <- colSums(indel_counts[13:24,])
two_bp_deletion_with_microhomology_ <- indel_counts[73,]
two_bp_deletion_together <- colSums(indel_counts[c(25:30, 73),])
one_bp_delete_ratio <- C_T_deletion_together / (C_T_deletion_together + two_bp_deletion_together)
two_bp_delete_ratio <- two_bp_deletion_together / (C_T_deletion_together + two_bp_deletion_together)
# one_bp_delete_ratio_age_match <- one_bp_delete_ratio[Cell_ID_list_age_match]
# two_bp_delete_ratio_age_match <- two_bp_delete_ratio[Cell_ID_list_age_match]
# t.test(one_bp_delete_ratio, two_bp_delete_ratio)
# t.test(one_bp_delete_ratio_age_match, two_bp_delete_ratio_age_match)
# t.test(one_bp_delete_ratio[control_range_age_match], two_bp_delete_ratio[control_range_age_match])
# t.test(one_bp_delete_ratio[disease_range_age_match], two_bp_delete_ratio[disease_range_age_match])
# t.test(one_bp_delete_ratio_age_match[control_range_age_match], one_bp_delete_ratio_age_match[disease_range_age_match])
# t.test(two_bp_delete_ratio_age_match[control_range_age_match], two_bp_delete_ratio_age_match[disease_range_age_match])
# t.test(C_T_deletion_together[control_range_age_match], C_T_deletion_together[disease_range_age_match])
# t.test(two_bp_deletion_together[control_range_age_match], two_bp_deletion_together[disease_range_age_match])
# t.test(two_bp_deletion_with_microhomology[1:12], two_bp_deletion_with_microhomology[13:22])
selected_indel_type <- as.data.frame(t(rbind(condensed_indel_type, C_T_deletion_together, C_T_insertion_together, 
                                             two_bp_deletion_with_microhomology_, one_bp_delete_ratio, two_bp_delete_ratio, two_bp_deletion_together)))
selected_indel_type$condition <- metadata_df$Condition
selected_indel_type$Case_ID <- metadata_df$Case_ID
selected_indel_type$age <- metadata_df$Age
selected_indel_type_relevel <- selected_indel_type %>% mutate(condition = relevel(factor(condition), ref = "Normal"))
selected_indel_type_Disease <- selected_indel_type_relevel[selected_indel_type_relevel$condition == "Disease",]
selected_indel_type_Control <- selected_indel_type_relevel[selected_indel_type_relevel$condition == "Normal",]
selected_indel_type_relevel$two_bp_deletion_with_microhomology
##### mix linear model fitting
indel_type_list <- c(type_list, "C_T_deletion_together", "C_T_insertion_together", 
                     "two_bp_deletion_with_microhomology_", "one_bp_delete_ratio", "two_bp_delete_ratio", "two_bp_deletion_together")
# plot_list <- c(1:21)
plot_list <- c(1:22)
all_SBS_burden <- metadata_df$Cell_ID
all_SBS_burden_age_match <- metadata_df_age_match$Cell_ID
all_mix_effect_model_p_value <- metadata_df[1,1]

pdf(paste0(sindel_plot_dir, "/8-top3_point_type_with_age_condition.pdf"), width = 12, height = 8)
for (i in plot_list){
  ##### all cell all age
  absolute_sig_contri <- selected_indel_type_relevel[c(indel_type_list[i],'age','condition','Case_ID')]
  colnames(absolute_sig_contri)[1] <- 'SBS'
  all_SBS_burden <- cbind(all_SBS_burden, absolute_sig_contri[1])
  absolute_sig_contri_Disease <- absolute_sig_contri[absolute_sig_contri$condition == "Disease",]
  absolute_sig_contri_Control <- absolute_sig_contri[absolute_sig_contri$condition == "Normal",]
  
  burden_age_model_02 <- lmer(SBS ~ age + condition + (1|Case_ID), absolute_sig_contri)
  burden_age_model_Disease_02 <- lmer(SBS ~ age + (1|Case_ID), absolute_sig_contri_Disease)
  burden_age_model_Control_02 <- lmer(SBS ~ age + (1|Case_ID), absolute_sig_contri_Control)
  summary_lm_all_02 <- summary(burden_age_model_02)
  summary_lm_Disease_02 <- summary(burden_age_model_Disease_02)
  summary_lm_Control_02 <- summary(burden_age_model_Control_02)
  
  anova_pvalue_02 <- anova(burden_age_model_02)$"Pr(>F)"[2]
  r.squaredGLMM(burden_age_model_Disease_02)[1]
  r.squaredGLMM(burden_age_model_Control_02)[1]
  all_mix_effect_model_p_value <- cbind(all_mix_effect_model_p_value, anova_pvalue_02)
  
  ## manually calculate the fitting lines
  Disease_fitting_x = range(absolute_sig_contri_Disease$age)
  Disease_fitting_y = Disease_fitting_x * fixef(burden_age_model_Disease_02)[2] + fixef(burden_age_model_Disease_02)[1]
  Control_fitting_x = range(absolute_sig_contri_Control$age)
  Control_fitting_y = Control_fitting_x * fixef(burden_age_model_Control_02)[2] + fixef(burden_age_model_Control_02)[1]
  
  # if (!all(absolute_sig_contri_Control$SBS == 0)){
  #   print(all(is.na(absolute_sig_contri_Control$SBS)))
  
  absolute_sig_contri <- absolute_sig_contri %>% mutate(condition = relevel(factor(condition), ref = "Normal"))
  absolute_sig_contri_plot <- absolute_sig_contri %>% group_by(factor(condition, levels = c(disease_name, control_name))) %>% arrange(age, .by_group = TRUE)
  absolute_sig_contri_plot$Case_IDs <- factor(absolute_sig_contri_plot$Case_ID, levels = unique(absolute_sig_contri_plot$Case_ID))
  legend_data_04 <- absolute_sig_contri_plot[c(7,30), c("age", "SBS", "condition")]
  
  p_SNV_burden_per_sig <- ggplot(absolute_sig_contri, aes(x = age, y = SBS, color = condition, fill = condition)) +
    geom_point(pch = 21, data = legend_data_04, aes(x = age, y = SBS, color = condition, fill = condition), size = 5) +
    geom_point(pch = 21, fill = c(disease_color, control_color), size = 5, show.legend = FALSE) +
    geom_segment(aes(x = Disease_fitting_x[1], xend = Disease_fitting_x[2], 
                     y = Disease_fitting_y[1], yend = Disease_fitting_y[2]), colour = "#FC4E07") + 
    geom_segment(aes(x = Control_fitting_x[1], xend = Control_fitting_x[2], 
                     y = Control_fitting_y[1], yend = Control_fitting_y[2]), colour = "dodgerblue3") + 
    annotate("text", size = 7, x = 0.0 * max(absolute_sig_contri$age), y = 1.0 * max(absolute_sig_contri$SBS),
             # label = paste("P = ", format(anova_pvalue_02, scientific = TRUE, digits = 3)), hjust = 0) +
             label = paste("P = ", sprintf("%.3f",anova_pvalue_02)), hjust = 0) +
    scale_color_manual(values = c("Disease" = "black", "Normal" = "black"), guide = "none") +
    scale_fill_manual(values = c("Normal" = dis_ctrl_color[1], "Disease" = dis_ctrl_color[2]), guide = "legend") +
    labs(x = "Age", y = paste0(indel_type_list[i], " \n (sindel rate per cell)")) + theme_classic() +
    theme(text = element_text(size=24), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position="right")
  print(p_SNV_burden_per_sig)
  
  ##### age matched cells
  # absolute_sig_contri_age_match <- (top10_contri_mat_age_match[c(refined_selected_sigs_list[i],'age','condition','Case_ID')])
  absolute_sig_contri_age_match <- absolute_sig_contri[rownames(absolute_sig_contri) %in% Cell_ID_list_age_match, ]
  colnames(absolute_sig_contri_age_match)[1] <- 'SBS'
  absolute_sig_contri_Disease_age_match <- absolute_sig_contri_age_match[absolute_sig_contri_age_match$condition == "Disease",]
  absolute_sig_contri_Control_age_match <- absolute_sig_contri_age_match[absolute_sig_contri_age_match$condition == "Normal",]
  absolute_sig_contri_Disease_age_match$label <- "Disease"
  absolute_sig_contri_Control_age_match$label <- "Age-matched control"
  
  burden_age_model_03 <- lmer(SBS ~ age + condition + (1|Case_ID), absolute_sig_contri_age_match)
  burden_age_model_Disease_03 <- lmer(SBS ~ age + (1|Case_ID), absolute_sig_contri_Disease_age_match)
  burden_age_model_Control_03 <- lmer(SBS ~ age + (1|Case_ID), absolute_sig_contri_Control_age_match)
  summary_lm_all_03 <- summary(burden_age_model_03)
  summary_lm_Disease_03 <- summary(burden_age_model_Disease_03)
  summary_lm_Control_03 <- summary(burden_age_model_Control_03)
  
  anova_pvalue_03 <- anova(burden_age_model_03)$"Pr(>F)"[2]
  r.squaredGLMM(burden_age_model_Disease_03)[1]
  r.squaredGLMM(burden_age_model_Control_03)[1]
  
  absolute_sig_contri_Control_age_match$SBS_remove_age <- absolute_sig_contri_Control_age_match$SBS - 
    fixef(burden_age_model_Control_02)[2] * absolute_sig_contri_Control_age_match$age
  burden_age_model_Control_remove_age_03 <- lmer(SBS_remove_age ~ age + (1|Case_ID), absolute_sig_contri_Control_age_match)
  summary_lm_Control_remove_age_03 <- summary(burden_age_model_Control_remove_age_03)
  
  absolute_sig_contri_Disease_age_match$SBS_remove_age <- absolute_sig_contri_Disease_age_match$SBS - 
    fixef(burden_age_model_Control_02)[2] * absolute_sig_contri_Disease_age_match$age
  burden_age_model_Disease_remove_age_03 <- lmer(SBS_remove_age ~ age + (1|Case_ID), absolute_sig_contri_Disease_age_match)
  summary_lm_Disease_remove_age_03 <- summary(burden_age_model_Disease_remove_age_03)
  
  absolute_sig_contri_age_match$SBS_remove_age <-
    c(absolute_sig_contri_Disease_age_match$SBS_remove_age, absolute_sig_contri_Control_age_match$SBS_remove_age)
  absolute_sig_contri_age_match$label <- c(absolute_sig_contri_Disease_age_match$label, absolute_sig_contri_Control_age_match$label)
  mean_absolute_sig_contri_disease_remove <- mean(absolute_sig_contri_Disease_age_match$SBS_remove_age)
  mean_absolute_sig_contri_control_remove <- mean(absolute_sig_contri_Control_age_match$SBS_remove_age)
  
  all_SBS_burden_age_match <- cbind(all_SBS_burden_age_match, absolute_sig_contri_age_match[5])
  
  t_test_result_02 <- t.test(absolute_sig_contri_Disease_age_match$SBS_remove_age, absolute_sig_contri_Control_age_match$SBS_remove_age)
  wilcox.test_result_02 <- wilcox.test(absolute_sig_contri_Disease_age_match$SBS_remove_age, absolute_sig_contri_Control_age_match$SBS_remove_age, alternative = "two.sided")
  
  p_SNV_burden_per_sig_age_match <-ggplot(absolute_sig_contri_age_match, aes(label, SBS_remove_age, color = label, fill = label)) + 
    geom_boxplot(color = c("#599CB4", "#E69191"), fill = c("#C7DFF0", "#FBDFE2")) +
    geom_jitter(pch=21, color = "black", fill = c(disease_color_palette_02(12), control_color_palette_02(10)), size = 5, na.rm = FALSE) +
    # geom_segment(aes(x = 0.65,xend = 1.35,y = mean_absolute_sig_contri_control_remove, 
    #                  yend = mean_absolute_sig_contri_control_remove), colour = "dodgerblue3", size = 1) + 
    # geom_segment(aes(x = 1.65,xend = 2.35,y = mean_absolute_sig_contri_disease_remove, 
    #                  yend = mean_absolute_sig_contri_disease_remove), colour = "#FC4E07", size = 1) + 
    annotate("text", size = 7, x = 1.32, y = 1 * max(absolute_sig_contri_age_match$SBS_remove_age),
             # label = paste("P = ", format(wilcox.test_result_02$p.value, scientific = TRUE, digits = 3)), hjust = 0) +
             label = paste("P = ", round(wilcox.test_result_02$p.value, digits = 3)), hjust = 0) +
    labs(x = "Age", y = paste0(indel_type_list[i], "\n (sindel rate per cell)")) + 
    scale_color_manual(values = dis_ctrl_color) +
    stat_summary(fun=mean, colour="black", geom="point", shape=18, size=3, show.legend=FALSE) + 
    stat_compare_means(comparisons = list(c("Age-matched control", "Disease")), 
                       label.x = 1.6, label.y = 1 * max(absolute_sig_contri_age_match$SBS_remove_age), 
                       label = "p.signif", bracket.size = 1, size = 7) + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
          axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), legend.position = "none", panel.background = element_rect(fill = "white"), 
          axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank())
  print(p_SNV_burden_per_sig_age_match)
}
dev.off()

pdf(paste0(sindel_plot_dir, "/8-heatmap_top3_indel_type.pdf"), width = 30, height = 8)
all_SBS_burden <- all_SBS_burden[,-1]
all_SBS_burden$condition <- metadata_df$Condition
all_SBS_burden_conditional <- aggregate(all_SBS_burden[,plot_list], by=list(condition=all_SBS_burden$condition), FUN=mean)
rownames(all_SBS_burden_conditional) <- all_SBS_burden_conditional$condition
all_SBS_burden_conditional_matrix <- t(as.matrix(all_SBS_burden_conditional[,-1]))

colnames(all_SBS_burden)[plot_list] <- indel_type_list[plot_list]
row.names(all_SBS_burden_conditional_matrix) <- indel_type_list[plot_list]

all_SBS_burden_age_match <- all_SBS_burden_age_match[,-1]
all_SBS_burden_age_match$condition <- metadata_df_age_match$Condition
all_SBS_burden_age_match_conditional <- aggregate(all_SBS_burden_age_match[,plot_list], by=list(condition=all_SBS_burden_age_match$condition), FUN=mean)
rownames(all_SBS_burden_age_match_conditional) <- all_SBS_burden_age_match_conditional$condition
all_SBS_burden_age_match_conditional_matrix <- t(as.matrix(all_SBS_burden_age_match_conditional[,-1]))
colnames(all_SBS_burden_age_match)[plot_list] <- indel_type_list[plot_list]
row.names(all_SBS_burden_age_match_conditional_matrix) <- indel_type_list[plot_list]


print(plot_contribution_heatmap_cus(all_SBS_burden_conditional_matrix, cluster_sigs = FALSE, cluster_samples = FALSE) + 
        scale_fill_gradient(low = "white", high = "dodgerblue4") + geom_tile(colour = "black") + 
        geom_text(aes(label = round(Reshape(all_SBS_burden_conditional_matrix,nrow(all_SBS_burden_conditional_matrix)*2,1), 0)), size = 6) +
        theme(text = element_text(size=24), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.5)))

print(plot_contribution_heatmap_cus(all_SBS_burden_age_match_conditional_matrix, cluster_sigs = FALSE, cluster_samples = FALSE) + 
        scale_fill_gradient(low = "white", high = "dodgerblue4") + geom_tile(colour = "black") + 
        geom_text(aes(label = round(Reshape(all_SBS_burden_age_match_conditional_matrix,nrow(all_SBS_burden_age_match_conditional_matrix)*2,1), 0)), size = 6) +
        theme(text = element_text(size=24), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.5)))

melt_all_SBS_burden <- melt(all_SBS_burden)
colnames(melt_all_SBS_burden) <- c("condition", "sig", "no_mutation")
melt_all_SBS_burden$conditions <- factor(melt_all_SBS_burden$condition, level = c(control_name, disease_name))
barplot_melt_all_SBS_burden <- ggbarplot(melt_all_SBS_burden, x = "sig", y = "no_mutation", color = "conditions",
                                         label = TRUE, lab.size = 6, lab.nb.digits = 0, add = c("mean_se", "jitter"), position = position_dodge(0.9), 
                                         palette = dis_ctrl_color) +
  stat_compare_means(aes(group = condition, label = sprintf("p = %1.5f", as.numeric(..p.format..))), label.y = 1.04 * max(melt_all_SBS_burden$no_mutation), size = 6) + 
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.15 * max(melt_all_SBS_burden$no_mutation), size = 6) +
  labs(x = "Point mutation type", y = "number of sSNV") + 
  # scale_y_break(c(3500, 5800, 5920, 9800), scale=0.2) +
  scale_y_continuous(limits=c(0,83)) +
  theme_classic(base_size = 24) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white"), 
        axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank())
# facet_wrap(factor(mut_type, c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) ~ .)
print(barplot_melt_all_SBS_burden)

melt_all_SBS_burden_age_match <- melt(all_SBS_burden_age_match)
colnames(melt_all_SBS_burden_age_match) <- c("condition", "sig", "no_mutation")
melt_all_SBS_burden_age_match$conditions <- factor(melt_all_SBS_burden_age_match$condition, level = c(control_name, disease_name))
barplot_melt_all_SBS_burden_age_match <- ggbarplot(melt_all_SBS_burden_age_match, x = "sig", y = "no_mutation", color = "conditions",
                                                   label = TRUE, lab.size = 6, lab.nb.digits = 0, add = c("mean_se", "jitter"), position = position_dodge(0.9), 
                                                   palette = dis_ctrl_color) +
  stat_compare_means(aes(group = condition, label = sprintf("p = %1.3f", as.numeric(..p.format..))), label.y = 1.04 * max(melt_all_SBS_burden_age_match$no_mutation), size = 6) + 
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.18 * max(melt_all_SBS_burden_age_match$no_mutation), size = 6) +
  labs(x = "Point mutation type", y = "number of sSNV") + 
  # scale_y_break(c(3500, 5800, 5920, 9800), scale=0.2) +
  scale_y_continuous(limits=c(0,80)) +
  theme_classic(base_size = 24) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white"), 
        axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank())
# facet_wrap(factor(mut_type, c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) ~ .)
print(barplot_melt_all_SBS_burden_age_match)

all_mix_effect_model_p_value <- all_mix_effect_model_p_value[-1]
all_mix_effect_model_p_value <- round(as.numeric(all_mix_effect_model_p_value), digits = 3)
all_mix_effect_model_p_value <- as.data.frame(all_mix_effect_model_p_value)
rownames(all_mix_effect_model_p_value) <- indel_type_list[plot_list]
# all_mix_effect_model_p_value <- t(all_mix_effect_model_p_value)
colnames(all_mix_effect_model_p_value) <- "p_value"

stars.pval <- function(x){
  stars <- c("***", "**", "*", "ns")
  var <- c(0, 0.01, 0.05, 0.10, 1)
  i <- findInterval(x, var, left.open = T, rightmost.closed = T)
  stars[i]
}

a <- transform(all_mix_effect_model_p_value, stars = stars.pval(all_mix_effect_model_p_value[,1]))
t_a <- t(a)

dev.off()
