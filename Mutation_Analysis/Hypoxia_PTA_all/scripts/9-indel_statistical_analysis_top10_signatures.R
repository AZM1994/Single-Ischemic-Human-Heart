##### statistical_analysis_top10_signatures
##### statistical_analysis_top10_signatures
## raw contribution
top10_contri_mat <- as.data.frame(t(refined_new_fit_res[-nrow(refined_new_fit_res),]))
## normalize total contribution for each cell to 100
top10_contri_mat <- 1 / rowSums(top10_contri_mat) * top10_contri_mat * SCAN2_df_age_match$indel.burden
# rowSums(top10_contri_mat)
# rowSums(top10_contri_mat[,-ncol(top10_contri_mat)])

control_sample_list <- metadata_df$Cell_ID[metadata_df$Condition == control_name]
disease_sample_list <- metadata_df$Cell_ID[metadata_df$Condition == disease_name]
control_condition <- grepl(paste(control_sample_list, collapse = "|"), rownames(top10_contri_mat))
disease_condition <- grepl(paste(disease_sample_list, collapse = "|"), rownames(top10_contri_mat))
top10_contri_mat$condition <- with(top10_contri_mat, ifelse(control_condition, control_name, 
                                                            ifelse(disease_condition, disease_name, 'unknown')))

top10_contri_mat_temp <- top10_contri_mat

melt_matrix_statistics <- melt(top10_contri_mat)
colnames(melt_matrix_statistics)[2:3] <- c("sigs", "contribution")

pdf(paste0(sindel_plot_dir, "/9-indel_statistical_analysis_top10_signatures.pdf"), width = 12, height = 6)

top10_contri_mat_conditional <- top10_contri_mat[,-ncol(top10_contri_mat)]
Control <- colSums(top10_contri_mat_conditional[control_range_age_match,]) / length(control_range_age_match)
top10_contri_mat_conditional <- rbind(top10_contri_mat_conditional, Control)
Disease <- colSums(top10_contri_mat_conditional[disease_range_age_match,]) / length(disease_range_age_match)
top10_contri_mat_conditional <- rbind(top10_contri_mat_conditional, Disease)
Net_change <- Disease - Control
# Net_change[Net_change<0] <- NA
top10_contri_mat_conditional <- rbind(top10_contri_mat_conditional, Net_change)
rownames(top10_contri_mat_conditional)[(nrow(top10_contri_mat_conditional)-2):nrow(top10_contri_mat_conditional)] <- c(control_name, disease_name, 'Net Change')

haha <- t(top10_contri_mat_conditional[(nrow(top10_contri_mat_conditional) - 2):(nrow(top10_contri_mat_conditional) - 1),])
haha <- haha[haha[,1] + haha[,2] > 0 ,]
melt_matrix_statistics <- melt_matrix_statistics[melt_matrix_statistics$sigs %in% rownames(haha),]

melt_matrix_statistics$conditions <- factor(melt_matrix_statistics$condition, level=c(control_name, disease_name))
boxplot_statistical_analysis_top10_sigs_conditional <- ggboxplot(melt_matrix_statistics, x = "sigs", y = "contribution", color = "conditions") +
  # geom_jitter(aes(colour = condition), shape = 21, size = 1, stroke = 1, width = 0.2) +
  stat_compare_means(aes(group = condition), label = "p.format") + 
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 950) +
  labs(x = "top10 extracted signatures", y = "Absolute Contribution") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20),
        axis.title.x = element_text(size = 25, hjust = 0.5), axis.title.y = element_text(size = 25, hjust = 0.5)) + 
  scale_fill_manual(values=c("darkorange", "cyan3")) + theme_gray(base_size = 20)
print(boxplot_statistical_analysis_top10_sigs_conditional)

# round(Reshape(refined_new_fit_res_conditional_matrix,nrow(refined_new_fit_res_conditional_matrix),1), 2)
print(plot_contribution_heatmap_cus(haha, cluster_sigs = FALSE, cluster_samples = FALSE) + 
        scale_fill_gradient(low = "white", high = "dodgerblue4") + geom_tile(colour = "black") + 
        # geom_text(aes(label = round(Reshape(haha,nrow(haha)*2,1), 0)), size = 6) +
        theme(text = element_text(size=15), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5)))

dev.off()

######### mix linear model for each signature
top10_contri_mat_all_cell <- as.data.frame(t(fit_res_strict_conditional_02_all[(row.names(fit_res_strict_conditional_02_all) %in% refined_selected_sigs_list),]))
## normalize total contribution for each cell to 1
top10_contri_mat_all_cell <- 1 / rowSums(top10_contri_mat_all_cell) * top10_contri_mat_all_cell * SCAN2_df$indel.rate.per.gb

control_condition <- grepl(paste(control_sample_list, collapse = "|"), rownames(top10_contri_mat_all_cell))
disease_condition <- grepl(paste(disease_sample_list, collapse = "|"), rownames(top10_contri_mat_all_cell))
top10_contri_mat_all_cell$condition <- with(top10_contri_mat_all_cell, ifelse(control_condition, control_name, 
                                                                              ifelse(disease_condition, disease_name, 'unknown')))

top10_contri_mat_all_cell$age <- metadata_df$Age
top10_contri_mat_all_cell$conditions <- metadata_df$Condition
top10_contri_mat_all_cell$gender <- metadata_df$Gender
top10_contri_mat_all_cell$Case_ID <- metadata_df$Case_ID
top10_contri_mat_all_cell <- top10_contri_mat_all_cell %>% mutate(condition = relevel(factor(condition), ref = "Normal"))
top10_contri_mat_age_match <- top10_contri_mat_all_cell[Cell_ID_list_age_match,]
# top10_contri_mat_all_cell <- top10_contri_mat_all_cell[-9,]
# top10_contri_mat_age_match <- top10_contri_mat_age_match[-9,]

# refined_selected_sigs_list_02 <- refined_selected_sigs_list[c(1,3,5,6)]
refined_selected_sigs_list_02 <- refined_selected_sigs_list
plot_list <- c(1:2)
all_SBS_burden <- metadata_df$Cell_ID
all_SBS_burden_age_match <- metadata_df_age_match$Cell_ID

all_mix_effect_model_p_value <- metadata_df[1,1]

pdf(paste0(sindel_plot_dir, "/9-indel_top10_signatures_contri_with_age_condition.pdf"), width = 12, height = 8)
for (i in plot_list){
  ##### all cell all age
  absolute_sig_contri <- top10_contri_mat_all_cell[c(refined_selected_sigs_list_02[i],'age','condition','Case_ID')]
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
  legend_data_02 <- absolute_sig_contri_plot[c(7,30), c("age", "SBS", "condition")]
  p_SNV_burden_per_sig <- ggplot(absolute_sig_contri_plot, aes(x = age, y = SBS, color = Case_IDs)) +
    geom_point(pch = 21, data = legend_data_02, aes(x = age, y = SBS, color = condition, fill = condition), size = 5) +
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
    # scale_y_break(c(1500, 3400), scale=0.25) + scale_y_continuous(limits=c(0,3500), breaks=yticks, labels=yticks) +
    # ggtitle("Estimated sSNV rate per GB for each sample") +
    labs(x = "Age", y = paste0(refined_selected_sigs_list_02[i], " Contribution \n (sindel rate per GB)")) + theme_classic() +
    theme(text = element_text(size=24), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position="right")
  print(p_SNV_burden_per_sig)
  
  ##### age matched cells
  absolute_sig_contri_age_match <- (top10_contri_mat_age_match[c(refined_selected_sigs_list_02[i],'age','condition','Case_ID')])
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
  
  p_SNV_burden_per_sig_age_match <-ggplot(absolute_sig_contri_age_match, aes(label, SBS_remove_age, color = label)) + 
    geom_boxplot(color = c("#599CB4", "#E69191"), fill = c("#C7DFF0", "#FBDFE2")) +
    geom_jitter(pch=21, color = "black", fill = c(disease_color_palette_02(12), control_color_palette_02(10)), size = 5, na.rm = FALSE) +
    # geom_segment(aes(x = 0.65,xend = 1.35,y = mean_absolute_sig_contri_control_remove, 
    #                  yend = mean_absolute_sig_contri_control_remove), colour = "dodgerblue3", size = 1) + 
    # geom_segment(aes(x = 1.65,xend = 2.35,y = mean_absolute_sig_contri_disease_remove, 
    #                  yend = mean_absolute_sig_contri_disease_remove), colour = "#FC4E07", size = 1) + 
    annotate("text", size = 7, x = 1.32, y = 1 * max(absolute_sig_contri_age_match$SBS_remove_age),
             # label = paste("P = ", format(wilcox.test_result_02$p.value, scientific = TRUE, digits = 3)), hjust = 0) +
             label = paste("P = ", round(wilcox.test_result_02$p.value, digits = 3)), hjust = 0) +
    labs(x = "Age", y = paste0(refined_selected_sigs_list_02[i], " Contribution \n (Age-corrected sindel rate per GB)")) + 
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

pdf(paste0(sindel_plot_dir, "/9-heatmap_top8_signatures.pdf"), width = 15, height = 8)
all_SBS_burden <- all_SBS_burden[,-1]
all_SBS_burden$condition <- metadata_df$Condition
all_SBS_burden_conditional <- aggregate(all_SBS_burden[,plot_list], by=list(condition=all_SBS_burden$condition), FUN=mean)
rownames(all_SBS_burden_conditional) <- all_SBS_burden_conditional$condition
all_SBS_burden_conditional_matrix <- t(as.matrix(all_SBS_burden_conditional[,-1]))

colnames(all_SBS_burden)[plot_list] <- refined_selected_sigs_list_02[plot_list]
row.names(all_SBS_burden_conditional_matrix) <- refined_selected_sigs_list_02[plot_list]

all_SBS_burden_age_match <- all_SBS_burden_age_match[,-1]
all_SBS_burden_age_match$condition <- metadata_df_age_match$Condition
all_SBS_burden_age_match_conditional <- aggregate(all_SBS_burden_age_match[,plot_list], by=list(condition=all_SBS_burden_age_match$condition), FUN=mean)
rownames(all_SBS_burden_age_match_conditional) <- all_SBS_burden_age_match_conditional$condition
all_SBS_burden_age_match_conditional_matrix <- t(as.matrix(all_SBS_burden_age_match_conditional[,-1]))
colnames(all_SBS_burden_age_match)[plot_list] <- refined_selected_sigs_list_02[plot_list]
row.names(all_SBS_burden_age_match_conditional_matrix) <- refined_selected_sigs_list_02[plot_list]


print(plot_contribution_heatmap_cus(all_SBS_burden_conditional_matrix, cluster_sigs = FALSE, cluster_samples = FALSE) + 
        scale_fill_gradient(low = "white", high = "dodgerblue4") + geom_tile(colour = "black") + 
        geom_text(aes(label = round(Reshape(all_SBS_burden_conditional_matrix,nrow(all_SBS_burden_conditional_matrix)*2,1), 0)), size = 6) +
        theme(text = element_text(size=24), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5)))

print(plot_contribution_heatmap_cus(all_SBS_burden_age_match_conditional_matrix, cluster_sigs = FALSE, cluster_samples = FALSE) + 
        scale_fill_gradient(low = "white", high = "dodgerblue4") + geom_tile(colour = "black") + 
        geom_text(aes(label = round(Reshape(all_SBS_burden_age_match_conditional_matrix,nrow(all_SBS_burden_age_match_conditional_matrix)*2,1), 0)), size = 6) +
        theme(text = element_text(size=24), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5)))

melt_all_SBS_burden <- melt(all_SBS_burden)
colnames(melt_all_SBS_burden) <- c("condition", "sig", "no_mutation")
melt_all_SBS_burden$conditions <- factor(melt_all_SBS_burden$condition, level = c(control_name, disease_name))
barplot_melt_all_SBS_burden <- ggbarplot(melt_all_SBS_burden, x = "sig", y = "no_mutation", color = "conditions",
                                         label = TRUE, lab.size = 6, lab.nb.digits = 0, add = c("mean_se", "jitter"), position = position_dodge(0.9), 
                                         palette = dis_ctrl_color) +
  stat_compare_means(aes(group = condition, label = sprintf("p = %1.5f", as.numeric(..p.format..))), label.y = 1.04 * max(melt_all_SBS_burden$no_mutation), size = 6) + 
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.08 * max(melt_all_SBS_burden$no_mutation), size = 6) +
  labs(x = "Point mutation type", y = "number of sSNV") + 
  # scale_y_break(c(3500, 5800, 5920, 9800), scale=0.2) +
  # scale_y_continuous(limits=c(0,11500), breaks=yticks_03, labels=yticks_03) +
  theme_classic(base_size = 24) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white"), 
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
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.08 * max(melt_all_SBS_burden_age_match$no_mutation), size = 6) +
  labs(x = "Point mutation type", y = "number of sSNV") + 
  # scale_y_break(c(3500, 5800, 5920, 9800), scale=0.2) +
  # scale_y_continuous(limits=c(0,11500), breaks=yticks_03, labels=yticks_03) +
  theme_classic(base_size = 24) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white"), 
        axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank())
# facet_wrap(factor(mut_type, c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) ~ .)
print(barplot_melt_all_SBS_burden_age_match)

all_mix_effect_model_p_value <- all_mix_effect_model_p_value[-1]
all_mix_effect_model_p_value <- round(as.numeric(all_mix_effect_model_p_value), digits = 3)
all_mix_effect_model_p_value <- as.data.frame(all_mix_effect_model_p_value)
rownames(all_mix_effect_model_p_value) <- refined_selected_sigs_list_02[plot_list]
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