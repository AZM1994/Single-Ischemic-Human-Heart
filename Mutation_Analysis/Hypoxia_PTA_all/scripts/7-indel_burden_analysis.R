##### sindel burden visualization
sindel_plot_dir <- paste0(project_dir, "/sindel_plots")
dir.create(sindel_plot_dir)

pdf(paste0(sindel_plot_dir, "/7-sindel_burden_per_GB.pdf"), width = 10, height = 6)

p_indel_burden_per_GB <- ggplot(SCAN2_df, aes(x = factor(cell_ID, level = Cell_ID_list), y = indel.rate.per.gb, fill = condition)) + 
  geom_point(size = 6, pch = 21, na.rm = FALSE, show.legend = FALSE) + 
  scale_fill_manual(values = c("Normal" = dis_ctrl_color[1], "Disease" = dis_ctrl_color[2]), guide = "legend") +
  labs(x = "Sample ID", y = "sindel rate per GB") + theme_classic() +
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0), 
        axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right") + 
  facet_grid(. ~  factor(condition, level = c(control_name, disease_name)), scale = "free_x")
print(p_indel_burden_per_GB)

burden_df_indel <- SCAN2_df[, c("indel.rate.per.gb", "age", "condition", "Case_ID")]
burden_df_age_match_indel <- SCAN2_df_age_match[, c("indel.rate.per.gb", "age", "condition", "Case_ID")]
Mix_Effects_Model_Func_indel(burden_df_indel, burden_df_age_match_indel)

##### mix linear model fitting
# burden_age_model_indel <- lmer(indel.rate.per.gb ~ age + condition + (1|donor), SCAN2_df_relevel)
# burden_age_model_Disease_indel <- lmer(indel.rate.per.gb ~ age + (1|donor), SCAN2_df_relevel_Disease)
# burden_age_model_Control_indel <- lmer(indel.rate.per.gb ~ age + (1|donor), SCAN2_df_relevel_Control)
# summary_lm_all_indel <- summary(burden_age_model_indel)
# summary_lm_Disease_indel <- summary(burden_age_model_Disease_indel)
# summary_lm_Control_indel <- summary(burden_age_model_Control_indel)
# 
# anova_pvalue_indel <- anova(burden_age_model_indel)$"Pr(>F)"[2]
# r.squaredGLMM(burden_age_model_Disease_indel)[1]
# r.squaredGLMM(burden_age_model_Control_indel)[1]
# 
# ## manually calculate the fitting lines
# # Disease_fitting_x = range(SCAN2_df_relevel_Disease$age)
# Disease_fitting_y = Disease_fitting_x * fixef(burden_age_model_Disease_indel)[2] + fixef(burden_age_model_Disease_indel)[1]
# # Control_fitting_x = range(SCAN2_df_relevel_Control$age)
# Control_fitting_y = Control_fitting_x * fixef(burden_age_model_Control_indel)[2] + fixef(burden_age_model_Control_indel)[1]
# 
# # colfunc1<-colorRampPalette(c("pink","firebrick"))
# # colfunc2<-colorRampPalette(c("skyblue","dodgerblue4"))
# # disease_color <- rep(colfunc1(4), c(3,3,3,3))
# # control_color <- rep(colfunc2(9), c(2,1,3,3,3,1,4,2,2))
# # SCAN2_df_plot <- SCAN2_df_relevel %>% group_by(factor(condition, levels = c(disease_name, control_name))) %>% arrange(age, .by_group = TRUE)
# # SCAN2_df_plot$donors <- factor(SCAN2_df_plot$donor, levels = unique(SCAN2_df_plot$donor))
# 
# # yticks <- c(0, 500, 1000, 1500, 3400, 3450, 3500)
# legend_data_03 <- SCAN2_df_plot[c(7,30), c("age", "indel.rate.per.gb", "condition")]
# p_indel_burden_per_GB_age_origin <- ggplot(SCAN2_df_plot, aes(x = age, y = indel.rate.per.gb, color = donors)) +
#   geom_point(pch = 21, data = legend_data_03, aes(x = age, y = indel.rate.per.gb, color = condition, fill = condition), size = 5) +
#   geom_point(pch = 21, fill = c(disease_color, control_color), size = 5) +
#   geom_segment(aes(x = Disease_fitting_x[1], xend = Disease_fitting_x[2], 
#                    y = Disease_fitting_y[1], yend = Disease_fitting_y[2]), colour = "#FC4E07") + 
#   geom_segment(aes(x = Control_fitting_x[1], xend = Control_fitting_x[2], 
#                    y = Control_fitting_y[1], yend = Control_fitting_y[2]), colour = "dodgerblue3") + 
#   annotate("text", size = 7, x = 0.0 * max(SCAN2_df_relevel$age), y = 1.0 * max(SCAN2_df_relevel$indel.rate.per.gb),
#            label = paste("P = ", format(anova_pvalue_indel, scientific = FALSE, digits = 3)), hjust = 0) +
#   scale_color_manual(values = c("Disease" = "black", "Normal" = "black"), guide = "none") +
#   scale_fill_manual(values = c("Normal" = dis_ctrl_color[1], "Disease" = dis_ctrl_color[2]), guide = "legend") +
#   labs(x = "Age", y = "sSNV rate per GB") + theme_classic() +
#   theme(text = element_text(size=24), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
#         panel.background = element_rect(fill = "white"), legend.position="right")
# print(p_indel_burden_per_GB_age_origin)
# # ggsave(paste0(sSNV_plot_dir, "/high_quality/p_SNV_burden_per_GB_age_origin.tiff"), plot = p_SNV_burden_per_GB_age_origin, width = 10, height = 6, dpi = 600)
# # ggsave(paste0(sSNV_plot_dir, "/high_quality/p_SNV_burden_per_GB_age_origin.pdf"), plot = p_SNV_burden_per_GB_age_origin, width = 10, height = 6, dpi = 600)
# 
# ##### age matched mix linear model fitting
# burden_age_model_age_indel <- lmer(indel.rate.per.gb ~ age + condition + (1|donor), SCAN2_df_relevel_age_match)
# burden_age_model_Disease_age_indel <- lmer(indel.rate.per.gb ~ age + (1|donor), SCAN2_df_relevel_age_match_Disease)
# burden_age_model_Control_age_indel <- lmer(indel.rate.per.gb ~ age + (1|donor), SCAN2_df_relevel_age_match_Control)
# summary_lm_all_age_indel <- summary(burden_age_model_age_indel)
# summary_lm_Disease_age_indel <- summary(burden_age_model_Disease_age_indel)
# summary_lm_Control_age_indel <- summary(burden_age_model_Control_age_indel)
# 
# anova_pvalue_age_indel <- anova(burden_age_model_age_indel)$"Pr(>F)"[2]
# r.squaredGLMM(burden_age_model_Disease_age_indel)[1]
# r.squaredGLMM(burden_age_model_Control_age_indel)[1]
# 
# ##### remove the effect of age
# SCAN2_df_relevel_age_match_Control$indel.rate.per.gb_remove_age <- SCAN2_df_relevel_age_match_Control$indel.rate.per.gb - 
#   fixef(burden_age_model_Control_indel)[2] * SCAN2_df_relevel_age_match_Control$age
# burden_age_model_Control_remove_age <- lmer(indel.rate.per.gb_remove_age ~ age + (1|donor), SCAN2_df_relevel_age_match_Control)
# summary_lm_Control_remove_age <- summary(burden_age_model_Control_remove_age)
# 
# SCAN2_df_relevel_age_match_Disease$indel.rate.per.gb_remove_age <- SCAN2_df_relevel_age_match_Disease$indel.rate.per.gb - 
#   fixef(burden_age_model_Control_indel)[2] * SCAN2_df_relevel_age_match_Disease$age
# burden_age_model_Disease_remove_age <- lmer(indel.rate.per.gb_remove_age ~ age + (1|donor), SCAN2_df_relevel_age_match_Disease)
# summary_lm_Disease_remove_age <- summary(burden_age_model_Disease_remove_age)
# 
# SCAN2_df_relevel_age_match$indel.rate.per.gb_remove_age <-
#   c(SCAN2_df_relevel_age_match_Disease$indel.rate.per.gb_remove_age, SCAN2_df_relevel_age_match_Control$indel.rate.per.gb_remove_age)
# SCAN2_df_relevel_age_match$label <- c(SCAN2_df_relevel_age_match_Disease$label, SCAN2_df_relevel_age_match_Control$label)
# mean_indel_disease_remove <- mean(SCAN2_df_relevel_age_match_Disease$indel.rate.per.gb_remove_age)
# mean_indel_control_remove <- mean(SCAN2_df_relevel_age_match_Control$indel.rate.per.gb_remove_age)
# 
# t_test_result_indel <- t.test(SCAN2_df_relevel_age_match_Disease$indel.rate.per.gb_remove_age, SCAN2_df_relevel_age_match_Control$indel.rate.per.gb_remove_age)
# wilcox.test_result_indel <- wilcox.test(SCAN2_df_relevel_age_match_Disease$indel.rate.per.gb_remove_age, 
#                                   SCAN2_df_relevel_age_match_Control$indel.rate.per.gb_remove_age, alternative = "two.sided")
# # yticks_02 <- c(0, 500, 1000, 1500, 3250, 3300, 3350)
# # colfunc3<-colorRampPalette(c("firebrick2","firebrick2"))
# # colfunc4<-colorRampPalette(c("dodgerblue3","dodgerblue3"))
# 
# p32 <- ggplot(SCAN2_df_relevel_age_match, aes(label, indel.rate.per.gb_remove_age, color = label)) +
#   geom_boxplot(color = c("#599CB4", "#E69191"), fill = c("#C7DFF0", "#FBDFE2")) +
#   geom_jitter(pch=21, color = "black", fill = c(disease_color_palette_02(12), control_color_palette_02(10)), size = 5, na.rm = FALSE) +
#   # geom_segment(aes(x = 0.65,xend = 1.35,y = mean_indel_control_remove, yend = mean_indel_control_remove), colour = "dodgerblue3", size = 1) + 
#   # geom_segment(aes(x = 1.65,xend = 2.35,y = mean_indel_disease_remove, yend = mean_indel_disease_remove), colour = "#FC4E07", size = 1) + 
#   annotate("text", size = 7, x = 1.35, y = 125, label = paste("P = ", format(wilcox.test_result_indel$p.value, scientific = FALSE, digits = 3)), hjust = 0) +
#   labs(x = "", y = "Age-corrected sindel rate per GB") + 
#   scale_color_manual(values = dis_ctrl_color) + 
#   stat_summary(fun=mean, colour="black", geom="point", shape=18, size=3, show.legend=FALSE) + 
#   stat_compare_means(comparisons = list(c("Age-matched control (n = 10)", "Disease (n = 12)")), label.x = 1.6, label.y = 125,
#                      label = "p.signif", bracket.size = 1, size = 7) + theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
#         axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), legend.position = "none", panel.background = element_rect(fill = "white"), 
#         axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank())
# print(p32)

dev.off()

# means <- aggregate(indel.rate.per.gb ~  condition, SCAN2_df, FUN = function(x) {round(mean(x), digits=0)})
# boxplot_sindel_burden <- ggplot(SCAN2_df, aes(x=factor(condition, level=c(control_name, disease_name)), y=indel.rate.per.gb, fill=condition)) +
#   geom_boxplot() + scale_fill_manual(values=c("darkorange", "cyan3"), limits = c(control_name, disease_name)) +
#   geom_jitter(shape = 21, colour = "black", fill = "white", size = 2, stroke = 1, width = 0.2) +
#   stat_summary(fun=mean, colour="darkred", geom="point", 
#                shape=18, size=3, show.legend=FALSE) + 
#   geom_text(data = means, aes(label = indel.rate.per.gb, y = indel.rate.per.gb, hjust = 1.5)) +
#   labs(x = "condition", y = "sindel rate per GB") + 
#   theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), 
#         axis.title.x = element_text(size = 25, hjust = 0.5), axis.title.y = element_text(size = 25, hjust = 0.5)) + 
#   stat_compare_means(method = "t.test") + 
#   ggtitle("Boxplot of estimated sindel/GB")
# ggsave(paste0(sindel_plot_dir, "/boxplot_sindel_burden.png"), plot = boxplot_sindel_burden, width = 10, height = 8)
