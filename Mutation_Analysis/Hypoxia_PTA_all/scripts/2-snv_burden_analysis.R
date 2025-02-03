#############################################################################
######################### sSNV burden visualization #########################
#############################################################################

yticks_all <- c(0, 100, 400, 900, 1600, 2500, 3600, 4900, 6400, 8100)
y_limit_all <- c(0, 8100)
yticks_all_02 <- c(-400, -100, 0, 100, 400, 900, 1600, 2500, 3600, 4900, 6400, 8100)
y_limit_all_02 <- c(-400, 8100)

yticks_control <- c(0, 100, 400, 900, 1600, 2500, 3600)
y_limit_control <- c(0, 3600)
yticks_control_02 <- c(-400, -100, 0, 100, 400, 900, 1600, 2500)
y_limit_control_02 <- c(-400, 2500)

#############################################################################
######################### Mix effects modeling ##############################
#############################################################################
sum(SCAN2_df$snv.burden[1:50])
sum(SCAN2_df$snv.burden[51:83])
sum(SCAN2_df$snv.burden)

##### input a matrix with columns: burden.per.gb, Case_ID, MAPD, CoV, Age, Condition
burden_df <- SCAN2_df[, c("snv.rate.per.gb", "Case_ID", "MAPD", "CoV", "Coverage", "Depth", "Age", "Condition", "Gender", "Color")] %>%
  mutate(Condition = relevel(factor(Condition), ref = "Control")) %>% group_by(factor(Condition, levels = c(control_name, disease_name))) %>% arrange(Age, .by_group = TRUE)
burden_df_age_match <- SCAN2_df_age_match[, c("snv.rate.per.gb", "Case_ID", "MAPD", "CoV", "Coverage", "Depth", "Age", "Condition", "Gender")]
burden_df_control <- burden_df[burden_df$Condition == "Control", ]
burden_df_disease <- burden_df[burden_df$Condition == "IHD", ]

#############################################################################
##### mixed linear modeling
sig_digits <- 2
# for (plot_type in c("uncorrected", "MAPD_corrected", "CoV_corrected", "Coverage_corrected", "Depth_corrected")){
for (plot_type in c("uncorrected")){
# for (plot_type in c("Coverage_corrected")){
  y_limit_plot <- y_limit_all
  yticks_plot <- yticks_all
  y_limit_control_plot <- y_limit_control
  yticks_control_plot <- yticks_control
  if (plot_type == "uncorrected"){
    burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + (1|Case_ID), burden_df, REML = FALSE)
    burden_age_model_control <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_control, REML = FALSE)
    burden_age_model_disease <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_disease, REML = FALSE)
    burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb
    ### for gender testing
    # burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + gender + (1|Case_ID), burden_df, REML = FALSE)
    # burden_age_model_control <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_control, REML = FALSE)
    # burden_age_model_disease <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_disease, REML = FALSE)
    # burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb
    } else if (plot_type == "MAPD_corrected"){
      burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + MAPD + (1|Case_ID), burden_df, REML = FALSE)
      burden_age_model_control <- lmer(snv.rate.per.gb ~ Age + MAPD + (1|Case_ID), burden_df_control, REML = FALSE)
      burden_age_model_disease <- lmer(snv.rate.per.gb ~ Age + MAPD + (1|Case_ID), burden_df_disease, REML = FALSE)
      burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb - fixef(burden_age_model_control)[3] * burden_df$MAPD
      } else if (plot_type == "CoV_corrected"){
        burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + CoV + (1|Case_ID), burden_df, REML = FALSE)
        burden_age_model_control <- lmer(snv.rate.per.gb ~ Age + CoV + (1|Case_ID), burden_df_control, REML = FALSE)
        burden_age_model_disease <- lmer(snv.rate.per.gb ~ Age + CoV + (1|Case_ID), burden_df_disease, REML = FALSE)
        burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb - fixef(burden_age_model_control)[3] * burden_df$CoV
        } else if (plot_type == "Coverage_corrected"){
          burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + Coverage + (1|Case_ID), burden_df, REML = FALSE)
          burden_age_model_control <- lmer(snv.rate.per.gb ~ Age + Coverage + (1|Case_ID), burden_df_control, REML = FALSE)
          burden_age_model_disease <- lmer(snv.rate.per.gb ~ Age + Coverage + (1|Case_ID), burden_df_disease, REML = FALSE)
          burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb - fixef(burden_age_model_control)[3] * burden_df$Coverage
          y_limit_plot <- y_limit_all_02
          yticks_plot <- yticks_all_02
          y_limit_control_plot <- y_limit_control_02
          yticks_control_plot <- yticks_control_02
          } else if (plot_type == "Depth_corrected"){
            burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + Depth + (1|Case_ID), burden_df, REML = FALSE)
            burden_age_model_control <- lmer(snv.rate.per.gb ~ Age + Depth + (1|Case_ID), burden_df_control, REML = FALSE)
            burden_age_model_disease <- lmer(snv.rate.per.gb ~ Age + Depth + (1|Case_ID), burden_df_disease, REML = FALSE)
            burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb - fixef(burden_age_model_control)[3] * burden_df$Depth
            } else {
              print("unknown correction type")
              }
  
  burden_df_control <- burden_df[burden_df$Condition == "Control",]
  burden_df_disease <- burden_df[burden_df$Condition == "IHD",]
  burden_age_model <- lmer(snv.rate.per.gb_cor ~ Age + Condition + (1|Case_ID), burden_df, REML = FALSE)
  # burden_age_model <- lmer(snv.rate.per.gb_cor ~ Age + Condition + gender + (1|Case_ID), burden_df, REML = FALSE)
  burden_age_model_control <- lmer(snv.rate.per.gb_cor ~ Age + (1|Case_ID), burden_df_control, REML = FALSE)
  burden_age_model_disease <- lmer(snv.rate.per.gb_cor ~ Age + (1|Case_ID), burden_df_disease, REML = FALSE)
  
  # aging_rate <- formatC(signif(fixef(burden_age_model_control)[2], digits = 2), digits = 2, format="fg", flag="#")
  aging_rate <- format(round(fixef(burden_age_model_control)[2], digits = sig_digits), nsmall = sig_digits)
  
  ## check model parameters
  summary_lm_all <- summary(burden_age_model)
  summary_lm_control <- summary(burden_age_model_control)
  summary_lm_disease <- summary(burden_age_model_disease)
  
  ## compute confident interval
  CI_age_model <- confint(burden_age_model, oldNames = FALSE)
  CI_age_model_control <- confint(burden_age_model_control, method = "Wald", oldNames = FALSE)
  
  ## check p value
  anova_pvalue <- anova(burden_age_model)$"Pr(>F)"[2]
  anova_pvalue_print <- formatC(signif(anova_pvalue, digits = 2), digits = 2, format="fg", flag="#")
  
  anova_pvalue_control <- anova(burden_age_model_control)$"Pr(>F)"[1]
  anova_pvalue_print_control <- formatC(signif(anova_pvalue_control, digits = 2), digits = 2, format="fg", flag="#")
  
  ## manually calculate the fitting lines
  geom_line_data <- burden_df_control %>% {
    Control_fitting_x <- range(.$Age)
    Control_fitting_y <- Control_fitting_x * fixef(burden_age_model_control)[2] + fixef(burden_age_model_control)[1]
    geom_line_data <- rbind(Control_fitting_x, Control_fitting_y)
    as.data.frame(t(geom_line_data))}
  
  ##### plot
  legend_data <- burden_df[c(7, 60), c("Age", "snv.rate.per.gb_cor", "Condition")]
  max_x <- max(burden_df$Age)
  max_y <- max(burden_df$snv.rate.per.gb_cor)
  
  p_SNV_burden_lme <- ggplot(burden_df, aes(x = Age, y = snv.rate.per.gb_cor)) + 
    geom_point(pch = 21, data = legend_data, aes(x = Age, y = snv.rate.per.gb_cor, color = Condition, fill = Condition), size = 5) + 
    geom_point(pch = 21, color = "black", fill = burden_df$Color, size = 5) + 
    geom_line(aes(x = Control_fitting_x, y = Control_fitting_y), color = "dodgerblue3", data = geom_line_data) + 
    annotate("text", size = 7, x = 0, y = 8000, label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB"), hjust = 0) +
    annotate("text", size = 7, x = 0, y = 7000, label = paste("P = ", anova_pvalue_print), hjust = 0) + 
    # annotate("text", size = 5, x = 0.0 * max_x, y = 5500, label = paste0("Age accum. rate = ", format(round(fixef(burden_age_model_control)[2], digits = sig_digits), nsmall = sig_digits), "/(GB.year), ", "95% CI = [", format(round(CI_age_model_control[4, 1], digits = sig_digits), nsmall = sig_digits), ", ", format(round(CI_age_model_control[4, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) + 
    # annotate("text", size = 5, x = 0.0 * max_x, y = 4500, label = paste0("IHD excess SNVs = ", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), ", ", "95% CI = [", format(round(CI_age_model[5, 1], digits = sig_digits), nsmall = sig_digits), ", ", format(round(CI_age_model[5, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) + 
    scale_color_manual(values = c("Control" = "black", "IHD" = "black"), guide = "legend") + scale_fill_manual(values = c("Control" = ctrl_dis_color[1], "IHD" = ctrl_dis_color[2]), guide = "legend") +
    scale_y_continuous(limits = y_limit_plot, breaks = yticks_plot, labels = yticks_plot, trans = scales::modulus_trans(0.5)) + labs(x = "Age (years)", y = "sSNV rate per GB") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
          panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
  ggsave(paste0(main_figure_dir, "/2-SNV_burden_lme_", plot_type, ".pdf"), plot = p_SNV_burden_lme, width = 9.5, height = 6, dpi = 600)
  # ggsave(paste0(suppl_figure_dir, "/2-SNV_burden_lme_", plot_type, ".pdf"), plot = p_SNV_burden_lme, width = 9.5, height = 6, dpi = 600)
  # ggsave(paste0(other_figure_dir, "/2-SNV_burden_lme_", plot_type, ".pdf"), plot = p_SNV_burden_lme, width = 9.5, height = 6, dpi = 600)
  
  p_SNV_burden_lme_control <- ggplot(burden_df_control, aes(x = Age, y = snv.rate.per.gb_cor)) + 
    geom_point(pch = 21, color = "black", fill = burden_df_control$Color, size = 5) + 
    geom_line(aes(x = Control_fitting_x, y = Control_fitting_y), color = "dodgerblue3", data = geom_line_data) + 
    annotate("text", size = 7, x = 0, y = 3500, label = paste("aging effect:", aging_rate, "sSNVs/(GB·year)"), hjust = 0) + 
    annotate("text", size = 7, x = 0, y = 3000, label = paste("P = ", anova_pvalue_print_control), hjust = 0) + 
    scale_y_continuous(limits = y_limit_control_plot, breaks = yticks_control_plot, labels = yticks_control_plot, trans = scales::modulus_trans(0.5)) + labs(x = "Age (years)", y = "sSNV rate per GB") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
          panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
  ggsave(paste0(main_figure_dir, "/2-SNV_burden_lme_", plot_type, "_control.pdf"), plot = p_SNV_burden_lme_control, width = 7, height = 6, dpi = 600)
  # ggsave(paste0(suppl_figure_dir, "/2-SNV_burden_lme_", plot_type, "_control.pdf"), plot = p_SNV_burden_lme_control, width = 9.5, height = 6, dpi = 600)
  }
