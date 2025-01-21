#############################################################################
######################### sSNV burden visualization #########################
#############################################################################

##### customized parameters for nice plot
yticks <- c(0, 300, 600, 900, 1200, 3450)
disease_color_palette<-colorRampPalette(c("pink1","firebrick3"))
control_color_palette<-colorRampPalette(c("skyblue1","dodgerblue4"))
dis_ctrl_color <- c(control_color_palette(9)[7], disease_color_palette(4)[3])
disease_color <- rep(disease_color_palette(5), c(3,3,3,3,3))
control_color <- rep(control_color_palette(8), c(2,1,3,3,1,4,2,2))
control_color_age_match <- rep(control_color_palette(3), c(1,4,2))

yticks_02 <- c(0, 500, 1000, 1500, 3000, 3200)
control_color_palette_02 <- colorRampPalette(c(dis_ctrl_color[1],dis_ctrl_color[1]))
disease_color_palette_02 <- colorRampPalette(c(dis_ctrl_color[2],dis_ctrl_color[2]))

scale_y_break <- c(1200, 3400)
y_limit <- c(0, 3500)
scale_y_break_02 <- c(1500, 2900)
y_limit_02 <- c(-70,3500)
y_pvalue <- 3250
y_psig <- 3200

#############################################################################
######################### Mix effects modeling ##############################
#############################################################################
sum(SCAN2_df$snv.burden[1:18])
sum(SCAN2_df$snv.burden[19:33])
sum(SCAN2_df$snv.burden)
##### input a matrix with columns: burden.per.gb, Case_ID, MAPD, CoV, age, condition
burden_df <- SCAN2_df[, c("snv.rate.per.gb", "Case_ID", "MAPD", "CoV", "Coverage", "Depth", "age", "condition")]
burden_df_age_match <- SCAN2_df_age_match[, c("snv.rate.per.gb", "Case_ID", "MAPD", "CoV", "Coverage", "Depth", "age", "condition")]

burden_df <- burden_df %>%
  mutate(condition = relevel(factor(condition), ref = "Control")) %>%
  group_by(factor(condition, levels = c(control_name, disease_name))) %>% arrange(age, .by_group = TRUE)
burden_df_control <- burden_df[burden_df$condition == "Control",]
burden_df_disease <- burden_df[burden_df$condition == "IHD",]

#############################################################################
##### mixed linear modeling (uncorrected)
sig_digits <- 2
for (plot_type in c("uncorrected", "MAPD_corrected", "CoV_corrected", "Coverage_corrected", "Depth_corrected")){
# for (plot_type in c("uncorrected")){
  if (plot_type == "uncorrected"){
    burden_age_model <- lmer(snv.rate.per.gb ~ age + condition + (1|Case_ID), burden_df, REML = FALSE)
    burden_age_model_control <- lmer(snv.rate.per.gb ~ age + (1|Case_ID), burden_df_control, REML = FALSE)
    burden_age_model_disease <- lmer(snv.rate.per.gb ~ age + (1|Case_ID), burden_df_disease, REML = FALSE)
    burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb
    } else if (plot_type == "MAPD_corrected"){
      burden_age_model <- lmer(snv.rate.per.gb ~ age + condition + MAPD + (1|Case_ID), burden_df, REML = FALSE)
      burden_age_model_control <- lmer(snv.rate.per.gb ~ age + MAPD + (1|Case_ID), burden_df_control, REML = FALSE)
      burden_age_model_disease <- lmer(snv.rate.per.gb ~ age + MAPD + (1|Case_ID), burden_df_disease, REML = FALSE)
      burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb - fixef(burden_age_model_control)[3] * burden_df$MAPD
      } else if (plot_type == "CoV_corrected"){
        burden_age_model <- lmer(snv.rate.per.gb ~ age + condition + CoV + (1|Case_ID), burden_df, REML = FALSE)
        burden_age_model_control <- lmer(snv.rate.per.gb ~ age + CoV + (1|Case_ID), burden_df_control, REML = FALSE)
        burden_age_model_disease <- lmer(snv.rate.per.gb ~ age + CoV + (1|Case_ID), burden_df_disease, REML = FALSE)
        burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb - fixef(burden_age_model_control)[3] * burden_df$CoV
        } else if (plot_type == "Coverage_corrected"){
          burden_age_model <- lmer(snv.rate.per.gb ~ age + condition + Coverage + (1|Case_ID), burden_df, REML = FALSE)
          burden_age_model_control <- lmer(snv.rate.per.gb ~ age + Coverage + (1|Case_ID), burden_df_control, REML = FALSE)
          burden_age_model_disease <- lmer(snv.rate.per.gb ~ age + Coverage + (1|Case_ID), burden_df_disease, REML = FALSE)
          burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb - fixef(burden_age_model_control)[3] * burden_df$Coverage
          } else if (plot_type == "Depth_corrected"){
            burden_age_model <- lmer(snv.rate.per.gb ~ age + condition + Depth + (1|Case_ID), burden_df, REML = FALSE)
            burden_age_model_control <- lmer(snv.rate.per.gb ~ age + Depth + (1|Case_ID), burden_df_control, REML = FALSE)
            burden_age_model_disease <- lmer(snv.rate.per.gb ~ age + Depth + (1|Case_ID), burden_df_disease, REML = FALSE)
            burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb - fixef(burden_age_model_control)[3] * burden_df$Depth
            } else {
              print("unknown correction type")
              }
  
  burden_df_control <- burden_df[burden_df$condition == "Control",]
  burden_df_disease <- burden_df[burden_df$condition == "IHD",]
  burden_age_model <- lmer(snv.rate.per.gb_cor ~ age + condition + (1|Case_ID), burden_df, REML = FALSE)
  burden_age_model_control <- lmer(snv.rate.per.gb_cor ~ age + (1|Case_ID), burden_df_control, REML = FALSE)
  burden_age_model_disease <- lmer(snv.rate.per.gb_cor ~ age + (1|Case_ID), burden_df_disease, REML = FALSE)
  
  # aging_rate <- formatC(signif(fixef(burden_age_model_control)[2], digits = 2), digits = 2, format="fg", flag="#")
  aging_rate <- format(round(fixef(burden_age_model_control)[2], digits = sig_digits), nsmall = sig_digits)
  
  ## check model parameters
  summary_lm_all <- summary(burden_age_model)
  summary_lm_control <- summary(burden_age_model_control)
  summary_lm_disease <- summary(burden_age_model_disease)
  
  ## compute confident interval
  CI_age_model <- confint(burden_age_model, oldNames=FALSE)
  CI_age_model_control <- confint(burden_age_model_control, oldNames=FALSE)
  
  ## check p value
  anova_pvalue <- anova(burden_age_model)$"Pr(>F)"[2]
  anova_pvalue_print <- formatC(signif(anova_pvalue, digits = 2), digits = 2, format="fg", flag="#")
  
  anova_pvalue_control <- anova(burden_age_model_control)$"Pr(>F)"[1]
  anova_pvalue_print_control <- formatC(signif(anova_pvalue_control, digits = 2), digits = 2, format="e", flag="#")
  mantissa <- format(as.numeric(sub("e.*", "", anova_pvalue_print_control)), nsmall = 2)
  exponent <- as.numeric(sub(".*e([+-]?\\d+)", "\\1", anova_pvalue_print_control))
  
  
  ## manually calculate the fitting lines
  Control_fitting_x = range(burden_df_control$age)
  Control_fitting_y = Control_fitting_x * fixef(burden_age_model_control)[2] + fixef(burden_age_model_control)[1]
  Disease_fitting_x = range(burden_df_disease$age)
  Disease_fitting_y = Disease_fitting_x * fixef(burden_age_model_disease)[2] + fixef(burden_age_model_disease)[1]
  
  ##### plot
  legend_data <- burden_df[c(7, 30), c("age", "snv.rate.per.gb_cor", "condition")]
  max_x <- max(burden_df$age)
  max_y <- max(burden_df$snv.rate.per.gb_cor)
  
  p_SNV_burden_lme <- ggplot(burden_df, aes(x = age, y = snv.rate.per.gb_cor), color = donors) + 
    geom_point(pch = 21, data = legend_data, aes(x = age, y = snv.rate.per.gb_cor, color = condition, fill = condition), size = 5) + 
    geom_point(pch = 21, color = "black", fill = c(control_color, disease_color), size = 5) + 
    geom_segment(aes(x = Control_fitting_x[1], xend = Control_fitting_x[2], y = Control_fitting_y[1], yend = Control_fitting_y[2]), colour = "dodgerblue3") + 
    # annotate("text", size = 7, x = 0.002 * max_x, y = 3470, label = bquote("s = " ~ .(aging_rate) ~ GB^{phantom()-1}*year^{phantom()-1}), hjust = 0, color = "black") + 
    annotate("text", size = 7, x = 0.0 * max_x, y = 3435, label = paste("p = ", anova_pvalue_print), hjust = 0, color = dis_ctrl_color[2]) + 
    annotate("text", size = 6, x = 0.0 * max_x, y = 3470, label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB"), hjust = 0, color = dis_ctrl_color[2]) +
    # annotate("text", size = 6, x = 0.0 * max_x, y = 3470, label = bquote("IHD effect: " ~ .(format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits)) ~ "sSNVs/GB"), hjust = 0, color = dis_ctrl_color[2]) +
    # annotate("text", size = 6, x = 0.0 * max_x, y = 3470, label = paste0("IHD excess ", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), " sSNVs"), hjust = 0, color = dis_ctrl_color[2]) +
    # annotate("text", size = 6, x = 0.0 * max_x, y = 1200, label = paste0("IHD excess SNVs = ", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), ", ",
    #                         "95% CI = [", format(round(CI_age_model[5, 1], digits = sig_digits), nsmall = sig_digits), ", ", format(round(CI_age_model[5, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) +
    # annotate("text", size = 6, x = 0.0 * max_x, y = 1130, label = paste0("Age accum. rate = ", format(round(fixef(burden_age_model_control)[2], digits = sig_digits), nsmall = sig_digits),
    #                         "/(GB.year), ", "95% CI = [", format(round(CI_age_model_control[4, 1], digits = sig_digits), nsmall = sig_digits), ", ", format(round(CI_age_model_control[4, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) +
    scale_color_manual(values = c("IHD" = "black", "Control" = "black"), guide = "legend") + 
    scale_fill_manual(values = c("Control" = dis_ctrl_color[1], "IHD" = dis_ctrl_color[2]), guide = "legend") + scale_y_break(scale_y_break, scale=0.25) + 
    scale_y_continuous(limits=y_limit, breaks=yticks, labels=yticks) + labs(x = "Age (years)", y = "sSNV rate per GB") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
          panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0))
  # ggsave(paste0(main_figure_dir, "/2-SNV_burden_lme_", plot_type, ".pdf"), plot = p_SNV_burden_lme, width = 9.5, height = 6, dpi = 600)
  ggsave(paste0(suppl_figure_dir, "/2-SNV_burden_lme_", plot_type, ".pdf"), plot = p_SNV_burden_lme, width = 9.5, height = 6, dpi = 600)
  # ggsave(paste0(other_figure_dir, "/2-SNV_burden_lme_", plot_type, ".pdf"), plot = p_SNV_burden_lme, width = 9.5, height = 6, dpi = 600)
  
  p_SNV_burden_lme_control <- ggplot(burden_df_control, aes(x = age, y = snv.rate.per.gb_cor), color = donors) + 
    geom_point(pch = 21, color = "black", fill = control_color, size = 5) +
    geom_segment(aes(x = Control_fitting_x[1], xend = Control_fitting_x[2], y = Control_fitting_y[1], yend = Control_fitting_y[2]), colour = "dodgerblue3") + 
    # annotate("text", size = 7, x = 0.002 * max_x, y = 500, label = bquote("aging effect: " ~ .(aging_rate) ~ "sSNVs" ~ GB^{phantom()-1}*year^{phantom()-1}), hjust = 0, color = dis_ctrl_color[1]) + 
    annotate("text", size = 7, x = 0.002 * max_x, y = 500, label = paste("aging effect:", aging_rate, "sSNVs/(GB·year)"), hjust = 0, color = dis_ctrl_color[1]) + 
    # annotate("text", size = 7, x = 0.002 * max_x, y = 500, label = bquote("s = " ~ .(aging_rate) ~ GB^{phantom()-1}*year^{phantom()-1}), hjust = 0, color = dis_ctrl_color[1]) + 
    annotate("text", size = 7, x = 0.0 * max_x, y = 460, label = bquote("p =" ~ .(mantissa) ~ "×" ~ 10^{.(exponent)}), hjust = 0, color = dis_ctrl_color[1]) + 
    # annotate("text", size = 6, x = 0.0 * max_x, y = 1200, label = paste0("IHD excess SNVs = ", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), ", ",
    #                         "95% CI = [", format(round(CI_age_model[5, 1], digits = sig_digits), nsmall = sig_digits), ", ", format(round(CI_age_model[5, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) +
    # annotate("text", size = 6, x = 0.0 * max_x, y = 1130, label = paste0("Age accum. rate = ", format(round(fixef(burden_age_model_control)[2], digits = sig_digits), nsmall = sig_digits),
    #                         "/(GB.year), ", "95% CI = [", format(round(CI_age_model_control[4, 1], digits = sig_digits), nsmall = sig_digits), ", ", format(round(CI_age_model_control[4, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) +
    scale_color_manual(values = c("IHD" = "black", "Control" = "black"), guide = "legend") + 
    scale_fill_manual(values = c("Control" = dis_ctrl_color[1], "IHD" = dis_ctrl_color[2]), guide = "legend") + 
    # scale_y_break(scale_y_break, scale=0.25) + 
    # scale_y_continuous(limits=y_limit, breaks=yticks, labels=yticks) + 
    labs(x = "Age (years)", y = "sSNV rate per GB") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
          panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0))
  p_SNV_burden_lme_control
  # ggsave(paste0(main_figure_dir, "/2-SNV_burden_lme_", plot_type, "_control.pdf"), plot = p_SNV_burden_lme_control, width = 7, height = 6, dpi = 600)
  ggsave(paste0(suppl_figure_dir, "/2-SNV_burden_lme_", plot_type, "_control.pdf"), plot = p_SNV_burden_lme_control, width = 9.5, height = 6, dpi = 600)
  }
