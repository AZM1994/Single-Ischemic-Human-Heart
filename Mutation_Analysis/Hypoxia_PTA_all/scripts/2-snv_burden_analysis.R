#############################################################################
##################### sSNV burden mixed-effects modeling ####################
#############################################################################
sum(SCAN2_df$snv.burden[ctrl_range])
sum(SCAN2_df$snv.burden[dis_range])
sum(SCAN2_df$snv.burden)

##### input a matrix with columns: burden.per.gb, Case_ID, MAPD, CoV, Age, Condition
burden_df <- SCAN2_df[, c("snv.rate.per.gb", "Case_ID", "MAPD", "CoV", "Coverage", "Depth", "Age", "Condition", "Gender", "Color")] %>%
  mutate(Condition = relevel(factor(Condition), ref = "Control")) %>% group_by(factor(Condition, levels = c(ctrl_name, dis_name))) %>% arrange(Age, .by_group = TRUE)
burden_df_ctrl <- burden_df %>% filter(Condition == "Control")
burden_df_dis <- burden_df %>% filter(Condition == "IHD")

#############################################################################
##### mixed linear modeling
sig_digits <- 2
# for (plot_type in c("uncorrected", "MAPD_corrected", "CoV_corrected", "Coverage_corrected", "Depth_corrected")){
for (plot_type in c("uncorrected")){
  if (plot_type == "uncorrected"){
    burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + (1|Case_ID), burden_df, REML = FALSE)
    burden_age_model_ctrl <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_ctrl, REML = FALSE)
    burden_age_model_dis <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_dis, REML = FALSE)
    burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb
    ### for gender testing: P = 0.1334
    # burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + Gender + (1|Case_ID), burden_df, REML = FALSE)
    # burden_age_model_ctrl <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_ctrl, REML = FALSE)
    # burden_age_model_dis <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_dis, REML = FALSE)
    # burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb
    } else if (plot_type == "MAPD_corrected"){
      burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + MAPD + (1|Case_ID), burden_df, REML = FALSE)
      burden_age_model_ctrl <- lmer(snv.rate.per.gb ~ Age + MAPD + (1|Case_ID), burden_df_ctrl, REML = FALSE)
      burden_age_model_dis <- lmer(snv.rate.per.gb ~ Age + MAPD + (1|Case_ID), burden_df_dis, REML = FALSE)
      burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb - fixef(burden_age_model_ctrl)[3] * burden_df$MAPD
      } else if (plot_type == "CoV_corrected"){
        burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + CoV + (1|Case_ID), burden_df, REML = FALSE)
        burden_age_model_ctrl <- lmer(snv.rate.per.gb ~ Age + CoV + (1|Case_ID), burden_df_ctrl, REML = FALSE)
        burden_age_model_dis <- lmer(snv.rate.per.gb ~ Age + CoV + (1|Case_ID), burden_df_dis, REML = FALSE)
        burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb - fixef(burden_age_model_ctrl)[3] * burden_df$CoV
        } else if (plot_type == "Coverage_corrected"){
          burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + Coverage + (1|Case_ID), burden_df, REML = FALSE)
          burden_age_model_ctrl <- lmer(snv.rate.per.gb ~ Age + Coverage + (1|Case_ID), burden_df_ctrl, REML = FALSE)
          burden_age_model_dis <- lmer(snv.rate.per.gb ~ Age + Coverage + (1|Case_ID), burden_df_dis, REML = FALSE)
          burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb - fixef(burden_age_model_ctrl)[3] * burden_df$Coverage
          } else if (plot_type == "Depth_corrected"){
            burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + Depth + (1|Case_ID), burden_df, REML = FALSE)
            burden_age_model_ctrl <- lmer(snv.rate.per.gb ~ Age + Depth + (1|Case_ID), burden_df_ctrl, REML = FALSE)
            burden_age_model_dis <- lmer(snv.rate.per.gb ~ Age + Depth + (1|Case_ID), burden_df_dis, REML = FALSE)
            burden_df$snv.rate.per.gb_cor <- burden_df$snv.rate.per.gb - fixef(burden_age_model_ctrl)[3] * burden_df$Depth
            } else {
              print("unknown correction type")
              }
  
  burden_df_ctrl <- burden_df %>% filter(Condition == "Control")
  burden_df_dis <- burden_df %>% filter(Condition == "IHD")
  burden_age_model <- lmer(snv.rate.per.gb_cor ~ Age + Condition + (1|Case_ID), burden_df, REML = FALSE)
  # burden_age_model <- lmer(snv.rate.per.gb_cor ~ Age + Condition + Gender + (1|Case_ID), burden_df, REML = FALSE)
  burden_age_model_ctrl <- lmer(snv.rate.per.gb_cor ~ Age + (1|Case_ID), burden_df_ctrl, REML = FALSE)
  burden_age_model_dis <- lmer(snv.rate.per.gb_cor ~ Age + (1|Case_ID), burden_df_dis, REML = FALSE)
  
  # aging_rate <- formatC(signif(fixef(burden_age_model_ctrl)[2], digits = 2), digits = 2, format="fg", flag="#")
  aging_rate <- format(round(fixef(burden_age_model_ctrl)[2], digits = sig_digits), nsmall = sig_digits)
  
  ## check model parameters
  summary_lm_all <- summary(burden_age_model)
  summary_lm_ctrl <- summary(burden_age_model_ctrl)
  summary_lm_dis <- summary(burden_age_model_dis)
  
  ## compute confident interval
  CI_age_model <- confint(burden_age_model, oldNames = FALSE)
  CI_age_model_ctrl <- confint(burden_age_model_ctrl, method = "Wald", oldNames = FALSE)
  
  ## check p value
  anova_pval <- anova(burden_age_model)$"Pr(>F)"[2]
  anova_pval_print <- formatC(signif(anova_pval, digits = 2), digits = 2, format="fg", flag="#")
  
  anova_pval_ctrl <- anova(burden_age_model_ctrl)$"Pr(>F)"[1]
  anova_pval_print_ctrl <- formatC(signif(anova_pval_ctrl, digits = 2), digits = 2, format="fg", flag="#")
  
  ## manually calculate the fitting lines
  geom_line_data <- tibble(x_fit = range(burden_df_ctrl$Age), y_fit = x_fit * fixef(burden_age_model_ctrl)[2] + fixef(burden_age_model_ctrl)[1])
  
  ##### plot
  legend_data <- burden_df[c(7, 60), c("Age", "snv.rate.per.gb_cor", "Condition")]

  generate_y_breaks <- function(y_data) {
    y_min <- min(y_data, na.rm = TRUE)
    y_max <- max(y_data, na.rm = TRUE)
    if (y_min < 0) {
      min_root <- ceiling(sqrt(abs(y_min) / 100))
      neg_breaks <- -((min_root:1)^2 * 100)
    } else {
      neg_breaks <- NULL
    }
    max_root <- ceiling(sqrt(y_max / 100))
    pos_breaks <- (0:max_root)^2 * 100
    return(c(neg_breaks, pos_breaks))
  }

  yticks_plot <- generate_y_breaks(burden_df$snv.rate.per.gb_cor)
  ylim_plot <- range(yticks_plot, na.rm = TRUE)
  p_SNV_burden_lme <- ggplot(burden_df, aes(x = Age, y = snv.rate.per.gb_cor)) + 
    geom_point(pch = 21, data = legend_data, aes(x = Age, y = snv.rate.per.gb_cor, color = Condition, fill = Condition), size = 5) + 
    geom_point(pch = 21, color = "black", fill = burden_df$Color, size = 5) + 
    geom_line(aes(x = x_fit, y = y_fit), color = "dodgerblue3", data = geom_line_data) + 
    annotate("text", size = 7, x = 0, y = 0.95 * ylim_plot[2], label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB"), hjust = 0) +
    annotate("text", size = 7, x = 0, y = 0.80 * ylim_plot[2], label = paste("P = ", anova_pval_print), hjust = 0) + 
    # annotate("text", size = 6, x = 0, y = 0.65 * ylim_plot[2], label = paste0("Age accum. rate = ", format(round(fixef(burden_age_model_ctrl)[2], digits = sig_digits), nsmall = sig_digits), "/(GB.year), ", "95% CI = [", format(round(CI_age_model_ctrl[4, 1], digits = sig_digits), nsmall = sig_digits), ", ", format(round(CI_age_model_ctrl[4, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) +
    # annotate("text", size = 6, x = 0, y = 0.50 * ylim_plot[2], label = paste0("IHD excess SNVs = ", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), ", ", "95% CI = [", format(round(CI_age_model[5, 1], digits = sig_digits), nsmall = sig_digits), ", ", format(round(CI_age_model[5, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) +
    scale_color_manual(values = c("Control" = "black", "IHD" = "black"), guide = "legend") + scale_fill_manual(values = c("Control" = ctrl_dis_color[1], "IHD" = ctrl_dis_color[2]), guide = "legend") +
    scale_y_continuous(limits = ylim_plot, breaks = yticks_plot, labels = yticks_plot, trans = scales::modulus_trans(0.5)) + labs(x = "Age (years)", y = "sSNV rate per GB") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
  ggsave(paste0(main_figure_dir, "/2-SNV_burden_lme_", plot_type, ".pdf"), plot = p_SNV_burden_lme, width = 9.5, height = 6, dpi = 600)
  # ggsave(paste0(suppl_figure_dir, "/2-SNV_burden_lme_", plot_type, ".pdf"), plot = p_SNV_burden_lme, width = 9.5, height = 6, dpi = 600)
  # ggsave(paste0(other_figure_dir, "/2-SNV_burden_lme_", plot_type, ".pdf"), plot = p_SNV_burden_lme, width = 9.5, height = 6, dpi = 600)
  
  yticks_ctrl_plot <- generate_y_breaks(burden_df_ctrl$snv.rate.per.gb_cor)
  ylimit_ctrl_plot <- range(yticks_ctrl_plot, na.rm = TRUE)
  p_SNV_burden_lme_control <- ggplot(burden_df_ctrl, aes(x = Age, y = snv.rate.per.gb_cor)) + 
    geom_point(pch = 21, color = "black", fill = burden_df_ctrl$Color, size = 5) + 
    geom_line(aes(x = x_fit, y = y_fit), color = "dodgerblue3", data = geom_line_data) + 
    annotate("text", size = 7, x = 0, y = 0.95 * ylimit_ctrl_plot[2], label = paste("aging effect:", aging_rate, "sSNVs/(GB·year)"), hjust = 0) + 
    annotate("text", size = 7, x = 0, y = 0.80 * ylimit_ctrl_plot[2], label = paste("P = ", anova_pval_print_ctrl), hjust = 0) + 
    scale_y_continuous(limits = ylimit_ctrl_plot, breaks = yticks_ctrl_plot, labels = yticks_ctrl_plot, trans = scales::modulus_trans(0.5)) + labs(x = "Age (years)", y = "sSNV rate per GB") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
  ggsave(paste0(main_figure_dir, "/2-SNV_burden_lme_", plot_type, "_control.pdf"), plot = p_SNV_burden_lme_control, width = 8, height = 6, dpi = 600)
  # ggsave(paste0(other_figure_dir, "/2-SNV_burden_lme_", plot_type, "_control.pdf"), plot = p_SNV_burden_lme_control, width = 8, height = 6, dpi = 600)
  }