#############################################################################
##################### sSNV burden mixed-effects modeling ####################
#############################################################################
# sum(sSNV_SCAN2_df$snv.burden[ctrl_range])
# sum(sSNV_SCAN2_df$snv.burden[dis_range])
# sum(sSNV_SCAN2_df$snv.burden)
y_break_range <- c(4000, 6800)

##### input a matrix with columns: burden.per.gb, Case_ID, MAPD, CoV, Age, Condition
burden_df <- sSNV_SCAN2_df[, c("snv.rate.per.gb", "Case_ID", "MAPD", "CoV", "Coverage", "Depth", "Age", "Condition", "Gender", "PMI", "Color", "Outline_Color")] %>% 
  mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% 
  mutate(Condition = relevel(Condition, ref = ctrl_name)) %>% 
  group_by(Condition) %>% arrange(Age, .by_group = TRUE)
burden_df_ctrl <- burden_df %>% filter(Condition == "Control")
burden_df_dis <- burden_df %>% filter(Condition == "IHD")

#############################################################################
##### mixed-effects modeling
sig_digits <- 2
CI_P_list <- list()

lme_corr_switch_func <- function(burden_df, burden_df_ctrl, plot_type = "Uncorr") {
  covariate_map <- c("MAPD_corr" = "MAPD", "CoV_corr" = "CoV", "Coverage_corr" = "Coverage", "Depth_corr" = "Depth", "PMI_corr" = "PMI")
  if (plot_type == "Uncorr") {
    formula_ctrl_str <- "snv.rate.per.gb ~ Age + (1|Case_ID)"
  } else if (plot_type %in% names(covariate_map)) {
    covariate <- covariate_map[[plot_type]]
    formula_ctrl_str <- paste0("snv.rate.per.gb ~ Age + ", covariate, " + (1|Case_ID)")
  } else {
    stop("Unknown plot type")
  }
  
  burden_age_model_ctrl <- lmer(as.formula(formula_ctrl_str), burden_df_ctrl, REML = FALSE)
  
  if (plot_type == "Uncorr") {
    burden_df$snv.rate.per.gb_corr <- burden_df$snv.rate.per.gb
  } else {
    burden_df$snv.rate.per.gb_corr <- burden_df$snv.rate.per.gb - fixef(burden_age_model_ctrl)[3] * burden_df[[covariate]]
  }
  
  return(list(df = burden_df))
}

# for (plot_type in c("Uncorr")){
for (plot_type in c("Uncorr", "MAPD_corr", "CoV_corr", "Coverage_corr", "Depth_corr", "PMI_corr")){
  cat("lmer modeling:", plot_type, "\n")
  burden_df_corr <- lme_corr_switch_func(burden_df, burden_df_ctrl, plot_type = plot_type)$df
  
  burden_df_ctrl <- burden_df_corr %>% filter(Condition == "Control")
  burden_age_model <- lmer(snv.rate.per.gb_corr ~ Age + Condition + (1|Case_ID), burden_df_corr, REML = FALSE)
  # burden_age_model <- lmer(snv.rate.per.gb_corr ~ Age + Condition + Gender + (1|Case_ID), burden_df_corr, REML = FALSE)
  burden_age_model_ctrl <- lmer(snv.rate.per.gb_corr ~ Age + (1|Case_ID), burden_df_ctrl, REML = FALSE)
  
  aging_rate <- format(round(fixef(burden_age_model_ctrl)[2], digits = sig_digits), nsmall = sig_digits)
  
  ## check model parameters
  summary_lm_all <- summary(burden_age_model)
  summary_lm_ctrl <- summary(burden_age_model_ctrl)
  
  ## compute confident interval
  CI_age_model <- confint(burden_age_model, method = "Wald", oldNames = FALSE)
  CI_age_model_ctrl <- confint(burden_age_model_ctrl, method = "Wald", oldNames = FALSE)
  
  ## check p value
  anova_pval <- anova(burden_age_model)$"Pr(>F)"[2]
  anova_pval_print <- formatC(signif(anova_pval, digits = 2), digits = 2, format="fg", flag="#")
  anova_pval_ctrl <- anova(burden_age_model_ctrl)$"Pr(>F)"[1]
  anova_pval_print_ctrl <- formatC(signif(anova_pval_ctrl, digits = 2), digits = 2, format="e", flag="#")
  
  ## save to csv
  Aging_effect <- paste0(plot_type, " Aging effect: ", format(round(fixef(burden_age_model_ctrl)[2], digits = sig_digits), nsmall = sig_digits), 
                         " sSNVs/(GB.year), 95% CI = [", format(round(CI_age_model_ctrl[4, 1], digits = sig_digits), nsmall = sig_digits), ", ", 
                         format(round(CI_age_model_ctrl[4, 2], digits = sig_digits), nsmall = sig_digits), "]", ", P = ", anova_pval_print_ctrl)
  IHD_effect <- paste0(plot_type, " IHD effect: ", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), 
                       " sSNVs/GB, 95% CI = [", format(round(CI_age_model[5, 1], digits = sig_digits), nsmall = sig_digits), ", ", 
                       format(round(CI_age_model[5, 2], digits = sig_digits), nsmall = sig_digits), "]", ", P = ", anova_pval_print)
  CI_P_list[[plot_type]] <- data.frame(Correction_type = plot_type, Measure = c("Age accumulation rate", "IHD excess sSNVs"), Description = c(Aging_effect, IHD_effect))
  
  ## manually calculate the fitting lines
  geom_line_data <- tibble(x_fit = range(burden_df_ctrl$Age), y_fit = x_fit * fixef(burden_age_model_ctrl)[2] + fixef(burden_age_model_ctrl)[1])
  legend_data <- burden_df_corr[c(1, nrow(burden_df_corr)), c("Age", "snv.rate.per.gb_corr", "Condition")]
  
  generate_y_breaks <- function(y_data) {
    y_min <- min(y_data, na.rm = TRUE)
    y_max <- max(y_data, na.rm = TRUE)
    start_val <- if (y_min < 0) y_min else 0
    range_span <- y_max - start_val
    step_size <- signif(range_span / 6, 1)
    breaks <- seq(0, ceiling(y_max / step_size) * step_size, by = step_size)
    return(breaks)
  }
  
  yticks_plot <- generate_y_breaks(c(burden_df_corr$snv.rate.per.gb_corr, outliers_df$snv.rate.per.gb, 8000))
  ylim_plot <- range(c(min(burden_df_corr$snv.rate.per.gb_corr), yticks_plot), na.rm = TRUE)
  p_sSNV_lme <- ggplot(burden_df_corr, aes(x = Age, y = snv.rate.per.gb_corr)) + 
    geom_point(pch = 21, data = legend_data, aes(x = Age, y = snv.rate.per.gb_corr, color = Condition, fill = Condition), size = 5) + 
    geom_text(data = outliers_df, aes(x = Age, y = snv.rate.per.gb), label = "×", size = 8, color = "black") + 
    geom_point(pch = 21, color = burden_df_corr$Outline_Color, fill = burden_df_corr$Color, size = 5) + 
    geom_line(aes(x = x_fit, y = y_fit), color = "dodgerblue2", linetype = "22", data = geom_line_data) + 
    geom_hline(yintercept = 0, color = "gray40", linetype = "22") + 
    annotate("text", size = 6, x = 0, y = 0.93 * ylim_plot[2], 
             label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB", ", P = ", anova_pval_print), hjust = 0) +
    annotate("text", size = 6, x = 0, y = 0.98 * ylim_plot[2], 
             label = paste("Aging effect:", aging_rate, "sSNVs/(GB·year)", ", P = ", anova_pval_print_ctrl), hjust = 0) + 
    scale_color_manual(values = c("Control" = "black", "IHD" = "black"), guide = "legend") + 
    scale_fill_manual(values = c("Control" = ctrl_dis_color[1], "IHD" = ctrl_dis_color[2]), guide = "legend") + 
    scale_y_break(y_break_range, scale = 0.25) +
    scale_y_continuous(limits = ylim_plot, breaks = yticks_plot, labels = yticks_plot) + labs(x = "Age (years)", y = "sSNV rate per GB") + theme_linedraw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.40, vjust = 0), 
          axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank())
  ggsave(paste0(suppl_figure_dir, "/2-sSNV_lme_", plot_type, ".pdf"), plot = p_sSNV_lme, width = 9, height = 6, dpi = 600)
  
  yticks_ctrl_plot <- generate_y_breaks(burden_df_ctrl$snv.rate.per.gb_corr)
  ylimit_ctrl_plot <- range(yticks_ctrl_plot, na.rm = TRUE)
  p_sSNV_lme_ctrl <- ggplot(burden_df_ctrl, aes(x = Age, y = snv.rate.per.gb_corr)) + 
    geom_point(pch = 21, color = burden_df_ctrl$Outline_Color, fill = burden_df_ctrl$Color, size = 5) + 
    # geom_text(data = outliers_df, aes(x = Age, y = snv.rate.per.gb), label = "×", size = 8, color = "black") +
    geom_line(aes(x = x_fit, y = y_fit), color = "dodgerblue2", linetype = "22", data = geom_line_data) + 
    geom_hline(yintercept = 0, color = "gray40", linetype = "22") + 
    annotate("text", size = 6, x = 0, y = 0.98 * ylimit_ctrl_plot[2], 
             label = paste("Aging effect:", aging_rate, "sSNVs/(GB·year)"), hjust = 0) + 
    annotate("text", size = 6, x = 0, y = 0.88 * ylimit_ctrl_plot[2], 
             label = paste("P = ", anova_pval_print_ctrl), hjust = 0) + 
    scale_y_continuous(limits = ylimit_ctrl_plot, breaks = yticks_ctrl_plot, labels = yticks_ctrl_plot) + 
    labs(x = "Age (years)", y = "sSNV rate per GB") + theme_linedraw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), 
          text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
  if (plot_type == "Uncorr"){
    ggsave(paste0(main_figure_dir, "/2-sSNV_lme_", plot_type, ".pdf"), plot = p_sSNV_lme, width = 9, height = 6, dpi = 600)
    ggsave(paste0(main_figure_dir, "/2-sSNV_lme_", plot_type, "_ctrl.pdf"), plot = p_sSNV_lme_ctrl, width = 7, height = 6, dpi = 600)
  }
}

CI_P_df <- do.call(rbind, CI_P_list)
write.csv(CI_P_df, paste0(table_dir, "/2-sSNV_lme_modeling.csv"), row.names = FALSE)

## lme test for smoking categories
burden_df_dis_02 <- burden_df_dis %>% 
  mutate(log10snv.rate.per.gb = log10(snv.rate.per.gb), 
         Smoking_category = c(rep("Smoker", 5), rep("Smoker", 8), rep("Smoker", 8), rep("Smoker", 4), rep("Hevay Smoker", 8)), 
         Smoking_category_num = c(rep(0, 25), rep(1, 8)))
fit_lmm <- lmer(snv.rate.per.gb ~ Age + Smoking_category_num + (1 | Case_ID), data = burden_df_dis_02)
summary(fit_lmm)

## lme test for smoking categories 1363, 1673, 604, 1743, 1113
burden_df_dis_smoking <- burden_df_dis %>% 
  mutate(log10snv.rate.per.gb = log10(snv.rate.per.gb), 
         Smoking_category = c(rep("Smoker", 5), rep("Smoker", 8), rep("Smoker", 8), rep("Smoker", 4), rep("Hevay Smoker", 8)), 
         Smoking_category_num = c(rep(0, 25), rep(1, 8)))
fit_lmm_smoking <- lmer(snv.rate.per.gb ~ Age + Smoking_category_num + (1 | Case_ID), data = burden_df_dis_smoking)
summary(fit_lmm_smoking)

## lme test for calcification 1363, 1673, 604, 1743, 1113
# left anterior:  60%, 100%, 80%, 75%, 90%
# right anterior: 80%, 75%,  20%, 40%, 60%

burden_df_dis_left_ant <- burden_df_dis %>% 
  mutate(log10snv.rate.per.gb = log10(snv.rate.per.gb), 
         left_cal = c(rep(0.6, 5), rep(1.0, 8), rep(0.8, 8), rep(0.75, 4), rep(0.9, 8)))
fit_lmm_left_ant <- lmer(snv.rate.per.gb ~ Age + left_cal + (1 | Case_ID), data = burden_df_dis_left_ant)
summary(fit_lmm_left_ant)

burden_df_dis_right_ant <- burden_df_dis %>% 
  mutate(log10snv.rate.per.gb = log10(snv.rate.per.gb), 
         right_cal = c(rep(0.8, 5), rep(0.75, 8), rep(0.2, 8), rep(0.4, 4), rep(0.6, 8)))
fit_lmm_right_ant <- lmer(snv.rate.per.gb ~ Age + right_cal + (1 | Case_ID), data = burden_df_dis_right_ant)
summary(fit_lmm_right_ant)
