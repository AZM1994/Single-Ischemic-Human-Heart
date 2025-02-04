###############################################################################
##### 1. determine the highest contribution signatures            #############
##### 2. find signatures that are significantly different between #############
#####    the normal and disease conditions                        #############
###############################################################################

##### determine the highest contribution signatures
sig_contri <- as.data.frame(t(SigNet_contri))
sig_contri_age_match <- sig_contri[Cell_ID_list_age_match, ]
contri_for_each_sig <- colSums(sig_contri)
top_contri_for_each_sig <- sort(contri_for_each_sig, decreasing = TRUE)[1:20]
top_sigs <- row.names(data.frame(top_contri_for_each_sig))

##### identify point mutation type in top signatures
## plot the top signature spectrum: COSMIC_v3.4_sigs, COSMIC_v3.3.1_sigs
top_sigs_spectrum <- COSMIC_v3.4_sigs %>% .[, top_sigs] %>% {rownames(.) <- rownames(mut_mat_age_match_conditional); .}
p_top_sigs_spectrum <- plot_96_profile_abs(top_sigs_spectrum)
# ggsave(paste0(other_figure_dir, "/p_top_sigs_spectrum.pdf"), plot = p_top_sigs_spectrum, width = 20, height = 15, dpi = 300)

## plot the top signature spectrum by Condition
top_sigs_spectrum_by_condition <- as.data.frame(top_sigs_spectrum) %>% 
  mutate(Control = top_sigs_spectrum %*% colMeans(sig_contri[control_range, top_sigs])) %>% 
  mutate(IHD = top_sigs_spectrum %*% colMeans(sig_contri[disease_range, top_sigs])) %>% 
  mutate(Net_change = IHD - Control) %>% dplyr::select(control_name, disease_name, "Net_change")
p_top_sigs_spectrum_by_condition <- plot_96_profile_abs(top_sigs_spectrum_by_condition[ , c(control_name, disease_name, "Net_change")])
# ggsave(paste0(other_figure_dir, "/p_top_sigs_spectrum_by_condition.pdf"), plot = p_top_sigs_spectrum_by_condition, width = 20, height = 10, dpi = 300)

## find differences between top signature spectrum and ground truth
# ground_truth <- mut_mat_raw_cond[c(control_name, disease_name), ] %>% t() %>% {1 / colSums(.) * .}
# diff_spectrum <- top_sigs_spectrum_by_condition - ground_truth
# p_diff_spectrum <- plot_96_profile_abs(diff_spectrum)
# ggsave(paste0(other_figure_dir, "/p_diff_spectrum.pdf"), plot = p_diff_spectrum, width = 20, height = 10, dpi = 300)

##### top signatures contribution
Reshape <- function(a, n, m){if (missing(m)) m <- length(a)%/%n
  if (length(a) != n * m) stop("Matrix 'a' does not have n*m elements")
  dim(a) <- c(n, m)
  return(a)}

## raw SigNet results
top_sig_contri_cond <- sig_contri %>% .[, colnames(.) %in% top_sigs] %>% 
  mutate(others = rowSums(sig_contri[, !(colnames(sig_contri) %in% top_sigs)])) %>% t() %>% as.data.frame() %>% 
  mutate(Control = rowMeans(.[, control_range_age_match]), IHD = rowMeans(.[, disease_range_age_match])) %>% dplyr::select(control_name, disease_name)
sorted_labels <- as.data.frame(row.names(top_sig_contri_cond)) %>% mutate(A = .[[1]]) %>% t() %>% as.matrix() %>% Reshape(nrow(top_sig_contri_cond) * 2, 1)
dim(sorted_labels) <- c(nrow(top_sig_contri_cond) * 2,1)

label_size <- 6 * rbind(top_sig_contri_cond[,1] / top_sig_contri_cond[,1], top_sig_contri_cond[,2] / top_sig_contri_cond[,2])
label_size <- Reshape(label_size,ncol(label_size)*2,1)
label_size[is.na(label_size)] <- 0

##### plot raw contribution for Control and IHD: normalized to 1
p_sSNV_top_sigs_contri_raw <- plot_contribution(top_sig_contri_cond, coord_flip = FALSE, mode = "relative") + aes(label = sorted_labels) + 
  geom_text(size = label_size, position = position_stack(vjust = 0.5), col = "white", fontface = "bold") + 
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), 
                           panel.grid.minor = element_blank(), panel.border = element_rect(size = 0.5), text = element_text(size = 12))
ggsave(paste0(other_figure_dir, "/4-sSNV_top_sigs_contri_raw.pdf"), plot = p_sSNV_top_sigs_contri_raw, width = 12, height = 18, dpi = 300)

##########################################################################
##### normalize total contribution for each cell to estimated burden
sig_contri_normalized <- 1 / rowSums(sig_contri) * sig_contri * SCAN2_df$snv.rate.per.gb
sig_contri_normalized_age_match <- sig_contri_normalized[Cell_ID_list_age_match, ]
contri_for_each_sig <- colSums(sig_contri_normalized)
top_contri_for_each_sig <- sort(contri_for_each_sig, decreasing = TRUE)[1:20]
top_sigs <- row.names(data.frame(top_contri_for_each_sig))

##### select signatures that are above certain percentage in either Condition
sig_percentage_thres <- 10
selected_sigs <- sig_contri_normalized[control_range, ] %>% colSums() %>% as.data.frame() %>% setNames(c("Control_sig_contri")) %>% mutate(Control_Pct = Control_sig_contri / sum(Control_sig_contri) * 100) %>% 
  bind_cols(sig_contri_normalized[disease_range, ] %>% colSums() %>% as.data.frame() %>% setNames(c("IHD_sig_contri")) %>% mutate(IHD_Pct = IHD_sig_contri / sum(IHD_sig_contri) * 100)) %>% 
  filter(Control_Pct >= sig_percentage_thres | IHD_Pct >= sig_percentage_thres)

sig_contri_normalized_age_match_cond <- sig_contri_normalized_age_match %>% .[, colnames(.) %in% top_sigs] %>% 
  mutate(others = rowSums(.[, !(colnames(.) %in% top_sigs)])) %>% t() %>% as.data.frame() %>% 
  mutate(Control = rowMeans(.[, control_range_age_match]), IHD = rowMeans(.[, disease_range_age_match])) %>% dplyr::select(control_name, disease_name) %>% 
  mutate(Control_Pct = Control / sum(Control) * 100, IHD_Pct = IHD / sum(IHD) * 100, 
         Control_label = paste0(rownames(.), " (", round(Control_Pct, 1), "%)"), IHD_label = paste0(rownames(.), " (", round(IHD_Pct, 1), "%)"), 
         Control_label = ifelse(Control_Pct < 0.7, NA, Control_label), IHD_label = ifelse(IHD_Pct < 0.7, NA, IHD_label))

sorted_labels <- as.data.frame(row.names(sig_contri_normalized_age_match_cond)) %>%
  mutate(A = .[[1]]) %>% t() %>% as.matrix() %>% Reshape(nrow(sig_contri_normalized_age_match_cond) * 2,1)

label_size <- 6 * rbind(sig_contri_normalized_age_match_cond[,1] / sig_contri_normalized_age_match_cond[,1], sig_contri_normalized_age_match_cond[,2] / sig_contri_normalized_age_match_cond[,2])
label_size <- Reshape(label_size,ncol(label_size)*2,1)
label_size[is.na(label_size)] <- 0

p_sSNV_top_sigs_burden <- plot_contribution(sig_contri_normalized_age_match_cond[,1:2], coord_flip = F, mode = "absolute") +
  # aes(label = sorted_labels) + geom_text(size = label_size, position = position_stack(vjust = 0.5), col = "white", fontface = "bold") +
  theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), 
                           panel.grid.minor = element_blank(), panel.border = element_rect(size = 0.5), text = element_text(size = 12))
ggsave(paste0(main_figure_dir, "/4-sSNV_top_sigs_burden.pdf"), plot = p_sSNV_top_sigs_burden, width = 4, height = 6, dpi = 600)

##############################################################################
### add metadata and remove the signatures with all zeros in one Condition
sig_contri_normalized <- sig_contri_normalized %>% 
  mutate(Age = metadata_df$Age, Condition = relevel(factor(metadata_df$Condition), ref = "Control"), 
         Gender = metadata_df$Gender, Case_ID = metadata_df$Case_ID, Color = SCAN2_df$Color) %>% 
  filter(row.names(.) %in% c(filtered_mut_mat_est_Control, filtered_mut_mat_est_IHD)) %>% {
  sig_contri_normalized_control <- .[.$Condition == control_name, ]
  zero_col_control <- which(colSums(sig_contri_normalized_control[, 1:(ncol(sig_contri_normalized_control) - 5)]) == 0)
  . <- .[, -zero_col_control]
  sig_contri_normalized_disease <- .[.$Condition == disease_name, ]
  zero_col_disease <- which(colSums(sig_contri_normalized_disease[, 1:(ncol(sig_contri_normalized_disease) - 5)]) == 0)
  .[, -zero_col_disease]}

refined_selected_sigs_list <- colnames(sig_contri_normalized)[1 : (ncol(sig_contri_normalized) - 5)]
sig_digits <- 2

##### Mix effects model & Age matched mix effects model
pdf(paste0(other_figure_dir, "/4-lmer_top_sigs_all.pdf"), width = 12, height = 8)
  # pdf(paste0(other_figure_dir, "/4-lmer_top_sigs_reported.pdf"), width = 12, height = 8)
  # for (i in c(1,4,5,9,21,22,29,34,36)){
  # for (i in c(1)){
  for (i in c(1 : (ncol(sig_contri_normalized) - 5))){
    cat("index", i, ":", refined_selected_sigs_list[i])
    ## burden_df: the input matrix with columns: burden, Condition, Case_ID, Age, Color
    burden_df <- sig_contri_normalized[c(refined_selected_sigs_list[i], "Age", "Condition", "Case_ID", "Color")] %>% setNames(c("snv.rate.per.gb", "Age", "Condition", "Case_ID", "Color")) %>% 
      mutate(Condition = relevel(factor(Condition), ref = "Control")) %>% group_by(factor(Condition, levels = c(disease_name, control_name))) %>% arrange(Age, .by_group = TRUE)
    burden_df_control <- burden_df %>% filter(Condition == "Control")
    burden_df_disease <- burden_df %>% filter(Condition == "IHD")
    
    if (sum(burden_df_control$snv.rate.per.gb) > 6 & sum(burden_df_disease$snv.rate.per.gb) > 6){
      ## mixed effects modeling
      burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + (1|Case_ID), burden_df, REML = FALSE)
      burden_age_model_control <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_control, REML = FALSE)
      # burden_age_model_disease <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_disease, REML = FALSE)

      ## compute confident interval
      # CI_age_model <- suppressMessages(confint(burden_age_model, oldNames = FALSE))
      # CI_age_model_control <- suppressMessages(confint(burden_age_model_control, oldNames = FALSE))
      
      ## check p value and R^2
      anova_pvalue <- anova(burden_age_model)$"Pr(>F)"[2]
      anova_pvalue_print <- formatC(signif(anova_pvalue, digits = sig_digits), digits = sig_digits, format="fg", flag="#")
      anova_pvalue_control <- anova(burden_age_model_control)$"Pr(>F)"[1]
      anova_pvalue_print_control <- formatC(signif(anova_pvalue_control, digits = sig_digits), digits = sig_digits, format="fg", flag="#")
      
      ## manually calculate the fitting lines
      geom_line_data <- burden_df_control %>% {
        Control_fitting_x <- range(.$Age)
        Control_fitting_y <- Control_fitting_x * fixef(burden_age_model_control)[2] + fixef(burden_age_model_control)[1]
        geom_line_data <- rbind(Control_fitting_x, Control_fitting_y)
        as.data.frame(t(geom_line_data))}
      
      aging_rate <- format(round(fixef(burden_age_model_control)[2], digits = sig_digits), nsmall = sig_digits)
      legend_data <- burden_df[c(7, 60), c("Age", "snv.rate.per.gb", "Condition")]
      max_y <- max(burden_df$snv.rate.per.gb)
      
      p_SNV_burden_lme <- ggplot(burden_df, aes(x = Age, y = snv.rate.per.gb)) + 
        geom_point(pch = 21, data = legend_data, aes(x = Age, y = snv.rate.per.gb, color = Condition, fill = Condition), size = 5) + 
        geom_point(pch = 21, color = "black", fill = burden_df$Color, size = 5) + 
        geom_line(aes(x = Control_fitting_x, y = Control_fitting_y), colour = "dodgerblue3", data = geom_line_data) + 
        annotate("text", size = 6, x = 0, y = 0.95 * max_y, label = paste("aging effect:", aging_rate, "sSNVs/(GB·year),", "P =", anova_pvalue_print_control), hjust = 0) +
        annotate("text", size = 6, x = 0, y = 0.87 * max_y, label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB,", "P =", anova_pvalue_print), hjust = 0) +
        # annotate("text", size = 6, x = 0, y = 0.87 * max_y, label = paste0("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "/GB, ", 
        #                                                                    "95% CI = [", format(round(CI_age_model[5, 1], digits = sig_digits), nsmall = sig_digits), ", ", 
        #                                                                    format(round(CI_age_model[5, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) + 
        # annotate("text", size = 6, x = 0, y = 0.95 * max_y, label = paste0("aging effect:", format(round(fixef(burden_age_model_control)[2], digits = sig_digits), nsmall = sig_digits), "/(GB.year), ", 
        #                                                                    "95% CI = [", format(round(CI_age_model_control[4, 1], digits = sig_digits), nsmall = sig_digits), ", ", 
        #                                                                    format(round(CI_age_model_control[4, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) +
        scale_color_manual(values = c("IHD" = "black", "Control" = "black"), guide = "legend") + scale_fill_manual(values = c("Control" = ctrl_dis_color[1], "IHD" = ctrl_dis_color[2]), guide = "legend") + 
        labs(x = "Age (years)", y = paste0(refined_selected_sigs_list[i], " Contribution \n (sSNV rate per GB)")) + theme_linedraw() + 
        theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
              panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0))
      print(p_SNV_burden_lme)
      if (anova_pvalue <= 0.05){
        ggsave(paste0(main_figure_dir, "/4-", refined_selected_sigs_list[i], "_SNV_burden_lme.pdf"), plot = p_SNV_burden_lme, width = 9, height = 6, dpi = 600)
      }
      
    } else {
      print("low_num")
      ## mixed linear modeling
      burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + (1|Case_ID), burden_df, REML = FALSE)
      burden_age_model_control <- lm(snv.rate.per.gb ~ Age, burden_df_control)
      
      ## compute confident interval
      # CI_age_model <- confint(burden_age_model, oldNames = FALSE)
      # CI_age_model_control <- confint(burden_age_model_control, oldNames = FALSE)
      
      ## check p value and R^2
      anova_pvalue <- anova(burden_age_model)$"Pr(>F)"[2]
      anova_pvalue_print <- formatC(signif(anova_pvalue, digits = 2), digits = 2, format="fg", flag="#")
      anova_pvalue_control <- anova(burden_age_model_control)$"Pr(>F)"[1]
      anova_pvalue_print_control <- formatC(signif(anova_pvalue_control, digits = 2), digits = 2, format="fg", flag="#")
      
      ## manually calculate the fitting lines
      geom_line_data <- burden_df_control %>% {
        Control_fitting_x <- range(.$Age)
        Control_fitting_y <- Control_fitting_x * coef(burden_age_model_control)[2] + coef(burden_age_model_control)[1]
        geom_line_data <- rbind(Control_fitting_x, Control_fitting_y)
        as.data.frame(t(geom_line_data))}
      
      aging_rate <- format(round(coef(burden_age_model_control)[2], digits = sig_digits), nsmall = sig_digits)
      legend_data <- burden_df[c(7, 60), c("Age", "snv.rate.per.gb", "Condition")]
      max_y <- max(burden_df$snv.rate.per.gb)
      
      p_SNV_burden_lme <- ggplot(burden_df, aes(x = Age, y = snv.rate.per.gb)) + 
        geom_point(pch = 21, data = legend_data, aes(x = Age, y = snv.rate.per.gb, color = Condition, fill = Condition), size = 5) + 
        geom_point(pch = 21, color = "black", fill = burden_df$Color, size = 5) +
        geom_line(aes(x = Control_fitting_x, y = Control_fitting_y), colour = "dodgerblue3", data = geom_line_data) + 
        annotate("text", size = 6, x = 0, y = 0.95 * max_y, label = paste("aging effect:", aging_rate, "sSNVs/(GB·year),", "P =", anova_pvalue_print_control), hjust = 0) +
        annotate("text", size = 6, x = 0, y = 0.87 * max_y, label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB,", "P =", anova_pvalue_print), hjust = 0) +
        # annotate("text", size = 6, x = 0, y = 0.87 * max_y, label = paste0("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "/GB, ", 
        #                                                                    "95% CI = [", format(round(CI_age_model[5, 1], digits = sig_digits), nsmall = sig_digits), ", ", 
        #                                                                    format(round(CI_age_model[5, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) + 
        # annotate("text", size = 6, x = 0, y = 0.95 * max_y, label = paste0("aging effect:", format(round(fixef(burden_age_model_control)[2], digits = sig_digits), nsmall = sig_digits), "/(GB.year), ", 
        #                                                                    "95% CI = [", format(round(CI_age_model_control[4, 1], digits = sig_digits), nsmall = sig_digits), ", ", 
        #                                                                    format(round(CI_age_model_control[4, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) + 
        scale_color_manual(values = c("IHD" = "black", "Control" = "black"), guide = "legend") + scale_fill_manual(values = c("Control" = ctrl_dis_color[1], "IHD" = ctrl_dis_color[2]), guide = "legend") + 
        labs(x = "Age (years)", y = paste0(refined_selected_sigs_list[i], " Contribution \n (sSNV rate per GB)")) + theme_linedraw() + 
        theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
              panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0))
      print(p_SNV_burden_lme)
      if (anova_pvalue <= 0.05){
        ggsave(paste0(main_figure_dir, "/4-", refined_selected_sigs_list[i], "_SNV_burden_lme.pdf"), plot = p_SNV_burden_lme, width = 9, height = 6, dpi = 600)
      }
    }
  }
dev.off()
