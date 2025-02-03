##### Mutational signatures extraction
# mut_mat_age_match <- mut_mat_age_match + 0.0001
# mut_mat_est_denovo <- t(mut_mat_est) + 0.0001
mut_mat_est_denovo <- t(mut_mat_est_filtered) + 0.0001
# signatures = get_known_signatures()
COSMIC_v3.4_sigs <- read.table("main/COSMIC/COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE) %>% {rownames(.) <- .$Type; .} %>% select(-Type) %>% as.matrix()

#################### De Novo signature extraction with NMF
# estimate <- nmf(mut_mat_est_denovo, rank = 2:10, method = "brunet", nrun = 30, seed = 123456, .opt = "v-p")
# print(plot(estimate))

################################################################################
############## choose rank 4
nmf_res_4 <- extract_signatures(mut_mat_est_denovo, rank = 4, nrun = 30, single_core = TRUE)
nmf_res_4_contri <- nmf_res_4$contribution
nmf_res_4_sigs <- nmf_res_4$signatures
nmf_res_4_sigs <- nmf_res_4_sigs[, c(3, 2, 4, 1)]
colnames(nmf_res_4_sigs) <- c("Signature N4", "Signature N3", "Signature N2", "Signature N1")
rownames(nmf_res_4_contri) <- c("Signature N4", "Signature N3", "Signature N2", "Signature N1")
nmf_res_4_sigs <- apply(nmf_res_4_sigs, 2, function(x) x / sum(x))
write.csv(t(nmf_res_4_sigs), file = paste0(table_dir, "/denovo_sigs_rank4.csv"))

pdf(paste0(other_figure_dir, "/4-Visualizing_nmf_res_rank4.pdf"), width = 14, height = 10)
  p_96_denovo_rank4 <- plot_96_profile(nmf_res_4_sigs) + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.15), panel.grid.minor = element_blank(), 
          panel.border = element_rect(size = 0.5), text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(paste0(suppl_figure_dir, "/4-denovo_sigs_spectrum_rank4.pdf"), plot = p_96_denovo_rank4, width = 12, height = 6, dpi = 600)
  
  nmf_res_4_contri_cond <- tibble::tibble(
    SigN4_control = sum(nmf_res_4_contri[1, control_range_age_match]), SigN4_disease = sum(nmf_res_4_contri[1, disease_range_age_match]), 
    SigN3_control = sum(nmf_res_4_contri[2, control_range_age_match]), SigN3_disease = sum(nmf_res_4_contri[2, disease_range_age_match]),
    SigN2_control = sum(nmf_res_4_contri[3, control_range_age_match]), SigN2_disease = sum(nmf_res_4_contri[3, disease_range_age_match]), 
    SigN1_control = sum(nmf_res_4_contri[4, control_range_age_match]), SigN1_disease = sum(nmf_res_4_contri[4, disease_range_age_match])) %>%
    t() %>% matrix(nrow = 4, ncol = 2, byrow = TRUE) %>% `rownames<-`(c("SigN4", "SigN3", "SigN2", "SigN1")) %>% `colnames<-`(c("Control", "IHD"))
  
  plot_contribution(nmf_res_4_contri_cond, nmf_res_4_sigs, mode = "relative") + theme(text = element_text(size=20))
  plot_contribution_heatmap(nmf_res_4_contri_cond, cluster_samples = FALSE, cluster_sigs = FALSE) + theme(text = element_text(size = 20))
  plot_contribution(nmf_res_4_contri, nmf_res_4_sigs, mode = "relative") + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0))
  plot_contribution_heatmap(nmf_res_4_contri, cluster_samples = FALSE, cluster_sigs = FALSE) + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5)) + geom_tile(colour = "white")

  cos_sim_samples_signatures_4 <- cos_sim_matrix(COSMIC_v3.4_sigs, nmf_res_4_sigs)
  plot_cosine_heatmap(t(cos_sim_samples_signatures_4), cluster_rows = F, cluster_cols = FALSE, method = "complete", plot_values = TRUE) + theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.0)) + geom_tile(colour = "white")
dev.off()

###### signature decomposition from Signet
### check the denovo sigs similar to which COSMIC sigs (use 61 COSMIC in SigNet)
SigNet_denovo_contri_rank4 <- read.csv(paste0(project_dir, "/data/SigNet/PTA_all_denovo_rank4/weight_guesses.csv"), row.names = 1) %>% t()
p_signet_denovo_rank4 <- plot_contribution_heatmap(SigNet_denovo_contri_rank4, cluster_sigs = FALSE, cluster_samples = FALSE, method = "complete", plot_values = FALSE) + scale_fill_gradient(low = "white", high = "dodgerblue3") + geom_tile(colour = "black") + theme_linedraw() + ggtitle(NULL) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.0))
ggsave(paste0(suppl_figure_dir, "/4-denovo_signet_contri_rank4.pdf"), plot = p_signet_denovo_rank4, width = 15, height = 5, dpi = 600)

################################################################################
################################################################################
############## choose rank 5
nmf_res_5 <- extract_signatures(mut_mat_est_denovo, rank = 5, nrun = 30, single_core = TRUE)
nmf_res_5_contri <- nmf_res_5$contribution
nmf_res_5_sigs <- nmf_res_5$signatures
nmf_res_5_sigs <- nmf_res_5_sigs[, c(2, 3, 4, 5, 1)] 
colnames(nmf_res_5_sigs) <- c("Signature N5", "Signature N4", "Signature N3", "Signature N2", "Signature N1")
rownames(nmf_res_5_contri) <- c("Signature N5", "Signature N4", "Signature N3", "Signature N2", "Signature N1")
nmf_res_5_sigs <- apply(nmf_res_5_sigs, 2, function(x) x / sum(x))
write.csv(t(nmf_res_5_sigs), file = paste0(table_dir, "/denovo_sigs_rank5.csv"))

pdf(paste0(other_figure_dir, "/4-Visualizing_nmf_res_rank5.pdf"), width = 14, height = 10)
  p_96_denovo_rank5 <- plot_96_profile(nmf_res_5_sigs) + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.15), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(paste0(suppl_figure_dir, "/4-denovo_sigs_spectrum_rank5.pdf"), plot = p_96_denovo_rank5, width = 12, height = 6, dpi = 600)

  nmf_res_5_contri_cond <- tibble::tibble(
    SigN5_control = sum(nmf_res_5_contri[1, control_range_age_match]), SigN5_disease = sum(nmf_res_5_contri[1, disease_range_age_match]), 
    SigN4_control = sum(nmf_res_5_contri[2, control_range_age_match]), SigN4_disease = sum(nmf_res_5_contri[2, disease_range_age_match]), 
    SigN3_control = sum(nmf_res_5_contri[3, control_range_age_match]), SigN3_disease = sum(nmf_res_5_contri[3, disease_range_age_match]),
    SigN2_control = sum(nmf_res_5_contri[4, control_range_age_match]), SigN2_disease = sum(nmf_res_5_contri[4, disease_range_age_match]), 
    SigN1_control = sum(nmf_res_5_contri[5, control_range_age_match]), SigN1_disease = sum(nmf_res_5_contri[5, disease_range_age_match])) %>%
    t() %>% matrix(nrow = 5, ncol = 2, byrow = TRUE) %>% `rownames<-`(c("SigN5", "SigN4", "SigN3", "SigN2", "SigN1")) %>% `colnames<-`(c("Control", "IHD"))

  plot_contribution(nmf_res_5_contri_cond, nmf_res_5_sigs, mode = "relative") + theme(text = element_text(size=20))
  plot_contribution_heatmap(nmf_res_5_contri_cond, cluster_samples = FALSE, cluster_sigs = FALSE) + theme(text = element_text(size = 20))
  plot_contribution(nmf_res_5_contri, nmf_res_5_sigs, mode = "relative") + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0))
  plot_contribution_heatmap(nmf_res_5_contri, cluster_samples = FALSE, cluster_sigs = FALSE) + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5)) + geom_tile(colour = "white")
  
  cos_sim_samples_signatures_5 <- cos_sim_matrix(COSMIC_v3.4_sigs, nmf_res_5_sigs)
  plot_cosine_heatmap(t(cos_sim_samples_signatures_5), cluster_rows = F, cluster_cols = FALSE, method = "complete", plot_values = TRUE) + theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.0)) + geom_tile(colour = "white")
dev.off()

###### signature decomposition from Signet
### check the denovo sigs similar to which COSMIC sigs (use 61 COSMIC in SigNet)
SigNet_denovo_contri_rank5 <- read.csv(paste0(project_dir, "/data/SigNet/PTA_all_denovo_rank5/weight_guesses.csv"), row.names = 1) %>% t()
p_signet_denovo_rank5 <- plot_contribution_heatmap(SigNet_denovo_contri_rank5, cluster_sigs = FALSE, cluster_samples = FALSE, method = "complete", plot_values = FALSE) + 
  scale_fill_gradient(low = "white", high = "dodgerblue3") + geom_tile(colour = "black") + theme_linedraw() + ggtitle(NULL) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), 
        panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.0))
ggsave(paste0(suppl_figure_dir, "/4-denovo_signet_contri_rank5.pdf"), plot = p_signet_denovo_rank5, width = 15, height = 5, dpi = 600)

################################################################################
### calculate the denovo sigs contribution in mut_mat_estimated_PTA_Cases (use 4 denovo in SigNet)
# SigNet_denovo_mut_mat_contri <- read.csv(paste0(main_figure_dir, "/SigNet/denovo_mut_mat/weight_guesses.csv"), row.names = 1) %>% t()
# sigNet_denovo_contri <- as.data.frame(t(SigNet_denovo_mut_mat_contri))
fit_res_loose_denovo <- fit_to_signatures(t(mut_mat_est), nmf_res_5_sigs)
loose_contri_denovo <- fit_res_loose_denovo$contribution
sigNet_denovo_contri_normalized <- as.data.frame(t(loose_contri_denovo)) %>% mutate(norm_factor = 1 / rowSums(.)) %>% 
  mutate(across(-norm_factor, ~ . * norm_factor * SCAN2_df$snv.rate.per.gb)) %>% select(-norm_factor) %>% 
  mutate(Age = metadata_df$Age, Condition = relevel(factor(metadata_df$Condition), ref = "Control"), Gender = metadata_df$Gender, Case_ID = metadata_df$Case_ID, Color = SCAN2_df$Color) %>% 
  filter(row.names(.) %in% c(filtered_mut_mat_est_Control, filtered_mut_mat_est_IHD))

refined_selected_sigs_list <- colnames(sigNet_denovo_contri_normalized)[1 : (ncol(sigNet_denovo_contri_normalized) - 5)]
sig_digits <- 2
pdf(paste0(other_figure_dir, "/4-lmer_denovo_sigs.pdf"), width = 12, height = 8)
# pdf(paste0(suppl_figure_dir, "/4-lmer_denovo_sigs.pdf"), width = 10, height = 8)
  for (i in c(1 : (ncol(sigNet_denovo_contri_normalized) - 5))){
    cat("index", i, ":", refined_selected_sigs_list[i])
    ##### Mix effects model & age matched mix effects model
    ## burden_df: the input matrix with columns: burden, Condition, Case_ID, Age, Color
    burden_df <- sigNet_denovo_contri_normalized[c(refined_selected_sigs_list[i], "Age", "Condition", "Case_ID", "Color")] %>% setNames(c("snv.rate.per.gb", "Age", "Condition", "Case_ID", "Color")) %>% 
      mutate(Condition = relevel(factor(Condition), ref = "Control")) %>% group_by(factor(Condition, levels = c(disease_name, control_name))) %>% arrange(Age, .by_group = TRUE)
    burden_df_control <- burden_df %>% filter(Condition == "Control")
    burden_df_disease <- burden_df %>% filter(Condition == "IHD")
    
    ## mixed linear modeling
    burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + (1|Case_ID), burden_df, REML = FALSE)
    burden_age_model_control <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_control, REML = FALSE)
    # burden_age_model_disease <- lmer(snv.rate.per.gb ~ age + (1|Case_ID), burden_df_disease, REML = FALSE)
    
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
      geom_point(pch = 21, fill = burden_df$Color, size = 5) + 
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
    ggsave(paste0(suppl_figure_dir, "/4-", refined_selected_sigs_list[i], "_SNV_burden_lme.pdf"), plot = p_SNV_burden_lme, width = 9, height = 6, dpi = 600)
  }
dev.off()
