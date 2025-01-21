##### Mutational signatures extraction
# mut_mat_age_match <- mut_mat_age_match + 0.0001
mut_mat_estimiated_denovo <- mut_mat_estimiated + 0.0001
# signatures = get_known_signatures()
COSMIC_v3.4_sigs <- read.table("main/COSMIC/COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE) %>%
  {rownames(.) <- .$Type; .} %>% select(-Type) %>% as.matrix()

#################### De Novo signature extraction with NMF
# estimate_2 <- nmf(mut_mat_estimiated_denovo, rank = 2:7, method = "brunet", nrun = 10, seed = 123456, .opt = "v-p")
# print(plot(estimate_2))

############## choose rank 4
nmf_res_4 <- extract_signatures(mut_mat_estimiated_denovo, rank = 4, nrun = 10, single_core = TRUE)
nmf_res_4_contri <- nmf_res_4$contribution
nmf_res_4_sigs <- nmf_res_4$signatures
nmf_res_4_sigs <- nmf_res_4_sigs[, c(4, 3, 1, 2)]
colnames(nmf_res_4_sigs) <- c("Signature N4", "Signature N3", "Signature N2", "Signature N1")
rownames(nmf_res_4_contri) <- c("Signature N4", "Signature N3", "Signature N2", "Signature N1")
# nmf_res_4 <- rename_nmf_signatures(nmf_res_4, signatures, cutoff = 0.88)
# nmf_res_4 <- rename_nmf_signatures(nmf_res_4, COSMIC_v3.4_sigs, cutoff = 0.88)
colnames(nmf_res_4_sigs)
nmf_res_4_sigs <- apply(nmf_res_4_sigs, 2, function(x) x / sum(x))
write.csv(t(nmf_res_4_sigs), file = paste0(sSNV_figure_dir, "/denovo_sigs.csv"))
write.csv(nmf_res_4_sigs, file = paste0(sSNV_figure_dir, "/denovo_sigs_02.csv"))

pdf(paste0(other_figure_dir, "/4-Visualizing_nmf_res_rank4.pdf"), width = 14, height = 10)
  p_96_denovo <- plot_96_profile(nmf_res_4_sigs) + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.15), panel.grid.minor = element_blank(), 
          panel.border = element_rect(size = 0.5), text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(paste0(suppl_figure_dir, "/4-denovo_sigs_spectrum.pdf"), plot = p_96_denovo, width = 12, height = 6, dpi = 600)
  
  SBSN4_control <- sum(nmf_res_4_contri[1, control_range_age_match])
  SBSN4_disease <- sum(nmf_res_4_contri[1, disease_range_age_match])
  SBSN3_control <- sum(nmf_res_4_contri[2, control_range_age_match])
  SBSN3_disease <- sum(nmf_res_4_contri[2, disease_range_age_match])
  SBSN2_control <- sum(nmf_res_4_contri[3, control_range_age_match])
  SBSN2_disease <- sum(nmf_res_4_contri[3, disease_range_age_match])
  SBSN1_control <- sum(nmf_res_4_contri[4, control_range_age_match])
  SBSN1_disease <- sum(nmf_res_4_contri[4, disease_range_age_match])

  nmf_res_4_contri_02 <- matrix(c(SBSN4_control,SBSN3_control,SBSN2_control,SBSN1_control,SBSN4_disease,SBSN3_disease,SBSN2_disease,SBSN1_disease),4,2)
  rownames(nmf_res_4_contri_02) = c("SigN4", "SigN3", "SigN2", "SigN1")
  colnames(nmf_res_4_contri_02) = c("Control", "IHD")
  plot_contribution(nmf_res_4_contri_02, nmf_res_4_sigs, mode = "relative") + theme(text = element_text(size=20))
  plot_contribution_heatmap(nmf_res_4_contri_02, cluster_samples = FALSE, cluster_sigs = FALSE) +
    theme(text = element_text(size=20), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5)) + geom_tile(colour = "white")
  
  plot_contribution(nmf_res_4_contri, nmf_res_4_sigs, mode = "relative") +
    theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0))

  plot_contribution_heatmap(nmf_res_4_contri, cluster_samples = FALSE, cluster_sigs = FALSE) + 
    theme(text = element_text(size=20), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5)) + geom_tile(colour = "white")

  # cos_sim_samples_signatures_4 <- cos_sim_matrix(COSMIC_v3.4_sigs, nmf_res_4_sigs)
  # p_cosine <- plot_cosine_heatmap(t(cos_sim_samples_signatures_4), cluster_rows = F, cluster_cols = FALSE, method = "complete", plot_values = TRUE) +
  #   theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.0)) + geom_tile(colour = "white")
  # print(p_cosine)
  # ggsave(paste0(sSNV_figure_dir, "/4-denovo-cosine.png"), plot = p_cosine, width = 12, height = 5, dpi = 600)
dev.off()

###### signature decomposition from Signet
### check the denovo sigs similar to which COSMIC sigs (use 61 COSMIC in SigNet)
SigNet_denovo_contri <- read.csv(paste0(main_figure_dir, "/SigNet/denovo_COSMIC/weight_guesses.csv"), row.names = 1) %>% t()
p_signet_denovo <- plot_contribution_heatmap(SigNet_denovo_contri, cluster_sigs = FALSE, cluster_samples = FALSE, method = "complete", plot_values = FALSE) + 
        scale_fill_gradient(low = "white", high = "dodgerblue3") + geom_tile(colour = "black") + theme_linedraw() + ggtitle(NULL) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), 
        panel.grid.minor = element_blank(), text = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.0))

print(p_signet_denovo)
ggsave(paste0(suppl_figure_dir, "/4-denovo-signet-contri.pdf"), plot = p_signet_denovo, width = 15, height = 5, dpi = 600)

### calculate the denovo sigs contribution in mut_mat_estimated_PTA_Cases (use 4 denovo in SigNet)
# SigNet_denovo_mut_mat_contri <- read.csv(paste0(main_figure_dir, "/SigNet/denovo_mut_mat/weight_guesses.csv"), row.names = 1) %>% t()
# sigNet_denovo_contri <- as.data.frame(t(SigNet_denovo_mut_mat_contri))
# sigNet_denovo_contri <- as.data.frame(t(nmf_res_4_contri))
fit_res_loose_denovo <- fit_to_signatures(mut_mat_estimiated, nmf_res_4_sigs)
loose_contri_denovo <- fit_res_loose_denovo$contribution
sigNet_denovo_contri <- as.data.frame(t(loose_contri_denovo))
sigNet_denovo_contri_normalized <- 1 / rowSums(sigNet_denovo_contri) * sigNet_denovo_contri * SCAN2_df$snv.rate.per.gb
sigNet_denovo_contri_normalized <- slice(sigNet_denovo_contri_normalized, match(Cell_ID_list, rownames(sigNet_denovo_contri_normalized)))

sigNet_denovo_contri_normalized_sum <- sigNet_denovo_contri_normalized
# signature_cols_N1 <- grep("Signature.N1", names(sigNet_denovo_contri_normalized))
# signature_cols_N2 <- grep("Signature.N2", names(sigNet_denovo_contri_normalized))
# signature_cols_N3 <- grep("Signature.N3", names(sigNet_denovo_contri_normalized))
# signature_cols_N4 <- grep("Signature.N4", names(sigNet_denovo_contri_normalized))

# sigNet_denovo_contri_normalized_sum <- data.frame(
#   SignatureN1 = rowSums(sigNet_denovo_contri_normalized[, signature_cols_N1]),
#   SignatureN2 = rowSums(sigNet_denovo_contri_normalized[, signature_cols_N2]),
#   SignatureN3 = rowSums(sigNet_denovo_contri_normalized[, signature_cols_N3]),
#   SignatureN4 = rowSums(sigNet_denovo_contri_normalized[, signature_cols_N4])
# )

# sigNet_denovo_contri_normalized_sum <- sigNet_denovo_contri_normalized_sum[,1:4]
sigNet_denovo_contri_normalized_sum <- sigNet_denovo_contri_normalized_sum %>% 
  mutate(age = metadata_df$Age) %>%
  mutate(condition = metadata_df$Condition) %>%
  mutate(gender = metadata_df$Gender) %>%
  mutate(Case_ID = metadata_df$Case_ID) %>% 
  mutate(condition = relevel(factor(condition), ref = "Control"))

refined_selected_sigs_list <- colnames(sigNet_denovo_contri_normalized_sum)[1 : (ncol(sigNet_denovo_contri_normalized_sum) - 4)]

# pdf(paste0(other_figure_dir, "/4-lmer_denovo_sigs.pdf"), width = 10, height = 8)
pdf(paste0(suppl_figure_dir, "/4-lmer_denovo_sigs.pdf"), width = 10, height = 8)
  for (i in c(1,2,3,4)){
  print(i)
  print(refined_selected_sigs_list[i])
  signature_burden_df <- sigNet_denovo_contri_normalized_sum[c(refined_selected_sigs_list[i],'age','condition','Case_ID')]
  colnames(signature_burden_df)[1] <- 'snv.rate.per.gb'
  sig_digits = 2
  ##### Mix effects model & age matched mix effects model
  ## burden_df: the input matrix with columns: burden, condition, Case_ID, age
  burden_df <- signature_burden_df
  burden_df <- burden_df %>%
    mutate(condition = relevel(factor(condition), ref = "Control")) %>%
    group_by(factor(condition, levels = c(disease_name, control_name))) %>% arrange(age, .by_group = TRUE)
  burden_df_control <- burden_df[burden_df$condition == "Control",]
  burden_df_disease <- burden_df[burden_df$condition == "IHD",]
  
  ## mixed linear modeling
    burden_age_model <- lmer(snv.rate.per.gb ~ age + condition + (1|Case_ID), burden_df, REML = FALSE)
    burden_age_model_control <- lmer(snv.rate.per.gb ~ age + (1|Case_ID), burden_df_control, REML = FALSE)
    burden_age_model_disease <- lmer(snv.rate.per.gb ~ age + (1|Case_ID), burden_df_disease, REML = FALSE)
    
    ## check model parameters
    summary_lm_all <- summary(burden_age_model)
    summary_lm_control <- summary(burden_age_model_control)
    # confint(burden_age_model_control, oldNames=FALSE)
    summary_lm_disease <- summary(burden_age_model_disease)
    
    ## compute confident interval
    CI_age_model <- confint(burden_age_model, oldNames = FALSE)
    CI_age_model_control <- confint(burden_age_model_control, oldNames = FALSE)
    
    ## check p value and R^2
    anova_pvalue <- anova(burden_age_model)$"Pr(>F)"[2]
    if (refined_selected_sigs_list[i] %in% c("SBS32", "SBS89")){
      # if (anova_pvalue < 0.0099){
      # anova_pvalue_print <- format(anova_pvalue, scientific = TRUE, digits = 2)
      # anova_pvalue_print <- formatC(signif(anova_pvalue, digits = 2), digits = 1, format="e", flag="#")
      anova_pvalue_print <- formatC(signif(anova_pvalue, digits = 2), digits = 2, format="fg", flag="#")
    } else{
      anova_pvalue_print <- formatC(signif(anova_pvalue, digits = 2), digits = 2, format="fg", flag="#")
    }
    anova_pvalue_control <- anova(burden_age_model_control)$"Pr(>F)"[1]
    anova_pvalue_print_control <- formatC(signif(anova_pvalue_control, digits = 2), digits = 2, format="e", flag="#")
    mantissa <- format(as.numeric(sub("e.*", "", anova_pvalue_print_control)), nsmall = 2)
    exponent <- as.numeric(sub(".*e([+-]?\\d+)", "\\1", anova_pvalue_print_control))
    
    # all_mix_effect_model_p_value <- cbind(all_mix_effect_model_p_value, anova_pvalue)
    
    ## manually calculate the fitting lines
    Control_fitting_x = range(burden_df_control$age)
    Control_fitting_y = Control_fitting_x * fixef(burden_age_model_control)[2] + fixef(burden_age_model_control)[1]
    geom_line_data <- rbind(Control_fitting_x, Control_fitting_y)
    geom_line_data <- as.data.frame(t(geom_line_data))
    
    # aging_rate <- formatC(signif(fixef(burden_age_model_control)[2], digits = sig_digits), digits = sig_digits, format="fg", flag="#")
    aging_rate <- format(round(fixef(burden_age_model_control)[2], digits = sig_digits), nsmall = sig_digits)
    legend_data <- burden_df[c(7,20), c("age", "snv.rate.per.gb", "condition")]
    p_SNV_burden_lme <- ggplot(burden_df, aes(x = age, y = snv.rate.per.gb), color = donors) + 
      geom_point(pch = 21, data = legend_data, aes(x = age, y = snv.rate.per.gb, color = condition, fill = condition), size = 5) + 
      geom_point(pch = 21, fill = c(disease_color, control_color), size = 5) + 
      geom_line(aes(x = Control_fitting_x, y = Control_fitting_y), colour = "dodgerblue3", data = geom_line_data) + 
      annotate("text", size = 7, x = 0.0 * max(burden_df$age), y = 1.0 * max(burden_df$snv.rate.per.gb), label = paste("aging effect:", aging_rate, "sSNVs/(GB·year)"), hjust = 0, color = dis_ctrl_color[1]) + 
      annotate("text", size = 7, x = 0.002 * max(burden_df$age), y = 0.93 * max(burden_df$snv.rate.per.gb), label = bquote("p =" ~ .(mantissa) ~ "×" ~ 10^{.(exponent)}), hjust = 0, color = dis_ctrl_color[1]) + 
      annotate("text", size = 7, x = 0.002 * max(burden_df$age), y = 0.87 * max(burden_df$snv.rate.per.gb), label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB"), hjust = 0, color = dis_ctrl_color[2]) + 
      annotate("text", size = 7, x = 0.0 * max(burden_df$age), y = 0.8 * max(burden_df$snv.rate.per.gb), label = paste("p = ", anova_pvalue_print), hjust = 0, color = dis_ctrl_color[2]) + 
      # annotate("text", size = 7, x = 0.0 * max(burden_df$age), y = 0.93 * max(burden_df$snv.rate.per.gb), label = paste("p = ", anova_pvalue_print), hjust = 0) + 
      # annotate("text", size = 6, x = 0.0 * max(burden_df$age), y = 0.85 * max(burden_df$snv.rate.per.gb),
      #          label = paste0("IHD excess SNVs = ", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "/GB, ",
      #                        "95% CI = [", format(round(CI_age_model[5, 1], digits = sig_digits), nsmall = sig_digits), ", ",
      #                                               format(round(CI_age_model[5, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) +
      # annotate("text", size = 6, x = 0.0 * max(burden_df$age), y = 0.90 * max(burden_df$snv.rate.per.gb),
      #          label = paste0("Age accum. rate = ", format(round(fixef(burden_age_model_control)[2], digits = sig_digits), nsmall = sig_digits), "/(GB.year), ",
      #                         "95% CI = [", format(round(CI_age_model_control[4, 1], digits = sig_digits), nsmall = sig_digits), ", ",
      #                         format(round(CI_age_model_control[4, 2], digits = sig_digits), nsmall = sig_digits), "]"), hjust = 0) +
      # annotate("text", size = 7, x = 0.002 * max(burden_df$age), y = 1.0 * max(burden_df$snv.rate.per.gb), label = bquote("s = " ~ .(aging_rate) ~ GB^{phantom()-1}*year^{phantom()-1}), hjust = 0) + 
      scale_color_manual(values = c("IHD" = "black", "Control" = "black"), guide = "legend") + scale_fill_manual(values = c("Control" = dis_ctrl_color[1], "IHD" = dis_ctrl_color[2]), guide = "legend") + 
      labs(x = "Age (years)", y = paste0(refined_selected_sigs_list[i], " Contribution \n (sSNV rate per GB)")) + theme_linedraw() + 
      theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
            panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0))
    print(p_SNV_burden_lme)
    ggsave(paste0(suppl_figure_dir, "/4-", refined_selected_sigs_list[i], "_SNV_burden_lme.pdf"), plot = p_SNV_burden_lme, width = 9, height = 6, dpi = 600)
}
dev.off()
