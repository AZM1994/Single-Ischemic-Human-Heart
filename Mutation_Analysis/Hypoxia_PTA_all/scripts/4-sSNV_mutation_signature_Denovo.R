##### Mutational signatures extraction
mut_mat_denovo <- t(mut_mat_est_AMG) + 0.0001
COSMIC_v3.4_sigs <- read.table("main/COSMIC/COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE) %>% {rownames(.) <- .$Type; .} %>% select(-Type) %>% as.matrix()

#################### De Novo signature extraction with NMF
# estimate <- nmf(mut_mat_denovo, rank = 2:8, method = "brunet", nrun = 20, seed = 123456, .opt = "v-p")
# print(plot(estimate))

################################################################################
############## choose rank 4
nmf_res_4 <- extract_signatures(mut_mat_denovo, rank = 4, nrun = 20, single_core = TRUE)
nmf_res_4_contri <- nmf_res_4$contribution
nmf_res_4_sigs <- nmf_res_4$signatures
# nmf_res_4_sigs <- nmf_res_4_sigs[, c(1, 2, 4, 3)]
nmf_res_4_sigs <- nmf_res_4_sigs[, c(2, 3, 1, 4)]
colnames(nmf_res_4_sigs) <- c("Signature_N4", "Signature_N3", "Signature_N2", "Signature_N1")
rownames(nmf_res_4_contri) <- c("Signature_N4", "Signature_N3", "Signature_N2", "Signature_N1")
nmf_res_4_sigs <- apply(nmf_res_4_sigs, 2, function(x) x / sum(x))
write.csv(t(nmf_res_4_sigs), file = paste0(table_dir, "/4-denovo_sigs_rank4.csv"))

p_96_denovo_rank4 <- plot_96_profile(nmf_res_4_sigs, ymax = 0.12) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.15), panel.grid.minor = element_blank(), 
        panel.border = element_rect(linewidth = 0.5), text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(suppl_figure_dir, "/4-denovo_sigs_spectrum_rank4.pdf"), plot = p_96_denovo_rank4, width = 12, height = 6, dpi = 600)

nmf_res_4_contri_cond <- tibble::tibble(
  SigN4_control = sum(nmf_res_4_contri[1, ctrl_range_AMG]), SigN4_disease = sum(nmf_res_4_contri[1, dis_range_AMG]), 
  SigN3_control = sum(nmf_res_4_contri[2, ctrl_range_AMG]), SigN3_disease = sum(nmf_res_4_contri[2, dis_range_AMG]),
  SigN2_control = sum(nmf_res_4_contri[3, ctrl_range_AMG]), SigN2_disease = sum(nmf_res_4_contri[3, dis_range_AMG]), 
  SigN1_control = sum(nmf_res_4_contri[4, ctrl_range_AMG]), SigN1_disease = sum(nmf_res_4_contri[4, dis_range_AMG])) %>%
  t() %>% matrix(nrow = 4, ncol = 2, byrow = TRUE) %>% `rownames<-`(c("Signature_N4", "Signature_N3", "Signature_N2", "Signature_N1")) %>% `colnames<-`(c("Control", "IHD"))

sig_colors_denovo <- c(Signature_N4 = "#91a64f", Signature_N3 = "#fca85c", Signature_N2 = "#f2b8c0", Signature_N1 = "#8882bd")

pdf(paste0(other_figure_dir, "/4-sSNV_nmf_res_rank4_barplot.pdf"), width = 14, height = 10)
  plot_contribution(nmf_res_4_contri_cond, mode = "relative") + 
    scale_fill_manual(name = "Signature", values = sig_colors_denovo) + theme_linedraw() + 
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12))
  plot_contribution(nmf_res_4_contri, mode = "relative") + 
    scale_fill_manual(name = "Signature", values = sig_colors_denovo) + theme_linedraw() + 
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

###### decompose denovo signature to COSMIC
SigNet_denovo_contri_rank4 <- read.csv(paste0(project_dir, "/data/SigNet/PTA_all_denovo_rank4/weight_guesses.csv"), row.names = 1) %>% t()
p_signet_denovo_rank4 <- plot_contribution_heatmap(SigNet_denovo_contri_rank4, cluster_sigs = FALSE, cluster_samples = FALSE, method = "complete", plot_values = FALSE) + 
  scale_fill_gradientn(colors = c("white", "#cce6ff", "dodgerblue1", "dodgerblue2", "dodgerblue3", "dodgerblue4")) + 
  geom_tile(colour = "grey80") + theme_linedraw() + ggtitle(NULL) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.0))
ggsave(paste0(suppl_figure_dir, "/4-denovo_signet_contri_rank4.pdf"), plot = p_signet_denovo_rank4, width = 15, height = 5, dpi = 600)

################################################################################
### calculate the denovo sigs contribution in mut_mat_estimated_PTA_Cases (use 4 denovo in SigNet)
fit_res_loose_denovo <- fit_to_signatures(t(mut_mat_est), nmf_res_4_sigs)
loose_contri_denovo <- fit_res_loose_denovo$contribution
sigNet_denovo_contri_nor <- as.data.frame(t(loose_contri_denovo)) %>% mutate(norm_factor = 1 / rowSums(.)) %>% 
  mutate(across(-norm_factor, ~ . * norm_factor * sSNV_SCAN2_df$snv.rate.per.gb)) %>% select(-norm_factor) %>% 
  mutate(Age = metadata_df$Age, Condition = relevel(factor(metadata_df$Condition), ref = "Control"), Gender = metadata_df$Gender, 
         Case_ID = metadata_df$Case_ID, Color = sSNV_SCAN2_df$Color, Outline_Color = sSNV_SCAN2_df$Outline_Color)

refined_selected_sigs_list <- colnames(sigNet_denovo_contri_nor)[1 : (ncol(sigNet_denovo_contri_nor) - 6)]
sig_digits <- 2

for (i in c(1 : (ncol(sigNet_denovo_contri_nor) - 6))){
  cat("index", i, ":", refined_selected_sigs_list[i], "\n")
  ##### Mix effects model & age matched mix effects model
  ## burden_df: the input matrix with columns: burden, Condition, Case_ID, Age, Color
  burden_df <- sigNet_denovo_contri_nor[c(refined_selected_sigs_list[i], "Age", "Condition", "Case_ID", "Color", "Outline_Color")] %>% 
    setNames(c("snv.rate.per.gb", "Age", "Condition", "Case_ID", "Color", "Outline_Color")) %>% 
    # filter(snv.rate.per.gb < 5000) %>% 
    mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% 
    mutate(Condition = relevel(Condition, ref = ctrl_name)) %>% 
    group_by(Condition) %>% arrange(Age, .by_group = TRUE)
  burden_df_ctrl <- burden_df %>% filter(Condition == "Control")
  
  ## mixed linear modeling
  burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + (1|Case_ID), burden_df, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
  burden_age_model_ctrl <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_ctrl, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
  
  ## compute confident interval
  # CI_age_model <- suppressMessages(confint(burden_age_model, method = "Wald", oldNames = FALSE))
  # CI_age_model_control <- suppressMessages(confint(burden_age_model_ctrl, method = "Wald", oldNames = FALSE))
  
  ## check p value and R^2
  anova_pvalue <- anova(burden_age_model)$"Pr(>F)"[2]
  anova_pvalue_print <- formatC(signif(anova_pvalue, digits = sig_digits), digits = sig_digits, format="fg", flag="#")
  anova_pvalue_ctrl <- anova(burden_age_model_ctrl)$"Pr(>F)"[1]
  anova_pvalue_print_ctrl <- formatC(signif(anova_pvalue_ctrl, digits = sig_digits), digits = sig_digits, format="e", flag="#")
  
  ## manually calculate the fitting lines
  geom_line_data <- burden_df_ctrl %>% {
    ctrl_fitting_x <- range(.$Age)
    ctrl_fitting_y <- ctrl_fitting_x * fixef(burden_age_model_ctrl)[2] + fixef(burden_age_model_ctrl)[1]
    geom_line_data <- rbind(ctrl_fitting_x, ctrl_fitting_y)
    as.data.frame(t(geom_line_data))}
  
  aging_rate <- format(round(fixef(burden_age_model_ctrl)[2], digits = sig_digits), nsmall = sig_digits)
  legend_data <- burden_df[c(1, nrow(burden_df)), c("Age", "snv.rate.per.gb", "Condition")]
  
  generate_y_breaks <- function(y_data) {
    y_min <- min(y_data, na.rm = TRUE)
    y_max <- max(y_data, na.rm = TRUE)
    range_span <- y_max - y_min
    step_size <- signif(range_span / 4, 1)
    breaks <- seq(floor(y_min / step_size) * step_size, ceiling(y_max / step_size) * step_size, by = step_size)
    return(breaks)
  }
  
  yticks_plot <- generate_y_breaks(burden_df$snv.rate.per.gb)
  ylim_plot <- range(yticks_plot, na.rm = TRUE)
  p_SNV_burden_lme <- ggplot(burden_df, aes(x = Age, y = snv.rate.per.gb)) + 
    geom_point(pch = 21, data = legend_data, aes(x = Age, y = snv.rate.per.gb, color = Condition, fill = Condition), size = 5) + 
    geom_point(pch = 21, color = burden_df$Outline_Color, fill = burden_df$Color, size = 5) + 
    geom_line(aes(x = ctrl_fitting_x, y = ctrl_fitting_y), color = "dodgerblue2", linetype = "22", data = geom_line_data) + 
    geom_hline(yintercept = 0, color = "gray40", linetype = "22") + 
    annotate("text", size = 6, x = 0, y = 0.95 * ylim_plot[2], 
             label = paste("Aging effect:", aging_rate, "sSNVs/(GB·year),", 
                           "P =", anova_pvalue_print_ctrl), hjust = 0) + 
    annotate("text", size = 6, x = 0, y = 0.87 * ylim_plot[2], 
             label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB,", 
                           "P =", anova_pvalue_print), hjust = 0) +
    scale_color_manual(values = c("IHD" = "black", "Control" = "black"), guide = "legend") + 
    scale_fill_manual(values = c("Control" = ctrl_dis_color[1], "IHD" = ctrl_dis_color[2]), guide = "legend") + 
    scale_y_continuous(limits = c(min(min(geom_line_data[, 2]), 0), ylim_plot[2]), breaks = yticks_plot, labels = yticks_plot) +
    labs(x = "Age (years)", y = paste0(refined_selected_sigs_list[i], " Contribution \n (sSNV rate per GB)")) + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
  print(p_SNV_burden_lme)
  ggsave(paste0(suppl_figure_dir, "/4-", refined_selected_sigs_list[i], "_SNV_burden_lme.pdf"), plot = p_SNV_burden_lme, width = 9, height = 6, dpi = 600)
}
