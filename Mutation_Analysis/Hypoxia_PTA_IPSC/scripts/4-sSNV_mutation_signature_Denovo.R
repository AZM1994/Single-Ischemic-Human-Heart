##### Mutational signatures extraction
mut_mat_denovo <- t(mut_mat_est) + 0.0001
COSMIC_v3.4_sigs <- read.table("main/COSMIC/COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE) %>% {rownames(.) <- .$Type; .} %>% select(-Type) %>% as.matrix()

#################### De Novo signature extraction with NMF
# estimate <- nmf(mut_mat_denovo, rank = 2:8, method = "brunet", nrun = 10, seed = 123456, .opt = "v-p")
# print(plot(estimate))

################################################################################
############## choose rank 3
nmf_res_2 <- extract_signatures(mut_mat_denovo, rank = 2, nrun = 20, single_core = TRUE)
nmf_res_2_contri <- nmf_res_2$contribution
nmf_res_2_sigs <- nmf_res_2$signatures
# nmf_res_2_sigs <- nmf_res_2_sigs[, c(1, 2, 4, 3)]
colnames(nmf_res_2_sigs) <- c("Signature_N2", "Signature_N1")
rownames(nmf_res_2_contri) <- c("Signature_N2", "Signature_N1")
nmf_res_2_sigs <- apply(nmf_res_2_sigs, 2, function(x) x / sum(x))
write.csv(t(nmf_res_2_sigs), file = paste0(table_dir, "/4-denovo_sigs_rank2.csv"))

p_96_denovo_rank2 <- plot_96_profile(nmf_res_2_sigs, ymax = 0.12) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.15), panel.grid.minor = element_blank(), 
        panel.border = element_rect(linewidth = 0.5), text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(suppl_figure_dir, "/4-denovo_sigs_spectrum_rank2.pdf"), plot = p_96_denovo_rank2, width = 12, height = 6, dpi = 600)

nmf_res_2_contri_cond <- tibble::tibble(
  SigN2_control = sum(nmf_res_2_contri[1, ctrl_range]), SigN2_disease = sum(nmf_res_2_contri[1, dis_range]), 
  SigN1_control = sum(nmf_res_2_contri[2, ctrl_range]), SigN1_disease = sum(nmf_res_2_contri[2, dis_range])) %>%
  t() %>% matrix(nrow = 2, ncol = 2, byrow = TRUE) %>% `rownames<-`(c("Signature_N2", "Signature_N1")) %>% `colnames<-`(c(ctrl_name, dis_name))

sig_colors_denovo <- c(Signature_N4 = "#91a64f", Signature_N3 = "#fca85c", Signature_N2 = "#f2b8c0", Signature_N1 = "#8882bd")

pdf(paste0(other_figure_dir, "/4-sSNV_nmf_res_rank2_barplot.pdf"), width = 14, height = 10)
  plot_contribution(nmf_res_2_contri_cond, mode = "relative") + 
    scale_fill_manual(name = "Signature", values = sig_colors_denovo) + theme_linedraw() + 
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12))
  plot_contribution(nmf_res_2_contri, mode = "relative") + 
    scale_fill_manual(name = "Signature", values = sig_colors_denovo) + theme_linedraw() + 
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

###### decompose denovo signature to COSMIC
SigNet_denovo_contri_rank2 <- read.csv(paste0(project_dir, "/data/SigNet/PTA_IPSC_denovo_rank2/weight_guesses.csv"), row.names = 1) %>% t()
p_signet_denovo_rank2 <- plot_contribution_heatmap(SigNet_denovo_contri_rank2, cluster_sigs = FALSE, cluster_samples = FALSE, method = "complete", plot_values = FALSE) + 
  scale_fill_gradientn(colors = c("white", "#cce6ff", "dodgerblue1", "dodgerblue2", "dodgerblue3", "dodgerblue4")) + 
  geom_tile(colour = "grey80") + theme_linedraw() + ggtitle(NULL) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.0))
ggsave(paste0(suppl_figure_dir, "/4-denovo_signet_contri_rank2.pdf"), plot = p_signet_denovo_rank2, width = 15, height = 5, dpi = 600)

################################################################################
### calculate the denovo sigs contribution in mut_mat_estimated_PTA_Cases (use 2 denovo in SigNet)
fit_res_loose_denovo <- fit_to_signatures(t(mut_mat_est), nmf_res_2_sigs)
loose_contri_denovo <- fit_res_loose_denovo$contribution
sigNet_denovo_contri_nor <- as.data.frame(t(loose_contri_denovo)) %>% mutate(norm_factor = 1 / rowSums(.)) %>% 
  mutate(across(-norm_factor, ~ . * norm_factor * sSNV_SCAN2_df$snv.rate.per.gb)) %>% select(-norm_factor) %>% 
  mutate(Age = metadata_df$Age, Condition = relevel(factor(metadata_df$Condition), ref = ctrl_name), Gender = metadata_df$Gender, 
         Color = sSNV_SCAN2_df$Color, Outline_Color = sSNV_SCAN2_df$Outline_Color)

refined_selected_sigs_list <- colnames(sigNet_denovo_contri_nor)[1 : (ncol(sigNet_denovo_contri_nor) - 5)]
sig_digits <- 2

pdf(paste0(other_figure_dir, "/4-lmer_top_sigs_denovo.pdf"), width = 8, height = 6)
  for (i in c(1 : (ncol(sigNet_denovo_contri_nor) - 5))){
    cat("index", i, ":", refined_selected_sigs_list[i], "\n")
    ## burden_df: the input matrix with columns: burden, Condition, Case_ID, Age, Color
    burden_df <- sigNet_denovo_contri_nor[c(refined_selected_sigs_list[i], "Age", "Condition", "Color", "Outline_Color")] %>% 
      setNames(c("snv.rate.per.gb", "Age", "Condition", "Color", "Outline_Color")) %>% 
      mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% 
      mutate(Condition = relevel(Condition, ref = ctrl_name)) %>% 
      mutate(label = c(rep("Normoxia (n = 9)", 9), rep("Hypoxia (n = 9)", 9))) %>% 
      group_by(Condition) %>% arrange(Age, .by_group = TRUE)
    burden_df_ctrl <- burden_df %>% filter(Condition == ctrl_name)
    burden_df_dis <- burden_df %>% filter(Condition == dis_name)
    
    wilcox.test_result_sig <- wilcox.test(burden_df_ctrl$snv.rate.per.gb, burden_df_dis$snv.rate.per.gb, alternative = "two.sided")$p.value
    if (wilcox.test_result_sig < 0.0099){
      wilcox.test_result_sig_print <- format(wilcox.test_result_sig, scientific = TRUE, digits = 2)
    } else{
      wilcox.test_result_sig_print <- formatC(signif(wilcox.test_result_sig, digits = 2), digits = 2, format="fg", flag="#")
    }
    
    boxplot_sSNV_burden_sig <- ggplot(burden_df, aes(x = factor(label, level = c("Normoxia (n = 9)", "Hypoxia (n = 9)")), y = snv.rate.per.gb, fill = label)) + 
      geom_boxplot(color = c("dodgerblue3", "firebrick3"), fill = c("#C7DFF0", "#FBDFE2")) + 
      geom_quasirandom(pch = 21, fill = burden_df$Color, color = burden_df$Outline_Color, size = 5) + 
      annotate("text", size = 7, x = 1.32, y = max(burden_df$snv.rate.per.gb), label = paste("P = ", wilcox.test_result_sig_print), hjust = 0) + 
      labs(x = "", y = paste0(refined_selected_sigs_list[i], " Contribution \n (sSNV rate per GB)")) + 
      stat_summary(fun = mean, colour = "black", geom = "point", shape = 18, size = 3, show.legend = FALSE) + 
      stat_compare_means(comparisons = list(c("Normoxia (n = 9)", "Hypoxia (n = 9)")), label.x = 1.6, label.y = max(burden_df$snv.rate.per.gb), 
                         label = "p.signif", bracket.size = 1, size = 7, method = "wilcox.test") + 
      scale_y_continuous(labels = number_format(accuracy = 0.1)) + theme_linedraw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), 
            text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
    print(boxplot_sSNV_burden_sig)
    if (wilcox.test_result_sig <= 0.05){
      cat("Significant \n")
      ggsave(paste0(suppl_figure_dir, "/2-boxplot_sSNV_burden_", refined_selected_sigs_list[i], ".pdf"), plot = boxplot_sSNV_burden_sig, width = 9, height = 8)
    }
  }
dev.off()