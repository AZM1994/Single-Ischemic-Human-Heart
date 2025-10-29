##### 96 mutational profile
mut_mat_raw <- mut_matrix(vcf_list = ssnv_grl, ref_genome = ref_genome)
mut_mat_raw <- t(mut_mat_raw[, Cell_ID_list])
# mut_num_summary <- as.data.frame(sSNV_SCAN2_df) %>% dplyr::select(Cell_ID, Condition, snv.burden) %>% mutate(true_call = rowSums(mut_mat_raw))

##### save mut_mat for SigNet (raw)
write.csv(mut_mat_raw, file = paste0(table_dir, "/3-mut_mat_raw_PTA_IPSC.csv"))
mut_mat_raw <- mut_mat_raw[Cell_ID_list, ]
mut_mat_raw_cond <- mut_mat_raw %>% t() %>% as.data.frame() %>% 
  mutate(Normoxia = rowMeans(.[, ctrl_range]), Hypoxia = rowMeans(.[, dis_range]), 
         Net_change = pmax(Hypoxia - Normoxia, 0, na.rm = TRUE))
write.csv(mut_mat_raw_cond, file = paste0(table_dir, "/3-mut_mat_raw_AMG_cond_PTA_IPSC.csv"))

##### save mut_mat for SigNet (normalize to est snv.burden)
mut_mat_est <- 1 / rowSums(mut_mat_raw) * mut_mat_raw * sSNV_SCAN2_df$snv.rate.per.gb
write.csv(mut_mat_est, file = paste0(table_dir, "/3-mut_mat_est_PTA_IPSC.csv"))
mut_mat_est <- mut_mat_est[Cell_ID_list, ]
mut_mat_est_cond <- mut_mat_est %>% t() %>% as.data.frame() %>% 
  mutate(Normoxia = rowMeans(.[, ctrl_range]), Hypoxia = rowMeans(.[, dis_range]), 
         Net_change = pmax(Hypoxia - Normoxia, 0, na.rm = TRUE))
write.csv(mut_mat_est_cond, file = paste0(table_dir, "/3-mut_mat_est_cond_PTA_IPSC.csv"))

p_96_raw_rel <- plot_96_profile(mut_mat_raw_cond[, c(ctrl_name, dis_name, "Net_change")], ymax = 0.10) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(linewidth = 0.5), text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(other_figure_dir, "/3-sSNV_mutational_profile_raw_relative.pdf"), plot = p_96_raw_rel, width = 10, height = 6, dpi = 600)

p_96_est_abs <- plot_96_profile_abs(mut_mat_est_cond[, c(ctrl_name, dis_name, "Net_change")]) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(linewidth = 0.5), text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(main_figure_dir, "/3-sSNV_mutational_profile_est_abs.pdf"), plot = p_96_est_abs, width = 10, height = 6, dpi = 600)
