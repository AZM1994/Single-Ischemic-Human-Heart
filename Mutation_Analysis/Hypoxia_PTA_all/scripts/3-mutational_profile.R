##### 96 mutational profile
mut_mat <- mut_matrix(vcf_list = snv_grl, ref_genome = ref_genome)
mut_mat <- t(mut_mat[, Cell_ID_list])
mut_num_summary <- as.data.frame(SCAN2_df) %>% dplyr::select(Cell_ID, Condition, snv.burden) %>% mutate(true_call = rowSums(mut_mat))

##### generate mut_mat for SigNet (raw)
mut_mat_raw <- mut_mat
write.csv(mut_mat_raw, file = paste0(table_dir, "/mut_mat_raw_PTA_all.csv"))
mut_mat_raw_AMG <- mut_mat_raw[Cell_ID_list_AMG, ]
write.csv(mut_mat_raw_AMG, file = paste0(table_dir, "/mut_mat_raw_age_match_PTA_all.csv"))
# mut_mat_raw_cond <- mut_mat_raw %>% t() %>% as.data.frame() %>% mutate(Control = rowSums(.[, ctrl_range]), IHD = rowSums(.[, dis_range])) %>% dplyr::select(Control, IHD) %>% t()
# write.csv(mut_mat_raw_cond, file = paste0(table_dir, "/mut_mat_raw_conditional_PTA_all.csv"))

##### generate mut_mat for SigNet (normalize mut_mat to est snv.burden (per cell))
mut_mat_est <- 1 / rowSums(mut_mat) * mut_mat * SCAN2_df$snv.rate.per.gb
write.csv(mut_mat_est, file = paste0(table_dir, "/mut_mat_est_PTA_all.csv"))
# mut_mat_est_cond <- mut_mat_est %>% t() %>% as.data.frame() %>% mutate(Control = rowSums(.[, ctrl_range]), IHD = rowSums(.[, dis_range])) %>% dplyr::select(Control, IHD) %>% t()
# write.csv(mut_mat_est_cond, file = paste0(table_dir, "/mut_mat_est_conditional_PTA_all.csv"))

##### age-match: normalize mut_mat to est snv.burden (per cell)
mut_mat_est_AMG <- mut_mat_est[Cell_ID_list_AMG, ]
write.csv(mut_mat_est_AMG, file = paste0(table_dir, "/mut_mat_est_age_match_PTA_all.csv"))
mut_mat_est_AMG_cond <- mut_mat_est_AMG %>% t() %>% as.data.frame() %>% mutate(Control = rowMeans(.[, ctrl_range_AMG]), IHD = rowMeans(.[, dis_range_AMG]), Net_change = pmax(IHD - Control, 0, na.rm = TRUE))

p_96 <- plot_96_profile_abs(mut_mat_est_AMG_cond[, c(ctrl_name, dis_name, "Net_change")]) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.15), panel.grid.minor = element_blank(), 
        panel.border = element_rect(linewidth = 0.5), text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(main_figure_dir, "/3-sSNV_mutational_profile.pdf"), plot = p_96, width = 10, height = 6, dpi = 600)
# write.csv(mut_mat_est_AMG_cond[, c(ctrl_name, dis_name, "Net_change")], file = paste0(table_dir, "/mut_mat_age_match_conditional_PTA_all.csv"))
