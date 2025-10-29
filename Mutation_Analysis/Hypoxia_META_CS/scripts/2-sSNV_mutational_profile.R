##### 96 mutational profile
mut_mat_raw <- mut_matrix(vcf_list = ssnv_grl_DS, ref_genome = ref_genome)
mut_mat_raw <- t(mut_mat_raw[, Cell_ID_list])
# mut_num_summary <- as.data.frame(sSNV_SCAN2_df) %>% dplyr::select(Cell_ID, Condition, snv.burden) %>% mutate(true_call = rowSums(mut_mat_raw))

##### save mut_mat for SigNet (raw)
write.csv(mut_mat_raw, file = paste0(table_dir, "/3-mut_mat_raw_META_CS.csv"))
mut_mat_raw_AMG <- mut_mat_raw[Cell_ID_list_AMG, ]
mut_mat_raw_AMG_cond <- mut_mat_raw_AMG %>% t() %>% as.data.frame() %>% 
  mutate(Control = rowMeans(.[, ctrl_range_AMG]), IHD = rowMeans(.[, dis_range_AMG]), 
         Net_change = pmax(IHD - Control, 0, na.rm = TRUE))
write.csv(mut_mat_raw_AMG_cond, file = paste0(table_dir, "/3-mut_mat_raw_AMG_cond_META_CS.csv"))

# p_96_raw_rel <- plot_96_profile(mut_mat_raw_AMG_cond[, c(ctrl_name, dis_name, "Net_change")], ymax = 0.10) + theme_linedraw() +
p_96_raw_rel <- plot_96_profile(mut_mat_raw_AMG_cond[, c(ctrl_name, dis_name)], ymax = 0.10) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.15), panel.grid.minor = element_blank(), 
        panel.border = element_rect(linewidth = 0.5), text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(main_figure_dir, "/3-sSNV_mutational_profile_raw_relative.pdf"), plot = p_96_raw_rel, width = 10, height = 4.25, dpi = 600)
