##### 96 mutational profile
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
mut_mat <- t(mut_mat[, Cell_ID_list])
mut_num_summary <- as.data.frame(SCAN2_df) |> base::`[`(c("cell_ID", "condition", "snv.burden")) %>% 
  mutate(true_call = rowSums(mut_mat))

## generate mut_mat for SigNet (raw)
mut_mat_raw <- mut_mat
write.csv(mut_mat_raw, file = paste0(sSNV_figure_dir, "/mut_mat_raw_PTA_Cases.csv"))
mut_mat_raw <- t(mut_mat_raw)
mut_mat_raw_cond <- as.data.frame(mut_mat_raw) %>% 
  mutate(Control = rowSums(mut_mat_raw[, control_range])) %>% 
  mutate(IHD = rowSums(mut_mat_raw[, disease_range])) %>% 
  mutate(Net_change = IHD - Control) %>% 
  mutate(Net_change = ifelse(Net_change < 0, NA, Net_change)) |> base::`[`(c("Control", "IHD")) %>% t()
write.csv(mut_mat_raw_cond, file = paste0(sSNV_figure_dir, "/mut_mat_raw_conditional_PTA_Cases.csv"))

## generate mut_mat for SigNet (normalize mut_mat to est snv.burden (per cell))
mut_mat_estimiated <- 1 / rowSums(mut_mat) * mut_mat * SCAN2_df$snv.rate.per.gb
write.csv(mut_mat_estimiated, file = paste0(sSNV_figure_dir, "/mut_mat_estimiated_PTA_Cases.csv"))
mut_mat_estimiated <- t(mut_mat_estimiated)
mut_mat_estimiated_cond <- as.data.frame(mut_mat_estimiated) %>% 
  mutate(Control = rowSums(mut_mat_estimiated[, control_range])) %>% 
  mutate(IHD = rowSums(mut_mat_estimiated[, disease_range])) %>% 
  mutate(Net_change = IHD - Control) %>% 
  mutate(Net_change = ifelse(Net_change < 0, NA, Net_change)) |> base::`[`(c("Control", "IHD")) %>% t()
write.csv(mut_mat_estimiated_cond, file = paste0(sSNV_figure_dir, "/mut_mat_estimiated_conditional_PTA_Cases.csv"))

## age-match: normalize mut_mat to est snv.burden (per cell)
mut_mat <- t(1 / rowSums(mut_mat) * mut_mat * SCAN2_df$snv.rate.per.gb)
mut_mat_age_match <- mut_mat[, Cell_ID_list_age_match]
mut_mat_age_match_conditional <- as.data.frame(mut_mat_age_match) %>% 
  mutate(Control = rowMeans(mut_mat_age_match[, control_range_age_match])) %>% mutate(IHD = rowMeans(mut_mat_age_match[, disease_range_age_match])) %>% 
  mutate(Net_change = IHD - Control) %>% mutate(Net_change = ifelse(Net_change < 0, NA, Net_change))

p_96 <- plot_96_profile_abs(mut_mat_age_match_conditional[, c(control_name, disease_name, "Net_change")]) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.15), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  # theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24), 
  #       axis.title.x = element_text(size = 24, hjust = 0.5), axis.title.y = element_text(size = 24, hjust = 0.5))
ggsave(paste0(main_figure_dir, "/3-sSNV_mutational_profile.pdf"), plot = p_96, width = 10, height = 6, dpi = 600)
write.csv(mut_mat_age_match_conditional[, c(control_name, disease_name, "Net_change")], file = paste0(sSNV_figure_dir, "/mut_mat_age_match_conditional_PTA_Cases.csv"))

# single_signature <- as.matrix(COSMIC_v3.4_sigs[,5])
# colnames(single_signature) <- "SBS5"
# p_COSMIC_SBS5 <- plot_96_profile(single_signature, ymax = 0.05) + 
#   theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.1), panel.grid.minor = element_blank(), 
#         panel.border = element_rect(size = 0.1), text = element_text(size = 4))
# ggsave(paste0(main_figure_dir, "/3-COSMIC_SBS5.pdf"), plot = p_COSMIC_SBS5, width = 6, height = 3, dpi = 600)