##### 96 mutational profile
mut_mat <- mut_matrix(vcf_list = snv_grl, ref_genome = ref_genome)
mut_mat <- t(mut_mat[, Cell_ID_list])
mut_num_summary <- as.data.frame(SCAN2_df) %>% dplyr::select(Cell_ID, Condition, snv.burden) %>% mutate(true_call = rowSums(mut_mat))

##### generate mut_mat for SigNet (raw)
mut_mat_raw <- mut_mat
write.csv(mut_mat_raw, file = paste0(table_dir, "/mut_mat_raw_PTA_all.csv"))
# mut_mat_raw_cond <- mut_mat_raw %>% t() %>% as.data.frame() %>% mutate(Control = rowSums(.[, control_range]), IHD = rowSums(.[, disease_range])) %>% dplyr::select(Control, IHD) %>% t()
# write.csv(mut_mat_raw_cond, file = paste0(table_dir, "/mut_mat_raw_conditional_PTA_all.csv"))

##### generate mut_mat for SigNet (normalize mut_mat to est snv.burden (per cell))
mut_mat_est <- 1 / rowSums(mut_mat) * mut_mat * SCAN2_df$snv.rate.per.gb
write.csv(mut_mat_est, file = paste0(table_dir, "/mut_mat_est_PTA_all.csv"))
# mut_mat_est_cond <- mut_mat_est %>% t() %>% as.data.frame() %>% mutate(Control = rowSums(.[, control_range]), IHD = rowSums(.[, disease_range])) %>% dplyr::select(Control, IHD) %>% t()
# write.csv(mut_mat_est_cond, file = paste0(table_dir, "/mut_mat_est_conditional_PTA_all.csv"))
################################################################################
## Function to filter out outlier Cell_ID from given mut_mat matrix
filter_outliers <- function(data, threshold = 3.5) {
  dists <- as.matrix(dist(data)) # Compute Euclidean distances between samples
  avg_dists <- apply(dists, 1, function(x) mean(x[x != 0])) # For each sample, calculate its average distance to all other samples
  med <- median(avg_dists) # Calculate median and standard deviation of the average distances
  sd_val <- sd(avg_dists)
  outliers <- names(avg_dists)[avg_dists > med + threshold * sd_val] # Identify outlier samples: those whose average distance exceeds (median + threshold*SD)
  cat("outliers:", outliers)
  setdiff(rownames(data), outliers) # Return the names of samples that are not outliers
}

## Filter out outliers in each condition
filtered_mut_mat_est_Control <- filter_outliers(mut_mat_est[SCAN2_df$Cell_ID[SCAN2_df$Condition == "Control"], ], threshold = 3.5)
filtered_mut_mat_est_IHD <- filter_outliers(mut_mat_est[SCAN2_df$Cell_ID[SCAN2_df$Condition == "IHD"], ], threshold = 3.5)
mut_mat_est_filtered <- mut_mat_est[c(filtered_mut_mat_est_Control, filtered_mut_mat_est_IHD), ]
write.csv(mut_mat_est_filtered, file = paste0(table_dir, "/mut_mat_est_PTA_all_filtered.csv"))
mut_mat_raw_filtered <- mut_mat_raw[c(filtered_mut_mat_est_Control, filtered_mut_mat_est_IHD), ]
write.csv(mut_mat_raw_filtered, file = paste0(table_dir, "/mut_mat_raw_PTA_all_filtered.csv"))
filtered_samples <- c(filtered_mut_mat_est_Control, filtered_mut_mat_est_IHD)

dists_all <- dist(mut_mat_est, method = "euclidean")
mds_fit <- cmdscale(dists_all, k = 2)  # Classical multidimensional scaling in two dimensions
mds_df <- data.frame(MDS1 = mds_fit[,1], MDS2 = mds_fit[,2], Cell_ID = rownames(mut_mat_est), Condition = SCAN2_df$Condition[match(rownames(mut_mat_est), SCAN2_df$Cell_ID)])
mds_df <- mds_df %>% mutate(Filtered = ifelse(mds_df$Cell_ID %in% filtered_samples, "Kept", "Filtered Out")) %>% mutate(Filtered = factor(Filtered, level = c("Kept", "Filtered Out")))
p_mds <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = Condition, shape = Filtered)) + geom_point(size = 3) + labs(x = "Dimension 1", y = "Dimension 2") + theme_linedraw() + 
  scale_color_manual(values = c("Control" = ctrl_dis_color[1], "IHD" = ctrl_dis_color[2]), guide = "legend") + scale_shape_manual(values = c(16, 17)) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
ggsave(paste0(suppl_figure_dir, "/3-MDS_mut_mat_est.pdf"), plot = p_mds, width = 8, height = 5.5, dpi = 600)

dist_mat <- as.matrix(dists_all)
annotation <- data.frame(Condition = SCAN2_df$Condition[match(rownames(mut_mat_est), SCAN2_df$Cell_ID)])
rownames(annotation) <- rownames(dist_mat)
p_heatmap_dist <- pheatmap(dist_mat, annotation_row = annotation, main = "Heatmap of Euclidean Distances Between Samples", cluster_rows = F, cluster_cols = F)
ggsave(paste0(suppl_figure_dir, "/3-Heatmap_dist_mut_mat_est.pdf"), plot = p_heatmap_dist, width = 15, height = 12, dpi = 600)

################################################################################

## age-match: normalize mut_mat to est snv.burden (per cell)
mut_mat_age_match <- t(1 / rowSums(mut_mat) * mut_mat * SCAN2_df$snv.rate.per.gb)
mut_mat_age_match <- mut_mat_age_match[, Cell_ID_list_age_match]
mut_mat_age_match_conditional <- as.data.frame(mut_mat_age_match) %>% 
  mutate(Control = rowMeans(mut_mat_age_match[, control_range_age_match])) %>% mutate(IHD = rowMeans(mut_mat_age_match[, disease_range_age_match])) %>% 
  mutate(Net_change = IHD - Control) %>% mutate(Net_change = ifelse(Net_change < 0, NA, Net_change))

p_96 <- plot_96_profile_abs(mut_mat_age_match_conditional[, c(control_name, disease_name, "Net_change")]) + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.15), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(main_figure_dir, "/3-sSNV_mutational_profile.pdf"), plot = p_96, width = 10, height = 6, dpi = 600)
write.csv(mut_mat_age_match_conditional[, c(control_name, disease_name, "Net_change")], file = paste0(table_dir, "/mut_mat_age_match_conditional_PTA_all.csv"))

# single_signature <- as.matrix(COSMIC_v3.4_sigs[,5])
# colnames(single_signature) <- "SBS5"
# p_COSMIC_SBS5 <- plot_96_profile(single_signature, ymax = 0.05) + 
#   theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.1), panel.grid.minor = element_blank(), 
#         panel.border = element_rect(size = 0.1), text = element_text(size = 4))
# ggsave(paste0(main_figure_dir, "/3-COSMIC_SBS5.pdf"), plot = p_COSMIC_SBS5, width = 6, height = 3, dpi = 600)