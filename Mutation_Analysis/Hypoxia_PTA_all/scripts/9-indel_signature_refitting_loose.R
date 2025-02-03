##### indel signature_refitting_loose
# indel_counts_02 <- indel_counts_est[,c(control_range_age_match, disease_range_age_match)]
# indel_counts_02 <- indel_counts[ ,sample_name_list_age_match]
indel_counts_02 <- indel_counts
indel_signatures = get_known_signatures(muttype = c("indel"))
fit_res_loose <- fit_to_signatures(indel_counts_02, indel_signatures)
loose_sig_names <- unique(names(which(fit_res_loose$contribution != 0, arr.ind = T)[,1]))

loose_contrilbution <- fit_res_loose$contribution
loose_contrilbution <- loose_contrilbution[,Cell_ID_list_age_match]

##### the contribution bar plot
pdf(paste0(sindel_plot_dir, "/9-indel_Signature_refitting_loose.pdf"), width = 16, height = 12)
print(plot_contribution(loose_contrilbution, coord_flip = FALSE, mode = "absolute") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)) +
  guides(fill=guide_legend(ncol=2)))
##### plot the cosine similarity between the original and loose reconstructed
print(plot_original_vs_reconstructed(indel_counts_02, fit_res_loose$reconstructed, y_intercept = 0.95) +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)) + 
  geom_hline(aes(yintercept = 0.95), color="red", size = 1))

print(plot_contribution_heatmap(loose_contrilbution, cluster_sigs = FALSE, cluster_samples = FALSE) + 
  scale_fill_gradient(low = "white", high = "blue") + geom_tile(colour = "black")  + 
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)))

############## by conditions
fit_res_loose_conditional <- loose_contrilbution
Control <- rowSums(fit_res_loose_conditional[,control_range_age_match])
fit_res_loose_conditional <- cbind(fit_res_loose_conditional, Control)
Disease <- rowSums(fit_res_loose_conditional[,disease_range_age_match])
fit_res_loose_conditional <- cbind(fit_res_loose_conditional, Disease)
sum_range <- (ncol(fit_res_loose_conditional)-1):ncol(fit_res_loose_conditional)
colnames(fit_res_loose_conditional)[sum_range] <- c(control_name, disease_name)

print(plot_contribution(fit_res_loose_conditional[,sum_range], coord_flip = FALSE, mode = "absolute") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)) + 
  guides(fill=guide_legend(ncol=2)))

indel_counts_integrated <- cbind(indel_counts_02, rowSums(indel_counts_02[,control_range_age_match]), rowSums(indel_counts_02[,disease_range_age_match]))
colnames(indel_counts_integrated)[sum_range] <- c(control_name, disease_name)
fit_res_loose_recon_integrated <- cbind(fit_res_loose$reconstructed, 
                                        rowSums(fit_res_loose$reconstructed[,control_range_age_match]), rowSums(fit_res_loose$reconstructed[,disease_range_age_match]))
colnames(fit_res_loose_recon_integrated)[sum_range] <- c(control_name, disease_name)
print(plot_original_vs_reconstructed(indel_counts_integrated[,sum_range], fit_res_loose_recon_integrated[,sum_range], y_intercept = 0.95) +
  scale_x_discrete(limits = c(control_name, disease_name)) + 
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)) + 
  geom_hline(aes(yintercept = 0.95), color="red", size = 1))

print(plot_contribution_heatmap(fit_res_loose_conditional[,sum_range], cluster_sigs = FALSE, cluster_samples = FALSE) + 
  scale_fill_gradient(low = "white", high = "blue") + geom_tile(colour = "black") + 
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)))

dev.off()
