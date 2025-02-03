##### Signature refitting strict
# strict_refit <- fit_to_signatures_strict(indel_counts_02, indel_signatures, max_delta = 0.21)
strict_refit <- fit_to_signatures_strict(indel_counts, indel_signatures, max_delta = 0.21)
# strict_refit <- fit_to_signatures_strict(mut_mat, signatures[,1], max_delta = 0.04, method = "best_subset")

########## find selected signatures
fit_res_strict <- strict_refit$fit_res
strict_sig_names <- unique(names(which(fit_res_strict$contribution != 0, arr.ind = T)[,1]))

########## put all the rest to others
new_fit_res_conditional <- fit_res_strict$contribution
Control <- rowSums(new_fit_res_conditional[,control_range_age_match]) / length(control_range_age_match)
new_fit_res_conditional <- cbind(new_fit_res_conditional, Control)
Disease <- rowSums(new_fit_res_conditional[,disease_range_age_match]) / length(disease_range_age_match)
new_fit_res_conditional <- cbind(new_fit_res_conditional, Disease)
sum_range <- (ncol(new_fit_res_conditional)-1):ncol(new_fit_res_conditional)
colnames(new_fit_res_conditional)[sum_range] <- c(control_name, disease_name)

pdf(paste0(sindel_plot_dir, "/9-indel_Signature_refitting_strict.pdf"), width = 16, height = 12)
print(plot_contribution(new_fit_res_conditional[,-sum_range], coord_flip = FALSE, mode = "absolute") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)))

print(plot_original_vs_reconstructed(indel_counts_02, fit_res_strict$reconstructed, y_intercept = 0.95) +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)) + 
  geom_hline(aes(yintercept = 0.95), color="red", size = 1) + 
  geom_hline(aes(yintercept = 0.85), color="green", size = 1))

print(plot_contribution_heatmap(new_fit_res_conditional, cluster_sigs = FALSE, cluster_samples = FALSE) + 
  scale_fill_gradient(low = "white", high = "blue") + geom_tile(colour = "black") + 
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)))

sorted_labels <- str_sort(strict_sig_names, numeric = TRUE)
sorted_labels <- rbind(sorted_labels, sorted_labels)
sorted_labels <- Reshape(sorted_labels,ncol(sorted_labels)*2,1)
new_fit_res_conditional <- new_fit_res_conditional[rowSums(new_fit_res_conditional[])>0,]
label_size <- 8 * rbind(new_fit_res_conditional[,sum_range[1]] / new_fit_res_conditional[,sum_range[1]], new_fit_res_conditional[,sum_range[2]] / new_fit_res_conditional[,sum_range[2]])
label_size <- Reshape(label_size,ncol(label_size)*2,1)
label_size[is.na(label_size)] <- 0

print(plot_contribution(new_fit_res_conditional[,sum_range], coord_flip = FALSE, mode = "absolute") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)) + 
  aes(label = sorted_labels) +
  geom_text(size = label_size, position = position_stack(vjust = 0.5)))

print(plot_contribution(new_fit_res_conditional[,sum_range], coord_flip = FALSE, mode = "relative") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)))

fit_res_strict_recon_integrated <- cbind(fit_res_strict$reconstructed, 
                                         rowSums(fit_res_strict$reconstructed[,control_range_age_match]), rowSums(fit_res_strict$reconstructed[,disease_range_age_match]))
colnames(fit_res_strict_recon_integrated)[sum_range] <- c(control_name, disease_name)
# print(plot_original_vs_reconstructed(indel_counts_integrated[,sum_range], fit_res_strict_recon_integrated[,sum_range], y_intercept = 0.95) +
#   theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)) + 
#   geom_hline(aes(yintercept = 0.95), color="red", size = 1))

print(plot_contribution_heatmap(new_fit_res_conditional[,sum_range], cluster_sigs = FALSE, cluster_samples = FALSE) + 
  scale_fill_gradient(low = "white", high = "blue") + geom_tile(colour = "black") + 
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)))

dev.off()