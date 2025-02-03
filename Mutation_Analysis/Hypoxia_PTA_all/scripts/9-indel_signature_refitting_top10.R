##### Signature refitting selected list of signatures
########## refined selected signatures
# strict_refit_02 <- fit_to_signatures_strict(indel_counts_02, indel_signatures, max_delta = 0.1)
strict_refit_02 <- fit_to_signatures_strict(indel_counts, indel_signatures, max_delta = 0.1)
fit_res_strict_02 <- strict_refit_02$fit_res
strict_sig_names <- unique(names(which(fit_res_strict_02$contribution != 0, arr.ind = T)[,1]))
# refined_selected_sigs_list <- strict_sig_names[c(1:9,11)]
refined_selected_sigs_list <- strict_sig_names
refined_selected_sigs_list
########## put all the rest to others
# fit_res_strict_conditional_02 <- fit_res_strict$contribution
# fit_res_strict_conditional_02 <- fit_res_loose$contribution
fit_res_strict_conditional_02_all <- fit_res_strict$contribution # strict refitting, all cell/age 60*33
fit_res_strict_conditional_02 <- fit_res_strict_conditional_02_all[,Cell_ID_list_age_match] # 60*22

refined_selected_sigs <- fit_res_strict_conditional_02[(row.names(fit_res_strict_conditional_02) %in% refined_selected_sigs_list),]
others <- colSums(fit_res_strict_conditional_02[!(row.names(fit_res_strict_conditional_02) %in% refined_selected_sigs_list),])
refined_new_fit_res <- rbind(refined_selected_sigs, others)

########## conditional new
# refined_new_fit_res_conditional <- refined_new_fit_res
refined_new_fit_res_conditional <- t(refined_new_fit_res)
refined_new_fit_res_conditional <- t(refined_new_fit_res_conditional / rowSums(refined_new_fit_res_conditional) * SCAN2_df_age_match$indel.burden)

Control <- rowSums(refined_new_fit_res_conditional[,control_range_age_match]) / length(control_range_age_match)
refined_new_fit_res_conditional <- cbind(refined_new_fit_res_conditional, Control)
Disease <- rowSums(refined_new_fit_res_conditional[,disease_range_age_match]) / length(disease_range_age_match)
refined_new_fit_res_conditional <- cbind(refined_new_fit_res_conditional, Disease)
sum_range <- (ncol(refined_new_fit_res_conditional)-1):ncol(refined_new_fit_res_conditional)
colnames(refined_new_fit_res_conditional)[sum_range] <- c(control_name, disease_name)

for (name_temp in refined_selected_sigs_list) {
  control_temp <- refined_new_fit_res_conditional[name_temp,control_range_age_match]
  disease_temp <- refined_new_fit_res_conditional[name_temp,disease_range_age_match]
  print(c(name_temp, round(t.test(control_temp, disease_temp)$p.value,2)))
}
refined_selected_sigs_list

sorted_labels <- str_sort(refined_selected_sigs_list, numeric = TRUE)
sorted_labels <- c(sorted_labels, "others")
sorted_labels <- rbind(sorted_labels, sorted_labels)
sorted_labels <- Reshape(sorted_labels,ncol(sorted_labels)*2,1)
new_fit_res_conditional <- new_fit_res_conditional[rowSums(new_fit_res_conditional[])>0,]

pdf(paste0(sindel_plot_dir, "/9-indel_Signature_refitting_strict_refined.pdf"), width = 10, height = 6)
print(plot_contribution(refined_new_fit_res[1:nrow(refined_new_fit_res_conditional),], coord_flip = FALSE, mode = "absolute") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)))

print(plot_contribution_heatmap(refined_new_fit_res[1:nrow(refined_new_fit_res_conditional),], cluster_sigs = FALSE, cluster_samples = FALSE) + 
  scale_fill_gradient(low = "white", high = "blue") + geom_tile(colour = "black") + 
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)))


refined_new_fit_res_conditional_matrix <- refined_new_fit_res_conditional[1:nrow(refined_new_fit_res_conditional),sum_range]
label_size <- 5 * rbind(refined_new_fit_res_conditional_matrix[,1] / refined_new_fit_res_conditional_matrix[,1], refined_new_fit_res_conditional_matrix[,2] / refined_new_fit_res_conditional_matrix[,2])
label_size <- Reshape(label_size,ncol(label_size)*2,1)
label_size[is.na(label_size)] <- 0
print(plot_contribution(refined_new_fit_res_conditional_matrix, coord_flip = FALSE, mode = "absolute") +
        theme(text = element_text(size=20), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5)) + 
        aes(label = sorted_labels) +
        geom_text(size = label_size, position = position_stack(vjust = 0.5)))

print(plot_contribution(refined_new_fit_res_conditional_matrix, coord_flip = FALSE, mode = "relative") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5)))

# a <- refined_new_fit_res_conditional_matrix/rowSums(t(refined_new_fit_res_conditional_matrix))
refined_new_fit_res_conditional_matrix[,1] <- refined_new_fit_res_conditional_matrix[,1]/rowSums(t(refined_new_fit_res_conditional_matrix))[1]
refined_new_fit_res_conditional_matrix[,2] <- refined_new_fit_res_conditional_matrix[,2]/rowSums(t(refined_new_fit_res_conditional_matrix))[2]
# round(Reshape(refined_new_fit_res_conditional_matrix,nrow(refined_new_fit_res_conditional_matrix),1), 2)
print(plot_contribution_heatmap(refined_new_fit_res_conditional_matrix, cluster_sigs = FALSE, cluster_samples = FALSE) + 
  scale_fill_gradient(low = "white", high = "blue") + geom_tile(colour = "black") + 
  geom_text(aes(label = round(Reshape(refined_new_fit_res_conditional_matrix,nrow(refined_new_fit_res_conditional_matrix)*2,1), 2)), size = 6) + 
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5)))

dev.off()