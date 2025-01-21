################################################################################
############################# Strand bias test #################################
################################################################################

################################################################################
##### 1. count SNV for all mutations types and strands
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
strand <- mut_strand(grl[[1]], genes_hg19)
mut_mat_s <- mut_matrix_stranded(grl, ref_genome, genes_hg19)

################################################################################
##### 2. reorder, normalize by estimated SNV burden, select age matched samples
mut_mat_s <- mut_mat_s[ , Cell_ID_list]
mut_mat_s <- t(mut_mat_s)
mut_mat_s <- 1 / rowSums(mut_mat_s) * mut_mat_s * SCAN2_df$snv.burden
mut_mat_s <- round(t(mut_mat_s), digits = 0)
mut_mat_s_age_match <- mut_mat_s[ , Cell_ID_list_age_match]

################################################################################
##### 3. strand_counts_age_match analysis by cell
strand_counts_age_match <- strand_occurrences(mut_mat_s_age_match, by = metadata_df_age_match$Condition)
strand_bias_age_match <- strand_bias_test(strand_counts_age_match)

strand_counts_age_match_by_cell <- strand_occurrences(mut_mat_s_age_match[, 1])[,1:3]
for (cell in Cell_ID_list_age_match){
  strand_counts_age_match_by_cell <- strand_counts_age_match_by_cell %>% 
    mutate(!!cell := strand_occurrences(mut_mat_s_age_match[, cell])$no_mutations)
}

melt_strand_counts_age_match_by_cell <- as.data.frame(t(strand_counts_age_match_by_cell[Cell_ID_list_age_match])) %>% 
  mutate(condition = metadata_df_age_match$Condition) %>% 
  setNames(paste0(strand_counts_age_match_by_cell$type, ' (', strand_counts_age_match_by_cell$strand, ")")) %>% 
  melt() %>% 
  setNames(c("condition", "type", "no_mutation")) %>% 
  mutate(condition = factor(condition, level = c(control_name, disease_name)))

barplot_melt_strand_counts_age_match <- ggbarplot(melt_strand_counts_age_match_by_cell, x = "type", y = "no_mutation", color = "condition", label = TRUE, lab.nb.digits = 0, add = c("mean_se", "jitter"), 
                                                  position = position_dodge(0.9), palette = c(dis_ctrl_color[1], dis_ctrl_color[2])) +
  stat_compare_means(aes(group = condition), label = "p.format", label.y = 1.02 * max(melt_strand_counts_age_match_by_cell$no_mutation)) + 
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.04 * max(melt_strand_counts_age_match_by_cell$no_mutation)) +
  labs(x = "mutation type and strand", y = "Absolute Contribution") + scale_y_break(c(3000, 5500)) + ylim(0, 6200) + theme_linedraw() + scale_fill_manual(values = dis_ctrl_color) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.5))
ggsave(paste0(suppl_figure_dir, "/5-strand_counts_age_match.pdf"), plot = barplot_melt_strand_counts_age_match, width = 15, height = 10, dpi = 600)

################################################################################
##### 4. strand ratio between two strands
melt_strand_ratio_age_match_by_cell <- as.data.frame(t(strand_counts_age_match_by_cell[Cell_ID_list_age_match])) %>% 
  mutate("C_A" = V1 / V2) %>% mutate("C>A" = replace(C_A, (is.infinite(C_A) | is.na(C_A)), 0)) %>% mutate("C_G" = V3 / V4) %>% mutate("C>G" = replace(C_G, (is.infinite(C_G) | is.na(C_G)), 0)) %>% 
  mutate("C_T" = V5 / V6) %>% mutate("C>T" = replace(C_T, (is.infinite(C_T) | is.na(C_T)), 0)) %>% mutate("T_A" = V7 / V8) %>% mutate("T>A" = replace(T_A, (is.infinite(T_A) | is.na(T_A)), 0)) %>% 
  mutate("T_C" = V9 / V10) %>% mutate("T>C" = replace(T_C, (is.infinite(T_C) | is.na(T_C)), 0)) %>% mutate("T_G" = V11 / V12) %>% 
  mutate("T>G" = replace(T_G, (is.infinite(T_G) | is.na(T_G)), 0)) |> base::`[`(unique(strand_counts_age_match_by_cell$type)) %>% 
  mutate(condition = metadata_df_age_match$Condition) %>% melt() %>% setNames(c("condition", "type", "ratio")) %>% mutate(condition = factor(condition, level = c(control_name, disease_name)))

barplot_strand_ratio_age_match <- ggbarplot(melt_strand_ratio_age_match_by_cell, x = "type", y = "ratio", color = "condition", label = TRUE, lab.nb.digits = 2, add = c("mean_se", "jitter"), 
                                            position = position_dodge(0.9), palette = c(dis_ctrl_color[1], dis_ctrl_color[2])) +
  stat_compare_means(aes(group = condition), label = "p.format", label.y = 1.02 * max(melt_strand_ratio_age_match_by_cell$ratio, na.rm = TRUE)) + 
  stat_compare_means(aes(group = condition), label = "p.signif", label.y = 1.06 * max(melt_strand_ratio_age_match_by_cell$ratio, na.rm = TRUE)) + 
  labs(x = "mutation type and strand", y = expression(paste("ratio of sSNV (", frac(transcribed, untranscribed), ')'))) + theme_linedraw() + scale_fill_manual(values = dis_ctrl_color) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5))
ggsave(paste0(suppl_figure_dir, "/5-strand_ratio_age_match.pdf"), plot = barplot_strand_ratio_age_match, width = 10, height = 6, dpi = 600)

################################################################################
##### 5. strand_counts_age_match analysis by condition
mut_mat_s_age_match_conditional <- as.data.frame(mut_mat_s_age_match) %>% 
  mutate(Control = round(rowMeans(mut_mat_s_age_match[,control_range_age_match]), digits = 0)) %>% 
  mutate(IHD = round(rowMeans(mut_mat_s_age_match[,disease_range_age_match]), digits = 0)) %>% 
  mutate(Net_change = IHD - Control) %>% mutate(Net_change = ifelse(Net_change < 0, NA, Net_change)) |> base::`[`(c(control_name, disease_name, "Net_change"))

##### age matched by condition all cell together
strand_counts_age_match_conditional <- strand_occurrences(mut_mat_s_age_match_conditional[, c(control_name, disease_name)], by = c(control_name, disease_name))
strand_bias_age_match_conditional <- strand_bias_test(strand_counts_age_match_conditional)

##### summarize 
summary_strand_counts_age_match_conditional <- strand_counts_age_match_conditional %>%
  group_by(group, strand) %>% summarize(mean_value = mean(no_mutations))

################################################################################
##### 6. standard MutationalPattern strand bias analysis
ps1 <- plot_strand(strand_counts_age_match, mode = "absolute") + 
  facet_wrap(factor(group, c(control_name, disease_name)) ~ .) +
  stat_compare_means(aes(label = ..p.signif..))
ggsave(paste0(suppl_figure_dir, "/5-plot_strand_bias_absolute.pdf"), plot = ps1, width = 10, height = 6, dpi = 600)

ps2 <- plot_strand_bias(strand_bias_age_match, sig_type = "p") + 
  facet_wrap(factor(group, c(control_name, disease_name)) ~ .)
ggsave(paste0(suppl_figure_dir, "/5-plot_strand_bias_ratio.pdf"), plot = ps2, width = 10, height = 6, dpi = 600)

strand_bias_notstrict <- strand_bias_test(strand_counts_age_match, p_cutoffs = c(0.5, 0.1, 0.05), fdr_cutoffs = 0.5)
ps3 <- plot_strand_bias(strand_bias_notstrict, sig_type = "p") + facet_wrap(factor(group, c(control_name, disease_name)) ~ .)
ggsave(paste0(suppl_figure_dir, "/5-plot_strand_bias_ratio_test.pdf"), plot = ps3, width = 10, height = 6, dpi = 600)

################################################################################
##### 6. strand bias analysis, condition separately
perform_wilcox_test <- function(mat1, mat2) {
  # wilcox_test_results <- vector("list", length = ncol(mat1))
  wilcox_test_results <- data.frame(matrix(NA, nrow = ncol(mat1), ncol = 1))
  for (i in 1 : ncol(mat1)) {
    wilcox_test_results[i,1] <- wilcox.test(mat1[, i], mat2[, i], alternative = "two.sided")$p.value
  }
  return(wilcox_test_results)
}

stars.pval <- function(x){
  # stars <- c("***", "**", "*", "ns")
  stars <- c("***", "**", "*", " ")
  var <- c(0, 0.001, 0.01, 0.05, 0.10, 1)
  i <- findInterval(x, var, left.open = T, rightmost.closed = T)
  stars[i]
}

strand_counts_02 <- as.data.frame(t(strand_counts_age_match_by_cell[Cell_ID_list_age_match])) %>% setNames(strand_counts_age_match_by_cell$strand)

wilcox_test_control <- perform_wilcox_test(strand_counts_02[control_range_age_match, colnames(strand_counts_02) == "transcribed"], 
                                           strand_counts_02[control_range_age_match, colnames(strand_counts_02) == "untranscribed"]) %>% 
  setNames("p_value") %>% mutate(p_value = round(p_value, 3)) %>% mutate(sig_level = stars.pval(p_value))
  
wilcox_test_disease <- perform_wilcox_test(strand_counts_02[disease_range_age_match, colnames(strand_counts_02) == "transcribed"], 
                                           strand_counts_02[disease_range_age_match, colnames(strand_counts_02) == "untranscribed"]) %>% 
  setNames("p_value") %>% mutate(p_value = round(p_value, 3)) %>% mutate(sig_level = stars.pval(p_value))

strand_counts_age_match_conditional$label <- c(rep("",24))
strand_counts_age_match_conditional$label[10] <- wilcox_test_control$sig_level[5]

p_strand_bias_analysis_condition_separately <- plot_strand(strand_counts_age_match_conditional, mode = "absolute") + 
  geom_text(aes(x = 4.6, y = 215, label = strand_counts_age_match_conditional$label), vjust = 0, hjust = 0, color = "black", size = 8, angle = 0, show.legend = FALSE, alpha = 1) + 
  labs(y = "sSNV rate per cell") + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5))
ggsave(paste0(main_figure_dir, "/5-strand_bias_analysis_condition_separately.pdf"), plot = p_strand_bias_analysis_condition_separately, width = 8, height = 5, dpi = 600)

ps4 <- plot_strand_bias(strand_bias_age_match_conditional, sig_type = "p") + 
  facet_wrap(factor(group, c(control_name, disease_name)) ~ .) + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white"))
print(ps4)
