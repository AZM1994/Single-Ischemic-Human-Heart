################################################################################
############################# Strand bias test #################################
################################################################################

################################################################################
##### 1. count SNV for all mutations types and strands
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
strand <- mut_strand(snv_grl[[1]], genes_hg19)
mut_mat_s <- mut_matrix_stranded(snv_grl, ref_genome, genes_hg19)

################################################################################
##### 2. reorder, normalize by estimated SNV burden, select age matched samples
mut_mat_s <- mut_mat_s[, Cell_ID_list] %>% t() %>% { . / rowSums(.) * SCAN2_df$snv.burden } %>% t() %>% round(digits = 0)
mut_mat_s_AMG <- mut_mat_s[, Cell_ID_list_AMG]
# mut_mat_s_AMG <- mut_mat_s[ , Cell_ID_list_AMG]

################################################################################
##### 3. strand_counts_AMG analysis by cell
strand_counts_AMG <- strand_occurrences(mut_mat_s_AMG, by = metadata_df_AMG$Condition)
strand_bias_AMG <- strand_bias_test(strand_counts_AMG)

strand_counts_AMG_by_cell <- strand_occurrences(mut_mat_s_AMG[, 1])[,1:3]
for (cell in Cell_ID_list_AMG){
  strand_counts_AMG_by_cell <- strand_counts_AMG_by_cell %>% mutate(!!cell := strand_occurrences(mut_mat_s_AMG[, cell])$no_mutations)
}

melt_strand_counts_AMG_by_cell <- as.data.frame(t(strand_counts_AMG_by_cell[Cell_ID_list_AMG])) %>% 
  mutate(Condition = metadata_df_AMG$Condition) %>% setNames(paste0(strand_counts_AMG_by_cell$type, ' (', strand_counts_AMG_by_cell$strand, ")")) %>% 
  melt() %>% setNames(c("Condition", "type", "no_mutation")) %>% mutate(Condition = factor(Condition, level = c(ctrl_name, dis_name)))

sc_yticks <- generate_y_breaks(melt_strand_counts_AMG_by_cell$no_mutation)
sc_ylim <- range(sc_yticks, na.rm = TRUE)
barplot_sc_AMG <- ggbarplot(melt_strand_counts_AMG_by_cell, x = "type", y = "no_mutation", color = "Condition", label = TRUE, lab.nb.digits = 0, 
                            add = c("mean_se", "jitter"), position = position_dodge(0.9), palette = c(ctrl_dis_color[1], ctrl_dis_color[2])) + 
  stat_compare_means(aes(group = Condition, label = sprintf("p = %1.2e", as.numeric(..p.format..))), label.y = 1.08 * sqrt(max(melt_strand_counts_AMG_by_cell$no_mutation))) + 
  stat_compare_means(aes(group = Condition), label = "p.signif", label.y = 1.05 * sqrt(max(melt_strand_counts_AMG_by_cell$no_mutation))) +
  labs(x = "mutation type and strand", y = "Absolute Contribution") + scale_y_continuous(limits = sc_ylim, breaks = sc_yticks, labels = sc_yticks, trans = "sqrt") + 
  theme_linedraw() + scale_fill_manual(values = ctrl_dis_color) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0), axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0))
ggsave(paste0(other_figure_dir, "/5-strand_counts_AMG.pdf"), plot = barplot_sc_AMG, width = 15, height = 10, dpi = 600)

################################################################################
##### 4. strand ratio between two strands
melt_strand_ratio_AMG_by_cell <- as.data.frame(t(strand_counts_AMG_by_cell[Cell_ID_list_AMG])) %>% 
  mutate("C_A" = V1 / V2) %>% mutate("C>A" = replace(C_A, (is.infinite(C_A) | is.na(C_A)), 0)) %>% mutate("C_G" = V3 / V4) %>% mutate("C>G" = replace(C_G, (is.infinite(C_G) | is.na(C_G)), 0)) %>% 
  mutate("C_T" = V5 / V6) %>% mutate("C>T" = replace(C_T, (is.infinite(C_T) | is.na(C_T)), 0)) %>% mutate("T_A" = V7 / V8) %>% mutate("T>A" = replace(T_A, (is.infinite(T_A) | is.na(T_A)), 0)) %>% 
  mutate("T_C" = V9 / V10) %>% mutate("T>C" = replace(T_C, (is.infinite(T_C) | is.na(T_C)), 0)) %>% mutate("T_G" = V11 / V12) %>% 
  mutate("T>G" = replace(T_G, (is.infinite(T_G) | is.na(T_G)), 0)) |> base::`[`(unique(strand_counts_AMG_by_cell$type)) %>% 
  mutate(Condition = metadata_df_AMG$Condition) %>% melt() %>% setNames(c("Condition", "type", "ratio")) %>% mutate(Condition = factor(Condition, level = c(ctrl_name, dis_name)))

barplot_sr_AMG <- ggbarplot(melt_strand_ratio_AMG_by_cell, x = "type", y = "ratio", color = "Condition", label = TRUE, lab.nb.digits = 2, 
                            add = c("mean_se", "jitter"), position = position_dodge(0.9), palette = c(ctrl_dis_color[1], ctrl_dis_color[2])) +
  stat_compare_means(aes(group = Condition), label = "p.format", label.y = 1.02 * max(melt_strand_ratio_AMG_by_cell$ratio, na.rm = TRUE)) + 
  stat_compare_means(aes(group = Condition), label = "p.signif", label.y = 1.06 * max(melt_strand_ratio_AMG_by_cell$ratio, na.rm = TRUE)) + 
  labs(x = "mutation type and strand", y = expression(paste("ratio of sSNV (", frac(transcribed, untranscribed), ')'))) + theme_linedraw() + scale_fill_manual(values = ctrl_dis_color) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5))
ggsave(paste0(other_figure_dir, "/5-strand_ratio_AMG.pdf"), plot = barplot_sr_AMG, width = 10, height = 6, dpi = 600)

################################################################################
##### 5. strand_counts_AMG analysis by Condition
mut_mat_s_AMG_cond <- as.data.frame(mut_mat_s_AMG) %>% 
  mutate(Control = round(rowMeans(mut_mat_s_AMG[, ctrl_range_AMG]), digits = 0)) %>% 
  mutate(IHD = round(rowMeans(mut_mat_s_AMG[, dis_range_AMG]), digits = 0)) %>% 
  mutate(Net_change = IHD - Control) %>% mutate(Net_change = ifelse(Net_change < 0, NA, Net_change)) |> base::`[`(c(ctrl_name, dis_name, "Net_change"))

##### age matched by Condition all cell together
strand_counts_AMG_cond <- strand_occurrences(mut_mat_s_AMG_cond[, c(ctrl_name, dis_name)], by = c(ctrl_name, dis_name))
strand_bias_AMG_cond <- strand_bias_test(strand_counts_AMG_cond)
##### summarize 
summary_strand_counts_AMG_cond <- strand_counts_AMG_cond %>% group_by(group, strand) %>% summarize(mean_value = mean(no_mutations))

################################################################################
##### 6. standard MutationalPattern strand bias analysis
ps1 <- plot_strand(strand_counts_AMG, mode = "absolute") + facet_wrap(factor(group, c(ctrl_name, dis_name)) ~ .) + stat_compare_means(aes(label = ..p.signif..))
ggsave(paste0(other_figure_dir, "/5-strand_bias_absolute.pdf"), plot = ps1, width = 10, height = 6, dpi = 600)

ps2 <- plot_strand_bias(strand_bias_AMG, sig_type = "p") + facet_wrap(factor(group, c(ctrl_name, dis_name)) ~ .)
ggsave(paste0(other_figure_dir, "/5-strand_bias_ratio.pdf"), plot = ps2, width = 10, height = 6, dpi = 600)

strand_bias_notstrict <- strand_bias_test(strand_counts_AMG, p_cutoffs = c(0.5, 0.1, 0.05), fdr_cutoffs = 0.5)
ps3 <- plot_strand_bias(strand_bias_notstrict, sig_type = "p") + facet_wrap(factor(group, c(ctrl_name, dis_name)) ~ .)
ggsave(paste0(other_figure_dir, "/5-strand_bias_ratio_test.pdf"), plot = ps3, width = 10, height = 6, dpi = 600)

################################################################################
##### 6. strand bias analysis, Condition separately
perform_wilcox_test <- function(mat1, mat2) {
  # wilcox_test_results <- vector("list", length = ncol(mat1))
  wilcox_test_results <- data.frame(matrix(NA, nrow = ncol(mat1), ncol = 1))
  for (i in 1 : ncol(mat1)) {
    wilcox_test_results[i,1] <- wilcox.test(mat1[, i], mat2[, i], alternative = "two.sided", exact = FALSE)$p.value
  }
  return(wilcox_test_results)
}

stars.pval <- function(x){
  # stars <- c("***", "**", "*", "ns")
  stars <- c("***", "**", "*", "NA")
  var <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, var, left.open = T, rightmost.closed = T)
  stars[i]
}

strand_counts_02 <- as.data.frame(t(strand_counts_AMG_by_cell[Cell_ID_list_AMG])) %>% setNames(strand_counts_AMG_by_cell$strand)

wilcox_test_control <- perform_wilcox_test(strand_counts_02[ctrl_range_AMG, colnames(strand_counts_02) == "transcribed"], 
                                           strand_counts_02[ctrl_range_AMG, colnames(strand_counts_02) == "untranscribed"]) %>% 
  setNames("p_value") %>% mutate(p_value = round(p_value, 3)) %>% mutate(sig_level = stars.pval(p_value))
  
wilcox_test_disease <- perform_wilcox_test(strand_counts_02[dis_range_AMG, colnames(strand_counts_02) == "transcribed"], 
                                           strand_counts_02[dis_range_AMG, colnames(strand_counts_02) == "untranscribed"]) %>% 
  setNames("p_value") %>% mutate(p_value = round(p_value, 3)) %>% mutate(sig_level = stars.pval(p_value))

strand_counts_AMG_cond$label <- c(rep("",24))
strand_counts_AMG_cond$label[22] <- wilcox_test_disease$sig_level[5]

p_strand_bias_cond <- plot_strand(strand_counts_AMG_cond, mode = "absolute") + geom_text(aes(x = 5, y = 630, label = strand_counts_AMG_cond$label), color = "black", size = 8, show.legend = FALSE, alpha = 1) + 
  labs(y = "sSNV rate per cell") + theme_linedraw() + theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
                                                            panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.4, vjust = 0), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5))
ggsave(paste0(main_figure_dir, "/5-strand_bias_cond.pdf"), plot = p_strand_bias_cond, width = 8, height = 5, dpi = 600)

ps4 <- plot_strand_bias(strand_bias_AMG_cond, sig_type = "p") + facet_wrap(factor(group, c(ctrl_name, dis_name)) ~ .) + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), panel.background = element_rect(fill = "white"))
