############################################################################# 
######################### sSNV burden visualization #########################
############################################################################# 

p_SNV_burden_per_GB <- ggplot(sSNV_SCAN2_df, aes(x = factor(Cell_ID, level = Cell_ID_list), y = snv.rate.per.gb, color = Condition)) + 
  geom_point(size = 5, na.rm = FALSE, show.legend = FALSE) + scale_color_manual(values = c("Normoxia" = ctrl_dis_color[1], "Hypoxia" = ctrl_dis_color[2]), guide = "legend") + 
  labs(x = "Sample ID", y = "sSNV rate per GB") + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", size = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.5), text = element_text(size = 24), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + 
  facet_grid(. ~  factor(Condition, level = c(ctrl_name, dis_name)), scale = "free_x")
ggsave(paste0(other_figure_dir, "/2-sSNV_burden_per_GB.pdf"), plot = p_SNV_burden_per_GB, width = 9, height = 6, dpi = 600)

wilcox.test_result <- wilcox.test(sSNV_SCAN2_df$snv.rate.per.gb[ctrl_range], sSNV_SCAN2_df$snv.rate.per.gb[dis_range], alternative = "two.sided")
sSNV_SCAN2_df$label <- c(rep("Normoxia (n = 9)", 9), rep("Hypoxia (n = 9)", 9))

boxplot_sSNV_burden <- ggplot(sSNV_SCAN2_df, aes(x = factor(label, level = c("Normoxia (n = 9)", "Hypoxia (n = 9)")), y = snv.rate.per.gb, fill = label)) + 
  geom_boxplot(color = c("dodgerblue3", "firebrick3"), fill = c("#C7DFF0", "#FBDFE2"), na.rm = TRUE, outlier.shape = NA) + 
  geom_jitter(pch = 21, fill = sSNV_SCAN2_df$Color, color = sSNV_SCAN2_df$Outline_Color, size = 5, na.rm = TRUE) + 
  annotate("text", size = 7, x = 1.32, y = 220, label = paste("P = ", format(wilcox.test_result$p.value, scientific = FALSE, digits = 2)), hjust = 0) + 
  labs(x = "", y = "sSNV rate per GB") + stat_summary(fun = mean, colour = "black", geom = "point", shape = 18, size = 3, show.legend = FALSE) + 
  stat_compare_means(comparisons = list(c("Normoxia (n = 9)", "Hypoxia (n = 9)")), label.x = 1.6, label.y = 220,
                     label = "p.signif", bracket.size = 1, size = 7, method = "wilcox.test") + theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), 
        text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
ggsave(paste0(main_figure_dir, "/2-boxplot_sSNV_burden.pdf"), plot = boxplot_sSNV_burden, width = 8, height = 8)
