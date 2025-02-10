################################################################################
############################### circos plot ####################################
################################################################################
Condition_list <- c("Control", "IHD")
genomic_context_colnames <- c("Cell_ID", "Case_ID", "Condition", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene")
genomic_context <- c()
for (condition_temp in Condition_list){
  heart_PTA_Cases_vcf_temp <- read.table(paste0(project_dir, "/data/vcfs/heart_PTA_all.all_age.", condition_temp, "_ssnv.vcf"), sep = "\t")
  genomic_context_temp <- read.csv(paste0(project_dir, "/data/vcfs/heart_PTA_all.all_age.", condition_temp, ".hg19_multianno.csv"), header = TRUE) %>%
    mutate(Cell_ID = heart_PTA_Cases_vcf_temp$V8) %>% mutate(Case_ID = str_extract(Cell_ID, "[^_]+")) %>% 
    mutate(Condition = condition_temp) |> base::`[`(genomic_context_colnames)
  
  genomic_context_temp <- genomic_context_temp |> base::`[`(c("Chr", "Start", "Cell_ID", "Case_ID", "Condition")) %>% 
    rename_with(~ c("Chr", "Postion", "Cell_ID", "Case_ID", "Condition"))
  
  genomic_context <- rbind(genomic_context, genomic_context_temp)
}

genomic_context <- genomic_context %>% filter(Cell_ID %in% Cell_ID_list)
case_list <- unique(genomic_context$Case_ID)
condition_list <- rev(unique(genomic_context$Condition))
control_color_palette <- colorRampPalette(c("skyblue","dodgerblue4"))
disease_color_palette <- colorRampPalette(c("pink","firebrick"))
color_list <- c(control_color_palette(4), disease_color_palette(5))

################################################################################
##### 1. by donor all together
# pdf(paste0(sSNV_figure_dir, "/6-ssnv_circos_plot_by_case.pdf"), width = 10, height = 6)
# circos.par(start.degree=90)
# circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr",1:22), 
#                               plotType = c("ideogram", "labels"), ideogram.height = 0.04)
# for(i in 1:length(case_list)){
#   df_site = genomic_context[genomic_context$Case_ID == case_list[i], 1:2]
#   df2_site = data.frame(Chr = sprintf("chr%s", df_site[, 1]), Start = df_site[, 2], End = df_site[, 2], Value = 1)
#   df3_site = df2_site[!df2_site[, 1] %in% c("chrX", "chrY"), ]
#   # circos.genomicTrack(df3_site,panel.fun=function(region,value,...){circos.genomicLines(region,value,type="h",col=color_list[i])},ylim=c(0,1),bg.border="#CCCCCC",track.height=0.06,track.margin=c(0.005,0.003))
#   circos.genomicDensity(df3_site, col = color_list[i], track.height = 0.05)
# }
# legend(x = 1, y = 0.5, lty = 1, col = color_list, legend = case_list, cex = 1, title = "Donor ID", box.lty = 0)
# circos.clear()
# dev.off()

################################################################################
##### 2. by condition
pdf(paste0(main_figure_dir, "/6-sSNV_circos_plot_by_condition.pdf"), width = 5, height = 3)
  circos.par(start.degree = 90)
  circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", 1:22), plotType = c("ideogram", "labels"), ideogram.height = 0.04)
  for(i in 1:length(condition_list)){
    df_site = genomic_context[genomic_context$Condition == condition_list[i], 1:2]
    df2_site = data.frame(Chr = sprintf("chr%s", df_site[, 1]), Start = df_site[, 2], End = df_site[, 2], Value = 1)
    df3_site = df2_site[!df2_site[, 1] %in% c("chrX", "chrY"), ]
    circos.genomicDensity(df3_site, col = rev(ctrl_dis_color)[i], track.height = 0.15, lwd = 0.1)
  }
  legend(x = 0.6, y = 1.1, lty = 1, col = rev(ctrl_dis_color), legend = condition_list, cex = 1, title = "Condition", box.lty = 0)
  circos.clear()
dev.off()

################################################################################
##### 3. compute summary dataframe to test if one chromosome is different from the rest
vcf_df_normal <- genomic_context[genomic_context$Case_ID %in% unique(metadata_df_AMG$Case_ID[metadata_df_AMG$Condition == ctrl_name]), ]
vcf_df_disease <- genomic_context[genomic_context$Case_ID %in% unique(metadata_df_AMG$Case_ID[metadata_df_AMG$Condition == dis_name]), ]

summary_normal_df <- vcf_df_normal %>% 
  group_by(Cell_ID, Chr) %>% 
  summarise(count = n()) %>% 
  mutate(percentage = count / sum(count)) %>% 
  ungroup() %>% 
  complete(Cell_ID, Chr, fill = list(count = 0)) %>% 
  arrange(match(Cell_ID, SCAN2_df_AMG$Cell_ID)) %>% 
  mutate(percentage = replace(percentage, is.na(percentage), 0)) %>% 
  mutate(est_snv_burden = rep(SCAN2_df_AMG$snv.rate.per.gb[SCAN2_df_AMG$Condition == ctrl_name], each = 22)) %>% 
  mutate(chr_seqlength = rep(seqlengths_list, length(unique(vcf_df_normal$Cell_ID))) / 1000000000) %>% 
  mutate(Chr_burden = percentage * est_snv_burden, condition = ctrl_name) %>% 
  mutate(est_chr_burden_density = Chr_burden / chr_seqlength)

summary_disease_df <- vcf_df_disease %>% 
  group_by(Cell_ID, Chr) %>% 
  summarise(count = n()) %>% 
  mutate(percentage = count / sum(count)) %>% 
  ungroup() %>% 
  complete(Cell_ID, Chr, fill = list(count = 0)) %>% 
  arrange(match(Cell_ID, SCAN2_df_AMG$Cell_ID)) %>% 
  mutate(percentage = replace(percentage, is.na(percentage), 0)) %>% 
  mutate(est_snv_burden = rep(SCAN2_df_AMG$snv.rate.per.gb[SCAN2_df_AMG$Condition == dis_name], each = 22)) %>% 
  mutate(chr_seqlength = rep(seqlengths_list, length(unique(vcf_df_disease$Cell_ID))) / 1000000000) %>% 
  mutate(Chr_burden = percentage * est_snv_burden, condition = dis_name) %>% 
  mutate(est_chr_burden_density = Chr_burden / chr_seqlength)

summary_all_chr_df <- rbind(summary_normal_df, summary_disease_df) |> base::`[`(c("Cell_ID", "Chr", "est_chr_burden_density", "condition", "percentage")) %>% 
  mutate(condition = factor(condition, level = c(ctrl_name, dis_name)))

percentage_p.value_list <- numeric()
for (chr_number in 1:22){
  chr_wilcox.test_result <- wilcox.test(summary_normal_df$percentage[summary_normal_df$Chr == chr_number], 
                                   summary_disease_df$percentage[summary_disease_df$Chr == chr_number], "two.sided", exact = FALSE)
  percentage_p.value_list <- c(percentage_p.value_list, round(chr_wilcox.test_result$p.value, digits = 3))
}

avg_df <- summary_all_chr_df %>% group_by(condition, Chr) %>% 
  summarize(avg_burden_density = mean(est_chr_burden_density), 
            avg_percentage = mean(percentage), 
            se = sd(percentage) / sqrt(n()), 
            CI_lower = t.test(percentage)$conf.int[1], 
            CI_upper = t.test(percentage)$conf.int[2])

avg_ctrl_df <- avg_df[avg_df$condition == ctrl_name, ]
avg_dis_df <- avg_df[avg_df$condition == dis_name, ]
merged_avg_df <- merge(avg_ctrl_df, avg_dis_df, by = "Chr") %>% 
  mutate(percentage.p.value = percentage_p.value_list, 
         fdr = p.adjust(percentage.p.value, method = "fdr"), 
         burden_dot_size = 10/log10(abs(avg_burden_density.x - avg_burden_density.y)), 
         percentage_dot_size = 2 - 10/log10(6*abs(avg_percentage.x - avg_percentage.y)), 
         significance = case_when(percentage.p.value < 0.001 ~ "***", percentage.p.value < 0.01 ~ "**", percentage.p.value < 0.05 ~ "*", TRUE ~ ""))

p_chr_percentage <- ggplot(merged_avg_df, aes(x = avg_percentage.x, y = avg_percentage.y, color = percentage_p.value_list)) +
  geom_segment(aes(x = 0, xend = 0.1, y = 0, yend = 0.1), colour = "grey20", linetype = "dashed", size = 1.5) + geom_point(size = 10) + 
  geom_errorbar(aes(ymin = avg_percentage.y - se.y, ymax = avg_percentage.y + se.y)) + geom_errorbarh(aes(xmin = avg_percentage.x - se.x, xmax = avg_percentage.x + se.x)) + 
  geom_text(aes(label = Chr), hjust = 0.5, vjust = 0.5, color = "white", size = 5) + 
  # geom_text(aes(label = significance), hjust = 0.5, vjust = 1.5, size = 15) + 
  annotate("text", size = 15, x = 0.0325, y = 0.0225, label = "*", hjust = 0) + annotate("text", size = 15, x = 0.026, y = 0.0175, label = "*", hjust = 0) +
  scale_x_continuous("avg. chr SNV percentage in Control") + scale_y_continuous("avg. chr SNV percentage in IHD") + 
  scale_color_gradient(low = "red", high = "blue", limits = c(0.01, 1), breaks = c(0.01, 0.05, 1), name = "P value", transform = "log10") + theme_linedraw() + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "grey80", linetype = "dashed", linewidth = 0.25), panel.grid.minor = element_blank(), 
        panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
ggsave(paste0(other_figure_dir, "/6-chr_percentage_by_condition.pdf"), plot = p_chr_percentage, width = 10, height = 8, dpi = 600)