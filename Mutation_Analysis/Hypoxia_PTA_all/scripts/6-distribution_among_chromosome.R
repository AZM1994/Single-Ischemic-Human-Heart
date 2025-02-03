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
  # genomic_SCAN2_df_temp <- merge(genomic_context_temp, metadata_df[c("Cell_ID", "Case_ID", "Age", "Gender", "Condition")]) %>% 
  #   rename_with(~ c("Cell_ID", "Case_ID", "Condition", "chr", "start", "end", "ref", "alt", "region", "gene"), .cols = 1:10) %>% 
  #   # filter(age >= 40 & age < 80) %>% 
  #   mutate(Condition = as.factor(Condition))
}
# vcf_df <- genomic_context %>% 
#   mutate(Condition = ifelse(grepl(paste(metadata_df$Case_ID[metadata_df$Condition == control_name], collapse = "|"), Case_ID), control_name,
#                             ifelse(metadata_df$Case_ID[metadata_df$Condition == disease_name], disease_name, 'unknown')))
vcf_df <- genomic_context

vcf_df2 <- read.csv(paste0(project_dir, "/data/vcfs/ssnv_all_vcf_age_matched.csv")) %>%
  mutate(Case_ID = str_extract(Cell_ID, "[^_]+")) %>%
  mutate(condition = ifelse(grepl(paste(metadata_df$Case_ID[metadata_df$Condition == control_name], collapse = "|"), Case_ID), control_name,
                                         ifelse(metadata_df$Case_ID[metadata_df$Condition == disease_name], disease_name, 'unknown')))

case_list <- unique(vcf_df$Case_ID)
condition_list <- rev(unique(vcf_df$Condition))
control_color_palette <- colorRampPalette(c("skyblue","dodgerblue4"))
disease_color_palette <- colorRampPalette(c("pink","firebrick"))
dis_ctrl_color <- c(control_color_palette(9)[7], disease_color_palette(4)[3])
color_list <- c(control_color_palette(3), disease_color_palette(5))

################################################################################
##### 1. by donor all together
# pdf(paste0(sSNV_figure_dir, "/6-ssnv_circos_plot_by_case.pdf"), width = 10, height = 6)
# circos.par(start.degree=90)
# circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr",1:22), 
#                               plotType = c("ideogram", "labels"), ideogram.height = 0.04)
# for(i in 1:length(case_list)){
#   df_site = vcf_df[vcf_df$Case_ID == case_list[i], 1:2]
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
pdf(paste0(main_figure_dir, "/6-ssnv_circos_plot_by_condition.pdf"), width = 5, height = 3)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", 1:22), plotType = c("ideogram", "labels"), ideogram.height = 0.04)
for(i in 1:length(condition_list)){
  df_site = vcf_df[vcf_df$Condition == condition_list[i], 1:2]
  df2_site = data.frame(Chr = sprintf("chr%s", df_site[, 1]), Start = df_site[, 2], End = df_site[, 2], Value = 1)
  df3_site = df2_site[!df2_site[, 1] %in% c("chrX", "chrY"), ]
  circos.genomicDensity(df3_site, col = rev(dis_ctrl_color)[i], track.height = 0.15, lwd = 0.1)
}
legend(x = 0.6, y = 1.1, lty = 1, col = rev(dis_ctrl_color), legend = condition_list, cex = 1, title = "Condition", box.lty = 0)
circos.clear()
dev.off()

################################################################################
##### 3. compute summary dataframe to test if one chromsome is different from the rest
vcf_df_normal <- vcf_df[vcf_df$Case_ID %in% c("1039", "5919", "5828"),]
vcf_df_disease <- vcf_df[vcf_df$Case_ID %in% c("604", "1673", "1113", "1363", "1743"),]
# vcf_df_normal <- vcf_df[vcf_df$condition == "Control", ]
# vcf_df_disease <- vcf_df[vcf_df$condition == "IHD", ]

# normal_reorder_cellID <- c("1039", "5919_1_C4-2n", "5919_1_D2-2n", "5919_1_E3-2n", "5919_1_F6-2n", "5828_CM_C2_2n", "5828_CM_G2_2n")
# disease_reorder_cellID <- c("1363_A4", "1363_D4", "1363_H4", "1673_CM_A2_2n", "1673_CM_A3_2n", "1673_CM_D2_2n", 
#                             "604_B2", "604_B3", "604_B6", "1743_A3", "1743_C3", "1743_F3", "1113_D1", "1113_E1", "1113_F1")

summary_normal_df <- vcf_df_normal %>% 
  group_by(Cell_ID, Chr) %>% 
  summarise(count = n()) %>% 
  mutate(percentage = count / sum(count)) %>% 
  ungroup() %>%
  complete(Cell_ID, Chr, fill = list(count = 0)) %>% 
  arrange(match(Cell_ID, SCAN2_df$cell_ID)) %>% 
  mutate(percentage = replace(percentage, is.na(percentage), 0)) %>% 
  # mutate(estimated_snv_burden = rep(SCAN2_df$snv.burden[SCAN2_df$condition == control_name], each = 22)) %>%
  mutate(estimated_snv_burden = rep(SCAN2_df_age_match$snv.burden[SCAN2_df_age_match$condition == control_name], each = 22)) %>%
  mutate(chr_seqlength = rep(seqlengths_list, length(unique(vcf_df_normal$Cell_ID))) / (1024*1024*1024)) %>%
  mutate(Chr_burden = percentage * estimated_snv_burden, condition = control_name) %>% 
  mutate(estimated_chr_burden_density = Chr_burden / chr_seqlength)

summary_disease_df <- vcf_df_disease %>% 
  group_by(Cell_ID, Chr) %>%
  summarise(count = n()) %>% 
  mutate(percentage = count / sum(count)) %>% 
  ungroup() %>% 
  complete(Cell_ID, Chr, fill = list(count = 0)) %>% 
  arrange(match(Cell_ID, SCAN2_df$cell_ID)) %>% 
  mutate(percentage = replace(percentage, is.na(percentage), 0)) %>% 
  # mutate(estimated_snv_burden = rep(SCAN2_df$snv.burden[SCAN2_df$condition == disease_name], each = 22)) %>%
  mutate(estimated_snv_burden = rep(SCAN2_df_age_match$snv.burden[SCAN2_df_age_match$condition == disease_name], each = 22)) %>%
  mutate(chr_seqlength = rep(seqlengths_list, length(unique(vcf_df_disease$Cell_ID))) / (1024*1024*1024)) %>%
  mutate(Chr_burden = percentage * estimated_snv_burden, condition = disease_name) %>% 
  mutate(estimated_chr_burden_density = Chr_burden / chr_seqlength)

summary_all_chr_df <- rbind(summary_normal_df, summary_disease_df) |> base::`[`(c("Cell_ID", "Chr", "estimated_chr_burden_density", "condition", "percentage")) %>% 
  mutate(condition = factor(condition, level = c(control_name, disease_name)))

estimated_chr_burden_density_p.value_list <- numeric()
for (chr_number in 1:22){
  chr_wilcox.test_result <- t.test(summary_normal_df$estimated_chr_burden_density[summary_normal_df$Chr == chr_number], 
                              summary_disease_df$estimated_chr_burden_density[summary_disease_df$Chr == chr_number], "two.sided")
  estimated_chr_burden_density_p.value_list <- c(estimated_chr_burden_density_p.value_list, round(chr_wilcox.test_result$p.value, digits = 3))
}

percentage_p.value_list <- numeric()
for (chr_number in 1:22){
  chr_wilcox.test_result <- wilcox.test(summary_normal_df$percentage[summary_normal_df$Chr == chr_number], 
                                   summary_disease_df$percentage[summary_disease_df$Chr == chr_number], "two.sided")
  percentage_p.value_list <- c(percentage_p.value_list, round(chr_wilcox.test_result$p.value, digits = 3))
}

averages_df <- summary_all_chr_df %>%
  group_by(condition, Chr) %>%
  summarize(average_burden_density = mean(estimated_chr_burden_density), 
            average_percentage = mean(percentage),
            CI_lower = t.test(percentage)$conf.int[1],
            CI_upper = t.test(percentage)$conf.int[2])

averages_normal_df <- averages_df[averages_df$condition == control_name, ]
averages_disease_df <- averages_df[averages_df$condition == disease_name, ]
merged_average_df <- merge(averages_normal_df, averages_disease_df, by = "Chr") %>% 
  mutate(burden.p.value = estimated_chr_burden_density_p.value_list) %>% 
  mutate(percentage.p.value = percentage_p.value_list) %>% 
  mutate(burden_dot_size = 10/log10(abs(average_burden_density.x - average_burden_density.y))) %>% 
  # mutate(percentage_dot_size = (abs(average_percentage.x - average_percentage.y)))
  mutate(percentage_dot_size = 2 - 10/log10(6*abs(average_percentage.x - average_percentage.y)))

pdf(paste0(sSNV_figure_dir, "/6-chr_percentage_by_condition.pdf"), width = 9, height = 8)
scatter_plot_chr_percentage <- ggplot(merged_average_df, aes(x = average_percentage.x, y = average_percentage.y, color = percentage_p.value_list)) +
  geom_segment(aes(x = 0, xend = 0.1, y = 0, yend = 0.1), colour = "grey20", linetype = "dashed", size = 1.5) +
  geom_point(size = 10) + 
  geom_errorbar(aes(ymin = CI_lower.y, ymax = CI_upper.y)) + 
  geom_errorbarh(aes(xmin = CI_lower.x, xmax = CI_upper.x)) + 
  geom_text(aes(label = Chr), vjust = 0.5, color = "white", size = 5) + 
  # annotate("text", size = 15, x = 0.030, y = 0.043, label = "*", hjust = 0) + 
  annotate("text", size = 15, x = 0.037, y = 0.050, label = "*", hjust = 0) +
  scale_x_continuous("avg. chr SNV percentage in Control") + 
  scale_y_continuous("avg. chr SNV percentage in IHD") + 
  theme_classic() + 
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), 
        panel.background = element_rect(fill = "white"), legend.position = c(0.9,0.3)) + 
  scale_color_gradient(low = "red", high = "blue", limits = c(0.01, 1), breaks = c(0.01, 0.05, 1), name = "pvalue", transform = "log10")
print(scatter_plot_chr_percentage)
ggsave(paste0(figures_for_manuscript_dir, "/6-chr_percentage_by_condition.png"),
       plot = scatter_plot_chr_percentage, width = 9, height = 8, dpi = 600)
dev.off()