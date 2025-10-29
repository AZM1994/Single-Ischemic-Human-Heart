################################################################################
############################### circos plot ####################################
################################################################################
vcf_df <- read.csv(paste0(project_dir, "/data/vcfs/vcfs_by_condition/ssnv_all_vcf.csv")) %>% 
  mutate(Cell_ID = replace(Cell_ID, Cell_ID == "Normoxia2n_IPSC", "A2n")) %>% 
  mutate(Cell_ID = replace(Cell_ID, Cell_ID == "Hypoxia2n_IPSC", "B2n")) %>% 
  mutate(Case_ID = gsub(control_name, "", Cell_ID)) %>% 
  mutate(Case_ID = gsub(disease_name, "", Case_ID)) %>% 
  mutate(Case_ID = gsub("iPSC", "", Case_ID)) %>% 
  mutate(condition = ifelse(Case_ID %in% Cell_ID_list[control_range], control_name,
                            ifelse(Case_ID %in% Cell_ID_list[disease_range], disease_name, 'unknown')))

case_list <- unique(vcf_df$Case_ID)
condition_list <- unique(vcf_df$condition)
disease_color_palette<-colorRampPalette(c("pink","firebrick"))
control_color_palette<-colorRampPalette(c("skyblue","dodgerblue4"))
dis_ctrl_color <- c(control_color_palette(9)[7], disease_color_palette(4)[3])
color_list <- c(control_color_palette(9), disease_color_palette(9))

################################################################################
##### 1. by cell by condition
pdf(paste0(sSNV_figure_dir, "/6-ssnv_circos_plot_by_cell_normoxia.pdf"), width = 10, height = 6)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", 1:22), 
                              plotType = c("ideogram", "labels"), ideogram.height = 0.04)

for(i in control_range){
  df_site = vcf_df[vcf_df$Case_ID == case_list[i], 1:2]
  df2_site = data.frame(Chr = sprintf("chr%s", df_site[, 1]), Start = df_site[, 2], End = df_site[, 2], Value = 1)
  df3_site = df2_site[!df2_site[, 1] %in% c("chrX", "chrY"), ]
  # circos.genomicTrack(df3_site, panel.fun = function(region, value, ...){
  #   circos.genomicLines(region, value, type = "h", col = color_list[i])
  #   }, 
  #   ylim = c(0, 1), bg.border = "#CCCCCC", track.height = 0.06, track.margin = c(0.005, 0.003))
  circos.genomicDensity(df3_site, col = color_list[i], track.height = 0.05)
}
legend(x = 1, y = 0.5, lty = 1, col = color_list[control_range], legend = case_list[control_range], cex = 1, title = "Cell ID", box.lty = 0)
circos.clear()
dev.off()

pdf(paste0(sSNV_figure_dir, "/6-ssnv_circos_plot_by_cell_hypoxia.pdf"), width = 10, height = 6)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", 1:22), 
                              plotType = c("ideogram", "labels"), ideogram.height = 0.04)

for(i in disease_range){
  df_site = vcf_df[vcf_df$Case_ID == case_list[i], 1:2]
  df2_site = data.frame(Chr = sprintf("chr%s", df_site[, 1]), Start = df_site[, 2], End = df_site[, 2], Value = 1)
  df3_site = df2_site[!df2_site[, 1] %in% c("chrX", "chrY"), ]
  # circos.genomicTrack(df3_site, panel.fun = function(region, value, ...){
  #   circos.genomicLines(region, value, type = "h", col = color_list[i])
  # }, 
  # ylim = c(0, 1), bg.border = "#CCCCCC", track.height = 0.06, track.margin = c(0.005, 0.003))
  circos.genomicDensity(df3_site, col = color_list[i], track.height = 0.05)
}
legend(x = 1, y = 0.5, lty = 1, col = color_list[disease_range], legend = case_list[disease_range], cex = 1, title = "Cell ID", box.lty = 0)
circos.clear()
dev.off()

################################################################################
##### 2. by condition
pdf(paste0(main_figure_dir, "/6-ssnv_circos_plot_by_condition.pdf"), width = 5, height = 3)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", 1:22), plotType = c("ideogram", "labels"), ideogram.height = 0.04)
for(i in 1:length(condition_list)){
  df_site=vcf_df[vcf_df$condition==condition_list[i],1:2]
  df2_site=data.frame(Chr=sprintf("chr%s",df_site[,1]),Start=df_site[,2],End=df_site[,2],Value=1)
  df3_site=df2_site[!df2_site[,1] %in% c("chrX","chrY"),]
  # circos.genomicTrack(df3_site,panel.fun=function(region,value,...){circos.genomicLines(region,value,type="h",col=dis_ctrl_color[i])},ylim=c(0,1),bg.border="#CCCCCC",track.height=0.06,track.margin=c(0.005,0.003))
  circos.genomicDensity(df3_site, col = dis_ctrl_color[i], track.height = 0.15)
  # circos.genomicDensity(df3_site, col = dis_ctrl_color[i], count_by = "number", track.height = 0.1)
}
legend(x = 0.6, y = 1.1, lty = 1, col = dis_ctrl_color, legend = condition_list, cex = 1, title = "condition", box.lty = 0)
circos.clear()
dev.off()

################################################################################
##### 3. compute summary dataframe to test if one chromsome is different from the rest
vcf_df_normoxia <- vcf_df[vcf_df$condition == "Normoxia", ]
vcf_df_hypoxia <- vcf_df[vcf_df$condition == "Hypoxia", ]

normoxia_reorder_CaseID <- c("A1", "A2", "A3", "A6", "A7", "A8", "A9", "A11", "A2n")
hypoxia_reorder_CaseID <- c("B3", "B4", "B7", "B8", "B9", "B10", "B11", "B12", "B2n")

summary_normoxia_df <- vcf_df_normoxia %>% 
  mutate(Case_ID = factor(Case_ID, levels = normoxia_reorder_CaseID)) %>%
  group_by(Case_ID, Chr) %>% 
  summarise(count = n()) %>% 
  mutate(percentage = count / sum(count)) %>% 
  ungroup() %>% 
  complete(Case_ID, Chr, fill = list(count = 0)) %>% 
  mutate(percentage = replace(percentage, is.na(percentage), 0)) %>%
  mutate(estimated_snv_burden = rep(SCAN2_df$snv.burden[SCAN2_df$condition == control_name], each = 22)) %>%
  mutate(chr_seqlength = rep(seqlengths_list, length(normoxia_reorder_CaseID)) / (1024*1024*1024)) %>%
  mutate(Chr_burden = percentage * estimated_snv_burden, condition = control_name) %>% 
  mutate(estimated_chr_burden_density = Chr_burden / chr_seqlength)

summary_hypoxia_df <- vcf_df_hypoxia %>% 
  mutate(Case_ID = factor(Case_ID, levels = hypoxia_reorder_CaseID)) %>%
  group_by(Case_ID, Chr) %>%
  summarise(count = n()) %>% 
  mutate(percentage = count / sum(count)) %>% 
  ungroup() %>% 
  complete(Case_ID, Chr, fill = list(count = 0)) %>% 
  mutate(percentage = replace(percentage, is.na(percentage), 0)) %>%
  mutate(estimated_snv_burden = rep(SCAN2_df$snv.burden[SCAN2_df$condition == disease_name], each = 22)) %>%
  mutate(chr_seqlength = rep(seqlengths_list, length(hypoxia_reorder_CaseID)) / (1024*1024*1024)) %>%
  mutate(Chr_burden = percentage * estimated_snv_burden, condition = disease_name) %>% 
  mutate(estimated_chr_burden_density = Chr_burden / chr_seqlength)

summary_all_chr_df <- rbind(summary_normoxia_df, summary_hypoxia_df) |> base::`[`(c("Case_ID", "Chr", "estimated_chr_burden_density", "condition", "percentage")) %>% 
  mutate(condition = factor(condition, level = c(control_name, disease_name)))

estimated_chr_burden_density_p.value_list <- numeric()
for (chr_number in 1:22){
  chr_wilcox.test_result <- wilcox.test(summary_normoxia_df$estimated_chr_burden_density[summary_normoxia_df$Chr == chr_number], 
                              summary_hypoxia_df$estimated_chr_burden_density[summary_hypoxia_df$Chr == chr_number], "two.sided")
  estimated_chr_burden_density_p.value_list <- c(estimated_chr_burden_density_p.value_list, round(chr_wilcox.test_result$p.value, digits = 3))
}

percentage_p.value_list <- numeric()
for (chr_number in 1:22){
  chr_wilcox.test_result <- wilcox.test(summary_normoxia_df$percentage[summary_normoxia_df$Chr == chr_number], 
                              summary_hypoxia_df$percentage[summary_hypoxia_df$Chr == chr_number], "two.sided")
  percentage_p.value_list <- c(percentage_p.value_list, round(chr_wilcox.test_result$p.value, digits = 3))
}

averages_df <- summary_all_chr_df %>%
  group_by(condition, Chr) %>%
  summarize(average_burden_density = mean(estimated_chr_burden_density), 
            average_percentage = mean(percentage),
            CI_lower = t.test(percentage)$conf.int[1],
            CI_upper = t.test(percentage)$conf.int[2])

averages_normoxia_df <- averages_df[averages_df$condition == control_name, ]
averages_hypoxia_df <- averages_df[averages_df$condition == disease_name, ]
merged_average_df <- merge(averages_normoxia_df, averages_hypoxia_df, by = "Chr") %>% 
  mutate(burden.p.value = estimated_chr_burden_density_p.value_list) %>% 
  mutate(percentage.p.value = percentage_p.value_list) %>% 
  mutate(burden_dot_size = 10/log10(abs(average_burden_density.x - average_burden_density.y))) %>% 
  # mutate(percentage_dot_size = (abs(average_percentage.x - average_percentage.y)))
  mutate(percentage_dot_size = 2 - 10/log10(6*abs(average_percentage.x - average_percentage.y)))

pdf(paste0(sSNV_figure_dir, "/6-chr_percentage_by_condition.pdf"), width = 9, height = 8)
scatter_plot_chr_percentage <- ggplot(merged_average_df, aes(x = average_percentage.x, y = average_percentage.y, color = percentage.p.value)) +
  geom_segment(aes(x = 0, xend = 0.1, y = 0, yend = 0.1), colour = "grey20", linetype = "dashed", size = 1.5) +
  geom_point(size = 10) + 
  geom_errorbar(aes(ymin = CI_lower.y, ymax = CI_upper.y)) +
  geom_errorbarh(aes(xmin = CI_lower.x, xmax = CI_upper.x)) +
  geom_text(aes(label = Chr), vjust = 0.5, color = "white", size = 5) + 
  annotate("text", size = 12, x = 0.025, y = 0.041, label = "*", hjust = 0) + 
  annotate("text", size = 12, x = 0.0225, y = 0.0125, label = "*", hjust = 0) + 
  scale_x_continuous("avg. chr sSNV percentage in normoxia", limits = c(0, 0.1)) + 
  scale_y_continuous("avg. chr sSNV percentage in hypoxia", limits = c(0, 0.1)) + 
  theme_classic() + 
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5), 
        panel.background = element_rect(fill = "white"), legend.position = c(0.9,0.3)) + 
  scale_color_gradient(low = "red", high = "blue", limits = c(0.01, 1), breaks = c(0.01, 0.05, 1), name = "pvalue", transform = "log10")
print(scatter_plot_chr_percentage)
ggsave(paste0(figures_for_manuscript_dir, "/6-chr_percentage_by_condition.png"),
       plot = scatter_plot_chr_percentage, width = 9, height = 8, dpi = 600)
dev.off()
