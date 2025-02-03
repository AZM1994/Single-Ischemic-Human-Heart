##### circos plot
# vcf_file_list <- list.files(path = paste0(project_dir, "/all_vcfs_in_one"), pattern = ".vcf", full.names = TRUE)
# vcf_df <- rbindlist(sapply(vcf_file_list, fread, simplify = FALSE))
vcf_df <- read.csv(paste0(project_dir, "/all_vcfs_in_one/ssnv_all_vcf_age_matched.csv"))
vcf_df$Case_ID <- str_extract(vcf_df$Cell_ID, "[^_]+")
control_sample_list <- metadata_df$Case_ID[metadata_df$Condition == control_name]
disease_sample_list <- metadata_df$Case_ID[metadata_df$Condition == disease_name]
control_condition <- grepl(paste(control_sample_list, collapse = "|"), vcf_df$Case_ID)
disease_condition <- grepl(paste(disease_sample_list, collapse = "|"), vcf_df$Case_ID)
vcf_df$condition <- with(vcf_df, ifelse(control_condition, control_name,
                                            ifelse(disease_condition, disease_name, 'unknown')))

case_list <- unique(vcf_df$Case_ID)
condition_list <- unique(vcf_df$condition)
disease_color_palette<-colorRampPalette(c("pink","firebrick"))
control_color_palette<-colorRampPalette(c("skyblue","dodgerblue4"))
dis_ctrl_color <- c(control_color_palette(9)[7], disease_color_palette(4)[3])
color_list <- c(control_color_palette(3), disease_color_palette(5))
sSNV_plot_dir <- paste0(project_dir, "/sSNV_plots")

##### by donor all together
pdf(paste0(sSNV_plot_dir, "/6-ssnv_circos_plot_by_case.pdf"), width = 10, height = 6)
circos.par(start.degree=90)
circos.initializeWithIdeogram(species="hg19",chromosome.index=paste0("chr",1:22),plotType=c("ideogram","labels"),ideogram.height=0.04)
for(i in 1:length(case_list)){
  df_site=vcf_df[vcf_df$Case_ID==case_list[i],1:2]
  df2_site=data.frame(Chr=sprintf("chr%s",df_site[,1]),Start=df_site[,2],End=df_site[,2],Value=1)
  df3_site=df2_site[!df2_site[,1] %in% c("chrX","chrY"),]
  # circos.genomicTrack(df3_site,panel.fun=function(region,value,...){circos.genomicLines(region,value,type="h",col=color_list[i])},ylim=c(0,1),bg.border="#CCCCCC",track.height=0.06,track.margin=c(0.005,0.003))
  circos.genomicDensity(df3_site, col = color_list[i], track.height = 0.05)
}
legend(x='topright',lty=1,col=color_list,legend=case_list,cex=1,title="Case ID",box.lty=0)
circos.clear()
dev.off()

##### by condition
pdf(paste0(sSNV_plot_dir, "/6-ssnv_circos_plot_by_condition.pdf"), width = 10, height = 6)
circos.par(start.degree=90)
circos.initializeWithIdeogram(species="hg19",chromosome.index=paste0("chr",1:22),plotType=c("ideogram","labels"),ideogram.height=0.04)
for(i in 1:length(condition_list)){
  # pdf(paste0(sSNV_plot_dir, '/ssnv_', condition_list[i], "_circos_plot.pdf"), width = 10, height = 6)
  df_site=vcf_df[vcf_df$condition==condition_list[i],1:2]
  df2_site=data.frame(Chr=sprintf("chr%s",df_site[,1]),Start=df_site[,2],End=df_site[,2],Value=1)
  df3_site=df2_site[!df2_site[,1] %in% c("chrX","chrY"),]
  # circos.genomicTrack(df3_site,panel.fun=function(region,value,...){circos.genomicLines(region,value,type="h",col=dis_ctrl_color[i])},ylim=c(0,1),bg.border="#CCCCCC",track.height=0.06,track.margin=c(0.005,0.003))
  circos.genomicDensity(df3_site, col = dis_ctrl_color[i], track.height = 0.15)
  # circos.genomicDensity(df3_site, col = dis_ctrl_color[i], count_by = "number", track.height = 0.1)
  # dev.off()
}
legend(x='topright',lty=1,col=dis_ctrl_color,legend=condition_list,cex=1,title="Case ID",box.lty=0)
circos.clear()
dev.off()

##### compute summary dataframe to test if one chromsome is different from the rest
vcf_df_disease <- vcf_df[vcf_df$Case_ID %in% c("604", "1673", "1113", "1363", "1743"),]
vcf_df_normal <- vcf_df[vcf_df$Case_ID %in% c("1039", "5919", "5828"),]

post_hoc_multiple_comparison(vcf_df_disease)
input_vcf_df <- vcf_df_disease
input_vcf_df <- vcf_df_normal
input_vcf_df <- vcf_df

post_hoc_multiple_comparison <- function(input_vcf_df, rep_number) {
  summary_df <- input_vcf_df %>% group_by(Chr, Cell_ID) %>% summarise(count = n())
  summary_df$seqlength <- rep(seqlengths_list, each = rep_number) / (1024*1024*1024)
  summary_df$density <- summary_df$count / summary_df$seqlength
  
  remove_outliers_by_group <- function(df) {
    df %>%
      group_by(Chr) %>%
      filter(density >= quantile(density, 0.25) - 1.5 * IQR(density) &
               density <= quantile(density, 0.75) + 1.5 * IQR(density))}
  summary_df <- remove_outliers_by_group(summary_df)
  mean_density <- mean(summary_df$density)
    
  # Performing ANOVA
  # model <- aov(count ~ factor(Chr), data = summary_df)
  model <- aov(density ~ factor(Chr), data = summary_df)
  # tukey_result <- TukeyHSD(model)
  
  pvalue_list <- list()
  weight_ori <- c(1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21, 1/21)
  for(i in 1:length(unique(input_vcf_df$Chr))){
    weight_chr <- weight_ori
    weight_chr[i] <- -1
    com_index <- paste(i, "vs Rest")
    contrast_spec <- list(com_index = weight_chr)
    # Perform contrast analysis
    contrast_results <- contrast(emmeans(model, "Chr"), contrast_spec)
    # View the results
    # print(summary(contrast_results)$p.value)
    pvalue_list[[i]] <- summary(contrast_results)$p.value
  }
  pvalue_df <- data.frame(iteration = 1:22, result = unlist(pvalue_list))
  colnames(pvalue_df) <- c("chr", "p_value")
  pvalue_df$p_value <- round(pvalue_df$p_value, 3)
  a <- transform(pvalue_df, sig. = stars.pval(pvalue_df[,2]))
  t_a <- t(a)
  
  summary_df$label <- factor(summary_df$Chr)
  p_multiple_comparison <- ggplot(summary_df, aes(x = label, y = density, color = label)) +
    geom_boxplot(outlier.shape = NA) +
    # annotate("text", x = 1:22, y = rep(590,22), label = a$sig., size = 5, color = "red") +
    # scale_y_continuous(limits = c(0,600)) +
    annotate("text", x = 1:22, y = rep(490,22), label = a$sig., size = 5, color = "red") +
    scale_y_continuous(limits = c(0,500)) +
    geom_segment(x = 0, xend = 22.5, y = mean_density, yend = mean_density, colour = "grey50", linetype = "dashed") +
    labs(x = "chromosome", y = "sSNV/GB") + theme_classic() +
    theme(text = element_text(size=24), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position="none")
  print(p_multiple_comparison)
  # ggsave(paste0(figures_for_manuscript_dir, "/6-chromsome_density_disease.png"),
  #        plot = p_multiple_comparison, width = 12, height = 6, dpi = 300)
  # ggsave(paste0(figures_for_manuscript_dir, "/6-chromsome_density_normal.png"),
  #        plot = p_multiple_comparison, width = 12, height = 6, dpi = 300)
  ggsave(paste0(figures_for_manuscript_dir, "/6-chromsome_density_all.png"),
         plot = p_multiple_comparison, width = 12, height = 6, dpi = 300)
}

############ for indel
vcf_df <- read.csv(paste0(project_dir, "/all_vcfs_in_one/sindel_all_vcf_age_matched.csv"))
vcf_df$Case_ID <- str_extract(vcf_df$Cell_ID, "[^_]+")
control_sample_list <- metadata_df$Case_ID[metadata_df$Condition == control_name]
disease_sample_list <- metadata_df$Case_ID[metadata_df$Condition == disease_name]
control_condition <- grepl(paste(control_sample_list, collapse = "|"), vcf_df$Case_ID)
disease_condition <- grepl(paste(disease_sample_list, collapse = "|"), vcf_df$Case_ID)
vcf_df$condition <- with(vcf_df, ifelse(control_condition, control_name,
                                        ifelse(disease_condition, disease_name, 'unknown')))

case_list <- unique(vcf_df$Case_ID)
condition_list <- unique(vcf_df$condition)
disease_color_palette<-colorRampPalette(c("pink","firebrick"))
control_color_palette<-colorRampPalette(c("skyblue","dodgerblue4"))
dis_ctrl_color <- c(control_color_palette(9)[7], disease_color_palette(4)[3])
color_list <- c(control_color_palette(3), disease_color_palette(5))
sindel_plot_dir <- paste0(project_dir, "/sindel_plots")

##### by cell all together
pdf(paste0(sindel_plot_dir, "/6-sindel_circos_plot.pdf"), width = 10, height = 6)
circos.par(start.degree=90)
circos.initializeWithIdeogram(species="hg19",chromosome.index=paste0("chr",1:22),plotType=c("ideogram","labels"),ideogram.height=0.04)
for(i in 1:length(case_list)){
  df_site=vcf_df[vcf_df$Case_ID==case_list[i],1:2]
  df2_site=data.frame(Chr=sprintf("chr%s",df_site[,1]),Start=df_site[,2],End=df_site[,2],Value=1)
  df3_site=df2_site[!df2_site[,1] %in% c("chrX","chrY"),]
  # circos.genomicTrack(df3_site,panel.fun=function(region,value,...){circos.genomicLines(region,value,type="h",col=color_list[i])},ylim=c(0,1),bg.border="#CCCCCC",track.height=0.06,track.margin=c(0.005,0.003))
  circos.genomicDensity(df3_site, col = color_list[i], track.height = 0.05)
}
circos.clear()
legend(x='topright',lty=1,col=color_list,legend=case_list,cex=1,title="Case ID",box.lty=0)
# dev.off()

circos.clear()
# circos.par(start.degree=90)
# circos.initializeWithIdeogram(species="hg19",chromosome.index=paste0("chr",1:22),plotType=c("ideogram","labels"),ideogram.height=0.04)

##### by condition separate
circos.par(start.degree=90)
circos.initializeWithIdeogram(species="hg19",chromosome.index=paste0("chr",1:22),plotType=c("ideogram","labels"),ideogram.height=0.04)
for(i in 1:length(condition_list)){
  # pdf(paste0(sindel_plot_dir, '/sindel_', condition_list[i], "_circos_plot.pdf"), width = 10, height = 6)
  df_site=vcf_df[vcf_df$condition==condition_list[i],1:2]
  df2_site=data.frame(Chr=sprintf("chr%s",df_site[,1]),Start=df_site[,2],End=df_site[,2],Value=1)
  df3_site=df2_site[!df2_site[,1] %in% c("chrX","chrY"),]
  # circos.genomicTrack(df3_site,panel.fun=function(region,value,...){circos.genomicLines(region,value,type="h",col=dis_ctrl_color[i])},ylim=c(0,1),bg.border="#CCCCCC",track.height=0.06,track.margin=c(0.005,0.003))
  circos.genomicDensity(df3_site, col = dis_ctrl_color[i], track.height = 0.15)
  # circos.genomicDensity(df3_site, col = dis_ctrl_color[i], count_by = "number", track.height = 0.1)
  # dev.off()
}
circos.clear()
legend(x='topright',lty=1,col=dis_ctrl_color,legend=condition_list,cex=1,title="Case ID",box.lty=0)
dev.off()

unique(vcf_df$Case_ID)
vcf_df_disease <- vcf_df[vcf_df$Case_ID %in% c("604","1113","1363","1743"),]
disease_chromosome <- data.frame(table(vcf_df_disease$Chr))
disease_chromosome$condition <- "disease"
disease_chromosome$seqlength <- seqlengths_list/(1024*1024*1024)
colnames(disease_chromosome) <- c('chromosome', "sindel")

vcf_df_normal <- vcf_df[vcf_df$Case_ID %in% c("1673","1039","5919","5828"),]
normal_chromosome <- data.frame(table(vcf_df_normal$Chr))
normal_chromosome$condition <- "normal"
normal_chromosome$seqlength <- seqlengths_list/(1024*1024*1024)
colnames(normal_chromosome) <- c('chromosome', "sindel")

chromosome_summary <- rbind(normal_chromosome, disease_chromosome)
colnames(chromosome_summary)[3:4] <- c("condition", "seq_length")

chromosome_summary$ssnv_per_GB <- chromosome_summary$sindel / chromosome_summary$seq_length

dis_nor_ratio <- chromosome_summary$ssnv_per_GB[chromosome_summary$condition == "disease"] / chromosome_summary$ssnv_per_GB[chromosome_summary$condition == "normal"]
chromosome_summary$dis_nor_ratio <- rep(dis_nor_ratio,2)
# max_ssnv <- max(chromosome_summary$normal_ssnv_per_GB, chromosome_summary$disease_ssnv_per_GB)

pdf(paste0(sindel_plot_dir, "/sindel_chromosome_distribution.pdf"), width = 10, height = 6)
p1 <- ggplot(chromosome_summary, aes(x = chromosome, y = ssnv_per_GB, fill = condition)) +
  geom_point(pch = 21, size = 5, fill = c(rep(control_color_palette(9)[7], 22), rep(disease_color_palette(4)[3], 22))) +
  labs(x = "chromosome", y = "sindels/GB") + theme_classic() +
  theme(text = element_text(size=24), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right")
# print(p1)

p2 <- ggplot(chromosome_summary, aes(x = chromosome, y = dis_nor_ratio)) +
  geom_point(pch = 21, size = 5, fill = 'black') +
  labs(x = "chromosome", y = "disease / normal") + theme_classic() +
  theme(text = element_text(size=24), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right")
# print(p1,p2)
grid.arrange(p1, p2, ncol=1, nrow =2)
dev.off()