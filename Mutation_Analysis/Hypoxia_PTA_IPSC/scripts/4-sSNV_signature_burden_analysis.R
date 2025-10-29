################################################################################
######### 1. Get signature contribution from SigNet                   ##########
######### 2. Determine the top contribution signatures.               ##########
######### 3. Find signatures that are significantly different between ##########
#########    the control and disease conditions                       ##########
################################################################################
library(ComplexHeatmap)
library(circlize)
library(ggbeeswarm)
library(scales)
##### get COSMIC signature spectrum
COSMIC_v3.4_sigs <- read.table("main/COSMIC/COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE) %>% 
  {rownames(.) <- .$Type; .} %>% dplyr::select(-Type) %>% as.matrix()
SigNet_contri <- read.csv(paste0(project_dir, "/data/SigNet/PTA_IPSC_est/weight_guesses.csv"), row.names = 1) %>% as.data.frame()

func_Reshape <- function(a, n, m){
  if (missing(m)) m <- length(a)%/%n 
  if (length(a) != n * m) 
    stop("Matrix 'a' does not have n*m elements")
  dim(a) <- c(n, m)
return(a)
  }

##### determine top contribution signatures in age-matched samples
sSNV_sig_contri_rel <- 1 / rowSums(SigNet_contri) * SigNet_contri
sSNV_sig_contri_est <- sSNV_sig_contri_rel * sSNV_SCAN2_df$snv.rate.per.gb

sig_colors <- c(SBS1 = "#acdb88", SBS2 = "#91a64f", SBS3 = "#f6dc87", SBS4 = "#fca85c", SBS5 = "#f2b8c0", SBS7a = "#d98284", SBS8 = "#b94e4e", SBS10a = "#853232", 
                SBS11 = "#e3a6cd", SBS12 = "#cb77b2", SBS16 = "#a94d98", SBS18 = "#c0dcf0", SBS19 = "#cde6f6", SBS24 = "#a4e4e0", SBS29 = "#40bfc1", SBS30 = "#6fb0da", SBS31 = "#a4e4e5", 
                SBS32 = "#3f91c2", SBS36 = "#1f6aa5", SBS37 = "#e2dff0", SBS39 = "#cbcce1", SBS40a = "#b5b5d2", SBS44 = "#8882bd", SBS84 = "#d8bda3", SBS86 = "#c6a589", SBS87 = "#a97f63", 
                SBS89 = "#8a6a4f", SBS92 = "#6b4f39", others = "grey")

plot_type_list <- c("relative", "absolute")
# plot_type_list <- c("relative")
for (plot_type in plot_type_list) {
  if (plot_type == "relative") {
    sSNV_sig_contri <- sSNV_sig_contri_rel[Cell_ID_list, ]
  }
  else if (plot_type == "absolute") {
    sSNV_sig_contri <- sSNV_sig_contri_est[Cell_ID_list, ]
  }
  top_sSNV_sig_contri <- sort(colSums(sSNV_sig_contri), decreasing = TRUE)[1 : 20]
  top_sigs <- row.names(data.frame(top_sSNV_sig_contri))
  top_sigs_reordered <- top_sigs[order(as.numeric(gsub("[^0-9]", "", top_sigs)))]
  
  ## by condition
  top_sSNV_sig_contri_cond <- sSNV_sig_contri %>% 
    mutate(others = rowSums(.[, !(colnames(.) %in% top_sigs)])) %>% t() %>% as.data.frame() %>% 
    mutate(Normoxia = rowMeans(.[, ctrl_range]), Hypoxia = rowMeans(.[, dis_range])) %>% 
    dplyr::select(all_of(ctrl_name), all_of(dis_name)) %>% 
    mutate(Normoxia_Pct = Normoxia / sum(Normoxia) * 100, Hypoxia_Pct = Hypoxia / sum(Hypoxia) * 100, 
           Normoxia_label = paste0(rownames(.), " (", round(Normoxia_Pct, 1), "%)"), Hypoxia_label = paste0(rownames(.), " (", round(Hypoxia_Pct, 1), "%)"), 
           Normoxia_label = ifelse(Normoxia_Pct < 0.7, NA, Normoxia_label), Hypoxia_label = ifelse(Hypoxia_Pct < 0.7, NA, Hypoxia_label)) %>% 
    filter(rownames(.) %in% c(top_sigs, "others"))
  sorted_labels <- as.data.frame(row.names(top_sSNV_sig_contri_cond)) %>% mutate(A = .[[1]]) %>% t() %>% as.matrix() %>% func_Reshape(nrow(top_sSNV_sig_contri_cond) * 2, 1)
  dim(sorted_labels) <- c(nrow(top_sSNV_sig_contri_cond) * 2,1)
  label_size <- 4 * rbind(top_sSNV_sig_contri_cond[,1] / top_sSNV_sig_contri_cond[,1], top_sSNV_sig_contri_cond[,2] / top_sSNV_sig_contri_cond[,2])
  label_size <- func_Reshape(label_size,ncol(label_size)*2,1)
  label_size[is.na(label_size)] <- 0
  
  p_sSNV_top_sigs_contri <- plot_contribution(top_sSNV_sig_contri_cond[, 1:2], coord_flip = FALSE, mode = plot_type) + 
    scale_fill_manual(name = "Signature", values = sig_colors) + geom_bar(stat = "identity", colour = NA, linewidth = 0) + theme_linedraw() + 
    # aes(label = sorted_labels) + geom_text(size = label_size, position = position_stack(vjust = 0.5), col = "white", fontface = "bold") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12))
  
  if (plot_type == "relative") {
    write.csv(sSNV_sig_contri, file = paste0(table_dir, "/4-sSNV_sig_contri_rel_AMG_2n_PTA_all.csv"))
    top_sSNV_sig_contri_cond_save <- sSNV_sig_contri %>% 
      mutate(others = rowSums(.[, !(colnames(.) %in% colnames(sSNV_sig_contri))])) %>% t() %>% as.data.frame() %>% 
      mutate(Normoxia = rowMeans(.[, ctrl_range]), Hypoxia = rowMeans(.[, dis_range])) %>% 
      dplyr::select(all_of(ctrl_name), all_of(dis_name)) %>% 
      mutate(Normoxia_Pct = Normoxia / sum(Normoxia) * 100, Hypoxia_Pct = Hypoxia / sum(Hypoxia) * 100, 
             Normoxia_label = paste0(rownames(.), " (", round(Normoxia_Pct, 1), "%)"), Hypoxia_label = paste0(rownames(.), " (", round(Hypoxia_Pct, 1), "%)"), 
             Normoxia_label = ifelse(Normoxia_Pct < 0.7, NA, Normoxia_label), Hypoxia_label = ifelse(Hypoxia_Pct < 0.7, NA, Hypoxia_label)) %>% 
      filter(rownames(.) %in% c(colnames(sSNV_sig_contri), "others"))
    write.csv(top_sSNV_sig_contri_cond_save, file = paste0(table_dir, "/4-top_sSNV_sig_contri_cond_2n_PTA_all.csv"))
    ggsave(paste0(suppl_figure_dir, "/4-sSNV_top_sigs_contri_", plot_type, ".pdf"), plot = p_sSNV_top_sigs_contri, width = 3.5, height = 7, dpi = 300)
    ## by cell
    p_sSNV_top_sigs_contri_bycell_heatmap <- plot_contribution_heatmap(t(sSNV_sig_contri[, top_sigs_reordered]), cluster_sigs = FALSE, cluster_samples = FALSE) + 
      scale_y_discrete(limits = rev(Cell_ID_list)) + scale_fill_gradient(low = "white", high = "dodgerblue2") + geom_tile(colour = "black") + theme_linedraw() + 
      theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12))
    ggsave(paste0(other_figure_dir, "/4-sSNV_top_sigs_contri_bycell_", plot_type, "_heatmap.pdf"), plot = p_sSNV_top_sigs_contri_bycell_heatmap, width = 12, height = 10, dpi = 300)
    
    p_sSNV_top_sigs_contri_bycell_barplot <- plot_contribution(t(sSNV_sig_contri[, top_sigs_reordered]), coord_flip = FALSE, mode = plot_type) + 
      scale_x_discrete(limits = Cell_ID_list) + guides(fill = guide_legend(ncol = 2)) + scale_fill_manual(name = "Signature", values = sig_colors) + theme_linedraw() + 
      theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.0))
    ggsave(paste0(other_figure_dir, "/4-sSNV_top_sigs_contri_bycell_", plot_type, "_barplot.pdf"), plot = p_sSNV_top_sigs_contri_bycell_barplot, width = 12, height = 8, dpi = 300)
  }
  else if (plot_type == "absolute") {
    ggsave(paste0(main_figure_dir, "/4-sSNV_top_sigs_contri_", plot_type, ".pdf"), plot = p_sSNV_top_sigs_contri, width = 5, height = 7.5, dpi = 300)
  }
}

##############################################################################
### add metadata and remove the signatures with all zeros in one Condition
sSNV_sig_contri_est <- sSNV_sig_contri_est %>% 
  mutate(Age = metadata_df$Age, Condition = relevel(factor(metadata_df$Condition), ref = ctrl_name), 
         Gender = metadata_df$Gender, Color = sSNV_SCAN2_df$Color, Outline_Color = sSNV_SCAN2_df$Outline_Color) %>% 
  {
    numeric_mat <- .[, 1:(ncol(.) - 5)]
    zero_col <- which(colSums(numeric_mat) == 0)
    if (length(zero_col) > 0) {
      .[, -zero_col]
    } else {
      .
    }
    }

refined_selected_sigs_list <- colnames(sSNV_sig_contri_est)[1 : (ncol(sSNV_sig_contri_est) - 5)]
sig_digits <- 2

##### Mix effects model & Age matched mix effects model
pdf(paste0(other_figure_dir, "/4-lmer_top_sigs_all.pdf"), width = 8, height = 6)
  for (i in c(1 : (ncol(sSNV_sig_contri_est) - 5))){
    cat("index", i, ":", refined_selected_sigs_list[i], "\n")
    ## burden_df: the input matrix with columns: burden, Condition, Case_ID, Age, Color
    burden_df <- sSNV_sig_contri_est[c(refined_selected_sigs_list[i], "Age", "Condition", "Color", "Outline_Color")] %>% 
      setNames(c("snv.rate.per.gb", "Age", "Condition", "Color", "Outline_Color")) %>% 
      mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% 
      mutate(Condition = relevel(Condition, ref = ctrl_name)) %>% 
      mutate(label = c(rep("Normoxia (n = 9)", 9), rep("Hypoxia (n = 9)", 9))) %>% 
      group_by(Condition) %>% arrange(Age, .by_group = TRUE)
    burden_df_ctrl <- burden_df %>% filter(Condition == ctrl_name)
    burden_df_dis <- burden_df %>% filter(Condition == dis_name)
    
    wilcox.test_result_sig <- wilcox.test(burden_df_ctrl$snv.rate.per.gb, burden_df_dis$snv.rate.per.gb, alternative = "two.sided")$p.value
    if (wilcox.test_result_sig < 0.0099){
      wilcox.test_result_sig_print <- format(wilcox.test_result_sig, scientific = TRUE, digits = 2)
    } else{
      wilcox.test_result_sig_print <- formatC(signif(wilcox.test_result_sig, digits = 2), digits = 2, format="fg", flag="#")
    }
    
    boxplot_sSNV_burden_sig <- ggplot(burden_df, aes(x = factor(label, level = c("Normoxia (n = 9)", "Hypoxia (n = 9)")), y = snv.rate.per.gb, fill = label)) + 
      geom_boxplot(color = c("dodgerblue3", "firebrick3"), fill = c("#C7DFF0", "#FBDFE2")) + 
      geom_quasirandom(pch = 21, fill = burden_df$Color, color = burden_df$Outline_Color, size = 5) + 
      annotate("text", size = 7, x = 1.32, y = max(burden_df$snv.rate.per.gb), label = paste("P = ", wilcox.test_result_sig_print), hjust = 0) + 
      labs(x = "", y = paste0(refined_selected_sigs_list[i], " Contribution \n (sSNV rate per GB)")) + 
      stat_summary(fun = mean, colour = "black", geom = "point", shape = 18, size = 3, show.legend = FALSE) + 
      stat_compare_means(comparisons = list(c("Normoxia (n = 9)", "Hypoxia (n = 9)")), label.x = 1.6, label.y = max(burden_df$snv.rate.per.gb), 
                         label = "p.signif", bracket.size = 1, size = 7, method = "wilcox.test") + 
      scale_y_continuous(labels = number_format(accuracy = 0.1)) + theme_linedraw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), 
            text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
    print(boxplot_sSNV_burden_sig)
    if (wilcox.test_result_sig <= 0.05){
      cat("Significant \n")
      ggsave(paste0(main_figure_dir, "/2-boxplot_sSNV_burden_", refined_selected_sigs_list[i], ".pdf"), plot = boxplot_sSNV_burden_sig, width = 9, height = 8)
    }
  }
dev.off()
