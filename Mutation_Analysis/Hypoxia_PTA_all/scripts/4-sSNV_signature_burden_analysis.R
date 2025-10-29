################################################################################
######### 1. Get signature contribution from SigNet                   ##########
######### 2. Determine the top contribution signatures.               ##########
######### 3. Find signatures that are significantly different between ##########
#########    the control and disease conditions                       ##########
################################################################################

sig_colors <- c(SBS1 = "#acdb88", SBS2 = "#91a64f", SBS3 = "#f6dc87", SBS4 = "#fca85c", SBS5 = "#f2b8c0", SBS8 = "#b94e4e", 
                SBS11 = "#e3a6cd", SBS12 = "#cb77b2", SBS16 = "#a94d98", SBS19 = "#cde6f6", SBS24 = "#a4e4e0", SBS29 = "#40bfc1", SBS30 = "#6fb0da", 
                SBS32 = "#3f91c2", SBS36 = "#1f6aa5", SBS37 = "#e2dff0", SBS40a = "#b5b5d2", SBS44 = "#8882bd", SBS84 = "#d8bda3", SBS86 = "#c6a589", SBS87 = "#a97f63", 
                SBS89 = "#8a6a4f", SBS92 = "#6b4f39", others = "grey")

################################################################################
################ calculate all types of contribution dataframe #################
################################################################################
## get COSMIC signature spectrum
COSMIC_v3.4_sigs <- read.table("main/COSMIC/COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE) %>% 
  {rownames(.) <- .$Type; .} %>% dplyr::select(-Type) %>% as.matrix()

## get SigNet signature decomposition results
SigNet_contri <- read.csv(paste0(project_dir, "/data/SigNet/PTA_all_est/weight_guesses.csv"), row.names = 1) %>% as.data.frame()

## calculate relative and absolute contribution by cell
sSNV_sig_contri_rel <- (1 / rowSums(SigNet_contri) * SigNet_contri) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Cell_ID") %>% 
  inner_join(sSNV_SCAN2_df %>% select(Cell_ID, Case_ID, Condition, Age, snv.rate.per.gb), by = "Cell_ID") %>% 
  relocate(Cell_ID, Case_ID, Condition, Age, snv.rate.per.gb, .before = 1) %>% 
  column_to_rownames(var = "Cell_ID")
sSNV_sig_contri_est <- sSNV_sig_contri_rel %>% mutate(across(starts_with("SBS"), ~ .x * snv.rate.per.gb))
write.csv(sSNV_sig_contri_rel, file = paste0(table_dir, "/4-sSNV_sig_contri_rel_by_cell_PTA.csv"))
write.csv(sSNV_sig_contri_est, file = paste0(table_dir, "/4-sSNV_sig_contri_est_by_cell_PTA.csv"))

## calculate relative and absolute contribution by donor
sSNV_sig_contri_rel_bydonor <- sSNV_sig_contri_rel %>% 
  mutate(Case_ID = factor(Case_ID, levels = unique(Case_ID))) %>% 
  group_by(Case_ID) %>% 
  summarise(Condition = first(Condition), Age = first(Age), snv.rate.per.gb = mean(snv.rate.per.gb), 
            across(starts_with("SBS"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>% 
  arrange(Case_ID) %>% column_to_rownames(var = "Case_ID")
sSNV_sig_contri_est_bydonor <- sSNV_sig_contri_rel_bydonor %>% mutate(across(starts_with("SBS"), ~ .x * snv.rate.per.gb))
write.csv(sSNV_sig_contri_rel_bydonor, file = paste0(table_dir, "/4-sSNV_sig_contri_rel_by_donor_PTA.csv"))
write.csv(sSNV_sig_contri_est_bydonor, file = paste0(table_dir, "/4-sSNV_sig_contri_est_by_donor_PTA.csv"))

##### determine top contribution signatures in age-matched samples
##### plot contribution by cell, by donor, by condition (age-matched)
plot_type_list <- c("absolute")
# plot_type_list <- c("relative")
for (plot_type in plot_type_list) {
  if (plot_type == "relative") {
    sSNV_sig_contri_AMG <- sSNV_sig_contri_rel[Cell_ID_list_AMG, ]
    sSNV_sig_contri_AMG_bydonor <- sSNV_sig_contri_rel_bydonor[Case_ID_list_AMG, ]
  }
  else if (plot_type == "absolute") {
    sSNV_sig_contri_AMG <- sSNV_sig_contri_est[Cell_ID_list_AMG, ]
    sSNV_sig_contri_AMG_bydonor <- sSNV_sig_contri_est_bydonor[Case_ID_list_AMG, ]
  }
  
  ## find top signatures from age-matched samples
  top_sSNV_sig_contri <- sort(colSums(sSNV_sig_contri_AMG %>% dplyr::select(starts_with("SBS")), na.rm = TRUE), decreasing = TRUE)[1:16]
  top_sigs <- row.names(data.frame(top_sSNV_sig_contri))
  top_sigs_reordered <- top_sigs[order(as.numeric(gsub("[^0-9]", "", top_sigs)))]
  
  ## by condition (non-top go to "others")
  top_sSNV_sig_contri_cond <- sSNV_sig_contri_AMG %>% 
    dplyr::select(starts_with("SBS")) %>% 
    mutate(others = rowSums(.[, !(colnames(.) %in% top_sigs)])) %>% t() %>% as.data.frame() %>% 
    mutate(Control = rowMeans(.[, ctrl_range_AMG]), IHD = rowMeans(.[, dis_range_AMG])) %>% 
    dplyr::select(all_of(ctrl_name), all_of(dis_name)) %>% 
    mutate(Control_Pct = Control / sum(Control) * 100, IHD_Pct = IHD / sum(IHD) * 100) %>% 
    filter(rownames(.) %in% c(top_sigs, "others"))
  
  ## by condition keep all
  sSNV_sig_contri_cond_save <- sSNV_sig_contri_AMG %>% 
    dplyr::select(starts_with("SBS")) %>% 
    mutate(others = rowSums(.[, !(colnames(.) %in% colnames(sSNV_sig_contri_AMG))])) %>% t() %>% as.data.frame() %>% 
    mutate(Control = rowMeans(.[, ctrl_range_AMG]), IHD = rowMeans(.[, dis_range_AMG])) %>% 
    dplyr::select(all_of(ctrl_name), all_of(dis_name)) %>% 
    mutate(Control_Pct = Control / sum(Control) * 100, IHD_Pct = IHD / sum(IHD) * 100) %>% 
    filter(rownames(.) %in% c(colnames(sSNV_sig_contri_AMG), "others"))
  
  p_top_sSNV_sig_contri_cond <- plot_contribution(top_sSNV_sig_contri_cond[, 1:2], coord_flip = FALSE, mode = plot_type) + 
    scale_fill_manual(name = "Signature", values = sig_colors) + geom_bar(stat = "identity", colour = NA, linewidth = 0) + theme_linedraw() + 
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12))

  ## by donor 
  top_sSNV_sig_contri_bydonor <- sSNV_sig_contri_AMG_bydonor %>% 
    dplyr::select(starts_with("SBS")) %>% 
    mutate(others = rowSums(.[, !(colnames(.) %in% top_sigs)])) %>% t() %>% as.data.frame() %>% 
    filter(rownames(.) %in% c(top_sigs, "others")) %>% as.matrix()
  
  p_top_sSNV_sig_bydonor_barplot <- plot_contribution(top_sSNV_sig_contri_bydonor, coord_flip = FALSE, mode = plot_type) + 
    scale_x_discrete(limits = Case_ID_list_AMG) + guides(fill = guide_legend(ncol = 2)) + scale_fill_manual(name = "Signature", values = sig_colors) + 
    geom_bar(stat = "identity", colour = NA, linewidth = 0) + theme_linedraw() + 
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.0))
  
  ## by cell
  top_sSNV_sig_contri_bycell <- sSNV_sig_contri_AMG %>% 
    dplyr::select(starts_with("SBS")) %>% 
    mutate(others = rowSums(.[, !(colnames(.) %in% top_sigs)])) %>% t() %>% as.data.frame() %>% 
    filter(rownames(.) %in% c(top_sigs, "others")) %>% as.matrix()
  
  p_top_sSNV_sig_bycell_barplot <- plot_contribution(top_sSNV_sig_contri_bycell, coord_flip = FALSE, mode = plot_type) + 
    scale_x_discrete(limits = Cell_ID_list_AMG) + guides(fill = guide_legend(ncol = 2)) + scale_fill_manual(name = "Signature", values = sig_colors) + 
    geom_bar(stat = "identity", colour = NA, linewidth = 0) + theme_linedraw() + 
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.0))
  
  if (plot_type == "relative") {
    ## by condition 
    write.csv(sSNV_sig_contri_cond_save, file = paste0(table_dir, "/4-sSNV_sig_contri_", plot_type, "_by_cond_AMG_PTA_all.csv"))
    write.csv(top_sSNV_sig_contri_cond, file = paste0(table_dir, "/4-sSNV_sig_contri_", plot_type, "_by_cond_AMG_PTA_top.csv"))
    ggsave(paste0(suppl_figure_dir, "/4-sSNV_top_sigs_contri_", plot_type, ".pdf"), plot = p_top_sSNV_sig_contri_cond, width = 3.5, height = 7, dpi = 300)
    ## by cell
    ggsave(paste0(other_figure_dir, "/4-sSNV_top_sigs_contri_bycell_", plot_type, "_barplot.pdf"), plot = p_top_sSNV_sig_bycell_barplot, width = 12, height = 8, dpi = 300)
    ## by donor
    ggsave(paste0(other_figure_dir, "/4-sSNV_top_sigs_contri_bydonor_", plot_type, "_barplot.pdf"), plot = p_top_sSNV_sig_bydonor_barplot, width = 12, height = 8, dpi = 300)
  }
  else if (plot_type == "absolute") {
    ## by condition
    write.csv(sSNV_sig_contri_cond_save, file = paste0(table_dir, "/4-sSNV_sig_contri_", plot_type, "_by_cond_AMG_PTA_all.csv"))
    write.csv(top_sSNV_sig_contri_cond, file = paste0(table_dir, "/4-sSNV_sig_contri_", plot_type, "_by_cond_AMG_PTA_top.csv"))
    ggsave(paste0(main_figure_dir, "/4-sSNV_top_sigs_contri_", plot_type, ".pdf"), plot = p_top_sSNV_sig_contri_cond, width = 4, height = 7, dpi = 300)
    ## by cell
    ggsave(paste0(other_figure_dir, "/4-sSNV_top_sigs_contri_bycell_", plot_type, "_barplot.pdf"), plot = p_top_sSNV_sig_bycell_barplot, width = 12, height = 8, dpi = 300)
    ## by donor
    ggsave(paste0(other_figure_dir, "/4-sSNV_top_sigs_contri_bydonor_", plot_type, "_barplot.pdf"), plot = p_top_sSNV_sig_bydonor_barplot, width = 12, height = 8, dpi = 300)
  }
}

##############################################################################
### add metadata and remove the signatures with all zeros in one Condition
sSNV_sig_contri_est <- sSNV_sig_contri_est %>% 
  dplyr::select(starts_with("SBS")) %>% 
  mutate(Age = metadata_df$Age, Condition = relevel(factor(metadata_df$Condition), ref = "Control"), 
         Gender = metadata_df$Gender, Case_ID = metadata_df$Case_ID, Color = sSNV_SCAN2_df$Color, Outline_Color = sSNV_SCAN2_df$Outline_Color) %>% 
  {
  # sSNV_sig_contri_est_control <- .[.$Condition == ctrl_name, ]
  # zero_col_control <- which(colSums(sSNV_sig_contri_est_control[, 1:(ncol(sSNV_sig_contri_est_control) - 6)]) == 0)
  # . <- .[, -zero_col_control]
  sSNV_sig_contri_est_disease <- .[.$Condition == dis_name, ]
  zero_col_disease <- which(colSums(sSNV_sig_contri_est_disease[, 1:(ncol(sSNV_sig_contri_est_disease) - 6)]) == 0)
  .[, -zero_col_disease]}

refined_selected_sigs_list <- colnames(sSNV_sig_contri_est)[1 : (ncol(sSNV_sig_contri_est) - 5)]
sig_digits <- 2

##### Mix effects model & Age matched mix effects model
pdf(paste0(other_figure_dir, "/4-lmer_top_sigs_all.pdf"), width = 9, height = 6)
  for (i in c(1 : (ncol(sSNV_sig_contri_est) - 6))){
    cat("index", i, ":", refined_selected_sigs_list[i])
    ## burden_df: the input matrix with columns: burden, Condition, Case_ID, Age, Color
    burden_df <- sSNV_sig_contri_est[c(refined_selected_sigs_list[i], "Age", "Condition", "Case_ID", "Color", "Outline_Color")] %>% 
      setNames(c("snv.rate.per.gb", "Age", "Condition", "Case_ID", "Color", "Outline_Color")) %>% 
      # filter(snv.rate.per.gb < 5000) %>% 
      mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% 
      mutate(Condition = relevel(Condition, ref = ctrl_name)) %>% 
      group_by(Condition) %>% arrange(Age, .by_group = TRUE)
    burden_df_ctrl <- burden_df %>% filter(Condition == "Control")
    burden_df_dis <- burden_df %>% filter(Condition == "IHD")
    
    generate_line_data <- function(df_ctrl, coef_vec) {
      x <- range(df_ctrl$Age)
      y <- x * coef_vec[2] + coef_vec[1]
      data.frame(ctrl_fitting_x = x, ctrl_fitting_y = y)
    }
    
    generate_y_breaks <- function(y_data) {
      y_min <- min(y_data, na.rm = TRUE)
      y_max <- max(y_data, na.rm = TRUE)
      range_span <- y_max - y_min
      step_size <- signif(range_span / 4, 1)
      breaks <- seq(floor(y_min / step_size) * step_size, ceiling(y_max / step_size) * step_size, by = step_size)
      return(breaks)
    }
    
    ## Check sums and choose model
    check_sums <- sum(burden_df_ctrl$snv.rate.per.gb) > 10 & sum(burden_df_dis$snv.rate.per.gb) > 10
    burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + (1|Case_ID), burden_df, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
    
    if (check_sums) {
      cat(", Normal\n")
      burden_age_model_ctrl <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_ctrl, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
      coef_ctrl <- fixef(burden_age_model_ctrl)
    } else {
      cat(", Low number\n")
      burden_age_model_ctrl <- lm(snv.rate.per.gb ~ Age, burden_df_ctrl)
      coef_ctrl <- coef(burden_age_model_ctrl)
    }
    
    aging_rate <- format(round(coef_ctrl[2], digits = sig_digits), nsmall = sig_digits)
    ## compute confident interval
    # CI_age_model <- confint(burden_age_model, method = "Wald", oldNames = FALSE)
    # CI_age_model_ctrl <- confint(burden_age_model_ctrl, method = "Wald", oldNames = FALSE)
    
    ## check p value and R^2
    anova_pval <- anova(burden_age_model)$"Pr(>F)"[2]
    anova_pval_print <- formatC(signif(anova_pval, digits = sig_digits), digits = sig_digits, format="fg", flag="#")
    anova_pval_ctrl <- anova(burden_age_model_ctrl)$"Pr(>F)"[1]
    anova_pval_print_ctrl <- formatC(signif(anova_pval_ctrl, digits = sig_digits), digits = sig_digits, format="e", flag="#")
    
    geom_line_data <- generate_line_data(burden_df_ctrl, coef_ctrl)
    legend_data <- burden_df[c(1, nrow(burden_df)), c("Age", "snv.rate.per.gb", "Condition")]
    
    yticks_plot <- generate_y_breaks(burden_df$snv.rate.per.gb)
    ylim_plot <- range(yticks_plot, na.rm = TRUE)
    p_SNV_burden_lme <- ggplot(burden_df, aes(x = Age, y = snv.rate.per.gb)) + 
      geom_point(pch = 21, data = legend_data, aes(x = Age, y = snv.rate.per.gb, color = Condition, fill = Condition), size = 5) + 
      geom_point(pch = 21, color = burden_df$Outline_Color, fill = burden_df$Color, size = 5) + 
      geom_line(aes(x = ctrl_fitting_x, y = ctrl_fitting_y), color = "dodgerblue2", linetype = "22", data = geom_line_data) + 
      geom_hline(yintercept = 0, color = "gray40", linetype = "22") + 
      annotate("text", size = 6, x = 0, y = 0.95 * ylim_plot[2], 
               label = paste("Aging effect:", aging_rate, "sSNVs/(GB·year),", "P =", anova_pval_print_ctrl), hjust = 0) +
      annotate("text", size = 6, x = 0, y = 0.87 * ylim_plot[2], 
               label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB,", "P =", anova_pval_print), hjust = 0) +
      scale_color_manual(values = c("IHD" = "black", "Control" = "black"), guide = "legend") + 
      scale_fill_manual(values = c("Control" = ctrl_dis_color[1], "IHD" = ctrl_dis_color[2]), guide = "legend") + 
      scale_y_continuous(limits = c(min(min(geom_line_data[, 2]), 0), ylim_plot[2]), breaks = yticks_plot, labels = yticks_plot) + 
      labs(x = "Age (years)", y = paste0(refined_selected_sigs_list[i], " Contribution \n (sSNV rate per GB)")) + theme_linedraw() + 
      theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
    print(p_SNV_burden_lme)
    if (anova_pval < 0.05 & (mean(burden_df_ctrl$snv.rate.per.gb) > 10 | mean(burden_df_dis$snv.rate.per.gb) > 10) | 
        refined_selected_sigs_list[i] %in% c("SBS1", "SBS4", "SBS5", "SBS92")) {
      ggsave(paste0(main_figure_dir, "/4-", refined_selected_sigs_list[i], "_SNV_burden_lme.pdf"),
             plot = p_SNV_burden_lme, width = 9, height = 6, dpi = 600)
    }
  }
dev.off()
