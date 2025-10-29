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
SigNet_contri <- read.csv(paste0(project_dir, "/data/SigNet/2n/DS/weight_guesses.csv"), row.names = 1) %>% as.data.frame()

## calculate relative and absolute contribution by cell
sSNV_sig_contri_rel <- (1 / rowSums(SigNet_contri) * SigNet_contri) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Cell_ID") %>% 
  inner_join(META_CS_burden_df %>% select(Cell_ID, Case_ID, Condition, Age), by = "Cell_ID") %>% 
  relocate(Cell_ID, Case_ID, Condition, Age, .before = 1) %>% 
  column_to_rownames(var = "Cell_ID")
# sSNV_sig_contri_est <- sSNV_sig_contri_rel %>% mutate(across(starts_with("SBS"), ~ .x * snv.rate.per.gb))
write.csv(sSNV_sig_contri_rel, file = paste0(table_dir, "/4-sSNV_sig_contri_rel_by_cell_META_CS.csv"))
# write.csv(sSNV_sig_contri_est, file = paste0(table_dir, "/4-sSNV_sig_contri_est_by_cell_META_CS.csv"))

## calculate relative and absolute contribution by donor
sSNV_sig_contri_rel_bydonor <- sSNV_sig_contri_rel %>% 
  mutate(Case_ID = factor(Case_ID, levels = unique(Case_ID))) %>% 
  group_by(Case_ID) %>% 
  summarise(Condition = first(Condition), Age = first(Age), across(starts_with("SBS"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>% 
  arrange(Case_ID) %>% column_to_rownames(var = "Case_ID")
# sSNV_sig_contri_est_bydonor <- sSNV_sig_contri_rel_bydonor %>% mutate(across(starts_with("SBS"), ~ .x * snv.rate.per.gb))
write.csv(sSNV_sig_contri_rel_bydonor, file = paste0(table_dir, "/4-sSNV_sig_contri_rel_by_donor_META_CS.csv"))
# write.csv(sSNV_sig_contri_est_bydonor, file = paste0(table_dir, "/4-sSNV_sig_contri_est_by_donor_META_CS.csv"))

##### determine top contribution signatures in age-matched samples
##### plot contribution by cell, by donor, by condition (age-matched)
plot_type_list <- c("relative")
for (plot_type in plot_type_list) {
  if (plot_type == "relative") {
    sSNV_sig_contri_AMG <- sSNV_sig_contri_rel[Cell_ID_list_AMG, ]
    sSNV_sig_contri_AMG_bydonor <- sSNV_sig_contri_rel_bydonor[Case_ID_list_AMG, ]
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
    write.csv(sSNV_sig_contri_cond_save, file = paste0(table_dir, "/4-sSNV_sig_contri_", plot_type, "_by_cond_AMG_META_CS_all.csv"))
    write.csv(top_sSNV_sig_contri_cond, file = paste0(table_dir, "/4-sSNV_sig_contri_", plot_type, "_by_cond_AMG_META_CS_top.csv"))
    ggsave(paste0(main_figure_dir, "/4-sSNV_top_sigs_contri_", plot_type, ".pdf"), plot = p_top_sSNV_sig_contri_cond, width = 3.5, height = 7, dpi = 300)
    ## by cell
    ggsave(paste0(main_figure_dir, "/4-sSNV_top_sigs_contri_bycell_", plot_type, "_barplot.pdf"), plot = p_top_sSNV_sig_bycell_barplot, width = 12, height = 8, dpi = 300)
    ## by donor
    ggsave(paste0(main_figure_dir, "/4-sSNV_top_sigs_contri_bydonor_", plot_type, "_barplot.pdf"), plot = p_top_sSNV_sig_bydonor_barplot, width = 12, height = 8, dpi = 300)
  }
}
