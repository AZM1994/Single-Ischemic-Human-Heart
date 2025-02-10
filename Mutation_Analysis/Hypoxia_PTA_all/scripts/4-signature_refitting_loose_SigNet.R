################################################################################
########################## Signature refitting loose ###########################

##### get COSMIC signature spectrum
# signatures = get_known_signatures()
COSMIC_v3.4_sigs <- read.table("main/COSMIC/COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE) %>% {rownames(.) <- .$Type; .} %>% dplyr::select(-Type) %>% as.matrix()
COSMIC_v3.3.1_sigs <- read.table("main/COSMIC/COSMIC_v3.3.1_SBS_GRCh37.txt", header = TRUE) %>% {rownames(.) <- .$Type; .} %>% dplyr::select(-Type) %>% as.matrix()

SigNet_contri <- read.csv(paste0(project_dir, "/data/SigNet/PTA_all_est/weight_guesses.csv"), row.names = 1) %>% t()
# SigNet_contri <- read.csv(paste0(project_dir, "/data/SigNet/PTA_all_raw/weight_guesses.csv"), row.names = 1) %>% t()
SigNet_contri_AMG <- SigNet_contri[ , Cell_ID_list_AMG]

##### the contribution bar plot
pdf(paste0(other_figure_dir, "/4-Sig_refitting_SigNet.pdf"), width = 16, height = 12)
  ## plot absolute signature contribution: barplot
  print(plot_contribution(SigNet_contri_AMG, coord_flip = FALSE, mode = "absolute") + scale_x_discrete(limits = Cell_ID_list_AMG) + guides(fill = guide_legend(ncol = 2))) + 
    theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0))

  ## plot absolute signature contribution: heatmap
  print(plot_contribution_heatmap(SigNet_contri_AMG, cluster_sigs = FALSE, cluster_samples = FALSE) + scale_y_discrete(limits = rev(Cell_ID_list_AMG)) + 
    scale_fill_gradient(low = "white", high = "blue") + geom_tile(colour = "black") + theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)))

  ## by conditions
  SigNet_contri_AMG_cond <- SigNet_contri_AMG %>% as.data.frame() %>% mutate(Control = rowMeans(.[, ctrl_range_AMG]), IHD = rowMeans(.[, dis_range_AMG]))
  
  ## plot absolute signature contribution: barplot
  print(plot_contribution(SigNet_contri_AMG_cond[ , c(ctrl_name, dis_name)], coord_flip = FALSE, mode = "absolute") + guides(fill = guide_legend(ncol = 2))) + 
    theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0))
  
  ## plot absolute signature contribution: heatmap
  print(plot_contribution_heatmap(as.matrix(SigNet_contri_AMG_cond[ , c(ctrl_name, dis_name)]), cluster_sigs = FALSE, cluster_samples = FALSE) + 
    scale_fill_gradient(low = "white", high = "blue") + geom_tile(colour = "black") + theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)))
dev.off()