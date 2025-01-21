################################################################################
########################## Signature refitting loose ###########################

##### get COSMIC signature spectrum
# signatures = get_known_signatures()
# COSMIC_v3.4_sigs <- read.table("main/COSMIC/COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE) %>% write.xlsx(file = "COSMIC_v3.4_SBS_GRCh37.xlsx")
# COSMIC_v3.4_sigs <- read.table("main/COSMIC/COSMIC_v3.3.1_SBS_GRCh37.txt", header = TRUE) %>% write.xlsx(file = "COSMIC_v3.3.1_SBS_GRCh37.xlsx")

COSMIC_v3.4_sigs <- read.table("main/COSMIC/COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE) %>%
  {rownames(.) <- .$Type; .} %>%
  select(-Type) %>%
  as.matrix()

COSMIC_v3.3.1_sigs <- read.table("main/COSMIC/COSMIC_v3.3.1_SBS_GRCh37.txt", header = TRUE) %>% 
  {rownames(.) <- .$Type; .} %>%
  select(-Type) %>% 
  as.matrix()

SigNet_contri <- read.csv(paste0(main_figure_dir, "/SigNet/signature_weights/weight_guesses.csv"), row.names = 1) %>% t()
SigNet_contri_meta <- as.data.frame(t(SigNet_contri)) %>%
  mutate(Case_ID_decile = rownames(.)) %>%
  mutate(decile = sub(".*(decile.*)", "\\1", Case_ID_decile)) %>%
  mutate(Case_ID = sub("(.*)_decile.*", "\\1", Case_ID_decile))

SigNet_contri_by_case <- SigNet_contri_meta %>% group_by(decile) %>% 
  summarise(across(where(is.numeric), sum)) %>% select(-1) %>% t()

SigNet_contri_age_match <- SigNet_contri[ , Cell_ID_list_age_match]

##### the contribution bar plot
pdf(paste0(other_figure_dir, "/4-Sig_refitting_SigNet.pdf"), width = 16, height = 12)
  ## plot absolute signature contribution: barplot
  print(plot_contribution(SigNet_contri_age_match, coord_flip = FALSE, mode = "absolute") +
    theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)) +
    scale_x_discrete(limits = Cell_ID_list_age_match) + guides(fill=guide_legend(ncol=2)))

  ## plot absolute signature contribution: heatmap
  print(plot_contribution_heatmap(SigNet_contri_age_match, cluster_sigs = FALSE, cluster_samples = FALSE) + 
    scale_y_discrete(limits = rev(Cell_ID_list_age_match)) + 
    scale_fill_gradient(low = "white", high = "blue") + geom_tile(colour = "black")  + 
    theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)))

  ############## by conditions
  SigNet_contri_age_match_conditional <- as.data.frame(SigNet_contri_age_match) %>% 
    mutate(Control = rowMeans(SigNet_contri_age_match[ , control_range_age_match])) %>% 
    mutate(IHD = rowMeans(SigNet_contri_age_match[ , disease_range_age_match]))
  
  ## plot absolute signature contribution: barplot
  print(plot_contribution(SigNet_contri_age_match_conditional[ , c(control_name, disease_name)], coord_flip = FALSE, mode = "absolute") +
    theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)) + 
    guides(fill=guide_legend(ncol=2)))
  
  ## plot absolute signature contribution: heatmap
  print(plot_contribution_heatmap(as.matrix(SigNet_contri_age_match_conditional[ , c(control_name, disease_name)]), cluster_sigs = FALSE, cluster_samples = FALSE) + 
    scale_fill_gradient(low = "white", high = "blue") + geom_tile(colour = "black") + 
    theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0)))

dev.off()
