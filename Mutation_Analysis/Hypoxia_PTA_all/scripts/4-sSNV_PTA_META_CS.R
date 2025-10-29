library(ComplexHeatmap)
################################################################################
## similarity between the PTA and META-CS control and IHD top signatures
ctrl_color_palette <- colorRampPalette(c("skyblue1","dodgerblue4"))
dis_color_palette <- colorRampPalette(c("pink1","firebrick3"))
ctrl_dis_color <- c(ctrl_color_palette(9)[7], dis_color_palette(4)[3])
sig_colors <- c(SBS6 = "#acdb88", SBS14 = "#fca85c", SBS15 = "#f2b8c0", SBS20 = "#cb77b2", SBS21 = "#a4e4e0", SBS26 = "#c6a589", SBS44 = "#8882bd")

## compare PTA and METS-CS raw relative signature contribution (top15)
PTA_contri_raw <- read.csv(file = paste0(table_dir, "/4-sSNV_sig_contri_relative_by_cond_AMG_PTA_all.csv"), row.names = 1) %>% 
  dplyr::select(ctrl_name, dis_name) %>% t() %>% as.data.frame() %>% { `rownames<-`(., paste0("PTA_", rownames(.))) }
PTA_sig_contri_rel_AMG <- read.csv(file = paste0(table_dir, "/4-sSNV_sig_contri_rel_by_cell_PTA.csv"), row.names = 1) %>% filter(Age >= 30 & Age <= 70)
top_sig_contri_PTA <- sort(colSums(PTA_sig_contri_rel_AMG %>% dplyr::select(starts_with("SBS"))), decreasing = TRUE)[1 : 15]
top_sigs_PTA <- row.names(data.frame(top_sig_contri_PTA))

METS_CS_contri_raw <- read.csv(file = paste0(table_dir, "/4-sSNV_sig_contri_relative_by_cond_AMG_META_CS_all.csv"), row.names = 1) %>%
  dplyr::select(ctrl_name, dis_name) %>% t() %>% as.data.frame() %>% { `rownames<-`(., paste0("META_CS_", rownames(.))) }
METS_CS_sig_contri_rel_AMG <- read.csv(file = paste0(table_dir, "/4-sSNV_sig_contri_rel_by_cell_META_CS.csv"), row.names = 1) %>% filter(Age >= 30 & Age <= 70)
top_sig_contri_META_CS <- sort(colSums(METS_CS_sig_contri_rel_AMG %>% dplyr::select(starts_with("SBS"))), decreasing = TRUE)[1 : 15]
top_sigs_META_CS <- row.names(data.frame(top_sig_contri_META_CS))

combind_contri_raw <- rbind(PTA_contri_raw, METS_CS_contri_raw) %>% dplyr::select(union(top_sigs_PTA, top_sigs_META_CS)) %>% as.matrix()
sig_names <- colnames(combind_contri_raw)
sig_df <- data.frame(sig = sig_names) %>% mutate(num = as.numeric(gsub("[^0-9]", "", sig)), suffix = gsub("^[^0-9]*[0-9]+", "", sig))
ordered_names <- sig_df %>% arrange(num, suffix) %>% pull(sig)
combind_contri_raw_ordered <- combind_contri_raw[, ordered_names]
combind_contri_raw_log <- log10(combind_contri_raw_ordered + 1e-3)
breaks <- c(1e-3, 0.01, 0.05, 0.2, 0.7)
colors <- c("grey80", "#fff7b2", "#a7e5c6", "#66c2a5", "#3a8f70")
col_fun <- colorRamp2(log10(breaks), colors)
heatmap_legend_param <- list(title = "Contribution", at = log10(c(1e-3, 0.01, 0.05, 0.2, 0.7)), labels = c("0", "0.01", "0.05", "0.2", "0.7"))
dend = as.dendrogram(hclust(dist(combind_contri_raw_log)))
desired_order <- c("PTA_Control", "META_CS_Control", "PTA_IHD", "META_CS_IHD")
dend <- dendextend::rotate(dend, order = desired_order)
dend = dendextend::color_branches(dend, k = 2, col = ctrl_dis_color)
pdf(paste0(suppl_figure_dir, "/4-top_sig_contri_heatmap.pdf"), width = 12, height = 3.5)
  Heatmap(combind_contri_raw_log, col = col_fun, cluster_columns = FALSE,
          heatmap_legend_param = heatmap_legend_param, column_names_rot = 90, cluster_rows = dend, row_split = 2,
          left_annotation = rowAnnotation(
            foo = anno_block(gp = gpar(fill = ctrl_dis_color), labels = c("Control", "IHD"), labels_gp = gpar(col = "white", fontsize = 10))))
dev.off()

################################################################################
## compare PTA and METS-CS MMR-relative signature contribution
sigs_reordered <- c("SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44")
PTA_sig_contri_rel_AMG_MMR <- PTA_sig_contri_rel_AMG %>% 
  mutate(MMR_sig = rowSums(across(all_of(sigs_reordered)), na.rm = TRUE), Technique = "PTA") %>% 
  select(Case_ID, Condition, Age, Technique, MMR_sig, all_of(sigs_reordered))

METS_CS_sig_contri_rel_AMG_MMR <- METS_CS_sig_contri_rel_AMG %>% 
  mutate(MMR_sig = rowSums(across(all_of(sigs_reordered)), na.rm = TRUE), Technique = "META-CS") %>% 
  select(Case_ID, Condition, Age, Technique, MMR_sig, all_of(sigs_reordered))

MMR_long <- bind_rows(PTA_sig_contri_rel_AMG_MMR, METS_CS_sig_contri_rel_AMG_MMR) %>% 
  group_by(Technique, Condition) %>% 
  summarise(across(all_of(sigs_reordered), mean, na.rm = TRUE), .groups = "drop") %>% 
  pivot_longer(cols = starts_with("SBS"), names_to = "Signature", values_to = "Mean_contribution") %>% 
  mutate(Signature = factor(Signature, levels = sigs_reordered))

p_MMR <- ggplot(MMR_long, aes(x = Condition, y = Mean_contribution, fill = Signature)) + 
  geom_col() + facet_wrap(~ Technique, scales = "fixed") + theme_linedraw() + scale_fill_manual(name = "Signature", values = sig_colors) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "", y = "MMR-related signature contribution")
ggsave(paste0(suppl_figure_dir, "/4-sSNV_MMR_Contribution_PTA_META_CS_AMG_by_Sig.pdf"), plot = p_MMR, width = 6, height = 6, dpi = 600)

## overall by condition and tech
MMR_total <- bind_rows(PTA_sig_contri_rel_AMG_MMR, METS_CS_sig_contri_rel_AMG_MMR) %>% 
  group_by(Technique, Condition) %>% 
  summarise(mean_MMR = mean(MMR_sig, na.rm = TRUE), .groups = "drop")

p_MMR_total <- ggplot(MMR_total, aes(x = Condition, y = mean_MMR, fill = Condition)) + 
  geom_col() + facet_wrap(~ Technique, scales = "fixed") + 
  scale_fill_manual(values = c("Control" = ctrl_dis_color[1], "IHD" = ctrl_dis_color[2])) + 
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "", y = "MMR-related signature contribution")
ggsave(paste0(suppl_figure_dir, "/4-sSNV_MMR_Contribution_PTA_META_CS_AMG.pdf"), plot = p_MMR_total, width = 6, height = 6, dpi = 600)

################################################################################
## similarity between the PTA and META-CS control and IHD mutational profiles
library(proxy)
library(circlize)
PTA_mut_mat_raw_AMG_cond <- read.csv(file = paste0(table_dir, "/3-mut_mat_raw_AMG_cond_PTA_all.csv"), row.names = 1) %>% 
  select(Control, IHD) %>% apply(2, function(x) x / sum(x, na.rm = TRUE)) %>% as.data.frame()
META_CS_mut_mat_raw_AMG_cond <- read.csv(file = paste0(table_dir, "/3-mut_mat_raw_AMG_cond_META_CS.csv"), row.names = 1) %>% 
  select(Control, IHD) %>% apply(2, function(x) x / sum(x, na.rm = TRUE)) %>% as.data.frame()

PTA_Control <- PTA_mut_mat_raw_AMG_cond$Control
PTA_IHD <- PTA_mut_mat_raw_AMG_cond$IHD
META_CS_Control <- META_CS_mut_mat_raw_AMG_cond$Control
META_CS_IHD <- META_CS_mut_mat_raw_AMG_cond$IHD
all_spectra <- cbind(PTA_Control = PTA_Control, PTA_IHD = PTA_IHD, META_Control = META_CS_Control, META_IHD = META_CS_IHD)
cosine_mat <- as.matrix(simil(t(all_spectra), method = "cosine"))
diag(cosine_mat) <- 1

col_fun <- colorRamp2(c(0.7, 1.0), c("white", "dodgerblue4"))
Heatmap(cosine_mat, name = "CosSimilarity", cluster_rows = TRUE, cluster_columns = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.2f", cosine_mat[i, j]), x, y, gp = gpar(fontsize = 20))},
        col = col_fun, column_title = "Cosine similarity of 96-trinucleotide spectra")
