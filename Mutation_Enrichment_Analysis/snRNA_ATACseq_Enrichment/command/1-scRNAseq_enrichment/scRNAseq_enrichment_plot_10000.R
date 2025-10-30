library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)
library(DescTools)

setwd("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration/Mutation_Enrichment_Analysis/snRNA_ATACseq_Enrichment")
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])
group_num <- 8
permutation_num <- 10000

##### read in metadata
Hypoxia_PTA_all_SCAN2 <- readRDS("./data/1-sSNV_SCAN2_df_filtered.rds") %>% dplyr::select(c("Cell_ID", "Age", "Gender", "Case_ID", "Condition", "snv.burden", "snv.rate.per.gb"))
Hypoxia_PTA_all_SCAN2_collapsed <- Hypoxia_PTA_all_SCAN2 %>% distinct(Case_ID, .keep_all = TRUE)
# Cell_ID_list <- Hypoxia_PTA_all_SCAN2$Cell_ID
Case_ID_list <- Hypoxia_PTA_all_SCAN2_collapsed$Case_ID
Condition_list <- unique(Hypoxia_PTA_all_SCAN2$Condition)

analysis_ID <- paste0(permutation_num, "_perms_", group_num, "G")
fig_save_dir <- paste0("./results/1-scRNAseq_analysis_by_case/", analysis_ID, "_final/")
suppressWarnings(dir.create(fig_save_dir, recursive = TRUE))

##### read in all mutmat
mut_mat_dir <- paste0("./results/1-scRNAseq_analysis_by_case/", analysis_ID, "_mut_mat/")
mutation_mut_mat_all <- c()
permutation_mut_mat_all <- c()
for (Case_ID_temp in Case_ID_list){
  mutation_mut_mat_temp <- readRDS(paste0(mut_mat_dir, Case_ID_temp, "/mutation_mut_mat_by_case_decile_", Case_ID_temp, ".rds"))
  mutation_mut_mat_all <- rbind(mutation_mut_mat_all, mutation_mut_mat_temp)
  permutation_mut_mat_temp <- readRDS(paste0(mut_mat_dir, Case_ID_temp, "/permutation_mut_mat_by_case_decile_", Case_ID_temp, ".rds"))
  permutation_mut_mat_all <- rbind(permutation_mut_mat_all, permutation_mut_mat_temp)
}

## normalize mut_mat to est snv.burden (per Case_ID)
mutation_avg_by_case <- Hypoxia_PTA_all_SCAN2 %>% group_by(Case_ID) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% 
  arrange(match(Case_ID, Hypoxia_PTA_all_SCAN2_collapsed$Case_ID)) %>% 
  mutate(Condition = Hypoxia_PTA_all_SCAN2_collapsed$Condition)

## calculate burden for each Case_ID and decile
mutation_mut_mat_all_burden_decile <- mutation_mut_mat_all %>% 
  mutate(row_sum_decile = rowSums(dplyr::select(., 1:96))) %>% 
  group_by(Case_ID) %>% 
  mutate(percentage_decile = row_sum_decile / sum(row_sum_decile) * 100) %>% 
  ungroup() %>% 
  mutate(burden_individual = rep(mutation_avg_by_case$snv.burden, each = group_num)) %>% 
  mutate(burden_decile = percentage_decile * burden_individual / 100)

permutation_mut_mat_all_burden_decile <- permutation_mut_mat_all %>%
  mutate(row_sum_decile = rowSums(dplyr::select(., 1:96))) %>% 
  group_by(perm.id, Case_ID) %>%
  mutate(percentage_decile = row_sum_decile / sum(row_sum_decile) * 100) %>%
  ungroup() %>%
  mutate(burden_individual = rep(mutation_avg_by_case$snv.burden, each = group_num * permutation_num)) %>%
  mutate(burden_decile = percentage_decile * burden_individual / 100)

mutation_mut_mat_all_raw <- mutation_mut_mat_all[, 1:96]
permutation_mut_mat_all_raw <- permutation_mut_mat_all[, 1:96]
mutation_mut_mat_all_est <- 1 / rowSums(mutation_mut_mat_all_raw) * mutation_mut_mat_all_raw * mutation_mut_mat_all_burden_decile$burden_decile
permutation_mut_mat_all_est <- 1 / rowSums(permutation_mut_mat_all_raw) * permutation_mut_mat_all_raw * permutation_mut_mat_all_burden_decile$burden_decile

## replace all 0s and NaNs in mutation and permutation results by 1
zero_sum_rows_mut <- rowSums(mutation_mut_mat_all_est) %in% c(0, NaN)
sum(zero_sum_rows_mut)
mutation_mut_mat_all_est[zero_sum_rows_mut, ] <- 1
mutation_mut_mat_all_raw[zero_sum_rows_mut, ] <- 1

zero_sum_rows_permut <- rowSums(permutation_mut_mat_all_est) %in% c(0, NaN)
sum(zero_sum_rows_permut)
permutation_mut_mat_all_est[zero_sum_rows_permut, ] <- 1
permutation_mut_mat_all_raw[zero_sum_rows_permut, ] <- 1

# write.csv(mutation_mut_mat_all_raw, paste0(fig_save_dir, "mutation_mut_mat_all_raw.csv"))
# write.csv(permutation_mut_mat_all_raw, paste0(fig_save_dir, "permutation_mut_mat_all_raw.csv"))
# write.csv(mutation_mut_mat_all_est, paste0(fig_save_dir, "mutation_mut_mat_all_est.csv"))
# write.csv(permutation_mut_mat_all_est, paste0(fig_save_dir, "permutation_mut_mat_all_est.csv"))

# Number of rows per split
# n_split <- ceiling(nrow(permutation_mut_mat_all_raw) / 10)
# for (i in 1:10) {
#   start_row <- (i - 1) * n_split + 1
#   end_row <- min(i * n_split, nrow(permutation_mut_mat_all_raw))
#   chunk <- permutation_mut_mat_all_raw[start_row:end_row, ]
#   suppressWarnings(dir.create(paste0(fig_save_dir, "1-permutation_mut_mat_all_raw_split"), recursive = TRUE))
#   write.csv(chunk, file = paste0(fig_save_dir, "1-permutation_mut_mat_all_raw_split/permutation_mut_mat_all_raw_part", i, ".csv"))
#   }

# n_split <- ceiling(nrow(permutation_mut_mat_all_est) / 10)
# for (i in 1:10) {
#   start_row <- (i - 1) * n_split + 1
#   end_row <- min(i * n_split, nrow(permutation_mut_mat_all_est))
#   chunk <- permutation_mut_mat_all_est[start_row:end_row, ]
#   suppressWarnings(dir.create(paste0(fig_save_dir, "1-permutation_mut_mat_all_est_split"), recursive = TRUE))
#   write.csv(chunk, file = paste0(fig_save_dir, "1-permutation_mut_mat_all_est_split/permutation_mut_mat_all_est_part", i, ".csv"))
# }

##### mutation mut_mat summary
mutation_mut_mat_summary <- mutation_mut_mat_all %>% 
  mutate(mut_sum = rowSums(dplyr::select(., 1:96))) %>% 
  dplyr::select(c("Case_ID", "mut_sum", "decile", "Condition")) %>% 
  group_by(Case_ID, decile) %>% 
  summarise(Value = mut_sum, .groups = "drop") %>% 
  pivot_wider(names_from = decile, values_from = Value) %>% 
  mutate(mut_num_percase = rowSums(dplyr::select(., 2:(2 + group_num - 1)))) %>% 
  arrange(match(Case_ID, Hypoxia_PTA_all_SCAN2$Case_ID)) %>% 
  mutate(Condition = Hypoxia_PTA_all_SCAN2_collapsed$Condition, 
         Case_ID = Hypoxia_PTA_all_SCAN2_collapsed$Case_ID, 
         Age = Hypoxia_PTA_all_SCAN2_collapsed$Age)

mutation_mut_mat_summary_ctrl <- mutation_mut_mat_summary %>% filter(Condition == "Control") %>% 
  summarise(across(2:(2 + group_num), sum, .names = "{.col}")) %>%
  bind_rows(mutation_mut_mat_summary[mutation_mut_mat_summary$Condition == "Control", ], .)

mutation_mut_mat_summary_dis <- mutation_mut_mat_summary %>% filter(Condition == "IHD") %>% 
  summarise(across(2:(2 + group_num), sum, .names = "{.col}")) %>%
  bind_rows(mutation_mut_mat_summary[mutation_mut_mat_summary$Condition == "IHD", ], .)

##### permutation mut_mat summary
permutation_mut_mat_summary <- permutation_mut_mat_all %>% 
  mutate(mut_sum = rowSums(dplyr::select(., 1:96))) %>%
  dplyr::select(c("perm.id", "Case_ID", "mut_sum", "decile", "Condition")) %>% 
  group_by(Case_ID, decile) %>% 
  summarise(mut_sum = mean(mut_sum), .groups = "drop") %>% 
  group_by(Case_ID, decile) %>% 
  summarise(Value = mut_sum, .groups = "drop") %>% 
  pivot_wider(names_from = decile, values_from = Value) %>% 
  mutate(mut_num_percase = rowSums(dplyr::select(., 2:(2 + group_num - 1)))) %>% 
  arrange(match(Case_ID, Hypoxia_PTA_all_SCAN2$Case_ID)) %>% 
  mutate(Condition = Hypoxia_PTA_all_SCAN2_collapsed$Condition, 
         Case_ID = Hypoxia_PTA_all_SCAN2_collapsed$Case_ID, 
         Age = Hypoxia_PTA_all_SCAN2_collapsed$Age)

permutation_mut_mat_summary_ctrl <- permutation_mut_mat_summary %>% filter(Condition == "Control") %>% 
  summarise(across(2:(2 + group_num), sum, .names = "{.col}")) %>%
  bind_rows(permutation_mut_mat_summary[permutation_mut_mat_summary$Condition == "Control", ], .)

permutation_mut_mat_summary_dis <- permutation_mut_mat_summary %>% filter(Condition == "IHD") %>% 
  summarise(across(2:(2 + group_num), sum, .names = "{.col}")) %>%
  bind_rows(permutation_mut_mat_summary[permutation_mut_mat_summary$Condition == "IHD", ], .)

##### plot 
mutation_mut_mat_summary_plot <- mutation_mut_mat_all %>% mutate(mut_num = rowSums(dplyr::select(., 1:96))) %>%
  dplyr::select(c("Case_ID", "mut_num", "decile", "Condition"))
permutation_mut_mat_summary_plot <- permutation_mut_mat_all %>% mutate(permut_sum = rowSums(dplyr::select(., 1:96))) %>%
  dplyr::select(c("perm.id", "Case_ID", "permut_sum", "decile", "Condition"))

merged_summary_plot <- merge(mutation_mut_mat_summary_plot, permutation_mut_mat_summary_plot) %>% 
  merge(Hypoxia_PTA_all_SCAN2_collapsed[, c("Case_ID", "Age")]) %>% 
  mutate(enrichment_ratio = mut_num / permut_sum) %>% 
  filter(!is.na(enrichment_ratio) & is.finite(enrichment_ratio)) %>% 
  filter(enrichment_ratio != 0) %>%
  group_by(Condition, decile) %>% 
  mutate(enrichment_ratio = Winsorize(enrichment_ratio, probs = c(0.05, 0.95))) %>%
  # filter(enrichment_ratio >= quantile(enrichment_ratio, 0.05) & enrichment_ratio <= quantile(enrichment_ratio, 0.95)) %>%
  summarise(mean_ER = mean(enrichment_ratio, na.rm = TRUE, .groups = "drop"), sd_ER = sd(enrichment_ratio, na.rm = TRUE)) %>% 
  mutate(Condition = factor(Condition, level = c("Control", "IHD")), 
         decile = factor(decile, level = seq(1:group_num)))

overall_p <- wilcox.test(mean_ER ~ Condition, data = merged_summary_plot, alternative = c("two.sided"))$p.value
overall_star <- case_when(overall_p <= 0.001 ~ "***", overall_p <= 0.01  ~ "**", overall_p <= 0.05  ~ "*", TRUE ~ "ns")
overall_label <- paste0("Wilcoxon test Control v.s. IHD, P = ", signif(overall_p, 2))

p_SNV_RNAseq_enrichment <- ggplot(merged_summary_plot, aes(x = decile, y = mean_ER, group = Condition, color = Condition)) + 
  geom_hline(yintercept = 1, color = "black", linewidth = 0.6) + geom_line(position = position_dodge(width = 0.1), linewidth = 0.7) + 
  geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE, label.x.npc = "right", hjust = 1) +
  geom_errorbar(aes(ymin = mean_ER - sd_ER, ymax = mean_ER + sd_ER), width = 0.2, position = position_dodge(width = 0.1)) +
  geom_smooth(data = merged_summary_plot, aes(x = decile, y = mean_ER, color = Condition, fill = Condition, group = Condition),
              linetype = "22", method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.7) + 
  annotate(geom = "polygon", x = c(1, group_num, group_num), y = c(0.4, 0.4, 0.48), fill = "orange") + 
  scale_color_manual(values = c("Control" = color_set[1], "IHD" = color_set[2]), guide = "legend") + 
  scale_fill_manual(values = c("Control" = color_set[1], "IHD" = color_set[2]), guide = "legend") + theme_linedraw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24)) + 
  scale_y_continuous(breaks = seq(0.4, 1.8, by = 0.2), limits = c(0.4, 1.85)) +
  labs(x = "Gene expression level", y = "sSNV enrichment ratio \n (obs/exp)", color = "Condition", subtitle = overall_label)
ggsave(paste0(fig_save_dir, "/sSNV_snRNAseq_Enrichment.pdf"), plot = p_SNV_RNAseq_enrichment, width = 9, height = 5, dpi = 600)

################################################################################
##### plot top signature contribution from SigNet for mutation and permutation
SigNet_contri_mutation <- read.csv(paste0(fig_save_dir, "mutation_mut_mat_all_est/weight_guesses.csv"), row.names = 1)
# SigNet_contri_permutation <- read.csv(paste0(fig_save_dir, "permutation_mut_mat_all_est/weight_guesses.csv"), row.names = 1)
# SigNet_contri_mutation <- read.csv(paste0(fig_save_dir, "mutation_mut_mat_all_raw/weight_guesses.csv"), row.names = 1)
# SigNet_contri_permutation <- read.csv(paste0(fig_save_dir, "permutation_mut_mat_all_raw/weight_guesses.csv"), row.names = 1)

SigNet_contri_permutation <- c()
for (weight_i in 1:10) {
  cat(weight_i)
  permut_weight_i <- read.csv(paste0(fig_save_dir, "permutation_mut_mat_all_est_split_SigNet/permutation_mut_mat_all_est_part", weight_i, "/weight_guesses.csv"), row.names = 1)
  SigNet_contri_permutation <- rbind(SigNet_contri_permutation, permut_weight_i)
}

# SigNet_contri_permutation <- c()
# for (weight_i in 1:10) {
#   cat(weight_i)
#   permut_weight_i <- read.csv(paste0(fig_save_dir, "permutation_mut_mat_all_raw_split_SigNet/permutation_mut_mat_all_raw_part", weight_i, "/weight_guesses.csv"), row.names = 1)
#   SigNet_contri_permutation <- rbind(SigNet_contri_permutation, permut_weight_i)
# }

## normalize total contribution for each individual to estimated burden
SigNet_contri_mutation_est <- 1 / rowSums(SigNet_contri_mutation) * SigNet_contri_mutation * mutation_mut_mat_all_burden_decile$burden_decile
SigNet_contri_permutation_est <- 1 / rowSums(SigNet_contri_permutation) * SigNet_contri_permutation * permutation_mut_mat_all_burden_decile$burden_decile
SigNet_contri_mutation_est[zero_sum_rows_mut, ] <- 0
SigNet_contri_permutation_est[zero_sum_rows_permut, ] <- 0

# signature_list <- c("SBS5", "SBS8", "SBS30", "SBS32", "SBS40a", "SBS44")
signature_list <- c("SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS8", "SBS19", "SBS30", "SBS32", "SBS40a", "SBS44")
SigNet_contri_mutation_plot <- as.data.frame(SigNet_contri_mutation_est) %>%
  mutate(decile = sub(".*_", "", rownames(.)), Case_ID = sub("_.*", "", rownames(.)), 
         Condition = ifelse(Case_ID %in% Hypoxia_PTA_all_SCAN2_collapsed$Case_ID[Hypoxia_PTA_all_SCAN2_collapsed$Condition == "Control"], "Control", 
                            ifelse(Case_ID %in% Hypoxia_PTA_all_SCAN2_collapsed$Case_ID[Hypoxia_PTA_all_SCAN2_collapsed$Condition == "IHD"], "IHD", "Unknown"))) %>% 
  dplyr::select(c(signature_list, "decile", "Case_ID", "Condition")) %>% 
  setNames(c(paste0("mut_", signature_list), "decile", "Case_ID", "Condition"))

SigNet_contri_permutation_plot <- as.data.frame(SigNet_contri_permutation_est) %>% 
  mutate(decile = permutation_mut_mat_all_burden_decile$decile, perm.id = permutation_mut_mat_all_burden_decile$perm.id, Case_ID = sub("_.*", "", rownames(.)), 
         Condition = ifelse(Case_ID %in% Hypoxia_PTA_all_SCAN2_collapsed$Case_ID[Hypoxia_PTA_all_SCAN2_collapsed$Condition == "Control"], "Control", 
    ifelse(Case_ID %in% Hypoxia_PTA_all_SCAN2_collapsed$Case_ID[Hypoxia_PTA_all_SCAN2_collapsed$Condition == "IHD"], "IHD", "Unknown"))) %>% 
  dplyr::select(c(signature_list, "perm.id", "decile", "Case_ID", "Condition")) %>% 
  setNames(c(paste0("permut_", signature_list), "perm.id", "decile", "Case_ID", "Condition"))

merged_SigNet_contri_plot <- merge(SigNet_contri_mutation_plot, SigNet_contri_permutation_plot) %>% 
  merge(Hypoxia_PTA_all_SCAN2_collapsed[, c("Case_ID", "Age")])

## select signatures that are above certain percentage in either condition
SigNet_contri_mutation_threshold <- as.data.frame(SigNet_contri_mutation) %>%
  mutate(decile = sub(".*_", "", rownames(.)), Case_ID = sub("_.*", "", rownames(.)), 
         Condition = ifelse(Case_ID %in% Hypoxia_PTA_all_SCAN2_collapsed$Case_ID[Hypoxia_PTA_all_SCAN2_collapsed$Condition == "Control"], "Control", 
                            ifelse(Case_ID %in% Hypoxia_PTA_all_SCAN2_collapsed$Case_ID[Hypoxia_PTA_all_SCAN2_collapsed$Condition == "IHD"], "IHD", "Unknown"))) %>%
  dplyr::select(c(signature_list, "decile", "Case_ID", "Condition")) %>%
  setNames(c(signature_list, "decile", "Case_ID", "Condition"))

sig_pct_thr <- 5
sig_pct_sum_thr <- 10
all_sig_contri <- SigNet_contri_mutation_threshold %>% dplyr::select(all_of(signature_list), Condition) %>% 
  group_by(Condition) %>% summarise(across(all_of(signature_list), ~ sum(.x, na.rm = TRUE))) %>% 
  column_to_rownames("Condition") %>% t() %>% as.data.frame() %>% 
  mutate(Ctrl_pct = round(Control / sum(Control) * 100), Dis_pct = round(IHD / sum(IHD) * 100), 
         Pct_sum = Ctrl_pct + Dis_pct) %>% dplyr::select(Ctrl_pct, Dis_pct, Pct_sum)

selected_sigs <- all_sig_contri %>% filter(Ctrl_pct >= sig_pct_thr | Dis_pct >= sig_pct_thr) %>% filter(Pct_sum >= sig_pct_sum_thr)
# selected_sigs_list <- rownames(selected_sigs)
selected_sigs_list <- c("SBS1", "SBS5")

# for (sig_id_index in 1:length(selected_sigs_list)){
for (sig_id_index in 2){
  sig_id <- selected_sigs_list[sig_id_index]
  SigNet_contri_temp <- merged_SigNet_contri_plot %>% dplyr::select("perm.id", "decile", "Case_ID", "Condition", "Age", paste0("mut_", sig_id), paste0("permut_", sig_id)) %>% 
    setNames(c("perm.id", "decile", "Case_ID", "Condition", "Age", "mut_sig", "permut_sig")) %>% 
    mutate(EnR_sig = mut_sig / permut_sig) %>% 
    filter(!is.na(EnR_sig) & is.finite(EnR_sig)) %>% 
    filter(EnR_sig != 0) %>% 
    group_by(Condition, decile) %>% 
    mutate(EnR_sig = Winsorize(EnR_sig, probs = c(0.05, 0.95))) %>% 
    # filter(EnR_sig >= quantile(EnR_sig, 0.05) & EnR_sig <= quantile(EnR_sig, 0.95)) %>%
    summarise(mean_EnR = mean(EnR_sig, na.rm = TRUE), sd_EnR = sd(EnR_sig, na.rm = TRUE)) %>% 
    mutate(Condition = factor(Condition, level = c("Control", "IHD")), decile = factor(decile, level = seq(1:group_num)))
  
  overall_p <- wilcox.test(mean_EnR ~ Condition, data = SigNet_contri_temp, alternative = c("two.sided"))$p.value
  overall_star <- case_when(overall_p <= 0.001 ~ "***", overall_p <= 0.01  ~ "**", overall_p <= 0.05  ~ "*", TRUE ~ "ns")
  overall_label <- paste0("Wilcoxon test Control v.s. IHD, P = ", signif(overall_p, 2))

  p_enrichment_SigNet <- ggplot(SigNet_contri_temp, aes(x = decile, y = mean_EnR, color = Condition, group = Condition)) + 
    geom_hline(yintercept = 1, color = "black", linewidth = 0.6) + geom_line(position = position_dodge(width = 0.1), linewidth = 0.7) + 
    geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE, label.x.npc = "right", hjust = 1) + 
    geom_errorbar(aes(ymin = mean_EnR - sd_EnR, ymax = mean_EnR + sd_EnR), width = 0.2, position = position_dodge(width = 0.1)) +
    geom_smooth(data = SigNet_contri_temp, aes(x = decile, y = mean_EnR, color = Condition, fill = Condition, group = Condition), 
                linetype = "22", method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.7) + 
    scale_color_manual(values = c("Control" = color_set[1], "IHD" = color_set[2]), guide = "legend") + 
    scale_fill_manual(values = c("Control" = color_set[1], "IHD" = color_set[2]), guide = "legend") + theme_linedraw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24)) + 
    # annotate(geom = "polygon", x = c(1, group_num, group_num), y = c(-0.7, -0.7, -0.7 + 0.4), fill = "orange") +
    # scale_y_continuous(breaks = seq(0, 7, by = 1), limits = c(-0.7, 7)) +
    annotate(geom = "polygon", x = c(1, group_num, group_num), y = c(0.0, 0.0, 0.0 + 0.25), fill = "orange") +
    scale_y_continuous(breaks = seq(0.0, 4.0, by = 1.0), limits = c(0.0, 4.5)) +
    labs(x = "Gene expression level", y = paste0(sig_id, " enrichment ratio \n (obs/exp)"), color = "Condition", subtitle = overall_label)
  ggsave(paste0(fig_save_dir, sig_id, "_snRNAseq_Enrichment.pdf"), plot = p_enrichment_SigNet, width = 8.8, height = 5)
}
