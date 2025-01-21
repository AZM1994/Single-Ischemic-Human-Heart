library(ggplot2)
library(stringr)
library(readxl)
library(ggsci)
library(dplyr)
library(tidyr)
library(tibble)
library(reshape2)
library(Seurat)
library(pheatmap)
library(readxl)
library(ggpubr)
ref_genome="BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = T)
library(MutationalPatterns)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

setwd("/Users/zhemingan/Documents/BCH_research/Gene_Expression_Analysis")
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])
group_num <- 9
batch_size <- 1
permutation_round <- 1000 / batch_size
perm_round_select <- 1 : permutation_round

##### read in metadata
Hypoxia_PTA_Cases_metadata <- readRDS("./data/SCAN2_df.rds") %>%
  as.data.frame() |> base::`[`(c("cell_ID", "age", "gender", "Case_ID", "condition", "snv.burden", "snv.rate.per.gb")) %>%
  rename_with(~ c("Cell_ID", "Condition"), .cols = c(1, 5))
Hypoxia_PTA_Cases_metadata_collapsed <- Hypoxia_PTA_Cases_metadata %>% distinct(Case_ID, .keep_all = TRUE)
Cell_ID_list <- Hypoxia_PTA_Cases_metadata$Cell_ID
Case_ID_list <- Hypoxia_PTA_Cases_metadata_collapsed$Case_ID
Condition_list <- unique(Hypoxia_PTA_Cases_metadata$Condition)
genomic_context_colnames <- c("Cell_ID", "Case_ID", "Condition", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene")

analysis_ID <- paste0(permutation_round, "P_", group_num, "G_3sigs_O2")
fig_save_dir <- paste0("./figures/ATAC_enrichment_analysis/", analysis_ID, "/")
dir.create(fig_save_dir)

### read in all mutmat O2
O2_results_dir <- paste0("/Users/zhemingan/Documents/BCH_research/Gene_Expression_Analysis/O2_ATAC/1000P_", group_num, "G/")
mutation_mut_mat_all <- c()
mutation_mut_mat_all <- readRDS(paste0(O2_results_dir, "/1000_perms_", group_num, "G_run_1/mutation_mut_mat_all_by_case_decile.rds"))
permutation_mut_mat_all <- c()
for (O2_batch in seq(1 : 10)){
  permutation_mut_mat_batch <- readRDS(paste0(O2_results_dir, "/1000_perms_", group_num, "G_run_", O2_batch, "/permutation_mut_mat_all_by_case_decile.rds"))
  permutation_mut_mat_all <- rbind(permutation_mut_mat_all, permutation_mut_mat_batch)
}

# saveRDS(mutation_mut_mat_all, paste0(fig_save_dir, "/mutation_mut_mat_all_by_case_decile.rds"))
# saveRDS(permutation_mut_mat_all, paste0(fig_save_dir, "/permutation_mut_mat_all_by_case_decile.rds"))
# write.csv(mutation_mut_mat_all[, 1:96], paste0(fig_save_dir, "mutation_mut_mat_all_by_case_decile.csv"))
# write.csv(permutation_mut_mat_all[, 1:96], paste0(fig_save_dir, "permutation_mut_mat_all_by_case_decile.csv"))

# ### replace all 0s in permutation results by 1
# mutation_mut_mat_all_modified <- mutation_mut_mat_all[, 1:96]
# zero_sum_rows_mut <- rowSums(mutation_mut_mat_all_modified) == 0
# mutation_mut_mat_all_modified[zero_sum_rows_mut, ] <- 1
# 
# permutation_mut_mat_all_modified <- permutation_mut_mat_all[, 1:96]
# zero_sum_rows_permut <- rowSums(permutation_mut_mat_all_modified) == 0
# permutation_mut_mat_all_modified[zero_sum_rows_permut, ] <- 1

# write.csv(mutation_mut_mat_all_modified, paste0(fig_save_dir, "mutation_mut_mat_all_modified.csv"))
# write.csv(permutation_mut_mat_all_modified, paste0(fig_save_dir, "permutation_mut_mat_all_modified.csv"))

## normalize mut_mat to eat snv.burden (per individual)
mutation_avg_snv_burden_by_individual <- Hypoxia_PTA_Cases_metadata %>% 
  # filter(Cell_ID != "1743_F3") %>% 
  group_by(Case_ID) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% 
  arrange(match(Case_ID, Hypoxia_PTA_Cases_metadata_collapsed$Case_ID)) %>% 
  mutate(Condition = Hypoxia_PTA_Cases_metadata_collapsed$Condition)

## calculate burden for each individual and decile
mutation_mut_mat_all_burden_decile <- mutation_mut_mat_all %>% 
  rowwise() %>%
  mutate(row_sum_decile = sum(c_across(1:96))) %>%
  ungroup() %>% 
  group_by(Case_ID) %>% 
  mutate(percentage_decile = row_sum_decile / sum(row_sum_decile) * 100) %>% 
  ungroup() %>% 
  mutate(burden_individual = rep(mutation_avg_snv_burden_by_individual$snv.burden, each = group_num)) %>% 
  mutate(burden_decile = percentage_decile * rep(mutation_avg_snv_burden_by_individual$snv.burden, each = group_num) / 100)

# snv.burden_normal <- as.matrix(mutation_avg_snv_burden_by_individual$snv.burden[mutation_avg_snv_burden_by_individual$Condition == "Normal"])
# snv.burden_disease <- as.matrix(mutation_avg_snv_burden_by_individual$snv.burden[mutation_avg_snv_burden_by_individual$Condition == "Disease"])
# permutation_avg_snv_burden_by_individual <- rbind(as.matrix(rep(rep(snv.burden_normal, each = group_num), permutation_round)), 
#                                                   as.matrix(rep(rep(snv.burden_disease, each = group_num), permutation_round)))
snv.burden_normal <- as.matrix(mutation_mut_mat_all_burden_decile$burden_decile[mutation_mut_mat_all_burden_decile$condition == "Normal"])
snv.burden_disease <- as.matrix(mutation_mut_mat_all_burden_decile$burden_decile[mutation_mut_mat_all_burden_decile$condition == "Disease"])
permutation_avg_snv_burden_by_individual <- rbind(as.matrix(rep(snv.burden_normal, permutation_round)), 
                                                  as.matrix(rep(snv.burden_disease, permutation_round)))

# mutation_mut_mat_all_modified_est <- 1 / rowSums(mutation_mut_mat_all_modified) * mutation_mut_mat_all_modified * 
#   rep(mutation_avg_snv_burden_by_individual$snv.burden, each = group_num)
mutation_mut_mat_all_modified_est <- 1 / rowSums(mutation_mut_mat_all[, 1:96]) * mutation_mut_mat_all[, 1:96] * 
  mutation_mut_mat_all_burden_decile$burden_decile
permutation_mut_mat_all_modified_est <- 1 / rowSums(permutation_mut_mat_all[, 1:96]) * permutation_mut_mat_all[, 1:96] * 
  permutation_avg_snv_burden_by_individual

### replace all 0s and NaNs in mutation and permutation results by 1
zero_sum_rows_mut <- rowSums(mutation_mut_mat_all_modified_est) %in% c(0, NaN)
mutation_mut_mat_all_modified_est[zero_sum_rows_mut, ] <- 1

# permutation_mut_mat_all_modified <- permutation_mut_mat_all[, 1:96]
zero_sum_rows_permut <- rowSums(permutation_mut_mat_all_modified_est) %in% c(0, NaN)
permutation_mut_mat_all_modified_est[zero_sum_rows_permut, ] <- 1

write.csv(mutation_mut_mat_all_modified_est, paste0(fig_save_dir, "mutation_mut_mat_all_modified_est.csv"))
write.csv(permutation_mut_mat_all_modified_est, paste0(fig_save_dir, "permutation_mut_mat_all_modified_est.csv"))

##### mutation mut_mat summary
mutation_mut_mat_summary <- mutation_mut_mat_all %>% 
  mutate(mut_sum = rowSums(select(., 1:96))) |> base::`[`(c("Case_ID", "mut_sum", "decile", "condition")) %>% 
  group_by(Case_ID, decile) %>% 
  summarise(Value = mut_sum, .groups = "drop") %>% 
  pivot_wider(names_from = decile, values_from = Value) %>% 
  mutate(mut_num_percase = rowSums(select(., 2:(2 + group_num - 1)))) %>% 
  arrange(match(Case_ID, Hypoxia_PTA_Cases_metadata$Case_ID)) %>% 
  mutate(Condition = Hypoxia_PTA_Cases_metadata_collapsed$Condition) %>% 
  mutate(Case_ID = Hypoxia_PTA_Cases_metadata_collapsed$Case_ID) %>% 
  mutate(Age = Hypoxia_PTA_Cases_metadata_collapsed$age)

mutation_mut_mat_summary_Normal <- mutation_mut_mat_summary[mutation_mut_mat_summary$Condition == "Normal", ] %>% 
  summarise(across(2:(2 + group_num), sum, .names = "{.col}")) %>%
  bind_rows(mutation_mut_mat_summary[mutation_mut_mat_summary$Condition == "Normal", ], .)

mutation_mut_mat_summary_IHD <- mutation_mut_mat_summary[mutation_mut_mat_summary$Condition == "Disease", ] %>% 
  summarise(across(2:(2 + group_num), sum, .names = "{.col}")) %>%
  bind_rows(mutation_mut_mat_summary[mutation_mut_mat_summary$Condition == "Disease", ], .)

# write.csv(mut_num_genic_merged, paste0(fig_save_dir, "mut_num_genic_merged.csv"))
# write.csv(mutation_mut_mat_summary_Normal, paste0(fig_save_dir, "mutation_mut_mat_summary_Normal.csv"))
# write.csv(mutation_mut_mat_summary_IHD, paste0(fig_save_dir, "mutation_mut_mat_summary_IHD.csv"))
# write.csv(mutation_mut_mat_summary, paste0(fig_save_dir, "mutation_mut_mat_summary.csv"))

##### permutation mut_mat summary
permutation_mut_mat_summary <- permutation_mut_mat_all %>% 
  mutate(mut_sum = rowSums(select(., 1:96))) |> base::`[`(c("Case_ID", "mut_sum", "decile", "condition", "perm.id_batch")) %>% 
  group_by(Case_ID, decile) %>% 
  summarise(mut_sum = mean(mut_sum), .groups = "drop") %>% 
  group_by(Case_ID, decile) %>% 
  summarise(Value = mut_sum, .groups = "drop") %>% 
  # reframe(Value = mut_sum) %>% 
  pivot_wider(names_from = decile, values_from = Value) %>% 
  mutate(mut_num_percase = rowSums(select(., 2:(2 + group_num - 1)))) %>% 
  arrange(match(Case_ID, Hypoxia_PTA_Cases_metadata$Case_ID)) %>% 
  mutate(Condition = Hypoxia_PTA_Cases_metadata_collapsed$Condition) %>% 
  mutate(Case_ID = Hypoxia_PTA_Cases_metadata_collapsed$Case_ID) %>% 
  mutate(Age = Hypoxia_PTA_Cases_metadata_collapsed$age)

permutation_mut_mat_summary_Normal <- permutation_mut_mat_summary[permutation_mut_mat_summary$Condition == "Normal", ] %>% 
  summarise(across(2:(2 + group_num), sum, .names = "{.col}")) %>%
  bind_rows(permutation_mut_mat_summary[permutation_mut_mat_summary$Condition == "Normal", ], .)

permutation_mut_mat_summary_IHD <- permutation_mut_mat_summary[permutation_mut_mat_summary$Condition == "Disease", ] %>% 
  summarise(across(2:(2 + group_num), sum, .names = "{.col}")) %>%
  bind_rows(permutation_mut_mat_summary[permutation_mut_mat_summary$Condition == "Disease", ], .)

# write.csv(permutation_mut_mat_summary_Normal, paste0(fig_save_dir, "permutation_mut_mat_summary_Normal.csv"))
# write.csv(permutation_mut_mat_summary_IHD, paste0(fig_save_dir, "permutation_mut_mat_summary_IHD.csv"))
# write.csv(permutation_mut_mat_summary, paste0(fig_save_dir, "permutation_mut_mat_summary.csv"))

##### plot 
mutation_mut_mat_summary_plot0 <- mutation_mut_mat_all %>% 
  mutate(mut_num = rowSums(select(., 1:96))) |> base::`[`(c("Case_ID", "mut_num", "decile", "condition"))

permutation_mut_mat_summary_plot0 <- permutation_mut_mat_all %>% 
  mutate(permut_sum = rowSums(select(., 1:96))) |> base::`[`(c("Case_ID", "permut_sum", "decile", "condition", "perm.id_batch")) %>% 
  group_by(Case_ID, decile) %>% 
  summarise(permut_sum = mean(permut_sum))

merged_summary_plot <- merge(mutation_mut_mat_summary_plot0, permutation_mut_mat_summary_plot0) %>% 
  merge(Hypoxia_PTA_Cases_metadata_collapsed[, c("Case_ID", "age")]) %>% 
  # filter(age >= 0.5 & age <= 90) %>%
  mutate(enrichment_ratio = mut_num / permut_sum) %>% 
  filter(!is.na(enrichment_ratio) & is.finite(enrichment_ratio)) %>% 
  # filter(enrichment_ratio != 0) %>%
  # filter(!is.na(enrichment_ratio)) %>% 
  # mutate(enrichment_ratio = ifelse(is.infinite(enrichment_ratio), NA, enrichment_ratio)) %>% 
  # mutate(enrichment_ratio = Winsorize(enrichment_ratio, probs = c(0.05, 0.95))) %>%
  filter(enrichment_ratio >= quantile(enrichment_ratio, 0.25) - 1.5 * IQR(enrichment_ratio) &
           enrichment_ratio <= quantile(enrichment_ratio, 0.75) + 1.5 * IQR(enrichment_ratio)) %>%
  # filter(enrichment_ratio >= quantile(enrichment_ratio, 0.05) & enrichment_ratio <= quantile(enrichment_ratio, 0.95)) %>%
  group_by(condition, decile) %>% 
  summarise(mean_ER = mean(enrichment_ratio, na.rm = TRUE), 
            sd_ER = sd(enrichment_ratio, na.rm = TRUE)) %>% 
  mutate(Condition = ifelse(condition == "Normal", "Control", ifelse(condition == "Disease", "IHD", "Unknown"))) %>% 
  mutate(Condition = factor(Condition, level = c("Control", "IHD"))) %>% 
  mutate(decile = factor(decile, level = seq(1:group_num)))

ggplot(merged_summary_plot, aes(x = decile, y = mean_ER, group = Condition, color = Condition)) + 
  geom_hline(yintercept = 1, color = "black", linewidth = 0.6) + 
  geom_line(position = position_dodge(width = 0.1), size = 1) + 
  geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE) + 
  geom_errorbar(aes(ymin = mean_ER - sd_ER, ymax = mean_ER + sd_ER), width = 0.2, position = position_dodge(width = 0.1)) +
  geom_smooth(data = merged_summary_plot, aes(x = decile, y = mean_ER, color = Condition, fill = Condition, group = Condition),
              method = "lm", se = TRUE, alpha = 0.2, linewidth = 1, linetype = "dashed") +
  theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right") +
  theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
  ylim(c(0.5, 1.6)) + 
  labs(x = "Gene expression levels", y = "Somatic Mutation enrichment ratio \n (obs/exp)", color = "Condition", title = "")
ggsave(paste0(fig_save_dir, "/total_snv_filtered.pdf"), width = 8, height = 5)
# ggsave(paste0(fig_save_dir, "/total_snv.pdf"), width = 8, height = 5)



################################################################################
################################################################################
##### plot top signature contribution from SigNet for mutation and permutation
################################################################################
################################################################################
SigNet_contri_mutation <- read.csv(paste0(fig_save_dir, "mutation_mut_mat_all_by_case_decile_est/weight_guesses.csv"), row.names = 1)
SigNet_contri_permutation <- read.csv(paste0(fig_save_dir, "permutation_mut_mat_all_by_case_decile_est/weight_guesses.csv"), row.names = 1)
# SigNet_contri_mutation <- read.csv("/Users/zhemingan/Documents/BCH_research/SigNet/mutation_mut_mat_all_by_case_decile_est/weight_guesses.csv", row.names = 1)
# SigNet_contri_permutation <- read.csv("/Users/zhemingan/Documents/BCH_research/SigNet/permutation_mut_mat_all_by_case_decile_est/weight_guesses.csv", row.names = 1)
# SigNet_contri_mutation_score <- read.csv("/Users/zhemingan/Documents/BCH_research/SigNet/mutation_mut_mat_all_by_case_decile_est/classification_guesses.csv", row.names = 1)
# SigNet_contri_permutation_score <- read.csv("/Users/zhemingan/Documents/BCH_research/SigNet/permutation_mut_mat_all_by_case_decile_est/classification_guesses.csv", row.names = 1)

## normalize total contribution for each individual to estimated burden
SigNet_contri_mutation_est <- 1 / rowSums(SigNet_contri_mutation) * SigNet_contri_mutation *
  mutation_mut_mat_all_burden_decile$burden_decile
SigNet_contri_permutation_est <- 1 / rowSums(SigNet_contri_permutation) * SigNet_contri_permutation * permutation_avg_snv_burden_by_individual
# SigNet_contri_mutation_est <- SigNet_contri_mutation
# SigNet_contri_permutation_est <- SigNet_contri_permutation

# write.csv(SigNet_contri_mutation_est, paste0(fig_save_dir, "weight_guesses_SigNet_mutation_est.csv"))
# write.csv(SigNet_contri_permutation_est, paste0(fig_save_dir, "weight_guesses_SigNet_permutation_est.csv"))

# low_mut_num_index <- rowSums(mutation_mut_mat_all_modified_est) <= quantile(rowSums(mutation_mut_mat_all_modified_est), 0.10)
low_mut_num_index <- rowSums(mutation_mut_mat_all_modified_est) <= 15
sum(low_mut_num_index)
SigNet_contri_mutation_est[zero_sum_rows_mut, ] <- NA
# SigNet_contri_mutation_est[SigNet_contri_mutation_score < 0.1, ] <- NA
SigNet_contri_mutation_est[low_mut_num_index, ] <- NA
SigNet_contri_permutation_est[zero_sum_rows_permut, ] <- NA
# SigNet_contri_permutation_est[SigNet_contri_permutation_score < 0.1, ] <- NA
# signature_list <- c("SBS5", "SBS30", "SBS19", "SBS4", "SBS1")
# signature_list <- c("SBS5", "SBS30", "SBS19", "SBS4", "SBS1", "SBS2", "SBS92", "SBS32", "SBS8", "SBS40")
signature_list <- c("SBS5", "SBS30", "SBS19", "SBS4", "SBS1", "SBS2", "SBS92", "SBS32", "SBS8", "SBS40",
                    "SBS44", "SBS29", "SBS7b", "SBS89", "SBS36")
# signature_list <- c("SBS5", "SBS30", "SBS19", "SBS4", "SBS1", "SBS2", "SBS92", "SBS32", "SBS8", "SBS40",
#                     "SBS44", "SBS29", "SBS7b", "SBS89", "SBS36", "SBS3", "SBS22", "SBS16", "SBS39", "SBS18")
# signature_list <- colnames(SigNet_contri_mutation_est)
# signature_list <- c("SBS5", "SBS30", "SBS19", "SBS4", "SBS1", "SBS2", "SBS92", "SBS32", "SBS8", "SBS40",
#                     "SBS44", "SBS29", "SBS7b", "SBS89", "SBS36", "SBS3")

SigNet_contri_mutation_plot <- as.data.frame(SigNet_contri_mutation_est) %>% 
  mutate(decile = sub(".*_", "", rownames(.))) %>% 
  mutate(Case_ID = sub("_.*", "", rownames(.))) %>% 
  mutate(Condition = ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
    Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Normal", "Case_ID"], "Control", 
    ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
      Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Disease", "Case_ID"], "IHD", "Unknown"))) |> 
  base::`[`(c(signature_list, "decile", "Case_ID", "Condition")) %>% 
  setNames(c(paste0("mut_", signature_list), "decile", "Case_ID", "Condition"))

## select signatures that are above certain percentage in either condition
SigNet_contri_mutation_threshold <- as.data.frame(SigNet_contri_mutation) %>%
  mutate(decile = sub(".*_", "", rownames(.))) %>%
  mutate(Case_ID = sub("_.*", "", rownames(.))) %>%
  mutate(Condition = ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
    Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Normal", "Case_ID"], "Control",
    ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
      Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Disease", "Case_ID"], "IHD", "Unknown"))) |>
  base::`[`(c(signature_list, "decile", "Case_ID", "Condition")) %>%
  setNames(c(paste0("mut_", signature_list), "decile", "Case_ID", "Condition"))
# SigNet_contri_mutation_plot <- SigNet_contri_mutation_threshold
sig_percentage_threshold <- 8
# sig_percentage_sum_threshold <- 10

control_sig_contri <- as.data.frame(colSums(SigNet_contri_mutation_threshold
                                            [SigNet_contri_mutation_threshold$Condition == "Control", 1 : length(signature_list)], na.rm = TRUE)) %>% 
  setNames(c("Control_sig_contri")) %>% mutate(Control_Percentage = round(Control_sig_contri / sum(Control_sig_contri) * 100))
disease_sig_contri <- as.data.frame(colSums(SigNet_contri_mutation_threshold
                                            [SigNet_contri_mutation_threshold$Condition == "IHD", 1 : length(signature_list)], na.rm = TRUE)) %>% 
  setNames(c("IHD_sig_contri")) %>% mutate(IHD_Percentage = round(IHD_sig_contri / sum(IHD_sig_contri) * 100))

all_sig_contri <- cbind(control_sig_contri, disease_sig_contri) %>% 
  mutate(Percentage_Sum = Control_Percentage + IHD_Percentage) |> base::`[`(c("Control_Percentage", "IHD_Percentage", "Percentage_Sum"))

# selected_sigs <- all_sig_contri[all_sig_contri$Control_Percentage >= sig_percentage_threshold | 
#                                   all_sig_contri$IHD_Percentage >= sig_percentage_threshold & 
#                                   all_sig_contri$Percentage_Sum >= sig_percentage_sum_threshold, ]
selected_sigs <- all_sig_contri[all_sig_contri$Control_Percentage >= sig_percentage_threshold |
                                  all_sig_contri$IHD_Percentage >= sig_percentage_threshold, ]
# cat(rownames(selected_sigs))
## determine the signatures to show
# haha <- SigNet_contri_mutation_plot[SigNet_contri_mutation_plot$Condition == "Control", ] %>% 
#   filter(all_sig_contri$Control_Percentage < sig_percentage_threshold)

df_A <- SigNet_contri_mutation_plot %>% filter(Condition == 'Control')
df_B <- SigNet_contri_mutation_plot %>% filter(Condition == 'IHD')
df_A_repeated <- df_A[rep(seq_len(nrow(df_A)), permutation_round), ]
df_B_repeated <- df_B[rep(seq_len(nrow(df_B)), permutation_round), ]
df_combined <- rbind(df_A_repeated, df_B_repeated)

SigNet_contri_permutation_plot <- as.data.frame(SigNet_contri_permutation_est) %>% 
  mutate(decile = permutation_mut_mat_all$decile) %>% 
  mutate(Case_ID = sub("_.*", "", rownames(.))) |> 
  mutate(Condition = ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
    Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Normal", "Case_ID"], "Control", 
    ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
      Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Disease", "Case_ID"], "IHD", "Unknown"))) |> 
  base::`[`(c(signature_list, "decile", "Case_ID", "Condition")) %>% 
  setNames(c(paste0("permut_", signature_list), "decile", "Case_ID", "Condition"))

# sig_id = "SBS1"
# b <- SigNet_contri_all[SigNet_contri_all$Case_ID == "1363", ]
# a <- SigNet_contri_all[SigNet_contri_all$decile == 2 & SigNet_contri_all$condition == "Control", ]
# a <- SigNet_contri_all[SigNet_contri_all$decile == 1 & SigNet_contri_all$Case_ID == "1039",]
# a <- SigNet_contri_all[SigNet_contri_all$condition == "Control" & SigNet_contri_all$decile == 1,]
# a <- SigNet_contri_all[SigNet_contri_all$condition == "IHD",]
# a_1 <- a %>% group_by(Case_ID) %>% summarise(mean_EnR = mean(EnR_sig))
# for (sig_id in signature_list){
selected_sigs_control <- selected_sigs$Control_Percentage >= sig_percentage_threshold
selected_sigs_IHD <- selected_sigs$IHD_Percentage >= sig_percentage_threshold
selected_sigs_list <- gsub("mut_", "", rownames(selected_sigs))
for (sig_id_index in 1:length(selected_sigs_list)){
# for (sig_id_index in 1){
  sig_id <- selected_sigs_list[sig_id_index]
  SigNet_contri_all <- cbind(df_combined[c(paste0("mut_", sig_id), "decile", "Case_ID", "Condition")], 
                             SigNet_contri_permutation_plot[c(paste0("permut_", sig_id))]) %>% 
    setNames(c("mut_sig", "decile", "Case_ID", "Condition", "permut_sig")) %>% 
    merge(Hypoxia_PTA_Cases_metadata_collapsed[, c("Case_ID", "age")]) %>% 
    # filter(!is.na(mut_sig) & is.finite(mut_sig)) %>% 
    # filter(mut_sig >= quantile(mut_sig, 0.25) - 1.5 * IQR(mut_sig) & mut_sig <= quantile(mut_sig, 0.75) + 1.5 * IQR(mut_sig)) %>%
    # filter(!is.na(permut_sig) & is.finite(permut_sig)) %>% 
    # filter(permut_sig >= quantile(permut_sig, 0.25) - 1.5 * IQR(permut_sig) & permut_sig <= quantile(permut_sig, 0.75) + 1.5 * IQR(permut_sig)) %>%
    # filter(age != 44) %>%
    # filter(age > 20 & age < 70) %>%
    mutate(EnR_sig = mut_sig / permut_sig) %>% 
    # filter(EnR_sig != 0) %>%
    # filter(!is.na(EnR_sig)) %>% 
    # mutate(EnR_sig = ifelse(is.infinite(EnR_sig), 1, EnR_sig)) %>% 
    filter(!is.na(EnR_sig) & is.finite(EnR_sig)) %>% 
    # mutate(EnR_sig = DescTools::Winsorize(EnR_sig, probs = c(0.05, 0.95))) %>%
    # filter(EnR_sig >= quantile(EnR_sig, 0.25) - 1.5 * IQR(EnR_sig) & EnR_sig <= quantile(EnR_sig, 0.75) + 1.5 * IQR(EnR_sig)) %>%
    filter(EnR_sig >= quantile(EnR_sig, 0.05) & EnR_sig <= quantile(EnR_sig, 0.95)) %>%
    # group_by(Condition, decile) %>%
    # filter(EnR_sig >= quantile(EnR_sig, 0.25) - 1.5 * IQR(EnR_sig) & EnR_sig <= quantile(EnR_sig, 0.75) + 1.5 * IQR(EnR_sig)) %>%
    # filter(EnR_sig >= quantile(EnR_sig, 0.05) & EnR_sig <= quantile(EnR_sig, 0.95)) %>%
    # ungroup() %>%
    group_by(Case_ID, decile) %>% 
    # filter(!is.na(EnR_sig) & is.finite(EnR_sig)) %>% 
    # filter(EnR_sig >= quantile(EnR_sig, 0.25) - 1.5 * IQR(EnR_sig) & EnR_sig <= quantile(EnR_sig, 0.75) + 1.5 * IQR(EnR_sig)) %>%
    summarise(mean_EnR = mean(EnR_sig), sd_EnR = sd(EnR_sig)) %>% 
    mutate(condition = ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
      Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Normal", "Case_ID"], "Control", 
      ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
        Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Disease", "Case_ID"], "IHD", "Unknown"))) %>% 
    mutate(decile = factor(decile, level = seq(1:group_num)))
  
  SigNet_contri_all_final <- SigNet_contri_all %>% 
    group_by(decile, condition) %>% 
    # filter(!is.na(mean_EnR) & is.finite(mean_EnR)) %>% 
    # filter(mean_EnR >= quantile(mean_EnR, 0.25) - 1.5 * IQR(mean_EnR) & mean_EnR <= quantile(mean_EnR, 0.75) + 1.5 * IQR(mean_EnR)) %>%
    # filter(mean_EnR >= quantile(mean_EnR, 0.05) & mean_EnR <= quantile(mean_EnR, 0.95)) %>%
    summarise(mean_EnR_final = mean(mean_EnR), sd_EnR_final = sd(mean_EnR)) %>% 
    mutate(condition = factor(condition, level = c("Control", "IHD"))) %>% 
  filter(condition == c("Control", "IHD")[c(selected_sigs_control[sig_id_index], selected_sigs_IHD[sig_id_index])])
  
  ggplot(SigNet_contri_all_final, aes(x = decile, y = mean_EnR_final, color = condition, group = condition)) + 
    geom_hline(yintercept = 1, color = "black", linewidth = 0.6) + 
    geom_line(position = position_dodge(width = 0.1), size = 1) + 
    geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE) + 
    geom_errorbar(aes(ymin = mean_EnR_final - sd_EnR_final, ymax = mean_EnR_final + sd_EnR_final), width = 0.2, position = position_dodge(width = 0.1)) +
    geom_smooth(data = SigNet_contri_all_final, aes(x = decile, y = mean_EnR_final, color = condition, fill = condition, group = condition),
                method = "lm", se = TRUE, alpha = 0.2, linewidth = 1, linetype = "dashed") +
    theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position="right") + 
    theme_classic() + scale_color_manual(values = color_set[c(selected_sigs_control[sig_id_index], selected_sigs_IHD[sig_id_index])]) + 
    scale_fill_manual(values = color_set[c(selected_sigs_control[sig_id_index], selected_sigs_IHD[sig_id_index])]) + 
    # ylim(c(0.0, 3.2)) +
    labs(x = "Gene expression levels", y = paste0(sig_id, " enrichment ratio \n (obs/exp)"), color = "condition")
  ggsave(paste0(fig_save_dir, sig_id, "_enrichment_filtered.pdf"), width = 8, height = 5)
}
selected_sigs