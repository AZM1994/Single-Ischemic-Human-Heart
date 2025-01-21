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
group_num <- 8
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

analysis_ID <- paste0(permutation_round, "P_", group_num, "G_10sigs_O2")
fig_save_dir <- paste0("./results/O2/ATACseq_results/", analysis_ID, "/")
dir.create(fig_save_dir, recursive = TRUE)

### read in all mutmat O2
O2_results_dir <- paste0("./data/ATACseq_all_permutation_results/1000P_", group_num, "G/")
mutation_mut_mat_all <- c()
mutation_mut_mat_all <- readRDS(paste0(O2_results_dir, "/1000_perms_", group_num, "G_run_1/mutation_mut_mat_all_by_case_decile.rds"))
permutation_mut_mat_all <- c()
for (O2_batch in seq(1 : 10)){
  permutation_mut_mat_batch <- readRDS(paste0(O2_results_dir, "/1000_perms_", group_num, "G_run_", O2_batch, "/permutation_mut_mat_all_by_case_decile.rds"))
  permutation_mut_mat_all <- rbind(permutation_mut_mat_all, permutation_mut_mat_batch)
}

## normalize mut_mat to eat snv.burden (per individual)
mutation_avg_snv_burden_by_individual <- Hypoxia_PTA_Cases_metadata %>% 
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

snv.burden_normal <- as.matrix(mutation_mut_mat_all_burden_decile$burden_decile[mutation_mut_mat_all_burden_decile$condition == "Normal"])
snv.burden_disease <- as.matrix(mutation_mut_mat_all_burden_decile$burden_decile[mutation_mut_mat_all_burden_decile$condition == "Disease"])
permutation_avg_snv_burden_by_individual <- rbind(as.matrix(rep(snv.burden_normal, permutation_round)), 
                                                  as.matrix(rep(snv.burden_disease, permutation_round)))

mutation_mut_mat_all_modified_est <- 1 / rowSums(mutation_mut_mat_all[, 1:96]) * mutation_mut_mat_all[, 1:96] * 
  mutation_mut_mat_all_burden_decile$burden_decile
permutation_mut_mat_all_modified_est <- 1 / rowSums(permutation_mut_mat_all[, 1:96]) * permutation_mut_mat_all[, 1:96] * 
  permutation_avg_snv_burden_by_individual

### replace all 0s and NaNs in mutation and permutation results by 1
zero_sum_rows_mut <- rowSums(mutation_mut_mat_all_modified_est) %in% c(0, NaN)
mutation_mut_mat_all_modified_est[zero_sum_rows_mut, ] <- 1

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

##### permutation mut_mat summary
permutation_mut_mat_summary <- permutation_mut_mat_all %>% 
  mutate(mut_sum = rowSums(select(., 1:96))) |> base::`[`(c("Case_ID", "mut_sum", "decile", "condition", "perm.id_batch")) %>% 
  group_by(Case_ID, decile) %>% 
  summarise(mut_sum = mean(mut_sum), .groups = "drop") %>% 
  group_by(Case_ID, decile) %>% 
  summarise(Value = mut_sum, .groups = "drop") %>% 
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

##### plot 
mutation_mut_mat_summary_plot0 <- mutation_mut_mat_all %>% 
  mutate(mut_num = rowSums(select(., 1:96))) |> base::`[`(c("Case_ID", "mut_num", "decile", "condition"))

permutation_mut_mat_summary_plot0 <- permutation_mut_mat_all %>% 
  mutate(permut_sum = rowSums(select(., 1:96))) |> base::`[`(c("Case_ID", "permut_sum", "decile", "condition", "perm.id_batch")) %>% 
  group_by(Case_ID, decile) %>% 
  summarise(permut_sum = mean(permut_sum))

merged_summary_plot <- merge(mutation_mut_mat_summary_plot0, permutation_mut_mat_summary_plot0) %>% 
  merge(Hypoxia_PTA_Cases_metadata_collapsed[, c("Case_ID", "age")]) %>% 
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
  geom_point(position = position_dodge(width = 0.1), size = 2) + 
  stat_cor(size = 6, show.legend = FALSE, label.x.npc = "right", hjust = 1) + 
  geom_errorbar(aes(ymin = mean_ER - sd_ER, ymax = mean_ER + sd_ER), width = 0.2, position = position_dodge(width = 0.1)) +
  geom_smooth(data = merged_summary_plot, aes(x = decile, y = mean_ER, color = Condition, fill = Condition, group = Condition),
              method = "lm", se = TRUE, alpha = 0.2, linewidth = 1, linetype = "dashed") +
  theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position="right") +
  theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
  ylim(c(0.8, 1.3)) + 
  labs(x = "Chromosome accessibility level", y = "Somatic Mutation enrichment ratio \n (obs/exp)", color = "Condition", title = "")
ggsave(paste0(fig_save_dir, "/enrichment_ATACseq_snv_filtered.pdf"), width = 8, height = 5)