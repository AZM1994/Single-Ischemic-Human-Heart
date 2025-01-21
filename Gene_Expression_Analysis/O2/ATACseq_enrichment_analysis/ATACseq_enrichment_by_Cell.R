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
batch_size <- 10
permutation_round <- 100 / batch_size
perm_round_select <- 1 : permutation_round
# perm_round_select <- 1:3

##### read in metadata
Hypoxia_PTA_Cases_metadata <- readRDS("./data/SCAN2_df.rds") %>%
  as.data.frame() |> base::`[`(c("cell_ID", "age", "gender", "Case_ID", "condition", "snv.burden", "snv.rate.per.gb")) %>%
  rename_with(~ c("Cell_ID", "Condition"), .cols = c(1, 5))
Hypoxia_PTA_Cases_metadata_collapsed <- Hypoxia_PTA_Cases_metadata %>% distinct(Case_ID, .keep_all = TRUE)
Cell_ID_list <- Hypoxia_PTA_Cases_metadata$Cell_ID
Case_ID_list <- Hypoxia_PTA_Cases_metadata_collapsed$Case_ID
Condition_list <- unique(Hypoxia_PTA_Cases_metadata$Condition)
genomic_context_colnames <- c("Cell_ID", "Case_ID", "Condition", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene")

##### read in scRNA-ATACseq data
atac_normal <- c()
atac_disease <- c()
atac_4638_2n <- read.table("./data/scATACseq/ATAC_cov_average_hg19_1kb_bin.4638_2n.bed") %>% setNames(c("chr", "start", "end", "id", "score_4638_2n"))
atac_5828_2n <- read.table("./data/scATACseq/ATAC_cov_average_hg19_1kb_bin.5828_2n.bed") %>% setNames(c("chr", "start", "end", "id", "score_5828_2n")) |> base::`[`(c("score_5828_2n"))
atac_5919_all <- read.table("./data/scATACseq/ATAC_cov_average_hg19_1kb_bin.5919_all.bed") %>% setNames(c("chr", "start", "end", "id", "score_5919_all")) |> base::`[`(c("score_5919_all"))
atac_normal <- cbind(atac_4638_2n, atac_5828_2n, atac_5919_all) %>% 
  mutate(score = rowMeans(select(., score_4638_2n, score_5828_2n, score_5919_all), na.rm = TRUE)) |> base::`[`(c("chr", "start", "end", "id", "score")) %>% 
  mutate(group = ntile(score, n = group_num)) %>% 
  mutate(group = as.factor(group)) %>% 
  mutate(Condition = "Normal")
atac_grange_normal <- GRanges(seqnames = atac_normal$chr, 
                       ranges = IRanges(start = atac_normal$start, 
                                        end = atac_normal$end),
                       group = atac_normal$group) 

atac_604_all <- read.table("./data/scATACseq/ATAC_cov_average_hg19_1kb_bin.604_all.bed") %>% setNames(c("chr", "start", "end", "id", "score_604_all"))
atac_disease <- atac_604_all %>% 
  setNames(c("chr", "start", "end", "id", "score")) %>% 
  mutate(group = ntile(score, n = group_num)) %>% 
  mutate(group = as.factor(group)) %>% 
  mutate(Condition = "Disease")
atac_grange_disease <- GRanges(seqnames = atac_disease$chr, 
                              ranges = IRanges(start = atac_disease$start, 
                                               end = atac_disease$end),
                              group = atac_disease$group) 

# atac <- read.table("./data/scATACseq/ATAC_cov_average_hg19_1kb_bin.bed")
# colnames(atac) <- c("chr", "start", "end", "id", "score")
# atac$group <- ntile(atac$score, n = group_num)
# atac_grange <- GRanges(seqnames = atac$chr, 
#                        ranges = IRanges(start = atac$start, 
#                                         end = atac$end),
#                        group = atac$group) 

##### read in permutation results for each cell
permutation_normal <- c()
permutation_disease <- c()
for (Cell_ID_temp in Cell_ID_list){
  permutation_path <- paste0("./data/permutation_by_cell/", Cell_ID_temp, "/permutation_snv.avinput.variant_function")
  case_temp0 <- Hypoxia_PTA_Cases_metadata$Case_ID[Hypoxia_PTA_Cases_metadata$Cell_ID == Cell_ID_temp]
  condition_temp0 <- Hypoxia_PTA_Cases_metadata$Condition[Hypoxia_PTA_Cases_metadata$Cell_ID == Cell_ID_temp]
  cat("Cell:", Cell_ID_temp, "Case:", case_temp0)
  if (condition_temp0 == "Normal"){
    cat(" Condition:", condition_temp0, "\n")
    permutation_temp <- read.table(permutation_path, header = FALSE) |> base::`[`(c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V15")) %>% 
      setNames(c("region", "gene", "chr", "start", "end", "ref", "alt", "perm.id")) %>% 
      mutate(perm.id_batch = ceiling(perm.id / batch_size)) %>% 
      mutate(Cell_ID = Cell_ID_temp) %>% 
      mutate(Case_ID = case_temp0) %>% 
      mutate(Condition = condition_temp0)
    permutation_normal <- rbind(permutation_normal, permutation_temp)
    }else if (condition_temp0 == "Disease"){
      cat(" Condition:", condition_temp0, "\n")
      permutation_temp <- read.table(permutation_path, header = FALSE) |> base::`[`(c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V15")) %>%
        setNames(c("region", "gene", "chr", "start", "end", "ref", "alt", "perm.id")) %>%
        mutate(perm.id_batch = ceiling(perm.id / batch_size)) %>%
        mutate(Cell_ID = Cell_ID_temp) %>% 
        mutate(Case_ID = case_temp0) %>% 
        mutate(Condition = condition_temp0)
      permutation_disease <- rbind(permutation_disease, permutation_temp)
    }
}

permutation_mut_mat_all <- c()
mutation_mut_mat_all <- c()
mut_num_genic_merged <- c()
mutation_num <- c()
permutation_num <- c()
# condition_temp = Condition_list[1]
for (condition_temp in Condition_list){
  ###########################################################################
  ##################### raw SCAN2 call mutation analysis ####################
  ###########################################################################
  cat("##### Raw SCAN2 call mutation analysis:", condition_temp, "...\n")
  cat("Get genomic context for", condition_temp, "...\n")
  Case_ID_order <- Hypoxia_PTA_Cases_metadata_collapsed[Hypoxia_PTA_Cases_metadata_collapsed$Condition == condition_temp, "Case_ID"]
  heart_PTA_Cases_vcf_temp <- read.table(paste0("data/heart_PTA_Cases.all_age.", condition_temp, "_ssnv.vcf"), sep = "\t")
  genomic_context_temp <- read.csv(paste0("data/heart_PTA_Cases.all_age.", condition_temp, ".annotation.csv"), header = TRUE) %>%
    mutate(Cell_ID = heart_PTA_Cases_vcf_temp$V8) %>% mutate(Case_ID = str_extract(Cell_ID, "[^_]+")) %>% 
    mutate(Condition = condition_temp) |> base::`[`(genomic_context_colnames)
  
  genomic_context_temp <- genomic_context_temp %>%
    mutate(Cell_ID = ifelse(Cell_ID == "1864_M-E3-2n", "1864_E3", ifelse(Cell_ID == "4402_1_A1-2n", "4402_A1", ifelse(Cell_ID == "4402_1_A3-2n", "4402_A3",
    ifelse(Cell_ID == "4638_1_A1-2n", "4638_A1", ifelse(Cell_ID == "4638_1_A4-2n", "4638_A4", ifelse(Cell_ID == "4638_1_B2-2n", "4638_B2",
    ifelse(Cell_ID == "5657_M-D1-2n", "5657_D1", ifelse(Cell_ID == "5657_M-H1-2n", "5657_H1", ifelse(Cell_ID == "5828_CM_C2_2n", "5828_C2",
    ifelse(Cell_ID == "5828_CM_G2_2n", "5828_G2", ifelse(Cell_ID == "5919_1_C4-2n", "5919_C4", ifelse(Cell_ID == "5919_1_D2-2n", "5919_D2",
    ifelse(Cell_ID == "5919_1_E3-2n", "5919_E3", ifelse(Cell_ID == "5919_1_F6-2n", "5919_F6", ifelse(Cell_ID == "6032_1_A1-2n", "6032_A1",
    ifelse(Cell_ID == "6032_M-C7-2n", "6032_C7", ifelse(Cell_ID == "6032_M-E7-2n", "6032_E7", ifelse(Cell_ID == "1673_CM_A2_2n", "1673_A2",
    ifelse(Cell_ID == "1673_CM_A3_2n", "1673_A3", ifelse(Cell_ID == "1673_CM_D2_2n", "1673_D2", Cell_ID)))))))))))))))))))))
  
  genomic_SCAN2_df_temp <- merge(genomic_context_temp, Hypoxia_PTA_Cases_metadata[c("Cell_ID", "Case_ID", "age", "gender", "Condition")]) %>% 
    rename_with(~ c("Cell_ID", "Case_ID", "Condition", "chr", "start", "end", "ref", "alt", "region", "gene"), .cols = 1:10) %>% 
    mutate(Condition = as.factor(Condition))

  if (condition_temp == "Normal"){
    atac_grange_temp <- atac_grange_normal
  }else if (condition_temp == "Disease"){
    atac_grange_temp <- atac_grange_disease
  }
  ##### raw SCAN2 call signature analysis
  cat("Raw SCAN2 call mutation signature analysis for", condition_temp, "...\n")
  mutation_decile_grange_list_temp <- GRangesList()
  mutation_num_temp <- c()
  for (case_index in Case_ID_order){
    genomic_SCAN2_df_by_case_temp <- genomic_SCAN2_df_temp[genomic_SCAN2_df_temp$Case_ID == case_index, ]
    mutation_grange_by_case_temp <- GRanges(seqnames = paste0("chr", genomic_SCAN2_df_by_case_temp$chr),
                               ranges = IRanges(start = genomic_SCAN2_df_by_case_temp$start,
                                                end = genomic_SCAN2_df_by_case_temp$end),
                               ref = genomic_SCAN2_df_by_case_temp$ref,
                               alt = genomic_SCAN2_df_by_case_temp$alt)
    mutation_num_by_case_temp <- c()
    for (i in 1:group_num){
      mutation_intersection <- findOverlaps(mutation_grange_by_case_temp, atac_grange_temp[atac_grange_temp$group == i,])
      mutation_intersected_ranges <- mutation_grange_by_case_temp[queryHits(mutation_intersection)]
      mutation_num_by_case_temp[i] <- length(mutation_intersected_ranges)
      mutation_decile_grange_list_temp[[paste0(case_index, "_", i)]] <- mutation_intersected_ranges
    }
    mutation_num_temp <- rbind(mutation_num_temp, mutation_num_by_case_temp)
  }
  
  mutation_num <- rbind(mutation_num, mutation_num_temp)
  
  chr_length <- seqlengths(Hsapiens)
  seqlengths(mutation_decile_grange_list_temp) <- chr_length[names(seqlengths(mutation_decile_grange_list_temp))]
  seqlevels(mutation_decile_grange_list_temp) <- seqlevels(mutation_decile_grange_list_temp)[order(factor(seqlevels(mutation_decile_grange_list_temp), levels = chr_orders))]
  genome(mutation_decile_grange_list_temp) = 'hg19'
  
  ##### calculate 96 types of snvs
  mutation_mut_mat_temp <- mut_matrix(mutation_decile_grange_list_temp, ref_genome = ref_genome)
  mutation_mut_mat_temp_save <- as.data.frame(t(mutation_mut_mat_temp)) %>% 
    mutate(condition = condition_temp) %>% 
    mutate(decile = sub(".*_", "", rownames(.))) %>% 
    mutate(Case_ID = sub("_.*", "", rownames(.)))

  mutation_mut_mat_all <- rbind(mutation_mut_mat_all, mutation_mut_mat_temp_save)
  
  ###########################################################################
  ###################### Permutation mutation analysis ######################
  ###########################################################################
  cat("##### Permutation mutation analysis:", condition_temp, "...\n")
  cat("Get permutation data for", condition_temp, "...\n")
  if (condition_temp == "Normal"){
    permutation_temp <- permutation_normal %>% filter(perm.id_batch %in% perm_round_select)
  }else if (condition_temp == "Disease"){
    permutation_temp <- permutation_disease %>% filter(perm.id_batch %in% perm_round_select)
  }
  
  ##### Permutation signature analysis
  cat("Permutation signature analysis for", condition_temp, "...\n")
  permutation_mut_mat_temp <- c()
  permutation_num_temp <- c()
  for (perm_index in perm_round_select){
    permutation_decile_grange_list_temp_perm_i <- GRangesList()
    permutation_num_temp_perm_i <- c()
    for (case_index in Case_ID_order){
      permutation_by_case_temp_perm_i <- permutation_temp[permutation_temp$perm.id_batch == perm_index & 
                                                            permutation_temp$Case_ID == case_index, ]
      permutation_grange_by_case_temp_perm_i <- GRanges(seqnames = permutation_by_case_temp_perm_i$chr,
                                      ranges = IRanges(start = permutation_by_case_temp_perm_i$start,
                                                       end = permutation_by_case_temp_perm_i$end),
                                      ref = permutation_by_case_temp_perm_i$ref,
                                      alt = permutation_by_case_temp_perm_i$alt)
      permutation_num_by_case_temp_perm_i <- c()
      for (i in 1:group_num){
        permutation_intersection <- findOverlaps(permutation_grange_by_case_temp_perm_i, atac_grange_temp[atac_grange_temp$group == i,])
        permutation_intersected_ranges <- permutation_grange_by_case_temp_perm_i[queryHits(permutation_intersection)]
        permutation_num_by_case_temp_perm_i[i] <- length(permutation_intersected_ranges)
        permutation_decile_grange_list_temp_perm_i[[paste0(case_index, "_", i)]] <- permutation_intersected_ranges
      }
      permutation_num_temp_perm_i <- rbind(permutation_num_temp_perm_i, permutation_num_by_case_temp_perm_i)
    }
    permutation_num_temp <- rbind(permutation_num_temp, permutation_num_temp_perm_i)

    chr_length <- seqlengths(Hsapiens)
    seqlengths(permutation_decile_grange_list_temp_perm_i) <- chr_length[names(seqlengths(permutation_decile_grange_list_temp_perm_i))]
    seqlevels(permutation_decile_grange_list_temp_perm_i) <- seqlevels(permutation_decile_grange_list_temp_perm_i)[order(factor(seqlevels(permutation_decile_grange_list_temp_perm_i), levels = chr_orders))]
    genome(permutation_decile_grange_list_temp_perm_i) = 'hg19'
    
    # calculate 96 type of snvs
    permutation_mut_mat_temp_perm_i <- mut_matrix(permutation_decile_grange_list_temp_perm_i, ref_genome = ref_genome) / batch_size
    # permutation_mut_mat_temp_perm_i <- mut_matrix(permutation_decile_grange_list_temp_perm_i, ref_genome = ref_genome)
    permutation_mut_mat_temp_perm_i_save <- as.data.frame(t(permutation_mut_mat_temp_perm_i)) %>% 
      mutate(perm.id_batch = perm_index) %>% 
      mutate(condition = condition_temp) %>% 
      mutate(decile = sub(".*_", "", rownames(.))) %>% 
      mutate(Case_ID = sub("_.*", "", rownames(.)))
    
    permutation_mut_mat_temp <- rbind(permutation_mut_mat_temp, permutation_mut_mat_temp_perm_i_save)

    if(perm_index %% 1 == 0){
      cat("Condition:", condition_temp, "Permutation index:", perm_index, "\n")
    }
  }
  permutation_num <- rbind(permutation_num, permutation_num_temp)
  permutation_mut_mat_all <- rbind(permutation_mut_mat_all, permutation_mut_mat_temp)
}

analysis_ID <- paste0(permutation_round, "P_", group_num, "G_3sigs_O2")
# analysis_ID <- "5_perms_group_10_3sigs"
fig_save_dir <- paste0("./figures/ATAC_enrichment_analysis/", analysis_ID, "/")
dir.create(fig_save_dir)

saveRDS(mutation_mut_mat_all, paste0(fig_save_dir, "/mutation_mut_mat_all_by_case_decile.rds"))
saveRDS(permutation_mut_mat_all, paste0(fig_save_dir, "/permutation_mut_mat_all_by_case_decile.rds"))
write.csv(mutation_mut_mat_all[, 1:96], paste0(fig_save_dir, "mutation_mut_mat_all_by_case_decile.csv"))
write.csv(permutation_mut_mat_all[, 1:96], paste0(fig_save_dir, "permutation_mut_mat_all_by_case_decile.csv"))

# ### replace all 0s in permutation resutls by 1
# mutation_mut_mat_all_modified <- mutation_mut_mat_all[, 1:96]
# zero_sum_rows_mut <- rowSums(mutation_mut_mat_all_modified) == 0
# mutation_mut_mat_all_modified[zero_sum_rows_mut, ] <- 1
# 
# permutation_mut_mat_all_modified <- permutation_mut_mat_all[, 1:96]
# zero_sum_rows_permut <- rowSums(permutation_mut_mat_all_modified) == 0
# permutation_mut_mat_all_modified[zero_sum_rows_permut, ] <- 1
# 
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
  # ylim(c(0.5, 1.6)) + 
  labs(x = "Gene expression levels", y = "Somatic Mutation enrichment ratio \n (obs/exp)", color = "Condition", title = "")
ggsave(paste0(fig_save_dir, "/total_snv_filtered.pdf"), width = 8, height = 5)
# ggsave(paste0(fig_save_dir, "/total_snv.pdf"), width = 8, height = 5)



################################################################################
################################################################################
##### plot top signature contribution from SigNet for mutation and permutation
################################################################################
################################################################################
# SigNet_contri_mutation <- read.csv(paste0(fig_save_dir, "weight_guesses_SigNet_mutation.csv"), row.names = 1)
# SigNet_contri_permutation <- read.csv(paste0(fig_save_dir, "weight_guesses_SigNet_permutation.csv"), row.names = 1)
SigNet_contri_mutation <- read.csv("/Users/zhemingan/Documents/BCH_research/SigNet/mutation_mut_mat_all_by_case_decile_est/weight_guesses.csv", row.names = 1)
SigNet_contri_permutation <- read.csv("/Users/zhemingan/Documents/BCH_research/SigNet/permutation_mut_mat_all_by_case_decile_est/weight_guesses.csv", row.names = 1)
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
low_mut_num_index <- rowSums(mutation_mut_mat_all_modified_est) <= 0
sum(low_mut_num_index)
SigNet_contri_mutation_est[zero_sum_rows_mut, ] <- NA
# SigNet_contri_mutation_est[SigNet_contri_mutation_score < 0.1, ] <- NA
SigNet_contri_mutation_est[low_mut_num_index, ] <- NA
SigNet_contri_permutation_est[zero_sum_rows_permut, ] <- NA
# SigNet_contri_permutation_est[SigNet_contri_permutation_score < 0.1, ] <- NA
# signature_list <- c("SBS5", "SBS30", "SBS19", "SBS4", "SBS1")
# signature_list <- c("SBS5", "SBS30", "SBS19", "SBS4", "SBS1", "SBS2", "SBS92", "SBS32", "SBS8", "SBS40")
# signature_list <- c("SBS5", "SBS30", "SBS19", "SBS4", "SBS1", "SBS2", "SBS92", "SBS32", "SBS8", "SBS40",
#                     "SBS44", "SBS29", "SBS7b", "SBS89", "SBS36")
# signature_list <- c("SBS5", "SBS30", "SBS19", "SBS4", "SBS1", "SBS2", "SBS92", "SBS32", "SBS8", "SBS40",
#                     "SBS44", "SBS29", "SBS7b", "SBS89", "SBS36", "SBS3", "SBS22", "SBS16", "SBS39", "SBS18")
# signature_list <- colnames(SigNet_contri_mutation_est)
signature_list <- c("SBS5", "SBS30", "SBS19", "SBS4", "SBS1", "SBS2", "SBS92", "SBS32", "SBS8", "SBS40",
                    "SBS44", "SBS29", "SBS7b", "SBS89", "SBS36", "SBS3")

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
sig_percentage_threshold <- 9
# sig_percentage_sum_threshold <- 10

control_sig_contri <- as.data.frame(colSums(SigNet_contri_mutation_plot
                                            [SigNet_contri_mutation_plot$Condition == "Control", 1 : length(signature_list)], na.rm = TRUE)) %>% 
  setNames(c("Control_sig_contri")) %>% mutate(Control_Percentage = Control_sig_contri / sum(Control_sig_contri) * 100)
disease_sig_contri <- as.data.frame(colSums(SigNet_contri_mutation_plot
                                            [SigNet_contri_mutation_plot$Condition == "IHD", 1 : length(signature_list)], na.rm = TRUE)) %>% 
  setNames(c("IHD_sig_contri")) %>% mutate(IHD_Percentage = IHD_sig_contri / sum(IHD_sig_contri) * 100)

all_sig_contri <- cbind(control_sig_contri, disease_sig_contri) %>% 
  mutate(Percentage_Sum = Control_Percentage + IHD_Percentage) |> base::`[`(c("Control_Percentage", "IHD_Percentage", "Percentage_Sum"))

# selected_sigs <- all_sig_contri[all_sig_contri$Control_Percentage >= sig_percentage_threshold | 
#                                   all_sig_contri$IHD_Percentage >= sig_percentage_threshold & 
#                                   all_sig_contri$Percentage_Sum >= sig_percentage_sum_threshold, ]
selected_sigs <- all_sig_contri[all_sig_contri$Control_Percentage >= sig_percentage_threshold |
                                  all_sig_contri$IHD_Percentage >= sig_percentage_threshold, ]
cat(rownames(selected_sigs))
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
# a <- SigNet_contri_all[SigNet_contri_all$condition == "IHD" & SigNet_contri_all$decile == 7,]
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
    mutate(EnR_sig = mut_sig / permut_sig) %>% 
    # filter(EnR_sig != 0) %>%
    # filter(!is.na(EnR_sig)) %>% 
    # mutate(EnR_sig = ifelse(is.infinite(EnR_sig), 1, EnR_sig)) %>% 
    filter(!is.na(EnR_sig) & is.finite(EnR_sig)) %>%
    # mutate(EnR_sig = DescTools::Winsorize(EnR_sig, probs = c(0.05, 0.95))) %>%
    filter(EnR_sig >= quantile(EnR_sig, 0.25) - 1.5 * IQR(EnR_sig) & EnR_sig <= quantile(EnR_sig, 0.75) + 1.5 * IQR(EnR_sig)) %>%
    # filter(EnR_sig >= quantile(EnR_sig, 0.05) & EnR_sig <= quantile(EnR_sig, 0.95)) %>%
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
    summarise(mean_EnR = mean(mean_EnR), sd_EnR = sd(sd_EnR)) %>% 
    mutate(condition = factor(condition, level = c("Control", "IHD")))
  # filter(condition == c("Control", "IHD")[c(selected_sigs_control[sig_id_index], selected_sigs_IHD[sig_id_index])])
  
  ggplot(SigNet_contri_all_final, aes(x = decile, y = mean_EnR, color = condition, group = condition)) + 
    geom_hline(yintercept = 1, color = "black", linewidth = 0.6) + 
    geom_line(position = position_dodge(width = 0.1), size = 1) + 
    geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE) + 
    geom_errorbar(aes(ymin = mean_EnR - sd_EnR, ymax = mean_EnR + sd_EnR), width = 0.2, position = position_dodge(width = 0.1)) +
    geom_smooth(data = SigNet_contri_all_final, aes(x = decile, y = mean_EnR, color = condition, fill = condition, group = condition),
                method = "lm", se = TRUE, alpha = 0.2, linewidth = 1, linetype = "dashed") +
    theme(text = element_text(size = 20), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position="right") + 
    theme_classic() + scale_color_manual(values = color_set) + scale_fill_manual(values = color_set) + 
    # ylim(c(0.5, 1.4)) + 
    labs(x = "Gene expression levels", y = paste0(sig_id, " enrichment ratio \n (obs/exp)"), color = "condition")
  ggsave(paste0(fig_save_dir, sig_id, "_enrichment_filtered.pdf"), width = 8, height = 5)
}