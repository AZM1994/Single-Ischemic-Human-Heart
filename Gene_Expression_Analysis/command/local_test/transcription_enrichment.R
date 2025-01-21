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
perm_round_select <- 1 : 3
# perm_round_select <- sample(1:1000, 10)

##### read in metadata
Hypoxia_PTA_Cases_metadata <- readRDS("./data/SCAN2_df.rds") %>%
  as.data.frame() |> base::`[`(c("cell_ID", "age", "gender", "Case_ID", "condition", "snv.burden", "snv.rate.per.gb")) %>%
  rename_with(~ c("Cell_ID", "Condition"), .cols = c(1, 5))
Hypoxia_PTA_Cases_metadata_collapsed <- Hypoxia_PTA_Cases_metadata %>% distinct(Case_ID, .keep_all = TRUE)
Cell_ID_list <- Hypoxia_PTA_Cases_metadata$Cell_ID
Case_ID_list <- Hypoxia_PTA_Cases_metadata_collapsed$Case_ID
Condition_list <- unique(Hypoxia_PTA_Cases_metadata$Condition)
genomic_context_colnames <- c("Cell_ID", "Case_ID", "Condition", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene")

##### read in scRNA-seq data
all_celltype_RNAseq <- readRDS("./data/Seurat.obj_with_annotation.RDS")
CM_cells <- subset(all_celltype_RNAseq, subset = annotated_clusters == "Cardiomyocytes")

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
# condition_temp = Condition_list[1]
for (condition_temp in Condition_list){
  ##### get transcription data
  cat("##### Get transcription data for:", condition_temp, "...\n")
  expr_level_temp <- data.frame(AverageExpression(CM_cells, group.by = "condition", slot = "data")$RNA) %>% 
    setNames(Condition_list) |> base::`[`(condition_temp) %>% 
    mutate(gene = row.names(.)) %>% 
    setNames(c("average_expr_level", "gene")) %>% 
    mutate(decile = ntile(average_expr_level, n = group_num)) %>% mutate(decile = as.factor(decile))
  
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

  # genic_mutation_temp <- genomic_SCAN2_df_temp
  genic_mutation_temp <- genomic_SCAN2_df_temp[genomic_SCAN2_df_temp$region %in% c("exonic", "exonic;splicing", "intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3"), ] %>%
    mutate(gene = str_remove(gene, "\\(.*\\)$")) %>% filter(!str_detect(gene, ","))
  # table(genomic_SCAN2_df_temp$Case_ID)
  # table(genic_mutation_temp$Case_ID)
  mutation_num_temp <- data.frame(table(genic_mutation_temp$gene, genic_mutation_temp$Case_ID)) %>% 
    setNames(c("gene", "Case_ID", "mut_number"))
  
  expr_level_mutation_temp <- inner_join(expr_level_temp, mutation_num_temp, by = "gene")
  
  cat("Summarize raw call mut_num for", condition_temp, "...\n")
  mutation_mut_num_temp <- expr_level_mutation_temp %>% 
    group_by(decile, Case_ID) %>% 
    summarise("mutation_number" = sum(mut_number)) %>% 
    mutate(condition = condition_temp) %>% 
    arrange(match(Case_ID, Case_ID_order))
  
  ###########################################################################
  ###################### Permutation mutation analysis ######################
  ###########################################################################
  cat("##### Permutation mutation analysis:", condition_temp, "...\n")
  cat("Get permutation data for", condition_temp, "...\n")
  if (condition_temp == "Normal"){
    # permutation_temp <- permutation_normal %>% filter(perm.id_batch %in% seq(1:permutation_round))
    permutation_temp <- permutation_normal %>% filter(perm.id_batch %in% perm_round_select)
    # permutation_temp <- permutation_normal
  }else if (condition_temp == "Disease"){
    # permutation_temp <- permutation_disease %>% filter(perm.id_batch %in% seq(1:permutation_round))
    permutation_temp <- permutation_disease %>% filter(perm.id_batch %in% perm_round_select)
    # permutation_temp <- permutation_disease
  }
  
  # genic_permutation_temp <- permutation_temp
  genic_permutation_temp <- permutation_temp[permutation_temp$region %in% c("exonic", "exonic;splicing", "intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3"),] %>%
    mutate(gene = str_remove(gene, "\\(.*\\)$")) %>% filter(!str_detect(gene, ","))
  permutation_num_temp <- data.frame(table(genic_permutation_temp$gene, genic_permutation_temp$Case_ID, genic_permutation_temp$perm.id_batch)) %>%
    setNames(c("gene", "Case_ID", "perm.id_batch", "permutation_number"))
  # table(permutation_temp$Case_ID)
  # table(genic_permutation_temp$Case_ID)
  
  expr_level_permutation_temp <- inner_join(expr_level_temp, permutation_num_temp, by = "gene")

  permutation_mut_num_temp <- expr_level_permutation_temp %>% group_by(decile, Case_ID, perm.id_batch) %>%
    summarise("permutation_number" = sum(permutation_number) / batch_size) %>%
    mutate(condition = condition_temp) %>% 
    arrange(match(Case_ID, Case_ID_order))
  
  cat("Summarize permutation mut_num for", condition_temp, "...\n")
  mut_num_genic_merged_temp <- merge(mutation_mut_num_temp, permutation_mut_num_temp) %>% group_by(Case_ID) %>% 
    summarise("mutation_number" = sum(mutation_number), "permutation_number" = sum(permutation_number) / batch_size) %>% 
    mutate(enrichment_ratio = mutation_number/permutation_number) %>% 
    mutate(condition = as.factor(condition_temp)) %>% 
    arrange(match(Case_ID, Case_ID_order))
  
  mut_num_genic_merged <- rbind(mut_num_genic_merged, mut_num_genic_merged_temp) 
  
  ###########################################################################
  ############################ signature analysis ###########################
  ###########################################################################
  ##### raw SCAN2 call signature analysis
  cat("Raw SCAN2 call mutation signature analysis for", condition_temp, "...\n")
  mutation_decile_grange_list_temp <- GRangesList()
  for (case_index in Case_ID_order){
    for (i in 1:group_num){
      ## select mutation in genes belonging to each decile for each case
      expr_level_mutation_temp_percase <- expr_level_mutation_temp[expr_level_mutation_temp$decile == i & expr_level_mutation_temp$Case_ID == case_index, ]
      gene_name_temp <- expr_level_mutation_temp_percase$gene[expr_level_mutation_temp_percase$mut_number != 0]
      # gene_name_temp <- expr_level_mutation_temp[expr_level_mutation_temp$decile == i & expr_level_mutation_temp$Case_ID == case_index, "gene"]
      genic_mutation_decile_case_temp <- genic_mutation_temp[genic_mutation_temp$gene %in% gene_name_temp & genic_mutation_temp$Case_ID == case_index, ]
      
      ## change the order
      genic_mutation_decile_case_temp <- genic_mutation_decile_case_temp[order(genic_mutation_decile_case_temp$start),]
      genic_mutation_decile_case_temp <- genic_mutation_decile_case_temp[order(genic_mutation_decile_case_temp$chr),]
      
      if (dim(genic_mutation_decile_case_temp)[1] == 0){
        mutation_decile_grange_list_temp[[paste0(case_index, "_", i)]] <- GRanges()
      }else{
        mutation_decile_grange_list_temp[[paste0(case_index, "_", i)]] <- 
          GRanges(seqnames = paste0("chr", genic_mutation_decile_case_temp$chr),
                  ranges = IRanges(start = genic_mutation_decile_case_temp$start,
                                   end = genic_mutation_decile_case_temp$end),
                  ref = genic_mutation_decile_case_temp$ref,
                  alt = genic_mutation_decile_case_temp$alt)
      }
    }
  }
  
  chr_length <- seqlengths(Hsapiens)
  seqlengths(mutation_decile_grange_list_temp) <- chr_length[names(seqlengths(mutation_decile_grange_list_temp))]
  seqlevels(mutation_decile_grange_list_temp) <- seqlevels(mutation_decile_grange_list_temp)[order(factor(seqlevels(mutation_decile_grange_list_temp), levels = chr_orders))]
  genome(mutation_decile_grange_list_temp) = 'hg19'
  
  ##### calculate 96 types of snvs
  mutation_mut_mat_temp <- mut_matrix(mutation_decile_grange_list_temp, ref_genome = ref_genome)
  # mutation_mut_mat_temp <- 1 / colSums(mutation_mut_mat_temp) * mutation_mut_mat_temp
  mutation_mut_mat_temp_save <- as.data.frame(t(mutation_mut_mat_temp)) %>% 
    mutate(condition = condition_temp) %>% 
    mutate(decile = sub(".*_", "", rownames(.))) %>% 
    mutate(Case_ID = sub("_.*", "", rownames(.)))

  # row.names(mutation_mut_mat_temp_save) <- seq(1, group_num)
  mutation_mut_mat_all <- rbind(mutation_mut_mat_all, mutation_mut_mat_temp_save)
  
  ##### calculate signature
  # selected_sig_list <- read.table("./data/COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE) %>%
  #   column_to_rownames(var = "Type") |> base::`[`(signature_list) %>% as.matrix()
  # mutation_sig_fit_temp <- fit_to_signatures(mutation_mut_mat_temp, selected_sig_list)
  # mutation_sig_contri_temp <- apply(mutation_sig_fit_temp$contribution, 2, function(x){x/sum(x)})
  
  # for (sbs_sig_index in 1:length(signature_list)){
  #   mutation_mut_num_temp[paste0(signature_list[sbs_sig_index], "_mutation")] <- rep(mutation_sig_contri_temp[sbs_sig_index, ], each = permutation_round)
  # }
  
  ##### Permutation signature analysis
  cat("Permutation signature analysis for", condition_temp, "...\n")
  # permutation_sig_contri_temp <- c()
  permutation_mut_mat_temp <- c()
  # for (perm_index in 1:permutation_round){
  for (perm_index in perm_round_select){
    permutation_decile_grange_list_temp_perm_i <- GRangesList()
    for (case_index in Case_ID_order){
      for (i in 1:group_num){
        # cat("Condition:", condition_temp, "Permutation index:", perm_index, "decile", i, "\n")
        # select mutation in genes belonging to each decile
        expr_level_permutation_temp_percase_perm_i <- expr_level_permutation_temp[expr_level_permutation_temp$decile == i & expr_level_permutation_temp$Case_ID == case_index, ]
        gene_name_temp <- expr_level_permutation_temp_percase_perm_i$gene[expr_level_permutation_temp_percase_perm_i$permutation_number != 0]
        
        genic_permutation_decile_case_temp_perm_i <- genic_permutation_temp[(genic_permutation_temp$gene %in% gene_name_temp) &
                                                                  (genic_permutation_temp$perm.id_batch == perm_index) & 
                                                                    (genic_permutation_temp$Case_ID == case_index), ]
        
        # change the order
        genic_permutation_decile_case_temp_perm_i <- genic_permutation_decile_case_temp_perm_i[order(genic_permutation_decile_case_temp_perm_i$start), ]
        genic_permutation_decile_case_temp_perm_i <- genic_permutation_decile_case_temp_perm_i[order(genic_permutation_decile_case_temp_perm_i$chr), ]
  
        if (dim(genic_permutation_decile_case_temp_perm_i)[1] == 0){
          permutation_decile_grange_list_temp_perm_i[[paste0(case_index, "_", i)]] <- GRanges()
        }else{
          permutation_decile_grange_list_temp_perm_i[[paste0(case_index, "_", i)]] <-
            GRanges(seqnames = genic_permutation_decile_case_temp_perm_i$chr,
                    ranges = IRanges(start = genic_permutation_decile_case_temp_perm_i$start,
                                     end = genic_permutation_decile_case_temp_perm_i$end),
                    ref = genic_permutation_decile_case_temp_perm_i$ref,
                    alt = genic_permutation_decile_case_temp_perm_i$alt)
        }
      }
      }

    chr_length <- seqlengths(Hsapiens)
    seqlengths(permutation_decile_grange_list_temp_perm_i) <- chr_length[names(seqlengths(permutation_decile_grange_list_temp_perm_i))]
    seqlevels(permutation_decile_grange_list_temp_perm_i) <- seqlevels(permutation_decile_grange_list_temp_perm_i)[order(factor(seqlevels(permutation_decile_grange_list_temp_perm_i), levels = chr_orders))]
    genome(permutation_decile_grange_list_temp_perm_i) = 'hg19'
    
    # calculate 96 type of snvs
    permutation_mut_mat_temp_perm_i <- mut_matrix(permutation_decile_grange_list_temp_perm_i, ref_genome = ref_genome) / batch_size
    # permutation_mut_mat_temp_perm_i <- 1 / colSums(permutation_mut_mat_temp_perm_i) * permutation_mut_mat_temp_perm_i
    permutation_mut_mat_temp_perm_i_save <- as.data.frame(t(permutation_mut_mat_temp_perm_i)) %>% 
      mutate(perm.id_batch = perm_index) %>% 
      mutate(condition = condition_temp) %>% 
      mutate(decile = sub(".*_", "", rownames(.))) %>% 
      mutate(Case_ID = sub("_.*", "", rownames(.)))
    
    permutation_mut_mat_temp <- rbind(permutation_mut_mat_temp, permutation_mut_mat_temp_perm_i_save)

    # permutation_sig_fit_temp_perm_i <- fit_to_signatures(permutation_mut_mat_temp_perm_i, selected_sig_list)
    # permutation_sig_contri_temp_perm_i <- apply(permutation_sig_fit_temp_perm_i$contribution, 2, function(x){x/sum(x)})

    # permutation_sig_contri_temp_perm_i <- data.frame(t(permutation_sig_contri_temp_perm_i))
    # row.names(permutation_sig_contri_temp_perm_i) <- seq(1, group_num)
    # colnames(permutation_sig_contri_temp_perm_i) <- str_c(colnames(permutation_sig_contri_temp_perm_i), "_permutation")
    # permutation_sig_contri_temp_perm_i <- permutation_sig_contri_temp_perm_i %>%
    #   mutate(perm.id_batch = perm_index) %>%
    #   mutate(decile = seq(1, group_num)) %>%
    #   mutate(condition = condition_temp)

    # permutation_sig_contri_temp <- rbind(permutation_sig_contri_temp, permutation_sig_contri_temp_perm_i)
    if(perm_index %% 1 == 0){
      # print(paste("perm_index is a multiple of 10:", perm_index))
      cat("Condition:", condition_temp, "Permutation index:", perm_index, "\n")
    }
  }
  permutation_mut_mat_all <- rbind(permutation_mut_mat_all, permutation_mut_mat_temp)
  # permutation_sig_contri_temp$decile <- as.factor(permutation_sig_contri_temp$decile)
  # mutation_mut_num_temp$perm.id_batch <- as.integer(as.character(mutation_mut_num_temp$perm.id_batch))
  # mut_num_mutation_temp_combine <- left_join(mutation_mut_num_temp, permutation_sig_contri_temp, by = c("decile", "perm.id_batch", "condition"))
  # mut_num_mutation_temp_combine <- mutation_mut_num_temp

  # signature_list <- colnames(selected_sig_list)
  # for (sbs_sig in signature_list){
  #   mut_num_mutation_temp_combine[paste0(sbs_sig, "_enrichment_ratio")] <-
  #     (mut_num_mutation_temp_combine[,paste0(sbs_sig, "_mutation")] * mut_num_mutation_temp_combine$mutation_number) /
  #     (mut_num_mutation_temp_combine[,paste0(sbs_sig, "_permutation")] * mut_num_mutation_temp_combine$permutation_number)
  # }
  # average_mut_num_all <- rbind(average_mut_num_all, mut_num_mutation_temp_combine)
  # average_mut_num_all <- rbind(average_mut_num_all, mutation_mut_num_temp)
}

# analysis_ID <- paste0(permutation_round, "P_", group_num, "G")
analysis_ID <- "test"
fig_save_dir <- paste0("./results/local_test/RNAseq_results/", analysis_ID, "/")
dir.create(fig_save_dir, recursive = TRUE)

# saveRDS(mutation_mut_mat_all, paste0(fig_save_dir, "/mutation_mut_mat_all_by_case_decile.rds"))
saveRDS(mutation_mut_mat_all, paste0(fig_save_dir, "/mutation_mut_mat_all_by_case_decile.rds"))
saveRDS(permutation_mut_mat_all, paste0(fig_save_dir, "/permutation_mut_mat_all_by_case_decile.rds"))
write.csv(mutation_mut_mat_all[, 1:96], paste0(fig_save_dir, "mutation_mut_mat_all_by_case_decile.csv"))
write.csv(permutation_mut_mat_all[, 1:96], paste0(fig_save_dir, "permutation_mut_mat_all_by_case_decile.csv"))

# mutation_mut_mat_all <- readRDS(paste0(fig_save_dir, "/mutation_mut_mat_all_by_case_decile.rds"))
# permutation_mut_mat_all <- readRDS(paste0(fig_save_dir, "/permutation_mut_mat_all_by_case_decile.rds"))

### replace all 0s in permutation resutls by 1
mutation_mut_mat_all_modified <- mutation_mut_mat_all[, 1:96]
zero_sum_rows_mut <- rowSums(mutation_mut_mat_all_modified) == 0
mutation_mut_mat_all_modified[zero_sum_rows_mut, ] <- 1

permutation_mut_mat_all_modified <- permutation_mut_mat_all[, 1:96]
zero_sum_rows_permut <- rowSums(permutation_mut_mat_all_modified) == 0
permutation_mut_mat_all_modified[zero_sum_rows_permut, ] <- 1

write.csv(mutation_mut_mat_all_modified, paste0(fig_save_dir, "mutation_mut_mat_all_modified.csv"))
write.csv(permutation_mut_mat_all_modified, paste0(fig_save_dir, "permutation_mut_mat_all_modified.csv"))

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

write.csv(mut_num_genic_merged, paste0(fig_save_dir, "mut_num_genic_merged.csv"))
write.csv(mutation_mut_mat_summary_Normal, paste0(fig_save_dir, "mutation_mut_mat_summary_Normal.csv"))
write.csv(mutation_mut_mat_summary_IHD, paste0(fig_save_dir, "mutation_mut_mat_summary_IHD.csv"))
write.csv(mutation_mut_mat_summary, paste0(fig_save_dir, "mutation_mut_mat_summary.csv"))

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

write.csv(permutation_mut_mat_summary_Normal, paste0(fig_save_dir, "permutation_mut_mat_summary_Normal.csv"))
write.csv(permutation_mut_mat_summary_IHD, paste0(fig_save_dir, "permutation_mut_mat_summary_IHD.csv"))
write.csv(permutation_mut_mat_summary, paste0(fig_save_dir, "permutation_mut_mat_summary.csv"))

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
  # filter(!is.na(enrichment_ratio)) %>% 
  # mutate(enrichment_ratio = ifelse(is.infinite(enrichment_ratio), NA, enrichment_ratio)) %>% 
  # mutate(enrichment_ratio = Winsorize(enrichment_ratio, probs = c(0.05, 0.95))) %>%
  filter(enrichment_ratio >= quantile(enrichment_ratio, 0.25) - 1.5 * IQR(enrichment_ratio) &
           enrichment_ratio <= quantile(enrichment_ratio, 0.75) + 1.5 * IQR(enrichment_ratio)) %>%
  group_by(condition, decile) %>% 
  summarise(mean_ER = mean(enrichment_ratio, na.rm = TRUE), 
            sd_ER = sd(enrichment_ratio, na.rm = TRUE)) %>% 
  # filter(!is.na(mean_ER)) %>% 
  # filter(mean_ER >= quantile(mean_ER, 0.25) - 1.5 * IQR(mean_ER) & mean_ER <= quantile(mean_ER, 0.75) + 1.5 * IQR(mean_ER)) %>% 
  # mutate(mean_ER = ifelse(is.infinite(mean_ER), 1, mean_ER)) %>% 
  mutate(Condition = ifelse(condition == "Normal", "Control", ifelse(condition == "Disease", "IHD", "Unknown"))) %>% 
  mutate(Condition = factor(Condition, level = c("Control", "IHD"))) %>% 
  mutate(decile = factor(decile, level = seq(1:group_num)))

merged_summary_plot <- merge(mutation_mut_mat_summary_plot0, permutation_mut_mat_summary_plot0) %>% 
  merge(Hypoxia_PTA_Cases_metadata_collapsed[, c("Case_ID", "age")]) %>% 
  # filter(age >= 0.5 & age <= 90) %>%
  mutate(enrichment_ratio = mut_num / permut_sum) %>% 
  filter(!is.na(enrichment_ratio) & is.finite(enrichment_ratio)) %>% 
  # filter(!is.na(enrichment_ratio)) %>% 
  # mutate(enrichment_ratio = ifelse(is.infinite(enrichment_ratio), NA, enrichment_ratio)) %>% 
  # mutate(enrichment_ratio = Winsorize(enrichment_ratio, probs = c(0.05, 0.95))) %>%
  # filter(enrichment_ratio >= quantile(enrichment_ratio, 0.25) - 1.5 * IQR(enrichment_ratio) &
  #          enrichment_ratio <= quantile(enrichment_ratio, 0.75) + 1.5 * IQR(enrichment_ratio)) %>% 
  filter(enrichment_ratio >= quantile(enrichment_ratio, 0.05) & enrichment_ratio <= quantile(enrichment_ratio, 0.95)) %>% 
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
  ylim(c(0.6, 1.6)) + 
  labs(x = "Gene expression levels", y = "Somatic Mutation enrichment ratio \n (obs/exp)", color = "Condition", title = "")
ggsave(paste0(fig_save_dir, "/total_snv_filtered.pdf"), width = 8, height = 5)
# ggsave(paste0(fig_save_dir, "/total_snv.pdf"), width = 8, height = 5)

##### plot top signature contribution from SigNet for mutation and permutation
SigNet_contri_mutation <- read.csv(paste0(fig_save_dir, "weight_guesses_SigNet_mutation.csv"), row.names = 1)
SigNet_contri_permutation <- read.csv(paste0(fig_save_dir, "weight_guesses_SigNet_permutation.csv"), row.names = 1)
SigNet_contri_mutation[zero_sum_rows_mut, ] <- NA
SigNet_contri_permutation[zero_sum_rows_permut, ] <- NA
signature_list <- c("SBS1", "SBS5", "SBS30")

SigNet_contri_mutation_plot <- as.data.frame(SigNet_contri_mutation) %>% 
  mutate(decile = sub(".*_", "", rownames(.))) %>% 
  mutate(Case_ID = sub("_.*", "", rownames(.))) %>% 
  mutate(Condition = ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
    Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Normal", "Case_ID"], "Control", 
    ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
      Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Disease", "Case_ID"], "IHD", "Unknown"))) |> 
  base::`[`(c(signature_list, "decile", "Case_ID", "Condition")) %>% 
  setNames(c(paste0("mut_", signature_list), "decile", "Case_ID", "Condition"))

df_A <- SigNet_contri_mutation_plot %>% filter(Condition == 'Control')
df_B <- SigNet_contri_mutation_plot %>% filter(Condition == 'IHD')
df_A_repeated <- df_A[rep(seq_len(nrow(df_A)), permutation_round), ]
df_B_repeated <- df_B[rep(seq_len(nrow(df_B)), permutation_round), ]
df_combined <- rbind(df_A_repeated, df_B_repeated)

SigNet_contri_permutation_plot <- as.data.frame(SigNet_contri_permutation) %>% 
  mutate(decile = permutation_mut_mat_all$decile) %>% 
  mutate(Case_ID = sub("_.*", "", rownames(.))) |> 
  mutate(Condition = ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
    Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Normal", "Case_ID"], "Control", 
    ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
      Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Disease", "Case_ID"], "IHD", "Unknown"))) |> 
  base::`[`(c(signature_list, "decile", "Case_ID", "Condition")) %>% 
  setNames(c(paste0("permut_", signature_list), "decile", "Case_ID", "Condition"))

# sig_id = "SBS5"
for (sig_id in signature_list){
  SigNet_contri_all <- cbind(df_combined[c(paste0("mut_", sig_id), "decile", "Case_ID", "Condition")], 
                             SigNet_contri_permutation_plot[c(paste0("permut_", sig_id))]) %>% 
    setNames(c("mut_sig", "decile", "Case_ID", "Condition", "permut_sig")) %>% 
    merge(Hypoxia_PTA_Cases_metadata_collapsed[, c("Case_ID", "age")]) %>% 
    # filter(age >= 0.5 & age <= 90) %>%
    mutate(EnR_sig = mut_sig / permut_sig) %>% 
    # filter(!is.na(EnR_sig)) %>% 
    # mutate(EnR_sig = ifelse(is.infinite(EnR_sig), 1, EnR_sig)) %>% 
    filter(!is.na(EnR_sig) & is.finite(EnR_sig)) %>%
    # mutate(EnR_sig = Winsorize(EnR_sig, probs = c(0.05, 0.95))) %>%
    # filter(EnR_sig >= quantile(EnR_sig, 0.25) - 1.5 * IQR(EnR_sig) & EnR_sig <= quantile(EnR_sig, 0.75) + 1.5 * IQR(EnR_sig)) %>%
    filter(EnR_sig >= quantile(EnR_sig, 0.05) & EnR_sig <= quantile(EnR_sig, 0.95)) %>% 
    group_by(Case_ID, decile) %>% 
    summarise(mean_EnR = mean(EnR_sig), sd_EnR = sd(EnR_sig)) %>% 
    mutate(condition = ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
      Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Normal", "Case_ID"], "Control", 
      ifelse(Case_ID %in% Hypoxia_PTA_Cases_metadata_collapsed[
        Hypoxia_PTA_Cases_metadata_collapsed$Condition == "Disease", "Case_ID"], "IHD", "Unknown"))) %>% 
    mutate(decile = factor(decile, level = seq(1:group_num)))
  
  SigNet_contri_all_final <- SigNet_contri_all %>% 
    group_by(decile, condition) %>% 
    summarise(mean_EnR = mean(mean_EnR), sd_EnR = sd(sd_EnR)) %>% 
    mutate(condition = factor(condition, level = c("Control", "IHD")))
  
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
    ylim(c(0.6, 1.3)) +
    labs(x = "Gene expression levels", y = paste0(sig_id, " enrichment ratio \n (obs/exp)"), color = "condition")
  ggsave(paste0(fig_save_dir, sig_id, "_enrichment_filtered.pdf"), width = 8, height = 5)
}