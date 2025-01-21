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
ref_genome="BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = T)
library(MutationalPatterns)

args <- commandArgs(trailingOnly = TRUE)

chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

setwd("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Enrichment_Analysis/transcription_analysis")
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])
group_num <- as.numeric(args[1])
batch_size <- 1
permutation_round <- 1000 / batch_size
run_index <- as.numeric(args[2])
perm_round_select <- ((run_index - 1) * 100 + 1) : (run_index * 100)

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
  
  genic_mutation_temp <- genomic_SCAN2_df_temp[genomic_SCAN2_df_temp$region %in% c("exonic", "exonic;splicing", "intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3"), ] %>%
    mutate(gene = str_remove(gene, "\\(.*\\)$")) %>% filter(!str_detect(gene, ","))
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
    permutation_temp <- permutation_normal %>% filter(perm.id_batch %in% perm_round_select)
  }else if (condition_temp == "Disease"){
    permutation_temp <- permutation_disease %>% filter(perm.id_batch %in% perm_round_select)
  }
  
  genic_permutation_temp <- permutation_temp[permutation_temp$region %in% c("exonic", "exonic;splicing", "intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3"),] %>%
    mutate(gene = str_remove(gene, "\\(.*\\)$")) %>% filter(!str_detect(gene, ","))
  permutation_num_temp <- data.frame(table(genic_permutation_temp$gene, genic_permutation_temp$Case_ID, genic_permutation_temp$perm.id_batch)) %>%
    setNames(c("gene", "Case_ID", "perm.id_batch", "permutation_number"))
  
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
  mutation_mut_mat_temp_save <- as.data.frame(t(mutation_mut_mat_temp)) %>% 
    mutate(condition = condition_temp) %>% 
    mutate(decile = sub(".*_", "", rownames(.))) %>% 
    mutate(Case_ID = sub("_.*", "", rownames(.)))
  
  mutation_mut_mat_all <- rbind(mutation_mut_mat_all, mutation_mut_mat_temp_save)
  
  ##### Permutation signature analysis
  cat("Permutation signature analysis for", condition_temp, "...\n")
  permutation_mut_mat_temp <- c()
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
  permutation_mut_mat_all <- rbind(permutation_mut_mat_all, permutation_mut_mat_temp)
}

analysis_ID <- paste0(permutation_round, "_perms_", group_num, "G_run_", run_index)
fig_save_dir <- paste0("./figures/transcription_analysis/", analysis_ID, "/")
dir.create(fig_save_dir)

saveRDS(mutation_mut_mat_all, paste0(fig_save_dir, "/mutation_mut_mat_all_by_case_decile.rds"))
saveRDS(permutation_mut_mat_all, paste0(fig_save_dir, "/permutation_mut_mat_all_by_case_decile.rds"))
write.csv(mutation_mut_mat_all[, 1:96], paste0(fig_save_dir, "mutation_mut_mat_all_by_case_decile.csv"))
write.csv(permutation_mut_mat_all[, 1:96], paste0(fig_save_dir, "permutation_mut_mat_all_by_case_decile.csv"))