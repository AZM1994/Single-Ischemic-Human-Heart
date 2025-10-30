# library(ggplot2)
library(stringr)
# library(purrr)
# library(readxl)
# library(ggsci)
library(dplyr)
library(tidyr)
# library(tibble)
# library(reshape2)
library(Seurat)
# library(pheatmap)
# library(readxl)
ref_genome="BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = T)
library(MutationalPatterns)
library(data.table)

setwd("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration/Mutation_Enrichment_Analysis/snRNA_ATACseq_Enrichment")
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
group_num <- 8
total_permutation <- 10000
permutation_num <- 1000
Case_selected <- 1864
# Case_selected <- c(4402, 1864, 6032, 4638, 1863, 1039, 1940, 5919, 5828, 5657, 1363, 1673, 604, 1743, 1113)

##### read in metadata
Hypoxia_PTA_all_SCAN2 <- readRDS("./data/1-sSNV_SCAN2_df_filtered.rds") %>% 
  dplyr::select(c("Cell_ID", "Age", "Gender", "Case_ID", "Condition", "snv.burden", "snv.rate.per.gb"))
Hypoxia_PTA_all_SCAN2_collapsed <- Hypoxia_PTA_all_SCAN2 %>% 
  distinct(Case_ID, .keep_all = TRUE)
# Cell_ID_list <- Hypoxia_PTA_all_SCAN2$Cell_ID
Case_ID_list <- Hypoxia_PTA_all_SCAN2_collapsed$Case_ID
Condition_list <- unique(Hypoxia_PTA_all_SCAN2$Condition)
genomic_context_colnames <- c("Cell_ID", "Case_ID", "Condition", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene")

##### read in CM scRNA-seq data
seurat_obj_CM <- readRDS("./data/Seurat.obj_sub_clustering_CM_only.RDS")

##### read in permutation results for selected Case_ID
for (Case_ID_temp in Case_selected){
  permutation_temp <- c()
  Cell_ID_list <- Hypoxia_PTA_all_SCAN2 %>% filter(Case_ID == Case_ID_temp) %>% pull(Cell_ID)
  cat("Cell_ID_list:", Cell_ID_list)
  
  for (Cell_ID_temp in Cell_ID_list){
    permutation_path <- paste0("./data/permutation_SCAN2/sSNV/", Cell_ID_temp, "/perms_snv.hg19_multianno.csv")
    case_temp <- Hypoxia_PTA_all_SCAN2$Case_ID[Hypoxia_PTA_all_SCAN2$Cell_ID == Cell_ID_temp]
    condition_temp <- Hypoxia_PTA_all_SCAN2$Condition[Hypoxia_PTA_all_SCAN2$Cell_ID == Cell_ID_temp]
    cat("Cell:", Cell_ID_temp, "Case:", case_temp, "Condition:", condition_temp, "\n")
    
    total_lines <- as.integer(system(paste("awk 'END{print NR}'", shQuote(permutation_path)), intern = TRUE))
    rows_to_read <- floor((total_lines - 1) * permutation_num / total_permutation)
    permutation_cell <- fread(permutation_path, nrows = rows_to_read) %>% as.data.table() %>% 
    # permutation_cell <- read.csv(permutation_path, header = TRUE) %>% as.data.table() %>% 
      .[, perm.id := rep(1:permutation_num, each = .N / permutation_num)] %>%
      .[, `:=`(Gene.refGene = {
        gene_lists <- strsplit(Gene.refGene, ";")
        dist_lists <- str_extract_all(GeneDetail.refGene, "\\d+") %>% lapply(as.integer)
        mapply(function(gs, ds) {
          if ("NONE" %in% gs && length(gs) > 1) gs <- gs[gs != "NONE"]
          if (length(ds) == length(gs) && length(ds) > 0) gs[which.min(ds)] else gs[1]
        }, gene_lists, dist_lists)},
        GeneDetail.refGene = {
          dist_lists <- str_extract_all(GeneDetail.refGene, "\\d+") %>% lapply(as.integer)
          sapply(dist_lists, function(ds) if (length(ds)) paste0("dist=", min(ds)) else NA)
        })] %>% 
      dplyr::select(c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "perm.id")) %>% 
      mutate(Cell_ID = Cell_ID_temp, Case_ID = case_temp, Condition = condition_temp)
    
    permutation_temp <- rbind(permutation_temp, permutation_cell)
  }

  ##### Get transcription data
  cat("Get transcription data for:", condition_temp, "...\n")
  expr_level_temp <- data.frame(AverageExpression(seurat_obj_CM, group.by = "condition", slot = "data")$RNA) %>% 
    setNames(Condition_list) %>% 
    dplyr::select(condition_temp) %>% 
    mutate(Gene.refGene = row.names(.)) %>% 
    setNames(c("average_expr_level", "Gene.refGene")) %>% 
    # filter(average_expr_level > 0) %>%
    # filter(between(average_expr_level, quantile(average_expr_level, 0.05), quantile(average_expr_level, 0.95))) %>%
    mutate(decile = ntile(average_expr_level, n = group_num), 
           decile = as.factor(decile))
  
  ##### Get SCAN2 genomic_context data with metadata
  cat("Get raw call genomic_context data for:", Case_ID_temp, "...\n")
  genomic_context_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene")
  heart_PTA_all_vcf_temp <- read.table(paste0("data/annotation_results/all_age_sSNV/", condition_temp, "/heart_PTA_all.all_age.", condition_temp, "_ssnv.vcf"), sep = "\t")
  genomic_context_temp <- read.csv(paste0("data/annotation_results/all_age_sSNV/", condition_temp, "/heart_PTA_all.all_age.", condition_temp, "_ssnv.hg19_multianno.csv"), header = TRUE) %>% 
    as.data.table() %>% 
    .[, `:=`(Gene.refGene = {
      gene_lists <- strsplit(Gene.refGene, ";")
      dist_lists <- str_extract_all(GeneDetail.refGene, "\\d+") %>% lapply(as.integer)
      mapply(function(gs, ds) {
        if ("NONE" %in% gs && length(gs) > 1) gs <- gs[gs != "NONE"]
        if (length(ds) == length(gs) && length(ds) > 0) gs[which.min(ds)] else gs[1]
      }, gene_lists, dist_lists)},
      GeneDetail.refGene = {
        dist_lists <- str_extract_all(GeneDetail.refGene, "\\d+") %>% lapply(as.integer)
        sapply(dist_lists, function(ds) if (length(ds)) paste0("dist=", min(ds)) else NA)
      })] %>% 
    dplyr::select(all_of(genomic_context_colnames)) %>% 
    mutate(Cell_ID = heart_PTA_all_vcf_temp$V8)
  
  genomic_SCAN2_df_temp <- genomic_context_temp %>% 
    inner_join(Hypoxia_PTA_all_SCAN2 %>% select(Cell_ID, Case_ID, Age, Gender, Condition), by = c("Cell_ID")) %>% 
    mutate(Condition = as.factor(Condition)) %>% filter(Case_ID == Case_ID_temp)
  genic_mutation_temp <- genomic_SCAN2_df_temp %>% 
    filter(Func.refGene %in% c("exonic", "exonic;splicing", "intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3"))
  mutation_num_temp <- data.frame(table(genic_mutation_temp$Gene.refGene)) %>% setNames(c("Gene.refGene", "mut_number"))
  
  expr_level_mutation_temp <- inner_join(expr_level_temp, mutation_num_temp, by = "Gene.refGene")
  
  cat("Summarize raw call mut_num for Case_ID:", Case_ID_temp, "...\n")
  mutation_mut_num_temp <- expr_level_mutation_temp %>% 
    group_by(decile) %>%
    summarise("mutation_number" = sum(mut_number)) %>% 
    mutate(Condition = condition_temp)
  
  ###########################################################################
  ##### Permutation mutation analysis
  cat("Get permutation genomic_context data for:", Case_ID_temp, "...\n")
  
  genic_permutation_temp <- permutation_temp %>% filter(Func.refGene %in% c("exonic", "exonic;splicing", "intronic", "splicing", "UTR3", "UTR5", "UTR5;UTR3"))
  permutation_num_temp <- data.frame(table(genic_permutation_temp$Gene.refGene, genic_permutation_temp$perm.id)) %>%
    setNames(c("Gene.refGene", "perm.id", "permutation_number"))
  
  expr_level_permutation_temp <- inner_join(expr_level_temp, permutation_num_temp, by = "Gene.refGene")
  
  cat("Summarize permutation mut_num for Case_ID:", Case_ID_temp, "...\n")
  permutation_mut_num_temp <- expr_level_permutation_temp %>% 
    group_by(decile, perm.id) %>% 
    summarise("permutation_number" = sum(permutation_number)) %>% 
    mutate(Condition = condition_temp)
  
  cat("Merge raw call and permutation mut_num for Case_ID:", Case_ID_temp, "...\n")
  mut_num_genic_merged_temp <- tibble(decile = as.factor(1:group_num)) %>% 
    left_join(permutation_mut_num_temp, by = "decile") %>% 
    left_join(mutation_mut_num_temp, by = c("decile")) %>% 
    group_by(decile) %>% 
    summarise("mutation_number" = sum(mutation_number), "permutation_number" = sum(permutation_number)) %>% 
    mutate(enrichment_ratio = mutation_number/permutation_number, Condition = as.factor(condition_temp)) %>% 
    replace_na(list(permutation_number = 0, mutation_number = 0, enrichment_ratio = 0))

  ###########################################################################
  ##### Calculate 96 trinuclutide context
  cat("Build raw call GRanges for Case_ID:", Case_ID_temp, "...\n")
  mutation_decile_grange_list_temp <- GRangesList()
  for (i in 1:group_num){
    ## select mutation in genes belonging to each decile for each case
    gene_name_temp <- expr_level_mutation_temp %>% filter(decile == i) %>% filter(mut_number != 0) %>% pull(Gene.refGene)
    genic_mutation_decile_case_temp <- genic_mutation_temp %>% filter(Gene.refGene %in% gene_name_temp) %>% arrange(Start, Chr)
    
    if (nrow(genic_mutation_decile_case_temp) == 0) {
      mutation_decile_grange_list_temp[[paste0(Case_ID_temp, "_", i)]] <- GRanges()
    } else {
      gr <- GRanges(
        seqnames = paste0("chr", genic_mutation_decile_case_temp$Chr),
        ranges = IRanges(
          start = genic_mutation_decile_case_temp$Start,
          end = genic_mutation_decile_case_temp$End
        )
      )
      mcols(gr)$Ref <- genic_mutation_decile_case_temp$Ref
      mcols(gr)$Alt <- genic_mutation_decile_case_temp$Alt
      mutation_decile_grange_list_temp[[paste0(Case_ID_temp, "_", i)]] <- gr
    }
  }
  
  chr_length <- seqlengths(Hsapiens)[1:22]
  seqlengths(mutation_decile_grange_list_temp) <- chr_length[names(seqlengths(mutation_decile_grange_list_temp))]
  seqlevels(mutation_decile_grange_list_temp) <- seqlevels(mutation_decile_grange_list_temp)[order(factor(seqlevels(mutation_decile_grange_list_temp), levels = chr_orders))]
  genome(mutation_decile_grange_list_temp) = 'hg19'
  
  ##### Calculate 96 trinuclutide context
  cat("Calculate raw call 96 trinuclutide context for Case_ID:", Case_ID_temp, "...\n")
  mutation_mut_mat_temp <- mut_matrix(mutation_decile_grange_list_temp, ref_genome = ref_genome) %>% 
    t() %>% as.data.frame() %>% 
    mutate(Condition = condition_temp, 
           decile = sub(".*_", "", rownames(.)), 
           Case_ID = sub("_.*", "", rownames(.)))

  ##### Permutation signature analysis
  cat("Build permutation GRanges for Case_ID:", Case_ID_temp, "...\n")
  cat("Calculate permutation 96 trinuclutide context for Case_ID:", Case_ID_temp, "...\n")
  genic_permutation_temp_split <- split(genic_permutation_temp, genic_permutation_temp$perm.id)
  expr_level_permutation_temp_split <- split(expr_level_permutation_temp, expr_level_permutation_temp$perm.id)
  permutation_mut_mat_list <- vector("list", permutation_num)

  for (perm_index in 1:permutation_num){
    permutation_decile_grange_list_temp_perm_i <- GRangesList()
    
    genic_permutation_temp_perm_i <- genic_permutation_temp_split[[as.character(perm_index)]]
    expr_level_permutation_temp_perm_i <- expr_level_permutation_temp_split[[as.character(perm_index)]]
    gene_by_decile <- expr_level_permutation_temp_perm_i %>% 
      filter(permutation_number != 0) %>% 
      split(.$decile) %>% 
      lapply(\(df) pull(df, Gene.refGene))
    for (i in 1:group_num){
      ## select mutation in genes belonging to each decile
      gene_name_temp <- gene_by_decile[[as.character(i)]]
      genic_permutation_decile_case_temp_perm_i <- genic_permutation_temp_perm_i %>% 
        filter(Gene.refGene %in% gene_name_temp) %>% arrange(Start, Chr)
      
      key_i <- paste0(Case_ID_temp, "_", i)
      if (nrow(genic_permutation_decile_case_temp_perm_i) == 0) {
        permutation_decile_grange_list_temp_perm_i[[key_i]] <- GRanges()
      } else {
        gr <- GRanges(
          seqnames = paste0("chr", genic_permutation_decile_case_temp_perm_i$Chr),
          ranges = IRanges(
            start = genic_permutation_decile_case_temp_perm_i$Start,
            end = genic_permutation_decile_case_temp_perm_i$End
          )
        )
        mcols(gr)$Ref <- genic_permutation_decile_case_temp_perm_i$Ref
        mcols(gr)$Alt <- genic_permutation_decile_case_temp_perm_i$Alt
        permutation_decile_grange_list_temp_perm_i[[key_i]] <- gr
      }
    }
    
    chr_length <- seqlengths(Hsapiens)
    seqlengths(permutation_decile_grange_list_temp_perm_i) <- chr_length[names(seqlengths(permutation_decile_grange_list_temp_perm_i))]
    seqlevels(permutation_decile_grange_list_temp_perm_i) <- seqlevels(permutation_decile_grange_list_temp_perm_i)[order(factor(seqlevels(permutation_decile_grange_list_temp_perm_i), levels = chr_orders))]
    genome(permutation_decile_grange_list_temp_perm_i) = 'hg19'
    
    ##### Calculate 96 trinuclutide context
    # cat("Calculate permutation 96 trinuclutide context for Case_ID:", Case_ID_temp, "...\n")
    permutation_mut_mat_temp_perm_i <- mut_matrix(permutation_decile_grange_list_temp_perm_i, ref_genome = ref_genome) %>% 
      t() %>% as.data.frame() %>% 
      mutate(perm.id = perm_index, 
             Condition = condition_temp, 
             decile = sub(".*_", "", rownames(.)), 
             Case_ID = sub("_.*", "", rownames(.)))

    permutation_mut_mat_list[[perm_index]] <- permutation_mut_mat_temp_perm_i

    if(perm_index %% 100 == 0){
      cat("Condition:", condition_temp, "Case_ID:", Case_ID_temp, "Permutation index:", perm_index, "\n")
    }
  }
  
  permutation_mut_mat_temp <- do.call(rbind, permutation_mut_mat_list)

  analysis_ID <- paste0(permutation_num, "_perms_", group_num, "G_mut_mat_")
  fig_save_dir <- paste0("./results/1-scRNAseq_analysis_by_case/", analysis_ID, "/", Case_ID_temp, "/")
  suppressWarnings(dir.create(fig_save_dir, recursive = TRUE))
  
  write.csv(mut_num_genic_merged_temp, paste0(fig_save_dir, "mut_num_genic_merged_by_case_decile_", Case_ID_temp, ".csv"))
  saveRDS(mutation_mut_mat_temp, paste0(fig_save_dir, "/mutation_mut_mat_by_case_decile_", Case_ID_temp, ".rds"))
  saveRDS(permutation_mut_mat_temp, paste0(fig_save_dir, "/permutation_mut_mat_by_case_decile_", Case_ID_temp, ".rds"))
  write.csv(mutation_mut_mat_temp[, 1:96], paste0(fig_save_dir, "mutation_mut_mat_by_case_decile_", Case_ID_temp, ".csv"))
  write.csv(permutation_mut_mat_temp[, 1:96], paste0(fig_save_dir, "permutation_mut_mat_by_case_decile_", Case_ID_temp, ".csv"))
}
