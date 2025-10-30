library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome)
library(GenomeInfoDb)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(data.table)
library(ggpubr)

setwd("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration/Mutation_Enrichment_Analysis")
replicates_num <- "Rep1" # Rap1, Rep2
results_dir <- paste0("Hi_C_Enrichment/", replicates_num)
color_set <- c(colorRampPalette(c("skyblue","dodgerblue4"))(9)[7], colorRampPalette(c("pink","firebrick"))(4)[3])
group_num <- 8

##### Load Hi-C PC1 data
Hi_C_data_path <- if (replicates_num == "Rep1") {
  "data/GSM2845454_CM_Rep1_500KB_PC1.bedGraph"
} else if (replicates_num == "Rep2") {
  "data/GSM2845455_CM_Rep2_500KB_PC1.bedGraph"
} else {
  stop("Invalid replicate number. Use 'Rep1' or 'Rep2'.")
}

pc1_df <- read_table(Hi_C_data_path, col_names = c("Chrom", "Start", "End", "PC1"), skip = 1, col_types = cols(.default = col_character())) %>% 
  mutate(Start = as.integer(Start), End = as.integer(End), PC1 = as.numeric(PC1)) %>% 
  filter(!is.na(PC1))

##### Convert to GRanges
pc1_gr <- GRanges(seqnames = pc1_df$Chrom,
                  ranges = IRanges(start = pc1_df$Start + 1, end = pc1_df$End), 
                  PC1 = pc1_df$PC1)
seqlevelsStyle(pc1_gr) <- "UCSC"

##### read in metadata, SCAN2 mutation grl
Hypoxia_PTA_all_SCAN2 <- readRDS("data/1-sSNV_SCAN2_df_filtered.rds") %>% dplyr::select(c("Cell_ID", "Age", "Gender", "Case_ID", "Condition", "snv.burden", "snv.rate.per.gb"))
Hypoxia_PTA_all_SCAN2_collapsed <- Hypoxia_PTA_all_SCAN2 %>% distinct(Case_ID, .keep_all = TRUE)
Cell_ID_list <- Hypoxia_PTA_all_SCAN2$Cell_ID
Case_ID_list <- Hypoxia_PTA_all_SCAN2_collapsed$Case_ID
Condition_list <- unique(Hypoxia_PTA_all_SCAN2$Condition)

##### Annotate mutations with PC1 value, flatten to single GRanges and annotate with sample name
ssnv_grl <- readRDS("data/1-ssnv_grl.rds")
ssnv_gr <- unlist(ssnv_grl, use.names = FALSE)
ssnv_gr$Cell_ID <- rep(names(ssnv_grl), elementNROWS(ssnv_grl))  # preserves sample info
seqlevelsStyle(ssnv_gr) <- "UCSC"

overlaps <- findOverlaps(ssnv_gr, pc1_gr)
ssnv_gr$PC1 <- NA_real_
ssnv_gr$PC1[queryHits(overlaps)] <- pc1_gr$PC1[subjectHits(overlaps)]

##### add groups, Case_ID, Condiditon
# ssnv_gr_df <- as.data.frame(mcols(ssnv_gr)) %>% 
#   filter(!is.na(PC1)) %>% 
#   mutate(PC1_decile = ntile(PC1, group_num), 
#          PC1_decile_label = case_when(PC1_decile <= 3 ~ "B-compartment (inactive)", PC1_decile >= 6 ~ "A-compartment (active)", TRUE ~ "Intermediate"), 
#          Compartment = ifelse(PC1 < 0, "B", "A")) %>% 
#   right_join(Hypoxia_PTA_all_SCAN2 %>% dplyr::select(Cell_ID, Case_ID, Condition), by = c("Cell_ID"))

# saveRDS(ssnv_gr_df, paste0(results_dir, "/ssnv_gr_df.rds"))
ssnv_gr_df <- readRDS(paste0(results_dir, "/ssnv_gr_df.rds"))

pdf(paste0(results_dir, "/PC1.pdf"), width = 5, height = 3)
  ggplot(ssnv_gr_df, aes(x = PC1)) + geom_histogram(aes(y = ..density..), bins = 100, fill = "steelblue", alpha = 0.7) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") + facet_wrap(. ~ Condition) + 
    labs(title = "Distribution of PC1 Values", x = "PC1 (Hi-C Eigenvector)", y = "Relative Frequency")
  
  ggplot(ssnv_gr_df, aes(x = Compartment), fill = Condition) + geom_bar() + facet_wrap(. ~ Condition) + 
    labs(title = "Mutation Counts by Hi-C Compartment", x = "Compartment", y = "Mutation Count")
dev.off()

##### read in all permutation results and annotate mutations with PC1 value
total_permutation <- 10000
batch_num <- 1
permutations_per_batch <- 2000
permutation_start <- (batch_num - 1) * permutations_per_batch + 1
permutation_end <- batch_num * permutations_per_batch
cat(permutation_start, permutation_end, "\n")

all_perm_cases <- list()
for (Case_ID_temp in Case_ID_list) {
  permutation_temp <- data.table()
  Cell_ID_list <- Hypoxia_PTA_all_SCAN2 %>% filter(Case_ID == Case_ID_temp) %>% pull(Cell_ID)

  cat("Case:", Case_ID_temp, Cell_ID_list, "\n")

  for (Cell_ID_temp in Cell_ID_list) {
    permutation_path <- paste0("data/permutation_SCAN2/sSNV/", Cell_ID_temp, "/perms_snv.hg19_multianno.csv")
    case_temp <- Case_ID_temp
    condition_temp <- Hypoxia_PTA_all_SCAN2$Condition[Hypoxia_PTA_all_SCAN2$Cell_ID == Cell_ID_temp]
    cat("Cell:", Cell_ID_temp, "Case:", case_temp, "Condition:", condition_temp, "\n")

    total_lines <- as.integer(system(paste("awk 'END{print NR}'", shQuote(permutation_path)), intern = TRUE))
    lines_to_read <- floor((total_lines - 1) * permutations_per_batch / total_permutation)
    skip_lines <- floor((total_lines - 1) * (permutation_start - 1) / total_permutation)
    cat(lines_to_read, skip_lines, "\n")

    header <- fread(permutation_path, nrows = 0)
    permutation_cell <- fread(permutation_path, nrows = lines_to_read, skip = skip_lines, col.names = colnames(header)) %>%
      as.data.table() %>%
      .[, perm.id := rep(permutation_start:permutation_end, each = .N / permutations_per_batch)] %>%
      .[, `:=`(
        Gene.refGene = {
          gene_lists <- strsplit(Gene.refGene, ";")
          dist_lists <- str_extract_all(GeneDetail.refGene, "\\d+") %>% lapply(as.integer)
          mapply(function(gs, ds) {
            if ("NONE" %in% gs && length(gs) > 1) gs <- gs[gs != "NONE"]
            if (length(ds) == length(gs) && length(ds) > 0) gs[which.min(ds)] else gs[1]
          }, gene_lists, dist_lists)
        },
        GeneDetail.refGene = {
          dist_lists <- str_extract_all(GeneDetail.refGene, "\\d+") %>% lapply(as.integer)
          sapply(dist_lists, function(ds) if (length(ds)) paste0("dist=", min(ds)) else NA)
        })] %>%
      dplyr::select(Chr, Start, End, Ref, Alt, Func.refGene, Gene.refGene, perm.id) %>%
      mutate(Cell_ID = Cell_ID_temp, Case_ID = case_temp, Condition = condition_temp)

    permutation_temp <- rbind(permutation_temp, permutation_cell)
  }

  ### Annotate PC1 for this case
  perm_gr <- GRanges(seqnames = permutation_temp$Chr,
                     ranges = IRanges(start = permutation_temp$Start, end = permutation_temp$End))
  seqlevelsStyle(perm_gr) <- "UCSC"
  seqlevelsStyle(pc1_gr) <- "UCSC"
  perm_hits <- findOverlaps(perm_gr, pc1_gr)
  permutation_temp$PC1 <- NA_real_
  permutation_temp$PC1[queryHits(perm_hits)] <- pc1_gr$PC1[subjectHits(perm_hits)]
  permutation_temp <- permutation_temp %>%
    filter(!is.na(PC1)) %>%
    mutate(PC1_decile = ntile(PC1, group_num),
           PC1_decile_label = case_when(PC1_decile <= 3 ~ "B-compartment (inactive)", PC1_decile >= 6 ~ "A-compartment (active)", TRUE ~ "Intermediate"),
           Compartment = ifelse(PC1 < 0, "B", "A")) %>%
    filter(!is.na(PC1_decile), !is.na(Condition)) %>%
    dplyr::select(c("perm.id", "Case_ID", "PC1_decile", "Condition")) %>%
    group_by(perm.id, PC1_decile, Condition, Case_ID) %>%
    summarize(permut_sum = n(), .groups = "drop")

  all_perm_cases[[as.character(Case_ID_temp)]] <- permutation_temp
  # saveRDS(all_perm_cases, paste0(results_dir, "/all_perm_cases", as.character(Case_ID_temp), ".rds"))
  saveRDS(all_perm_cases, paste0(results_dir, "/batch/all_perm_cases_", as.character(Case_ID_temp), "_batch", batch_num, ".rds"))
}

##### Combine all batches for 1363, 1673, 604, 1743, 1113
# case_id <- "1113"
# batch_files <- list.files(paste0(results_dir, "/batch"), pattern = paste0("all_perm_cases_", case_id, "_batch\\d+\\.rds$"), full.names = TRUE)
# all_batches_case <- lapply(batch_files, readRDS)
# combined_case <- bind_rows(all_batches_case)
# saveRDS(combined_case, paste0(results_dir, "/all_perm_cases", case_id, ".rds"))

##### read all permutation results and perform enrichment analysis
all_perm_cases <- list()
perm_rds_files <- list.files(results_dir, pattern = "^all_perm_cases.*\\.rds$", full.names = TRUE)
for (file in perm_rds_files) {
  Case_ID_temp <- stringr::str_extract(basename(file), "(?<=all_perm_cases).*(?=\\.rds)")
  all_perm_cases[[Case_ID_temp]] <- readRDS(file)
}

all_perm_cases_flat <- lapply(all_perm_cases, function(x) {
  if (is.list(x) && length(x) == 1 && is.data.frame(x[[1]])) {
    return(x[[1]])
  } else {
    return(x)
  }
})

mutation_summary_plot <- ssnv_gr_df %>% filter(!is.na(PC1_decile), !is.na(Condition)) %>% 
  dplyr::select(c("Case_ID", "PC1_decile", "Condition")) %>% 
  group_by(PC1_decile, Condition, Case_ID) %>% 
  summarize(mut_num = n(), .groups = "drop")

permutation_summary_plot <- bind_rows(all_perm_cases_flat) 

merged_summary_plot <- merge(mutation_summary_plot, permutation_summary_plot) %>% 
  merge(Hypoxia_PTA_all_SCAN2_collapsed[, c("Case_ID", "Age")]) %>% 
  mutate(enrichment_ratio = mut_num / permut_sum) %>% 
  filter(!is.na(enrichment_ratio) & is.finite(enrichment_ratio)) %>% 
  filter(enrichment_ratio != 0) %>% 
  group_by(Condition, PC1_decile) %>% 
  filter(enrichment_ratio >= quantile(enrichment_ratio, 0.05) & enrichment_ratio <= quantile(enrichment_ratio, 0.95)) %>%
  summarise(mean_ER = mean(enrichment_ratio, na.rm = TRUE, .groups = "drop"), sd_ER = sd(enrichment_ratio, na.rm = TRUE)) %>% 
  mutate(Condition = factor(Condition, level = c("Control", "IHD")), 
         PC1_decile = factor(PC1_decile, level = seq(1 : group_num)))

overall_p <- wilcox.test(mean_ER ~ Condition, data = merged_summary_plot, alternative = c("two.sided"))$p.value
overall_star <- case_when(overall_p <= 0.001 ~ "***", overall_p <= 0.01  ~ "**", overall_p <= 0.05  ~ "*", TRUE ~ "ns")
overall_label <- paste0("Wilcoxon test Control v.s. IHD, P = ", signif(overall_p, 2))

p_SNV_HiC_enrichment <- ggplot(merged_summary_plot, aes(x = PC1_decile, y = mean_ER, group = Condition, color = Condition)) + 
  geom_hline(yintercept = 1, color = "black", linewidth = 0.6) + geom_line(position = position_dodge(width = 0.1), linewidth = 0.7) + 
  geom_point(position = position_dodge(width = 0.1), size = 2) + stat_cor(size = 6, show.legend = FALSE, label.x.npc = "right", hjust = 1) +
  geom_errorbar(aes(ymin = mean_ER - sd_ER, ymax = mean_ER + sd_ER), width = 0.2, position = position_dodge(width = 0.1)) +
  geom_smooth(data = merged_summary_plot, aes(x = PC1_decile, y = mean_ER, color = Condition, fill = Condition, group = Condition),
              linetype = "22", method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.7) + 
  annotate(geom = "polygon", x = c(1, group_num, group_num), y = c(0.8, 0.8, 0.8 + 0.04), fill = "#36C9CB") + 
  scale_color_manual(values = c("Control" = color_set[1], "IHD" = color_set[2]), guide = "legend") + 
  scale_fill_manual(values = c("Control" = color_set[1], "IHD" = color_set[2]), guide = "legend") + theme_linedraw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24)) + 
  scale_y_continuous(breaks = seq(0.8, 1.4, by = 0.1), limits = c(0.8, 1.4)) +
  labs(x = "Chromatin activity level / replication timing", y = "sSNV enrichment ratio \n (obs/exp)", color = "Condition", subtitle = overall_label)
ggsave(paste0(results_dir, "/sSNV_HiC_Enrichment_", replicates_num, ".pdf"), plot = p_SNV_HiC_enrichment, width = 9, height = 5, dpi = 600)
