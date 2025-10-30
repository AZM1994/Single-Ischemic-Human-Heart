library(readr)
library(dplyr)

setwd("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration")

read_and_filter <- function(Donor_id, Condition) {
  path <- paste0("Mutation_Enrichment_Analysis/extract_germline_GATK/4-ANNOVAR_vcf_per_donor/", Condition, "/", Donor_id, ".snvs.pass.hg19_multianno.csv")
  
  selected_cols <- cols_only(Chr = col_character(), Start = col_double(), End = col_double(), Ref = col_character(), Alt = col_character(), 
                             Func.refGene = col_character(), Gene.refGene = col_character(), ExonicFunc.refGene = col_character(), AAChange.refGene = col_character(), 
                             AF_popmax = col_character(), AF = col_character(),  Polyphen2_HVAR_pred = col_character(), SIFT_pred = col_character(), 
                             CADD_phred = col_character(), REVEL_score = col_character(), MetaSVM_pred = col_character(), MPC_score = col_character())
  
  df <- read_csv(path, col_types = selected_cols, show_col_types = FALSE) %>% 
    mutate(AF_popmax = as.numeric(na_if(AF_popmax, ".")), AF = as.numeric(na_if(AF, ".")), CADD_phred = as.numeric(na_if(CADD_phred, ".")), 
           REVEL_score = as.numeric(na_if(REVEL_score, ".")), MPC_score = as.numeric(na_if(MPC_score, ".")))
    # filter(is.na(AF_popmax) | AF_popmax < 0.001) %>% 
    # mutate(AF_status = case_when(is.na(AF_popmax) ~ "Novel",
    #                              AF_popmax < 0.0001 ~ "Extremely rare (<0.01%)",
    #                              AF_popmax < 0.001  ~ "Ultra-rare (<0.1%)",
    #                              AF_popmax < 0.01   ~ "Rare (<1%)",
    #                              AF_popmax < 0.05   ~ "Low-freq (<5%)", TRUE ~ "Common"),
    #        donor = Donor_id, Condition = Condition)
  return(df)
}

donor_cell_map <- data.frame(
  Donor_id  = c("4402_liv_bulk", "1864-bulk-kid", "6032LivBulk", "4638_1024-pfc-bulk", "1863-bulk-kidney", 
                "1039-gDNA", "1940-bulk-kidney", "5919-bulk-kidney", "5828_liv_bulk", "5657-bulk-kid", 
                "1363_Brain_Bulk", "1673_liv_bulk", "604_liver_bulk", "1743_Brain_bulk", "1113-bulk-kidney"), 
  Condition = c("Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "IHD", "IHD", "IHD", "IHD", "IHD"), 
  stringsAsFactors = FALSE) %>% 
  mutate(Condition = factor(Condition, levels = c("Control", "IHD"))) %>% arrange(Condition)

# for (i in 1:nrow(donor_cell_map)) {
i=6
  Donor_id <- donor_cell_map$Donor_id[i]
  Condition <- donor_cell_map$Condition[i]
  message("Processing ", Condition, " donor: ", Donor_id)
  
  donor_variants <- read_and_filter(Donor_id, Condition)
  
  # out_path <- paste0("Mutation_Enrichment_Analysis/extract_germline_GATK/5-Germline_by_donor/", Donor_id, "_germline_filtered.csv")
  # write.csv(donor_variants, out_path, row.names = FALSE)
# }

  library(vcfR)
  
  # read VCF file
  vcf <- read.vcfR("Mutation_Enrichment_Analysis/extract_germline_GATK/4-ANNOVAR_vcf_per_donor/Control/5919_F2_germline_02.vcf")
  scan2_germline <- as.data.frame(vcf@fix) %>%
    dplyr::select(CHROM, POS, REF, ALT) %>%
    dplyr::mutate(
      CHROM = as.character(CHROM),
      POS   = as.numeric(POS),
      REF   = as.character(REF),
      ALT   = as.character(ALT)
    ) %>%
    dplyr::distinct()
  
  gatk_germline <- donor_variants %>%
    dplyr::select(Chr, Start, Ref, Alt, Gene.refGene, Func.refGene) %>%
    dplyr::mutate(
      Chr  = as.character(Chr),
      Start = as.numeric(Start),
      Ref  = as.character(Ref),
      Alt  = as.character(Alt)
    ) %>%
    dplyr::distinct()
  
  intersect_germline <- dplyr::inner_join(
    scan2_germline,
    gatk_germline,
    by = c("CHROM" = "Chr", "POS" = "Start", "REF" = "Ref", "ALT" = "Alt")
  )
  
  n_scan2  <- nrow(scan2_germline)
  n_gatk   <- nrow(gatk_germline)
  n_common <- nrow(intersect_germline)
  
  cat("SCAN2 germline:", n_scan2, "\n")
  cat("GATK germline:", n_gatk, "\n")
  cat("Overlap:", n_common, " (", round(100 * n_common / min(n_scan2, n_gatk), 2), "% of smaller set)\n")
  