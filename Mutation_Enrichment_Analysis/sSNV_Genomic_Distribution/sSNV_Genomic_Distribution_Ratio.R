###################################################################################################
##### 1. Read annotation of genomic location and functional categories                        #####
##### 2. Calculate the exonic-to-intronic ratio and nonsynonymous-to-synonymous ratio (dN/dS) #####
##### 3. Normalize the ratio by germline mutation (from SCAN2 resampled training variants)    #####
##### 4. Use LME modeling and emmeans test between conditions, and test the ratio vs 1        #####
###################################################################################################

library(ggplot2)
library(reshape2)
library(stringr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggbreak)
library(MuMIn)

setwd("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration/Mutation_Enrichment_Analysis")
exon_intron_dir <- "sSNV_Genomic_Distribution/1-exon_vs_intron"
dir.create(exon_intron_dir, recursive = TRUE, showWarnings = FALSE)
nonsyn_syn_dir <- "sSNV_Genomic_Distribution/2-nonsyn_vs_syn"
dir.create(nonsyn_syn_dir, recursive = TRUE, showWarnings = FALSE)
Ctrl_IHD_color <- c("dodgerblue3", "firebrick3")

genomic_factor <- 5.845
ctrl_name <- "Control"
dis_name <- "IHD"

##### genomic_context for Control, IHD, and germline
SCAN2_df <- readRDS("data/1-sSNV_SCAN2_df_filtered.rds") %>% 
  select(Cell_ID, Case_ID, snv.rate.per.gb, Condition, Age, Gender, Color, Outline_Color)

heart_PTA_all_Control_vcf <- read.table("data/annotation_results/all_age_sSNV/Control/heart_PTA_all.all_age.Control_ssnv.vcf", sep = "\t")
heart_PTA_all_IHD_vcf <- read.table("data/annotation_results/all_age_sSNV/IHD/heart_PTA_all.all_age.IHD_ssnv.vcf", sep = "\t")

genomic_context_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Cell_ID", "Func.refGene", "Gene.refGene")
genomic_context_ctrl <- read.csv("data/annotation_results/all_age_sSNV/Control/heart_PTA_all.all_age.Control_ssnv.hg19_multianno.csv", header = TRUE) %>% 
  mutate(Gene.refGene.split = strsplit(as.character(Gene.refGene), ";"), 
         dist_vals = str_extract_all(GeneDetail.refGene, "(?<=dist=)[0-9]+"), 
         dist_numeric = lapply(dist_vals, as.numeric)) %>% 
  rowwise() %>% 
  mutate(closest_idx = ifelse(length(dist_numeric) > 0, which.min(dist_numeric), NA_integer_), 
         Gene.refGene = ifelse(!is.na(closest_idx), Gene.refGene.split[[closest_idx]], Gene.refGene), 
         ClosestDist = ifelse(!is.na(closest_idx), dist_numeric[closest_idx], NA), 
         Func.refGene = case_when(Func.refGene == "upstream;downstream" & closest_idx == 1 ~ "upstream", 
                                  Func.refGene == "upstream;downstream" & closest_idx == 2 ~ "downstream", TRUE ~ Func.refGene), 
         Func.refGene = gsub("exonic;splicing", "exonic", Func.refGene), Func.refGene = gsub("UTR5;UTR3", "UTR5", Func.refGene)) %>% 
  ungroup() %>% 
  mutate(Cell_ID = heart_PTA_all_Control_vcf$V8) %>% dplyr::select(all_of(genomic_context_colnames))

genomic_context_dis <- read.csv("data/annotation_results/all_age_sSNV/IHD/heart_PTA_all.all_age.IHD_ssnv.hg19_multianno.csv", header = TRUE) %>% 
  mutate(Gene.refGene.split = strsplit(as.character(Gene.refGene), ";"), 
         dist_vals = str_extract_all(GeneDetail.refGene, "(?<=dist=)[0-9]+"), 
         dist_numeric = lapply(dist_vals, as.numeric)) %>% 
  rowwise() %>% 
  mutate(closest_idx = ifelse(length(dist_numeric) > 0, which.min(dist_numeric), NA_integer_), 
         Gene.refGene = ifelse(!is.na(closest_idx), Gene.refGene.split[[closest_idx]], Gene.refGene), 
         ClosestDist = ifelse(!is.na(closest_idx), dist_numeric[closest_idx], NA), 
         Func.refGene = case_when(Func.refGene == "upstream;downstream" & closest_idx == 1 ~ "upstream", 
                                  Func.refGene == "upstream;downstream" & closest_idx == 2 ~ "downstream", TRUE ~ Func.refGene), 
         Func.refGene = gsub("exonic;splicing", "exonic", Func.refGene), Func.refGene = gsub("UTR5;UTR3", "UTR5", Func.refGene)) %>% 
  ungroup() %>% 
  mutate(Cell_ID = heart_PTA_all_IHD_vcf$V8) %>% dplyr::select(all_of(genomic_context_colnames))

genomic_context <- rbind(genomic_context_ctrl, genomic_context_dis)
genomic_SCAN2_df <- genomic_context %>% 
  inner_join(SCAN2_df %>% select(Cell_ID, Case_ID, Age, Gender, Condition), by = "Cell_ID")

##### summarize all genetic regions
Func.refGene_oldnames <- c("intergenic", "upstream", "UTR5", "exonic", "UTR3", "downstream", "splicing", "intronic")
Func.refGene_newnames <- c("intergenic", "upstream", "5' UTR", "exonic", "3' UTR", "downstream", "splice site", "intronic")
Func.refGene_Group_df <- table(factor(genomic_SCAN2_df$Func.refGene, levels = Func.refGene_oldnames), genomic_SCAN2_df$Cell_ID) %>% `row.names<-`(Func.refGene_newnames)

## read germline exo_intro_ratio_df
germline_exo_intro_ratio_df <- read.csv("extract_germline_SCAN2/resampled_training_variants/exo_intro_ratio.csv") %>% rename(germline_ratio = exo_intro_germline_ratio)

Func.refGene_Group_df_summary <- data.frame(unclass(t(Func.refGene_Group_df[rownames(Func.refGene_Group_df) %in% c("exonic","intronic"), ])))
Func.refGene_Group_df_summary <- Func.refGene_Group_df_summary %>% 
  mutate(Cell_ID = factor(rownames(Func.refGene_Group_df_summary), levels = c(SCAN2_df$Cell_ID))) %>% 
  mutate(scRatio = exonic / intronic) %>% 
  inner_join(germline_exo_intro_ratio_df, by = "Cell_ID") %>% 
  mutate(normalized_ratio = scRatio / germline_ratio)

df_ratio_Func <- Func.refGene_Group_df_summary %>% inner_join(SCAN2_df, by = "Cell_ID") %>% 
  mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% 
  mutate(Condition = relevel(Condition, ref = ctrl_name)) %>% 
  group_by(Condition) %>% arrange(Age, .by_group = TRUE)

selected_func_list <- "normalized_ratio"
sig_digits <- 2
pval_df_Func <- data.frame(Func = character(), estimate = numeric(), SE = numeric(), 
                           df = numeric(), t = numeric(), pval = numeric(), stringsAsFactors = FALSE)

##### plot for exon_vs_intron ratio
for (i in c(1 : length(selected_func_list))){
  cat("index", i, ":", selected_func_list[i])
  ## burden_df: the input matrix with columns: burden, Condition, Case_ID, Age, Color
  burden_df <- df_ratio_Func[c(selected_func_list[i], "Age", "Condition", "Case_ID", "Color", "Outline_Color")] %>% 
    setNames(c("snv.rate.per.gb", "Age", "Condition", "Case_ID", "Color", "Outline_Color")) %>% 
    # filter(snv.rate.per.gb < 10) %>% 
    mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% 
    mutate(Condition = relevel(Condition, ref = ctrl_name)) %>% 
    group_by(Condition) %>% arrange(Age, .by_group = TRUE) %>% ungroup()
  burden_df_ctrl <- burden_df %>% filter(Condition == "Control")
  burden_df_dis <- burden_df %>% filter(Condition == "IHD")
  
  generate_line_data <- function(df_ctrl, coef_vec) {
    x <- range(df_ctrl$Age)
    y <- x * coef_vec[2] + coef_vec[1]
    data.frame(ctrl_fitting_x = x, ctrl_fitting_y = y)
  }
  
  generate_y_breaks <- function(y_data) {
    y_min <- min(y_data, na.rm = TRUE)
    y_max <- max(y_data, na.rm = TRUE)
    range_span <- y_max - y_min
    step_size <- signif(range_span / 4, 1)
    breaks <- seq(floor(y_min / step_size) * step_size, ceiling(y_max / step_size) * step_size, by = step_size)
    return(breaks)
  }
  
  ## Check sums and choose model
  check_sums <- sum(burden_df_ctrl$snv.rate.per.gb) > 0 & sum(burden_df_dis$snv.rate.per.gb) > 0
  burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + (1|Case_ID), burden_df, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
  
  if (check_sums) {
    cat(", Normal\n")
    burden_age_model_ctrl <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_ctrl, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
    coef_ctrl <- fixef(burden_age_model_ctrl)
  } else {
    cat(", Low number\n")
    burden_age_model_ctrl <- lm(snv.rate.per.gb ~ Age, burden_df_ctrl)
    coef_ctrl <- coef(burden_age_model_ctrl)
  }
  
  aging_rate <- format(round(coef_ctrl[2], digits = sig_digits), nsmall = sig_digits)
  
  ## check p value and R^2
  anova_pval <- anova(burden_age_model)$"Pr(>F)"[2]
  anova_pval_print <- formatC(signif(anova_pval, digits = sig_digits), digits = sig_digits, format="fg", flag="#")
  anova_pval_ctrl <- anova(burden_age_model_ctrl)$"Pr(>F)"[1]
  anova_pval_print_ctrl <- formatC(signif(anova_pval_ctrl, digits = sig_digits), digits = sig_digits, format="e", flag="#")
  
  aov_tab <- anova(burden_age_model)
  if ("Condition" %in% rownames(aov_tab)) {
    anova_pval <- aov_tab["Condition", "Pr(>F)"]
  } else if (any(grepl("^Condition", rownames(aov_tab)))) {
    anova_pval <- aov_tab[grep("^Condition", rownames(aov_tab))[1], "Pr(>F)"]
  } else {
    anova_pval <- NA_real_
  }
  ct <- coef(summary(burden_age_model))
  cond_row <- ct[grep("^Condition", rownames(ct)), , drop = FALSE]
  est  <- if (nrow(cond_row)) cond_row[1, "Estimate"] else NA_real_
  SE   <- if (nrow(cond_row)) cond_row[1, "Std. Error"] else NA_real_
  tval <- if (nrow(cond_row)) cond_row[1, "t value"] else NA_real_
  ddf  <- if (nrow(cond_row) && "df" %in% colnames(cond_row)) cond_row[1, "df"] else NA_real_
  
  ci <- suppressMessages(confint(burden_age_model, method = "Wald"))
  cond_row <- rownames(ci)[grepl("^Condition", rownames(ci))]
  if (length(cond_row)) {
    ci_vals <- ci[cond_row[1], ]
    ci_low  <- ci_vals[1]; ci_high <- ci_vals[2]
  } else { ci_low <- NA; ci_high <- NA }
  ci_low  <- as.numeric(ci_vals[1])
  ci_high <- as.numeric(ci_vals[2])
  
  pval_df_Func <- rbind(pval_df_Func, data.frame(Func = selected_func_list[i], estimate = est, SE = SE, 
                                                 ci_low=ci_low, ci_high=ci_high, pval = anova_pval))
  
  geom_line_data <- generate_line_data(burden_df_ctrl, coef_ctrl)
  legend_data <- burden_df[c(1, nrow(burden_df)), c("Age", "snv.rate.per.gb", "Condition")]
  
  yticks_plot <- generate_y_breaks(burden_df$snv.rate.per.gb)
  ylim_plot <- range(yticks_plot, na.rm = TRUE)
  p_SNV_burden_lme <- ggplot(burden_df, aes(x = Age, y = snv.rate.per.gb)) + 
    geom_point(pch = 21, data = legend_data, aes(x = Age, y = snv.rate.per.gb, color = Condition, fill = Condition), size = 5) + 
    geom_point(pch = 21, color = burden_df$Outline_Color, fill = burden_df$Color, size = 5) + 
    geom_line(aes(x = ctrl_fitting_x, y = ctrl_fitting_y), color = "dodgerblue2", linetype = "22", data = geom_line_data) + 
    geom_hline(yintercept = 0, color = "gray40", linetype = "22") + 
    # annotate("text", size = 6, x = 0, y = 0.95 * ylim_plot[2], 
    #          label = paste("Aging effect:", aging_rate, "sSNVs/(GB·year),", "P =", anova_pval_print_ctrl), hjust = 0) +
    annotate("text", size = 6, x = 0, y = 0.87 * ylim_plot[2], 
             label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB,", "P =", anova_pval_print), hjust = 0) +
    scale_color_manual(values = c("IHD" = "black", "Control" = "black"), guide = "legend") + 
    scale_fill_manual(values = c("Control" = Ctrl_IHD_color[1], "IHD" = Ctrl_IHD_color[2]), guide = "legend") + 
    scale_y_continuous(limits = c(min(min(geom_line_data[, 2]), 0), ylim_plot[2]), breaks = yticks_plot, labels = yticks_plot) + 
    labs(x = "Age (years)", y = "Normalized exonic-intronic ratio") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
  ggsave(paste0(exon_intron_dir, "/ratio/exon_vs_intron.ratio.pdf"), plot = p_SNV_burden_lme, width = 9, height = 6, dpi = 600)
}

################################################################################
##################### nonsynonymous vs synonymous analysis #####################
################################################################################

##### exonic_context for Control, IHD, and germline
exonic_context_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Cell_ID", "ExonicFunc.refGene", "AAChange.refGene")
exonic_context_ctrl <- read.csv("data/annotation_results/all_age_sSNV/Control/heart_PTA_all.all_age.Control_ssnv.hg19_multianno.csv", header = TRUE) %>% 
  mutate(Cell_ID = heart_PTA_all_Control_vcf$V8) %>% filter(ExonicFunc.refGene != ".") |> dplyr::select(all_of(exonic_context_colnames))
exonic_context_dis <- read.csv("data/annotation_results/all_age_sSNV/IHD/heart_PTA_all.all_age.IHD_ssnv.hg19_multianno.csv", header = TRUE) %>% 
  mutate(Cell_ID = heart_PTA_all_IHD_vcf$V8) %>% filter(ExonicFunc.refGene != ".") |> dplyr::select(all_of(exonic_context_colnames))
exonic_context_germline <- read.csv("data/GATK_PTA_Germline/1039-gDNA.snvs.pass.hg19_multianno_exonic_context.csv", header = TRUE) %>%
  mutate(Cell_ID = "germline_bulk") %>% filter(ExonicFunc.refGene != ".") |> dplyr::select(all_of(exonic_context_colnames))

exonic_context <- rbind(exonic_context_ctrl, exonic_context_dis, exonic_context_germline)
exonic_SCAN2_df <- SCAN2_df %>% select(Cell_ID, Case_ID, Age, Gender, Condition) %>% left_join(exonic_context, by = "Cell_ID")

################################################################################
##### summarize all mutation types
ExonicFunc.refGene_oldnames <- c("nonsynonymous SNV", "synonymous SNV", "stopgain", "stoploss")
ExonicFunc.refGene_newnames <- c("nonsynonymous", "synonymous", "stopgain", "stoploss")
ExonicFunc.refGene_Group_df <- table(factor(exonic_SCAN2_df$ExonicFunc.refGene, levels = ExonicFunc.refGene_oldnames), exonic_SCAN2_df$Cell_ID) %>% 
  `row.names<-`(ExonicFunc.refGene_newnames)

## read germline dnds_ratio_df
germline_dnds_ratio_df <- read.csv("extract_germline_SCAN2/resampled_training_variants/dnds.csv") %>% rename(germline_ratio = germline_dNdS)

ExonicFunc.refGene_Group_df_summary <- data.frame(unclass(t(ExonicFunc.refGene_Group_df[rownames(ExonicFunc.refGene_Group_df) %in% c("nonsynonymous","synonymous"), ])))
ExonicFunc.refGene_Group_df_summary <- ExonicFunc.refGene_Group_df_summary %>% 
  mutate(Cell_ID = factor(rownames(ExonicFunc.refGene_Group_df_summary), levels = c(SCAN2_df$Cell_ID))) %>% 
  mutate(scRatio = nonsynonymous / synonymous) %>% 
  inner_join(germline_dnds_ratio_df, by = "Cell_ID") %>% 
  mutate(normalized_ratio = scRatio / germline_ratio)

df_ratio_ExonicFunc <- ExonicFunc.refGene_Group_df_summary %>% inner_join(SCAN2_df, by = "Cell_ID") %>% 
  mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% 
  mutate(Condition = relevel(Condition, ref = ctrl_name)) %>% 
  group_by(Condition) %>% arrange(Age, .by_group = TRUE) %>% 
  filter(!is.na(normalized_ratio), is.finite(normalized_ratio))

selected_ExonicFunc_list <- "normalized_ratio"
sig_digits <- 2
pval_df_ExonicFunc <- data.frame(ExonicFunc = character(), estimate = numeric(), SE = numeric(), 
                                 df = numeric(), t = numeric(), pval = numeric(), stringsAsFactors = FALSE)

for (i in c(1 : length(selected_func_list))){
  cat("index", i, ":", selected_func_list[i])
  ## burden_df: the input matrix with columns: burden, Condition, Case_ID, Age, Color
  burden_df <- df_ratio_ExonicFunc[c(selected_func_list[i], "Age", "Condition", "Case_ID", "Color", "Outline_Color")] %>% 
    setNames(c("snv.rate.per.gb", "Age", "Condition", "Case_ID", "Color", "Outline_Color")) %>% 
    mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% 
    mutate(Condition = relevel(Condition, ref = ctrl_name)) %>% 
    group_by(Condition) %>% arrange(Age, .by_group = TRUE) %>% ungroup()
  burden_df_ctrl <- burden_df %>% filter(Condition == "Control")
  burden_df_dis <- burden_df %>% filter(Condition == "IHD")
  
  generate_line_data <- function(df_ctrl, coef_vec) {
    x <- range(df_ctrl$Age)
    y <- x * coef_vec[2] + coef_vec[1]
    data.frame(ctrl_fitting_x = x, ctrl_fitting_y = y)
  }
  
  generate_y_breaks <- function(y_data) {
    y_min <- min(y_data, na.rm = TRUE)
    y_max <- max(y_data, na.rm = TRUE)
    range_span <- y_max - y_min
    step_size <- signif(range_span / 4, 1)
    breaks <- seq(floor(y_min / step_size) * step_size, ceiling(y_max / step_size) * step_size, by = step_size)
    return(breaks)
  }
  
  ## Check sums and choose model
  check_sums <- sum(burden_df_ctrl$snv.rate.per.gb) > 10 & sum(burden_df_dis$snv.rate.per.gb) > 10
  burden_age_model <- lmer(snv.rate.per.gb ~ Age + Condition + (1|Case_ID), burden_df, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
  
  if (check_sums) {
    cat(", Normal\n")
    burden_age_model_ctrl <- lmer(snv.rate.per.gb ~ Age + (1|Case_ID), burden_df_ctrl, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
    coef_ctrl <- fixef(burden_age_model_ctrl)
  } else {
    cat(", Low number\n")
    burden_age_model_ctrl <- lm(snv.rate.per.gb ~ Age, burden_df_ctrl)
    coef_ctrl <- coef(burden_age_model_ctrl)
  }
  
  aging_rate <- format(round(coef_ctrl[2], digits = sig_digits), nsmall = sig_digits)
  
  ## check p value and R^2
  anova_pval <- anova(burden_age_model)$"Pr(>F)"[2]
  anova_pval_print <- formatC(signif(anova_pval, digits = sig_digits), digits = sig_digits, format="fg", flag="#")
  anova_pval_ctrl <- anova(burden_age_model_ctrl)$"Pr(>F)"[1]
  anova_pval_print_ctrl <- formatC(signif(anova_pval_ctrl, digits = sig_digits), digits = sig_digits, format="e", flag="#")
  
  aov_tab <- anova(burden_age_model)
  if ("Condition" %in% rownames(aov_tab)) {
    anova_pval <- aov_tab["Condition", "Pr(>F)"]
  } else if (any(grepl("^Condition", rownames(aov_tab)))) {
    anova_pval <- aov_tab[grep("^Condition", rownames(aov_tab))[1], "Pr(>F)"]
  } else {
    anova_pval <- NA_real_
  }
  ct <- coef(summary(burden_age_model))
  cond_row <- ct[grep("^Condition", rownames(ct)), , drop = FALSE]
  est  <- if (nrow(cond_row)) cond_row[1, "Estimate"] else NA_real_
  SE   <- if (nrow(cond_row)) cond_row[1, "Std. Error"] else NA_real_
  tval <- if (nrow(cond_row)) cond_row[1, "t value"] else NA_real_
  ddf  <- if (nrow(cond_row) && "df" %in% colnames(cond_row)) cond_row[1, "df"] else NA_real_
  
  ci <- suppressMessages(confint(burden_age_model, method = "Wald"))
  cond_row <- rownames(ci)[grepl("^Condition", rownames(ci))]
  if (length(cond_row)) {
    ci_vals <- ci[cond_row[1], ]
    ci_low  <- ci_vals[1]; ci_high <- ci_vals[2]
  } else { ci_low <- NA; ci_high <- NA }
  ci_low  <- as.numeric(ci_vals[1])
  ci_high <- as.numeric(ci_vals[2])
  
  pval_df_ExonicFunc <- rbind(pval_df_ExonicFunc, data.frame(ExonicFunc = selected_func_list[i], estimate = est, SE = SE, 
                                                             ci_low=ci_low, ci_high=ci_high, pval = anova_pval))
  
  geom_line_data <- generate_line_data(burden_df_ctrl, coef_ctrl)
  legend_data <- burden_df[c(1, nrow(burden_df)), c("Age", "snv.rate.per.gb", "Condition")]
  
  yticks_plot <- generate_y_breaks(burden_df$snv.rate.per.gb)
  ylim_plot <- range(yticks_plot, na.rm = TRUE)
  p_SNV_burden_lme <- ggplot(burden_df, aes(x = Age, y = snv.rate.per.gb)) + 
    geom_point(pch = 21, data = legend_data, aes(x = Age, y = snv.rate.per.gb, color = Condition, fill = Condition), size = 5) + 
    geom_point(pch = 21, color = burden_df$Outline_Color, fill = burden_df$Color, size = 5) + 
    geom_line(aes(x = ctrl_fitting_x, y = ctrl_fitting_y), color = "dodgerblue2", linetype = "22", data = geom_line_data) + 
    geom_hline(yintercept = 0, color = "gray40", linetype = "22") + 
    # annotate("text", size = 6, x = 0, y = 0.95 * ylim_plot[2], 
    #          label = paste("Aging effect:", aging_rate, "sSNVs/(GB·year),", "P =", anova_pval_print_ctrl), hjust = 0) +
    annotate("text", size = 6, x = 0, y = 0.87 * ylim_plot[2], 
             label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB,", "P =", anova_pval_print), hjust = 0) +
    scale_color_manual(values = c("IHD" = "black", "Control" = "black"), guide = "legend") + 
    scale_fill_manual(values = c("Control" = Ctrl_IHD_color[1], "IHD" = Ctrl_IHD_color[2]), guide = "legend") + 
    scale_y_continuous(limits = c(min(min(geom_line_data[, 2]), 0), ylim_plot[2]), breaks = yticks_plot, labels = yticks_plot) + 
    labs(x = "Age (years)", y = "Normalized nonsyn-syn ratio") + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
  ggsave(paste0(nonsyn_syn_dir, "/ratio/nonsyn_vs_syn.ratio.pdf"), plot = p_SNV_burden_lme, width = 9, height = 6, dpi = 600)
}

# write.csv(df_ratio_Func, file.path(exon_intron_dir, "2-normalized_ratio_Func.csv"), row.names = FALSE)
# write.csv(df_ratio_ExonicFunc, file.path(nonsyn_syn_dir, "2-normalized_ratio_ExonicFunc.csv"), row.names = FALSE)

################################################################################
##### mixed effects model
# df_input = df_ratio_Func
# p_traj_y_label = "Normalized exonic-intronic ratio"
# df_input = df_ratio_ExonicFunc
# p_traj_y_label = "Normalized nonsyn-syn ratio"

analyze_ratio <- function(SCAN2_df, df_input, p_traj_y_label) {
  # 1. Midpoint age centering
  mean_age <- round(mean(range(unique(SCAN2_df$Age))))
  df <- df_input %>% mutate(Age_cor = Age - mean_age)

  # 2. Fit models (with and without interaction)
  df_01 <- df %>% filter(normalized_ratio < 30)
  fit_simple   <- lmer(normalized_ratio ~ Age_cor + Condition + (1 | Case_ID), data = df_01, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
  fit_interact <- lmer(normalized_ratio ~ Age_cor * Condition + (1 | Case_ID), data = df_01, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
  
  # Test if interaction improves model
  lr_test <- anova(fit_simple, fit_interact)
  
  fit <- lmer(normalized_ratio ~ Age_cor + Condition + (1 | Case_ID), data = df, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
  
  # 3. Estimated marginal means at midpoint
  emm_mid <- emmeans(fit, ~ Condition, at = list(Age_cor = 0))
  emm_summary <- summary(emm_mid)
  emm_test <- test(emm_mid, null = 1, side = "two.sided")
  results <- left_join(as.data.frame(emm_summary), as.data.frame(emm_test) %>% select(Condition, p.value), by = "Condition") %>% 
    rename(Estimate = emmean, SE = SE, DF = df, LowerCL = lower.CL, UpperCL = upper.CL, P_value_vs1 = p.value)
  
  # 4. Predicted ratios across full age range
  age_range <- range(df_input$Age)
  age_seq <- seq(age_range[1], age_range[2], length.out = 100) - mean_age
  emm_age <- emmeans(fit, ~ Condition | Age_cor, at = list(Age_cor = age_seq), type = "response")
  emm_df <- as.data.frame(emm_age) %>% mutate(Age = Age_cor + mean_age)
  
  raw_df <- df %>% mutate(Age = Age_cor + mean_age)
  emm_obs_df <- emmeans(fit, ~ Condition | Age_cor, at = list(Age_cor = sort(unique(raw_df$Age_cor))), type = "response") %>% 
    as.data.frame() %>% mutate(Age = Age_cor + mean_age)
  obs_combos <- raw_df %>% distinct(Age, Condition)
  emm_obs_df <- inner_join(emm_obs_df, obs_combos, by = c("Age", "Condition"))
  
  p_traj <- ggplot(emm_df, aes(x = Age, y = emmean, color = Condition, fill = Condition)) + 
    geom_point(data = raw_df, aes(x = Age, y = normalized_ratio, color = Condition), alpha = 0.35, size = 5, shape = 18) + 
    geom_ribbon(data = emm_df, aes(x = Age, ymin = lower.CL, ymax = upper.CL, fill = Condition), alpha = 0.17, color = NA) +
    geom_line(data = emm_df, linetype = "solid", aes(x = Age, y = emmean, color = Condition), linewidth = 1.2) + 
    geom_errorbar(data = emm_obs_df, aes(x = Age, ymin = lower.CL, ymax = upper.CL, color = Condition), width = 1.0, linewidth = 1.0) + 
    geom_point(data = emm_obs_df, aes(x = Age, y = emmean, color = Condition), shape = 21, fill = "white", size = 2, stroke = 0.5, show.legend = FALSE) + 
    geom_hline(yintercept = 1, linetype = "11", linewidth = 1.0, color = "black") + 
    scale_fill_manual(values = c("Control" = Ctrl_IHD_color[1], "IHD" = Ctrl_IHD_color[2]), guide = "legend") + 
    scale_color_manual(values = c("Control" = Ctrl_IHD_color[1], "IHD" = Ctrl_IHD_color[2]), guide = "legend") + 
    labs(x = "Age (years)", y = p_traj_y_label) + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
  
  # Extract main Condition effect and Age slope
  ct <- coef(summary(fit))
  age_est <- ct["Age_cor", "Estimate"]
  age_p   <- ct["Age_cor", "Pr(>|t|)"]
  cond_est <- ct["ConditionIHD", "Estimate"]
  cond_p   <- ct["ConditionIHD", "Pr(>|t|)"]
  
  # Labels
  age_label <- paste0("Age slope = ", signif(age_est, 3), ", P = ", signif(age_p, 2))
  cond_label <- paste0("IHD effect = ", signif(cond_est, 3), ", P = ", signif(cond_p, 2))
  
  if (max(raw_df$normalized_ratio, na.rm = TRUE) > 30) {
    p_traj <- p_traj + scale_y_break(c(7, 42), scale = 0.2) + scale_y_continuous(limits = c(-2, 45), breaks = function(x) {
      pretty_vals <- scales::pretty_breaks()(x) 
      pretty_vals[pretty_vals %% 1 == 0]}) + 
      annotate("text", x = min(df$Age), y = max(raw_df$normalized_ratio) * 0.16, label = age_label, hjust = 0, size = 6) + 
      annotate("text", x = min(df$Age), y = max(raw_df$normalized_ratio)*0.14, label = cond_label, hjust = 0, size = 6)
  } else {
    p_traj <- p_traj + scale_y_continuous(breaks = function(x) {
      pretty_vals <- scales::pretty_breaks()(x) 
      pretty_vals[pretty_vals %% 1 == 0]}) + 
      annotate("text", x = min(df$Age), y = max(raw_df$normalized_ratio) * 0.99, label = age_label, hjust = 0, size = 6) + 
      annotate("text", x = min(df$Age), y = max(raw_df$normalized_ratio)*0.93, label = cond_label, hjust = 0, size = 6)
  }
  
  # Return results
  list(fit = fit, results_table = results, trajectory_plot = p_traj)
}

out_func <- analyze_ratio(SCAN2_df, df_ratio_Func, "Normalized exonic-intronic ratio")
out_exonic <- analyze_ratio(SCAN2_df, df_ratio_ExonicFunc, "Normalized nonsyn-syn ratio")

# View results table
# summary(out_func$fit)
# summary(out_exonic$fit)
out_func_results_table <- out_func$results_table
out_exonic_results_table <- out_exonic$results_table
write.csv(out_func_results_table, file.path(exon_intron_dir, "3-out_func_results_table.csv"), row.names = FALSE)
write.csv(out_exonic_results_table, file.path(nonsyn_syn_dir, "3-out_exonic_results_table.csv"), row.names = FALSE)

# Show trajectory plot
p_traj_plot_Func <- out_func$trajectory_plot
ggsave(paste0(exon_intron_dir, "/ratio/traj_plot_Func.pdf"), plot = p_traj_plot_Func, width = 16, height = 6, dpi = 600)
p_traj_plot_ExonicFunc <- out_exonic$trajectory_plot
ggsave(paste0(nonsyn_syn_dir, "/ratio/traj_plot_ExonicFunc.pdf"), plot = p_traj_plot_ExonicFunc, width = 16, height = 6, dpi = 600)