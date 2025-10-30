###################################################################################################
##### 1. Read annotation of genomic location and functional categories ############################
##### 2. Calculate the exonic-to-intronic ratio and nonsynonymous-to-synonymous ratio (dN/dS) #####
#####    for middle-aged control and IHD groups, normalize the ratio by germline mutation     #####
###################################################################################################

library(ggplot2)
library(reshape2)
library(stringr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(lme4)
library(lmerTest)

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
genomic_SCAN2_df <- genomic_context %>% inner_join(SCAN2_df %>% select(Cell_ID, Case_ID, Age, Gender, Condition), by = "Cell_ID")

##### summarize all genetic regions
Func.refGene_oldnames <- c("intergenic", "upstream", "UTR5", "exonic", "UTR3", "downstream", "splicing", "intronic")
Func.refGene_newnames <- c("intergenic", "upstream", "5' UTR", "exonic", "3' UTR", "downstream", "splice site", "intronic")
Func.refGene_Group_df <- table(factor(genomic_SCAN2_df$Func.refGene, levels = Func.refGene_oldnames), genomic_SCAN2_df$Cell_ID) %>% `row.names<-`(Func.refGene_newnames)

Func.refGene_Group_df_melt <- melt(Func.refGene_Group_df) %>% 
  setNames(c("Func.refGene","Cell_ID","Count")) %>% 
  mutate(Func.refGene = factor(Func.refGene, levels = Func.refGene_newnames)) %>% 
  mutate(Total = rep(table(genomic_SCAN2_df$Cell_ID), each = length(unique(Func.refGene)))) %>% 
  mutate(Prop = Count / Total)

##### plot sSNV burden by each type
sSNV_Func.refGene_CellID <- SCAN2_df %>% 
  left_join(Func.refGene_Group_df_melt %>% select(c("Func.refGene", "Cell_ID", "Prop")), by = "Cell_ID", relationship = "many-to-many") %>% 
  mutate(value = snv.rate.per.gb * Prop)

df_wide_Func <- sSNV_Func.refGene_CellID %>% 
  select(Age, Condition, Case_ID, Gender, Cell_ID, Func.refGene, value, Color, Outline_Color) %>% 
  pivot_wider(names_from = Func.refGene, values_from = value, values_fill = list(value = 0)) %>% 
  mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% 
  mutate(Condition = relevel(Condition, ref = ctrl_name)) %>% 
  group_by(Condition) %>% arrange(Age, .by_group = TRUE)

selected_func_list <- colnames(df_wide_Func)[8 : (ncol(df_wide_Func))]
sig_digits <- 2
pval_df_Func <- data.frame(Func = character(), estimate = numeric(), SE = numeric(), 
                           df = numeric(), t = numeric(), pval = numeric(), stringsAsFactors = FALSE)

##### Mix effects model & Age matched mix effects model
for (i in c(1 : length(selected_func_list))){
  cat("index", i, ":", selected_func_list[i])
  ## burden_df: the input matrix with columns: burden, Condition, Case_ID, Age, Color
  burden_df <- df_wide_Func[c(selected_func_list[i], "Age", "Condition", "Case_ID", "Color", "Outline_Color")] %>% 
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
  
  # Age effect p-value (from control-only model)
  aov_tab_ctrl <- anova(burden_age_model_ctrl)
  if ("Age" %in% rownames(aov_tab_ctrl)) {
    pval_age <- aov_tab_ctrl["Age", "Pr(>F)"]
  } else {
    pval_age <- NA_real_
  }
  
  pval_df_Func <- rbind(pval_df_Func, data.frame(Func = selected_func_list[i], estimate = est, SE = SE, 
                                                 ci_low=ci_low, ci_high=ci_high, pval_condition = anova_pval, pval_age = pval_age))
  
  geom_line_data <- generate_line_data(burden_df_ctrl, coef_ctrl)
  legend_data <- burden_df[c(1, nrow(burden_df)), c("Age", "snv.rate.per.gb", "Condition")]
  
  yticks_plot <- generate_y_breaks(burden_df$snv.rate.per.gb)
  ylim_plot <- range(yticks_plot, na.rm = TRUE)
  p_SNV_burden_lme <- ggplot(burden_df, aes(x = Age, y = snv.rate.per.gb)) + 
    geom_point(pch = 21, data = legend_data, aes(x = Age, y = snv.rate.per.gb, color = Condition, fill = Condition), size = 5) + 
    geom_point(pch = 21, color = burden_df$Outline_Color, fill = burden_df$Color, size = 5) + 
    geom_line(aes(x = ctrl_fitting_x, y = ctrl_fitting_y), color = "dodgerblue2", linetype = "22", data = geom_line_data) + 
    geom_hline(yintercept = 0, color = "gray40", linetype = "22") + 
    annotate("text", size = 6, x = 0, y = 0.95 * ylim_plot[2], 
             label = paste("Aging effect:", aging_rate, "sSNVs/(GB·year),", "P =", anova_pval_print_ctrl), hjust = 0) +
    annotate("text", size = 6, x = 0, y = 0.87 * ylim_plot[2], 
             label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB,", "P =", anova_pval_print), hjust = 0) +
    scale_color_manual(values = c("IHD" = "black", "Control" = "black"), guide = "legend") + 
    scale_fill_manual(values = c("Control" = Ctrl_IHD_color[1], "IHD" = Ctrl_IHD_color[2]), guide = "legend") + 
    scale_y_continuous(limits = c(min(min(geom_line_data[, 2]), 0), ylim_plot[2]), breaks = yticks_plot, labels = yticks_plot) + 
    labs(x = "Age (years)", y = paste0(selected_func_list[i], " Contribution \n (sSNV rate per GB)")) + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
  ggsave(paste0(exon_intron_dir, "/lme/", selected_func_list[i], "_SNV_burden_lme.pdf"), plot = p_SNV_burden_lme, width = 9, height = 5, dpi = 600)
}

pval_df_Func$pval_condition_BH <- p.adjust(pval_df_Func$pval_condition, method = "BH")
pval_df_Func$pval_age_BH <- p.adjust(pval_df_Func$pval_age, method = "BH")
write.csv(pval_df_Func, file.path(exon_intron_dir, "Func_refGene_LME_with_MTC.csv"), row.names = FALSE)

################################################################################
##################### nonsynonymous vs synonymous analysis #####################
################################################################################

##### exonic_context for Control, IHD, and germline
exonic_context_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Cell_ID", "ExonicFunc.refGene", "AAChange.refGene")
exonic_context_ctrl <- read.csv("data/annotation_results/all_age_sSNV/Control/heart_PTA_all.all_age.Control_ssnv.hg19_multianno.csv", header = TRUE) %>% 
  mutate(Cell_ID = heart_PTA_all_Control_vcf$V8) %>% filter(ExonicFunc.refGene != ".") |> dplyr::select(all_of(exonic_context_colnames))
exonic_context_dis <- read.csv("data/annotation_results/all_age_sSNV/IHD/heart_PTA_all.all_age.IHD_ssnv.hg19_multianno.csv", header = TRUE) %>% 
  mutate(Cell_ID = heart_PTA_all_IHD_vcf$V8) %>% filter(ExonicFunc.refGene != ".") |> dplyr::select(all_of(exonic_context_colnames))

exonic_context <- rbind(exonic_context_ctrl, exonic_context_dis)
exonic_SCAN2_df <- exonic_context %>% inner_join(SCAN2_df %>% select(Cell_ID, Case_ID, Age, Gender, Condition), by = "Cell_ID")

################################################################################
##### summarize all mutation types
ExonicFunc.refGene_oldnames <- c("nonsynonymous SNV", "synonymous SNV", "stopgain", "stoploss")
ExonicFunc.refGene_newnames <- c("nonsynonymous", "synonymous", "stopgain", "stoploss")
ExonicFunc.refGene_Group_df <- table(factor(exonic_SCAN2_df$ExonicFunc.refGene, levels = ExonicFunc.refGene_oldnames), exonic_SCAN2_df$Cell_ID) %>% 
  `row.names<-`(ExonicFunc.refGene_newnames)

ExonicFunc.refGene_Group_df_melt <- melt(ExonicFunc.refGene_Group_df) %>% 
  setNames(c("ExonicFunc.refGene","Cell_ID","Count")) %>% 
  mutate(ExonicFunc.refGene = factor(ExonicFunc.refGene, levels = ExonicFunc.refGene_newnames)) %>% 
  mutate(Total = rep(table(exonic_SCAN2_df$Cell_ID), each = length(unique(ExonicFunc.refGene)))) %>% 
  mutate(Prop = Count / Total)

p_exon <- Func.refGene_Group_df_melt %>% filter(Func.refGene == "exonic") %>% 
  select(Cell_ID, Prop) %>% setNames(c("Cell_ID", "exon_Prop"))

#####
SNV_ExonicFunc.refGene_CellID <- SCAN2_df %>% 
  left_join(ExonicFunc.refGene_Group_df_melt %>% select(c("ExonicFunc.refGene", "Cell_ID", "Prop")), by = "Cell_ID", relationship = "many-to-many") %>% 
  left_join(p_exon, by = "Cell_ID") %>% 
  mutate(value = snv.rate.per.gb * Prop * exon_Prop)

df_wide_ExonicFunc <- SNV_ExonicFunc.refGene_CellID %>% 
  select(Age, Condition, Case_ID, Gender, Cell_ID, ExonicFunc.refGene, value, Color, Outline_Color) %>% 
  pivot_wider(names_from = ExonicFunc.refGene, values_from = value, values_fill = list(value = 0)) %>% 
  mutate(Condition = factor(Condition, levels = c(ctrl_name, dis_name))) %>% 
  mutate(Condition = relevel(Condition, ref = ctrl_name)) %>% 
  group_by(Condition) %>% arrange(Age, .by_group = TRUE)

selected_func_list <- colnames(df_wide_ExonicFunc)[8 : (ncol(df_wide_ExonicFunc)-1)]
sig_digits <- 2

pval_df_ExonicFunc <- data.frame(
  ExonicFunc = character(),
  estimate = numeric(),
  SE = numeric(),
  df = numeric(),
  t = numeric(),
  pval = numeric(),
  stringsAsFactors = FALSE
)

##### Mix effects model & Age matched mix effects model
for (i in c(1 : length(selected_func_list))){
  cat("index", i, ":", selected_func_list[i])
  ## burden_df: the input matrix with columns: burden, Condition, Case_ID, Age, Color
  burden_df <- df_wide_ExonicFunc[c(selected_func_list[i], "Age", "Condition", "Case_ID", "Color", "Outline_Color")] %>% 
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
  
  # Age effect p-value (from control-only model)
  aov_tab_ctrl <- anova(burden_age_model_ctrl)
  if ("Age" %in% rownames(aov_tab_ctrl)) {
    pval_age <- aov_tab_ctrl["Age", "Pr(>F)"]
  } else {
    pval_age <- NA_real_
  }
  
  pval_df_ExonicFunc <- rbind(pval_df_ExonicFunc, data.frame(ExonicFunc = selected_func_list[i], estimate = est, SE = SE, 
                                                 ci_low=ci_low, ci_high=ci_high, pval_condition = anova_pval, pval_age = pval_age))
  
  geom_line_data <- generate_line_data(burden_df_ctrl, coef_ctrl)
  legend_data <- burden_df[c(1, nrow(burden_df)), c("Age", "snv.rate.per.gb", "Condition")]
  
  yticks_plot <- generate_y_breaks(burden_df$snv.rate.per.gb)
  ylim_plot <- range(yticks_plot, na.rm = TRUE)
  p_SNV_burden_lme <- ggplot(burden_df, aes(x = Age, y = snv.rate.per.gb)) + 
    geom_point(pch = 21, data = legend_data, aes(x = Age, y = snv.rate.per.gb, color = Condition, fill = Condition), size = 5) + 
    geom_point(pch = 21, color = burden_df$Outline_Color, fill = burden_df$Color, size = 5) + 
    geom_line(aes(x = ctrl_fitting_x, y = ctrl_fitting_y), color = "dodgerblue2", linetype = "22", data = geom_line_data) + 
    geom_hline(yintercept = 0, color = "gray40", linetype = "22") + 
    annotate("text", size = 6, x = 0, y = 0.95 * ylim_plot[2], 
             label = paste("Aging effect:", aging_rate, "sSNVs/(GB·year),", "P =", anova_pval_print_ctrl), hjust = 0) +
    annotate("text", size = 6, x = 0, y = 0.87 * ylim_plot[2], 
             label = paste("IHD effect:", format(round(fixef(burden_age_model)[3], digits = sig_digits), nsmall = sig_digits), "sSNVs/GB,", "P =", anova_pval_print), hjust = 0) +
    scale_color_manual(values = c("IHD" = "black", "Control" = "black"), guide = "legend") + 
    scale_fill_manual(values = c("Control" = Ctrl_IHD_color[1], "IHD" = Ctrl_IHD_color[2]), guide = "legend") + 
    scale_y_continuous(limits = c(min(min(geom_line_data[, 2]), 0), ylim_plot[2]), breaks = yticks_plot, labels = yticks_plot) + 
    labs(x = "Age (years)", y = paste0(selected_func_list[i], " Contribution \n (sSNV rate per GB)")) + theme_linedraw() + 
    theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 0.5), text = element_text(size = 24), axis.title.x = element_text(hjust = 0.5, vjust = 0))
  ggsave(paste0(nonsyn_syn_dir, "/lme/", selected_func_list[i], "_SNV_burden_lme.pdf"), plot = p_SNV_burden_lme, width = 9, height = 5, dpi = 600)
}

pval_df_ExonicFunc$pval_condition_BH <- p.adjust(pval_df_ExonicFunc$pval_condition, method = "BH")
pval_df_ExonicFunc$pval_age_BH <- p.adjust(pval_df_ExonicFunc$pval_age, method = "BH")
write.csv(pval_df_ExonicFunc, file.path(nonsyn_syn_dir, "ExonicFunc_refGene_LME_with_MTC.csv"), row.names = FALSE)
