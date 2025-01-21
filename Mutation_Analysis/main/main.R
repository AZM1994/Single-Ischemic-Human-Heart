# ##### load packages and customized functions
library(MutationalPatterns)
library(BSgenome)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(ggplot2)
library(ggpubr)
library(ggbreak)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(lme4)
library(lmerTest)
library(circlize)

##### set project folder and read vcf, metadata
# setwd("/Users/zhemingan/Documents/BCH_research/Mutation_Analysis")
source("main/modified_mutational_pattern_scripts/plot_spectrum_relative_diff.R")
source("main/modified_mutational_pattern_scripts/plot_spectrum_absolute_diff.R")
source("main/modified_mutational_pattern_scripts/plot_spectrum_absolute.R")
source("main/modified_mutational_pattern_scripts/plot_spectrum_customize.R")
source("main/modified_mutational_pattern_scripts/plot_96_profiles_absolute.R")

project_dir <- "Hypoxia_PTA_Cases"
# project_dir <- "Hypoxia_PTA_IPSC"
# project_dir <- "Hypoxia_PTA_all"
source(paste0(project_dir, "/scripts/1-setup_metadata.R"))

##### load SCAN2 results and add metadata
source(paste0(project_dir, "/scripts/1-load_scan2_results.R"))

##### sSNV burden analysis
source(paste0(project_dir, "/scripts/2-SNV_burden_analysis.R"))

##### Mutation spectrum
source(paste0(project_dir, "/scripts/3-mutation_spectrum.R"))

##### 96 mutational profile
source(paste0(project_dir, "/scripts/3-mutational_profile.R"))

##### De novo Mutational signatures extraction rank 2
# source(paste0(project_dir, "/scripts/4-mutation_signature_Denovo_rank4.R"))

##### Signature refitting loose all SBS signatures
source(paste0(project_dir, "/scripts/4-signature_refitting_loose_SigNet.R"))

##### Boxplot of statistical test of top 10 extracted signatures between two conditions
source(paste0(project_dir, "/scripts/4-statistical_analysis_top_signatures_SigNet.R"))

##### Strand bias test
source(paste0(project_dir, "/scripts/5-strand_bias_analysis.R"))

##### SNV distribution among chromosome
source(paste0(project_dir, "/scripts/6-distribution_among_chromosome.R"))

##### plot coverage
# source(paste0(project_dir, "/scripts/coverage_plot.R"))