library("dndscv")
library(dplyr)

setwd("/Users/zhemingan/Documents/annovar/dNdS/heart_PTA_Cases")
# data("dataset_simbreast", package="dndscv")
heart_PTA_Cases_Normal_vcf <- read.table("heart_PTA_Cases.all_age.normal_ssnv.vcf", sep = "\t") %>% 
  setNames(c("chr", "pos", "V3", "ref", "mut", "V6", "V7", "sampleID")) |> base::`[`(c("sampleID", "chr", "pos", "ref", "mut"))
heart_PTA_Cases_Disease_vcf <- read.table("heart_PTA_Cases.all_age.disease_ssnv.vcf", sep = "\t") %>% 
  setNames(c("chr", "pos", "V3", "ref", "mut", "V6", "V7", "sampleID")) |> base::`[`(c("sampleID", "chr", "pos", "ref", "mut"))

mutations = heart_PTA_Cases_Normal_vcf
# mutations = heart_PTA_Cases_Disease_vcf
dndsout = dndscv(mutations,max_muts_per_gene_per_sample=Inf,max_coding_muts_per_sample=Inf,outmats=T)
head(mutations)
sel_cv = dndsout$sel_cv
print(head(sel_cv), digits = 3)
# signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
# signif_genes = sel_cv[sel_cv$pallsubs_cv<0.1, c("gene_name","pallsubs_cv")]
rownames(signif_genes) = NULL
print(signif_genes)
print(dndsout$globaldnds)
head(dndsout$annotmuts)
print(dndsout$nbreg$theta)
signif_genes_localmodel = as.vector(dndsout$sel_loc$gene_name[dndsout$sel_loc$qall_loc<0.1])
print(signif_genes_localmodel)

# library("dndscv")
# data("dataset_normalskin", package="dndscv")
# data("dataset_normalskin_genes", package="dndscv")
# dndsskin = dndscv(m, gene_list=target_genes, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
# sel_cv = dndsskin$sel_cv
# print(sel_cv[sel_cv$qglobal_cv<0.1,c(1:10,19)], digits = 3)
# print(dndsskin$globaldnds, digits = 3)
# 
# library("dndscv")
# # 192 rates (used as default)
# data("submod_192r_3w", package="dndscv")
# colnames(substmodel) = c("syn","mis","non","spl")
# head(substmodel)
# 
# # 12 rates (no context-dependence)
# data("submod_12r_3w", package="dndscv")
# colnames(substmodel) = c("syn","mis","non","spl")
# head(substmodel)
# 
# # 2 rates (classic ts/tv model)
# data("submod_2r_3w", package="dndscv")
# colnames(substmodel) = c("syn","mis","non","spl")
# head(substmodel)
# 
# library("dndscv")
# data("dataset_normalskin", package="dndscv")
# data("dataset_normalskin_genes", package="dndscv")
# dndsskin_2r = dndscv(m, gene_list=target_genes, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, sm = "2r_3w")
# 
# print(dndsskin_2r$mle_submodel)
# 
# sel_cv = dndsskin_2r$sel_cv
# print(head(sel_cv[sel_cv$qglobal_cv<0.1, c(1:10,19)]), digits = 3)
# 
# AIC(dndsskin$poissmodel)
# AIC(dndsskin_2r$poissmodel)
# 
