library(maftools)
library(psych)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)

setwd("/Users/zhemingan/Documents/BCH_research/Gene_Expression_Analysis")
##### read in metadata
Hypoxia_PTA_Cases_metadata <- readRDS("./data/SCAN2_df.rds") %>% as.data.frame() |> 
  base::`[`(c("cell_ID", "age", "gender", "Case_ID", "condition", "snv.burden", "snv.rate.per.gb")) %>% 
  rename_with(~ c("Cell_ID", "Condition"), .cols = c(1, 5))
Condition_list <- unique(Hypoxia_PTA_Cases_metadata$Condition)
selected_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", 
               "ExonicFunc.refGene", "AAChange.refGene", "Cell_ID", "Case_ID", "Condition", "mut_type", "age")

genomic_SCAN2_df <- c()
for (condition_temp in Condition_list) {
  for (mutation_type in c("ssnv", "sindel")) {
    cat("Get genomic context for", condition_temp, mutation_type, "...\n")
    heart_PTA_Cases_vcf_temp <- read.table(paste0("data/heart_PTA_Cases_annovar/heart_PTA_Cases.all_age.", condition_temp, "_", mutation_type, ".vcf"), sep = "\t") %>% mutate(V8 = sub(";.*", "", V8))
    genomic_context_temp <- read.csv(paste0("data/heart_PTA_Cases_annovar/heart_PTA_Cases.all_age.", condition_temp, "_", mutation_type, ".csv"), header = TRUE) %>%
      mutate(Cell_ID = heart_PTA_Cases_vcf_temp$V8) %>% mutate(Case_ID = str_extract(Cell_ID, "[^_]+")) %>% 
      mutate(Condition = condition_temp) %>% mutate(mut_type = mutation_type)
    
    genomic_context_temp <- genomic_context_temp %>%
      mutate(Cell_ID = ifelse(Cell_ID == "1864_M-E3-2n", "1864_E3", ifelse(Cell_ID == "4402_1_A1-2n", "4402_A1", ifelse(Cell_ID == "4402_1_A3-2n", "4402_A3",
      ifelse(Cell_ID == "4638_1_A1-2n", "4638_A1", ifelse(Cell_ID == "4638_1_A4-2n", "4638_A4", ifelse(Cell_ID == "4638_1_B2-2n", "4638_B2",
      ifelse(Cell_ID == "5657_M-D1-2n", "5657_D1", ifelse(Cell_ID == "5657_M-H1-2n", "5657_H1", ifelse(Cell_ID == "5828_CM_C2_2n", "5828_C2",
      ifelse(Cell_ID == "5828_CM_G2_2n", "5828_G2", ifelse(Cell_ID == "5919_1_C4-2n", "5919_C4", ifelse(Cell_ID == "5919_1_D2-2n", "5919_D2",
      ifelse(Cell_ID == "5919_1_E3-2n", "5919_E3", ifelse(Cell_ID == "5919_1_F6-2n", "5919_F6", ifelse(Cell_ID == "6032_1_A1-2n", "6032_A1",
      ifelse(Cell_ID == "6032_M-C7-2n", "6032_C7", ifelse(Cell_ID == "6032_M-E7-2n", "6032_E7", ifelse(Cell_ID == "1673_CM_A2_2n", "1673_A2",
      ifelse(Cell_ID == "1673_CM_A3_2n", "1673_A3", ifelse(Cell_ID == "1673_CM_D2_2n", "1673_D2", Cell_ID)))))))))))))))))))))
    
    genomic_SCAN2_df_temp <- merge(genomic_context_temp, Hypoxia_PTA_Cases_metadata[c("Cell_ID", "Case_ID", "Condition", "age")]) |> base::`[`(selected_colnames) %>% 
      # filter(age >= 40 & age < 80) %>% 
      filter(Func.refGene %in% c("splicing", "exonic;splicing") | ExonicFunc.refGene %in% c("frameshift deletion", "frameshift insertion", "nonframeshift deletion", "nonsynonymous SNV", "stopgain"))
      # filter(Func.refGene %in% c("intronic", "splicing", "exonic;splicing") | 
      #          ExonicFunc.refGene %in% c("frameshift deletion", "frameshift insertion", "nonframeshift deletion", "nonsynonymous SNV", "stopgain")) %>% 
      # filter(!grepl(";", Gene.refGene))
    
    genomic_SCAN2_df <- rbind(genomic_SCAN2_df, genomic_SCAN2_df_temp)
  }
}

# genomic_SCAN2_df <- genomic_SCAN2_df %>% filter(Gene.refGene %in% list_haha)
PTA_Cases_annovar_output <- write.table(genomic_SCAN2_df, paste0("data/MAFTools/PTA_Cases_annovar_output.txt"), sep = "\t", row.names = FALSE)

PTA_Cases_maf <- annovarToMaf(
  annovar = paste0("data/MAFTools/PTA_Cases_annovar_output.txt"), 
  refBuild = "hg19",
  # tsbCol = "Cell_ID",
  tsbCol = "Case_ID",
  table = "refGene",
  MAFobj = T,
)

# Control_num <- 7
# IHD_num <- 15
# Control_num <- 3
# IHD_num <- 5
# Control_num <- 18
# IHD_num <- 15
Control_num <- 8
IHD_num <- 5
Control.maf <- subsetMaf(maf = PTA_Cases_maf, query = "Condition == 'Normal'", mafObj = TRUE)
Control.maf@summary$summary[3] <- Control_num
IHD.maf <- subsetMaf(maf = PTA_Cases_maf, query = "Condition == 'Disease'", mafObj = TRUE)
IHD.maf@summary$summary[3] <- IHD_num
my.color <- c("#FFDE17","#E21F26","#F57F20","#2179B4","#6B3F98","#009933","#66CBE4","#010101")
names(my.color) <- c("Missense_Mutation","Nonsense_Mutation","Splice_Site","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Multi_Hit")

### select the gene list
hypoxia_gene_list <- read.csv("data/signet_sup/supplementary_table_5.csv", header = TRUE) %>% 
  pivot_longer(cols = everything(), names_to = "gene_set", values_to = "genes") %>% filter(!duplicated(genes))

hot_gene_list <- read.csv("data/MAFTools/hot_gene_list_heart_pta_02.csv", header = FALSE) %>% filter(!duplicated(V1))
all.genes <- names(sort(table(genomic_SCAN2_df$Gene.refGene), decreasing = TRUE))
top.genes <- names(sort(table(genomic_SCAN2_df$Gene.refGene), decreasing = TRUE))[1:200]

# list_haha <- intersect(genomic_SCAN2_df$Gene.refGene, a$V1)
# a<- as.data.frame(hot_gene_list)
# b <- as.data.frame(genomic_SCAN2_df$Gene.refGene)
# genomic_SCAN2_df[genomic_SCAN2_df$Gene.refGene == "LOC400684", ]
# 
list_haha_2 <- intersect(genomic_SCAN2_df$Gene.refGene, hypoxia_gene_list$genes)
# list_haha_2 <- intersect(genomic_SCAN2_df$Gene.refGene, hot_gene_list)
genomic_SCAN2_df_haha <- genomic_SCAN2_df %>% filter(Gene.refGene %in% list_haha_2)
# sort(table(genomic_SCAN2_df_haha$Gene.refGene), decreasing = T)
selected_gene_list <- list_haha_2

pdf("results/local_test/MAFTools/PTA_Cases.MAFTools.summary.pdf", width = 8, height = 5)
plotmafSummary(maf = Control.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = IHD.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf("results/local_test/MAFTools/PTA_Cases.MAFTools.oncoplot.pdf", width = 10, height = 25)
oncoplot(maf = Control.maf, genes = selected_gene_list, keepGeneOrder = F, colors = my.color)
oncoplot(maf = IHD.maf, genes = selected_gene_list, keepGeneOrder = F, colors = my.color)
dev.off()

# pdf("results/local_test/MAFTools/PTA_Cases.MAFTools.oncoplot.topgenes.pdf", width = 10, height = 35)
# oncoplot(maf = Control.maf, genes = top.genes, keepGeneOrder = F, colors = my.color)
# oncoplot(maf = IHD.maf, genes = top.genes, keepGeneOrder = F, colors = my.color)
# dev.off()

# pdf("results/local_test/MAFTools/PTA_Cases.MAFTools.oncoplot.hypoxiagenes.pdf", width = 10, height = 25)
# oncoplot(maf = Control.maf, genes = hypoxia_gene_list$Elvidge, keepGeneOrder = TRUE, colors = my.color)
# oncoplot(maf = IHD.maf, genes = hypoxia_gene_list$Elvidge, keepGeneOrder = TRUE, colors = my.color)
# dev.off()

# pdf("results/local_test/MAFTools/PTA_Cases.MAFTools.oncoplot.hotgenes.pdf", width = 10, height = 25)
# oncoplot(maf = Control.maf, genes = list_haha_2, keepGeneOrder = TRUE, colors = my.color)
# oncoplot(maf = IHD.maf, genes = list_haha_2, keepGeneOrder = TRUE, colors = my.color)
# dev.off()


# pdf("AD_panel.both.loose.MAFTools.2.pdf",width=8,height=5)
# for(i in top.genes)
# {
#   print(lollipopPlot2(m1=AD.maf,m2=control.maf,gene=i,m1_name="AD",m2_name="Control",AACol1="AAChange.refGene",AACol2="AAChange.refGene",colors=my.color,showDomainLabel=F))
#   print(lollipopPlot2(m1=AD.maf,m2=chip.maf,gene=i,m1_name="AD",m2_name="CHIP",AACol1="AAChange.refGene",AACol2="AAChange.refGene",colors=my.color,showDomainLabel=F))
# }
# dev.off()
# 
