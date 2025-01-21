library("dndscv")
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)

setwd("/Users/zhemingan/Documents/annovar/dNdS/heart_PTA_Cases")

# heart_PTA_Cases_Normal_vcf <- read.table("heart_PTA_Cases.all_age.normal_ssnv.vcf", sep = "\t") %>% 
#   setNames(c("chr", "pos", "V3", "ref", "mut", "V6", "V7", "sampleID")) |> base::`[`(c("sampleID", "chr", "pos", "ref", "mut"))
# heart_PTA_Cases_Disease_vcf <- read.table("heart_PTA_Cases.all_age.disease_ssnv.vcf", sep = "\t") %>% 
#   setNames(c("chr", "pos", "V3", "ref", "mut", "V6", "V7", "sampleID")) |> base::`[`(c("sampleID", "chr", "pos", "ref", "mut"))

SCAN2_df <- readRDS("SCAN2_df.rds") %>% 
  rbind(list("germline_bulk", 2, 0, 0, 0, 0, 0, 0, 0, 0, "FALSE", 
             0, 0, 0, 0, 0, 0, 0, 0, "FALSE", 999, "bulk", "bulk", "bulk"))
colnames(SCAN2_df)[1] <- "Cell_ID"
age_line <- 5

################################################################################
##### genomic_context for normal, disease, and germline
heart_PTA_Cases_Normal_vcf <- read.table("heart_PTA_Cases.all_age.normal_ssnv.vcf", sep = "\t")
heart_PTA_Cases_Disease_vcf <- read.table("heart_PTA_Cases.all_age.disease_ssnv.vcf", sep = "\t")

genomic_context_colnames <- c("Chr","Start","End","Ref","Alt","Cell_ID","Func.refGene","Gene.refGene")
genomic_context_normal <- read.csv("heart_PTA_Cases.all_age.normal.annotation.csv", header = TRUE) %>% 
  mutate(Cell_ID = heart_PTA_Cases_Normal_vcf$V8) %>% 
  mutate(Case_ID = str_extract(Cell_ID, "[^_]+")) |> base::`[`(genomic_context_colnames)

genomic_context_disease <- read.csv("heart_PTA_Cases.all_age.disease.annotation.csv", header = TRUE) %>% 
  mutate(Cell_ID = heart_PTA_Cases_Disease_vcf$V8) %>% 
  mutate(Case_ID = str_extract(Cell_ID, "[^_]+")) |> base::`[`(genomic_context_colnames)

genomic_context <- rbind(genomic_context_normal, genomic_context_disease)
genomic_context <- genomic_context %>% 
  mutate(Cell_ID = ifelse(Cell_ID == "1864_M-E3-2n", "1864_E3", ifelse(Cell_ID == "4402_1_A1-2n", "4402_A1", ifelse(Cell_ID == "4402_1_A3-2n", "4402_A3", 
  ifelse(Cell_ID == "4638_1_A1-2n", "4638_A1", ifelse(Cell_ID == "4638_1_A4-2n", "4638_A4", ifelse(Cell_ID == "4638_1_B2-2n", "4638_B2", 
  ifelse(Cell_ID == "5657_M-D1-2n", "5657_D1", ifelse(Cell_ID == "5657_M-H1-2n", "5657_H1", ifelse(Cell_ID == "5828_CM_C2_2n", "5828_C2", 
  ifelse(Cell_ID == "5828_CM_G2_2n", "5828_G2", ifelse(Cell_ID == "5919_1_C4-2n", "5919_C4", ifelse(Cell_ID == "5919_1_D2-2n", "5919_D2", 
  ifelse(Cell_ID == "5919_1_E3-2n", "5919_E3", ifelse(Cell_ID == "5919_1_F6-2n", "5919_F6", ifelse(Cell_ID == "6032_1_A1-2n", "6032_A1", 
  ifelse(Cell_ID == "6032_M-C7-2n", "6032_C7", ifelse(Cell_ID == "6032_M-E7-2n", "6032_E7", ifelse(Cell_ID == "1673_CM_A2_2n", "1673_A2", 
  ifelse(Cell_ID == "1673_CM_A3_2n", "1673_A3", ifelse(Cell_ID == "1673_CM_D2_2n", "1673_D2", Cell_ID)))))))))))))))))))))

genomic_SCAN2_df <- merge(genomic_context, SCAN2_df[c("Cell_ID", "Case_ID", "age", "gender", "condition")]) %>% 
  mutate(Group = ifelse(condition == "Normal" & age < age_line, "normal_young", 
                        ifelse(condition == "Normal" & age > age_line, "normal_elder", 
                               ifelse(condition == "Disease", "disease", 
                                      ifelse(Cell_ID == "germline_bulk", "germline_mutation", NA))))) %>% 
  mutate(Group = factor(Group, levels = c("normal_young", "normal_elder", "disease", "germline_mutation"))) %>% 
  filter(!is.na(Group))

sorted_counts <- sort(table(genomic_SCAN2_df$Gene.refGene), decreasing = TRUE)[1:20]
genomic_SCAN2_df[genomic_SCAN2_df$Gene.refGene == "PTPRD", ]

SNV_normal <- genomic_SCAN2_df[genomic_SCAN2_df$condition == "Normal", ] |> base::`[`(c("Cell_ID", "Chr", "Start", "Ref", "Alt")) %>% 
  setNames(c("sampleID", "chr", "pos", "ref", "mut"))
SNV_disease <- genomic_SCAN2_df[genomic_SCAN2_df$condition == "Disease", ] |> base::`[`(c("Cell_ID", "Chr", "Start", "Ref", "Alt")) %>% 
  setNames(c("sampleID", "chr", "pos", "ref", "mut"))

hot_gene_list <- read.csv("hot_gene_list_heart_pta_02.csv", header = FALSE)
hot_gene_list <- hot_gene_list[!duplicated(hot_gene_list$V1), ]
hot_gene_list <- hot_gene_list[!hot_gene_list %in% c("LPAL2", "7q22", "9p21", "AK097927", "ACO74093.1", "C6orf105", "ABO", "LOC400684", "POM121L9P")]

dndsout_normal <- dndscv(SNV_normal, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
dndsout_disease <- dndscv(SNV_disease, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
# dndsout_normal <- dndscv(SNV_normal, gene_list = hot_gene_list, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
# dndsout_disease <- dndscv(SNV_disease, gene_list = hot_gene_list, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)

# dndsout_all <- dndsout_disease
dndsout_all <- rbind(data.frame(dndsout_normal$sel_cv[, c(1:8, 11, 14)], Condition = "Normal"), 
                     data.frame(dndsout_disease$sel_cv[, c(1:8, 11, 14)], Condition = "Disease"))
colnames(dndsout_all)=c("Gene","N_synonymous","N_missense","N_nonsense","N_splicing","Ratio_missense","Ratio_nonsense","Ratio_splicing","Pvalue","Qvalue","Condition")
dndsout_all$Gene=factor(dndsout_all$Gene,levels=rev(dndsout_disease$sel_cv[,1]))
dndsout_all$Condition=factor(dndsout_all$Condition,levels=c("Normal","Disease"))
dndsout_all$Pconvert=-log10(dndsout_all$Pvalue)
dndsout_all$Qconvert=-log10(dndsout_all$Qvalue)
dndsout_all$Genetype="Other"
# dndsout_all$Genetype[dndsout_all$Gene %in% cancer_census$Symbol[cancer_census$Role=="TSG"]]="TSG"
# dndsout_all$Genetype[dndsout_all$Gene %in% sig.genes]="Hotspot"
# dndsout_all$Genetype=factor(dndsout_all$Genetype,levels=c("Hotspot","TSG","Other"))
write.table(dndsout_all[dndsout_all$Gene %in% dndsout_all$Gene[dndsout_all$Pvalue<0.05],],file="dNdScv.gene.tsv",quote=F,sep="\t",row.names=F)

dndsout_all$sum <- rowSums(dndsout_all[, c("N_synonymous", "N_missense", "N_nonsense", "N_splicing")])
sorted_df <- dndsout_all[order(dndsout_all$sum, decreasing = TRUE), ]
write.csv(sorted_df,file="sorted_df.csv",quote=F,sep="\t",row.names=F)

pdf("dNdScv.1.pdf",width = 15, height = 7)
p1 <- ggplot(dndsout_all, aes(x = Gene, y = Pconvert, fill = Condition, color = Condition)) + 
  geom_point() + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_text_repel(data = dndsout_all[dndsout_all$Pvalue < 0.05, ], aes(label = Gene), max.overlaps = Inf) + 
  scale_fill_manual(values = c("royalblue1", "indianred3")) + scale_color_manual(values = c("royalblue1", "indianred3")) + 
  theme_classic() + theme(text = element_text(size = 15), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("-log10(P-value)")
p1
p2 <- ggplot(dndsout_all, aes(x = Gene, y = Qconvert, fill = Condition, color = Condition)) + 
  geom_point() + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_text_repel(data = dndsout_all[dndsout_all$Qvalue < 0.05, ], aes(label = Gene), max.overlaps = Inf) + 
  scale_fill_manual(values = c("royalblue1", "indianred3")) + scale_color_manual(values = c("royalblue1", "indianred3")) + 
  theme_classic() + theme(text = element_text(size = 15), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("-log10(adjusted P-value)")
p2
dev.off()

dndsout_all02 <- rbind(data.frame(dndsout_normal$globaldnds, Condition = "Normal"), 
                     data.frame(dndsout_disease$globaldnds, Condition = "Disease"))
colnames(dndsout_all02) <- c("Type","MLE","LowCI","HighCI","Condition")
dndsout_all02$Type=factor(c("Missense","Nonsense","Splicing","Truncating","All"),levels=c("All","Missense","Nonsense","Splicing","Truncating"))
dndsout_all02$Condition=factor(dndsout_all02$Condition,levels=c("Normal","Disease"))


pdf("dNdScv.2.pdf",width=5.5,height=2.5)
p3 <- ggplot(dndsout_all02[dndsout_all02$Type != "Truncating", ], aes(x = Type, y = MLE, fill = Condition, color = Condition))+ 
  geom_pointrange(aes(ymin = LowCI, ymax = HighCI), position = position_dodge(width = 0.5)) + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  scale_fill_manual(values = c("royalblue1", "indianred3")) + scale_color_manual(values = c("royalblue1", "indianred3")) + 
  theme_classic() + theme(text = element_text(size = 15)) + xlab("Mutation type") + ylab("dN/dS ratio")
p3
dev.off()

# Disease_gene <- geneci(dndsout_disease, gene_list = unique(dndsout_all$Gene[dndsout_all$Qvalue < 0.05]))
Disease_gene <- geneci(dndsout_disease, gene_list = unique(dndsout_all$Gene[dndsout_all$Pvalue < 0.05]))
# Disease_gene <- geneci(dndsout_disease, gene_list = unique(dndsout_all$Gene))
# Normal_gene <- geneci(dndsout_normal, gene_list = unique(dndsout_all$Gene[dndsout_all$Qvalue < 0.05]))
Normal_gene <- geneci(dndsout_normal, gene_list = unique(dndsout_all$Gene[dndsout_all$Pvalue < 0.05]))
# Normal_gene <- geneci(dndsout_normal, gene_list = unique(dndsout_all$Gene))
missense=data.frame(Disease_gene[,c(1,2,4,6)],Condition="Disease",Type="Missense")
# missense=rbind(data.frame(Disease_gene[,c(1,2,4,6)],Condition="AD",Type="Missense"),data.frame(control_gene[,c(1,2,4,6)],Condition="CTRL",Type="Missense"))
colnames(missense)=c("Gene","MLE","LowCI","HighCI","Condition","Type")
# truncating=rbind(data.frame(Disease_gene[,c(1,3,5,7)],Condition="AD",Type="Truncating"),data.frame(control_gene[,c(1,3,5,7)],Condition="CTRL",Type="Truncating"))
truncating=data.frame(Disease_gene[,c(1,3,5,7)],Condition="Disease",Type="Truncating")
colnames(truncating)=c("Gene","MLE","LowCI","HighCI","Condition","Type")
final3=rbind(missense,truncating)
final3$Type=factor(final3$Type,levels=c("Missense","Truncating"))
final3$Condition=factor(final3$Condition,levels=c("Normal","Disease"))

pdf("dNdScv.3.pdf",width = 10, height = 4)
p4 <- ggplot(final3[1:10,],aes(x=Type,y=MLE,fill=Condition,color=Condition)) + 
  geom_pointrange(aes(ymin=LowCI,ymax=HighCI),position=position_dodge(width=0.5)) + 
  geom_hline(yintercept=1,linetype="dashed") + 
  scale_fill_manual(values=c("royalblue1","indianred3")) + scale_color_manual(values=c("royalblue1","indianred3")) + 
  theme_classic() + theme(text=element_text(size=15)) + xlab("Mutation type") + ylab("dN/dS ratio") + 
  scale_y_log10() + facet_grid(. ~ Gene)
p4
dev.off()
