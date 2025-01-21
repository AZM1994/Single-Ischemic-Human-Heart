library(goseq)
library(ggplot2)
library(reshape2)
library(stringr)

setwd("/Users/zhemingan/Documents/BCH_research/annovar/heart_PTA_Cases/GOseq_Zheming")
set.seed(1000)

my.GOseq <- function(input_hit,input_length,output_file,permutation){
	pwf=nullp(input_hit,bias.data=input_length)
	GO=goseq(pwf,"hg19","geneSymbol")
	GO_filtered=GO[GO$numDEInCat>=numDEInCat_threshold&GO$numInCat<=numInCat_threshold,]
	GO_filtered$over_represented_pvalue_corrected=p.adjust(GO_filtered$over_represented_pvalue,method="fdr")
	if(!permutation){#print(head(GO_filtered))
		write.table(GO_filtered,file=output_file,quote=F,sep="\t",row.names=F)
	  }
	GO
	}

my.indexDiff <- function(go,input1,input2){
	a=which(input1$category==go)/nrow(input1)
	b=which(input2$category==go)/nrow(input2)
	if(length(a)==0){
		a=1
		}
	if(length(b)==0){
		b=1
		}
	a-b
	}

my.fisherTest <- function(a,b,c,d){
	if(a/c>=b/d){
		fisher.test(matrix(c(a,b,c-a,d-b),nrow=2),alternative="greater")$p.value
	  }
	else{
		fisher.test(matrix(c(a,b,c-a,d-b),nrow=2),alternative="less")$p.value
	  }
  }

permutation_round = 2
# numDEInCat_threshold = 10
# numInCat_threshold = 1000

numDEInCat_threshold = 3
numInCat_threshold = 300

##### GOseq for exonic only
type_list = c("exonic")
normal_raw_tsv_name = c("GO_results/normal.exonic_only.GO.tsv")
disease_raw_tsv_name = c("GO_results/disease.exonic_only.GO.tsv")
permutation_table_name = c("GO_results/disease_vs_normal.exonic_only.raw.tsv")
go_table_name = c("GO_results/disease_vs_normal.exonic_only.GO.tsv")
figure_name = c("GO_results/disease_vs_normal.exonic_only.GOseq.pdf")
figure_name_png = c("GO_results/disease_vs_normal.exonic_only.GOseq")
gene_length_type = "Exon_length"

##### GOseq for both exonic and intronic
# type_list = c("exonic", "intronic")
# normal_raw_tsv_name = c("GO_results/normal.exonic_and_intronic.GO.tsv")
# disease_raw_tsv_name = c("GO_results/disease.exonic_and_intronic.GO.tsv")
# permutation_table_name = c("GO_results/disease_vs_normal.exonic_and_intronic.raw.tsv")
# go_table_name = c("GO_results/disease_vs_normal.exonic_and_intronic.GO.tsv")
# figure_name = c("GO_results/disease_vs_normal.exonic_and_intronic.GOseq.pdf")
# figure_name_png = c("GO_results/disease_vs_normal.exonic_and_intronic.GOseq")
# gene_length_type = "Gene_length"

##### load disease and control data, add metadata
# genomic_context_normal <- read.csv("heart_PTA_Cases_Normal.SNV_gene_list.csv", header = TRUE)
# heart_PTA_Cases_Normal_vcf <- read.table("heart_PTA_Cases.all_normal_ssnv.vcf", sep = "\t")
# genomic_context_normal$Cell_ID <- heart_PTA_Cases_Normal_vcf$V8
# genomic_context_normal$Case_ID <- str_extract(genomic_context_normal$Cell_ID, "[^_]+")
# genomic_context_normal <- genomic_context_normal[c("Chr","Start","End","Ref","Alt","Cell_ID","Func.refGene","Gene.refGene")]
# colnames(genomic_context_normal)=c("Chr","Start","End","Ref","Alt","Cell_ID","Type","Gene_symbol")
# 
# genomic_context_disease <- read.csv("heart_PTA_Cases_Disease.SNV_gene_list.csv", header = TRUE)
# heart_PTA_Cases_Disease_vcf <- read.table("heart_PTA_Cases.all_disease_ssnv.vcf", sep = "\t")
# genomic_context_disease$Cell_ID <- heart_PTA_Cases_Disease_vcf$V8
# genomic_context_disease$Case_ID <- str_extract(genomic_context_disease$Cell_ID, "[^_]+")
# genomic_context_disease <- genomic_context_disease[c("Chr","Start","End","Ref","Alt","Cell_ID","Func.refGene","Gene.refGene")]
# colnames(genomic_context_disease)=c("Chr","Start","End","Ref","Alt","Cell_ID","Type","Gene_symbol")
### read gene list with metadata
Hypoxia_PTA_Cases_metadata <- readRDS("SCAN2_df.rds") %>% as.data.frame() |> 
  base::`[`(c("cell_ID", "age", "gender", "Case_ID", "condition", "snv.burden", "snv.rate.per.gb")) %>% 
  rename_with(~ c("Cell_ID", "Condition"), .cols = c(1, 5))
Condition_list <- unique(Hypoxia_PTA_Cases_metadata$Condition)
selected_colnames <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", 
                       "ExonicFunc.refGene", "AAChange.refGene", "Cell_ID", "Case_ID", "Condition", "mut_type", "age")

genomic_SCAN2_df <- c()
for (condition_temp in Condition_list) {
  for (mutation_type in c("ssnv", "sindel")) {
    cat("Get genomic context for", condition_temp, mutation_type, "...\n")
    heart_PTA_Cases_vcf_temp <- read.table(paste0("heart_PTA_Cases_annovar/heart_PTA_Cases.all_age.", condition_temp, "_", mutation_type, ".vcf"), sep = "\t") %>% mutate(V8 = sub(";.*", "", V8))
    genomic_context_temp <- read.csv(paste0("heart_PTA_Cases_annovar/heart_PTA_Cases.all_age.", condition_temp, "_", mutation_type, ".csv"), header = TRUE) %>%
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

genomic_context_normal <- genomic_SCAN2_df[genomic_SCAN2_df$Condition == "Normal", ]
genomic_context_disease <- genomic_SCAN2_df[genomic_SCAN2_df$Condition == "Disease", ]

##### read gene length
gene_length <- read.delim("hg19_refGene.length.tsv", header = F)
colnames(gene_length) <- c("Gene_symbol", "Transcript", "Gene_length", "Exon_length")
gene_length_deduped <- gene_length[!duplicated(gene_length$Gene_symbol),]

##### GOseq raw
normal_exon_ready <- gene_length_deduped$Gene_symbol %in% genomic_context_normal$Gene_symbol[genomic_context_normal$Type == type_list]
names(normal_exon_ready) <- gene_length_deduped$Gene_symbol
normal_raw <- my.GOseq(normal_exon_ready, gene_length_deduped[, gene_length_type], normal_raw_tsv_name, FALSE)

disease_exon_ready <- gene_length_deduped$Gene_symbol %in% genomic_context_disease$Gene_symbol[genomic_context_disease$Type == type_list]
names(disease_exon_ready) <- gene_length_deduped$Gene_symbol
disease_raw <- my.GOseq(disease_exon_ready,gene_length_deduped[, gene_length_type], disease_raw_tsv_name, FALSE)

##### GOseq filtered
normal_go_list=as.character(normal_raw$category[normal_raw$over_represented_pvalue<0.01&normal_raw$numDEInCat>=numDEInCat_threshold&normal_raw$numInCat<=numInCat_threshold])
disease_go_list=as.character(disease_raw$category[disease_raw$over_represented_pvalue<0.01&disease_raw$numDEInCat>=numDEInCat_threshold&disease_raw$numInCat<=numInCat_threshold])
go_list=unique(c(normal_go_list, disease_go_list))
gene_list=unique(c(names(normal_exon_ready)[normal_exon_ready],names(disease_exon_ready)[disease_exon_ready]))
gene_map=getgo(gene_list,"hg19","geneSymbol")
idx_diff=sapply(go_list,function(x){my.indexDiff(x, normal_raw[normal_raw$numDEInCat>=numDEInCat_threshold&normal_raw$numInCat<=numInCat_threshold,],
                                                 disease_raw[disease_raw$numDEInCat>=numDEInCat_threshold&disease_raw$numInCat<=numInCat_threshold,])})

##### GOseq permutation
idx_diff_permutation=numeric(0)
for(i in 1:permutation_round){
  print(i)
  normal_exon_ready <- gene_length_deduped$Gene_symbol %in% genomic_context_normal$Gene_symbol[genomic_context_normal$Type == type_list]
  normal_exon_ready <- sample(normal_exon_ready)
  names(normal_exon_ready) <- gene_length_deduped$Gene_symbol
  normal_permutation <- my.GOseq(normal_exon_ready, gene_length_deduped[, gene_length_type], NA, TRUE)
  
  disease_exon_ready <- gene_length_deduped$Gene_symbol %in% genomic_context_disease$Gene_symbol[genomic_context_disease$Type == type_list]
  disease_exon_ready <- sample(disease_exon_ready)
  names(disease_exon_ready) <- gene_length_deduped$Gene_symbol
  disease_permutation <- my.GOseq(disease_exon_ready, gene_length_deduped[, gene_length_type], NA, TRUE)
  
  idx_diff_permutation <- rbind(idx_diff_permutation,sapply(go_list,function(x){
    my.indexDiff(x, normal_permutation[normal_permutation$numDEInCat>=numDEInCat_threshold&normal_permutation$numInCat<=numInCat_threshold,],
                 disease_permutation[disease_permutation$numDEInCat>=numDEInCat_threshold&disease_permutation$numInCat<=numInCat_threshold,])}))
}

write.table(idx_diff_permutation, file= permutation_table_name, quote = F, sep = "\t", row.names = F)

##### build GO table
idx_diff_permutation <- read.delim(permutation_table_name, header = T)
go_table <- data.frame(go_list, 
                       sapply(go_list, function(x){sum(genomic_context_normal$Gene_symbol[which(genomic_context_normal$Type == type_list)] %in% names(gene_map)[grep(x, gene_map)], na.rm = T)}), 
                       sapply(go_list, function(x){sum(genomic_context_disease$Gene_symbol[which(genomic_context_disease$Type == type_list)] %in% names(gene_map)[grep(x, gene_map)], na.rm = T)}), 
                       sum(genomic_context_normal$Type == type_list, na.rm = T), 
                       sum(genomic_context_disease$Type == type_list, na.rm = T), 
                       sapply(go_list, function(x){normal_raw$numDEInCat[normal_raw$category == x]}), 
                       sapply(go_list, function(x){disease_raw$numDEInCat[disease_raw$category == x]}), 
                       sapply(go_list, function(x){disease_raw$numInCat[disease_raw$category == x]}), 
                       sapply(go_list, function(x){normal_raw$over_represented_pvalue[normal_raw$category == x]}), 
                       sapply(go_list, function(x){disease_raw$over_represented_pvalue[disease_raw$category == x]}), 
                       sapply(go_list, function(x){disease_raw$term[disease_raw$category == x]}), 
                       sapply(go_list, function(x){disease_raw$ontology[disease_raw$category == x]}))

colnames(go_table)=c("category", "hit_normal", "hit_disease", "total_normal", "total_disease",
                     "numDE_normal", "numDE_disease", "numInCat", "p_normal", "p_disease", "term", "ontology")

go_table$corrected_p_normal <- p.adjust(go_table$p_normal, method = "fdr")
go_table$corrected_p_disease <- p.adjust(go_table$p_disease, method = "fdr")
#go_table$p_fisher=sapply(1:nrow(go_table),function(x){my.fisherTest(go_table$hit_smoking[x],go_table$hit_control[x],go_table$total_smoking[x],go_table$total_control[x])})
go_table$p_permutation <- sapply(1 : nrow(go_table), function(x){
  if(idx_diff[x] < 0){
    sum(idx_diff[x] > idx_diff_permutation[, x]) / length(idx_diff_permutation[, x])
  }
  else{
    sum(idx_diff[x] < idx_diff_permutation[, x]) / length(idx_diff_permutation[, x])
  }
})
go_table$corrected_p_permutation <- p.adjust(go_table$p_permutation, method = "fdr")
write.table(go_table, file = go_table_name, quote = F, sep = "\t", row.names = F) 

##### enrichment plot
pdf(figure_name, height = 10, width = 8)
ready <- melt(go_table[go_table$corrected_p_disease < 0.05|go_table$corrected_p_normal < 0.05, c(1,11,13,14)])
colnames(ready) <- c("GO", "Term", "Variable", "Pvalue")
ready$Term <- factor(ready$Term, levels = rev(unique(ready$Term)))
ready$Condition <- factor(sapply(ready$Variable, function(x){
  if(grepl("normal", x)){
    "normal"
  }
  else if(grepl("disease", x)){
    "disease"
  }
}), levels = c("normal", "disease"))

ready$Group <- factor(sapply(ready$GO, function(x){
  if(x %in% normal_go_list){
    "normal-enriched"
  }else if(x %in% disease_go_list){
    "disease-enriched"
  }
}), levels = c("normal-enriched", "disease-enriched"))

ready$Pconvert <- -log10(ready$Pvalue)
p1 <- ggplot(ready, aes(x = Term, y = Pconvert, fill = Condition)) + 
  geom_col(position="dodge") + geom_hline(yintercept=-log10(0.05),linetype="dashed") + 
  scale_fill_discrete(name="") + 
  scale_fill_manual(values = c("normal" = "blue", "disease" = "red")) +
  xlab("") + ylab("-log10(FDR-adjusted P-value)") + 
  coord_flip() + theme(legend.position="bottom") + facet_grid(Group ~ .,space="free",scales="free")
print(p1)
ggsave(paste0(figure_name_png, "_01.png"),
       plot = p1, width = 8, height = 10, dpi = 300)

# new plot according to Mike Miller's request
#p <- ggplot(ready,aes(x=Term,y=Pconvert,fill=Tissue))
#p + geom_col(position="dodge",color="black") + geom_hline(yintercept=-log10(0.05),linetype="dashed") + scale_fill_manual(values=c("indianred3","royalblue1"),name="") + xlab("") + ylab("-log10(FDR-adjusted P-value)") + coord_flip() + theme_classic() + theme(legend.position="bottom",text=element_text(size=16))
#p <- ggplot(ready,aes(x=Term,y=Pconvert,fill=factor(Tissue,levels=c("control","smoking"))))
#p + geom_col(position="dodge",color="black") + geom_hline(yintercept=-log10(0.05),linetype="dashed") + scale_fill_manual(values=c("royalblue1","indianred3"),name="") + xlab("") + ylab("-log10(FDR-adjusted P-value)") + coord_flip() + theme_classic() + theme(legend.position="bottom",text=element_text(size=16))

ready <- melt(go_table[go_table$corrected_p_permutation < 0.05, c(1,11,13,14)])
colnames(ready) <- c("GO", "Term", "Variable", "Pvalue")
ready$Term <- factor(ready$Term, levels = rev(unique(ready$Term)))
ready$Condition <- factor(sapply(ready$Variable, function(x){
  if(grepl("normal", x)){
    "normal"
  }
  else if(grepl("disease", x)){
    "disease"
  }
}), levels = c("normal", "disease"))

ready$Group <- factor(sapply(ready$GO, function(x){
  if(x %in% normal_go_list){
    "normal-enriched"
  }else if(x %in% disease_go_list){
    "disease-enriched"
  }
}), levels = c("normal-enriched", "disease-enriched"))

ready$Pconvert <- -log10(ready$Pvalue)
p2 <- ggplot(ready, aes(x = Term, y = Pconvert, fill = Condition)) + 
  geom_col(position="dodge") + geom_hline(yintercept=-log10(0.05),linetype="dashed") + 
  scale_fill_discrete(name="") + 
  scale_fill_manual(values = c("normal" = "blue", "disease" = "red")) +
  xlab("") + ylab("-log10(FDR-adjusted P-value)") + 
  coord_flip() + theme(legend.position="bottom") + facet_grid(Group ~ .,space="free",scales="free")
print(p2)
ggsave(paste0(figure_name_png, "_02.png"),
       plot = p2, width = 8, height = 10, dpi = 300)

ready <- go_table[go_table$corrected_p_permutation < 0.05,c(1,11,16)]
colnames(ready)=c("GO","Term","Pvalue")
ready$Condition <- factor(sapply(ready$GO, function(x){
  if(x %in% normal_go_list){
    "normal-enriched"
  }else if(x %in% disease_go_list){
    "disease-enriched"
  }
}), levels = c("normal-enriched", "disease-enriched"))

ready$Pconvert <- -log10(ready$Pvalue) * sapply(ready$Condition, function(x){
  if(x == "normal-enriched"){
    -1
  }else if(x == "disease-enriched"){
    1
  }
})
ready$Pconvert[ready$Pconvert == Inf] <- log10(permutation_round)
ready$Pconvert[ready$Pconvert == -Inf] <- -log10(permutation_round)
ready$Term <- factor(ready$Term, levels = ready$Term[order(ready$Condition, -ready$Pvalue)])
p3 <- ggplot(ready, aes(x = Term, y = Pconvert, fill = Condition)) + 
  geom_col() + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  scale_fill_manual(values = c("normal-enriched" = "blue", "disease-enriched" = "red")) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") + xlab("") + ylab("-log10(FDR-adjusted P-value)") + 
  coord_flip() + theme(legend.position = "none") + facet_grid(Condition ~ ., space = "free", scales = "free")
print(p3)
ggsave(paste0(figure_name_png, "_03.png"),
       plot = p3, width = 10, height = 5, dpi = 300)

dev.off()
