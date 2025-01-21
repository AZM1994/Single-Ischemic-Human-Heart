library(ggplot2)
library(maftools)
library(psych)

my.sensitivity <- function(x,sensitivity_table)
{
	index=as.numeric(cut(x,sensitivity_table$MAF))
	return(sensitivity_table$Prop[index]+(x-sensitivity_table$MAF[index])/(sensitivity_table$MAF[index+1]-sensitivity_table$MAF[index])*(sensitivity_table$Prop[index+1]-sensitivity_table$Prop[index]))
}

AD_num=190
control_num=121

panel_design=read.delim("3277911_Covered.bed",header=F,stringsAsFactors=F)[-(1:2),]
colnames(panel_design)=c("Chr","Start","End","Symbol")

top.genes=c("DNMT3A","TET2","ASXL1","KMT2D","ATRX","BCR","CBL","TP53","MLH1","STAT3")

my.color=c("#FFDE17","#E21F26","#F57F20","#2179B4","#6B3F98","#009933","#66CBE4","#010101")
names(my.color)=c("Missense_Mutation","Nonsense_Mutation","Splice_Site","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Multi_Hit")

all.maf=annovarToMaf(annovar="AD_panel.both.loose.hg19_multianno.txt",refBuild="hg19",tsbCol="ID",table="refGene",MAFobj=T)

AD.maf=subsetMaf(maf=all.maf,query="Cogdx=='AD'",mafObj=T)
AD.maf@summary$summary[3]=AD_num
control.maf=subsetMaf(maf=all.maf,query="Cogdx=='Normal'",mafObj=T)
control.maf@summary$summary[3]=control_num
chip.maf=annovarToMaf(annovar="CHIP.hg38_multianno.txt",refBuild="hg38",tsbCol="ID",table="refGene",MAFobj=T)

pdf("AD_panel.both.loose.MAFTools.1.pdf",width=8,height=5)

oncoplot(maf=AD.maf,genes=top.genes,keepGeneOrder=T,colors=my.color)
oncoplot(maf=control.maf,genes=top.genes,keepGeneOrder=T,colors=my.color)

dev.off()

pdf("AD_panel.both.loose.MAFTools.2.pdf",width=8,height=5)

for(i in top.genes)
{
	print(lollipopPlot2(m1=AD.maf,m2=control.maf,gene=i,m1_name="AD",m2_name="Control",AACol1="AAChange.refGene",AACol2="AAChange.refGene",colors=my.color,showDomainLabel=F))
	print(lollipopPlot2(m1=AD.maf,m2=chip.maf,gene=i,m1_name="AD",m2_name="CHIP",AACol1="AAChange.refGene",AACol2="AAChange.refGene",colors=my.color,showDomainLabel=F))
}

dev.off()

AD.pathway=OncogenicPathways(maf=AD.maf)
control.pathway=OncogenicPathways(maf=control.maf)
merged.pathway=merge(AD.pathway,control.pathway,by="Pathway",all=T)

pathdb=system.file("extdata","oncogenic_sig_patwhays.tsv",package="maftools")
pathdb=data.table::fread(input=pathdb)

merged.pathway$Mutated_samples.x[is.na(merged.pathway$Mutated_samples.x)]=0
merged.pathway$Mutated_samples.y[is.na(merged.pathway$Mutated_samples.y)]=0
merged.pathway$P_value=sapply(1:nrow(merged.pathway),function(x){prop.test(c(merged.pathway$Mutated_samples.x[x],merged.pathway$Mutated_samples.y[x]),c(190,121),alternative="greater")$p.value})
merged.pathway[merged.pathway$P_value<0.1]

pi3k.list=intersect(panel_design$Symbol,pathdb$Gene[pathdb$Pathway=="PI3K"])
pi3k.nums.AD=unlist(sapply(pi3k.list,function(x){AD.maf@gene.summary$total[AD.maf@gene.summary$Hugo_Symbol==x]}))
pi3k.nums.control=unlist(sapply(pi3k.list,function(x){control.maf@gene.summary$total[control.maf@gene.summary$Hugo_Symbol==x]}))
pi3k.genes=unique(c(names(pi3k.nums.AD)[order(pi3k.nums.AD,decreasing=T)],names(pi3k.nums.control)[order(pi3k.nums.control,decreasing=T)]))

pdf("AD_panel.both.loose.MAFTools.3.pdf",width=8,height=5)

oncoplot(maf=AD.maf,genes=pi3k.genes,keepGeneOrder=T,colors=my.color)
oncoplot(maf=control.maf,genes=pi3k.genes,keepGeneOrder=T,colors=my.color)

PlotOncogenicPathways(maf=AD.maf,pathways="PI3K")
PlotOncogenicPathways(maf=control.maf,pathways="PI3K")

dev.off()

my.dgi=drugInteractions(genes=c(top.genes,pi3k.genes),drugs=T)
write.table(my.dgi[my.dgi$interaction_types!="",.(Gene,interaction_types,drug_name,drug_claim_name)],file="AD_panel.both.loose.dgi.tsv",sep="\t",quote=F,row.names=F)

input=read.delim("mixing_experiment.summary.tsv",header=F,stringsAsFactors=F)
colnames(input)=c("ID","Num_detected","Mean_detected","SD_detected","Num_all","Mean_all","SD_all")
input$ID=factor(input$ID,levels=c("AF0.2","AF0.5","AF1","AF2","AF5","AF10"))
input$MAF=c(0.002,0.005,0.01,0.02,0.05,0.1)
input$Prop=input$Num_detected/input$Num_all

ready=data.frame(MAF=c(0,input$MAF,1),Prop=c(0,input$Prop,1))

AD.top.maf=subsetMaf(maf=AD.maf,genes=top.genes,mafObj=T)
control.top.maf=subsetMaf(maf=control.maf,genes=top.genes,mafObj=T)
harmonic.mean(my.sensitivity(as.numeric(AD.top.maf@data$MAF),ready))
harmonic.mean(my.sensitivity(as.numeric(control.top.maf@data$MAF),ready))

AD_multi_num=17
control_multi_num=0
fisher.test(matrix(c(AD_multi_num,control_multi_num,AD_num-AD_multi_num,control_num-control_multi_num),ncol=2))

pdf("AD_panel.both.loose.multigene.pdf",width=5,height=5)

ready=data.frame(Cogdx=c("AD","Control","AD","Control"),Type=c("Multi","Multi","Single","Single"),Count=c(AD_multi_num,control_multi_num,AD_num-AD_multi_num,control_num-control_multi_num))
ready$Cogdx=factor(ready$Cogdx,levels=c("Control","AD"))
p <- ggplot(ready,aes(x=Cogdx,y=Count,fill=Type))
p + geom_bar(stat="identity",color="black",position="stack") + theme_classic() + theme(text=element_text(size=15)) + xlab("") + ylab("Number of samples")

dev.off()
