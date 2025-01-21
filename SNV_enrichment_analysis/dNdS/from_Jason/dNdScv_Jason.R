library(dndscv)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(gage)
library(gageData)

setwd("/Users/zhemingan/Documents/annovar/dNdS/from_Jason")
data(egSymb)
data(kegg.sets.hs)
data(sigmet.idx.hs)

AD_num=157
control_num=162

#sig.genes=c("TET2","ASXL1","KMT2D","ATRX","CBL")
sig.genes=c("CBL","DNMT3A","EZH2","WDR59")

my.clinical <- function()
{
        id_clinical=read.delim("AD_panel.sample.clinical.tsv",header=F,stringsAsFactors=F)
        colnames(id_clinical)=c("ID","ApoE","Age","Gender","MMSE","Cogdx")
        id_clinical$Gender[id_clinical$Gender=="F"]="Female"
        id_clinical$Gender[id_clinical$Gender=="M"]="Male"
        id_clinical$ApoE[id_clinical$ApoE=="ε3/ε3"]=33
        id_clinical$ApoE[id_clinical$ApoE=="ε3/ε4"]=34
        id_clinical$Age=as.numeric(id_clinical$Age)
        id_clinical$Gender=factor(id_clinical$Gender,levels=c("Female","Male"))
        id_clinical$ApoE=factor(id_clinical$ApoE,levels=sort(unique(id_clinical$ApoE)))
        id_clinical$Cogdx=factor(id_clinical$Cogdx,levels=c("CTRL","C","AD"))
        id_clinical=id_clinical[order(id_clinical$ID),]
        id_clinical
}

my.cancerCensus <- function()
{
	raw_list=read.delim("CancerGeneCensus_01212022.tsv",header=T,stringsAsFactors=F)
	final_list=data.frame(Symbol=raw_list$Gene.Symbol,Tier=raw_list$Tier)
	final_list$Hallmark=raw_list$Hallmark=="Yes"
	final_list$Role=sapply(1:nrow(raw_list),function(x){if(grepl("oncogene, TSG",raw_list$Role.in.Cancer[x])){return("both")}else if(grepl("oncogene",raw_list$Role.in.Cancer[x])){return("oncogene")}else if(grepl("TSG",raw_list$Role.in.Cancer[x])){return("TSG")}else{return("other")}})
	final_list
}

# my.summary <- function(input_header,output_header) 
# {
input_header = "AD_panel.MH.loose"
output_header = "AD_panel.MH.loose"
	input=read.delim(sprintf("%s.tsv",input_header),header=F,stringsAsFactors=F)[,c(1:4,6:7,9,10)]
	colnames(input)=c("Chr","Pos","Ref","Alt","Ref_num","Alt_num","ID","Gnomad")
	input$MAF=input$Alt_num/(input$Ref_num+input$Alt_num)
	
	id_clinical=my.clinical()
	# exclude problematic samples
	#id_excluded=read.delim("AD_panel.excluded.sample.list",header=F,stringsAsFactors=F)
	#id_clinical=id_clinical[!id_clinical$ID %in% id_excluded$V1[id_excluded$V2!="cross-contamination"],]

	panel_design=read.delim("3277911_Covered.bed",header=F,stringsAsFactors=F)[-(1:2),]
	colnames(panel_design)=c("Chr","Start","End","Symbol")
	panel_design$Length=panel_design$End-panel_design$Start
	
	id_coverage=read.delim("AD_panel.coverage.summary",header=F)
	colnames(id_coverage)=c("ID","Depth","Proportion")
	# the proportion of bases covered by 500 or more deduped reads
	id_coverage$Coverage=round(id_coverage$Proportion*sum(panel_design$Length))
	
	merged=unique(merge(merge(input,id_clinical,by="ID"),id_coverage,by="ID"))
	
	cancer_census=my.cancerCensus()
	
	# identify protential contamination samples
	count=merge(as.data.frame(table(merged$ID)),as.data.frame(table(merged$ID[which(merged$Gnomad<1e-4)])),by="Var1",all=T)
	id_contamination=count[which(count[,2]>=10&(is.na(count[,3])|count[,3]/count[,2]<=0.1)),1]
	
	#merged=merged[which(!merged$ID %in% id_contamination),]
	merged=merged[which(merged$Ref_num+merged$Alt_num>=100),]
	merged=merged[which(merged$Gnomad<1e-4),]
	
	AD_list=merged[which(merged$Cogdx=="AD"),1:5]
	# AD_list <- AD_list[AD_list$ID %in% unique(AD_list$ID)[1:15], ]
	control_list=merged[which(merged$Cogdx=="CTRL"),1:5]
	# control_list <- control_list[control_list$ID %in% unique(control_list$ID)[1:15], ]
	
	panel_gene=unique(panel_design$Symbol)
	panel_gene=panel_gene[!grepl("chr",panel_gene)]
	
	AD_dnds=dndscv(AD_list,max_muts_per_gene_per_sample=Inf,max_coding_muts_per_sample=Inf,outmats=T)
	control_dnds=dndscv(control_list,max_muts_per_gene_per_sample=Inf,max_coding_muts_per_sample=Inf,outmats=T)
	
	# AD_dnds=dndscv(AD_list,gene_list=panel_gene,max_muts_per_gene_per_sample=Inf,max_coding_muts_per_sample=Inf,outmats=T)
	# control_dnds=dndscv(control_list,gene_list=panel_gene,max_muts_per_gene_per_sample=Inf,max_coding_muts_per_sample=Inf,outmats=T)
	
	final=rbind(data.frame(AD_dnds$sel_cv[,c(1:8,11,14)],Cogdx="AD"),data.frame(control_dnds$sel_cv[,c(1:8,11,14)],Cogdx="CTRL"))
	colnames(final)=c("Gene","N_synonymous","N_missense","N_nonsense","N_splicing","Ratio_missense","Ratio_nonsense","Ratio_splicing","Pvalue","Qvalue","Cogdx")
	final$Gene=factor(final$Gene,levels=rev(AD_dnds$sel_cv[,1]))
	final$Cogdx=factor(final$Cogdx,levels=c("CTRL","AD"))
	final$Pconvert=-log10(final$Pvalue)
	final$Qconvert=-log10(final$Qvalue)
	final$Genetype="Other"
	final$Genetype[final$Gene %in% cancer_census$Symbol[cancer_census$Role=="TSG"]]="TSG"
	final$Genetype[final$Gene %in% sig.genes]="Hotspot"
	final$Genetype=factor(final$Genetype,levels=c("Hotspot","TSG","Other"))
	write.table(final[final$Gene %in% final$Gene[final$Pvalue<0.05],],file=sprintf("%s.dNdScv.gene.tsv",output_header),quote=F,sep="\t",row.names=F)
	
	pdf(sprintf("%s.dNdScv.1.pdf",output_header),width=8,height=3)
	p <- ggplot(final,aes(x=Gene,y=Pconvert,fill=Cogdx,color=Cogdx))
	print(p + geom_point() + geom_hline(yintercept=-log10(0.05),linetype="dashed") + geom_text_repel(data=final[final$Pvalue<0.05,],aes(label=Gene),max.overlaps=Inf) + scale_fill_manual(values=c("royalblue1","indianred3")) + scale_color_manual(values=c("royalblue1","indianred3")) + theme_classic() + theme(text=element_text(size=15),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ylab("-log10(P-value)"))
	p <- ggplot(final,aes(x=Gene,y=Qconvert,fill=Cogdx,color=Cogdx))
	print(p + geom_point() + geom_hline(yintercept=-log10(0.05),linetype="dashed") + geom_text_repel(data=final[final$Qvalue<0.05,],aes(label=Gene),max.overlaps=Inf) + scale_fill_manual(values=c("royalblue1","indianred3")) + scale_color_manual(values=c("royalblue1","indianred3")) + theme_classic() + theme(text=element_text(size=15),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ylab("-log10(adjusted P-value)"))
	dev.off()

	final2=rbind(data.frame(AD_dnds$globaldnds,Cogdx="AD"),data.frame(control_dnds$globaldnds,Cogdx="CTRL"))
	colnames(final2)=c("Type","MLE","LowCI","HighCI","Cogdx")
	final2$Type=factor(c("Missense","Nonsense","Splicing","Truncating","All"),levels=c("All","Missense","Nonsense","Splicing","Truncating"))
	final2$Cogdx=factor(final2$Cogdx,levels=c("CTRL","AD"))
	
	pdf(sprintf("%s.dNdScv.2.pdf",output_header),width=5.5,height=2.5)
	p <- ggplot(final2[final2$Type!="Truncating",],aes(x=Type,y=MLE,fill=Cogdx,color=Cogdx))
	print(p + geom_pointrange(aes(ymin=LowCI,ymax=HighCI),position=position_dodge(width=0.5)) + geom_hline(yintercept=1,linetype="dashed") + scale_fill_manual(values=c("royalblue1","indianred3")) + scale_color_manual(values=c("royalblue1","indianred3")) + theme_classic() + theme(text=element_text(size=15)) + xlab("Mutation type") + ylab("dN/dS ratio"))
	dev.off()
	
	AD_gene=geneci(AD_dnds,gene_list=unique(final$Gene[final$Qvalue<0.05]))
	control_gene=geneci(control_dnds,gene_list=unique(final$Gene[final$Qvalue<0.05]))
	missense=rbind(data.frame(AD_gene[,c(1,2,4,6)],Cogdx="AD",Type="Missense"),data.frame(control_gene[,c(1,2,4,6)],Cogdx="CTRL",Type="Missense"))
	colnames(missense)=c("Gene","MLE","LowCI","HighCI","Cogdx","Type")
	truncating=rbind(data.frame(AD_gene[,c(1,3,5,7)],Cogdx="AD",Type="Truncating"),data.frame(control_gene[,c(1,3,5,7)],Cogdx="CTRL",Type="Truncating"))
	colnames(truncating)=c("Gene","MLE","LowCI","HighCI","Cogdx","Type")
	final3=rbind(missense,truncating)
	final3$Type=factor(final3$Type,levels=c("Missense","Truncating"))
	final3$Cogdx=factor(final3$Cogdx,levels=c("CTRL","AD"))
	
	pdf(sprintf("%s.dNdScv.3.pdf",output_header),width=8,height=2.5)
	p <- ggplot(final3,aes(x=Type,y=MLE,fill=Cogdx,color=Cogdx))
	print(p + geom_pointrange(aes(ymin=LowCI,ymax=HighCI),position=position_dodge(width=0.5)) + geom_hline(yintercept=1,linetype="dashed") + scale_fill_manual(values=c("royalblue1","indianred3")) + scale_color_manual(values=c("royalblue1","indianred3")) + theme_classic() + theme(text=element_text(size=15)) + xlab("Mutation type") + ylab("dN/dS ratio") + scale_y_log10() + facet_grid(. ~ Gene))
	dev.off()
	
	mut_anno=unique(rbind(AD_dnds$annotmuts,control_dnds$annotmuts)[,c(2:6,15)])
	colnames(mut_anno)=c("Chr","Pos","Ref","Alt","Gene","Type")
	mut_merged=merge(merged,mut_anno,by=c("Chr","Pos","Ref","Alt"))
	
	AD_gene_sum=data.frame()
	control_gene_sum=data.frame()
	for(i in panel_gene)
	{
		AD_gene_sum=rbind(AD_gene_sum,unlist(c(AD_dnds$sel_cv[AD_dnds$sel_cv$gene_name==i,3:7],mean(mut_merged$MAF[mut_merged$Gene==i&mut_merged$Type=="Missense"&mut_merged$Cogdx=="AD"]),mean(mut_merged$MAF[mut_merged$Gene==i&mut_merged$Type %in% c("Nonsense","Essential_Splice")&mut_merged$Cogdx=="AD"]))))
		control_gene_sum=rbind(control_gene_sum,unlist(c(control_dnds$sel_cv[control_dnds$sel_cv$gene_name==i,3:7],mean(mut_merged$MAF[mut_merged$Gene==i&mut_merged$Type=="Missense"&mut_merged$Cogdx=="CTRL"]),mean(mut_merged$MAF[mut_merged$Gene==i&mut_merged$Type %in% c("Nonsense","Essential_Splice")&mut_merged$Cogdx=="CTRL"]))))
	}
	colnames(AD_gene_sum)=c("N_missense","N_nonsense","N_splicing","Ratio_missense","Ratio_truncating","MAF_missense","MAF_truncating")
	colnames(control_gene_sum)=c("N_missense","N_nonsense","N_splicing","Ratio_missense","Ratio_truncating","MAF_missense","MAF_truncating")
	
	AD_gene_sum$Gene=panel_gene
	AD_gene_sum[is.na(AD_gene_sum)]=0
	AD_gene_sum$Ratio_missense[AD_gene_sum$Ratio_missense<1]=1
	AD_gene_sum$Ratio_truncating[AD_gene_sum$Ratio_truncating<1]=1
	AD_gene_sum$Prop_missense=AD_gene_sum$N_missense*(AD_gene_sum$Ratio_missense-1)/AD_gene_sum$Ratio_missense*AD_gene_sum$MAF_missense*2
	AD_gene_sum$Prop_truncating=(AD_gene_sum$N_nonsense+AD_gene_sum$N_splicing)*(AD_gene_sum$Ratio_truncating-1)/AD_gene_sum$Ratio_truncating*AD_gene_sum$MAF_truncating*2
	AD_gene_sum$Prop_all=AD_gene_sum$Prop_missense+AD_gene_sum$Prop_truncating
	
	control_gene_sum$Gene=panel_gene
	control_gene_sum[is.na(control_gene_sum)]=0
	control_gene_sum$Ratio_missense[control_gene_sum$Ratio_missense<1]=1
	control_gene_sum$Ratio_truncating[control_gene_sum$Ratio_truncating<1]=1
	control_gene_sum$Prop_missense=control_gene_sum$N_missense*(control_gene_sum$Ratio_missense-1)/control_gene_sum$Ratio_missense*control_gene_sum$MAF_missense*2
	control_gene_sum$Prop_truncating=(control_gene_sum$N_nonsense+control_gene_sum$N_splicing)*(control_gene_sum$Ratio_truncating-1)/control_gene_sum$Ratio_truncating*control_gene_sum$MAF_truncating*2
	control_gene_sum$Prop_all=control_gene_sum$Prop_missense+control_gene_sum$Prop_truncating
	
	final4=rbind(data.frame(Prop=sum(AD_gene_sum$Prop_all)/AD_num,Cogdx="AD",Gene_list="All_genes"),
				 data.frame(Prop=sum(AD_gene_sum$Prop_all[AD_gene_sum$Gene %in% sig.genes])/AD_num,Cogdx="AD",Gene_list="Hotspot_genes"),
				 data.frame(Prop=sum(control_gene_sum$Prop_all)/control_num,Cogdx="CTRL",Gene_list="All_genes"),
				 data.frame(Prop=sum(control_gene_sum$Prop_all[control_gene_sum$Gene %in% sig.genes])/control_num,Cogdx="CTRL",Gene_list="Hotspot_genes"))
	final4$Prop_normalized=final4$Prop/final4$Prop[final4$Cogdx=="CTRL"]
	final4$Cogdx=factor(final4$Cogdx,levels=c("CTRL","AD"))
	
	pdf(sprintf("%s.dNdScv.4.pdf",output_header),width=5,height=3)
	p <- ggplot(final4,aes(x=Cogdx,y=Prop_normalized,fill=Cogdx))
	print(p + geom_col(color="black",width=0.7,position=position_dodge()) + scale_fill_manual(values=c("royalblue1","indianred3")) + theme_classic() + theme(text=element_text(size=15)) + xlab("") + ylab("Normalized number of\npositively selected cells") + facet_grid(. ~ Gene_list))
	dev.off()
	
	kegg.gs=kegg.sets.hs[sigmet.idx.hs]
	kegg.gs.sym=lapply(kegg.gs,eg2sym)
	
	AD_pathway_sum=numeric(0)
	control_pathway_sum=numeric(0)
	pathway_name=numeric(0)
	pathway_gene_count=numeric(0)
	for(i in names(kegg.gs.sym))
	{
		a=intersect(unlist(kegg.gs.sym[i]),panel_gene)
		if(length(a)>1)
		{
			if("CDKN2A" %in% a)
			{
				a=c(setdiff(a,"CDKN2A"),"CDKN2A.p16INK4a","CDKN2A.p14arf")
			}
			AD_pathway_sum=rbind(AD_pathway_sum,genesetdnds(AD_dnds,gene_list=a)$globaldnds_geneset[5,])
			control_pathway_sum=rbind(control_pathway_sum,genesetdnds(control_dnds,gene_list=a)$globaldnds_geneset[5,])
			pathway_name=c(pathway_name,i)
			pathway_gene_count=c(pathway_gene_count,length(a))
		}
	}
	AD_pathway_sum=data.frame(AD_pathway_sum)
	colnames(AD_pathway_sum)=c("MLE","LowCI","HighCI")
	AD_pathway_sum$Pathway=pathway_name
	AD_pathway_sum$Count=pathway_gene_count
	AD_pathway_sum$Cogdx="AD"
	control_pathway_sum=data.frame(control_pathway_sum)
	colnames(control_pathway_sum)=c("MLE","LowCI","HighCI")
	control_pathway_sum$Pathway=pathway_name
	control_pathway_sum$Count=pathway_gene_count
	control_pathway_sum$Cogdx="CTRL"
	
	final5=rbind(AD_pathway_sum,control_pathway_sum)
	write.table(final5[final5$Pathway %in% unique(final5$Pathway[final5$MLE>1&final5$LowCI>1]),],file=sprintf("%s.dNdScv.pathway.tsv",output_header),quote=F,sep="\t",row.names=F)
}

# my.summary("AD_panel.MH.stringent","AD_panel.MH.stringent")
my.summary("AD_panel.MH.loose","AD_panel.MH.loose")
