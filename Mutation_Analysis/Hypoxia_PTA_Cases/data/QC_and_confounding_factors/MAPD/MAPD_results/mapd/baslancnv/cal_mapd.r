library(ggplot2)
# library(DescTools)
library(ineq)
args = (commandArgs(TRUE));
pathfname = args[1];

cal_mapd <- function(ncnr){
	ncnr_1<-c(0,ncnr)
	n<-length(ncnr)
	b <- abs(log2(ncnr) - log2(ncnr_1[1:n]))
	#return (median(b[2:n]))
	mapd<-median(b[2:n])
	mean_mapd<-mean(b[2:n])
	std_mapd<-sqrt(var(b[2:n]))
	iqr<-IQR(b[2:n])
	result<-c(mapd,mean_mapd,std_mapd,iqr)
	return (result)
}

cal_stat_by_normBinCnt <- function(ncnr){
	n<-length(ncnr)
	mean_binReadCnt<-mean(ncnr)
	std_binReadCnt<-sqrt(var(ncnr))
	iqr<-IQR(ncnr)
	q25<-quantile(ncnr,0.25)
	q75<-quantile(ncnr,0.75)
	binFrac_0.1<-sum(ncnr<0.1)/n
	outlier_upper<-sum(ncnr > quantile(ncnr, 0.75) + 1.5 * IQR(ncnr))
	outlier_lower<-sum(ncnr < quantile(ncnr, 0.25) - 1.5 * IQR(ncnr))
	outlier_total<-sum(ncnr > quantile(ncnr, 0.75) + 1.5 * IQR(ncnr) | ncnr < quantile(ncnr, 0.25) - 1.5 * IQR(ncnr))
	gini<-Gini(ncnr)
	cov<-std_binReadCnt/mean_binReadCnt
	result<-c(mean_binReadCnt,std_binReadCnt,iqr,q25,q75,binFrac_0.1,outlier_total,outlier_upper,outlier_lower,gini,cov)
	return (result)
}

cal_mapd_plot <- function(ncnr,pathfname){
	ncnr_1<-c(0,ncnr)
	n<-length(ncnr)
	b <- abs(log2(ncnr) - log2(ncnr_1[1:n]))
	PD_raw <- log2(ncnr) - log2(ncnr_1[1:n])
	#return (median(b[2:n]))
	mapd<-median(b[2:n])
	mean_mapd<-mean(b[2:n])
	std_mapd<-sqrt(var(b[2:n]))
	iqr<-IQR(b[2:n])
	result<-c(mapd,mean_mapd,std_mapd,iqr)
	PD<-sort(PD_raw)
	df_pd<-data.frame(PD)
	bin<-1:n
	df_pd$bin<-bin
	p<-ggplot(df_pd, aes(x=bin,y=PD))+geom_bar(stat="identity",color="gray")+coord_cartesian(ylim=c(-10,10))+theme_classic()
	ggsave(file=gsub(".data.copynumber.txt", ".pd.dist.pdf",pathfname),p)
	return (result)
}

cal_stat_by_normBinCnt_plot <- function(ncnr,pathfname){
	n<-length(ncnr)
	mean_binReadCnt<-mean(ncnr)
	std_binReadCnt<-sqrt(var(ncnr))
	iqr<-IQR(ncnr)
	q25<-quantile(ncnr,0.25)
	q75<-quantile(ncnr,0.75)
	binFrac_0.1<-sum(ncnr<0.1)/n
	outlier_upper<-sum(ncnr > quantile(ncnr, 0.75) + 1.5 * IQR(ncnr))
	outlier_lower<-sum(ncnr < quantile(ncnr, 0.25) - 1.5 * IQR(ncnr))
	outlier_total<-sum(ncnr > quantile(ncnr, 0.75) + 1.5 * IQR(ncnr) | ncnr < quantile(ncnr, 0.25) - 1.5 * IQR(ncnr))
	gini<-Gini(ncnr)
	cov<-std_binReadCnt/mean_binReadCnt
	result<-c(mean_binReadCnt,std_binReadCnt,iqr,q25,q75,binFrac_0.1,outlier_total,outlier_upper,outlier_lower,gini,cov)
	normBinReadCnt<-sort(ncnr,decreasing=TRUE)
	df_bin<-data.frame(normBinReadCnt)
	bin<-1:n
	df_bin$bin<-bin
	#p<-ggplot(df_bin, aes(x=bin,y=normBinReadCnt))+geom_bar(stat="identity",color="gray")+coord_cartesian(ylim=c(0,100000))+theme_classic()
	p<-ggplot(df_bin, aes(x=bin,y=normBinReadCnt))+geom_bar(stat="identity",color="gray")+theme_classic()
	ggsave(file=gsub(".data.copynumber.txt", ".normBinCnt.dist.pdf",pathfname),p)

	return (result)
}




dat <- read.table(paste(pathfname,sep=""), header = T, sep="\t")
#cnr <- sort(log2(dat$cn.ratio))
#cnr <- log2(dat$cn.ratio)

chromlist <- c(1:22, "X", "Y")

cont <- c()
cont_normBin <- c()
#cont <- rbind(cont, c("chrom", "MAPD"))
cont <- rbind(cont, c("chrom", "MAPD","AVG_APD","STD_APD","IQR_APD"))
cont_normBin <- rbind(cont_normBin, c("chrom", "AVG_normBinCnt","STD_normBinCnt","IQR_normBinCnt","Q25","Q75","binFrac_0.1","outlierBinCnt_total","outlierBinCnt_upper","outlierBinCnt_lower","Gini","CoV"))
for (chrom in unique(dat$chrom)){
	idx = which(dat$chrom==chrom)
	mapd <- cal_mapd(dat[idx,]$cn.ratio)
	normBinStat <- cal_stat_by_normBinCnt(dat[idx,]$cn.ratio)
	#cat (mapd, "\n")
	cont <- rbind(cont, c(chromlist[chrom], mapd))
	cont_normBin <- rbind(cont_normBin, c(chromlist[chrom], normBinStat))
}
mapd <- cal_mapd_plot(dat$cn.ratio,pathfname)
normBinStat <- cal_stat_by_normBinCnt_plot(dat$cn.ratio,pathfname)
cont <- rbind(cont, c("Total",mapd))
cont_normBin <- rbind(cont_normBin, c("Total",normBinStat))

write.table(cont, file=gsub(".data.copynumber.txt", ".mapd.txt",pathfname), sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE, append=FALSE)
write.table(cont_normBin, file=gsub(".data.copynumber.txt", ".normBinStat.txt",pathfname), sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE, append=FALSE)


cat (mapd)


