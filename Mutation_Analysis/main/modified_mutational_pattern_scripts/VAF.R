# Rscript /n/no_backup2/bch/lee/zheming/heart_hypoxia/scan2_analysis/Scan2_R_Analysis/Extract_Scan2_Result.R 604_B2 604_Hypoxia_PTA
############# read scan2 raw results
# args <- c("604_B2", "604_Hypoxia_PTA")
args <- commandArgs(TRUE)
verbose  = TRUE
cell_id  = args[1]
case_dir = args[2]

wd <- setwd("/n/no_backup2/bch/lee/zheming/heart_hypoxia/scan2_analysis/")
library(scan2)
load(sprintf("%s/call_mutations/%s/scan2_object.rda",case_dir,cell_id))

####################################### extract snv results
#######################################
snv_summary_df = data.frame(results@call.mutations$snv.pass, results@mutburden$snv[2,], cell_id)
if(!is.na(snv_summary_df$burden)&snv_summary_df$burden > 10^8) {
	snv_summary_df$burden = snv_summary_df$burden*2/10^9
	snv_summary_df$somatic.sens = snv_summary_df$somatic.sens/2*10^9
}

snv_df = results@gatk
# snv_df = snv_df[snv_df$pass == TRUE&snv_df$muttype == "snv",]

############# define and write vcf file header
# new_header <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
# vcf_head <- paste0("##fileformat=VCFv4.1\n#",
                   # paste(new_header, collapse = "\t"))
# fileConn_ssnv <- file(sprintf("Scan2_R_Analysis/all_vcfs/ssnv_list.%s.vcf", cell_id))
# writeLines(vcf_head, fileConn_ssnv)
# close(fileConn_ssnv)

############# save snv vcf
if(nrow(snv_df) > 0) {
	snv_vcf_df = data.frame(snv_df$chr, snv_df$pos, '.', snv_df$refnt, snv_df$altnt, 
	                        snv_df$scref, snv_df$scalt, snv_df$muttype, snv_df$dp, snv_df$ab, 
	                        snv_df$pass, cell_id)
	write.table(snv_vcf_df, file = sprintf("Scan2_R_Analysis/all_vcfs/ssnv_list.%s.vcf", cell_id),
	            quote = F, sep = "\t", row.names = F, col.names = FALSE, append = T)
}

############# save snv summary vcf
# if(nrow(snv_summary_df)>0) {
# 	write.table(snv_summary_df, file = sprintf("Scan2_R_Analysis/all_vcfs/summary/ssnv_summary.%s.vcf",cell_id),
# 	            quote = F, sep = "\t", row.names = F)
# }

# ####################################### extract indel results
# #######################################
# indel_summary_df = data.frame(results@call.mutations$indel.pass, results@mutburden$indel[2,], cell_id)
# if(!is.na(indel_summary_df$burden)&indel_summary_df$burden > 10^8) {
# 	indel_summary_df$burden = indel_summary_df$burden*2/10^9
# 	indel_summary_df$somatic.sens = indel_summary_df$somatic.sens/2*10^9
# }
# 
# indel_df = results@gatk
# indel_df = indel_df[indel_df$pass == TRUE&indel_df$muttype == "indel",]
# 
# ############# define and write vcf file header
# fileConn_sindel <- file(sprintf("Scan2_R_Analysis/all_vcfs/sindel_list.%s.vcf", cell_id))
# writeLines(vcf_head, fileConn_sindel)
# close(fileConn_sindel)
# 
# ############# save indel vcf
# if(nrow(indel_df) > 0) {
#   indel_vcf_df = data.frame(indel_df$chr, indel_df$pos, '.', indel_df$refnt, indel_df$altnt, 
#                             indel_df$scref, indel_df$scalt, cell_id)
# 	write.table(indel_vcf_df,file = sprintf("Scan2_R_Analysis/all_vcfs/sindel_list.%s.vcf", cell_id),
# 	            quote = F, sep = "\t", row.names = F, col.names = FALSE, append = T)
# }
# 
# ############# save indel summary vcf
# if(nrow(indel_summary_df) > 0) {
# 	write.table(indel_summary_df, file = sprintf("Scan2_R_Analysis/all_vcfs/summary/sindel_summary.%s.vcf", cell_id),
# 	            quote = F, sep = "\t", row.names = F)
# }