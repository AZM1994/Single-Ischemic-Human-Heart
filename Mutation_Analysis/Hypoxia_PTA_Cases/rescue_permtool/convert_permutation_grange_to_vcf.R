library(stringr)

arg <- commandArgs(trailingOnly = TRUE)
# case_ID <- arg[1]
path <- arg[1]
# snv permutation results
snv_permutation_data_path <- paste0(path, "/perms_snv_pass.rda")
load(snv_permutation_data_path)

# convert GRange to vcf
mutation_type <- str_extract(as.character(zperml@unlistData$mutsig), pattern = "[ATCG]>[ATCG]")
vcf <- data.frame("CHROM" = paste0("chr", as.character(zperml@unlistData@seqnames)),
                         "POS" = zperml@unlistData@ranges@start,
                         "ID" = ".",
                         "REF" = str_sub(mutation_type, 1, 1),
                         "ALT" = str_sub(mutation_type, 3, 3),
                         "QUAL" = ".",
                         "FILTER" = "PASS",
                         "INFO" = zperml@unlistData$perm.id)
colnames(vcf)[1] <- "#CHROM"

rm(zperml, mutation_type)

# save vcf file
vcf_file <- paste0(path, "/permutation_snv.vcf")
write.table(vcf, vcf_file, quote = F, row.names = F, sep = "\t")
system(paste("sed -i '1 i\\##fileformat=VCF'", vcf_file))

# annotation
# working_dir <- paste0(path, "/permtool_", case_ID)
working_dir <- path
convert2annovar_path <- "/n/data1/bch/genetics/lee/shulin/softwares/annovar/convert2annovar.pl"
system(paste0(convert2annovar_path, " -format vcf4 ", vcf_file, " -includeinfo > ", working_dir, "/permutation_snv.avinput"))

annotate_variation_path <- "/n/data1/bch/genetics/lee/shulin/softwares/annovar/annotate_variation.pl"
system(paste0(annotate_variation_path, " -build hg19 ", working_dir, "/permutation_snv.avinput /n/data1/bch/genetics/lee/shulin/softwares/annovar/humandb/"))


# indel permutation results
# indel_permutation_data_path <- paste0(path, "/permtool_", case_ID, "/perms_indel_pass.rda")
indel_permutation_data_path <- paste0(path, "/perms_indel_pass.rda")
load(indel_permutation_data_path)

# convert GRange to vcf
vcf <- data.frame("CHROM" = paste0("chr", as.character(zperml@unlistData@seqnames)),
                  "POS" = zperml@unlistData@ranges@start,
                  "ID" = ".",
                  "REF" = ".",
                  "ALT" = ".",
                  "QUAL" = ".",
                  "FILTER" = "PASS",
                  "INFO" = zperml@unlistData$perm.id)
colnames(vcf)[1] <- "#CHROM"

rm(zperml)

# save vcf file
vcf_file <- paste0(path, "/permutation_indel.vcf")
write.table(vcf, vcf_file, quote = F, row.names = F, sep = "\t")
system(paste("sed -i '1 i\\##fileformat=VCF'", vcf_file))

# annotation
# working_dir <- paste0(path, "/permtool_", case_ID)
working_dir <- path
convert2annovar_path <- "/n/data1/bch/genetics/lee/shulin/softwares/annovar/convert2annovar.pl"
system(paste0(convert2annovar_path, " -format vcf4 ", vcf_file, " -includeinfo > ", working_dir, "/permutation_indel.avinput"))

annotate_variation_path <- "/n/data1/bch/genetics/lee/shulin/softwares/annovar/annotate_variation.pl"
system(paste0(annotate_variation_path, " -build hg19 ", working_dir, "/permutation_indel.avinput /n/data1/bch/genetics/lee/shulin/softwares/annovar/humandb/"))
