# Rscript /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Scan2_R_Analysis/extract_gatk_table.R 604_B2 604_Hypoxia_PTA

args <- commandArgs(TRUE)
verbose  = TRUE
cell_id  = args[1]
case_dir = args[2]

library(scan2)
setwd("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Scan2_R_Analysis/new_mut_burden_table/")

load(sprintf("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/%s/call_mutations/%s/scan2_object.rda",
             case_dir,cell_id))

object <- results

# sink(file = sprintf("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Scan2_R_Analysis/scan2_output_rds/scan2_output_%s.rds",cell_id))
# results
# sink(file = NULL)

# write.csv(results@mutburden, sprintf("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Scan2_R_Analysis/all_burden_table/%s.csv",cell_id), row.names=TRUE)

get.gbp.by.genome <- function(object) {
  if (object@genome.string == 'hs37d5') {
    # 93 contigs includes unplaced; see http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics.
    total <- 3137161264
    chrx <- 155270560
    chry <- 59373566
    chrm <- 16571
    return((total - chrx - chry - chrm)*2 / 1e9) # = 5.845001134
  } else if (object@genome.string == 'hg38') {
    # 455 contigs; see http://genomewiki.ucsc.edu/index.php/Hg38_100-way_Genome_size_statistics
    total <- 3209286105
    chrx <- 156040895
    chry <- 57227415
    chrm <- 16569
    return((total - chrx - chry - chrm)*2 / 1e9) # = 5.992002452
  } else if (object@genome.string == 'mm10') {
    # 66 contigs; see http://genomewiki.ucsc.edu/index.php/Hg38_100-way_Genome_size_statistics
    total <- 2730871774
    chrx <- 171031299
    chry <- 91744698
    chrm <- 16299
    return((total - chrx - chry - chrm)*2 / 1e9) # = 4.936158956
  } else {
    warn(paste('gbp not yet implemented for genome', object@genome.string))
    warn("the mutation burden for this analysis is a placeholder!")
    warn("DO NOT USE!")
    # hopefully returning a negative number will alert people that something
    # has gone wrong so they don't ignore the warning messages above
    return(-1)
  }
}

gbp.per.genome=get.gbp.by.genome(object)

# setGeneric("compute.mutburden", function(object, gbp.per.genome=get.gbp.by.genome(object), quiet=FALSE)
#         standardGeneric("compute.mutburden"))
# setMethod("compute.mutburden", "SCAN2", function(object, gbp.per.genome=get.gbp.by.genome(object), quiet=FALSE) {
#     check.slots(object, c('call.mutations', 'depth.profile'))

calculation_mod <- c('scan2', 'revised')
muttypes <- c('snv', 'indel')

for (mod_type in calculation_mod){
for (mt in muttypes){
# mt <- 'snv'
# mt <- 'indel'
# object@mutburden <- setNames(lapply(muttypes, function(mt) {
# [2] is the maximum burden; the minimum burden [1] is almost always ~0
print(c('mod = ', mod_type, 'mut = ', mt))
pre.geno.burden <- object@fdr.prior.data[[mt]]$burden[2]
sfp <- object@static.filter.params[[mt]]

# germline sites
if (mod_type == 'scan2') {
  g <- object@gatk[resampled.training.site == TRUE & muttype == mt] 
} else if (mod_type == 'revised') {
  g <- object@gatk[resampled.training.site == TRUE & muttype == mt & dp != 0]
} else {
  print('wrong mod')
}

# (single cell  x  bulk) depth table
dptab <- object@depth.profile$dptab
dptab <- dptab[1:min(max(g$dp)+1, nrow(dptab)),]

# these computations rely on there being a reasonably large number
# of germline sites tested. even 100 is very few; we expect more like
# 100,000.
if (nrow(g) < 100) {
  warning(paste('only', nrow(g), 'resampled germline', mt, 'sites were detected; aborting genome-wide extrapolation. Typical whole-genome experiments include ~10-100,000 germline sites'))
  ret <- data.frame(
    ncalls=NA,
    callable.sens=NA,
    callable.bp=NA
  )[c(1,1,1),]  # repeat row 1 3 times
} else {
  # somatic sites
  s <- object@gatk[pass == TRUE & muttype == mt]
  
  # Break data into 4 quantiles based on depth, use the middle 2 (i.e.,
  # middle 50%) to reduce noise caused by very low and very high depth.
  q=4
  qstouse <- c(1,2,4)
  qbreaks <- quantile(g$dp, prob=0:q/q)
  # qbreaks[2] <- 1
  # qbreaks[3] <- 2
  
  if (length(unique(qbreaks)) != q+1) {
    # if (length(unique(qbreaks)) > 100) {
    ncalls <- rep(NA, 3)
    callable.sens <- rep(NA, 3)
    rowqs <- rep(NA, 3)
    callable.bp <- rep(NA, 3)
    warning(paste('could not derive unique breakpoints for quartiles of sequencing depth at germline hSNPs.  this usually indicates that sequencing depth is heavily skewed toward low depths (typically DP=0)\ngot qbreaks = ', deparse(qbreaks)))
  } else {
    # s also uses g-based depth quantiles
    s$dpq <- cut(s$dp, qbreaks, include.lowest=T, labels=F)
    s$dpq[s$dpq==3] <- 2 # merge 25-75% into a single bin
    g$dpq <- cut(g$dp, qbreaks, include.lowest=T, labels=F)
    g$dpq[g$dpq==3] <- 2
    
    # select the subset of the depth profile passing the bulk depth requirement
    # cut down dptab to the max value in g$dp (+1 because 1 corresponds to dp=0)
    rowqs <- cut(0:(nrow(dptab)-1), qbreaks, include.lowest=T, labels=F)
    rowqs[rowqs==3] <- 2
    
    s <- s[dpq %in% qstouse]
    g <- g[dpq %in% qstouse]
    
    ncalls <- sapply(qstouse, function(q) sum(s[dpq == q]$pass, na.rm=TRUE))
    callable.bp <- sapply(split(dptab[,-(1:sfp$min.bulk.dp)], rowqs), sum)
    callable.sens <- sapply(qstouse, function(q) mean(g[bulk.dp >= sfp$min.bulk.dp & dpq == q & resampled.training.site == TRUE]$training.pass, na.rm=TRUE))  # only use resampled training sites to estimate sensitivity (hSNPs in general are closer to other hSNPs than somatic candidates are to hSNPs, leading to better local allele balance estimates and thus would overestimate sensitivity. Resampling the sites to match candidate-to-hSNP distances reduces this bias)
  }
  
  # this data.frame has 1 row for each quantile. the second row (=middle 50%)
  # is ultimately what we're interested in, but having the other calculations
  # around can also be interesting.
  ret <- data.frame(ncalls=ncalls, callable.sens=callable.sens, callable.bp=callable.bp)
}

# "callable" means:
# Sensitivity estimates only from germline training sites with the same
# depth cutoffs as somatic candidates. Detailed depth tables will be used
# to ensure extrapolation to the rest of the genome is equitable.
ret$callable.burden <- ret$ncalls / ret$callable.sens
# dividing by 2 makes it haploid gb
ret$rate.per.gb <- ret$callable.burden / ret$callable.bp * 1e9/2
ret$burden <- ret$rate.per.gb * gbp.per.genome
ret$somatic.sens <- ret$ncalls / ret$burden
ret$pre.genotyping.burden <- pre.geno.burden

ret$unsupported.filters <- sfp$max.bulk.alt > 0 | sfp$max.bulk.af > 0
if (any(ret$unsupported.filters)) {
  warning('mutation burdens must be extrapolated without clonal mutations (i.e., max.bulk.alt=0 and max.bulk.af=0)! burdens will be estimated, but they are invalid')
}
# ret
# saveRDS(ret, sprintf("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Scan2_R_Analysis/new_mut_burden_table/%s_gatk_table_%s.rds",mt,cell_id))
write.csv(ret, sprintf("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Scan2_R_Analysis/new_mut_burden_table/%s/%s_gatk_table_%s_%s.csv",mod_type,mt,cell_id,mod_type))
}
}
