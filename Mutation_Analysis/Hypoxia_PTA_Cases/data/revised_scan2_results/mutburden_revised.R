load("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/5657_Normal_PTA/call_mutations/5657_M-D1-2n/scan2_object.rda")
load("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/6032_Normal_PTA/call_mutations/6032_M-C7-2n/scan2_object.rda")
setwd("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Scan2_R_Analysis")
library(scan2)
g <- results@gatk[resampled.training.site == TRUE & muttype == mt]

saveRDS(g, "5657_D1.rds")
saveRDS(g, "6032_C7.rds")


##### local
# wd_dir <- setwd(getwd())
# object <- readRDS(paste0(wd_dir, "/Hypoxia_PTA_Cases/check_sth/5657_D1.rds"))
# object <- readRDS(paste0(wd_dir, "/Hypoxia_PTA_Cases/check_sth/5657_H_muttable.rds"))
object <- readRDS("/Users/zhemingan/My Drive/BCH/Huang Lab/Heart_hypoxia/Mutation_Analysis/Hypoxia_PTA_Cases/check_sth/5657_D1.rds")
object <- readRDS("/Users/zhemingan/My Drive/BCH/Huang Lab/Heart_hypoxia/Mutation_Analysis/Hypoxia_PTA_Cases/check_sth/5919_C4.rds")
object <- readRDS("/Users/zhemingan/My Drive/BCH/Huang Lab/Heart_hypoxia/Mutation_Analysis/Hypoxia_PTA_Cases/check_sth/6032_C7.rds")

g <- object[resampled.training.site == TRUE & muttype == 'snv']
# 
q=4
gdp <- g$dp
qbreaks <- quantile(gdp, prob=0:q/q)
ggplot(g, aes(x=dp)) + geom_histogram()
# 
# 
qbreaks <- quantile(g$dp, prob=0:q/q)
length(unique(qbreaks)) != q+1
# nrow(g[g$dp == 0])

g2 <- object[resampled.training.site == TRUE & muttype == 'snv' & dp != 0]
# 
q=4
gdp <- g2$dp
qbreaks <- quantile(gdp, prob=0:q/q)
ggplot(g2, aes(x=dp)) + geom_histogram()
length(unique(qbreaks)) != q+1

# Number of haploid basepairs (in billions) per genome. The default
# value of 5.845001134 corresponds to AUTOSOMES as determined by GRCh37
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


setGeneric("compute.mutburden", function(object, gbp.per.genome=get.gbp.by.genome(object), quiet=FALSE)
        standardGeneric("compute.mutburden"))
setMethod("compute.mutburden", "SCAN2", function(object, gbp.per.genome=get.gbp.by.genome(object), quiet=FALSE) {
    check.slots(object, c('call.mutations', 'depth.profile'))

    muttypes <- c('snv', 'indel')
    object@mutburden <- setNames(lapply(muttypes, function(mt) {
        # [2] is the maximum burden; the minimum burden [1] is almost always ~0
        pre.geno.burden <- object@fdr.prior.data[[mt]]$burden[2]
        sfp <- object@static.filter.params[[mt]]
    
        # germline sites
        g <- object@gatk[resampled.training.site == TRUE & muttype == mt]
        # g <- object@gatk[resampled.training.site == TRUE & muttype == mt & dp != 0]

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
    
            if (length(unique(qbreaks)) != q+1) {
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
        ret
    }), muttypes)

    object
})


setGeneric("mutburden", function(object, muttype=c('all', 'snv', 'indel'))
    standardGeneric("mutburden"))
setMethod("mutburden", "SCAN2", function(object, muttype=c('all', 'snv', 'indel')) {
    check.slots(object, 'mutburden')

    muttype <- match.arg(muttype)
    if (muttype == 'all')
        muttype <- c('snv', 'indel')

    # We use the middle 50% of the depth distribution to avoid biasing
    # the mutation burden by sensitivity differences in very low or very high
    # depth regions. Row 2 corresponds to the middle 50%; row 1 is the bottom
    # 25% and row 3 is the top 25%.
    sapply(muttype, function(mt) {
        ret <- object@mutburden[[mt]][2,]$burden

        if (object@mutburden[[mt]]$unsupported.filters[2]) {
            warning('unsupported static.filter.params were used, invalidating mutburden estimation (this is currently true for allowing mutation supporting reads in bulk)')
            ret <- NA
        }
        if (mt == 'indel' & (object@call.mutations$suppress.all.indels | object@call.mutations$suppress.shared.indels)) {
            warning("indel mutation burden estimates ARE NOT VALID (cross-sample panel insufficiency)!")
            ret <- NA
        }
        ret
    })
})
