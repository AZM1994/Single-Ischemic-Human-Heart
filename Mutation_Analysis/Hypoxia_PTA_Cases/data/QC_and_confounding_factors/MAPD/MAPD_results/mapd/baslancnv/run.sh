


python ~/bin/baslancnv/varbin.50k.sam.py a.sam a.varbin.50k.txt a.varbin.50k.stats.txt hg19.bin.boundaries.50k.bowtie.k50.sorted.txt
Rscript ~/bin/baslancnv/cbs.r a	~/bin/baslancnv/hg19.50k.k50.bad.bins.txt ~/bin/baslancnv/hg19.varbin.gc.content.50k.bowtie.k50.txt
Rscript ~/bin/baslancnv/copynumber.r a
Rscript ~/bin/baslancnv/cal_mapd.r a

python ~/bin/baslancnv/varbin.50k.sam.py a.sam a.varbin.6k.txt a.varbin.6k.stats.txt ~/bin/baslancnv/hg19.bin.boundaries.6k.bowtie.k50.sorted.txt ~/bin/baslancnv/hg19.chrom.sizes.txt
Rscript ~/bin/baslancnv/cbs.r a 6k ~/bin/baslancnv/hg19.50k.k50.bad.bins.txt ~/bin/baslancnv/hg19.varbin.gc.content.6k.bowtie.k50.txt
Rscript ~/bin/baslancnv/copynumber.r a.hg19.6k.k50.varbin
Rscript ~/bin/baslancnv/cal_mapd.r a.hg19.6k.k50.varbin.data.copynumber.txt

python ~/bin/baslancnv/varbin.50k.sam.py a.sam a.varbin.50k.txt a.varbin.50k.stats.txt ~/bin/baslancnv/hg19.bin.boundaries.50k.bowtie.k50.sorted.txt ~/bin/baslancnv/hg19.chrom.sizes.txt
Rscript ~/bin/baslancnv/cbs.r a 50k ~/bin/baslancnv/hg19.50k.k50.bad.bins.txt ~/bin/baslancnv/hg19.varbin.gc.content.50k.bowtie.k50.txt
Rscript ~/bin/baslancnv/copynumber.r a.hg19.50k.k50.varbin
Rscript ~/bin/baslancnv/cal_mapd.r a.hg19.50k.k50.varbin.data.copynumber.txt

### for bad bin
Rscript ~/bin/baslancnv/copynumber.r a.hg19.50k.k50.nobad.varbin
Rscript ~/bin/baslancnv/cal_mapd.r a.hg19.50k.k50.nobad.varbin.data.copynumber.txt



Rscript ~/bin/baslancnv/cbs.r S59 ~/bin/baslancnv/hg19.50k.k50.bad.bins.txt ~/bin/baslancnv/hg19.varbin.gc.content.50k.bowtie.k50.txt
Rscript ~/bin/baslancnv/copynumber.r S59 > S59.cnv.log
Rscript ~/bin/baslancnv/cal_mapd.r S59