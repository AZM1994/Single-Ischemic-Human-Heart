#!/bin/bash
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH -o log/MAPD_5919_1_F6-2n_step_1_%j.out
#SBATCH -e log/MAPD_5919_1_F6-2n_step_1_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/6.2.0
module load R/3.6.1
module load java/jdk-1.8u112
module load bwa/0.7.15
module load picard/2.8.0
module load samtools/1.3.1
module load pigz/2.3.4
module load python/2.7.12
module load openblas/0.2.19

cd /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD
mkdir -p 5919_1_F6-2n && cd 5919_1_F6-2n

samtools view -h /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/5919_1_F6-2n.bam | python /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/varbin.50k.bam.py 5919_1_F6-2n.varbin.50k.txt 5919_1_F6-2n.varbin.50k.stats.txt /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/hg19.bin.boundaries.50k.bowtie.k50.sorted.txt /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/hg19.chrom.sizes.txt
