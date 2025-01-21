#!/bin/bash
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH -t 0-3:00
#SBATCH -p short
#SBATCH -o casejob_genomic_coverage_%j.out
#SBATCH -e casejob_genomic_coverage_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/9.2.0
module load samtools/1.15.1

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/6032_M-C7-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_6032_C7_Cases.txt