#!/bin/bash
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH -o casejob_genomic_coverage2_%j.out
#SBATCH -e casejob_genomic_coverage2_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/9.2.0
module load samtools/1.15.1

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1743_A3.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1743_A3_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1743_C3.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1743_C3_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1743_F3.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1743_F3_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1363_A4.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1363_A4_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1363_D4.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1363_D4_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1363_H4.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1363_H4_Cases.txt