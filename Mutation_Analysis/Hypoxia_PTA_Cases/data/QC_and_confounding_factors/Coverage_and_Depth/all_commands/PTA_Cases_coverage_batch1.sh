#!/bin/bash
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH -o casejob_genomic_coverage1_%j.out
#SBATCH -e casejob_genomic_coverage1_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/9.2.0
module load samtools/1.15.1

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/604_B2.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_604_B2_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/604_B3.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_604_B3_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/604_B6.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_604_B6_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1113_D1.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1113_D1_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1113_E1.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1113_E1_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1113_F1.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1113_F1_Cases.txt