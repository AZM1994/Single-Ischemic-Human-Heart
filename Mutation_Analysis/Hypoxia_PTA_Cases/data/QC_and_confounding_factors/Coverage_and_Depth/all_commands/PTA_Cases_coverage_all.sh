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

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/604_B2.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_604_B2_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/604_B3.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_604_B3_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/604_B6.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_604_B6_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1113_D1.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1113_D1_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1113_E1.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1113_E1_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1113_F1.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1113_F1_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1743_A3.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1743_A3_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1743_C3.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1743_C3_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1743_F3.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1743_F3_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1363_A4.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1363_A4_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1363_D4.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1363_D4_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1363_H4.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1363_H4_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1039.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1039_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1864_M-E3-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1864_M-E3-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/4638_1_A1-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_4638_1_A1-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/4638_1_A4-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_4638_1_A4-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/4638_1_B2-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_4638_1_B2-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/5657_M-D1-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_5657_M-D1-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/5657_M-H1-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_5657_M-H1-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/5919_1_C4-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_5919_1_C4-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/5919_1_D2-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_5919_1_D2-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/5919_1_E3-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_5919_1_E3-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/5919_1_F6-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_5919_1_F6-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/6032_1_A1-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_6032_1_A1-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/6032_1_C7-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_6032_1_C7-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/6032_1_E7-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_6032_1_E7-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/4402_1_A1-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_4402_1_A1-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/4402_1_A3-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_4402_1_A3-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1673_CM_A2_2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1673_CM_A2_2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1673_CM_A3_2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1673_CM_A3_2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1673_CM_D2_2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1673_CM_D2_2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/5828_CM_C2_2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_5828_CM_C2_2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/5828_CM_G2_2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_5828_CM_G2_2n_Cases.txt




