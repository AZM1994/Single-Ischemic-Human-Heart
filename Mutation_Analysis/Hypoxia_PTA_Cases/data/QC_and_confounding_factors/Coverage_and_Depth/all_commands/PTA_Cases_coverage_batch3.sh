#!/bin/bash
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH -o casejob_genomic_coverage3_%j.out
#SBATCH -e casejob_genomic_coverage3_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/9.2.0
module load samtools/1.15.1

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1039.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1039_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/1864_M-E3-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_1864_M-E3-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/4638_1_A1-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_4638_1_A1-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/4638_1_A4-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_4638_1_A4-2n_Cases.txt

samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/4638_1_B2-2n.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_4638_1_B2-2n_Cases.txt
