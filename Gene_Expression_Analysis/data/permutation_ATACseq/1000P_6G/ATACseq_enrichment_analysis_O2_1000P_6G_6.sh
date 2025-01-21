#!/bin/bash
#SBATCH --mem=35G
#SBATCH -c 1
#SBATCH -t 0-3:00
#SBATCH -p short
#SBATCH -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Enrichment_Analysis/ATACseq_analysis/log/ATACseq_analysis_1000P_6G_6_%j.out
#SBATCH -e /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Enrichment_Analysis/ATACseq_analysis/log/ATACseq_analysis_1000P_6G_6_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/9.2.0
module load R/4.3.1

Rscript /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Enrichment_Analysis/ATACseq_analysis/commands/ATACseq_enrichment_by_Cell_O2.R 6 6