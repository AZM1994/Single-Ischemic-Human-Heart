#!/bin/bash
#SBATCH --mem=50G
#SBATCH -c 1
#SBATCH -t 3-00:00
#SBATCH -p medium
#SBATCH -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Enrichment_Analysis/transcription_analysis/log/transcription_analysis_1000P_8G_3_%j.out
#SBATCH -e /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Enrichment_Analysis/transcription_analysis/log/transcription_analysis_1000P_8G_3_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/9.2.0
module load R/4.3.1

Rscript /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Enrichment_Analysis/transcription_analysis/working_generate_mutmat_by_Cell_O2.R 8 3