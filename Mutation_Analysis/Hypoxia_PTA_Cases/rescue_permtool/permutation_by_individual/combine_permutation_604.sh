#!/bin/bash
#SBATCH --mem=2G
#SBATCH -c 4
#SBATCH -t 0-0:10
#SBATCH -p short
#SBATCH -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/rescue_log_file/permutation_by_individual/combine_permutation_604_%j.out
#SBATCH -e /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/rescue_log_file/permutation_by_individual/combine_permutation_604_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/6.2.0 slurm-drmaa

/home/zha255/miniconda3/envs/scan2/lib/scan2/combine_permutations.R hs37d5 ./perms_snv_pass.rda ./seedinfo_snv_pass.rda 4 ./604_B2/snv_pass.rda ./604_B3/snv_pass.rda ./604_B6/snv_pass.rda