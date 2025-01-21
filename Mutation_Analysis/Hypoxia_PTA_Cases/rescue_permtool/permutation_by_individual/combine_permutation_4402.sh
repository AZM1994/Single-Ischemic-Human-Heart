#!/bin/bash
#SBATCH --mem=5G
#SBATCH -c 4
#SBATCH -t 0-0:10
#SBATCH -p short
#SBATCH -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/rescue_log_file/permutation_by_individual/combine_permutation_4402_%j.out
#SBATCH -e /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/rescue_log_file/permutation_by_individual/combine_permutation_4402_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/6.2.0 slurm-drmaa

/home/zha255/miniconda3/envs/scan2/lib/scan2/combine_permutations.R hs37d5 /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/permutation_by_individual/4402/perms_snv_pass.rda /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/permutation_by_individual/4402/seedinfo_snv_pass.rda 4 /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/permtool_control/perms_by_sample/4402_1_A1-2n/snv_pass.rda /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/permtool_control/perms_by_sample/4402_1_A3-2n/snv_pass.rda