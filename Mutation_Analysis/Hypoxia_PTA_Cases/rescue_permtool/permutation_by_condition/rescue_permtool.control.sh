#!/bin/bash
#SBATCH --mem=80G
#SBATCH -c 1
#SBATCH -t 1-12:00
#SBATCH -p medium
#SBATCH -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/rescue_log_file/rescue_permtool.control_%j.out
#SBATCH -e /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/rescue_log_file/rescue_permtool.control_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/6.2.0 slurm-drmaa

scan2 -d rescue_control init

scan2 -d rescue_control config \
	--verbose \
	--analysis rescue \
	--rescue-target-fdr 0.01 \
	--scan2-object 1039 /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/1039_Hypoxia_PTA/call_mutations/1039/scan2_object.rda \
	--scan2-object 1864_M-E3-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/1864_Normal_PTA/call_mutations/1864_M-E3-2n/scan2_object.rda \
	--scan2-object 4638_1_A1-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/4638_Normal_PTA/call_mutations/4638_1_A1-2n/scan2_object.rda \
	--scan2-object 4638_1_A4-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/4638_Normal_PTA/call_mutations/4638_1_A4-2n/scan2_object.rda \
	--scan2-object 4638_1_B2-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/4638_Normal_PTA/call_mutations/4638_1_B2-2n/scan2_object.rda \
	--scan2-object 5657_M-D1-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/5657_Normal_PTA/call_mutations/5657_M-D1-2n/scan2_object.rda \
	--scan2-object 5657_M-H1-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/5657_Normal_PTA/call_mutations/5657_M-H1-2n/scan2_object.rda \
	--scan2-object 5919_1_C4-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/5919_Normal_PTA/call_mutations/5919_1_C4-2n/scan2_object.rda \
	--scan2-object 5919_1_D2-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/5919_Normal_PTA/call_mutations/5919_1_D2-2n/scan2_object.rda \
	--scan2-object 5919_1_E3-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/5919_Normal_PTA/call_mutations/5919_1_E3-2n/scan2_object.rda \
	--scan2-object 5919_1_F6-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/5919_Normal_PTA/call_mutations/5919_1_F6-2n/scan2_object.rda \
	--scan2-object 6032_1_A1-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/6032_Normal_PTA/call_mutations/6032_1_A1-2n/scan2_object.rda \
	--scan2-object 6032_M-C7-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/6032_Normal_PTA/call_mutations/6032_M-C7-2n/scan2_object.rda \
	--scan2-object 6032_M-E7-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/6032_Normal_PTA/call_mutations/6032_M-E7-2n/scan2_object.rda \
	--scan2-object 4402_1_A1-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/4402_Normal_PTA/call_mutations/4402_1_A1-2n/scan2_object.rda \
	--scan2-object 4402_1_A3-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/4402_Normal_PTA/call_mutations/4402_1_A3-2n/scan2_object.rda \
	--scan2-object 5828_CM_C2_2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/5828_Normal_PTA/call_mutations/5828_CM_C2_2n/scan2_object.rda \
	--scan2-object 5828_CM_G2_2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/5828_Normal_PTA/call_mutations/5828_CM_G2_2n/scan2_object.rda \

scan2 -d rescue_control validate

scan2 -d rescue_control rescue --joblimit 1

cd digest_calls
digest_calls.R \
	--muts ../rescue_control/rescued_muts.txt \
	--metadata /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/scan2.sample_to_subject_map.Hypoxia_PTA_Cases.csv \
	--individual-column donor \
	Control_PTA_snv_indel_pass_rescue.UNFILTERED.txt

awk -F, 'NR == 1 || $16 == "FALSE"' Control_PTA_snv_indel_pass_rescue.UNFILTERED.txt > Control_PTA_snv_indel_pass.FILTERED.txt

cd ..

scan2 -d permtool_control init

scan2 -d permtool_control config \
	--verbose \
	--analysis permtool \
	--permtool-muts digest_calls/Control_PTA_snv_indel_pass.FILTERED.txt \
	--permtool-bedtools-genome-file /home/yh174/reference/human_hg19_Broad_hs37d5/human_hg19_Broad_hs37d5.fasta.fai \
	--permtool-n-permutations 1000 \
	--permtool-sample 1039 1039_Hypoxia_PTA \
	--permtool-sample 1864_M-E3-2n 1864_Normal_PTA \
	--permtool-sample 4638_1_A1-2n 4638_Normal_PTA \
	--permtool-sample 4638_1_A4-2n 4638_Normal_PTA \
	--permtool-sample 4638_1_B2-2n 4638_Normal_PTA \
	--permtool-sample 5657_M-D1-2n 5657_Normal_PTA \
	--permtool-sample 5657_M-H1-2n 5657_Normal_PTA \
	--permtool-sample 5919_1_C4-2n 5919_Normal_PTA \
	--permtool-sample 5919_1_D2-2n 5919_Normal_PTA \
	--permtool-sample 5919_1_E3-2n 5919_Normal_PTA \
	--permtool-sample 5919_1_F6-2n 5919_Normal_PTA \
	--permtool-sample 6032_1_A1-2n 6032_Normal_PTA \
	--permtool-sample 6032_M-C7-2n 6032_Normal_PTA \
	--permtool-sample 6032_M-E7-2n 6032_Normal_PTA \
	--permtool-sample 4402_1_A1-2n 4402_Normal_PTA \
	--permtool-sample 4402_1_A3-2n 4402_Normal_PTA \
	--permtool-sample 5828_CM_C2_2n 5828_Normal_PTA \
	--permtool-sample 5828_CM_G2_2n 5828_Normal_PTA \

scan2 -d permtool_control validate

scan2 -d permtool_control permtool \
	--joblimit 300 \
	--cluster 'sbatch -p short  -c {threads} --mem={resources.mem_mb} -t 11:00:00 -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/rescue_log_file/permtool_control_slurm-%A.log' \
	--snakemake-args ' --max-threads 12 --rerun-triggers mtime'
