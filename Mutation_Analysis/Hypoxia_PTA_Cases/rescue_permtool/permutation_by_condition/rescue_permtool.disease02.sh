#!/bin/bash
#SBATCH --mem=100G
#SBATCH -c 1
#SBATCH -p medium
#SBATCH -t 1-00:00
#SBATCH -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/rescue_log_file/rescue_permtool.disease.%j.out
#SBATCH -e /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/rescue_log_file/rescue_permtool.disease.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/6.2.0 slurm-drmaa

scan2 -d rescue_disease init

scan2 -d rescue_disease config \
	--verbose \
	--analysis rescue \
	--rescue-target-fdr 0.01 \
	--scan2-object 604_B2 604_Hypoxia_PTA/call_mutations/604_B2/scan2_object.rda \
	--scan2-object 604_B3 604_Hypoxia_PTA/call_mutations/604_B3/scan2_object.rda \
	--scan2-object 604_B6 604_Hypoxia_PTA/call_mutations/604_B6/scan2_object.rda \
	--scan2-object 1113_D1 1113_Hypoxia_PTA/call_mutations/1113_D1/scan2_object.rda \
	--scan2-object 1113_E1 1113_Hypoxia_PTA/call_mutations/1113_E1/scan2_object.rda \
	--scan2-object 1113_F1 1113_Hypoxia_PTA/call_mutations/1113_F1/scan2_object.rda \
	--scan2-object 1743_A3 1743_Hypoxia_PTA/call_mutations/1743_A3/scan2_object.rda \
	--scan2-object 1743_C3 1743_Hypoxia_PTA/call_mutations/1743_C3/scan2_object.rda \
	--scan2-object 1743_F3 1743_Hypoxia_PTA/call_mutations/1743_F3/scan2_object.rda \
	--scan2-object 1363_A4 1363_Hypoxia_PTA/call_mutations/1363_A4/scan2_object.rda \
	--scan2-object 1363_D4 1363_Hypoxia_PTA/call_mutations/1363_D4/scan2_object.rda \
	--scan2-object 1363_H4 1363_Hypoxia_PTA/call_mutations/1363_H4/scan2_object.rda \
	--scan2-object 1673_CM_A2_2n 1673_Normal_PTA/call_mutations/1673_CM_A2_2n/scan2_object.rda \
	--scan2-object 1673_CM_A3_2n 1673_Normal_PTA/call_mutations/1673_CM_A3_2n/scan2_object.rda \
	--scan2-object 1673_CM_D2_2n 1673_Normal_PTA/call_mutations/1673_CM_D2_2n/scan2_object.rda \

scan2 -d rescue_disease validate

scan2 -d rescue_disease rescue --joblimit 1

cd digest_calls
digest_calls.R \
	--muts ../rescue_disease/rescued_muts.txt \
	--metadata /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/scan2.sample_to_subject_map.Hypoxia_PTA_Cases.csv \
	--individual-column donor \
	Disease_PTA_snv_indel_pass_rescue.UNFILTERED.txt

awk -F, 'NR == 1 || $16 == "FALSE"' Disease_PTA_snv_indel_pass_rescue.UNFILTERED.txt > Disease_PTA_snv_indel_pass.FILTERED.txt

cd ..

scan2 -d permtool_disease init

scan2 -d permtool_disease config \
	--verbose \
	--analysis permtool \
	--permtool-muts digest_calls/Disease_PTA_snv_indel_pass.FILTERED.txt \
	--permtool-bedtools-genome-file /home/yh174/reference/human_hg19_Broad_hs37d5/human_hg19_Broad_hs37d5.fasta.fai \
	--permtool-n-permutations 1000 \
	--permtool-sample 604_B2 604_Hypoxia_PTA \
	--permtool-sample 604_B3 604_Hypoxia_PTA \
	--permtool-sample 604_B6 604_Hypoxia_PTA \
	--permtool-sample 1113_D1 1113_Hypoxia_PTA \
	--permtool-sample 1113_E1 1113_Hypoxia_PTA \
	--permtool-sample 1113_F1 1113_Hypoxia_PTA \
	--permtool-sample 1743_A3 1743_Hypoxia_PTA \
	--permtool-sample 1743_C3 1743_Hypoxia_PTA \
	--permtool-sample 1743_F3 1743_Hypoxia_PTA \
	--permtool-sample 1363_A4 1363_Hypoxia_PTA \
	--permtool-sample 1363_D4 1363_Hypoxia_PTA \
	--permtool-sample 1363_H4 1363_Hypoxia_PTA \
	--permtool-sample 1673_CM_A2_2n 1673_Normal_PTA \
	--permtool-sample 1673_CM_A3_2n 1673_Normal_PTA \
	--permtool-sample 1673_CM_D2_2n 1673_Normal_PTA \

scan2 -d permtool_disease validate

scan2 -d permtool_disease permtool \
	--joblimit 300 \
	--cluster 'sbatch -p short  -c {threads} --mem={resources.mem_mb} -t 11:00:00 -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/rescue_log_file/slurm_log/permtool_disease_slurm-%A.log' \
	--snakemake-args ' --max-threads 12 --rerun-triggers mtime'
