cd /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD

module load gcc/9.2.0
module load python/3.8.12
python generate_sh_MAPD_step1.py
python generate_sh_MAPD_step2.py

python /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/submit_all_sh_MAPD.py /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/step1_all_commonds

python /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/submit_all_sh_MAPD.py /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/step2_all_commonds


