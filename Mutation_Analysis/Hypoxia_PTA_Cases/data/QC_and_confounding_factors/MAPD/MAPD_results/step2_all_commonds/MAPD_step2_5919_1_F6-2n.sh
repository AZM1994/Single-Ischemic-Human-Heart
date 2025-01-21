#!/bin/bash
#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH -t 0-00:30
#SBATCH -p short
#SBATCH -o log/MAPD_5919_1_F6-2n_step_2_%j.out
#SBATCH -e log/MAPD_5919_1_F6-2n_step_2_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/9.2.0
module load R/4.3.1

cd /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/5919_1_F6-2n

Rscript /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/cbs.r /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/5919_1_F6-2n 5919_1_F6-2n /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/5919_1_F6-2n/5919_1_F6-2n 50k /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/hg19.50k.k50.bad.bins.txt /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/hg19.varbin.gc.content.50k.bowtie.k50.txt

Rscript /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/copynumber.r /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/5919_1_F6-2n/5919_1_F6-2n.hg19.50k.k50.varbin > /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/5919_1_F6-2n/5919_1_F6-2n.hg19.50k.k50.varbin.copynumber.log

Rscript /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/cal_mapd.r /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/5919_1_F6-2n/5919_1_F6-2n.hg19.50k.k50.varbin.data.copynumber.txt > /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/5919_1_F6-2n/5919_1_F6-2n.hg19.50k.k50.mapd
