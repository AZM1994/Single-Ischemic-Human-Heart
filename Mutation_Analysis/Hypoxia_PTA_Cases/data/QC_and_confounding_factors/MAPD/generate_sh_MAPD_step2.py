#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 22:38:09 2024

@author: zhemingan
"""
import sys
import os

os.system("mkdir " + "step2_all_commonds")

# Number of .sh files to generate
list_of_cases = ["604_B2", "604_B3", "604_B6", "1039", "1113_D1", "1113_E1", "1113_F1", "1363_A4", "1363_D4", "1363_H4", 
                 "1673_CM_A2_2n", "1673_CM_A3_2n", "1673_CM_D2_2n", "1743_A3", "1743_C3", "1743_F3", "1864_M-E3-2n", 
                 "4402_1_A1-2n", "4402_1_A3-2n", "4638_1_A1-2n", "4638_1_A4-2n", "4638_1_B2-2n", 
                 "5657_M-D1-2n", "5657_M-H1-2n", "5828_CM_C2_2n", "5828_CM_G2_2n", 
                 "5919_1_C4-2n", "5919_1_D2-2n", "5919_1_E3-2n", "5919_1_F6-2n", "6032_1_A1-2n", "6032_M-C7-2n", "6032_M-E7-2n"]

# Generate .sh files
for i in list_of_cases:
    # File name
    file_name = f"step2_all_commonds/MAPD_step2_{i}.sh"
    # Content for each .sh file
    file_contents = f"""#!/bin/bash
#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH -t 0-00:30
#SBATCH -p short
#SBATCH -o log/MAPD_{i}_step_2_%j.out
#SBATCH -e log/MAPD_{i}_step_2_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/9.2.0
module load R/4.3.1

cd /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/{i}

Rscript /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/cbs.r /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/{i} {i} /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/{i}/{i} 50k /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/hg19.50k.k50.bad.bins.txt /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/hg19.varbin.gc.content.50k.bowtie.k50.txt

Rscript /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/copynumber.r /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/{i}/{i}.hg19.50k.k50.varbin > /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/{i}/{i}.hg19.50k.k50.varbin.copynumber.log

Rscript /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/cal_mapd.r /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/{i}/{i}.hg19.50k.k50.varbin.data.copynumber.txt > /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/{i}/{i}.hg19.50k.k50.mapd
"""
    
    # Write content to file
    with open(file_name, "w") as f:
        f.write(file_contents)

    # Print message
    print(f"Created file: {file_name}")
