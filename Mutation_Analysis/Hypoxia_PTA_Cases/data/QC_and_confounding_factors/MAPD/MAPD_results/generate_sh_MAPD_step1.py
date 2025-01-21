#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 22:38:09 2024

@author: zhemingan
"""
import sys
import os

os.system("mkdir " + "step1_all_commonds")

# Number of .sh files to generate
list_of_cases = ["604_B2", "604_B3", "604_B6", "1039", "1113_D1", "1113_E1", "1113_F1", "1363_A4", "1363_D4", "1363_H4", 
                 "1673_CM_A2_2n", "1673_CM_A3_2n", "1673_CM_D2_2n", "1743_A3", "1743_C3", "1743_F3", "1864_M-E3-2n", 
                 "4402_1_A1-2n", "4402_1_A3-2n", "4638_1_A1-2n", "4638_1_A4-2n", "4638_1_B2-2n", 
                 "5657_M-D1-2n", "5657_M-H1-2n", "5828_CM_C2_2n", "5828_CM_G2_2n", 
                 "5919_1_C4-2n", "5919_1_D2-2n", "5919_1_E3-2n", "5919_1_F6-2n", "6032_1_A1-2n", "6032_M-C7-2n", "6032_M-E7-2n"]

# Generate .sh files
for i in list_of_cases:
    # File name
    file_name = f"step1_all_commonds/MAPD_step1_{i}.sh"
    # Content for each .sh file
    file_contents = f"""#!/bin/bash
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH -o log/MAPD_{i}_step_1_%j.out
#SBATCH -e log/MAPD_{i}_step_1_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/6.2.0
module load R/3.6.1
module load java/jdk-1.8u112
module load bwa/0.7.15
module load picard/2.8.0
module load samtools/1.3.1
module load pigz/2.3.4
module load python/2.7.12
module load openblas/0.2.19

cd /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD
mkdir -p {i} && cd {i}

samtools view -h /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/{i}.bam | python /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/varbin.50k.bam.py {i}.varbin.50k.txt {i}.varbin.50k.stats.txt /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/hg19.bin.boundaries.50k.bowtie.k50.sorted.txt /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/QC/MAPD/mapd/baslancnv/hg19.chrom.sizes.txt
"""
    
    # Write content to file
    with open(file_name, "w") as f:
        f.write(file_contents)

    # Print message
    print(f"Created file: {file_name}")
