#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 22:38:09 2024

@author: zhemingan
"""
import sys

# Number of .sh files to generate
num_files = 60

if __name__ == "__main__":
    input_fastq_file = sys.argv[1]

# Generate .sh files
for i in range(1, num_files + 1):
    # File name
    file_name = f"all_sh/split_reads_{i}.sh"
    
    # Content for each .sh file
    file_contents = f"""#!/bin/bash

#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH -t 0-3:00
#SBATCH -p short
#SBATCH -o casejob_genomic_coverage_%j.out
#SBATCH -e casejob_genomic_coverage_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load gcc/9.2.0
module load samtools/1.15.1

cd /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/
samtools coverage /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/604_B3.bam -o /n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/all_bams/genome-coverage/coverage_604_B3_Cases.txt
python /n/no_backup2/bch/lee/zheming/AD_clonal/split_reads.py "{input_fastq_file}" {i}
"""
    
    # Write content to file
    with open(file_name, "w") as f:
        f.write(file_contents)

    # Print message
    print(f"Created file: {file_name}")
