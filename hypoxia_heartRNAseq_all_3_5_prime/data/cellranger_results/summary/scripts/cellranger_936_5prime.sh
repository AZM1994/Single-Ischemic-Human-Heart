#!/bin/bash
#SBATCH --mem=80G
#SBATCH -c 1
#SBATCH -t 3-00:00
#SBATCH -p medium
#SBATCH -o 936_5prime_%j.out
#SBATCH -e 936_5prime_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load cellranger/7.1.0

cd /n/no_backup2/bch/lee/zheming/single_cell_RNAseq/heart_RNAseq/cellranger/

cellranger count \
--id=936_5prime \
--transcriptome=/n/shared_db/GRCh38/uk/cellranger/7.0.0/7.0.0/refdata-gex-GRCh38-2020-A \
--fastqs=/n/no_backup2/bch/lee/zheming/single_cell_RNAseq/heart_RNAseq/fastq/936_5prime