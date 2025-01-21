#!/bin/bash
#SBATCH --mem=80G
#SBATCH -c 8
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH -o 5364_3prime_%j.out
#SBATCH -e 5364_3prime_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load cellranger/7.1.0

cd /n/scratch/users/z/zha255/heart_hypoxia/heart-scRNAseq/cellranger/

cellranger count \
--id=5364_3prime \
--transcriptome=/n/shared_db/GRCh38/uk/cellranger/7.0.0/7.0.0/refdata-gex-GRCh38-2020-A \
--fastqs=/n/scratch/users/z/zha255/heart_hypoxia/heart-scRNAseq/fastq/3-prime/5364_3prime/5364-3-IHD \
--localcores=8