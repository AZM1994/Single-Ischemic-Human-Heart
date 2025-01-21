#!/bin/bash
#SBATCH --mem=80G
#SBATCH -c 8
#SBATCH -t 0-16:00
#SBATCH -p medium
#SBATCH -o 5874_3prime_%j.out
#SBATCH -e 5874_3prime_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zheming.an@childrens.harvard.edu

module load cellranger/7.1.0

cd /n/no_backup2/bch/lee/zheming/test/

cellranger count \
--id=5874_3prime \
--transcriptome=/n/shared_db/GRCh38/uk/cellranger/7.0.0/7.0.0/refdata-gex-GRCh38-2020-A \
--fastqs=/n/no_backup2/bch/lee/zheming/test/5874-3-IHD \
--localcores=8