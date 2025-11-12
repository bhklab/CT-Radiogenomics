#!/bin/bash
#SBATCH --job-name=r2r_readii_negative_CPTAC-HNSCC
#SBATCH --mem=96G
#SBATCH -t 01:59:59
#SBATCH -c 32
#SBATCH --partition=veryhimem
#SBATCH --output="/cluster/home/t138199uhn/temp/logs/r2r_readii_negative_CPTAC-HNSCC_%j.out"

DISEASE_REGION="Abdomen"
DATA_SOURCE="TCIA"
DATASET="CPTAC-HNSCC"


USERNAME="t138199uhn"
source /cluster/home/$USERNAME/.bashrc

echo Processing $DISEASE_REGION/$DATA_SOURCE\_$DATASET

cd /cluster/home/$USERNAME/projects/readii_2_roqc

pixi run readii_negative $DATASET false false 10