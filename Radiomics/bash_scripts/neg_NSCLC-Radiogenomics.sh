#!/bin/bash
#SBATCH --job-name=r2r_NSCLC-Radiogenomics_negative
#SBATCH --mem=64G
#SBATCH -t 05:59:59
#SBATCH -c 16
#SBATCH --partition=himem
#SBATCH -D /cluster/home/t138199uhn/projects/readii_2_roqc
#SBATCH --output="/cluster/home/t138199uhn/temp/logs/r2r_readii_negative_CPTAC-PDA_%j.out"

DISEASE_REGION="Abdomen"
DATA_SOURCE="TCIA"
DATASET="NSCLC-Radiogenomics"


USERNAME="t138199uhn"
source /cluster/home/$USERNAME/.bashrc
source /cluster/home/$USERNAME/.zshrc

echo Processing $DISEASE_REGION/$DATA_SOURCE\_$DATASET

cd /cluster/home/$USERNAME/projects/readii_2_roqc

pixi run readii_negative $DATASET false false 10