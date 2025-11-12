#!/bin/bash
#SBATCH --job-name=ISPY2_r2r_negative
#SBATCH --mem=512G
#SBATCH -t 20:00:00
#SBATCH -c 42
#SBATCH --partition=veryhimem
#SBATCH -D /cluster/home/t138199uhn/projects/readii_2_roqc
#SBATCH --output="/cluster/home/t138199uhn/temp/logs/r2r_readii_negative_ISPY2_%j.out"

DISEASE_REGION="Breast"
DATA_SOURCE="TCIA"
DATASET="ISPY2"


USERNAME="t138199uhn"
source /cluster/home/$USERNAME/.bashrc

echo Processing $DISEASE_REGION/$DATA_SOURCE\_$DATASET

cd /cluster/home/$USERNAME/projects/readii_2_roqc

pixi run readii_negative $DATASET false false 10