#!/bin/bash
#SBATCH --job-name=r2r_readii_negative_CPTAC-CCRCC
#SBATCH --mem=64G
#SBATCH -t 05:59:59
#SBATCH -c 16
#SBATCH --partition=himem
#SBATCH --output="/cluster/home/t138199uhn/temp/logs/r2r_readii_negative_CPTAC-CCRCC_%j.out"

DISEASE_REGION="Abdomen"
DATA_SOURCE="TCIA"
DATASET="CPTAC-CCRCC"


USERNAME="t138199uhn"
source /cluster/home/$USERNAME/.bashrc

echo Processing $DISEASE_REGION/$DATA_SOURCE\_$DATASET

cd /cluster/home/$USERNAME/projects/readii_2_roqc

pixi run readii_negative $DATASET false true 60