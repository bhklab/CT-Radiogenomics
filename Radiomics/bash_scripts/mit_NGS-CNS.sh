#!/usr/bin/env bash
#SBATCH --job-name=r2r_readii_mit_NGS-CNS
#SBATCH --mem=96G
#SBATCH -t 05:45:00
#SBATCH -c 24
#SBATCH --partition=veryhimem
#SBATCH -D /cluster/home/t138199uhn/projects/med-imagetools
#SBATCH --output="/cluster/home/t138199uhn/slurm/feature_extraction/r2r_readii_mit_NGS-CNS_%j.out"

DATASET="NGS-CNS"
USERNAME="t138199uhn"

# LOAD YOUR ENVIRONMENT HERE 
source /cluster/home/$USERNAME/.bashrc
cd /cluster/home/$USERNAME/projects/med-imagetools

# Run the mit command (ensure imgtools is in PATH or loaded in environment)
pixi run imgtools autopipeline \
  "$RRAD/radiomics/PMCC_NGS-CNS" \
  "$PRAD/features/Brain/PMCC_NGS-CNS/images/mit_NGS-CNS" \
  --modalities all \
  --roi-strategy SEPARATE \
  -rmap "GTV:GTVp" \
  -rmap "BRAIN:Brain" \
  --filename-format '{PatientID}_{SampleNumber}/{Modality}_{SeriesInstanceUID}/{ImageID}.nii.gz' \
  --update-crawl \
  --existing-file-mode skip \
  --jobs 10