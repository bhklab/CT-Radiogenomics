#!/usr/bin/env bash
#SBATCH --job-name=r2r_readii_mit_ISPY
#SBATCH --mem=512G
#SBATCH -t 05:00:00
#SBATCH -c 38
#SBATCH -D /cluster/home/t138199uhn/projects/med-imagetools
#SBATCH --partition=veryhimem
#SBATCH --output="/cluster/home/t138199uhn/slurm/feature_extraction/r2r_readii_mit_ISPY_%j.out"

DATASET="ISPY2"
USERNAME="t138199uhn"

# LOAD YOUR ENVIRONMENT HERE (uncomment/adapt one of the following)
source /cluster/home/$USERNAME/.bashrc
cd /cluster/home/$USERNAME/projects/med-imagetools

# Run the mit command (ensure imgtools is in PATH or loaded in environment)
pixi run imgtools autopipeline \
  "/cluster/projects/radiomics/PublicDatasets/srcdata/Breast/TCIA_${DATASET}/images" \
  "$PDATA/Radiogenomics/features/Breast/TCIA_${DATASET}/images/mit_${DATASET}" \
  --modalities 'MR,SEG' \
  --roi-strategy SEPARATE \
  -rmap "ROI:.*" \
  --filename-format '{PatientID}_{SampleNumber}/{Modality}_{SeriesInstanceUID}/{ImageID}.nii.gz' \
  --spacing 1.0,1.0,1.0 \
  --update-crawl \
  --existing-file-mode skip