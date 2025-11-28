#!/usr/bin/env bash
#SBATCH --job-name=r2r_readii_mit_CPTAC-CCRCC
#SBATCH --mem=64G
#SBATCH -t 5:00:00
#SBATCH -c 16
#SBATCH --partition=himem
#SBATCH -D /cluster/home/t138199uhn/projects/med-imagetools
#SBATCH --output="/cluster/home/t138199uhn/slurm/feature_extraction/r2r_readii_mit_CPTAC-CCRCC_%j.out"

DATASET="CPTAC-CCRCC"
USERNAME="t138199uhn"

# LOAD YOUR ENVIRONMENT HERE 
source /cluster/home/$USERNAME/.bashrc
cd /cluster/home/$USERNAME/projects/med-imagetools

# Run the mit command (ensure imgtools is in PATH or loaded in environment)
pixi run imgtools autopipeline \
  "/cluster/projects/radiomics/PublicDatasets/srcdata/Abdomen/TCIA_${DATASET}/images" \
  "$PDATA/Radiogenomics/features/Abdomen/TCIA_${DATASET}/images/mit_${DATASET}" \
  --modalities 'CT,RTSTRUCT' \
  --roi-strategy SEPARATE \
  -rmap "GTV:.*" \
  --filename-format '{PatientID}_{SampleNumber}/{Modality}_{SeriesInstanceUID}/{ImageID}.nii.gz' \
  --spacing 1.0,1.0,1.0 \
  --update-crawl \
  --existing-file-mode skip