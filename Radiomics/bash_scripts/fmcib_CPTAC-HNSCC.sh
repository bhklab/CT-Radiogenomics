#!/usr/bin/env bash
#SBATCH --job-name=fmcib_HNSCC
#SBATCH --mem=64G
#SBATCH -t 5:00:00
#SBATCH -c 16
#SBATCH --partition=himem
#SBATCH -D /cluster/home/t138199uhn/projects/readii-fmcib
#SBATCH --output="/cluster/home/t138199uhn/slurm/feature_extraction/fmcib_CPTAC-HNSCC_%j.out"

DATASET="TCIA_CPTAC-HNSCC"
USERNAME="t138199uhn"

# LOAD YOUR ENVIRONMENT HERE 
source /cluster/home/$USERNAME/.bashrc
source /cluster/home/$USERNAME/.zshrc

# Run the fmcib command (ensure imgtools is in PATH or loaded in environment)
# inside your .sbatch script
pixi run -e fmcib -- python run_fmcib.py \
  --config-path "config/${DATASET}.yaml" \
  --csv-root "data/srcdata/${DATASET}/features/fmcib/cube_50_50_50" \
  --results-root "data/procdata"

pixi run python run_fmcib.py \
  --config-path "config/TCIA_CPTAC-HNSCC.yaml" \
  --csv-root "data/srcdata/TCIA_CPTAC-HNSCC/features/fmcib/cube_50_50_50" \
  --results-root "data/procdata"