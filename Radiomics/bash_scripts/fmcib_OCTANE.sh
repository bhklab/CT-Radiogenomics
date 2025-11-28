#!/usr/bin/env bash
#SBATCH --job-name=fmcib_OCTANE
#SBATCH --mem=64G
#SBATCH -t 2:00:00
#SBATCH -c 16
#SBATCH --partition=himem
#SBATCH -D /cluster/home/t138199uhn/projects/readii-fmcib
#SBATCH --output="/cluster/home/t138199uhn/slurm/feature_extraction/fmcib_OCTANE_%j.out"

DATASET="PMCC_OCTANE"
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