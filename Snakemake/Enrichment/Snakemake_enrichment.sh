#!/bin/bash
#SBATCH --job-name=snakemake_master_GSVA_run
#SBATCH --output=/cluster/home/t138199uhn/slurm/genomics_temp/GSVA_run/snakemakemaster_%j.out
#SBATCH --error=/cluster/home/t138199uhn/slurm/genomics_temp/GSVA_run/snakemakemaster_%j.err
#SBATCH -t 4:30:00
#SBATCH --mem=32G
#SBATCH --partition=himem
#SBATCH -D /cluster/home/t138199uhn/Scripts/genomics/snakemake
#SBATCH -c 16

set -euo pipefail

module load snakemake  # or activate your environment

# Ensure log directory exists for cluster job logs
mkdir -p /cluster/home/t138199uhn/slurm/genomics_temp/GSVA_run

snakemake \
  --snakefile Snakefile_Enrichment.snakefile \
  --configfile config/enrichment.yaml \
  --jobs 5 \
  --cluster "sbatch --job-name=smk_{rule} --output=/cluster/home/t138199uhn/slurm/genomics_temp/GSVA_run/{rule}.{wildcards}.out --error=/cluster/home/t138199uhn/slurm/genomics_temp/GSVA_run/{rule}.{wildcards}.err --time={resources.runtime} --mem={resources.mem_mb} --cpus-per-task=1" \
  --latency-wait 120 
