#!/bin/bash
#SBATCH --job-name=snakemake_master_GSVA_run
#SBATCH --output=/cluster/home/t138199uhn/slurm/slurm_logs/GSVA_run/snakemakemaster_%j.out
#SBATCH --error=/cluster/home/t138199uhn/slurm/slurm_logs/GSVA_run/snakemakemaster_%j.err
#SBATCH --time=30:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

cd /cluster/home/t138199uhn/Scripts

module load snakemake  # or activate your environment

snakemake --snakefile Snakefile_Enrichment \
  --jobs 300 \
  --configfile enrichment_config.yaml \
  --cluster "sbatch --job-name=smk_{rule} --output=/cluster/home/t138199uhn/slurm/slurm_logs/GSVA_run/{rule}.{wildcards}.out --error=/cluster/home/t138199uhn/slurm/slurm_logs/GSVA_run/{rule}.{wildcards}.err --time={resources.runtime} --mem={resources.mem_mb} --cpus-per-task=1" \
  --latency-wait 120
