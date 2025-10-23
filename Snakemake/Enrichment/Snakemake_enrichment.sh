#!/bin/bash
#SBATCH --job-name=snakemake_master_GSVA_run
#SBATCH --output=/cluster/home/t138199uhn/slurm/genomics_temp/GSVA_run/snakemakemaster_%j.out
#SBATCH --error=/cluster/home/t138199uhn/slurm/genomics_temp/GSVA_run/snakemakemaster_%j.err
#SBATCH -t 8:00:00
#SBATCH --mem=32G
#SBATCH --partition=himem
#SBATCH -D /cluster/home/t138199uhn/Scripts/genomics/snakemake
#SBATCH -c 16

module load snakemake  # or activate your environment

snakemake --snakefile Snakefile_Enrichment.snakefile \
  --jobs 5 \
  --configfile config/enrichment.yaml \
  --cluster "sbatch --job-name=smk_{rule} --output=/cluster/home/t138199uhn/slurm/genomics_temp/GSVA_run/{rule}.{wildcards}.out --error=/cluster/home/t138199uhn/slurm/genomics_temp/GSVA_run/{rule}.{wildcards}.err --time={resources.runtime} --mem={resources.mem_mb} --cpus-per-task=1" \
  --latency-wait 120
