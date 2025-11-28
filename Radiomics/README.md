# Radiomics: Feature Extraction

This directory contains tools and scripts for extracting radiomics features used in the Radiogenomics pipeline.

## READII → FMCIB Runner

We use a Python CLI to run FMCIB feature extraction over pre-cropped CT volumes described by CSV manifests.

- Script: `Radiomics/readii-fmcib/run_fmcib.py` (Click CLI)
- Inputs (default): `srcdata/<dataset>/features/fmcib/cube_50_50_50/*.csv`
- Outputs (default): `procdata/<dataset>/fmcib_features/`
- Weights: looks for `model_weights.torch` in the current working directory

Start here for detailed usage, data layout, CSV schema, and Slurm examples:

- `Radiomics/readii-fmcib/README.md`

## Bash helpers

Batch scripts to run radiomics tooling on HPC are under `Radiomics/bash_scripts/`. Adjust paths and resource requests to your environment before submitting.
