# READII → FMCIB Feature Extraction

## Overview

This tool runs FMCIB feature extraction over pre-cropped CT volumes using CSV manifests. It reads a dataset YAML to get `dataset_name`, discovers input CSVs, and writes features under a dataset-specific output folder.

To get and use the FMCIB tool, clone the following repo:

Clone the `aerts-example` branch of the `readii-fmcib` repo

```bash 
git clone https://github.com/bhklab/readii-fmcib.git -b aerts-example
```

Script: `Radiomics/readii-fmcib/run_fmcib.py` (Click CLI)

## Data layout

- Default input root (when `--csv-root` is omitted):
  - `srcdata/<dataset_name>/features/fmcib/cube_50_50_50/*.csv` (flat)
  - If no top-level CSVs exist, it will also scan one level of subfolders for legacy layouts.

- Default output root (when `--results-root` is omitted or set to `procdata`):
  - `procdata/<dataset_name>/fmcib_features/` (flat)
  - With legacy subfolders, features are written under matching subdirectories to avoid name collisions.

Symlinks are supported—`Path.glob` will follow them. A common pattern is to create repo-level links to your real storage trees:

```bash
ln -s /abs/storage/src/<dataset>  srcdata/<dataset>
ln -s /abs/storage/proc/<dataset> procdata/<dataset>
```

## CSV format

Each input CSV should minimally contain:

- `filepath` (required): absolute path to a pre-cropped CT volume (NIfTI: .nii/.nii.gz)
- Optional: `PatientID`, `SampleNumber` (ignored by FMCIB but useful for bookkeeping)

One row = one image.

Example:

```csv
filepath,PatientID,SampleNumber
/cluster/projects/.../CT_0001.nii.gz,P001,1
/cluster/projects/.../CT_0002.nii.gz,P002,1
```

## Weights

- The script looks for `model_weights.torch` in the current working directory.
- If present, it uses that file; otherwise it passes `None` to the underlying library.
- On clusters without internet egress, passing `None` often triggers a network error if the library tries to download defaults. Place the weights file locally to avoid this.

Tip: Run from `Radiomics/readii-fmcib` and put the weights file there, or run from another directory where the file exists.

## Usage

From repo root:

```bash
python Radiomics/readii-fmcib/run_fmcib.py \
  --config-path Radiomics/readii-fmcib/config/TCIA_CPTAC-CCRCC.yaml \
  --csv-root srcdata/TCIA_CPTAC-CCRCC/features/fmcib/cube_50_50_50 \
  --results-root procdata
```

From inside `Radiomics/readii-fmcib`:

```bash
python run_fmcib.py \
  --config-path config/TCIA_CPTAC-CCRCC.yaml \
  --csv-root ../../srcdata/TCIA_CPTAC-CCRCC/features/fmcib/cube_50_50_50 \
  --results-root ../../procdata
```

Notes:
- `--csv-root` is optional; if omitted, it uses `srcdata/<dataset>/features/fmcib/cube_50_50_50` from the YAML `dataset_name`.
- Features are written to `<results-root>/<dataset>/fmcib_features/…`.

## Slurm batch example

```bash
#!/usr/bin/env bash
#SBATCH --job-name=fmcib_CCRCC
#SBATCH -t 05:00:00
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -D /cluster/home/$USER/Radiogenomics/Radiomics/readii-fmcib
#SBATCH -o /cluster/home/$USER/slurm/fmcib_%j.out
#SBATCH -e /cluster/home/$USER/slurm/fmcib_%j.err

set -euo pipefail

DATASET="TCIA_CPTAC-CCRCC"
mkdir -p /cluster/home/$USER/slurm

python run_fmcib.py \
  --config-path "config/${DATASET}.yaml" \
  --csv-root "/cluster/projects/.../src/${DATASET}/features/fmcib/cube_50_50_50" \
  --results-root "/cluster/projects/.../proc"
```

## Troubleshooting

- Import error for `infer`: run the script from a directory where a `code/` folder with `infer.py` is importable (the script appends `code` to `sys.path`).
- Network is unreachable: provide a local `model_weights.torch` file; offline nodes cannot download weights.
- No CSVs found: confirm files exist under the expected `srcdata/<dataset>/features/fmcib/cube_50_50_50` or pass `--csv-root` explicitly.
