# Radiogenomics: Integrative Analysis Pipeline

## Overview

This repository contains an end-to-end pipeline to discover, visualize, and validate associations between quantitative imaging features (radiomics) and molecular profiles (genomics) across cancer types. It provides:

- Radiomics feature processing and feature-level clustering to select compact “medoid” feature sets
- Dimensionality reduction utilities for fast 2D/3D visualization (PCA, t-SNE, UMAP)
- Radiogenomic correlation analysis (Spearman and Pearson) with ID harmonization
- Survival modeling and enrichment analyses (legacy/optional; see project pivot below)

## Research hypothesis

Quantitative imaging features capture phenotypes driven by underlying molecular programs. Therefore, radiomic features should correlate with genomic pathway signatures. Note: clinical prediction (treatment response/survival) was investigated initially but is currently de‑emphasized after weak associations were observed (see pivot).

## Project history and pivot

- Original aim: evaluate whether radiomic features can serve as surrogate predictive markers for cancer diagnosis and survival outcomes, and whether radiomics/genomics complement each other for clinical prediction.
- Findings: across explored cohorts, associations between (i) radiomics and genomics, and (ii) either modality and clinical outcomes, were weak.
- Pivot: narrow the scope to characterizing radiomic–genomic relationships only, with emphasis on pathway-level signatures (MSigDB Hallmarks), and drop clinical analyses from the main line of inquiry.
- Feature extraction change: moved from hand-crafted PyRadiomics to deep learning features using FMCIB (Foundation Model for Cancer Imaging Biomarkers). Current analyses center on these deep features.

### Current objectives

- Identify robust, reproducible associations between deep learning FMCIB radiomic features and Hallmarks-of-cancer gene signatures.
- Reduce radiomic feature redundancy via spectral clustering and medoid selection to improve interpretability.
- Characterize cross-cohort consistency (e.g., NSCLC and CPTAC cohorts) of radiogenomic signals.

## Repository organization

```
Radiogenomics/
├── Data_Analysis/                           # Analysis notebooks, scripts and helper utilities
│   ├── Correlations/                        # Radiogenomic correlation & selection pipelines
│   │   ├── LASSO_corr.py                    # Per-signature LASSO -> test-set Spearman correlations + heatmaps
│   │   ├── PCA_correlation.py               # PCA -> correlation pipeline (PCA selection + correlations/heatmaps)
│   │   ├── lasso_per_signature.py           # helpers to run per-signature LASSO jobs
│   │   ├── radiogenomic_correlation.py      # Multi-task LASSO & correlation pipeline (residualization, variance filtering)
│   │   └── radiomics_self_correlation.R     # small R utilities for radiomics self-correlation
│   ├── Machine_Learning_Models/             # Linear/LASSO/Ridge/Elastic Net models (R implementations)
│   ├── Enrichment/                          # Pathway/gene set enrichment (R scripts)
│   ├── radiomics/                           # radiomics analysis helpers and feature-selection kernels
│   │   └── feature-selection/
│   │       ├── spectral_cluster.py
│   │       ├── hierarchical_medoid_clustering.py
│   │       ├── feature_PCA.py
│   │       ├── feature_vif_pca.py
│   │       ├── feature_volume_corr.py
│   │       └── feature_volume_sPCA.py
│   ├── misc_tasks.ipynb                     # utility notebook: heatmaps, PCA plotting, LASSO presence scans
│   └── README.md
├── Preprocessing/                           # Data harmonization & preprocessing scripts
│   ├── extract_tpm_by_symbols.py
│   ├── Data_sampler.R
│   └── archived/                            # older preprocessing scripts
├── Radiomics/                               # FMCIB feature extraction + run scripts and configs
│   ├── readii-fmcib/
	│   ├── run_fmcib.py
	│   └── config/                          # per-dataset FMCIB configs (TCIA, CPTAC, etc.)
│   ├── readii_2_roqc/
│   └── bash_scripts/                        # dataset-specific run scripts
├── Snakemake/                               # Workflow definitions and Snakefiles (Correlations, Enrichment, Clinical)
│   ├── Correlations/
│   ├── Enrichment/
│   └── clinical_associations/
├── Old_script_version/                      # legacy scripts kept for reference (includes older visualization kernels)
│   └── Visualization_scripts/
│       ├── UMAP_cancer.R
│       ├── clinical_correlation_viz.R
│       ├── correlative_heatmaps.R
│       └── signature_distribution.R
└── data/                                    # Raw and processed data (see data/procdata and data/rawdata)
```

If a folder has its own README, start there for details (e.g., Data_Analysis/Enrichment, Machine_Learning_Models).

## Data modalities and conventions

- Radiomics (current): deep learning features from FMCIB (columns prefixed with `pred_`), with `SampleID` formatted like `R01-###_####` where the suffix identifies an instance/region.
- Radiomics (historical): traditional PyRadiomics features were evaluated earlier but are not the current focus.
- Genomics: gene signatures (e.g., MSigDB Hallmarks) with samples as rows indexed by base IDs (e.g., `R01-###`).
- Clinical: outcomes and covariates for survival analyses.

Best-practice conventions:
- Keep radiomics rows filtered to permutation == original and region == full (where applicable).
- Use standardized feature matrices for downstream correlation to reduce redundancy.
- Maintain consistent sample ID schemes across modalities; derive base IDs by stripping the trailing `_####` from `SampleID` when aligning to genomics.

## Installation

You can use any Python 3.9+ environment and R ≥ 4.0. A minimal Python stack is listed below.

### R packages (typical)
- Bioconductor core: DESeq2, GSVA, biomaRt, org.Hs.eg.db
- Stats/ML: survival, survminer, glmnet, caret
- Viz: ggplot2, pheatmap, plotly, heatmaply

### Python packages (typical)
- numpy, pandas, scipy
- scikit-learn
- matplotlib, seaborn
- umap-learn (only for UMAP plots)
- snakemake (if running workflows)


## Datasets and signatures

- Datasets: CPTAC-CCRCC, CPTAC-PDA, CPTAC-HNSCC, TCGA-KIRC, NSCLC-Radiogenomics, OCTANE
- Gene signatures: MSigDB Hallmarks (CSV matrices with samples as rows and signatures as columns). Example files live under `data/procdata/gene_signatures/`.

Data syncing and curation are environment-specific. Keep signature matrices and radiomics matrices versioned and documented.

## Best practices

- Reproducibility: set seeds (e.g., `--seed 10`) and record package versions.
- ID hygiene: consistently maintain `R01-###_####` for radiomics and base IDs `R01-###` for genomics.
- Scaling: visualization utilities standardize features; avoid re-scaling inputs twice.
- Missing data: numeric coercion + median imputation is applied in visualization and clustering steps; for correlations, pairwise NaNs are skipped.
- Plots at scale: t-SNE perplexity must be < number of observations; the script will auto-adjust when needed.
- Storage: keep derived matrices (standardized, medoid) alongside cluster assignments to ease reproducibility.

## Legacy components

These parts of the repo reflect the original broader scope and are kept for reference:

- `Data_Analysis/CoxPH_models`, `Snakemake/clinical_associations`: clinical outcome modeling and workflows.
- `Old_script_version/` and `Preprocessing/archived/`: earlier enrichment and clinical preprocessing utilities.
- PyRadiomics-derived feature scripts/outputs (historical); current analyses prioritize FMCIB deep features.

## File formats

Inputs
- Genomics: CSV, samples as rows (index), signatures as columns
- Radiomics: CSV, rows=samples with `SampleID` and `pred_*` columns
- Clinical: CSV/TSV with outcomes and covariates

Outputs
- Correlations: CSV matrices and PNG heatmaps
- Clustering: cluster assignments, Gap-statistic tables, spectral/t-SNE embeddings, medoid matrices
- Models: RData objects and figures

## Contributing

1) Fork → 2) Branch → 3) Commit → 4) Push → 5) PR
—

Last updated: November 2025
