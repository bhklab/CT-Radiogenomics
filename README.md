# Radiogenomics: Integrative Analysis Pipeline

## Overview

This project looks focuses on identifying predictors of treatment response by leveraging and integrating various patient data modalities. Specifically, this project focuses on discovering associations between radiomic and genomic features predictive of patient outcomes.

We integrate two parallel initiatives: **radiomics**, which investigates radiological medical images to identify features that inform patient phenotypic attributes, and **genomics**, where we leverage patient molecular profiles to identify patterns (or "signatures") of gene expression associated with patient outcomes. At its core, this project investigates the association between radiomic features, genomic signatures, and patient outcomes.

## Research Hypothesis

**Quantitative imaging features (radiomics) correlate with underlying genomic signatures, and these correlations can predict patient treatment response and survival outcomes.**

This hypothesis is based on the understanding that:
1. Radiomic features capture tumor phenotypes resulting from underlying molecular processes
2. Gene expression patterns influence cellular behavior visible in imaging
3. Combined radiomic and genomic information provides more robust predictive models

## Objectives

### Primary Objectives
- **Discover Radiogenomic Associations**: Identify significant correlations between quantitative imaging features and genomic pathway signatures across cancer types
- **Clinical Outcome Prediction**: Develop integrated models predicting treatment response and survival outcomes
- **Biomarker Identification**: Identify radiomic-genomic feature pairs serving as predictive biomarkers

### Secondary Objectives
- **Cross-Cancer Validation**: Assess generalizability across cancer types and imaging modalities
- **Mechanistic Understanding**: Investigate biological pathways underlying radiogenomic correlations
- **Clinical Translation**: Develop clinically applicable models for routine oncological practice

## Pipeline Architecture

```
Radiogenomics/
├── Data_Analysis/           # Core analytical scripts
│   ├── Correlations/        # Radiogenomic correlation analysis
│   ├── CoxPH_models/        # Survival analysis models
│   ├── Enrichment/          # Gene set enrichment analysis
│   └── Machine_Learning_Models/  # Predictive modeling
├── Preprocessing/           # Data harmonization and quality control
├── Snakemake/              # Automated workflow management
├── Visualization_scripts/   # Results visualization
└── Old_script_version/     # Legacy code
```

## Methodology

### Data Integration
- **Radiomic Data**: Quantitative features extracted from medical image segmentation data using pyradiomics from TCIA (The Cancer Imaging Archive)
- **Genomic Data**: RNA-sequencing data from MSigDB (KEGG, Hallmark, Reactome, and BioCarta databases)
- **Clinical Data**: Treatment information, survival outcomes, and demographic data retrieved from NIH GDC (Genomic Data Commons)

### Analytical Framework
1. **Data Harmonization**: Sample identifier standardization and quality control
2. **Feature Engineering**: GSVA pathway scoring and radiomic feature preprocessing
3. **Association Discovery**: Spearman correlation analysis with FDR correction
4. **Predictive Modeling**: Cox regression and machine learning approaches
5. **Clinical Validation**: Survival analysis and outcome associations

### Key Methods
- **Pathway Analysis**: KEGG, Hallmark, Reactome, BioCarta databases
- **Correlation Analysis**: Spearman correlations between radiomic and genomic features
- **Survival Modeling**: Cox proportional hazards regression
- **Machine Learning**: LASSO, Ridge, Elastic Net regression with cross-validation

## Data Sources

- **The Cancer Genome Atlas (TCGA)**: Multi-cancer genomic data
- **Clinical Proteomic Tumor Analysis Consortium (CPTAC)**: CPTAC3 study (CCRCC, HNSCC and PDA)
- **MSigDB (Molecular Signatures Database)**: RNA-sequencing tumour gene expression data
- **TCIA (The Cancer Imaging Archive)**: Medical imaging data for radiomic feature extraction
- **NIH GDC (Genomic Data Commons)**: Clinical data including treatment information and survival outcomes for each cancer type
- **NSCLC Dataset**: Independent radiogenomics non-small cell lung cancer data source

### Data Processing Tools
- **Pyradiomics**: Feature extraction from medical image segmentation data for quantitative radiomic analysis
- **GSVA**: Gene Set Variation Analysis for pathway-level signature generation across KEGG, Hallmark, Reactome, and BioCarta databases

### Cancer Types
- **CCRCC**: Clear Cell Renal Cell Carcinoma
- **PDA**: Pancreatic Ductal Adenocarcinoma
- **HNSCC**: Head and Neck Squamous Cell Carcinoma
- **BRCA**: Breast Invasive Carcinoma
- **LGG**: Brain Lower Grade Glioma
- **GBM**: Glioblastoma Multiforme
- **KIRC**: Kidney Renal Clear Cell Carcinoma
- **NSCLC**: Non-Small Cell Lung Cancer (independent data source)

## Installation

### Prerequisites
- R ≥ 4.0 with bioinformatics packages
- Python ≥ 3.8 with Snakemake
- Sufficient computational resources for large-scale analysis

### R Dependencies
```r
# Bioinformatics packages
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "GSVA", "biomaRt", "org.Hs.eg.db"))

# Statistical and ML packages
install.packages(c("survival", "survminer", "glmnet", "caret"))

# Data manipulation and visualization
install.packages(c("data.table", "dplyr", "ggplot2", "pheatmap", "plotly", "heatmaply"))
```

### Python Dependencies
```bash
pip install snakemake
```

## Usage

### Quick Start
1. **Organize Data**: Place genomic, radiomic, and clinical data in appropriate directories
2. **Configure Paths**: Update file paths in configuration files
3. **Run Pipeline**: Execute Snakemake workflows for automated analysis
4. **Explore Results**: Use visualization scripts to interpret findings

### Workflow Execution
```bash
# Run complete pipeline
snakemake --configfile config.yaml --cores 8

# Run specific analysis
snakemake --configfile Snakemake/Correlations/correlative_config.yaml
```

### Manual Script Execution
```bash
# Preprocessing
Rscript Preprocessing/unique_ID_generator.R [args]

# Analysis
Rscript Data_Analysis/Correlations/correlative_analysis.R [args]

# Visualization
Rscript Visualization_scripts/correlative_heatmaps.R
```

## Output

### Analysis Results
- Correlation matrices between radiomic and genomic features
- Statistical significance testing with multiple comparison correction
- Survival analysis results and hazard ratios
- Machine learning model performance metrics

### Visualizations
- Interactive and static heatmaps of radiogenomic correlations
- Survival curves and forest plots
- Feature importance and correlation networks
- Clinical outcome associations

## File Formats

### Input Data
- **Genomic**: CSV files with samples as rows, genes as columns
- **Radiomic**: CSV files with samples as rows, features as columns
- **Clinical**: CSV/TSV files with patient metadata and outcomes

### Output Data
- **Correlations**: CSV matrices with correlation coefficients and p-values
- **Models**: RData objects containing fitted models
- **Visualizations**: PNG/HTML files with plots and interactive elements

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-analysis`)
3. Commit changes (`git commit -am 'Add new analysis method'`)
4. Push to branch (`git push origin feature/new-analysis`)
5. Create Pull Request


2. **Clinical Outcome Prediction**: Develop integrated models that predict treatment response and survival outcomes using combined radiomic and genomic data

3. **Biomarker Identification**: Identify specific radiomic-genomic feature pairs that serve as robust predictive biomarkers for clinical decision-making

### Secondary Objectives

1. **Cross-Cancer Validation**: Assess the generalizability of radiogenomic associations across different cancer types and imaging modalities

2. **Mechanistic Understanding**: Investigate the biological pathways that underlie radiogenomic correlations to understand the molecular basis of imaging phenotypes

3. **Clinical Translation**: Develop clinically applicable models and tools that can be integrated into routine oncological practice

## Methodology & Approach

### Multi-Modal Data Integration

Our approach integrates three critical data modalities:

#### 🖼️ **Radiomic Data**
- **Source**: Quantitative features extracted from medical images (CT, MRI, PET)
- **Content**: Shape, texture, intensity, and wavelet features characterizing tumor regions
- **Processing**: Feature standardization, correlation filtering, and quality control

#### 🧬 **Genomic Data**
- **Source**: RNA-sequencing data from tumor samples
- **Content**: Gene expression profiles transformed into pathway-level signatures
- **Processing**: Differential expression analysis, pathway enrichment scoring, and signature generation

#### 🏥 **Clinical Data**
- **Source**: Curated clincal data from trial patients
- **Content**: Treatment information, survival outcomes, demographic data
- **Processing**: Outcome harmonization, treatment filtering, and survival analysis preparation

### Analytical Framework

#### 1. **Data Harmonization**
- Sample identifier standardization across modalities
- Quality control and missing data handling
- Multi-institutional data integration protocols

#### 2. **Feature Filtering**
- Pathway-level gene expression signature generation using GSVA (Gene Set Variation Analysis)
- Radiomic feature preprocessing and correlation-based filtering
- Clinical outcome variable preparation for survival analysis

#### 3. **Association Analysis**
- Spearman correlation analysis between radiomic features and genomic signatures
- Multiple testing correction using False Discovery Rate (FDR)
- Clinical relevance filtering based on outcome associations

#### 4. **Predictive Modeling**
- Cox proportional hazards models for survival analysis
- Machine learning approaches (LASSO, Ridge, Elastic Net) for feature selection
- Cross-validation and performance assessment

#### 5. **Clinical Validation**
- Survival analysis integration throughout the pipeline
- Clinical metadata incorporation for biomarker validation
- Treatment-specific outcome associations

**Last Updated**: July 2025
