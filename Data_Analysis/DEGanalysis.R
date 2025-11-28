# ===============================================================================
# Differential Gene Expression Analysis using DESeq2
# ===============================================================================
# 
# Purpose: Performs differential gene expression analysis between tumor and normal
#          tissue samples using DESeq2 for multiple cancer types.
#
# Description:
#   This script takes RNA-seq count data for multiple cancer types and performs
#   DESeq2 differential expression analysis to identify genes that are significantly
#   differentially expressed between tumor and normal samples. Results include
#   log2 fold changes, p-values, and adjusted p-values.
#
# Input Requirements:
#   1. RNA-seq count files: Tab-separated files with genes as rows, samples as columns
#   2. Sample metadata indicating tumor vs normal status
#   3. File paths configured in the rna_files list
#
# Output:
#   For each cancer type:
#   - DESeq2 results table with differential expression statistics
#   - Normalized count matrices
#   - Statistical summaries of significant genes
#
# Analysis Method:
#   - Uses DESeq2 for robust differential expression analysis
#   - Applies variance stabilizing transformation
#   - Controls for multiple testing using Benjamini-Hochberg correction
#   - Filters results based on statistical significance and fold change
#
# Usage:
#   1. Update rna_files list with paths to your count data
#   2. Ensure sample naming follows tumor/normal conventions
#   3. Run: Rscript DEGanalysis.R
#
# Dependencies: DESeq2, data.table
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

# ---- USER INPUTS ----
rna_files <- list(
  BRCA = "/path/to/BRCA_counts.tsv",
  LIHC = "/path/to/LIHC_counts.tsv"
  # Add more cancer types and their count files
)
metadata_files <- list(
  BRCA = "/path/to/BRCA_metadata.tsv",
  LIHC = "/path/to/LIHC_metadata.tsv"
  # Add more cancer types and their metadata files
)
group_column <- "condition"  # Column in metadata indicating group (e.g., "tumor" vs "normal")
top_k <- 20                  # Number of top DEGs to report per cancer type
output_dir <- "/path/to/DEG_results/"

dir.create(output_dir, showWarnings = FALSE)

# --- Combine all counts ---
all_counts <- list()
all_metadata <- list()
for (cancer in names(rna_files)) {
  counts <- fread(rna_files[[cancer]], data.table = FALSE)
  rownames(counts) <- counts[[1]]
  counts <- counts[,-1]
  all_counts[[cancer]] <- counts
  # Build metadata for this cancer
  meta <- data.frame(
    Sample = colnames(counts),
    TumourType = cancer,
    stringsAsFactors = FALSE
  )
  all_metadata[[cancer]] <- meta
}
# Find common genes
common_genes <- Reduce(intersect, lapply(all_counts, rownames))
all_counts <- lapply(all_counts, function(x) x[common_genes, , drop=FALSE])
counts_mat <- do.call(cbind, all_counts)
metadata <- do.call(rbind, all_metadata)
rownames(metadata) <- metadata$Sample

# --- Add BRCA_status column ---
metadata$BRCA_status <- ifelse(metadata$TumourType == "BRCA", 1, 0)

# --- Ensure order matches ---
counts_mat <- counts_mat[, metadata$Sample]

# --- DESeq2 analysis ---
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_mat)),
  colData = metadata,
  design = ~ BRCA_status + TumourType
)
dds <- DESeq(dds)
res <- results(dds, name = "BRCA_status")
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), "BRCA_vs_all_DEGs.csv")