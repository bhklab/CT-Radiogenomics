# ===============================================================================
# Radiogenomic Correlative Analysis
# ===============================================================================
# 
# Purpose: Performs Spearman correlation analysis between genomic signatures 
#          (pathway enrichment scores) and radiomic features to identify 
#          radiogenomic associations.
#
# Description:
#   This script computes pairwise Spearman correlations between filtered genomic
#   signatures and radiomic features for a specific dataset. It focuses on 
#   original and transformed radiomic features while excluding diagnostic metadata.
#   Results include correlation coefficients, p-values, and FDR corrections.
#
# Input Requirements:
#   1. Genomic signatures file: CSV with samples as rows, pathway signatures as columns
#   2. Radiomic features file: CSV with samples as rows, radiomic features as columns
#   3. Both files must have matching sample IDs as row names
#
# Output:
#   CSV file containing correlation results:
#   - Genomic signature names
#   - Radiomic feature names  
#   - Spearman correlation coefficients
#   - P-values and FDR-adjusted p-values
#   - Sample sizes used for each correlation
#
# Analysis Method:
#   - Uses Spearman rank correlation (non-parametric)
#   - Filters radiomic features by allowed prefixes (original_, wavelet-, etc.)
#   - Applies False Discovery Rate (FDR) multiple testing correction
#   - Handles missing data by pairwise complete observations
#
# Usage:
#   Rscript correlative_analysis.R <genomic_signatures_file> <radiomic_features_file> <dataset_prefix> <output_file>
#
# Example:
#   Rscript correlative_analysis.R HNSCC_KEGG_filtered.csv HNSCC_radiomics_filtered.csv HNSCC HNSCC_KEGG_correlations.csv
#
# Dependencies: Standard R libraries (stats, utils)
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

# =========================
# Radiomic Feature Categories (for reference)
# =========================
# 1. diagnostics_ (Don't use just metadata etc.)
# 2. original_ (grab)
# 3. wavelet-(grab)
# 4. square_ (grab)
# 5. squareroot_
# 6. logarithm_
# 7. exponential_
# 8. gradient_

# =========================
# Load Libraries
# =========================
suppressPackageStartupMessages({
  library(data.table)
})

# =========================
# Parse Command Line Arguments
# =========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript correlative_analysis.R <input_genes.csv> <input_radiomics.csv> <prefix> <output_file>")
}
input_genes <- args[1]
input_radiomics <- args[2]
prefix <- args[3]
output_file <- args[4]
if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)

# =========================
# Read Data
# =========================
gene_signatures <- fread(input_genes, data.table = FALSE)
rownames(gene_signatures) <- gene_signatures[[1]]
gene_signatures <- gene_signatures[, -1, drop = FALSE]

radiomic_signatures <- fread(input_radiomics, data.table = FALSE)
rownames(radiomic_signatures) <- radiomic_signatures[[1]]
radiomic_signatures <- radiomic_signatures[, -1, drop = FALSE]

# =========================
# Filter Radiomics Features by Prefix
# =========================
feature_prefixes <- c("original_", "wavelet-", "square_", "squareroot_", "logarithm_", "exponential_", "gradient_")
feature_pattern <- paste0("^(", paste(feature_prefixes, collapse = "|"), ")")
radiomic_feature_cols <- grep(feature_pattern, colnames(radiomic_signatures), value = TRUE)
if (length(radiomic_feature_cols) == 0) {
  stop("No radiomics features found with specified prefixes.")
}
radiomic_signatures <- radiomic_signatures[, radiomic_feature_cols, drop = FALSE]

# =========================
# Sample Consistency Check & Subsetting to Same Order
# =========================
samples_to_keep <- intersect(rownames(gene_signatures), rownames(radiomic_signatures))
if (length(samples_to_keep) == 0) {
  stop("No overlapping samples between genomics and radiomics datasets.")
}
gene_signatures <- gene_signatures[samples_to_keep, , drop = FALSE]
radiomic_signatures <- radiomic_signatures[samples_to_keep, , drop = FALSE]

# =========================
# Correlative Analysis: All Genomics x Filtered Radiomics Features
# =========================
results <- data.frame()
total_perms <- length(colnames(gene_signatures)) * length(colnames(radiomic_signatures))
perm_counter <- 0

for (pathway in colnames(gene_signatures)) {
  for (feature in colnames(radiomic_signatures)) {
    perm_counter <- perm_counter + 1
    if (perm_counter %% 500 == 0) {
      cat(sprintf("  Running permutation %d of %d\n", perm_counter, total_perms))
    }
    x <- gene_signatures[, pathway]
    y <- radiomic_signatures[, feature]
    # Remove NA or empty values for both vectors
    valid_idx <- which(!is.na(x) & !is.na(y) & x != "" & y != "")
    x_valid <- as.numeric(x[valid_idx])
    y_valid <- as.numeric(y[valid_idx])
    # Check for constant values or too few valid pairs
    if (length(unique(x_valid)) < 2 || length(unique(y_valid)) < 2 || length(x_valid) < 3) {
      rho <- NA
      pval <- NA
    } else {
      cor_test <- suppressWarnings(cor.test(x_valid, y_valid, method = "spearman"))
      rho <- cor_test$estimate
      pval <- cor_test$p.value
    }
    results <- rbind(
      results,
      data.frame(
        GenomicFeature = pathway,
        RadiomicFeature = feature,
        SpearmanRho = rho,
        p.value = pval
      )
    )
  }
}
# Adjust p-values for Multiple Testing
results$adj.p.value <- p.adjust(results$p.value, method = "fdr")
write.csv(
  results,
  file = output_file,
  row.names = FALSE
)
cat("Results saved at", output_file, "\n")

