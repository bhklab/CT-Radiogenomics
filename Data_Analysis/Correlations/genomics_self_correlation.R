# ===============================================================================
# Genomic Signatures Self-Correlation Analysis
# ===============================================================================
# 
# Purpose: Computes Spearman correlation matrix among genomic signatures within
#          a single pathway database (KEGG, HALLMARK, REACTOME, or BIOCARTA)
#          to identify co-regulated pathway groups.
#
# Description:
#   This script calculates pairwise Spearman correlations between all genomic
#   signatures (pathway enrichment scores) within a dataset to identify which
#   pathways are co-expressed or co-regulated. High correlations suggest 
#   functional relationships between pathways.
#
# Input Requirements:
#   1. Genomic signatures file: CSV with samples as rows, pathway signatures as columns
#   2. File must contain "SampleID" as first column header
#   3. Pathway signatures should be from a single source (e.g., all KEGG pathways)
#
# Output:
#   CSV file containing self-correlation matrix:
#   - Square matrix with pathway signatures as both rows and columns
#   - Spearman correlation coefficients ranging from -1 to +1
#   - Diagonal values = 1 (perfect self-correlation)
#   - Used for identifying highly correlated pathway clusters
#
# Analysis Method:
#   - Uses Spearman rank correlation (robust to outliers)
#   - Computes pairwise complete observations
#   - Handles missing data appropriately
#   - Results used downstream for feature filtering
#
# Usage:
#   Rscript genomics_self_correlation.R <genomic_signatures_file> <output_directory> <dataset_prefix>
#
# Example:
#   Rscript genomics_self_correlation.R HNSCC_KEGG_enrichment.csv /outputs/ HNSCC
#   Output: HNSCC_KEGG_self_correlation.csv
#
# Dependencies: Standard R libraries (stats, utils)
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================
# The input file should contain enrichment data from one pathway source already separated.
#
# Input format: CSV file with samples as rows and pathways as columns (first column = sample IDs)
#               This should be an output file from enrichment_separator.R
# Output: One self-correlation matrix CSV file
#
# Usage: Rscript genomics_self_correlation.R <genomics_file> <output_dir> <prefix>
#

suppressPackageStartupMessages({
  library(data.table)
  library(tools)
})

cat("Starting genomics self-correlation script...\n")

# =========================
# Parse Command Line Arguments
# =========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript genomics_self_correlation.R <genomics_file> <output_dir> <prefix>")
}

genomics_file <- args[1]
output_dir <- args[2]
prefix <- args[3]

cat("Input genomics file:", genomics_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Prefix:", prefix, "\n")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# =========================
# Process Input File
# =========================
cat("\n=== Processing genomics file:", genomics_file, "===\n")

# Extract source name from filename
file_basename <- basename(tools::file_path_sans_ext(genomics_file))
# Try to extract source from filename (e.g., KEGG, HALLMARK, REACTOME, BIOCARTA)
source_name <- NULL
for (src in c("kegg", "hallmark", "reactome", "biocarta")) {
  if (grepl(src, file_basename, ignore.case = TRUE)) {
    source_name <- toupper(src)  # Convert to uppercase
    break
  }
}
if (is.null(source_name)) {
  source_name <- "GENOMICS"
}
cat("Identified source:", source_name, "\n")

# Read genomics data
cat("Reading genomics data from:", genomics_file, "\n")
genomic_matrix <- fread(genomics_file, data.table = FALSE)
cat("Genomics data dimensions (including sample ID column):", dim(genomic_matrix), "\n")

# The input files from enrichment_separator.R have samples as rows and pathways as columns
# First column should be sample IDs (with "SampleID" header), set as row names and remove it
rownames(genomic_matrix) <- genomic_matrix[[1]]
genomic_matrix <- genomic_matrix[, -1, drop = FALSE]
cat("Genomics matrix dimensions (samples x pathways):", dim(genomic_matrix), "\n")

# Ensure all values are numeric
cat("Converting all values to numeric...\n")
genomic_matrix[] <- lapply(genomic_matrix, function(x) as.numeric(as.character(x)))

# Compute correlation if enough features
if (ncol(genomic_matrix) >= 2) {
  cat("Computing Spearman self-correlation matrix for", source_name, "...\n")
  cor_mat <- cor(genomic_matrix, method = "spearman", use = "pairwise.complete.obs")
  cat("Correlation matrix dimensions:", dim(cor_mat), "\n")
  
  # Save correlation matrix
  out_file <- file.path(output_dir, paste0(prefix, "_", source_name, "_self_correlation.csv"))
  write.csv(cor_mat, file = out_file, row.names = TRUE)
  cat("Self-correlation matrix saved for", source_name, "at", out_file, "\n")
  
  # Print summary statistics
  cor_values <- cor_mat[upper.tri(cor_mat)]
  cat("Correlation summary for", source_name, ":\n")
  cat("  Number of correlation pairs:", length(cor_values), "\n")
  cat("  Mean correlation:", round(mean(cor_values, na.rm = TRUE), 3), "\n")
  cat("  Median correlation:", round(median(cor_values, na.rm = TRUE), 3), "\n")
  cat("  Min correlation:", round(min(cor_values, na.rm = TRUE), 3), "\n")
  cat("  Max correlation:", round(max(cor_values, na.rm = TRUE), 3), "\n")
  
} else {
  cat("Not enough pathways for correlation analysis in", source_name, "\n")
  cat("Number of pathways:", ncol(genomic_matrix), "(minimum 2 required)\n")
}
cat("Genomics self-correlation script complete.\n")
