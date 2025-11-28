#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages(library(data.table))

# Parse positional arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript Data_sampler.R <dataset> <ids.txt> <output.csv>")
}

dataset_file <- args[1]
ids_file <- args[2]
output_file <- args[3]

dataset_file <- "/Users/jackie-mac/Desktop/VSCode/data/NSCLC/NSCLC_norm_counts_FPKM.csv"
ids_file <- "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/NSCLC_patientids.txt"
output_file <- "/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/NSCLC_filtered_RNAseq.csv"


# Read patient IDs
patient_ids <- scan(ids_file, what=character(), quiet=TRUE)

# Read dataset (let fread auto-detect delimiter)
data <- fread(dataset_file, data.table=FALSE, check.names=FALSE)

# Extract columns: keep the first column (genes) and matching patient columns
gene_col <- colnames(data)[1]
cols_to_keep <- c(gene_col, intersect(patient_ids, colnames(data)))
if (length(cols_to_keep) <= 1) {
  stop("No matching patient IDs found in the dataset.")
}
filtered_data <- data[, cols_to_keep, drop=FALSE]

# Write to output CSV
write.csv(filtered_data, output_file, row.names=FALSE, quote=FALSE)
cat(sprintf("Filtered dataset saved to %s\n", output_file))