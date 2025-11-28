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

# Read patient IDs
patient_ids <- scan(ids_file, what=character(), quiet=TRUE)

# Read dataset (let fread auto-detect delimiter)
data <- fread(dataset_file, data.table=FALSE, check.names=FALSE)

# Extract columns: keep the first column (genes) and matching patient columns
gene_col <- colnames(data)[1]
kept_ids <- patient_ids[patient_ids %in% colnames(data)]  # preserve order from ids.txt
missing_ids <- setdiff(patient_ids, colnames(data))
if (length(missing_ids)) {
  message(sprintf("Skipping %d IDs not found in dataset: %s",
                  length(missing_ids), paste(head(missing_ids, 10), collapse=", ")))
  if (length(missing_ids) > 10) message("... (truncated)")
}
cols_to_keep <- c(gene_col, kept_ids)
if (length(cols_to_keep) <= 1) {
  stop("No matching patient IDs found in the dataset.")
}
filtered_data <- data[, cols_to_keep, drop=FALSE]

# Write to output CSV
write.csv(filtered_data, output_file, row.names=FALSE, quote=FALSE)
cat(sprintf("Filtered dataset saved to %s\n", output_file))

# Backup and update the ids file to include only the kept IDs
bak <- paste0(ids_file, ".bak")
invisible(try(file.copy(ids_file, bak, overwrite=TRUE), silent=TRUE))
writeLines(kept_ids, con = ids_file, sep = "\n")
cat(sprintf("Updated patient IDs written to %s (backup at %s)\n", ids_file, bak))