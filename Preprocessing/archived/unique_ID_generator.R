# ===============================================================================
# Unique Sample ID Generator for Radiogenomic Data Integration
# ===============================================================================
# 
# Purpose: Harmonizes sample identifiers between genomic and radiomic datasets
#          to enable integrated radiogenomic analysis by creating consistent
#          sample IDs across different data modalities.
#
# Description:
#   This script takes genomic (RNA-seq) and radiomic feature datasets and
#   standardizes sample identifiers to ensure proper matching between modalities.
#   It identifies common samples, resolves ID format differences, and outputs
#   harmonized datasets ready for downstream correlation analysis.
#
# Input Requirements:
#   1. Genomics file: CSV with samples as rows, genes/features as columns
#   2. Radiomics file: CSV with samples as rows, radiomic features as columns
#   3. Sample IDs may be in different formats (TCGA submitter IDs, etc.)
#
# Output:
#   1. Harmonized genomics file: Standardized sample IDs with consistent formatting
#   2. Harmonized radiomics file: Matched sample IDs with corresponding radiomic data
#   3. Only samples present in both datasets are retained
#
# Processing Steps:
#   - Standardizes sample ID formats (removes prefixes, handles TCGA conventions)
#   - Identifies overlapping samples between genomic and radiomic datasets
#   - Filters datasets to include only common samples
#   - Ensures consistent row/column structure for downstream analysis
#
# Usage:
#   Rscript unique_ID_generator.R <genomics_file> <radiomics_file> <dataset_prefix> <output_directory>
#
# Example:
#   Rscript unique_ID_generator.R HNSCC_RNAseq.csv HNSCC_radiomics.csv HNSCC /outputs/
#
# Dependencies: data.table
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

suppressPackageStartupMessages(library(data.table))

# ---- USER INPUTS ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript unique_ID_generator.R <genomics_file> <radiomics_file> <output_prefix> <output_dir>")
}
genomics_file <- args[1]
radiomics_file <- args[2]
output_prefix <- args[3]
output_dir <- args[4]
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- READ DATA ----
genomics <- fread(genomics_file, data.table = FALSE, quote = "\"", check.names = FALSE)

# ---- Remove any fully empty rows (just in case) ----
is_empty_row <- apply(genomics, 1, function(x) all(is.na(x) | x == ""))
genomics_clean <- genomics[!is_empty_row, , drop = FALSE]

cat("Genomics matrix cleaned: empty rows removed.\n")

# ---- READ RADIOMICS DATA ----
radiomics <- fread(radiomics_file, data.table = FALSE, quote = "\"", check.names = FALSE)

# ---- AUTO-DETECT RELEVANT COLUMNS ----
# For radiomics: patient_ID is always the first column, UID is auto-detected
radiomics_patient_col <- colnames(radiomics)[1]

# Priority search for segmentation UID column: seg_series_UID first, then SeriesInstanceUID_mask
radiomics_uid_col <- grep("seg_series_UID", colnames(radiomics), ignore.case = TRUE, value = TRUE)[1]
if (is.na(radiomics_uid_col)) {
  cat("seg_series_UID column not found, searching for SeriesInstanceUID_mask...\n")
  radiomics_uid_col <- grep("SeriesInstanceUID_Mask", colnames(radiomics), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(radiomics_uid_col)) {
    stop("Neither seg_series_UID nor SeriesInstanceUID_mask column found in radiomics file.")
  } else {
    cat("Using SeriesInstanceUID_Mask as segmentation UID column.\n")
  }
} else {
  cat("Using seg_series_UID as segmentation UID column.\n")
}

# ---- STORE GENOMICS SAMPLE IDS (COLUMN NAMES, EXCLUDING FIRST COLUMN) ----
genomics_sample_ids <- colnames(genomics_clean)[-1]

# ---- FILTER RADIOMICS AND STANDARDIZE SAMPLE IDs ----
# Use partial matching: check if genomics sample ID is contained in radiomics sample ID
# AND replace radiomics sample ID with the exact genomics sample ID
radiomics_sample_ids <- radiomics[[radiomics_patient_col]]
keep_rows <- rep(FALSE, length(radiomics_sample_ids))
standardized_radiomics_ids <- character(length(radiomics_sample_ids))

for (i in seq_along(radiomics_sample_ids)) {
  radio_id <- radiomics_sample_ids[i]
  
  # Find genomics sample IDs that are contained within this radiomics sample ID
  matching_gen_ids <- genomics_sample_ids[sapply(genomics_sample_ids, function(gen_id) grepl(gen_id, radio_id, fixed = TRUE))]
  
  if (length(matching_gen_ids) == 1) {
    # Found exactly one matching genomics ID - use it to replace the radiomics ID
    keep_rows[i] <- TRUE
    standardized_radiomics_ids[i] <- matching_gen_ids[1]
    cat("Matched radiomics ID:", radio_id, "-> standardized to genomics ID:", matching_gen_ids[1], "\n")
  } else if (length(matching_gen_ids) > 1) {
    # Multiple genomics IDs match - use the longest match (most specific)
    longest_match <- matching_gen_ids[which.max(nchar(matching_gen_ids))]
    keep_rows[i] <- TRUE
    standardized_radiomics_ids[i] <- longest_match
    cat("Multiple matches for radiomics ID:", radio_id, "- using longest genomics match:", longest_match, "\n")
  } else {
    # No match found - this sample will be excluded
    keep_rows[i] <- FALSE
    cat("Warning: No genomics sample ID found matching radiomics ID:", radio_id, "- excluding\n")
  }
}

# Filter radiomics to only matched samples
radiomics <- radiomics[keep_rows, , drop = FALSE]
standardized_radiomics_ids <- standardized_radiomics_ids[keep_rows]

# Replace radiomics sample IDs with standardized genomics sample IDs
radiomics[[radiomics_patient_col]] <- standardized_radiomics_ids

cat("Number of radiomics samples after filtering:", nrow(radiomics), "\n")
cat("Number of genomics samples available:", length(genomics_sample_ids), "\n")

# Check if any samples remain after filtering
if (nrow(radiomics) == 0) {
  stop("No matching samples found between genomics and radiomics datasets. Please check sample ID formats.")
}

# ---- DUPLICATE GENOMICS COLUMNS TO MATCH RADIOMICS ----
new_genomics_mat <- genomics_clean[, 1, drop = FALSE]  # keep the first column (e.g., gene/pathway names)
matched_samples <- 0

for (i in seq_len(nrow(radiomics))) {
  # Use the standardized radiomics sample ID (now matching genomics format)
  standardized_radio_id <- radiomics[[radiomics_patient_col]][i]
  
  # Find the corresponding genomics column (should be exact match now)
  col_idx <- which(colnames(genomics_clean) == standardized_radio_id)
  
  if (length(col_idx) == 1) {
    # Found exact matching genomics column
    new_genomics_mat <- cbind(new_genomics_mat, genomics_clean[, col_idx, drop = FALSE])
    matched_samples <- matched_samples + 1
  } else {
    cat("Warning: No exact genomics column found for standardized ID:", standardized_radio_id, "\n")
  }
}

cat("Successfully matched", matched_samples, "samples between datasets.\n")

# Check if any samples were successfully matched
if (matched_samples == 0) {
  stop("No samples could be matched between genomics and radiomics datasets. Please check sample ID formats.")
}

new_genomics_mat <- new_genomics_mat[, -1, drop = FALSE]

# ---- GENERATE NEW SAMPLE IDS ----
# Use the standardized radiomics sample IDs (now matching genomics format) + segmentation UIDs
new_ids <- paste0(radiomics[[radiomics_patient_col]], "_", radiomics[[radiomics_uid_col]])

# Only assign column names if we have columns to name
if (ncol(new_genomics_mat) > 0) {
  colnames(new_genomics_mat) <- new_ids[seq_len(ncol(new_genomics_mat))]
} else {
  stop("No genomics columns remain after matching. Check sample ID compatibility.")
}

# ---- REPLACE patient_ID COLUMN IN RADIOMICS ----
radiomics[[radiomics_patient_col]] <- new_ids

# ---- ADD BACK GENE/PATHWAY NAMES COLUMN TO GENOMICS ----
final_genomics <- cbind(genomics_clean[, 1, drop = FALSE], new_genomics_mat)

# ---- WRITE OUTPUT ----
fwrite(final_genomics, file.path(output_dir, paste0(output_prefix, "_radiogenomic_RNAseq.csv")), sep = ",", quote = TRUE, na = "NA")
fwrite(radiomics, file.path(output_dir, paste0(output_prefix, "_radiogenomic_features.csv")), sep = ",", quote = TRUE, na = "NA")

cat("New genomics and radiomics files with harmonized sample IDs have been written.\n")