# ===============================
# Clinical_filter_uniqueID.R
# -------------------------------
# Purpose:
#   Harmonizes clinical, radiomics, and multiple genomics pathway datasets by aligning sample IDs across all files.
#   Ensures all output files have consistent sample IDs for downstream multi-omics analysis.
#
# Inputs:
#   1. clinical_file: Path to the clinical data file (CSV).
#   2. radiomics_file: Path to the radiomics data file (CSV).
#   3. hallmark_genomics_file: Path to Hallmark genomics matrix (CSV).
#   4. kegg_genomics_file: Path to Kegg genomics matrix (CSV).
#   5. reactome_genomics_file: Path to Reactome genomics matrix (CSV).
#   6. biocarta_genomics_file: Path to Biocarta genomics matrix (CSV).
#   7. output_prefix: Prefix for output files.
#   8. output_dir: Directory to write harmonized outputs.
#
# Outputs:
#   - Harmonized clinical, radiomics, and genomics pathway files (all with matching sample IDs and orientation).
#   - Console messages summarizing harmonization steps.
#
# Main Steps:
#   - Reads all input files.
#   - Extracts and harmonizes sample IDs across datasets.
#   - Duplicates clinical rows to match radiomics segmentations if needed.
#   - Filters and outputs harmonized files for each data type.
# ===============================

suppressPackageStartupMessages(library(data.table))

# ---- SCRIPT DESCRIPTION ----
# This script harmonizes clinical, radiomics, and multiple genomics pathway datasets by ensuring consistent sample IDs across all files.
# It performs the following steps:
#
# 1. Reads clinical, radiomics, and four genomics pathway files (Hallmark, Kegg, Reactome, Biocarta) provided as input.
# 2. Extracts original sample IDs from radiomics and harmonizes all datasets to the same set of sample IDs.
# 3. Duplicates clinical rows to match radiomics segmentation duplicates, ensuring each segmentation ID has a corresponding clinical row.
# 4. Filters and harmonizes each genomics matrix to match the harmonized clinical sample IDs.
# 5. Writes harmonized clinical, radiomics, and all four genomics pathway datasets to separate CSV files in the specified output directory.
#
# The output files ensure that all datasets have consistent sample IDs, enabling downstream analyses that require alignment across clinical, radiomics, and genomics data.

# ---- USER INPUTS ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript Clinical_filter_uniqueID.R <clinical_file> <radiomics_file> <hallmark_genomics_file> <kegg_genomics_file> <reactome_genomics_file> <biocarta_genomics_file> <output_prefix> <output_dir>")
}
clinical_file <- args[1]
radiomics_file <- args[2]
hallmark_file <- args[3]
kegg_file <- args[4]
reactome_file <- args[5]
biocarta_file <- args[6]
output_prefix <- args[7]
output_dir <- args[8]
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- READ DATA ----
clinical <- fread(clinical_file, data.table = FALSE, quote = "\"", check.names = FALSE)
radiomics <- fread(radiomics_file, data.table = FALSE, quote = "\"", check.names = FALSE)
hallmark <- fread(hallmark_file, data.table = FALSE, quote = "\"", check.names = FALSE)
kegg <- fread(kegg_file, data.table = FALSE, quote = "\"", check.names = FALSE)
reactome <- fread(reactome_file, data.table = FALSE, quote = "\"", check.names = FALSE)
biocarta <- fread(biocarta_file, data.table = FALSE, quote = "\"", check.names = FALSE)

# ---- AUTO-DETECT RELEVANT COLUMNS ----
radiomics_patient_col <- colnames(radiomics)[1]  # First column is assumed to be sampleid_segid

# ---- EXTRACT ORIGINAL SAMPLE IDs FROM RADIOMICS ----
# Extract original sample IDs and seg IDs from sampleid_segid (format: sampleid_segid)
radiomics$original_sample_ID <- sub("_.*$", "", radiomics[[radiomics_patient_col]])
radiomics$seg_id <- sub("^[^_]*_", "", radiomics[[radiomics_patient_col]])

# ---- FILTER CLINICAL DATA ----
# Keep only clinical rows where the original sample ID matches the radiomics dataset
clinical <- clinical[clinical[[1]] %in% radiomics$original_sample_ID, , drop = FALSE]

# ---- FILTER RADIOMICS DATA ----
# Remove rows in radiomics where the original sample ID does not exist in clinical data
radiomics <- radiomics[radiomics$original_sample_ID %in% clinical[[1]], , drop = FALSE]

# ---- DUPLICATE CLINICAL DATA TO MATCH RADIOMICS DUPLICATES ----
duplicated_clinical <- data.frame()
for (i in seq_len(nrow(radiomics))) {
  sample_id <- radiomics$original_sample_ID[i]
  seg_id <- radiomics$seg_id[i]
  sampleid_segid <- radiomics[[radiomics_patient_col]][i]
  clinical_row <- clinical[clinical[[1]] == sample_id, , drop = FALSE]
  if (nrow(clinical_row) == 1) {
    clinical_row[[1]] <- sampleid_segid
    duplicated_clinical <- rbind(duplicated_clinical, clinical_row)
  }
}

# Replace clinical data with the updated duplicated clinical data
clinical <- duplicated_clinical

# ---- REMOVE UNMATCHED RADIOMICS SAMPLE ROWS ----
matched_radiomics_sample_ids <- intersect(radiomics[[radiomics_patient_col]], clinical[[1]])
unmatched_radiomics_sample_ids <- setdiff(radiomics[[radiomics_patient_col]], matched_radiomics_sample_ids)

# Remove unmatched rows
radiomics <- radiomics[radiomics[[radiomics_patient_col]] %in% matched_radiomics_sample_ids, , drop = FALSE]

# Debugging: Print dimensions of radiomics data after filtering
cat("Radiomics data dimensions after filtering:", dim(radiomics), "\n")
cat(sprintf("Removed %d unmatched radiomics sample rows.\n", length(unmatched_radiomics_sample_ids)))

# ---- FILTER EACH GENOMICS MATRIX ----
harmonize_genomics <- function(genomics, clinical_ids) {
  # Expect sample IDs in the first column (not as column names)
  sample_id_col <- colnames(genomics)[1]
  matched_rows <- genomics[[sample_id_col]] %in% clinical_ids
  genomics <- genomics[matched_rows, , drop = FALSE]
  return(genomics)
}
hallmark <- harmonize_genomics(hallmark, clinical[[1]])
kegg <- harmonize_genomics(kegg, clinical[[1]])
reactome <- harmonize_genomics(reactome, clinical[[1]])
biocarta <- harmonize_genomics(biocarta, clinical[[1]])

# ---- WRITE OUTPUT ----
# Write clinical dataset
fwrite(clinical, file.path(output_dir, paste0(output_prefix, "_harmonized_clinical.csv")), sep = ",", quote = TRUE, na = "NA")

# Write radiomics dataset
fwrite(radiomics, file.path(output_dir, paste0(output_prefix, "_harmonized_radiomics.csv")), sep = ",", quote = TRUE, na = "NA")

# Write genomic datasets in pathway x sample ID format
fwrite(hallmark, file.path(output_dir, paste0(output_prefix, "_harmonized_hallmark_genomics.csv")), sep = ",", quote = TRUE, na = "NA")
fwrite(kegg, file.path(output_dir, paste0(output_prefix, "_harmonized_kegg_genomics.csv")), sep = ",", quote = TRUE, na = "NA")
fwrite(reactome, file.path(output_dir, paste0(output_prefix, "_harmonized_reactome_genomics.csv")), sep = ",", quote = TRUE, na = "NA")
fwrite(biocarta, file.path(output_dir, paste0(output_prefix, "_harmonized_biocarta_genomics.csv")), sep = ",", quote = TRUE, na = "NA")

cat("Harmonized clinical, radiomics, and all genomics pathway files have been written.\n")