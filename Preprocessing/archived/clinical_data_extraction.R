# ===============================
# clinical_data_extraction.R
# -------------------------------
# Purpose:
#   Extracts and filters clinical data for a specific cancer type, matching a set of sample IDs and filtering for Chemotherapy or Radiation Therapy, NOS as treatment types.
#   Designed for use in radiogenomics pipelines where harmonized clinical data is required for downstream analysis.
#
# Inputs:
#   1. input_clinical_file: Path to the clinical data file (CSV or TSV).
#   2. sample_id_file: Path to a file containing sample IDs to extract (one per line).
#   3. output_file: Path to write the filtered clinical data.
#   4. prefix: Cancer type or cohort prefix (used for column renaming).
#
# Outputs:
#   - Filtered clinical data file with columns renamed to include the prefix.
#   - Console messages summarizing extraction and filtering steps.
#
# Main Steps:
#   - Reads clinical data and sample IDs.
#   - Filters for Chemotherapy or Radiation Therapy, NOS as treatment types.
#   - Extracts rows matching sample IDs, including duplicates.
#   - Renames columns with the cohort prefix.
#   - Writes the filtered data to output.
# ===============================

suppressPackageStartupMessages(library(data.table))
library(tools) # for file_path_sans_ext and basename

# ---- USER INPUTS ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript clinical_data_extraction.R <input_clinical_file> <sample_id_file> <output_file> <prefix>")
}
input_clinical_file <- args[1]
sample_id_file <- args[2]
output_file <- args[3]
prefix <- args[4]

columns_to_extract <- c(
  "cases.submitter_id",
  "cases.case_id",
  "treatments.treatment_outcome",
  "treatments.treatment_type",
  "treatments.treatment_dose",
  "treatments.treatment_dose_units",
  "treatments.treatment_intent_type",
  "diagnoses.ajcc_pathologic_stage",
  "demographic.gender",
  "diagnoses.age_at_diagnosis",
  "OS_event",
  "OS_days"  # Overall survival days
)

cat("Extracting clinical data for:", prefix, "\n")
sep <- if (grepl("\\.tsv$", input_clinical_file, ignore.case = TRUE)) "\t" else ","
cat("Reading dataset file:", input_clinical_file, "\n")
dataset <- fread(input_clinical_file, sep = sep, na.strings = c("NA", ""))

# Read sample IDs (strip whitespace)
sample_ids <- trimws(scan(sample_id_file, what = character(), quiet = TRUE))
cat("Total sample IDs to extract:", length(sample_ids), "\n")

# Hard-coded filter for Chemotherapy or Radiation Therapy, NOS
treatment_types_to_keep <- c("Chemotherapy", "Radiation Therapy, NOS")
cat("Filtering for treatment types:", paste(treatment_types_to_keep, collapse = ", "), "\n")
dataset <- dataset[treatments.treatment_type %in% treatment_types_to_keep]
cat("Rows remaining after treatment type filter:", nrow(dataset), "\n")

# Find matching and non-matching sample IDs
found_ids <- intersect(sample_ids, dataset[["cases.submitter_id"]])
not_found_ids <- setdiff(sample_ids, found_ids)
cat("Found", length(found_ids), "matching sample IDs in dataset.\n")
if (length(not_found_ids) > 0) {
  cat("Sample IDs not found in", prefix, ":", paste(not_found_ids, collapse = ", "), "\n")
}

# Extract all rows for found sample IDs (including duplicates)
subset_data <- dataset[cases.submitter_id %in% found_ids, ..columns_to_extract]
cat("Extracted", nrow(subset_data), "rows for", prefix, "\n")

# Rename columns with suffix (e.g., cases.submitter_id_BRCA)
setnames(subset_data, columns_to_extract, paste0(columns_to_extract, "_", prefix))

# Write the output file for this cancer type
fwrite(subset_data, output_file, na = "NA", quote = TRUE)
cat("Extraction complete for", prefix, ". Output written to:", output_file, "\n")