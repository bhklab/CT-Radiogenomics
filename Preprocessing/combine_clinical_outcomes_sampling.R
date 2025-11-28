suppressPackageStartupMessages(library(data.table))
library(tools) # for file_path_sans_ext and basename

# ---- USER INPUTS ----
# Existing cancer types
cancer_types <- c("CCRCC", "HNSCC", "PDA") 

# New cancer types with transposed survival data
new_cancer_types <- c("BRCA", "GBM", "KIRC", "LGG")

# File paths for clinical data, survival data, sample IDs, and output files
clinical_data_files <- list(
  CCRCC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/CCRCC_clinical_data_filtered.csv",
  HNSCC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/HNSCC_clinical_data_filtered.csv",
  PDA   = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/PDA_clinical_data_filtered.csv",
  BRCA  = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/BRCA_clinical_data_filtered.csv",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/KIRC_clinical_data_filtered.csv",
  LGG = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/LGG_clinical_data_filtered.csv",
  GBM = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/GBM_clinical_data_filtered.csv"
)
survival_data_files <- list(
  CCRCC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/survival/CCRCC_survival.tsv",
  HNSCC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/survival/HNSCC_survival.tsv",
  PDA   = "/Users/jackie-mac/Desktop/VSCode/data/clinical/survival/PDA_survival.tsv",
  BRCA  = "/Users/jackie-mac/Desktop/VSCode/data/clinical/survival/Human__TCGA_BRCA__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsv",
  LGG = "/Users/jackie-mac/Desktop/VSCode/data/clinical/survival/Human__TCGA_LGG__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsv",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/survival/Human__TCGA_KIRC__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsv",
  GBM = "/Users/jackie-mac/Desktop/VSCode/data/clinical/survival/Human__TCGA_GBM__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsv"
)
sample_id_files <- list(
  CCRCC = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/CCRCC_patientids.txt",
  HNSCC = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/HNSCC_patientids.txt",
  PDA   = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/PDA_patientids.txt",
  BRCA  = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/BRCA_patientids.txt",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/KIRC_patientids.txt",
  LGG = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/LGG_patientids.txt",
  GBM = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/GBM_patientids.txt"
)
output_files <- list(
  CCRCC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/combined_clinical/CCRCC_combined_clinical_data.csv",
  HNSCC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/combined_clinical/HNSCC_combined_clinical_data.csv",
  PDA   = "/Users/jackie-mac/Desktop/VSCode/data/clinical/combined_clinical/PDA_combined_clinical_data.csv",
  BRCA  = "/Users/jackie-mac/Desktop/VSCode/data/clinical/combined_clinical/BRCA_combined_clinical_data.csv",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/combined_clinical/KIRC_combined_clinical_data.csv",
  LGG = "/Users/jackie-mac/Desktop/VSCode/data/clinical/combined_clinical/LGG_combined_clinical_data.csv",
  GBM = "/Users/jackie-mac/Desktop/VSCode/data/clinical/combined_clinical/GBM_combined_clinical_data.csv"
)

columns_to_extract_clinical <- c(
  "cases.submitter_id",
  "cases.case_id",
  "treatments.treatment_outcome",
  "treatments.treatment_type",
  "treatments.treatment_dose",
  "treatments.treatment_dose_units",
  "treatments.treatment_intent_type",
  "diagnoses.ajcc_pathologic_stage",
  "demographic.gender",
  "diagnoses.age_at_diagnosis"
)
columns_to_extract_survival <- c("OS_days", "OS_event")
columns_to_extract_new_survival <- c("overall_survival", "status")

# ---- PROCESS EACH CANCER TYPE ----
for (cancer in c(cancer_types, new_cancer_types)) {
  cat("========================================\n")
  cat("Processing clinical and survival data for:", cancer, "\n")
  
  # ---- CLINICAL DATA EXTRACTION ----
  clinical_file <- clinical_data_files[[cancer]]
  sample_id_file <- sample_id_files[[cancer]]
  output_file <- output_files[[cancer]]
  
  sep <- if (grepl("\\.tsv$", clinical_file, ignore.case = TRUE)) "\t" else ","
  cat("Reading clinical data file:", clinical_file, "\n")
  clinical_data <- fread(clinical_file, sep = sep, na.strings = c("NA", ""))
  
  # Read sample IDs
  cat("Reading sample ID file:", sample_id_file, "\n")
  sample_ids <- trimws(scan(sample_id_file, what = character(), quiet = TRUE))
  cat("Total sample IDs to extract:", length(sample_ids), "\n")
  
  # Filter clinical data based on sample IDs
  found_ids_clinical <- intersect(sample_ids, clinical_data[["cases.submitter_id"]])
  not_found_ids_clinical <- setdiff(sample_ids, found_ids_clinical)
  cat("Found", length(found_ids_clinical), "matching sample IDs in clinical data.\n")
  if (length(not_found_ids_clinical) > 0) {
    cat("Sample IDs not found in clinical data:", paste(not_found_ids_clinical, collapse = ", "), "\n")
  }
  
  # Extract relevant columns from clinical data
  filtered_clinical_data <- clinical_data[cases.submitter_id %in% found_ids_clinical, ..columns_to_extract_clinical]
  cat("Extracted", nrow(filtered_clinical_data), "rows from clinical data.\n")
  
  # ---- SURVIVAL DATA EXTRACTION ----
  survival_file <- survival_data_files[[cancer]]
  sep_survival <- if (grepl("\\.tsv$", survival_file, ignore.case = TRUE)) "\t" else ","
  cat("Reading survival data file:", survival_file, "\n")
  survival_data <- fread(survival_file, sep = sep_survival, na.strings = c("NA", ""))
  
  if (cancer %in% new_cancer_types) {
    # Re-transpose survival data for new cancer types
    # Transpose and fix column names for new cancer types
    survival_data_t <- as.data.table(t(survival_data))
    # The first row now contains the clinical descriptors (column names)
    setnames(survival_data_t, as.character(unlist(survival_data_t[1, ])))
    # Add a column for sample names (which were the original column names)
    survival_data_t[, attrib_name := colnames(survival_data)]
    # Remove the first row (now used as column names)
    survival_data_t <- survival_data_t[-1]
    # Move attrib_name to the first column
    setcolorder(survival_data_t, c("attrib_name", setdiff(names(survival_data_t), "attrib_name")))
    cat("Re-transposed survival data for:", cancer, "\n")
    
    # Filter survival data based on sample IDs (using attrib_name column)
    found_ids_survival <- intersect(sample_ids, survival_data_t[["attrib_name"]])
    not_found_ids_survival <- setdiff(sample_ids, found_ids_survival)
    cat("Found", length(found_ids_survival), "matching sample IDs in survival data.\n")
    if (length(not_found_ids_survival) > 0) {
      cat("Sample IDs not found in survival data:", paste(not_found_ids_survival, collapse = ", "), "\n")
    }
    # Extract relevant columns from survival data
    filtered_survival_data <- survival_data_t[attrib_name %in% found_ids_survival, .SD, .SDcols = c("attrib_name", "overall_survival", "status")]
    cat("Extracted", nrow(filtered_survival_data), "rows from survival data.\n")
    # Rename attrib_name column to match clinical data column name (cases.submitter_id)
    setnames(filtered_survival_data, "attrib_name", "cases.submitter_id")
    # Rename survival columns to standard names for output
    if ("overall_survival" %in% colnames(filtered_survival_data)) {
      setnames(filtered_survival_data, "overall_survival", "OS_days")
    }
    if ("status" %in% colnames(filtered_survival_data)) {
      setnames(filtered_survival_data, "status", "OS_event")
    }
  } else {
    # Filter survival data for existing cancer types
    found_ids_survival <- intersect(sample_ids, survival_data[["case_id"]])
    not_found_ids_survival <- setdiff(sample_ids, found_ids_survival)
    cat("Found", length(found_ids_survival), "matching sample IDs in survival data.\n")
    if (length(not_found_ids_survival) > 0) {
      cat("Sample IDs not found in survival data:", paste(not_found_ids_survival, collapse = ", "), "\n")
    }
    # Extract relevant columns from survival data
    filtered_survival_data <- survival_data[case_id %in% found_ids_survival, c("case_id", ..columns_to_extract_survival), with = FALSE]
    cat("Extracted", nrow(filtered_survival_data), "rows from survival data.\n")
    # Rename case_id column to match clinical data column name (cases.submitter_id)
    setnames(filtered_survival_data, "case_id", "cases.submitter_id")
  }
  
  # Merge clinical and survival data
  combined_data <- merge(filtered_clinical_data, filtered_survival_data, by = "cases.submitter_id", all.x = TRUE)
  cat("Merged clinical and survival data for", cancer, ".\n")
  
  # Remove rows where survival column is NA or blank
  survival_column <- if (cancer %in% new_cancer_types) "OS_days" else "OS_days"
  combined_data <- combined_data[!is.na(get(survival_column)) & get(survival_column) != "", ]
  cat("Removed rows with NA or blank survival data. Remaining rows:", nrow(combined_data), "\n")
  
  # ---- WRITE OUTPUT ----
  fwrite(combined_data, output_file, na = "NA", quote = TRUE)
  cat("Combined clinical data saved to:", output_file, "\n")
  cat("========================================\n\n")
}