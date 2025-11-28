#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

# ---- USER INPUTS ----
cancer_types <- c("CCRCC", "HNSCC", "PDA", "BRCA", "LGG", "GBM", "KIRC") #NSCLC

id_files <- list(
  BRCA = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/BRCA_patientids.txt",
  LGG  = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/LGG_patientids.txt",
  GBM  = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/GBM_patientids.txt",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/KIRC_patientids.txt",
  #NSCLC = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/NSCLC_radiomic_patientids.txt",
  CCRCC = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/CCRCC_patientids.txt",
  HNSCC = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/HNSCC_patientids.txt",
  PDA = "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/PDA_patientids.txt"
)
data_files <- list(
  BRCA = "/Users/jackie-mac/Desktop/VSCode/data/clinical/BRCA_clinical.tsv",
  LGG  = "/Users/jackie-mac/Desktop/VSCode/data/clinical/LGG_clinical.tsv",
  GBM  = "/Users/jackie-mac/Desktop/VSCode/data/clinical/GBM_clinical.tsv",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/KIRC_clinical.tsv",
  NSCLC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/NSCLCR01Radiogenomic_DATA_LABELS_2018-05-22_1500-shifted.csv",
  CCRCC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/CCRCC_clinical.tsv",
  HNSCC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/HNSCC_clinical.tsv",
  PDA = "/Users/jackie-mac/Desktop/VSCode/data/clinical/PDA_clinical.tsv"
)
output_files <- list(
  BRCA = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/BRCA_clinical_data_filtered.csv",
  LGG  = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/LGG_clinical_data_filtered.csv",
  GBM  = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/GBM_clinical_data_filtered.csv",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/KIRC_clinical_data_filtered.csv",
  NSCLC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/NSCLC_clinical_data_filtered.csv",
  CCRCC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/CCRCC_clinical_data_filtered.csv",
  HNSCC = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/HNSCC_clinical_data_filtered.csv",
  PDA = "/Users/jackie-mac/Desktop/VSCode/data/clinical/Filtered_clinical/PDA_clinical_data_filtered.csv"
)

for (cancer in cancer_types) {
  cat("Processing", cancer, "...\n")
  ids <- scan(id_files[[cancer]], what = character(), quiet = TRUE)

  sep <- if (grepl("\\.tsv$", data_files[[cancer]], ignore.case = TRUE)) "\t" else
         if (grepl("\\.csv$", data_files[[cancer]], ignore.case = TRUE)) "," else "\t"
  data <- fread(data_files[[cancer]], data.table = FALSE, check.names = FALSE, sep = sep, quote = "\"")
  
  # Look for ID column 
  id_colname <- NULL
  if ("cases.submitter_id" %in% colnames(data)) {
    id_colname <- "cases.submitter_id"
  } else if ("Case ID" %in% colnames(data)) {
    id_colname <- "Case ID"
  } else {
    stop(sprintf("No ID column ('cases.submitter_id' or 'Case ID') found in %s", data_files[[cancer]]))
  }
  
  matched_data <- data[data[[id_colname]] %in% ids, , drop = FALSE]
  if (nrow(matched_data) == 0) {
    warning(sprintf("No matching sample IDs found in the data matrix for %s.", cancer))
  } else {
    write.csv(matched_data, output_files[[cancer]], row.names = FALSE, quote = TRUE)
    cat(sprintf("Filtered clinical data saved to %s\n", output_files[[cancer]]))
  }
}