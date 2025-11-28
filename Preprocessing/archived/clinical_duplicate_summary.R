library(data.table)

# ---- SCRIPT DESCRIPTION ----
# This script processes clinical intersect files for multiple cancer types to summarize treatment types 
# associated with each sample ID. It performs the following steps:
#
# 1. Reads clinical intersect files for specified cancer types.
# 2. Identifies relevant columns for sample IDs and treatment types.
# 3. Filters out rows with invalid or missing treatment types (e.g., NA, blank, or '--').
# 4. Collects all unique sample IDs and treatment types across datasets.
# 5. Maps each sample ID to its associated cancer type and treatment types.
# 6. Constructs a summary matrix indicating whether each sample received each treatment type.
# 7. Converts the summary matrix into a data frame and adds a column for cancer type.
# 8. Writes the summary data frame to a CSV file, ensuring treatment types with commas are quoted.
#
# The output file provides a comprehensive summary of treatment types for each sample ID across cancer
# types, facilitating downstream analysis and visualization.

# ---- USER INPUTS ----
input_files <- list(
  BRCA = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/BRCA_clinical_intersect.csv",
  GBM  = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/GBM_clinical_intersect.csv",
  LGG  = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/LGG_clinical_intersect.csv",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/KIRC_clinical_intersect.csv",
  CCRCC = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/CCRCC_clinical_intersect.csv",
  HNSCC = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/HNSCC_clinical_intersect.csv",
  PDA   = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/PDA_clinical_intersect.csv"
)
output_file <- "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/clinical_treatment_summary.csv"

# Store all unique sample IDs and treatment types
all_samples <- character()
all_treatments <- character()
sample_treatments <- list()

for (dataset in names(input_files)) {
  file <- input_files[[dataset]]
  cat("Processing:", dataset, "\n")
  df <- fread(file, data.table = FALSE, check.names = FALSE, stringsAsFactors = FALSE, quote = "\"")
  
  # Find relevant columns
  id_col <- grep("^sample_id", colnames(df), value = TRUE)[1]
  treat_col <- grep("^treatments.treatment_type", colnames(df), value = TRUE)[1]
  if (is.na(id_col) || is.na(treat_col)) {
    cat("  Columns not found for", dataset, "- skipping.\n")
    next
  }
  
  # Remove NA/blank/-- treatment types
  valid <- !(is.na(df[[treat_col]]) | df[[treat_col]] == "" | df[[treat_col]] == "'--")
  df <- df[valid, ]
  
  # Update master lists
  all_samples <- unique(c(all_samples, df[[id_col]]))
  all_treatments <- unique(c(all_treatments, unique(df[[treat_col]])))
  
  # For each sample, record treatments
  for (i in seq_len(nrow(df))) {
    sid <- df[[id_col]][i]
    ttype <- df[[treat_col]][i]
    if (!sid %in% names(sample_treatments)) {
      sample_treatments[[sid]] <- character()
    }
    sample_treatments[[sid]] <- unique(c(sample_treatments[[sid]], ttype))
  }
}

# Build a mapping from sample_id to cancer type prefix
sample_to_cancer <- list()
for (dataset in names(input_files)) {
  file <- input_files[[dataset]]
  df <- fread(file, data.table = FALSE, check.names = FALSE, stringsAsFactors = FALSE, quote = "\"")
  id_col <- grep("^sample_id", colnames(df), value = TRUE)[1]
  if (!is.na(id_col)) {
    for (sid in unique(df[[id_col]])) {
      # Only set if not already set (in case of duplicates across datasets)
      if (is.null(sample_to_cancer[[sid]])) {
        sample_to_cancer[[sid]] <- dataset
      }
    }
  }
}

# Build summary matrix
all_samples <- sort(unique(all_samples))
all_treatments <- sort(unique(all_treatments))
summary_mat <- matrix("No", nrow = length(all_samples), ncol = length(all_treatments),
                      dimnames = list(all_samples, all_treatments))

for (sid in names(sample_treatments)) {
  for (ttype in sample_treatments[[sid]]) {
    summary_mat[sid, ttype] <- "Yes"
  }
}

# Convert to data.frame for output and add cancer_type column at the end
summary_df <- data.frame(sample_id = rownames(summary_mat), summary_mat, row.names = NULL, check.names = FALSE)
summary_df$cancer_type <- unname(sample_to_cancer[summary_df$sample_id])

# Use quote=TRUE to ensure treatment types with commas are quoted in the output
fwrite(summary_df, output_file, row.names = FALSE, quote = TRUE)
cat("Summary written to", output_file, "\n")

