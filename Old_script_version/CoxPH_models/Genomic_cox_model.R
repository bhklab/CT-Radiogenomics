# ===============================
# Genomic_cox_model.R
# -------------------------------
# Purpose:
#   Performs Cox proportional hazards regression for each genomic feature (signature), adjusted for clinical covariates.
#   Outputs both binarized (median split) and continuous model results, including FDR-corrected p-values.
#
# Inputs:
#   1. genomic_features_file: Path to the genomics matrix (samples as rows, features as columns).
#   2. clinical_data_file: Path to harmonized clinical data (CSV).
#   3. output_dir: Directory to write results (CSV; suffixes _binary and _continuous will be added).
#   4. column_prefix: Prefix for clinical column matching.
#   5. output_prefix: Prefix for output file naming.
#
# Output:
#   CSV file containing Cox regression results for each radiomics feature:
#   - Feature name
#   - Hazard Ratio (HR)
#   - 95% Confidence Interval bounds
#   - P-value
#   - Concordance index (C-index)
#
# Main Steps:
#   - Loads and filters genomics and clinical data.
#   - Binarizes features by median for binary analysis.
#   - For each feature, fits a Cox model adjusted for age, gender, and stage.
#   - Outputs results with FDR correction for multiple testing.
# ===============================

library(survival)
library(survminer)

# ---- COMMAND LINE ARGUMENTS ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript Genomic_cox_model.R <genomic_features_file> <clinical_data_file> <output_dir> <column_prefix> <output_prefix>")
}

genomic_features_file <- args[1]
clinical_data_file <- args[2]
output_dir <- args[3]
column_prefix <- args[4]
output_prefix <- args[5]


# Ensure output_dir exists
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- LOAD DATA ----
genomic_features <- read.csv(genomic_features_file, row.names = 1)
clinical_data <- read.csv(clinical_data_file)

# ---- ENSURE CLINICAL DATA CONTAINS REQUIRED COLUMNS ----
# All columns have the prefix as a suffix, e.g., diagnoses.age_at_diagnosis_HNSCC
col_suffix <- paste0("_", column_prefix)
required_base_columns <- c("OS_days", "OS_event", "diagnoses.age_at_diagnosis", "demographic.gender", "diagnoses.ajcc_pathologic_stage")
required_columns <- paste0(required_base_columns, col_suffix)

if (!all(required_columns %in% colnames(clinical_data))) {
  stop(paste0("Clinical data is missing required columns: ", paste(required_columns, collapse=", ")))
}

# ---- MERGE DATA ----
# Use the prefixed column names for merging and modeling
# The sample ID column is also suffixed
sample_id_col <- paste0("cases.submitter_id", col_suffix)

# Filter and merge
filtered_genomic_features <- genomic_features[rownames(genomic_features) %in% clinical_data[[sample_id_col]], ]
cat(sprintf("Filtered genomic features based on clinical data. Remaining samples: %d\n", nrow(filtered_genomic_features)))

binarized_features <- filtered_genomic_features
for (feature in colnames(filtered_genomic_features)) {
  median_value <- median(filtered_genomic_features[[feature]], na.rm = TRUE)
  binarized_features[[feature]] <- ifelse(filtered_genomic_features[[feature]] > median_value, 1, 0)
}
cat("Genomic features have been binarized based on median expression scores.\n")

# Merge using the prefixed sample ID
merged_data <- merge(clinical_data, binarized_features, by.x = sample_id_col, by.y = "row.names", all.x = TRUE)
rownames(merged_data) <- merged_data[[sample_id_col]]

# Ensure no missing values in OS_days and OS_event
if (any(is.na(merged_data[[required_columns[1]]])) || any(is.na(merged_data[[required_columns[2]]]))) {
  stop("Missing values detected in OS_days or OS_event column after merging.")
}

# ---- PREPARE SURVIVAL OBJECT ----
surv_object <- Surv(time = merged_data[[required_columns[1]]], event = merged_data[[required_columns[2]]])

# ---- COX PROPORTIONAL HAZARDS MODEL (BINARY) ----
results_binary <- data.frame(
  Signature = character(),
  HR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p_value = numeric(),
  C_index = numeric(),
  stringsAsFactors = FALSE
)
for (feature in colnames(binarized_features)) {
  if (any(is.na(merged_data[[feature]]))) {
    cat("Skipping binary feature due to missing values:", feature, "\n")
    next
  }
  cox_model <- coxph(surv_object ~ merged_data[[feature]] + merged_data[[required_columns[2]]] + merged_data[[required_columns[3]]] + merged_data[[required_columns[4]]], data = merged_data)
  summary_cox <- summary(cox_model)
  hr <- summary_cox$coefficients[1, "exp(coef)"]
  ci <- summary_cox$conf.int[1, c("lower .95", "upper .95")]
  p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
  c_index <- summary_cox$concordance["C"]
  results_binary <- rbind(results_binary, data.frame(
    Signature = feature,
    HR = hr,
    CI_lower = ci[1],
    CI_upper = ci[2],
    p_value = p_value,
    C_index = c_index
  ))
}
# FDR correction for binary model
results_binary$FDR <- p.adjust(results_binary$p_value, method = "fdr")

# ---- COX PROPORTIONAL HAZARDS MODEL (CONTINUOUS) ----
results_continuous <- data.frame(
  Signature = character(),
  HR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p_value = numeric(),
  C_index = numeric(),
  stringsAsFactors = FALSE
)
for (feature in colnames(filtered_genomic_features)) {
  if (any(is.na(merged_data[[feature]]))) {
    cat("Skipping continuous feature due to missing values:", feature, "\n")
    next
  }
  cox_model <- coxph(surv_object ~ merged_data[[feature]] + merged_data[[required_columns[2]]] + merged_data[[required_columns[3]]] + merged_data[[required_columns[4]]], data = merged_data)
  summary_cox <- summary(cox_model)
  hr <- summary_cox$coefficients[1, "exp(coef)"]
  ci <- summary_cox$conf.int[1, c("lower .95", "upper .95")]
  p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
  # Correct way to extract concordance index
  c_index <- summary_cox$concordance["C"]
  results_continuous <- rbind(results_continuous, data.frame(
    Signature = feature,
    HR = hr,
    CI_lower = ci[1],
    CI_upper = ci[2],
    p_value = p_value,
    C_index = c_index
  ))
}
# FDR correction for continuous model
results_continuous$FDR <- p.adjust(results_continuous$p_value, method = "fdr")

# ---- WRITE OUTPUTS ----
if (nrow(results_binary) == 0) {
  cat("No binary Cox model results to write. Writing empty file.\n")
  write.csv(data.frame(), file.path(output_dir, paste0(output_prefix, "_binary.csv")), row.names = FALSE)
} else {
  write.csv(results_binary, file.path(output_dir, paste0(output_prefix, "_binary.csv")), row.names = FALSE)
  cat("Binary Cox model results saved to:", file.path(output_dir, paste0(output_prefix, "_binary.csv")), "\n")
}

if (nrow(results_continuous) == 0) {
  cat("No continuous Cox model results to write. Writing empty file.\n")
  write.csv(data.frame(), file.path(output_dir, paste0(output_prefix, "_continuous.csv")), row.names = FALSE)
} else {
  write.csv(results_continuous, file.path(output_dir, paste0(output_prefix, "_continuous.csv")), row.names = FALSE)
  cat("Continuous Cox model results saved to:", file.path(output_dir, paste0(output_prefix, "_continuous.csv")), "\n")
}