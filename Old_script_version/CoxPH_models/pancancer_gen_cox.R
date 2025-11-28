# ===============================
# pancancer_gen_cox.R
# -------------------------------
# Purpose:
#   Performs univariate Cox proportional hazards regression for each genomic signature across all cancer types.
#   Models are adjusted for cancer type, age, gender, and tumor stage.
#   Uses binarized (median split) features and outputs results with FDR-corrected p-values.
#
# Inputs:
#   1. genomic_features_file: Path to combined genomics matrix (samples as rows, features as columns)
#   2. clinical_data_file: Path to combined harmonized clinical data (CSV)
#   3. output_dir: Directory to write results (CSV with _binary suffix)
#   4. output_prefix: Prefix for output file naming
#
# Output:
#   CSV file containing Cox regression results for each genomic feature:
#   - Feature name
#   - Hazard Ratio (HR)
#   - 95% Confidence Interval bounds
#   - P-value
#   - Concordance index (C-index)
#   - FDR-corrected p-value
#
# Main Steps:
#   - Loads combined genomics and clinical data
#   - Binarizes features by median split
#   - For each feature, fits univariate Cox model adjusted for cancer type, age, gender, and stage
#   - Outputs results with FDR correction for multiple testing
# ===============================

library(survival)
library(survminer)
library(data.table)

# ---- COMMAND LINE ARGUMENTS ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript pancancer_gen_cox.R <genomic_features_file> <clinical_data_file> <output_dir> <output_prefix>")
}

genomic_features_file <- args[1]
clinical_data_file <- args[2]
output_dir <- args[3]
output_prefix <- args[4]

# Ensure output directory exists
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Starting pan-cancer univariate Cox modeling...\n")
cat("Genomic file:", genomic_features_file, "\n")
cat("Clinical file:", clinical_data_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Output prefix:", output_prefix, "\n")

# ---- LOAD DATA ----
genomic_features <- read.csv(genomic_features_file, row.names = 1)
clinical_data <- read.csv(clinical_data_file)

cat("Loaded genomic data:", nrow(genomic_features), "samples x", ncol(genomic_features), "features\n")
cat("Loaded clinical data:", nrow(clinical_data), "samples\n")

# ---- ENSURE CLINICAL DATA CONTAINS REQUIRED COLUMNS ----
# Since clinical data suffixes have been removed, use standard column names
required_columns <- c("OS_days", "OS_event", "diagnoses.age_at_diagnosis", "demographic.gender", "diagnoses.ajcc_pathologic_stage", "cancer_type")

if (!all(required_columns %in% colnames(clinical_data))) {
  missing_cols <- required_columns[!required_columns %in% colnames(clinical_data)]
  stop(paste0("Clinical data is missing required columns: ", paste(missing_cols, collapse=", ")))
}

# ---- MERGE DATA ----
# Use the standard sample ID column name (without suffix)
sample_id_col <- "cases.submitter_id"

if (!sample_id_col %in% colnames(clinical_data)) {
  stop(paste0("Clinical data is missing sample ID column: ", sample_id_col))
}

# Filter genomic features to only include samples present in clinical data
filtered_genomic_features <- genomic_features[rownames(genomic_features) %in% clinical_data[[sample_id_col]], ]
cat(sprintf("Filtered genomic features based on clinical data. Remaining samples: %d\n", nrow(filtered_genomic_features)))

# Binarize features by median for binary analysis
binarized_features <- filtered_genomic_features
for (feature in colnames(filtered_genomic_features)) {
  median_value <- median(filtered_genomic_features[[feature]], na.rm = TRUE)
  binarized_features[[feature]] <- ifelse(filtered_genomic_features[[feature]] > median_value, 1, 0)
}
cat("Genomic features have been binarized based on median expression scores.\n")

# Merge using the standard sample ID
merged_data <- merge(clinical_data, binarized_features, by.x = sample_id_col, by.y = "row.names", all.x = TRUE)
rownames(merged_data) <- merged_data[[sample_id_col]]

# Remove samples with missing survival data
merged_data <- merged_data[!is.na(merged_data$OS_days) & !is.na(merged_data$OS_event), ]
cat("Final merged dataset:", nrow(merged_data), "samples\n")

# ---- PREPARE SURVIVAL OBJECT ----
# Create survival object with both time and event status
surv_object <- Surv(time = merged_data$OS_days, event = merged_data$OS_event)

# ---- COX PROPORTIONAL HAZARDS MODEL (BINARY) ----
cat("Running univariate Cox models (binary features)...\n")
results_binary <- data.frame(
  Signature = character(),
  HR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p_value = numeric(),
  C_index = numeric(),
  FDR = numeric(),
  stringsAsFactors = FALSE
)

for (feature in colnames(binarized_features)) {
  if (any(is.na(merged_data[[feature]]))) {
    cat("Skipping binary feature due to missing values:", feature, "\n")
    next
  }
  
  # Univariate Cox model controlling for cancer type, age, gender, and stage
  cox_model <- coxph(surv_object ~ merged_data[[feature]] + 
                     merged_data$diagnoses.age_at_diagnosis + 
                     merged_data$demographic.gender + 
                     merged_data$diagnoses.ajcc_pathologic_stage + 
                     factor(merged_data$cancer_type), 
                     data = merged_data)
  
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
    C_index = c_index,
    FDR = NA  # Will be calculated after all features are processed
  ))
}

# FDR correction for binary model
results_binary$FDR <- p.adjust(results_binary$p_value, method = "fdr")

# ---- WRITE OUTPUTS ----
binary_output_file <- file.path(output_dir, paste0(output_prefix, "_binary.csv"))

if (nrow(results_binary) == 0) {
  cat("No binary Cox model results to write. Writing empty file.\n")
  write.csv(data.frame(), binary_output_file, row.names = FALSE)
} else {
  write.csv(results_binary, binary_output_file, row.names = FALSE)
  cat("Binary Cox model results saved to:", binary_output_file, "\n")
  cat("Number of significant binary features (p < 0.05):", sum(results_binary$p_value < 0.05), "\n")
  cat("Number of significant binary features (FDR < 0.05):", sum(results_binary$FDR < 0.05), "\n")
}

cat("Pan-cancer univariate Cox modeling completed successfully!\n")