# ===============================================================================
# Clinical Correlation Filtering Script
# ===============================================================================
# 
# Purpose: Filters clinical correlation results based on correlation strength and
#          statistical significance criteria.
#
# Input Requirements:
#   1. Clinical correlation results CSV files from clinical_correlations.R
#      (must contain columns: Signature, Correlation, P_value, P_value_adjusted)
#
# Output:
#   For each input file:
#     - Filtered correlation results CSV: Only correlations meeting specified criteria
#
# Filtering Criteria:
#   - Default: |correlation| > 0.3 AND p-value < 0.05
#   - Customizable via positional arguments
#
# Usage:
#   Rscript clinical_correlation_filter.R <cor_threshold> <p_threshold> <input_file1> <input_file2> ...
#   
#   Output files will be named: <input_basename>_filtered.csv
#
# ===============================================================================

suppressPackageStartupMessages(library(data.table))

filter_clinical_correlations <- function(input_files, cor_threshold = 0.35, p_threshold = 1) {
  
  cat("Starting correlation filtering with criteria:\n")
  cat("  Correlation threshold: |r| >", cor_threshold, "\n")
  cat("  P-value threshold: p <", p_threshold, "\n\n")
  
  for (input_file in input_files) {
    cat(rep("=", 60), "\n")
    cat("Processing file:", input_file, "\n")
    cat(rep("=", 60), "\n")
    
    # Check if file exists
    if (!file.exists(input_file)) {
      cat("Warning: File not found:", input_file, ". Skipping...\n")
      next
    }
    
    # Load correlation results
    tryCatch({
      correlation_data <- fread(input_file, data.table = FALSE, check.names = FALSE)
    }, error = function(e) {
      cat("Error reading file:", input_file, "\n")
      cat("Error message:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(correlation_data) || nrow(correlation_data) == 0) {
      cat("Warning: No data found in file:", input_file, ". Skipping...\n")
      next
    }
    
    # Validate required columns
    required_cols <- c("Signature", "Correlation", "P_value")
    missing_cols <- setdiff(required_cols, colnames(correlation_data))
    
    if (length(missing_cols) > 0) {
      cat("Error: Missing required columns in", input_file, ":", paste(missing_cols, collapse = ", "), "\n")
      cat("Available columns:", paste(colnames(correlation_data), collapse = ", "), "\n")
      next
    }
    
    # Display input data summary
    cat("Input data summary:\n")
    cat("  Total correlations:", nrow(correlation_data), "\n")
    cat("  Correlation range:", round(min(correlation_data$Correlation, na.rm = TRUE), 3), 
        "to", round(max(correlation_data$Correlation, na.rm = TRUE), 3), "\n")
    cat("  P-value range:", round(min(correlation_data$P_value, na.rm = TRUE), 6), 
        "to", round(max(correlation_data$P_value, na.rm = TRUE), 6), "\n")
    
    # Apply filtering criteria
    filtered_results <- correlation_data[
      abs(correlation_data$Correlation) > cor_threshold & 
      correlation_data$P_value < p_threshold, 
    ]
    
    # Generate output filename
    input_basename <- tools::file_path_sans_ext(basename(input_file))
    output_dir <- dirname(input_file)
    output_file <- file.path(output_dir, paste0(input_basename, "_filtered.csv"))
    
    # Save filtered results
    write.csv(filtered_results, output_file, row.names = FALSE, quote = FALSE)
    
    # Display filtering results
    cat("\nFiltering results:\n")
    cat("  Input correlations:", nrow(correlation_data), "\n")
    cat("  Correlations meeting criteria:", nrow(filtered_results), "\n")
    cat("  Filtering efficiency:", round(nrow(filtered_results) / nrow(correlation_data) * 100, 1), "%\n")
    cat("  Output saved to:", output_file, "\n")
    
    # Display filtered results preview if any meet the criteria
    if (nrow(filtered_results) > 0) {
      cat("\nFiltered Results Preview:\n")
      print(head(filtered_results, 10))
      
      # Additional statistics for filtered results
      cat("\nFiltered Results Statistics:\n")
      cat("  Strongest correlation:", round(max(abs(filtered_results$Correlation)), 3), "\n")
      cat("  Most significant p-value:", format(min(filtered_results$P_value), scientific = TRUE), "\n")
      
      if ("P_value_adjusted" %in% colnames(filtered_results)) {
        significant_adjusted <- sum(filtered_results$P_value_adjusted < 0.05, na.rm = TRUE)
        cat("  Significant after FDR correction (p_adj < 0.05):", significant_adjusted, "\n")
      }
      
    } else {
      cat("\nNo correlations met the filtering criteria (|r| >", cor_threshold, 
          "AND p <", p_threshold, ")\n")
    }
    
    cat("Completed filtering for:", input_file, "\n\n")
  }
  
  cat("All files processed successfully.\n")
}

# ===============================================================================
# Command Line Interface
# ===============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript clinical_correlation_filter.R <cor_threshold> <p_threshold> <input_file1> <input_file2> ...")
}

# Parse positional arguments
cor_threshold <- as.numeric(args[1])
p_threshold <- as.numeric(args[2])
input_files <- args[3:length(args)]

# Validate arguments
if (is.na(cor_threshold) || cor_threshold < 0 || cor_threshold > 1) {
  stop("Error: Correlation threshold must be a number between 0 and 1")
}

if (is.na(p_threshold) || p_threshold < 0 || p_threshold > 1) {
  stop("Error: P-value threshold must be a number between 0 and 1")
}

if (length(input_files) == 0) {
  stop("Error: No input files specified")
}

cat("Clinical Correlation Filtering\n")
cat("Input files:", length(input_files), "\n")
for (i in seq_along(input_files)) {
  cat("  ", i, ":", input_files[i], "\n")
}
cat("Correlation threshold: |r| >", cor_threshold, "\n")
cat("P-value threshold: p <", p_threshold, "\n\n")

# Execute filtering
filter_clinical_correlations(input_files, cor_threshold, p_threshold)

# ===============================================================================
# Example Usage:
# Rscript clinical_correlation_filter.R 0.3 0.05 HNSCC_KEGG_clinical_correlations.csv HNSCC_HALLMARK_clinical_correlations.csv
# Rscript clinical_correlation_filter.R 0.5 0.01 *.csv
# 
# This will create filtered files like:
# - HNSCC_KEGG_clinical_correlations_filtered.csv
# - HNSCC_HALLMARK_clinical_correlations_filtered.csv
# ===============================================================================
