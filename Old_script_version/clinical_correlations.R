# ===============================================================================
# Clinical Correlation Analysis Script
# ===============================================================================
# 
# Purpose: Performs Spearman correlation analysis between multiple genomic signature 
#          datasets and overall survival days (OS_days) from clinical data.
#
# Input Requirements:
#   1. Multiple genomic signatures files: CSV files with samples as rows, genomic signatures as columns
#   2. Clinical data file: CSV with samples as rows, clinical outcomes as columns
#      (must contain 'OS_days' column)
#
# Output:
#   For each genomic dataset:
#     1. Full correlation results CSV: All correlations with p-values and adjusted p-values
#     (Use clinical_correlation_filter.R for filtering by correlation strength and significance)
#
# Analysis Method:
#   - Uses Spearman rank correlation (non-parametric)
#   - Applies False Discovery Rate (FDR) multiple testing correction
#   - Filters out signatures with insufficient variance or data points
#   - Matches sample IDs between genomic and clinical datasets
#   - Processes multiple genomic datasets against single clinical dataset
#
# Usage:
#   Rscript clinical_correlations.R <clinical_file> <output_prefix> <genomic_file1> <genomic_file2> ... <genomic_fileN>
#   
#   Output files will be named: <output_prefix>_<source group>_clinical_correlations.csv
#
# ===============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

perform_clinical_correlations <- function(clinical_file, genomic_files, output_prefix) {
  # =============================================================================
  # Clinical Data Loading and Preprocessing
  # =============================================================================
  
  cat("Loading clinical dataset from:", clinical_file, "\n")
  clinical_data <- fread(clinical_file, data.table = FALSE, check.names = FALSE)
  
  # Set sample IDs as rownames and remove the ID column
  rownames(clinical_data) <- clinical_data[[1]]
  clinical_data[[1]] <- NULL
  
  # =============================================================================
  # Clinical Data Validation
  # =============================================================================
  
  # Verify that OS_days column exists in clinical data
  if (!"OS_days" %in% colnames(clinical_data)) {
    stop("Error: 'OS_days' column not found in clinical data. Available columns: ", 
         paste(colnames(clinical_data), collapse = ", "))
  }
  
  # Extract OS_days and identify samples with valid survival data
  os_days <- clinical_data$OS_days
  valid_clinical_indices <- !is.na(os_days) & is.finite(os_days)
  
  if (sum(valid_clinical_indices) == 0) {
    stop("Error: No valid OS_days values found (all are NA or infinite)")
  }
  
  cat("Number of clinical samples with valid OS_days:", sum(valid_clinical_indices), "\n")
  clinical_samples_valid <- rownames(clinical_data)[valid_clinical_indices]
  os_days_valid <- os_days[valid_clinical_indices]
  names(os_days_valid) <- clinical_samples_valid
  
  # =============================================================================
  # Process Each Genomic Dataset
  # =============================================================================
  
  for (genomic_file in genomic_files) {
    cat("\n" , rep("=", 60), "\n")
    cat("Processing genomic dataset:", genomic_file, "\n")
    cat(rep("=", 60), "\n")
    
    # Load genomic signature dataset
    genomic_data <- fread(genomic_file, data.table = FALSE, check.names = FALSE)
    
    # Debug: Check data structure after loading
    cat("Initial genomic data dimensions:", nrow(genomic_data), "x", ncol(genomic_data), "\n")
    if (ncol(genomic_data) > 0) {
      n_cols_to_show <- min(5, ncol(genomic_data))
      cat("First few column names:", paste(colnames(genomic_data)[seq_len(n_cols_to_show)], collapse = ", "), "\n")
    }
    cat("First column class:", class(genomic_data[[1]]), "\n")
    
    # Set sample IDs as rownames and remove the ID column
    rownames(genomic_data) <- genomic_data[[1]]
    genomic_data[[1]] <- NULL
    
    # Debug: Check data structure after removing ID column
    cat("Genomic data dimensions after ID removal:", nrow(genomic_data), "x", ncol(genomic_data), "\n")
    if (ncol(genomic_data) > 0) {
      n_cols_to_check <- min(5, ncol(genomic_data))
      cat("Column classes (first 5):", paste(sapply(genomic_data[seq_len(n_cols_to_check)], class), collapse = ", "), "\n")
    }
    
    # Check if all remaining columns can be converted to numeric
    numeric_check <- sapply(genomic_data, function(x) {
      test_numeric <- suppressWarnings(as.numeric(as.character(x)))
      !all(is.na(test_numeric)) || all(is.na(x))  # TRUE if convertible to numeric OR already all NA
    })
    cat("Number of columns that can be converted to numeric:", sum(numeric_check), "out of", ncol(genomic_data), "\n")
    
    if (sum(numeric_check) < ncol(genomic_data)) {
      cat("Warning: Some columns cannot be converted to numeric. Column types:\n")
      non_numeric_cols <- names(genomic_data)[!numeric_check]
      if (length(non_numeric_cols) > 0) {
        n_cols_to_show <- min(10, length(non_numeric_cols))
        cat("Non-numeric columns:", paste(non_numeric_cols[seq_len(n_cols_to_show)], collapse = ", "), "\n")
      }
    }
    
    # Extract dataset name from filename for labeling
    dataset_name <- tools::file_path_sans_ext(basename(genomic_file))
    
    # Extract source group (KEGG, HALLMARK, REACTOME, BIOCARTA, radiomics) from filename
    source_group <- NA
    if (grepl("KEGG", dataset_name, ignore.case = TRUE)) {
      source_group <- "KEGG"
    } else if (grepl("HALLMARK", dataset_name, ignore.case = TRUE)) {
      source_group <- "HALLMARK"
    } else if (grepl("REACTOME", dataset_name, ignore.case = TRUE)) {
      source_group <- "REACTOME"
    } else if (grepl("BIOCARTA", dataset_name, ignore.case = TRUE)) {
      source_group <- "BIOCARTA"
    } else if (grepl("radiomics", dataset_name, ignore.case = TRUE)) {
      source_group <- "radiomics"
    }
    
    if (is.na(source_group)) {
      cat("Warning: Could not identify source group from filename:", genomic_file, "\n")
      cat("Using full dataset name for output files.\n")
      source_group <- dataset_name
    }
    
    # =============================================================================
    # Sample Matching and Quality Control
    # =============================================================================
    
    # Find overlapping samples between current genomic dataset and valid clinical samples
    common_samples <- intersect(rownames(genomic_data), clinical_samples_valid)
    cat("Number of common samples with valid clinical data:", length(common_samples), "\n")
    
    if (length(common_samples) == 0) {
      cat("Warning: No common samples found between", dataset_name, "and clinical data. Skipping...\n")
      next
    }
    
    # Subset genomic data to common samples with valid clinical data
    genomic_subset <- genomic_data[common_samples, , drop = FALSE]
    os_days_subset <- os_days_valid[common_samples]
    
    # First, filter out any non-numeric columns
    numeric_columns <- sapply(genomic_subset, function(x) {
      test_numeric <- suppressWarnings(as.numeric(as.character(x)))
      !all(is.na(test_numeric)) || all(is.na(x))  # TRUE if convertible to numeric OR already all NA
    })
    
    if (sum(numeric_columns) == 0) {
      cat("Warning: No numeric columns found in", dataset_name, ". Skipping...\n")
      next
    }
    
    if (sum(numeric_columns) < ncol(genomic_subset)) {
      cat("Filtering out", ncol(genomic_subset) - sum(numeric_columns), "non-numeric columns\n")
      genomic_subset <- genomic_subset[, numeric_columns, drop = FALSE]
    }
    
    # Remove genomic signatures that have all missing values or no variance
    cat("Checking genomic signatures for validity...\n")
    cat("Genomic subset dimensions:", nrow(genomic_subset), "x", ncol(genomic_subset), "\n")
    
    if (ncol(genomic_subset) == 0) {
      cat("Warning: No genomic signatures found in", dataset_name, ". Skipping...\n")
      next
    }
    
    valid_signatures <- apply(genomic_subset, 2, function(x) {
      # Convert to numeric first, handling any non-numeric values
      x_numeric <- suppressWarnings(as.numeric(x))
      # Check if all values are NA after conversion
      if (all(is.na(x_numeric))) {
        return(FALSE)
      }
      # Check if there's sufficient variance
      variance <- var(x_numeric, na.rm = TRUE)
      if (is.na(variance) || variance == 0) {
        return(FALSE)
      }
      return(TRUE)
    })
    
    # Debug: Print information about valid signatures
    cat("Total signatures before filtering:", ncol(genomic_subset), "\n")
    cat("Valid signatures after filtering:", sum(valid_signatures), "\n")
    cat("Valid signatures logical vector length:", length(valid_signatures), "\n")
    cat("Names of valid signatures (first 10):", paste(names(valid_signatures)[valid_signatures][1:min(10, sum(valid_signatures))], collapse = ", "), "\n")
    
    # Check if any signatures are valid
    if (sum(valid_signatures) == 0) {
      cat("Warning: No valid genomic signatures found in", dataset_name, ". Skipping...\n")
      next
    }
    
    # Additional safety check - ensure we have valid column indices
    if (length(valid_signatures) != ncol(genomic_subset)) {
      cat("Error: Mismatch between valid_signatures length (", length(valid_signatures), 
          ") and genomic_subset columns (", ncol(genomic_subset), "). Skipping...\n")
      next
    }
    
    genomic_final <- genomic_subset[, valid_signatures, drop = FALSE]
    cat("Number of valid genomic signatures in", dataset_name, ":", ncol(genomic_final), "\n")
    
    if (ncol(genomic_final) == 0) {
      cat("Warning: No valid genomic signatures found in", dataset_name, ". Skipping...\n")
      next
    }
  
    # =============================================================================
    # Spearman Correlation Analysis for Current Dataset
    # =============================================================================
    
    cat("Performing Spearman correlation analysis for", dataset_name, "...\n")
    
    # Initialize results data frame for current dataset
    correlation_results <- data.frame(
      Signature = character(),
      Sample_Count = numeric(),
      Mean_OS_Days = numeric(),
      OS_Days_Range = character(),
      Correlation = numeric(),
      P_value = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Iterate through each genomic signature in current dataset
    for (i in seq_len(ncol(genomic_final))) {
      signature_name <- colnames(genomic_final)[i]
      signature_values <- genomic_final[, i]
      
      # Identify samples with valid data for both signature and OS_days
      valid_idx <- which(!is.na(signature_values) & !is.na(os_days_subset) & 
                         signature_values != "" & os_days_subset != "")
      
      # Convert to numeric vectors for correlation analysis
      x_valid <- as.numeric(signature_values[valid_idx])
      y_valid <- as.numeric(os_days_subset[valid_idx])
      
      # Quality control: check for sufficient data points and variance
      if (length(unique(x_valid)) < 2 || length(unique(y_valid)) < 2 || length(x_valid) < 3) {
        cat("Warning: Skipping", signature_name, "in", dataset_name, "- insufficient data points or no variance\n")
        next
      }
      
      # Perform Spearman rank correlation test
      cor_test <- suppressWarnings(cor.test(x_valid, y_valid, method = "spearman"))
      rho <- cor_test$estimate    # Correlation coefficient
      pval <- cor_test$p.value    # Statistical significance
      
      # Calculate summary statistics for the OS_days values used in this correlation
      sample_count <- length(y_valid)
      mean_os_days <- round(mean(y_valid, na.rm = TRUE), 2)
      min_os_days <- round(min(y_valid, na.rm = TRUE), 2)
      max_os_days <- round(max(y_valid, na.rm = TRUE), 2)
      os_range <- paste0(min_os_days, "-", max_os_days)
      
      # Store one correlation result per signature
      correlation_results <- rbind(correlation_results, data.frame(
        Signature = signature_name,
        Sample_Count = sample_count,
        Mean_OS_Days = mean_os_days,
        OS_Days_Range = os_range,
        Correlation = rho,
        P_value = pval,
        stringsAsFactors = FALSE
      ))
    }
    
    # =============================================================================
    # Statistical Correction and Results Processing for Current Dataset
    # =============================================================================
    
    if (nrow(correlation_results) == 0) {
      cat("Warning: No correlations computed for", dataset_name, ". Skipping output generation.\n")
      next
    }
    
    # Apply False Discovery Rate (FDR) multiple testing correction for current dataset
    correlation_results$P_value_adjusted <- p.adjust(correlation_results$P_value, method = "fdr")
    
    # Sort results by absolute correlation strength (strongest first)
    correlation_results <- correlation_results[order(abs(correlation_results$Correlation), decreasing = TRUE), ]
    
    # =============================================================================
    # Results Output and Filtering for Current Dataset
    # =============================================================================
    
    # Generate output filenames for current dataset using new naming convention
    # Extract prefix from output_prefix (remove the trailing part after last underscore)
    prefix_parts <- strsplit(output_prefix, "_")[[1]]
    if (length(prefix_parts) >= 1) {
      cancer_prefix <- prefix_parts[1]  # Take the first part as cancer type (e.g., HNSCC)
    } else {
      cancer_prefix <- output_prefix
    }
    
    dataset_output_file <- paste0(cancer_prefix, "_", source_group, "_clinical_correlations.csv")
    
    # Display debug information: top correlations for current dataset
    cat("Top 10 correlation results for", dataset_name, ":\n")
    print(head(correlation_results, 10))
    
    # Save complete correlation results for current dataset
    write.csv(correlation_results, dataset_output_file, row.names = FALSE, quote = FALSE)
    cat("Full correlation results for", dataset_name, "saved to:", dataset_output_file, "\n")
    
    # =============================================================================
    # Summary Statistics for Current Dataset
    # =============================================================================
    
    cat("\nSummary Statistics for", dataset_name, ":\n")
    cat("  Total signatures analyzed:", nrow(correlation_results), "\n")
    cat("  Significant correlations (adjusted p < 0.05):", sum(correlation_results$P_value_adjusted < 0.05), "\n")
    cat("  Strong correlations (|r| > 0.3):", sum(abs(correlation_results$Correlation) > 0.3), "\n")
    cat("  Very strong correlations (|r| > 0.7):", sum(abs(correlation_results$Correlation) > 0.7), "\n")
    if (nrow(correlation_results) > 0) {
      cat("  Range of correlations:", round(min(correlation_results$Correlation), 3), "to", 
          round(max(correlation_results$Correlation), 3), "\n")
      cat("  Sample count range:", min(correlation_results$Sample_Count), "to", 
          max(correlation_results$Sample_Count), "samples per signature\n")
    }
    
    # Display top correlations preview
    if (nrow(correlation_results) > 0) {
      cat("\nTop Correlation Results for", dataset_name, ":\n")
      print(head(correlation_results, 5))
    }
    
    cat("Completed analysis for", dataset_name, ": analyzed", ncol(genomic_final), "signatures\n")
    cat("Use clinical_correlation_filter.R to filter results by correlation strength and significance.\n")
  }
  
  cat("\nAll datasets processed successfully.\n")
}


# ===============================================================================
# Command Line Interface
# ===============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript clinical_correlations.R <clinical_file> <output_prefix> <genomic_file1> <genomic_file2> ... <genomic_fileN>")
}

clinical_file <- args[1]           # Path to clinical data CSV file
output_prefix <- args[2]           # Prefix for output correlation results CSV files
genomic_files <- args[3:length(args)]  # Paths to genomic/radiomics signatures CSV files

cat("Clinical data file:", clinical_file, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("Genomic datasets to process:", length(genomic_files), "\n")
for (i in seq_along(genomic_files)) {
  cat("  ", i, ":", genomic_files[i], "\n")
}

# Execute correlation analysis
cat("\nStarting multi-dataset clinical correlation analysis...\n")
perform_clinical_correlations(clinical_file, genomic_files, output_prefix)
cat("Analysis completed.\n")

# ===============================================================================
# Example Usage:
# Rscript clinical_correlations.R clinical_data.csv correlation_results kegg_signatures.csv hallmark_signatures.csv reactome_signatures.csv biocarta_signatures.csv
# 
# This will create output files like:
# - correlation_results_kegg_signatures.csv
# - correlation_results_hallmark_signatures.csv
# - etc.
# 
# Use clinical_correlation_filter.R to filter results:
# Rscript clinical_correlation_filter.R *.csv --cor_threshold 0.7 --p_threshold 0.05
# ===============================================================================
