# ===============================
# combine_datasets.R
# -------------------------------
# Purpose:
#   Combines clinical and genomic datasets across multiple cancer types for pan-cancer analysis.
#   Removes cancer type suffixes from clinical column headers to enable proper column matching.
#   Outputs 5 combined CSV files: 1 clinical file and 4 genomic files (one per pathway group).
#
# Inputs:
#   1. input_dir: Directory containing all the harmonized CSV files
#   2. output_dir: Directory to save combined files
#
# File Naming Convention:
#   Clinical files: {cancer_type}_harmonized_clinical.csv
#   Genomic files: {cancer_type}_harmonized_{pathway}_genomics.csv
#   Where pathway = kegg, hallmark, reactome, biocarta
#
# Output:
#   - COMBINED_harmonized_clinical.csv: Combined clinical data with cancer_type column
#   - COMBINED_harmonized_kegg_genomics.csv: Combined KEGG pathway data
#   - COMBINED_harmonized_hallmark_genomics.csv: Combined Hallmark pathway data  
#   - COMBINED_harmonized_reactome_genomics.csv: Combined Reactome pathway data
#   - COMBINED_harmonized_biocarta_genomics.csv: Combined BioCarta pathway data
#
# Usage:
#   Rscript combine_datasets.R "<input_dir>" "<output_dir>"
# ===============================

library(data.table)

# ---- FUNCTION: FIND FILES BY PATTERN ----
find_files_by_pattern <- function(input_dir, pattern) {
  # Find all files matching the pattern in the input directory and subdirectories
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) {
    cat("Warning: No files found matching pattern:", pattern, "in directory:", input_dir, "\n")
    cat("Searching recursively in subdirectories...\n")
  }
  return(files)
}

# ---- FUNCTION: EXTRACT CANCER TYPE FROM FILE NAME ----
extract_cancer_type_from_filename <- function(filename) {
  # Extract just the filename from the full path
  basename_file <- basename(filename)
  
  # Split by underscore and take the first part (cancer type)
  parts <- strsplit(basename_file, "_")[[1]]
  return(parts[1])
}

# ---- FUNCTION: REMOVE CANCER TYPE SUFFIX FROM COLUMN NAMES ----
remove_cancer_suffix <- function(column_names) {
  # Remove everything after the last underscore (including the underscore)
  sapply(column_names, function(col_name) {
    parts <- strsplit(col_name, "_")[[1]]
    if (length(parts) > 1) {
      # Remove the last part (cancer type suffix)
      paste(parts[-length(parts)], collapse = "_")
    } else {
      col_name  # No underscore found, return as is
    }
  }, USE.NAMES = FALSE)
}

# ---- FUNCTION: COMBINE DATASETS ACROSS CANCER TYPES ----
combine_datasets_across_cancers <- function(input_dir, output_dir) {
  # input_dir: directory containing all harmonized CSV files
  # output_dir: directory to save combined files
  
  cat("Starting dataset combination across cancer types...\n")
  cat("Input directory:", input_dir, "\n")
  cat("Output directory:", output_dir, "\n")
  
  # Find all clinical files
  clinical_files <- find_files_by_pattern(input_dir, ".*_harmonized_clinical\\.csv$")
  
  # Find genomic files for each pathway group
  kegg_files <- find_files_by_pattern(input_dir, ".*_harmonized_kegg_genomics\\.csv$")
  hallmark_files <- find_files_by_pattern(input_dir, ".*_harmonized_hallmark_genomics\\.csv$")
  reactome_files <- find_files_by_pattern(input_dir, ".*_harmonized_reactome_genomics\\.csv$")
  biocarta_files <- find_files_by_pattern(input_dir, ".*_harmonized_biocarta_genomics\\.csv$")
  
  # Create genomic files list
  genomic_files_list <- list(
    kegg = kegg_files,
    hallmark = hallmark_files,
    reactome = reactome_files,
    biocarta = biocarta_files
  )
  
  # Extract cancer types from clinical file names
  cancer_types <- sapply(clinical_files, extract_cancer_type_from_filename)
  cat("Detected cancer types:", paste(cancer_types, collapse = ", "), "\n")
  cat("Found", length(clinical_files), "clinical files\n")
  cat("Found", length(kegg_files), "KEGG files\n")
  cat("Found", length(hallmark_files), "Hallmark files\n")
  cat("Found", length(reactome_files), "Reactome files\n")
  cat("Found", length(biocarta_files), "BioCarta files\n")
  
  # ---- COMBINE CLINICAL DATA ----
  combined_clinical <- data.frame()
  
  for (i in seq_along(clinical_files)) {
    cancer_type <- cancer_types[i]
    clinical_file <- clinical_files[i]
    
    if (!file.exists(clinical_file)) {
      cat("Warning: Clinical file does not exist:", clinical_file, "\n")
      next
    }
    
    cat("Processing clinical data for:", cancer_type, "\n")
    clinical_data <- fread(clinical_file, data.table = FALSE)
    
    # Remove cancer type suffix from column names
    original_names <- colnames(clinical_data)
    cleaned_names <- remove_cancer_suffix(original_names)
    colnames(clinical_data) <- cleaned_names
    
    cat("  Removed suffixes from", length(original_names), "columns\n")
    cat("  Example: '", original_names[1], "' -> '", cleaned_names[1], "'\n")
    
    # Add cancer type column for stratification in Cox models
    clinical_data$cancer_type <- cancer_type
    
    # If this is the first dataset, initialize combined_clinical
    if (nrow(combined_clinical) == 0) {
      combined_clinical <- clinical_data
    } else {
      # Combine datasets, filling missing columns with NA
      all_cols <- unique(c(colnames(combined_clinical), colnames(clinical_data)))
      
      # Add missing columns to both datasets
      for (col in all_cols) {
        if (!col %in% colnames(combined_clinical)) {
          combined_clinical[[col]] <- NA
        }
        if (!col %in% colnames(clinical_data)) {
          clinical_data[[col]] <- NA
        }
      }
      
      # Reorder columns to match
      clinical_data <- clinical_data[, colnames(combined_clinical), drop = FALSE]
      combined_clinical <- rbind(combined_clinical, clinical_data)
    }
    
    cat("  Added", nrow(clinical_data), "samples from", cancer_type, "\n")
  }
  
  # Save combined clinical data
  clinical_output_file <- file.path(output_dir, "COMBINED_harmonized_clinical.csv")
  write.csv(combined_clinical, clinical_output_file, row.names = FALSE)
  cat("Combined clinical data saved to:", clinical_output_file, "\n")
  cat("Total clinical samples:", nrow(combined_clinical), "\n")
  
  # ---- COMBINE GENOMIC DATA FOR EACH PATHWAY GROUP ----
  genomic_output_files <- list()
  
  for (pathway_group in names(genomic_files_list)) {
    cat("\nProcessing pathway group:", pathway_group, "\n")
    genomic_files <- genomic_files_list[[pathway_group]]
    combined_genomic <- data.frame()
    
    for (i in seq_along(genomic_files)) {
      genomic_file <- genomic_files[i]
      cancer_type <- extract_cancer_type_from_filename(genomic_file)
      
      if (!file.exists(genomic_file)) {
        cat("Warning: Genomic file does not exist:", genomic_file, "\n")
        next
      }
      
      cat("  Processing", pathway_group, "data for:", cancer_type, "\n")
      genomic_data <- fread(genomic_file, data.table = FALSE)
      
      # Set row names from first column (sample IDs) and remove the ID column
      if (ncol(genomic_data) > 1) {
        rownames(genomic_data) <- genomic_data[, 1]
        genomic_data <- genomic_data[, -1, drop = FALSE]
      }
      
      # No need to add cancer type suffix since sample IDs should be unique across cancer types
      
      # If this is the first dataset for this pathway group, initialize combined_genomic
      if (nrow(combined_genomic) == 0) {
        combined_genomic <- genomic_data
      } else {
        # Get all unique pathway signatures (columns) across datasets
        all_pathways <- unique(c(colnames(combined_genomic), colnames(genomic_data)))
        
        # Add missing pathways as new columns to existing combined data (fill with 0)
        for (pathway in all_pathways) {
          if (!pathway %in% colnames(combined_genomic)) {
            combined_genomic[[pathway]] <- 0
            cat("    Added new pathway column to existing data:", pathway, "\n")
          }
          if (!pathway %in% colnames(genomic_data)) {
            genomic_data[[pathway]] <- 0
            cat("    Added missing pathway column to", cancer_type, "data:", pathway, "\n")
          }
        }
        
        # Reorder columns to match between datasets
        genomic_data <- genomic_data[, colnames(combined_genomic), drop = FALSE]
        
        # Add new samples (rows) to the combined dataset
        combined_genomic <- rbind(combined_genomic, genomic_data)
      }
      
      cat("    Added", nrow(genomic_data), "samples from", cancer_type, "\n")
      cat("    Current combined dataset:", nrow(combined_genomic), "samples x", ncol(combined_genomic), "pathways\n")
    }
    
    # Add sample IDs as first column for output
    combined_genomic_output <- data.frame(SampleID = rownames(combined_genomic), combined_genomic)
    
    # Save combined genomic data for this pathway group
    genomic_output_file <- file.path(output_dir, paste0("COMBINED_harmonized_", tolower(pathway_group), "_genomics.csv"))
    write.csv(combined_genomic_output, genomic_output_file, row.names = FALSE)
    cat("Combined", pathway_group, "genomic data saved to:", genomic_output_file, "\n")
    cat("Final dimensions:", nrow(combined_genomic_output), "samples x", ncol(combined_genomic_output)-1, "pathways\n")
    
    genomic_output_files[[pathway_group]] <- genomic_output_file
  }
  
  cat("\n=== DATASET COMBINATION SUMMARY ===\n")
  cat("Clinical file:", clinical_output_file, "\n")
  for (pathway_group in names(genomic_output_files)) {
    cat(pathway_group, "genomic file:", genomic_output_files[[pathway_group]], "\n")
  }
  cat("Dataset combination completed successfully!\n")
  
  return(list(
    clinical_file = clinical_output_file,
    genomic_files = genomic_output_files
  ))
}

# ---- COMMAND LINE ARGUMENTS ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript combine_datasets.R '<input_dir>' '<output_dir>'")
}

input_dir <- args[1]
output_dir <- args[2]

# Validate input directory exists
if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}

# Create output directory
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Combine datasets
result <- combine_datasets_across_cancers(
  input_dir = input_dir,
  output_dir = output_dir
)

cat("\n=== COMBINATION COMPLETE ===\n")
cat("All 5 output files have been generated successfully!\n")
