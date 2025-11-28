# ===============================================================================
# Pathway Enrichment Matrix Separator by Source Database
# ===============================================================================
# 
# Purpose: Separates a combined genomic enrichment matrix into individual pathway
#          source-specific matrices (KEGG, HALLMARK, REACTOME, BIOCARTA) for
#          downstream analysis and correlation studies.
#
# Description:
#   This script takes a single genomic enrichment matrix containing pathway
#   scores from multiple databases and splits it into separate matrices based
#   on pathway source. Each output matrix contains only pathways from one
#   database, enabling source-specific radiogenomic correlation analysis.
#
# Input Requirements:
#   1. Combined enrichment matrix: CSV with samples as rows, all pathway scores as columns
#   2. Pathway names should contain source identifiers (KEGG_, HALLMARK_, etc.)
#   3. First column should contain sample IDs
#
# Output:
#   Four separate enrichment matrices:
#   - {prefix}_KEGG_enrichment.csv: KEGG pathway scores only
#   - {prefix}_HALLMARK_enrichment.csv: MSigDB Hallmark pathway scores only
#   - {prefix}_REACTOME_enrichment.csv: Reactome pathway scores only
#   - {prefix}_BIOCARTA_enrichment.csv: BioCarta pathway scores only
#
# Processing Method:
#   - Identifies pathway columns by source-specific prefixes or patterns
#   - Transposes matrices to have samples as rows (standard format)
#   - Adds "SampleID" header to first column for consistency
#   - Maintains sample order across all output files
#
# Usage:
#   Rscript enrichment_separator.R <combined_enrichment_file> <output_directory> <prefix>
#
# Dependencies: Standard R libraries (utils, stats)
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================
# based on pathway source (KEGG, HALLMARK, BIOCARTA, REACTOME).
# Each output file contains pathways from one source with samples as rows and pathways as columns.
#
# Input format: CSV file with pathways as rows and samples as columns (first column = pathway names)
# Output: Four CSV files, one for each pathway source, transposed (samples x pathways)
#
# Usage: Rscript enrichment_separator.R <input_matrix.csv> <output_dir> <output_prefix>
#

suppressPackageStartupMessages({
  library(data.table)
})

cat("Starting enrichment matrix separator script...\n")

# =========================
# Parse Command Line Arguments
# =========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript enrichment_separator.R <input_matrix.csv> <output_dir> <output_prefix>")
}

input_file <- args[1]
output_dir <- args[2]
output_prefix <- args[3]

cat("Input file:", input_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Output prefix:", output_prefix, "\n")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# =========================
# Read Input Data
# =========================
cat("Reading input enrichment matrix...\n")
enrichment_matrix <- fread(input_file, data.table = FALSE)
cat("Input matrix dimensions:", dim(enrichment_matrix), "\n")

# Set row names from first column and remove it
pathway_names <- enrichment_matrix[[1]]
rownames(enrichment_matrix) <- pathway_names
enrichment_matrix <- enrichment_matrix[, -1, drop = FALSE]
cat("Matrix dimensions after removing pathway name column:", dim(enrichment_matrix), "\n")
cat("Number of pathways:", nrow(enrichment_matrix), "\n")
cat("Number of samples:", ncol(enrichment_matrix), "\n")

# =========================
# Define Pathway Sources
# =========================
sources <- c("KEGG", "HALLMARK", "BIOCARTA", "REACTOME")
cat("\nSearching for pathway sources:", paste(sources, collapse = ", "), "\n")

# =========================
# Separate Pathways by Source
# =========================
for (source in sources) {
  cat("\n=== Processing", source, "pathways ===\n")
  
  # Find pathways matching this source
  source_indices <- grepl(source, pathway_names, ignore.case = TRUE)
  num_pathways <- sum(source_indices)
  
  cat("Number of", source, "pathways found:", num_pathways, "\n")
  
  if (num_pathways == 0) {
    cat("No", source, "pathways found. Skipping...\n")
    next
  }
  
  # Extract pathways for this source
  source_matrix <- enrichment_matrix[source_indices, , drop = FALSE]
  cat("Source matrix dimensions (pathways x samples):", dim(source_matrix), "\n")
  
  # Show first few pathway names for verification
  source_pathway_names <- pathway_names[source_indices]
  cat("First few", source, "pathway names:\n")
  print(head(source_pathway_names, 3))
  
  # Transpose matrix: pathways x samples -> samples x pathways
  source_matrix_transposed <- t(source_matrix)
  source_matrix_transposed <- as.data.frame(source_matrix_transposed)
  cat("Transposed matrix dimensions (samples x pathways):", dim(source_matrix_transposed), "\n")
  
  # Ensure all values are numeric
  cat("Converting all values to numeric...\n")
  source_matrix_transposed[] <- lapply(source_matrix_transposed, function(x) as.numeric(as.character(x)))
  
  # Add SampleID column header and prepare final output
  final_output <- cbind(SampleID = rownames(source_matrix_transposed), source_matrix_transposed)
  
  # Create output filename
  output_file <- file.path(output_dir, paste0(output_prefix, "_", source, "_enrichment.csv"))
  
  # Write to CSV without row names (since SampleID is now a proper column)
  write.csv(final_output, file = output_file, row.names = FALSE)
  cat("Saved", source, "enrichment matrix to:", output_file, "\n")
  
  # Print summary statistics
  cat("Summary statistics for", source, ":\n")
  cat("  Sample names (first 3):", paste(head(rownames(source_matrix_transposed), 3), collapse = ", "), "\n")
  cat("  Pathway names (first 3):", paste(head(colnames(source_matrix_transposed), 3), collapse = ", "), "\n")
  cat("  Value range: [", round(min(source_matrix_transposed, na.rm = TRUE), 3), ", ", 
      round(max(source_matrix_transposed, na.rm = TRUE), 3), "]\n", sep = "")
  cat("  Mean value:", round(mean(as.matrix(source_matrix_transposed), na.rm = TRUE), 3), "\n")
}

# =========================
# Summary Report
# =========================
cat("\n=== SEPARATION SUMMARY ===\n")
total_pathways_processed <- 0
for (source in sources) {
  source_indices <- grepl(source, pathway_names, ignore.case = TRUE)
  num_pathways <- sum(source_indices)
  total_pathways_processed <- total_pathways_processed + num_pathways
  cat(sprintf("%-10s: %4d pathways\n", source, num_pathways))
}
cat(sprintf("%-10s: %4d pathways\n", "TOTAL", total_pathways_processed))
cat(sprintf("%-10s: %4d pathways\n", "INPUT", length(pathway_names)))

# Check for unclassified pathways
unclassified <- length(pathway_names) - total_pathways_processed
if (unclassified > 0) {
  cat("\nWARNING:", unclassified, "pathways were not classified into any source!\n")
  
  # Find unclassified pathways
  all_classified <- rep(FALSE, length(pathway_names))
  for (source in sources) {
    source_indices <- grepl(source, pathway_names, ignore.case = TRUE)
    all_classified <- all_classified | source_indices
  }
  unclassified_names <- pathway_names[!all_classified]
  
  cat("Unclassified pathway examples (first 5):\n")
  print(head(unclassified_names, 5))
}

cat("\nEnrichment matrix separation complete!\n")
cat("Output files saved in:", output_dir, "\n")
