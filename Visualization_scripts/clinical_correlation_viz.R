# ===============================================================================
# Clinical-Guided Radiogenomic Correlation Heatmap Visualization
# ===============================================================================
# 
# Purpose: Creates heatmaps of genomic-radiomic correlations for genomic features
#          that show significant clinical correlations, enabling visualization of
#          clinically relevant radiogenomic associations.
#
# Description:
#   This script takes clinical correlation results and genomic-radiomic correlation
#   results from four pathway collections (KEGG, HALLMARK, REACTOME, BIOCARTA).
#   It identifies genomic features with significant clinical correlations and
#   visualizes their correlations with radiomic features in heatmap format.
#
# Input Requirements:
#   1. Clinical correlation files (4): CSV files containing clinical associations
#      - Required columns: "Signature" (genomic feature names)
#   2. Genomic-radiomic correlation files (4): CSV files containing radiogenomic correlations
#      - Required columns: "GenomicFeature", "RadiomicFeature", "SpearmanRho"
#   3. Files should correspond by pathway collection (same order)
#
# Output:
#   - Static heatmap plots saved as PNG files for each pathway collection
#   - Interactive heatmap plots saved as HTML files for zooming and exploration
#   - Only includes genomic features present in clinical correlation data
#   - Correlation matrices saved as CSV files (optional)
#
# Usage:
#   Simply run: Rscript clinical_correlation_viz.R
#   (Modify the hardcoded file paths in the script as needed)
#
# Dependencies: ggplot2, reshape2, pheatmap, plotly, heatmaply
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(plotly)
  library(heatmaply)
})

# =========================
# Hardcoded File Paths - Modify these paths as needed
# =========================

# Clinical correlation files
clinical_files <- list(
  KEGG = "/Users/jackie-mac/Desktop/VSCode/outputs/correlations/PDA/PDA_KEGG_clinical_correlations_filtered.csv",
  HALLMARK = "/Users/jackie-mac/Desktop/VSCode/outputs/correlations/PDA/PDA_HALLMARK_clinical_correlations_filtered.csv", 
  REACTOME = "/Users/jackie-mac/Desktop/VSCode/outputs/correlations/PDA/PDA_REACTOME_clinical_correlations_filtered.csv",
  BIOCARTA = "/Users/jackie-mac/Desktop/VSCode/outputs/correlations/PDA/PDA_BIOCARTA_clinical_correlations_filtered.csv"
)

# Genomic x radiomic correlation files  
genomic_files <- list(
  KEGG = "/Users/jackie-mac/Desktop/VSCode/outputs/correlations/PDA/PDA_KEGG_correlative_analysis.csv",
  HALLMARK = "/Users/jackie-mac/Desktop/VSCode/outputs/correlations/PDA/PDA_HALLMARK_correlative_analysis.csv",
  REACTOME = "/Users/jackie-mac/Desktop/VSCode/outputs/correlations/PDA/PDA_REACTOME_correlative_analysis.csv", 
  BIOCARTA = "/Users/jackie-mac/Desktop/VSCode/outputs/correlations/PDA/PDA_BIOCARTA_correlative_analysis.csv"
)

# Output directory
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/plots/Heatmaps/PDA/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =========================
# Function to create correlation matrix from long format
# =========================
create_correlation_matrix <- function(correlation_data, genomic_features_to_include) {
  # Filter to only include specified genomic features
  filtered_data <- correlation_data[correlation_data$GenomicFeature %in% genomic_features_to_include, ]
  
  if (nrow(filtered_data) == 0) {
    cat("Warning: No matching genomic features found\n")
    return(NULL)
  }
  
  # Create matrix using reshape2
  correlation_matrix <- acast(filtered_data, GenomicFeature ~ RadiomicFeature, 
                             value.var = "SpearmanRho", fill = NA)
  
  return(correlation_matrix)
}

# =========================
# Function to create and save heatmap
# =========================
create_heatmap <- function(correlation_matrix, pathway_name, output_dir) {
  if (is.null(correlation_matrix) || nrow(correlation_matrix) == 0) {
    cat("Warning: No data to plot for", pathway_name, "\n")
    return(NULL)
  }
  
  # Remove columns (radiomic features) that are all NA
  valid_cols <- colSums(!is.na(correlation_matrix)) > 0
  correlation_matrix <- correlation_matrix[, valid_cols, drop = FALSE]
  
  # Remove rows (genomic features) that are all NA  
  valid_rows <- rowSums(!is.na(correlation_matrix)) > 0
  correlation_matrix <- correlation_matrix[valid_rows, , drop = FALSE]
  
  if (nrow(correlation_matrix) == 0 || ncol(correlation_matrix) == 0) {
    cat("Warning: No valid correlations to plot for", pathway_name, "\n")
    return(NULL)
  }
  
  cat("Creating heatmap for", pathway_name, "with", nrow(correlation_matrix), "genomic features and", 
      ncol(correlation_matrix), "radiomic features\n")
  
  # Clean genomic feature names - remove everything before first underscore
  cleaned_row_names <- sapply(rownames(correlation_matrix), function(name) {
    if (grepl("_", name)) {
      # Find first underscore and take everything after it
      sub("^[^_]*_", "", name)
    } else {
      # If no underscore, keep original name
      name
    }
  })
  rownames(correlation_matrix) <- cleaned_row_names
  
  # Set up output filename
  output_file <- file.path(output_dir, paste0(pathway_name, "_clinical_guided_heatmap.png"))
  interactive_file <- file.path(output_dir, paste0(pathway_name, "_clinical_guided_heatmap_interactive.html"))
  
  # Calculate symmetric color scale limits based on actual data
  max_abs_corr <- max(abs(correlation_matrix), na.rm = TRUE)
  color_limits <- c(-max_abs_corr, max_abs_corr)
  
  # Create color palette
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # Calculate dynamic dimensions for better visibility
  plot_width <- max(1200, ncol(correlation_matrix) * 40)  # Increased base width and scaling
  plot_height <- max(800, nrow(correlation_matrix) * 50)  # Increased base height and scaling
  
  # =========================
  # Create Static PNG Heatmap
  # =========================
  png(output_file, width = plot_width, height = plot_height, res = 300, bg = "white")
  
  pheatmap(
    correlation_matrix,
    color = color_palette,
    breaks = seq(color_limits[1], color_limits[2], length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,    # Show genomic feature names (cleaned)
    show_colnames = TRUE,    # Show radiomic feature names
    main = paste("Clinical-Guided Radiogenomic Correlations:", pathway_name),
    fontsize_row = max(6, min(12, 400/nrow(correlation_matrix))),   # Increased row font size
    fontsize_col = max(5, min(10, 250/ncol(correlation_matrix))),   # Increased column font size
    fontsize = 12,           # Title font size
    na_col = "grey90",
    border_color = NA,       # Remove cell borders for cleaner look
    cellwidth = max(8, min(20, 800/ncol(correlation_matrix))),     # Dynamic cell width
    cellheight = max(8, min(20, 600/nrow(correlation_matrix)))     # Dynamic cell height
  )
  
  dev.off()
  cat("Static heatmap saved:", output_file, "\n")
  
  # =========================
  # Create Interactive HTML Heatmap
  # =========================
  cat("Creating interactive heatmap...\n")
  
  # Create interactive heatmap with heatmaply
  interactive_heatmap <- heatmaply(
    correlation_matrix,
    colors = colorRampPalette(c("blue", "white", "red"))(100),
    limits = color_limits,
    na_value = "grey90",
    show_dendrogram = c(TRUE, TRUE),
    row_dend_left = TRUE,
    plot_method = "plotly",
    title = paste("Interactive Clinical-Guided Radiogenomic Correlations:", pathway_name),
    xlab = "Radiomic Features",
    ylab = "Genomic Features",
    margin = list(l = 200, r = 50, t = 100, b = 200),  # Increased margins for labels
    fontsize_row = max(8, min(14, 500/nrow(correlation_matrix))),
    fontsize_col = max(6, min(12, 300/ncol(correlation_matrix))),
    subplot_widths = c(0.85, 0.15),
    subplot_heights = c(0.15, 0.85),
    row_side_colors = NULL,
    col_side_colors = NULL,
    showticklabels = c(TRUE, TRUE),
    hide_colorbar = FALSE
  )
  
  # Save interactive heatmap
  htmlwidgets::saveWidget(
    interactive_heatmap,
    file = interactive_file,
    selfcontained = TRUE,
    title = paste("Interactive Heatmap:", pathway_name)
  )
  
  cat("Interactive heatmap saved:", interactive_file, "\n")
  
  # Also save correlation matrix as CSV
  matrix_file <- file.path(output_dir, paste0(pathway_name, "_correlation_matrix.csv"))
  write.csv(correlation_matrix, matrix_file, row.names = TRUE)
  cat("Correlation matrix saved:", matrix_file, "\n")
  cat("Both static PNG and interactive HTML heatmaps created successfully!\n")
  
  return(correlation_matrix)
}

# =========================
# Main Processing Loop
# =========================
cat("Starting clinical-guided radiogenomic visualization...\n\n")

pathway_names <- names(clinical_files)
summary_stats <- data.frame(
  Pathway = character(),
  Clinical_Features = integer(),
  Genomic_Features_Found = integer(),
  Radiomic_Features = integer(),
  Total_Correlations = integer(),
  Valid_Correlations = integer(),
  stringsAsFactors = FALSE
)

for (pathway in pathway_names) {
  cat(rep("=", 60), "\n")
  cat("Processing pathway:", pathway, "\n")
  cat(rep("=", 60), "\n")
  
  # Read clinical correlation data
  cat("Reading clinical correlations from:", clinical_files[[pathway]], "\n")
  clinical_data <- fread(clinical_files[[pathway]], data.table = FALSE)
  
  # Validate clinical data columns
  if (!"Signature" %in% colnames(clinical_data)) {
    cat("Error: 'Signature' column not found in clinical data for", pathway, "\n")
    cat("Available columns:", paste(colnames(clinical_data), collapse = ", "), "\n")
    next
  }
  
  # Extract genomic features from clinical data
  clinical_genomic_features <- unique(clinical_data$Signature)
  cat("Found", length(clinical_genomic_features), "unique genomic features in clinical data\n")
  
  # Read genomic x radiomic correlation data
  cat("Reading genomic-radiomic correlations from:", genomic_files[[pathway]], "\n")
  genomic_data <- fread(genomic_files[[pathway]], data.table = FALSE)
  
  # Validate genomic data columns
  required_cols <- c("GenomicFeature", "RadiomicFeature", "SpearmanRho")
  missing_cols <- setdiff(required_cols, colnames(genomic_data))
  if (length(missing_cols) > 0) {
    cat("Error: Missing required columns in genomic data for", pathway, ":", 
        paste(missing_cols, collapse = ", "), "\n")
    cat("Available columns:", paste(colnames(genomic_data), collapse = ", "), "\n")
    next
  }
  
  # Find overlap between clinical genomic features and genomic-radiomic data
  available_genomic_features <- unique(genomic_data$GenomicFeature)
  overlapping_features <- intersect(clinical_genomic_features, available_genomic_features)
  
  cat("Genomic features in clinical data:", length(clinical_genomic_features), "\n")
  cat("Genomic features in radiogenomic data:", length(available_genomic_features), "\n") 
  cat("Overlapping genomic features:", length(overlapping_features), "\n")
  
  if (length(overlapping_features) == 0) {
    cat("Warning: No overlapping genomic features found for", pathway, "\n")
    next
  }
  
  # Create correlation matrix for overlapping features
  correlation_matrix <- create_correlation_matrix(genomic_data, overlapping_features)
  
  if (!is.null(correlation_matrix)) {
    # Create and save heatmap
    result_matrix <- create_heatmap(correlation_matrix, pathway, output_dir)
    
    # Collect summary statistics
    if (!is.null(result_matrix)) {
      total_correlations <- nrow(result_matrix) * ncol(result_matrix)
      valid_correlations <- sum(!is.na(result_matrix))
      
      summary_stats <- rbind(summary_stats, data.frame(
        Pathway = pathway,
        Clinical_Features = length(clinical_genomic_features),
        Genomic_Features_Found = nrow(result_matrix),
        Radiomic_Features = ncol(result_matrix),
        Total_Correlations = total_correlations,
        Valid_Correlations = valid_correlations,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  cat("\n")
}

# =========================
# Summary Report
# =========================
cat("\n", rep("=", 60), "\n")
cat("SUMMARY REPORT\n")
cat(rep("=", 60), "\n")

if (nrow(summary_stats) > 0) {
  print(summary_stats)
  
  # Save summary statistics
  summary_file <- file.path(output_dir, "visualization_summary.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat("\nSummary statistics saved:", summary_file, "\n")
  
  # Overall statistics
  cat("\nOverall Statistics:\n")
  cat("Total pathways processed:", nrow(summary_stats), "\n")
  cat("Total clinical features across all pathways:", sum(summary_stats$Clinical_Features), "\n")
  cat("Total genomic features visualized:", sum(summary_stats$Genomic_Features_Found), "\n")
  cat("Total radiomic features included:", length(unique(unlist(sapply(pathway_names, function(p) {
    if (file.exists(genomic_files[[p]])) {
      data <- fread(genomic_files[[p]], data.table = FALSE)
      return(unique(data$RadiomicFeature))
    }
    return(character(0))
  })))), "\n")
  
} else {
  cat("No successful visualizations were created. Please check input files and data formats.\n")
}

cat("\nVisualization script completed successfully!\n")
cat("Output directory:", output_dir, "\n")
