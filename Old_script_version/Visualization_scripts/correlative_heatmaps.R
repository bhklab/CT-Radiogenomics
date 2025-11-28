# ===============================================================================
# Radiogenomic Correlation Heatmap Visualization Generator
# ===============================================================================
# 
# Purpose: Generates correlation heatmaps visualizing Spearman correlations between
#          genomic pathway signatures and radiomic features for multiple pathway
#          databases to identify radiogenomic association patterns.
#
# Description:
#   This script processes correlation analysis results from multiple pathway sources
#   (KEGG, HALLMARK, REACTOME, BIOCARTA) and creates publication-ready heatmaps
#   showing the strength and direction of correlations between genomic signatures
#   and radiomic features. Each heatmap provides a comprehensive view of
#   radiogenomic associations for a specific pathway database.
#
# Input Requirements:
#   1. Correlation analysis results: CSV files containing pairwise correlations
#      - Columns: GenomicFeature, RadiomicFeature, SpearmanRho, PValue, AdjPValue
#   2. Multiple pathway source files (KEGG, HALLMARK, REACTOME, BIOCARTA)
#   3. Results from correlative_analysis.R script output
#
# Output:
#   - High-resolution heatmap images (PNG format, 300 DPI)
#   - One heatmap per pathway database showing correlation matrix
#   - Genomic features as rows, radiomic features as columns
#   - Color-coded correlation coefficients (-1 to +1 scale)
#   - Hierarchical clustering for pattern identification
#
# Visualization Features:
#   - Red-Blue color scheme: Red (negative), White (no correlation), Blue (positive)
#   - Hierarchical clustering of rows and columns to group similar patterns
#   - Customizable dimensions and resolution for publication quality
#   - Feature names hidden for clarity (can be enabled if needed)
#
# Usage:
#   1. Configure input file paths for each pathway database
#   2. Set output directory and visualization parameters
#   3. Run: Rscript correlative_heatmaps.R
#
# Dependencies: data.table, pheatmap, RColorBrewer, reshape2
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(RColorBrewer)
  library(reshape2)
})

# ---- USER INPUTS ----
input_correlation_files <- list(
  KEGG = "/Users/jackie-mac/Desktop/VSCode/outputs/Correlative_analysis/HNSCC/HNSCC_KEGG_correlative_analysis.csv",
  HALLMARK = "/Users/jackie-mac/Desktop/VSCode/outputs/Correlative_analysis/HNSCC/HNSCC_HALLMARK_correlative_analysis.csv",
  REACTOME = "/Users/jackie-mac/Desktop/VSCode/outputs/Correlative_analysis/HNSCC/HNSCC_REACTOME_correlative_analysis.csv",
  BIOCARTA = "/Users/jackie-mac/Desktop/VSCode/outputs/Correlative_analysis/HNSCC/HNSCC_BIOCARTA_correlative_analysis.csv"
)
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/Visualizations/Heatmaps/HNSCC"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
heatmap_width <- 14  # in inches
heatmap_height <- 14 # in inches
cluster_rows_cols <- TRUE

for (matrix_name in names(input_correlation_files)) {
  input_file <- input_correlation_files[[matrix_name]]
  output_file <- file.path(output_dir, paste0(matrix_name, "_correlation_heatmap.png"))
  heatmap_title <- paste(matrix_name, "Genomic vs Radiomic Correlation (Spearman rho)")

  cat("Reading correlation summary from:", input_file, "\n")
  df <- fread(input_file, data.table = FALSE, check.names = FALSE)
  # Ensure column names are as expected
  colnames(df)[1:5] <- c("GenomicFeature", "RadiomicFeature", "SpearmanRho", "PValue", "AdjPValue")

  # Reshape to matrix: rows = Genomic, columns = Radiomic, values = SpearmanRho
  mat <- acast(df, GenomicFeature ~ RadiomicFeature, value.var = "SpearmanRho")

  cat("Generating heatmap for", matrix_name, "...\n")
  png(output_file, width = heatmap_width, height = heatmap_height, units = "in", res = 300)
  pheatmap(
    mat,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
    cluster_rows = cluster_rows_cols,
    cluster_cols = cluster_rows_cols,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = heatmap_title
  )
  dev.off()
  cat("Heatmap saved to:", output_file, "\n")
}
