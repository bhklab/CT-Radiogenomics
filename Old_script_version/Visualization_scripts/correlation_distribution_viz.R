# ===============================================================================
# Correlation Distribution Visualization for Clinically Relevant Features
# ===============================================================================
# 
# Purpose: Creates histograms showing the distribution of correlation values 
#          between genomic and radiomic features for genomic features that 
#          show significant clinical correlations.
#
# Description:
#   This script takes clinical correlation results and genomic-radiomic correlation
#   results from four pathway collections (KEGG, HALLMARK, REACTOME, BIOCARTA).
#   It filters genomic-radiomic correlations to include only genomic features
#   present in clinical correlation data, then visualizes the distribution of
#   correlation values as histograms.
#
# Input Requirements:
#   1. Clinical correlation files (4): CSV files containing clinical associations
#      - Required columns: "Signature" (genomic feature names)
#   2. Genomic-radiomic correlation files (4): CSV files containing radiogenomic correlations
#      - Required columns: "GenomicFeature", "RadiomicFeature", "SpearmanRho"
#   3. Files should correspond by pathway collection (same order)
#
# Output:
#   - Histogram plots saved as PNG files for each pathway collection
#   - Combined histogram plot showing all pathways together
#   - Summary statistics of correlation distributions
#
# Usage:
#   Simply run: Rscript correlation_distribution_viz.R
#   (Modify the hardcoded file paths in the script as needed)
#
# Dependencies: ggplot2, data.table, dplyr
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
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
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/plots/Histograms/correlative_analysis/PDA"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =========================
# Function to filter correlations by clinical relevance
# =========================
filter_correlations_by_clinical <- function(genomic_data, clinical_genomic_features) {
  # Filter genomic-radiomic correlations to include only clinically relevant features
  filtered_data <- genomic_data[genomic_data$GenomicFeature %in% clinical_genomic_features, ]
  
  if (nrow(filtered_data) == 0) {
    cat("Warning: No matching genomic features found\n")
    return(NULL)
  }
  
  # Remove rows with missing correlation values
  filtered_data <- filtered_data[!is.na(filtered_data$SpearmanRho), ]
  
  return(filtered_data)
}

# =========================
# Function to create correlation distribution histogram
# =========================
create_correlation_histogram <- function(correlation_data, pathway_name, output_dir) {
  if (is.null(correlation_data) || nrow(correlation_data) == 0) {
    cat("Warning: No data to plot for", pathway_name, "\n")
    return(NULL)
  }
  
  cat("Creating histogram for", pathway_name, "with", nrow(correlation_data), "correlations\n")
  
  # Calculate summary statistics
  corr_values <- correlation_data$SpearmanRho
  mean_corr <- mean(corr_values, na.rm = TRUE)
  median_corr <- median(corr_values, na.rm = TRUE)
  sd_corr <- sd(corr_values, na.rm = TRUE)
  min_corr <- min(corr_values, na.rm = TRUE)
  max_corr <- max(corr_values, na.rm = TRUE)
  
  # Create histogram plot
  p <- ggplot(correlation_data, aes(x = SpearmanRho)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black", size = 0.3) +
    geom_vline(xintercept = mean_corr, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = median_corr, color = "orange", linetype = "dashed", size = 1) +
    labs(
      title = paste("Correlation Distribution:", pathway_name),
      subtitle = paste("Clinical-Guided Radiogenomic Correlations (n =", nrow(correlation_data), ")"),
      x = "Spearman Correlation Coefficient",
      y = "Frequency",
      caption = paste("Mean =", round(mean_corr, 3), "| Median =", round(median_corr, 3), "| SD =", round(sd_corr, 3), 
                     "| Min =", round(min_corr, 3), "| Max =", round(max_corr, 3))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      plot.caption = element_text(size = 10, hjust = 0.5)
    ) +
    scale_x_continuous(breaks = seq(-1, 1, 0.2), limits = c(-1, 1))
  
  # Save histogram plot
  output_file <- file.path(output_dir, paste0(pathway_name, "_correlation_distribution.png"))
  ggsave(output_file, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
  cat("Histogram saved:", output_file, "\n")
  
  # Return summary statistics
  return(data.frame(
    Pathway = pathway_name,
    N_Correlations = nrow(correlation_data),
    Mean = mean_corr,
    Median = median_corr,
    SD = sd_corr,
    Min = min_corr,
    Max = max_corr,
    Q25 = quantile(corr_values, 0.25, na.rm = TRUE),
    Q75 = quantile(corr_values, 0.75, na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
}

# =========================
# Main Processing Loop
# =========================
cat("Starting correlation distribution visualization...\n\n")

pathway_names <- names(clinical_files)
all_correlations <- data.frame()
summary_stats <- data.frame()

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
  
  # Filter correlations by clinical relevance
  filtered_correlations <- filter_correlations_by_clinical(genomic_data, clinical_genomic_features)
  
  if (!is.null(filtered_correlations)) {
    # Add pathway information
    filtered_correlations$Pathway <- pathway
    
    # Append to combined dataset
    all_correlations <- rbind(all_correlations, filtered_correlations)
    
    # Create individual histogram
    pathway_stats <- create_correlation_histogram(filtered_correlations, pathway, output_dir)
    
    if (!is.null(pathway_stats)) {
      summary_stats <- rbind(summary_stats, pathway_stats)
    }
  }
  
  cat("\n")
}

# =========================
# Create Combined Histogram
# =========================
if (nrow(all_correlations) > 0) {
  cat("Creating combined histogram for all pathways...\n")
  
  # Create combined histogram with facets
  combined_plot <- ggplot(all_correlations, aes(x = SpearmanRho, fill = Pathway)) +
    geom_histogram(bins = 40, alpha = 0.7, color = "black", size = 0.2) +
    facet_wrap(~Pathway, scales = "free_y", ncol = 2) +
    labs(
      title = "Correlation Distributions Across Pathway Collections",
      subtitle = "Clinical-Guided Radiogenomic Correlations",
      x = "Spearman Correlation Coefficient",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 1))
  
  # Save combined histogram
  combined_file <- file.path(output_dir, "combined_correlation_distributions.png")
  ggsave(combined_file, plot = combined_plot, width = 12, height = 8, dpi = 300, bg = "white")
  cat("Combined histogram saved:", combined_file, "\n")
  
  # Create overall distribution plot
  overall_plot <- ggplot(all_correlations, aes(x = SpearmanRho)) +
    geom_histogram(bins = 60, fill = "steelblue", alpha = 0.7, color = "black", size = 0.3) +
    geom_vline(xintercept = mean(all_correlations$SpearmanRho, na.rm = TRUE), 
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = median(all_correlations$SpearmanRho, na.rm = TRUE), 
               color = "orange", linetype = "dashed", size = 1) +
    labs(
      title = "Overall Correlation Distribution",
      subtitle = paste("All Clinical-Guided Radiogenomic Correlations (n =", nrow(all_correlations), ")"),
      x = "Spearman Correlation Coefficient",
      y = "Frequency",
      caption = paste("Mean =", round(mean(all_correlations$SpearmanRho, na.rm = TRUE), 3), 
                     "| Median =", round(median(all_correlations$SpearmanRho, na.rm = TRUE), 3),
                     "| Min =", round(min(all_correlations$SpearmanRho, na.rm = TRUE), 3),
                     "| Max =", round(max(all_correlations$SpearmanRho, na.rm = TRUE), 3))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      plot.caption = element_text(size = 11, hjust = 0.5)
    ) +
    scale_x_continuous(breaks = seq(-1, 1, 0.2), limits = c(-1, 1))
  
  # Save overall histogram
  overall_file <- file.path(output_dir, "overall_correlation_distribution.png")
  ggsave(overall_file, plot = overall_plot, width = 10, height = 6, dpi = 300, bg = "white")
  cat("Overall histogram saved:", overall_file, "\n")
}

# =========================
# Summary Report
# =========================
cat("\n", rep("=", 60), "\n")
cat("CORRELATION DISTRIBUTION SUMMARY\n")
cat(rep("=", 60), "\n")

if (nrow(summary_stats) > 0) {
  print(summary_stats)
  
  # Save summary statistics
  summary_file <- file.path(output_dir, "correlation_distribution_summary.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat("\nSummary statistics saved:", summary_file, "\n")
  
  # Overall statistics
  cat("\nOverall Statistics:\n")
  cat("Total pathways processed:", nrow(summary_stats), "\n")
  cat("Total correlations analyzed:", sum(summary_stats$N_Correlations), "\n")
  cat("Mean correlation across all pathways:", round(mean(all_correlations$SpearmanRho, na.rm = TRUE), 3), "\n")
  cat("Median correlation across all pathways:", round(median(all_correlations$SpearmanRho, na.rm = TRUE), 3), "\n")
  cat("Standard deviation across all pathways:", round(sd(all_correlations$SpearmanRho, na.rm = TRUE), 3), "\n")
  
  # Correlation strength categories
  strong_pos <- sum(all_correlations$SpearmanRho > 0.5, na.rm = TRUE)
  moderate_pos <- sum(all_correlations$SpearmanRho > 0.3 & all_correlations$SpearmanRho <= 0.5, na.rm = TRUE)
  weak_pos <- sum(all_correlations$SpearmanRho > 0 & all_correlations$SpearmanRho <= 0.3, na.rm = TRUE)
  weak_neg <- sum(all_correlations$SpearmanRho < 0 & all_correlations$SpearmanRho >= -0.3, na.rm = TRUE)
  moderate_neg <- sum(all_correlations$SpearmanRho < -0.3 & all_correlations$SpearmanRho >= -0.5, na.rm = TRUE)
  strong_neg <- sum(all_correlations$SpearmanRho < -0.5, na.rm = TRUE)
  
  cat("\nCorrelation Strength Distribution:\n")
  cat("Strong positive (> 0.5):", strong_pos, "(", round(100*strong_pos/nrow(all_correlations), 1), "%)\n")
  cat("Moderate positive (0.3 to 0.5):", moderate_pos, "(", round(100*moderate_pos/nrow(all_correlations), 1), "%)\n")
  cat("Weak positive (0 to 0.3):", weak_pos, "(", round(100*weak_pos/nrow(all_correlations), 1), "%)\n")
  cat("Weak negative (-0.3 to 0):", weak_neg, "(", round(100*weak_neg/nrow(all_correlations), 1), "%)\n")
  cat("Moderate negative (-0.5 to -0.3):", moderate_neg, "(", round(100*moderate_neg/nrow(all_correlations), 1), "%)\n")
  cat("Strong negative (< -0.5):", strong_neg, "(", round(100*strong_neg/nrow(all_correlations), 1), "%)\n")
  
} else {
  cat("No correlation distributions were created. Please check input files and data formats.\n")
}

cat("\nCorrelation distribution visualization completed successfully!\n")
cat("Output directory:", output_dir, "\n")
