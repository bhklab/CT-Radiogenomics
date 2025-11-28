# ===============================================================================
# Forest Plot Generator for Cox Regression Results Visualization
# ===============================================================================
# 
# Purpose: Creates publication-ready forest plots to visualize Cox proportional
#          hazards regression results, showing hazard ratios with confidence
#          intervals for genomic signatures and radiomic features.
#
# Description:
#   This script processes Cox regression results and generates forest plots
#   displaying hazard ratios, confidence intervals, and p-values for the most
#   statistically significant features. It automatically selects top features
#   based on FDR-corrected p-values and creates professional visualizations.
#
# Input Requirements:
#   1. Cox regression results: CSV files containing hazard ratios, confidence intervals, p-values
#   2. Multiple input files can be processed in batch
#   3. Required columns: Feature names, HR, CI_lower, CI_upper, p-values
#
# Output:
#   - Forest plot visualizations (PDF/PNG format)
#   - Hazard ratios with 95% confidence intervals
#   - Statistical significance indicators
#   - Feature rankings by effect size and significance
#
# Visualization Features:
#   - Selects top 20 FDR-significant features per dataset
#   - Displays hazard ratios with confidence intervals
#   - Color-codes by statistical significance
#   - Professional formatting for publication
#
# Usage:
#   1. Configure input file paths and selection criteria
#   2. Run: Rscript forest_plot_generator.R
#   3. Review generated forest plots
#
# Dependencies: ggplot2, forestplot, dplyr
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================
# and generates a forest plot for each. The y-axis is the signature, x-axis is HR, box size is -log10(FDR),
# and error bars are standard error (SE). Input columns: 'Signature', 'HR', 'FDR', 'CI_lower', 'CI_upper'.
#
# Usage:
#   - Set input_files as a named vector or list of CSVs (e.g. KEGG, HALLMARK, etc)
#   - Set output_dir
#   - Run the script

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(dplyr)
})

# ---- USER INPUTS ----
input_files <- list(
  KEGG_bin = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/HNSCC/HNSCC_KEGG_cox_results_binary.csv",
  HALLMARK_bin = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/HNSCC/HNSCC_HALLMARK_cox_results_binary.csv",
  BIOCARTA_bin = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/HNSCC/HNSCC_BIOCARTA_cox_results_binary.csv",
  REACTOME_bin = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/HNSCC/HNSCC_REACTOME_cox_results_binary.csv",
  radiomics = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/HNSCC/HNSCC_radiomics_cox_results_binary.csv"
)
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/plots/forest_plots/HNSCC"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

top_n <- 20  # Number of top FDR-significant signatures to plot

plot_forest <- function(df, group_name, output_file) {
  # Check required columns
  required_cols <- c("Signature", "HR", "FDR", "CI_lower", "CI_upper")
  if (!all(required_cols %in% colnames(df))) {
    stop(paste("Missing columns in", group_name, ":", paste(setdiff(required_cols, colnames(df)), collapse=", ")))
  }
  # Calculate SE from confidence intervals (SE = (log(CI_upper) - log(CI_lower)) / (2*1.96))
  df$SE <- (log(df$CI_upper) - log(df$CI_lower)) / (2 * 1.96)
  # For error bars, use the actual CI limits, not HR ± SE
  # Filter for finite values
  df <- df %>% filter(is.finite(HR), is.finite(FDR), is.finite(SE))
  # Order by FDR, select top N
  df <- df %>% arrange(FDR) %>% head(top_n)
  # For plotting, order y-axis by HR or FDR
  # Wrap long signature names to two lines if needed (at ~40 characters)
  wrap_signature <- function(x, width = 40) {
    sapply(x, function(s) {
      if (nchar(s) > width) {
        # Try to break at the last space before width, else just break
        idx <- max(gregexpr(' ', substr(s, 1, width))[[1]])
        if (idx > 1 && idx < nchar(s)) {
          paste0(substr(s, 1, idx-1), "\n", substr(s, idx+1, nchar(s)))
        } else {
          paste(strwrap(s, width = width), collapse = "\n")
        }
      } else {
        s
      }
    }, USE.NAMES = FALSE)
  }
  df$Signature_wrapped <- wrap_signature(df$Signature)
  df$Signature_wrapped <- factor(df$Signature_wrapped, levels = rev(df$Signature_wrapped))
  # Calculate box size
  df$Size <- -log10(df$FDR + 1e-10)  # Add small value to avoid log(0)
  # Forest plot
  # Use patchwork to add a left-aligned title above the plot area
  suppressPackageStartupMessages(library(patchwork))
  p_main <- ggplot(df, aes(x = HR, y = Signature_wrapped)) +
    geom_point(aes(size = Size), shape = 22, fill = "#377EB8", color = "black", alpha = 0.8) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.3, color = "gray30") +
    scale_size_continuous(name = "-log10(FDR)", range = c(3, 10)) +
    labs(
      x = "Hazard Ratio (HR)",
      y = "Genomic Signature"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 11),
      axis.title.x = element_text(size = 11),
      plot.title = element_blank()
    )
  p_title <- ggplot() +
    theme_void() +
    labs(title = paste0("Forest Plot: Top ", top_n, " FDR-significant (", group_name, ")")) +
    theme(
      plot.title = element_text(hjust = 0, face = "plain", size = 12, margin = margin(b = 5), vjust = 0)
    )
  p_combined <- p_title / p_main + plot_layout(heights = c(0.06, 1))
  ggsave(output_file, plot = p_combined, width = 15, height = 8, dpi = 300, bg = "white")
  cat("Forest plot saved to", output_file, "\n")
}

for (group_name in names(input_files)) {
  input_file <- input_files[[group_name]]
  if (!file.exists(input_file)) {
    cat("Input file does not exist:", input_file, "\n")
    next
  }
  df <- fread(input_file, data.table = FALSE, check.names = FALSE)
  # Try to guess signature column if not present
  if (!"Signature" %in% colnames(df)) {
    sig_col <- colnames(df)[1]
    colnames(df)[colnames(df) == sig_col] <- "Signature"
  }
  output_file <- file.path(output_dir, paste0(group_name, "_forest_plot.png"))
  plot_forest(df, group_name, output_file)
}
