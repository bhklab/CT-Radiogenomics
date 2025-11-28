# ===============================================================================
# Cox Proportional Hazards Model Results Visualization Suite
# ===============================================================================
# 
# Purpose: Creates comprehensive visualizations for Cox regression analysis results
#          including hazard ratio distributions, top feature rankings, radiogenomic
#          correlations, and treatment outcome comparisons.
#
# Description:
#   This script generates multiple visualization types to explore and present
#   Cox proportional hazards model results. 
#   It includes:
#   - Histograms of hazard ratios for binary and continuous models
#   - Top significant features ranked by effect size and FDR
#   - Treatment outcome distribution plots
# Input Requirements:
#   1. Cox regression results: CSV files with hazard ratios, p-values, confidence intervals
#   2. Radiogenomic correlation data: Optional correlation matrices
#   3. Clinical outcome data: Treatment types and survival information
#
# Output Visualizations:
#   - Hazard ratio distribution histograms
#   - Top significant features ranked by effect size
#   - Volcano plots showing HR vs statistical significance
#   - Correlation heatmaps between genomic and radiomic features
#   - Treatment outcome distribution plots
#   - Combined multi-panel figures for comprehensive analysis
#
# Visualization Types:
#   - Statistical distributions and summaries
#   - Feature ranking and importance plots
#   - Correlation network visualizations
#   - Clinical outcome stratifications
#   - Publication-ready multi-panel figures
#
# Usage:
#   1. Configure input file paths and visualization parameters
#   2. Run: Rscript cox_visualization.R
#   3. Review generated plots and select for publication
#
# Dependencies: ggplot2, data.table, corrplot, cowplot
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

library(ggplot2)
library(data.table)
# library(pheatmap) # Uncomment if you want to use pheatmap for heatmaps
library(tools) # for file_path_sans_ext

# ---- USER INPUTS ----
# Assign filepaths to designated variables for each data type
binary_files <- list(
  kegg = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_KEGG_cox_results_binary.csv",
  hallmark = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_HALLMARK_cox_results_binary.csv",
  biocarta = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_BIOCARTA_cox_results_binary.csv",
  reactome = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_REACTOME_cox_results_binary.csv",
  radiomics = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_radiomics_cox_results_binary.csv"
)
continuous_files <- list(
  kegg = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_KEGG_cox_results_continuous.csv",
  hallmark = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_HALLMARK_cox_results_continuous.csv",
  biocarta = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_BIOCARTA_cox_results_continuous.csv",
  reactome = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_REACTOME_cox_results_continuous.csv",
  radiomics = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_radiomics_cox_results_continuous.csv"
)
merged_data_files <- list(
  kegg = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_extracted_clinical.csv",
  hallmark = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_extracted_clinical.csv",
  biocarta = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_extracted_clinical.csv",
  reactome = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_extracted_clinical.csv",
  radiomics = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_extracted_clinical.csv"
)              # Path to merged clinical+genomic data (for treatment outcome plot)
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/Visualizations/survival_outcomes"                        # Directory to save plots
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
top_n <- 15

# New input: file containing highly correlated pairs (genomic, radiomic)
highly_correlated_pairs_file <- c(
  kegg = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_KEGG_correlative_analysis.csv",
  hallmark = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_HALLMARK_correlative_analysis.csv",
  biocarta = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_BIOCARTA_correlative_analysis.csv",
  reactome = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_REACTOME_correlative_analysis.csv"
)
# ---- LOAD DATA ----
# results_binary <- fread(results_binary_file)
# results_continuous <- fread(results_continuous_file)

# ---- 1. HR DISTRIBUTION PLOTS FOR MULTIPLE FILES ----
for (name in names(binary_files)) {
  f <- binary_files[[name]]
  dat <- fread(f)
  label <- paste0(name, "_binary")
  p <- ggplot(dat, aes(x = HR)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    labs(title = paste0("Distribution of Hazard Ratios (Binary): ", name), x = "Hazard Ratio", y = "Count") +
    theme_minimal(base_size = 14) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_line(color = "grey95"))
  ggsave(file.path(output_dir, paste0(label, "_HR_distribution.png")), p, width = 6, height = 4, dpi = 300)
  dat_ordered <- dat[order(FDR, abs(HR - 1))]
  print(head(dat_ordered, top_n))
}

for (name in names(continuous_files)) {
  f <- continuous_files[[name]]
  dat <- fread(f)
  label <- paste0(name, "_continuous")
  p <- ggplot(dat, aes(x = HR)) +
    geom_histogram(bins = 30, fill = "orange", color = "black") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    labs(title = paste0("Distribution of Hazard Ratios (Continuous): ", name), x = "Hazard Ratio", y = "Count") +
    theme_minimal(base_size = 14) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_line(color = "grey95"))
  ggsave(file.path(output_dir, paste0(label, "_HR_distribution.png")), p, width = 6, height = 4, dpi = 300)
  dat_ordered <- dat[order(FDR, abs(HR - 1))]
  print(head(dat_ordered, top_n))
}

# ---- 2. TOP FEATURES BY FDR AND HR (ALL FILES) ----
top_n <- 15
all_top_genomic <- character(0)
all_top_radiomic <- character(0)

# Helper to extract top features from a file
top_features_from_file <- function(file, fdr_thresh = 0.05) {
  dat <- fread(file)
  dat_ordered <- dat[order(dat[['FDR']], abs(dat[['HR']] - 1))]
  # Genomic features: those not containing 'radiomic' in the name
  genomic_features <- dat_ordered[!grepl("radiomic", dat_ordered[['Signature']], ignore.case=TRUE) & dat_ordered[['FDR']] < fdr_thresh, ]
  top_genomic <- genomic_features[['Signature']]
  # Radiomic features: those containing 'radiomic' in the name
  radiomic_features <- dat_ordered[grepl("radiomic", dat_ordered[['Signature']], ignore.case=TRUE) & dat_ordered[['FDR']] < fdr_thresh, ]
  top_radiomic <- radiomic_features[['Signature']]
  list(top_genomic = top_genomic, top_radiomic = top_radiomic)
}

# Loop through all binary and continuous files
for (f in c(unlist(binary_files), unlist(continuous_files))) {
  if (file.exists(f)) {
    res <- top_features_from_file(f, fdr_thresh = 0.05)
    all_top_genomic <- unique(c(all_top_genomic, res$top_genomic))
    all_top_radiomic <- unique(c(all_top_radiomic, res$top_radiomic))
    # Optionally, save top features per file
    fwrite(data.table(GenomicFeature = res$top_genomic), file.path(output_dir, paste0(basename(f), "_top_genomic.csv")))
    fwrite(data.table(RadiomicFeature = res$top_radiomic), file.path(output_dir, paste0(basename(f), "_top_radiomic.csv")))
  }
}

# ---- 3. CHECK HIGHLY CORRELATED PAIRS BETWEEN ALL TOP GENOMIC AND RADIOMIC FEATURES ----
for (corr_name in names(highly_correlated_pairs_file)) {
  corr_file <- highly_correlated_pairs_file[[corr_name]]
  if (file.exists(corr_file)) {
    pairs_df <- fread(corr_file) # columns: GenomicFeature, RadiomicFeature
    pairs_to_check <- expand.grid(GenomicFeature = all_top_genomic, RadiomicFeature = all_top_radiomic, stringsAsFactors = FALSE)
    merged_pairs <- merge(pairs_to_check, pairs_df, by = c("GenomicFeature", "RadiomicFeature"))
    if (nrow(merged_pairs) > 0) {
      cat(paste0("Highly correlated top feature pairs found (", corr_name, "):\n"))
      print(merged_pairs)
      fwrite(merged_pairs, file.path(output_dir, paste0(corr_name, "_highly_correlated_top_feature_pairs.csv")))
    } else {
      cat(paste0("No highly correlated pairs among top genomic and radiomic features (", corr_name, ").\n"))
    }
  }
}

# ---- 4. OS_days DISTRIBUTION FOR MULTIPLE FILES ----
for (name in names(merged_data_files)) {
  f <- merged_data_files[[name]]
  if (file.exists(f)) {
    merged_data <- fread(f)
    # OS_days histogram only
    if ("OS_days_HNSCC" %in% colnames(merged_data)) {
      p_os <- ggplot(merged_data, aes(x = OS_days_HNSCC)) +
        geom_histogram(bins = 30, fill = "purple", color = "black") +
        labs(title = paste0("Distribution of Overall Survival (OS_days): ", name), x = "OS_days_HNSCC", y = "Count") +
        theme_minimal(base_size = 14) +
        theme(panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA),
              panel.grid.major = element_line(color = "grey90"),
              panel.grid.minor = element_line(color = "grey95"))
      ggsave(file.path(output_dir, paste0(name, "_OS_days_distribution.png")), p_os, width = 6, height = 4, dpi = 300)
    }
  }
}

top_n <- 50
# Extract top radiomic features from the single radiomics file
radiomics_file <- binary_files$radiomics
if (file.exists(radiomics_file)) {
  radiomic_res <- top_features_from_file(radiomics_file, fdr_thresh = 0.05)
  top_radiomic <- radiomic_res$top_radiomic
} else {
  top_radiomic <- character(0)
}

# For each genomic file, compare its significant genomic features to the significant radiomic features
# Define correlation matrix files for each pathway group
correlation_files <- list(
  kegg = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/kegg_correlation_matrix.csv",
  hallmark = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/hallmark_correlation_matrix.csv",
  biocarta = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/biocarta_correlation_matrix.csv",
  reactome = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/reactome_correlation_matrix.csv"
)

for (genomic_name in c("kegg", "hallmark", "biocarta", "reactome")) {
  genomic_file <- binary_files[[genomic_name]]
  correlation_file <- correlation_files[[genomic_name]]
  if (file.exists(genomic_file) && file.exists(correlation_file)) {
    genomic_res <- top_features_from_file(genomic_file, fdr_thresh = 0.05)
    top_genomic <- genomic_res$top_genomic
    # Save significant features per file
    fwrite(data.table(GenomicFeature = top_genomic), file.path(output_dir, paste0(genomic_name, "_significant_genomic.csv")))
    fwrite(data.table(RadiomicFeature = top_radiomic), file.path(output_dir, "radiomics_significant_radiomic.csv"))
    # Check highly correlated pairs using the pathway-specific correlation matrix
    pairs_df <- fread(correlation_file) # Must have 'GenomicFeature', 'RadiomicFeature', 'SpearmanRho', 'p.value', 'adj.p.value'
    pairs_to_check <- expand.grid(GenomicFeature = top_genomic, RadiomicFeature = top_radiomic, stringsAsFactors = FALSE)
    merged_pairs <- merge(pairs_to_check, pairs_df, by = c("GenomicFeature", "RadiomicFeature"))
    # Filter for abs(SpearmanRho) > 0.7
    filtered_pairs <- merged_pairs[abs(SpearmanRho) > 0.7, ]
    # Add FDR for genomic and radiomic features
    # Read cox results for radiomics
    radiomics_file <- binary_files$radiomics
    radiomics_cox <- fread(radiomics_file)
    genomic_cox <- fread(genomic_file)
    filtered_pairs <- merge(filtered_pairs, genomic_cox[, .(GenomicFeature = Signature, Genomic_FDR = FDR)], by = "GenomicFeature", all.x = TRUE)
    filtered_pairs <- merge(filtered_pairs, radiomics_cox[, .(RadiomicFeature = Signature, Radiomic_FDR = FDR)], by = "RadiomicFeature", all.x = TRUE)
    # Reorder columns for clarity
    out_cols <- c("GenomicFeature", "RadiomicFeature", "Genomic_FDR", "Radiomic_FDR", "SpearmanRho", "p.value", "adj.p.value")
    filtered_pairs <- filtered_pairs[, ..out_cols]
    out_file <- file.path(output_dir, paste0(genomic_name, "_vs_radiomics_highly_correlated_pairs.csv"))
    if (nrow(filtered_pairs) > 0) {
      cat(paste0("Highly correlated significant feature pairs found for ", genomic_name, " vs radiomics (|SpearmanRho| > 0.7):\n"))
      print(filtered_pairs)
      fwrite(filtered_pairs, out_file)
    } else {
      cat(paste0("No highly correlated pairs (|SpearmanRho| > 0.7) among significant genomic and radiomic features for ", genomic_name, " vs radiomics.\n"))
      fwrite(filtered_pairs, out_file) # Write empty file with headers
    }
  }
}
