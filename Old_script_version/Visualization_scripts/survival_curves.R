# ===============================================================================
# Kaplan-Meier Survival Curve Generator for Radiogenomic Signatures
# ===============================================================================
# 
# Purpose: Generates Kaplan-Meier survival curves comparing high vs low expression
#          groups for genomic signatures and radiomic features using median split
#          stratification to visualize survival differences.
#
# Description:
#   This script creates survival curve plots for radiogenomic features by dividing
#   patients into high and low expression groups based on median values. It generates
#   publication-ready Kaplan-Meier curves with log-rank test statistics to assess
#   the prognostic value of different signatures and features.
#
# Input Requirements:
#   1. Feature expression data: CSV with samples as rows, features as columns
#   2. Clinical survival data: Must include overall survival time and event status
#   3. Matching sample IDs between feature and clinical datasets
#
# Output:
#   - Kaplan-Meier survival curve plots (PDF/PNG format)
#   - Log-rank test p-values for survival differences
#   - Risk tables showing number at risk over time
#   - Summary statistics for high vs low groups
#
# Analysis Method:
#   - Median split stratification for high/low group assignment
#   - Kaplan-Meier survival estimation
#   - Log-rank test for statistical significance testing
#   - Customizable plotting parameters and aesthetics
#
# Usage:
#   1. Configure input file paths and parameters in the script
#   2. Run: Rscript survival_curves.R
#   3. Review generated survival plots and statistics
#
# Dependencies: survival, survminer, ggplot2
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

library(survival)
library(survminer)
library(data.table)
library(ggplot2)

# ---- USER INPUTS ----
# File lists (update paths as needed)
cox_files <- list(
  kegg = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/PANCAN/COMBINED_KEGG_pancancer_cox_results_binary.csv",
  hallmark = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/PANCAN/COMBINED_HALLMARK_pancancer_cox_results_binary.csv",
  biocarta = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/PANCAN/COMBINED_BIOCARTA_pancancer_cox_results_binary.csv",
  reactome = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/PANCAN/COMBINED_REACTOME_pancancer_cox_results_binary.csv"
)
# No longer using correlation files
# correlation_files <- list(...)
clinical_file <- "/Users/jackie-mac/Desktop/VSCode/procdata/PANCAN/COMBINED_harmonized_clinical.csv"
genomic_files <- list(
  kegg = "/Users/jackie-mac/Desktop/VSCode/procdata/PANCAN/COMBINED_harmonized_kegg_genomics.csv",
  hallmark = "/Users/jackie-mac/Desktop/VSCode/procdata/PANCAN/COMBINED_harmonized_hallmark_genomics.csv",
  biocarta = "/Users/jackie-mac/Desktop/VSCode/procdata/PANCAN/COMBINED_harmonized_biocarta_genomics.csv",
  reactome = "/Users/jackie-mac/Desktop/VSCode/procdata/PANCAN/COMBINED_harmonized_reactome_genomics.csv"
)
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/plots/survival_curves/PANCAN"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- 1. STANDARDIZE SAMPLE ID COLUMN NAMES ----
standardize_sampleid <- function(file, old_col) {
  dt <- fread(file)
  if (colnames(dt)[1] != "SampleID") {
    setnames(dt, old_col, "SampleID")
  }
  return(dt)
}
clinical <- standardize_sampleid(clinical_file, "cases.submitter_id")

for (g in names(genomic_files)) {
  dt <- fread(genomic_files[[g]])
  if (colnames(dt)[1] != "SampleID") setnames(dt, 1, "SampleID")
  assign(paste0(g, "_genomic"), dt)
}

# ---- 2. FILTER COX MODEL RESULTS BY FDR ----
filtered_features <- list()
cox_results <- list()  # Store full cox results for C-index extraction
for (g in names(cox_files)) {
  cox <- fread(cox_files[[g]])
  
  # Strictly filter by FDR < 0.05
  cox_filtered <- cox[!is.na(FDR) & FDR < 0.05]
  
  # Sort by FDR for easier interpretation (smallest to largest)
  cox_filtered <- cox_filtered[order(FDR, decreasing = FALSE)]
  
  # Log information about initial filtering
  cat(g, "pathway: Initially filtered from", nrow(cox), "to", nrow(cox_filtered), 
      "signatures with FDR < 0.05\n")
  
  # Take only top 20 features with lowest FDR values (if available)
  top_count <- min(20, nrow(cox_filtered))
  cox_top <- cox_filtered[1:top_count]
  
  # Store both filtered signatures and full filtered results
  filtered_features[[g]] <- cox_top$Signature
  cox_results[[g]] <- cox_top  # Store top results for later use
  
  # Log information about top feature selection
  cat(g, "pathway: Selected top", top_count, "signatures with lowest FDR values\n")
  
  if (top_count > 0) {
    cat("  Lowest FDR =", format(cox_top$FDR[1], scientific = TRUE, digits = 2),
        ", Highest FDR =", format(cox_top$FDR[top_count], scientific = TRUE, digits = 2), "\n")
  }
}
# ---- 3. USE GENOMIC FEATURES DIRECTLY (SKIP CORRELATION) ----
final_features <- list()
for (g in names(cox_files)) {
  # Skip if we have no features passing FDR filter
  if (length(filtered_features[[g]]) == 0) {
    cat(g, "pathway: No features passed FDR < 0.05 filter, skipping further analysis\n")
    final_features[[g]] <- data.table()
    next
  }
  
  # Create a data.table with the genomic features directly from Cox results
  features_dt <- data.table(
    GenomicFeature = filtered_features[[g]],
    RadiomicFeature = "None"  # Placeholder since we're not using correlations
  )
  
  final_features[[g]] <- features_dt
  cat(g, "pathway: Using", nrow(features_dt), "genomic features with FDR < 0.05\n")
}

# ---- 4. GENERATE SURVIVAL CURVES FOR EACH GENOMIC FEATURE ----
for (g in names(final_features)) {
  if (nrow(final_features[[g]]) == 0) next
  genomic_mat <- get(paste0(g, "_genomic"))
  cox_data <- cox_results[[g]]  # Get cox results for this pathway
  
  for (i in seq_len(nrow(final_features[[g]]))) {
    feature <- final_features[[g]]$GenomicFeature[i]
    # Not using radiomic features anymore
    
    # Extract C-index from Cox results (already calculated for binary high/low groups)
    cox_c_index <- cox_data[Signature == feature, C_index]
    if (length(cox_c_index) == 0 || is.na(cox_c_index)) {
      cox_c_index <- "N/A"
    } else {
      cox_c_index <- round(as.numeric(cox_c_index), 3)
    }
    
    # Merge clinical and feature
    merged <- merge(clinical, genomic_mat[, .(SampleID, value = get(feature))], by = "SampleID")
    
    # Remove rows with missing values - adjust column names to match your PANCAN dataset
    merged <- merged[!is.na(value) & !is.na(OS_days) & !is.na(OS_event)]
    
    if (nrow(merged) < 10) {
      cat("Warning: Too few samples with complete data for", feature, ". Skipping.\n")
      next
    }
    
    # Split by median (same as original Cox model)
    median_val <- median(merged$value, na.rm = TRUE)
    merged$group <- ifelse(merged$value > median_val, "High", "Low")
    
    # Get FDR value from Cox results
    cox_fdr <- cox_data[Signature == feature, FDR]
    if (length(cox_fdr) == 0 || is.na(cox_fdr)) {
      cat("Warning: FDR value not found for", feature, ". Skipping.\n")
      next
    }
    
    # Format FDR value
    fdr_formatted <- format(cox_fdr, scientific = TRUE, digits = 2)
    
    # Create annotation text with C-index and FDR from Cox model
    c_index_text <- paste0("C-index = ", cox_c_index, "\nFDR = ", fdr_formatted)
    
    # Survival curve
    surv_obj <- Surv(time = merged$OS_days, event = merged$OS_event)
    plot_title <- paste0("Survival Curve: High vs Low ", feature, " (", g, ")\nFDR = ", fdr_formatted)
    p <- ggsurvplot(
      survfit(surv_obj ~ group, data = merged),
      data = merged, pval = TRUE, risk.table = TRUE, conf.int = TRUE,
      title = plot_title,
      legend.title = feature, legend.labs = c("Low", "High"),
      palette = c("#377EB8", "#E41A1C"),
      ggtheme = theme_minimal(base_size = 18) +
        theme(panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA),
              legend.background = element_rect(fill = "white", color = NA),
              legend.key = element_rect(fill = "white", color = NA),
              plot.title = element_text(size = 9),  # Much smaller title font
              plot.subtitle = element_text(size = 9)) # Smaller subtitle font if present
    )
    
    # Add C-index annotation to the plot
    p$plot <- p$plot + 
      annotate("text", 
               x = Inf, 
               y = Inf, 
               label = c_index_text, 
               hjust = 1.05, 
               vjust = 2.5, 
               size = 5, 
               color = "black",
               fontface = "bold")
    
    # Clean feature name for filename
    clean_feature <- gsub("[^A-Za-z0-9_]", "_", feature)
    
    # Format FDR for filename (replace decimal point with p)
    fdr_str <- sprintf("FDR%.5f", cox_fdr)
    fdr_str <- gsub("\\.", "p", fdr_str)
    
    ggsave(file.path(output_dir, paste0("survival_curve_", g, "_", clean_feature, "_", fdr_str, ".png")), 
           p$plot, width = 14, height = 7, dpi = 300)
  }
}
cat("Survival curves generated for all selected feature pairs.\n")
