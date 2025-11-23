library(data.table)
library(ggplot2)

# ---- USER INPUTS ----
input_files <- list(
    HNSCC = "/Users/jackie-mac/Desktop/VSCode/procdata/RNA_enrichments/HNSCC_compiled_enrichment.csv",
    BRCA = "/Users/jackie-mac/Desktop/VSCode/procdata/RNA_enrichments/BRCA_compiled_enrichment.csv",
    KIRC = "/Users/jackie-mac/Desktop/VSCode/procdata/RNA_enrichments/KIRC_compiled_enrichment.csv",
    LGG = "/Users/jackie-mac/Desktop/VSCode/procdata/RNA_enrichments/LGG_compiled_enrichment.csv",
    GBM = "/Users/jackie-mac/Desktop/VSCode/procdata/RNA_enrichments/GBM_compiled_enrichment.csv",
    CCRCC = "/Users/jackie-mac/Desktop/VSCode/procdata/RNA_enrichments/CCRCC_compiled_enrichment.csv",
    PDA = "/Users/jackie-mac/Desktop/VSCode/procdata/RNA_enrichments/PDA_compiled_enrichment.csv"
)
radiomics_file <- "/Users/jackie-mac/Desktop/VSCode/outputs/Correlative_analysis/HNSCC_radiomics_features_filtered.csv"
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/plots/histograms/signature_distribution"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Pathway source patterns
sources <- c("KEGG_MEDICUS", "HALLMARK", "REACTOME", "BIOCARTA")

for (cancer in names(input_files)) {
  file <- input_files[[cancer]]
  cat("Processing:", cancer, "\n")
  df <- fread(file, data.table = FALSE, check.names = FALSE)
  
  # Assume first column is pathway name, rest are enrichment scores
  pathway_col <- 1
  score_cols <- 2:ncol(df)
  
  for (source in sources) {
    # Grep for pathway rows belonging to this source
    source_rows <- grepl(source, df[[pathway_col]], ignore.case = TRUE)
    if (!any(source_rows)) {
      cat("  No pathways found for", source, "in", cancer, "\n")
      next
    }
    # Get all enrichment scores for this source (flatten to vector)
    scores <- as.numeric(unlist(df[source_rows, score_cols, drop = FALSE]))
    scores <- scores[!is.na(scores)]
    if (length(scores) == 0) {
      cat("  No enrichment scores for", source, "in", cancer, "\n")
      next
    }
    # Plot histogram
    plot_df <- data.frame(EnrichmentScore = scores)
    p <- ggplot(plot_df, aes(x = EnrichmentScore)) +
      geom_histogram(binwidth = 0.05, fill = "#0072B2", color = "black") +
      labs(
        title = paste(cancer, "-", source, "Enrichment Score Distribution"),
        x = "Enrichment Score",
        y = "Count"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    out_file <- file.path(output_dir, paste0(cancer, "_", source, "_enrichment_histogram.png"))
    ggsave(out_file, plot = p, width = 7, height = 5, dpi = 300, bg = "white")
    cat("  Histogram saved for", source, "at", out_file, "\n")
  }
}

# ---- Radiomics Features Histogram ----
cat("Processing radiomics features...\n")
radiomics_df <- fread(radiomics_file, data.table = FALSE, check.names = FALSE)

# Feature prefixes to include
feature_prefixes <- c("original_", "wavelet-", "square_", "squareroot_", "logarithm_", "exponential_", "gradient_")
feature_pattern <- paste0("^(", paste(feature_prefixes, collapse = "|"), ")")
feature_cols <- grep(feature_pattern, colnames(radiomics_df), value = TRUE)

if (length(feature_cols) == 0) {
  cat("No radiomics features found with specified prefixes.\n")
} else {
  # Gather all values for these features
  feature_vals <- as.numeric(unlist(radiomics_df[, feature_cols, drop = FALSE]))
  feature_vals <- feature_vals[!is.na(feature_vals)]
  if (length(feature_vals) == 0) {
    cat("No valid numeric values found for radiomics features.\n")
  } else {
    plot_df_rad <- data.frame(RadiomicsValue = feature_vals)
    p_rad <- ggplot(plot_df_rad, aes(x = RadiomicsValue)) +
      geom_histogram(binwidth = 1000, fill = "#009E73", color = "black") +
      labs(
        title = "Radiomics Feature Value Distribution",
        subtitle = paste("Features:", length(feature_cols), "| Values:", length(feature_vals)),
        x = "Radiomics Feature Value",
        y = "Count"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    out_file_rad <- file.path(output_dir, "radiomics_features_histogram.png")
    ggsave(out_file_rad, plot = p_rad, width = 7, height = 5, dpi = 300, bg = "white")
    cat("  Histogram saved for radiomics features at", out_file_rad, "\n")
  }
}