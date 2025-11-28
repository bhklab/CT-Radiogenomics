# ===============================================================================
# Pie Chart Script
# This script generates pie charts for each treatment type in a dataset
# ===============================================================================

library(ggplot2)
library(data.table)

# Set output directory
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/piecharts/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Read your data
df <- read.csv("/Users/jackie-mac/Desktop/VSCode/data/clinical/Clinical Data Variables - Treatment type summary.csv", check.names = FALSE)

# Get treatment type columns (skip the first column, which is CancerType)
treatment_types <- colnames(df)[-1]

# Loop through each treatment type and generate pie chart for all treatment types
for (treatment_col in treatment_types) {
  pie_data <- data.frame(
    CancerType = df[[1]],
    Count = df[[treatment_col]]
  )
  total_samples <- sum(pie_data$Count, na.rm = TRUE)
  pie_data <- pie_data[pie_data$Count > 0, ]
  pie_data$Percent <- round(100 * pie_data$Count / total_samples, 1)
  pie_data$Label <- paste0(pie_data$CancerType, "\n", pie_data$Count, " (", pie_data$Percent, "%)")
  p <- ggplot(pie_data, aes(x = "", y = Count, fill = CancerType)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3) +
    labs(title = paste("Cancer Type Distribution for", treatment_col, "\n(Total samples:", total_samples, ")")) +
    theme_void()
  ggsave(file.path(output_dir, paste0("piechart_", treatment_col, ".png")), plot = p, width = 6, height = 6)
  print(paste("Pie chart saved for", treatment_col, "at", output_dir))
}
# ===============================================================================
# Gy filtering for clinical data
# ===============================================================================
# ===============================================================================
# Palliative/Tumour Control Sample Splitter Script
# This script loops through multiple CSV files and splits samples into
# "palliative control" and "tumour control" based on treatment dose.
# ===============================================================================
library(data.table)

# ---- USER INPUTS ----
input_files <- list(
  GBM  = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/GBM_updated_clinical_data.csv",
  BRCA = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/BRCA_updated_clinical_data.csv",
  LGG  = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/LGG_updated_clinical_data.csv",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/KIRC_updated_clinical_data.csv"
)
datasets <- names(input_files)

treatment_types_to_keep <- c(
  "Radiation, External Beam",
  "Radiation, Stereotactic/Gamma Knife/SRS",
  "Radiation, Intensity-Modulated Radiotherapy"
)

for (dataset in datasets) {
  cat("\nProcessing dataset:", dataset, "\n")
  file <- input_files[[dataset]]
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  
  dose_col <- paste0("treatments.treatment_dose_", dataset, "_patientids")
  id_col   <- paste0("cases.submitter_id_", dataset, "_patientids")
  type_col <- paste0("treatments.treatment_type_", dataset, "_patientids")
  
  if (!(dose_col %in% colnames(df)) || !(id_col %in% colnames(df)) || !(type_col %in% colnames(df))) {
    cat("  Columns not found in", file, "- skipping.\n")
    next
  }
  
  # Filter for specified treatment types only
  keep <- df[[type_col]] %in% treatment_types_to_keep
  df <- df[keep, , drop = FALSE]
  
  dose_raw <- df[[dose_col]]
  ids <- df[[id_col]]
  
  # Identify NA/--/blank cells
  is_no_data <- is.na(dose_raw) | dose_raw == "" | dose_raw == "'--"
  # Convert to numeric for valid entries, NA for others
  dose_vals <- suppressWarnings(as.numeric(ifelse(is_no_data, NA, dose_raw)))
  
  # Split into lists
  palliative_control <- ids[!is.na(dose_vals) & dose_vals > 20]
  tumour_control    <- ids[!is.na(dose_vals) & dose_vals <= 20]
  no_data           <- ids[is_no_data]
  
  cat("  Palliative control samples (dose > 20):\n")
  if (length(palliative_control) > 0) {
    for (id in palliative_control) cat("   -", id, "\n")
  } else {
    cat("    None\n")
  }
  
  cat("  Tumour control samples (dose <= 20):\n")
  if (length(tumour_control) > 0) {
    for (id in tumour_control) cat("   -", id, "\n")
  } else {
    cat("    None\n")
  }
  
  cat("  No data samples (dose is NA, blank, or '--'):\n")
  if (length(no_data) > 0) {
    for (id in no_data) cat("   -", id, "\n")
  } else {
    cat("    None\n")
  }
}
# ===============================================================================
# Histogram Generation for Radiation Gray values
# ===============================================================================
library(ggplot2)
library(data.table)

# ---- USER INPUTS ----
input_files <- list(
  BRCA = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/BRCA_clinical_intersect.csv",
  GBM  = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/GBM_clinical_intersect.csv",
  LGG  = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/LGG_clinical_intersect.csv",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/KIRC_clinical_intersect.csv"
)
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs/histograms/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Treatment type groups
group_ext <- c(
  "Radiation, External Beam",
  "Radiation, Stereotactic/Gamma Knife/SRS",
  "Radiation, Intensity-Modulated Radiotherapy"
)
group_nos <- c("Radiation Therapy, NOS")

# Collect all dose values across all datasets
all_dose_ext <- c()
all_dose_nos <- c()
total_ext <- 0
total_nos <- 0
n_no_data_ext <- 0
n_no_data_nos <- 0

for (dataset in names(input_files)) {
  file <- input_files[[dataset]]
  cat("\nProcessing:", dataset, "\n")
  df <- fread(file, data.table = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Find the correct columns for this dataset
  dose_col <- grep("treatments.treatment_dose", colnames(df), value = TRUE)[1]
  type_col <- grep("treatments.treatment_type", colnames(df), value = TRUE)[1]
  if (is.na(dose_col) || is.na(type_col)) {
    cat("  Columns not found for", dataset, "- skipping.\n")
    next
  }
  
  # External Beam + SRS + IMRT group
  df_ext <- df[df[[type_col]] %in% group_ext, ]
  total_ext <- total_ext + nrow(df_ext)
  is_no_data_ext <- is.na(df_ext[[dose_col]]) | df_ext[[dose_col]] == "" | df_ext[[dose_col]] == "'--"
  n_no_data_ext <- n_no_data_ext + sum(is_no_data_ext)
  dose_vals_ext <- suppressWarnings(as.numeric(df_ext[[dose_col]]))
  all_dose_ext <- c(all_dose_ext, dose_vals_ext[!is_no_data_ext])
  
  # NOS group
  df_nos <- df[df[[type_col]] %in% group_nos, ]
  total_nos <- total_nos + nrow(df_nos)
  is_no_data_nos <- is.na(df_nos[[dose_col]]) | df_nos[[dose_col]] == "" | df_nos[[dose_col]] == "'--"
  n_no_data_nos <- n_no_data_nos + sum(is_no_data_nos)
  dose_vals_nos <- suppressWarnings(as.numeric(df_nos[[dose_col]]))
  all_dose_nos <- c(all_dose_nos, dose_vals_nos[!is_no_data_nos])
}

# Plot External Beam + SRS + IMRT (all cancers combined)
if (length(all_dose_ext) > 0) {
  plot_df_ext <- data.frame(Dose = all_dose_ext)
  p_ext <- ggplot(plot_df_ext, aes(x = Dose)) +
    geom_histogram(binwidth = 1, fill = "#0072B2", color = "black", boundary = 0) +
    labs(
      title = "Gray Value Distribution (All Cancers)",
      subtitle = paste("External Beam/SRS/IMRT\nTotal:", total_ext, 
                       "| With data:", length(all_dose_ext), 
                       "| No data:", n_no_data_ext),
      x = "Gray Value",
      y = "Count"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 20))
  out_file_ext <- file.path(output_dir, "histogram_gray_ALL_ExternalBeam_SRS_IMRT.png")
  ggsave(out_file_ext, plot = p_ext, width = 10, height = 6, dpi = 300, bg = "white")
  cat("Histogram saved for External Beam/SRS/IMRT at", out_file_ext, "\n")
}

# Plot NOS (all cancers combined)
if (length(all_dose_nos) > 0) {
  plot_df_nos <- data.frame(Dose = all_dose_nos)
  p_nos <- ggplot(plot_df_nos, aes(x = Dose)) +
    geom_histogram(binwidth = 1, fill = "#D55E00", color = "black", boundary = 0) +
    labs(
      title = "Gray Value Distribution (All Cancers)",
      subtitle = paste("Radiation, NOS\nTotal:", total_nos, 
                       "| With data:", length(all_dose_nos), 
                       "| No data:", n_no_data_nos),
      x = "Gray Value",
      y = "Count"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 20))
  out_file_nos <- file.path(output_dir, "histogram_gray_ALL_NOS.png")
  ggsave(out_file_nos, plot = p_nos, width = 10, height = 6, dpi = 300, bg = "white")
  cat("Histogram saved for Radiation, NOS at", out_file_nos, "\n")
}
# =========================
# Batch Heatmap Generator for Correlation Matrices
# =========================

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
heatmap_width <- 12  # in inches
heatmap_height <- 12 # in inches
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

# =========================
# Histogram Generation for Correlation Distributions
# =========================

library(ggplot2)
library(data.table)

# ---- USER INPUTS ----
# Hardcoded input files for FPKM and TPM
fpkm_files <- list(
  KEGG = "/Users/jackie-mac/Desktop/VSCode/outputs/Correlative_analysis/HNSCC/HNSCC_KEGG_correlative_analysis.csv",
  HALLMARK = "/Users/jackie-mac/Desktop/VSCode/outputs/Correlative_analysis/HNSCC/HNSCC_HALLMARK_correlative_analysis.csv",
  REACTOME = "/Users/jackie-mac/Desktop/VSCode/outputs/Correlative_analysis/HNSCC/HNSCC_REACTOME_correlative_analysis.csv",
  BIOCARTA = "/Users/jackie-mac/Desktop/VSCode/outputs/Correlative_analysis/HNSCC/HNSCC_BIOCARTA_correlative_analysis_FPKM.csv"
)
tpm_files <- list(
  KEGG = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_KEGG_correlative_analysis.csv",
  HALLMARK = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_HALLMARK_correlative_analysis.csv",
  REACTOME = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_REACTOME_correlative_analysis.csv",
  BIOCARTA = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_BIOCARTA_correlative_analysis.csv"
)
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/Visualizations/histograms"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- HISTOGRAM FUNCTION ----
plot_corr_hist <- function(df, matrix_name, norm_type, color, output_dir) {
  # Bin width of 0.1, no forced x-axis limits
  p <- ggplot(df, aes(x = SpearmanRho)) +
    geom_histogram(binwidth = 0.1, fill = color, color = "black", boundary = -1, closed = "left") +
    labs(
      title = paste0(matrix_name, " Correlation Distribution (", norm_type, ")"),
      x = "Spearman Rho",
      y = "Count"
    ) +
    theme_minimal(base_size = 18) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  out_file <- file.path(output_dir, paste0("histogram_", matrix_name, "_", norm_type, ".png"))
  ggsave(out_file, plot = p, width = 14, height = 8, dpi = 300, bg = "white")
  cat("Histogram saved for", matrix_name, norm_type, "at", out_file, "\n")
}

# ---- FPKM HISTOGRAMS ----
for (matrix_name in names(fpkm_files)) {
  file <- fpkm_files[[matrix_name]]
  df <- fread(file, data.table = FALSE, check.names = FALSE)
  if (!"SpearmanRho" %in% colnames(df)) stop(paste("SpearmanRho column not found in", file))
  plot_corr_hist(df, matrix_name, "FPKM", "orange", output_dir)
}

# ---- TPM HISTOGRAMS ----
for (matrix_name in names(tpm_files)) {
  file <- tpm_files[[matrix_name]]
  df <- fread(file, data.table = FALSE, check.names = FALSE)
  if (!"SpearmanRho" %in% colnames(df)) stop(paste("SpearmanRho column not found in", file))
  plot_corr_hist(df, matrix_name, "TPM", "#377EB8", output_dir)
}

# =========================
# Density Plot Generator (R)
# =========================

# Usage:
# plot_density("input_file.csv", "column_name", "output_plot.png")
# Or run as a script:
# Rscript miscellaneous_scripts.R input_file.csv column_name output_plot.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

plot_density <- function(input_file, column_name, output_file, plot_title = NULL) {
  # Read data
  df <- fread(input_file, data.table = FALSE, check.names = FALSE)
  if (!(column_name %in% colnames(df))) {
    stop(paste("Column", column_name, "not found in", input_file))
  }
  values <- df[[column_name]]
  # Remove NA and non-numeric
  values <- suppressWarnings(as.numeric(values))
  values <- values[!is.na(values)]
  if (length(values) == 0) stop("No valid numeric data found in the specified column.")

  # Set output file if not provided
  if (is.null(output_file)) {
    output_file <- sub("\\.[^.]+$", paste0("_density_", column_name, ".png"), input_file)
  }
  # Set plot title if not provided
  if (is.null(plot_title)) {
    plot_title <- paste("Density Plot of", column_name, names(input_file))
  }

  # Plot
  p <- ggplot(data.frame(Value = values), aes(x = Value)) +
    geom_density(fill = "#377EB8", alpha = 0.7, color = "black") +
    labs(title = plot_title, x = column_name, y = "Density") +
    theme_minimal(base_size = 16) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  ggsave(output_file, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
  cat("Density plot saved to", output_file, "\n")
}

input_files <- c(KEGG_bin ="/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_KEGG_cox_results_binary.csv",
                 HALLMARK_bin = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_HALLMARK_cox_results_binary.csv",
                 BIOCARTA_bin = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_BIOCARTA_cox_results_binary.csv",
                 REACTOME_bin = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_REACTOME_cox_results_binary.csv",
                 KEGG_con = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_KEGG_cox_results_continuous.csv",
                 HALLMARK_con = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_HALLMARK_cox_results_continuous.csv",
                 BIOCARTA_con = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_BIOCARTA_cox_results_continuous.csv",
                 REACTOME_con = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_REACTOME_cox_results_continuous.csv",
                 radiomics_bin = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_radiomics_cox_results_binary.csv",
                 radiomics_cont = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_radiomics_cox_results_continuous.csv")  # Fill with your input file paths, e.g. c("/path/to/file1.csv", "/path/to/file2.csv")
column_name <- "HR"
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/Visualizations/density_plots"  # Set your output directory
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
for (input_name in names(input_files)) {
  input <- input_files[[input_name]]
  if (!file.exists(input)) {
    cat("Input file does not exist:", input, "\n")
    next
  }
  # Generate output file name based on input file and column
  base <- tools::file_path_sans_ext(basename(input))
  output_file <- file.path(output_dir, paste0(base, "_density_", column_name, ".png"))
  plot_title <- paste("Density Plot of", column_name, "for", input_name)
  plot_density(input, column_name, output_file, plot_title)
}