suppressPackageStartupMessages(library(data.table))

# ---- USER INPUTS ----
group1 <- c("BRCA", "GBM", "KIRC", "LGG")
group2 <- c("CCRCC", "HNSCC", "PDA")
cancer_types <- c(group1, group2)

# Hardcoded file paths
    outcome_file_group1 <- "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/combined_TCGA_clinical_outcomes.csv"
outcome_file_group2 <- "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/CPTAC_combined_OS_clinical_outcomes.csv"
treatment_type_files <- list(
  BRCA = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/BRCA_updated_clinical_data2.csv",
  GBM  = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/GBM_updated_clinical_data2.csv",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/KIRC_updated_clinical_data2.csv",
  LGG  = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/LGG_updated_clinical_data2.csv",
  CCRCC = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/CCRCC_updated_clinical_data2.csv",
  HNSCC = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/HNSCC_updated_clinical_data2.csv",
  PDA   = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/PDA_updated_clinical_data2.csv"
)
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/intersect_outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

is_na_or_blank <- function(x) is.na(x) | x == "" | x == "NA" | x == "--"

for (cancer in cancer_types) {
  cat("\nProcessing", cancer, "...\n")
  # Select correct outcome file
  if (cancer %in% group1) {
    outcome_mat <- fread(outcome_file_group1, data.table = FALSE, check.names = FALSE)
    id_col <- grep(paste0("^", cancer, "_SampleID$"), colnames(outcome_mat), value = TRUE)[1]
    outcome_col <- grep(paste0("^", cancer, "_overall_survival$"), colnames(outcome_mat), value = TRUE)[1]
    if (is.na(id_col) || is.na(outcome_col)) {
      warning(sprintf("Could not find all required columns for %s in outcome matrix.", cancer))
      next
    }
    outcome_sub <- outcome_mat[, c(id_col, outcome_col), drop = FALSE]
    colnames(outcome_sub) <- c("sample_id", "outcome")
  } else if (cancer %in% group2) {
    outcome_mat <- fread(outcome_file_group2, data.table = FALSE, check.names = FALSE)
    id_col <- grep(paste0("^", cancer, "_SampleID$"), colnames(outcome_mat), value = TRUE)[1]
    os_col <- grep(paste0("^", cancer, "_OS_days$"), colnames(outcome_mat), value = TRUE)[1]
    pfs_col <- grep(paste0("^", cancer, "_PFS_days$"), colnames(outcome_mat), value = TRUE)[1]
    if (is.na(id_col) || is.na(os_col) || is.na(pfs_col)) {
      warning(sprintf("Could not find all required columns for %s in outcome matrix.", cancer))
      next
    }
    outcome_sub <- outcome_mat[, c(id_col, os_col, pfs_col), drop = FALSE]
    colnames(outcome_sub) <- c("sample_id", "OS_days", "PFS_days")
  } else {
    warning(sprintf("Cancer type %s not recognized in any group.", cancer))
    next
  }

  # Read treatment type matrix for this cancer
  treat_file <- treatment_type_files[[cancer]]
  treat_mat <- fread(treat_file, data.table = FALSE, check.names = FALSE)
  treat_id_col <- grep("^cases.submitter_id", colnames(treat_mat), value = TRUE)[1]
  if (is.na(treat_id_col)) {
    warning(sprintf("Could not find required sample_id column for %s in treatment type matrix.", cancer))
    next
  }
  treat_sub <- treat_mat
  setnames(treat_sub, old = treat_id_col, new = "sample_id")

  # Merge on sample_id only (allowing for duplicates in treatment matrix)
  merged <- merge(
    treat_sub, outcome_sub,
    by = "sample_id",
    all.x = TRUE, suffixes = c("", ".outcome")
  )

  # Identify valid samples
  if (cancer %in% group1) {
    valid <- !is_na_or_blank(merged$outcome)
    filtered <- merged[valid, , drop = FALSE]
    # Move outcome column to the last column
    outcome_data <- filtered$outcome
    filtered$outcome <- NULL
    filtered$outcome <- outcome_data
  } else if (cancer %in% group2) {
    valid <- !(is_na_or_blank(merged$OS_days) | is_na_or_blank(merged$PFS_days))
    filtered <- merged[valid, , drop = FALSE]
    # Move OS_days and PFS_days columns to the last columns
    os_data <- filtered$OS_days
    pfs_data <- filtered$PFS_days
    filtered$OS_days <- NULL
    filtered$PFS_days <- NULL
    filtered$OS_days <- os_data
    filtered$PFS_days <- pfs_data
  }

  output_file <- file.path(output_dir, paste0(cancer, "_clinical_intersect.csv"))
  fwrite(filtered, output_file, row.names = FALSE, quote = TRUE)
  cat(sprintf("Filtered intersect for %s saved to %s\n", cancer, output_file))

  # Print missing samples
  if (cancer %in% group1) {
    missing <- merged[!valid, "sample_id", drop = FALSE]
    if (nrow(missing) > 0) {
      cat("Samples missing outcome data for", cancer, ":\n")
      for (i in seq_len(nrow(missing))) {
        cat(" -", missing$sample_id[i], "\n")
      }
    } else {
      cat("All samples for", cancer, "have outcome data.\n")
    }
  } else if (cancer %in% group2) {
    missing <- merged[!valid, "sample_id", drop = FALSE]
    if (nrow(missing) > 0) {
      cat("Samples missing OS_days or PFS_days data for", cancer, ":\n")
      for (i in seq_len(nrow(missing))) {
        cat(" -", missing$sample_id[i], "\n")
      }
    } else {
      cat("All samples for", cancer, "have OS_days and PFS_days data.\n")
    }
  }
}