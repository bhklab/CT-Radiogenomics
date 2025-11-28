# ===============================
# gen_clin_newID.R
# -------------------------------
# Purpose:
#   Harmonizes clinical data and compiled genomic enrichment scores by separating pathway groups,
#   transposing data format, and aligning sample IDs across all datasets.
#
# Inputs:
#   1. clinical_file: Path to the clinical data file (CSV).
#   2. compiled_genomics_file: Path to compiled genomics matrix (CSV) - pathways × samples format
#   3. output_dir: Directory to write harmonized outputs.
#
# Outputs:
#   - Harmonized clinical file
#   - Four separated genomics pathway files (Hallmark, KEGG, Reactome, BioCarta) in samples × pathways format
#   - All files have matching sample IDs for downstream analysis.
#
# Main Steps:
#   - Reads clinical and compiled genomics files
#   - Separates genomics data by pathway group (HALLMARK_, KEGG_, REACTOME_, BIOCARTA_)
#   - Transposes genomics data from pathways × samples to samples × pathways format
#   - Finds common sample IDs between clinical and genomic datasets
#   - Filters all datasets to matching sample IDs
#   - Outputs harmonized files for each data type
# ===============================

suppressPackageStartupMessages(library(data.table))

# ---- SCRIPT DESCRIPTION ----
# This script harmonizes clinical data and compiled genomic enrichment scores by separating pathway groups
# and ensuring consistent sample IDs across all datasets. It performs the following steps:
#
# 1. Reads clinical data and compiled genomics file (pathways × samples format).
# 2. Separates genomics data into four pathway groups based on pathway name prefixes:
#    - HALLMARK_ pathways
#    - KEGG_ pathways  
#    - REACTOME_ pathways
#    - BIOCARTA_ pathways
# 3. Transposes each pathway group from pathways × samples to samples × pathways format.
# 4. Identifies common sample IDs between clinical data and genomic datasets.
# 5. Filters all datasets to only include samples present in both clinical and genomic data.
# 6. Writes harmonized clinical and four genomics pathway datasets to separate CSV files.
#
# The output files ensure that clinical and genomic datasets have matching sample IDs for downstream analysis.

# ---- USER INPUTS ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript gen_clin_newID.R <clinical_file> <compiled_genomics_file> <output_dir>")
}
clinical_file <- args[1]
compiled_genomics_file <- args[2]
output_dir <- args[3]
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- READ DATA ----
clinical <- fread(clinical_file, data.table = FALSE, quote = "\"", check.names = FALSE)
compiled_genomics <- fread(compiled_genomics_file, data.table = FALSE, quote = "\"", check.names = FALSE)

cat("Initial data dimensions:\n")
cat("Clinical:", nrow(clinical), "samples\n")
cat("Compiled genomics:", nrow(compiled_genomics), "pathways x", ncol(compiled_genomics) - 1, "samples\n")

# ---- SEPARATE GENOMICS DATA BY PATHWAY GROUPS ----
# Assume first column contains pathway names, remaining columns are samples
pathway_names <- compiled_genomics[[1]]
sample_columns <- compiled_genomics[, -1, drop = FALSE]

cat("\nSeparating pathways by group...\n")

# Define pathway group patterns
pathway_groups <- list(
  HALLMARK = "^HALLMARK_",
  KEGG = "^KEGG_", 
  REACTOME = "^REACTOME_",
  BIOCARTA = "^BIOCARTA_"
)

# Separate pathways by group
separated_pathways <- list()
for (group_name in names(pathway_groups)) {
  pattern <- pathway_groups[[group_name]]
  group_indices <- grep(pattern, pathway_names, ignore.case = TRUE)
  
  if (length(group_indices) > 0) {
    # Extract pathways for this group (pathways × samples format)
    group_pathways <- pathway_names[group_indices]
    group_data <- sample_columns[group_indices, , drop = FALSE]
    rownames(group_data) <- group_pathways
    
    # Transpose to samples × pathways format
    group_data_t <- as.data.frame(t(group_data))
    
    # Add SampleID column (sample names become first column)
    group_data_t <- data.frame(SampleID = rownames(group_data_t), group_data_t, stringsAsFactors = FALSE)
    rownames(group_data_t) <- NULL
    
    separated_pathways[[group_name]] <- group_data_t
    cat(paste0(group_name, ": ", length(group_pathways), " pathways, ", nrow(group_data_t), " samples\n"))
  } else {
    cat(paste0("Warning: No pathways found for ", group_name, "\n"))
    # Create empty dataframe with SampleID column for consistency
    separated_pathways[[group_name]] <- data.frame(SampleID = character(0), stringsAsFactors = FALSE)
  }
}

# Extract the separated datasets for easier reference
hallmark <- separated_pathways[["HALLMARK"]]
kegg <- separated_pathways[["KEGG"]]
reactome <- separated_pathways[["REACTOME"]]
biocarta <- separated_pathways[["BIOCARTA"]]

# ---- IDENTIFY SAMPLE ID COLUMNS ----
# For clinical data, look for cases.submitter_id column, fallback to first column
clinical_id_col <- if("cases.submitter_id" %in% colnames(clinical)) {
  "cases.submitter_id"
} else {
  colnames(clinical)[1]
}

# For genomic data, use SampleID column (created during transposition)
genomic_id_col <- "SampleID"

cat("\nUsing sample ID columns:\n")
cat("Clinical:", clinical_id_col, "\n")
cat("Genomic datasets:", genomic_id_col, "\n")

# ---- HANDLE DUPLICATE SAMPLE IDs IN CLINICAL DATA ----
clinical_ids <- clinical[[clinical_id_col]]

# Check for duplicate sample IDs in clinical data
duplicate_ids <- clinical_ids[duplicated(clinical_ids) | duplicated(clinical_ids, fromLast = TRUE)]
has_duplicates <- length(duplicate_ids) > 0

if (has_duplicates) {
  cat("\nFound", length(duplicate_ids), "duplicate sample IDs in clinical data\n")
  cat("Examples of duplicated IDs:", paste(head(unique(duplicate_ids), 3), collapse=", "), "\n")
  
  # Create a mapping of original IDs to new unique IDs
  id_mapping <- data.frame(
    original_id = clinical_ids,
    stringsAsFactors = FALSE
  )
  
  # Add a counter for each ID to create unique identifiers
  id_counts <- table(clinical_ids)
  new_ids <- character(length(clinical_ids))
  
  for (i in seq_along(clinical_ids)) {
    id <- clinical_ids[i]
    # If this ID appears multiple times, add a suffix
    if (id_counts[id] > 1) {
      # Find how many times we've seen this ID so far
      prev_count <- sum(clinical_ids[1:i] == id)
      new_ids[i] <- paste0(id, "_dup", prev_count)
    } else {
      new_ids[i] <- id
    }
  }
  
  # Create the mapping
  id_mapping$new_id <- new_ids
  
  # Update the clinical data with the new IDs
  clinical[[clinical_id_col]] <- new_ids
  
  cat("\nCreated unique identifiers for duplicate clinical samples\n")
  cat("Example mappings:\n")
  for (dup_id in head(unique(duplicate_ids), 3)) {
    mappings <- id_mapping[id_mapping$original_id == dup_id, ]
    cat("  Original:", dup_id, "-> New:", paste(mappings$new_id, collapse=", "), "\n")
  }
  
  # Original clinical IDs are now the modified ones
  clinical_ids <- clinical[[clinical_id_col]]
} else {
  cat("\nNo duplicate sample IDs found in clinical data\n")
}

# Get sample IDs from each genomic dataset (only if they have data)
genomic_ids_list <- list()
for (group_name in c("HALLMARK", "KEGG", "REACTOME", "BIOCARTA")) {
  group_data <- separated_pathways[[group_name]]
  if (nrow(group_data) > 0 && genomic_id_col %in% colnames(group_data)) {
    genomic_ids_list[[group_name]] <- group_data[[genomic_id_col]]
  }
}

if (length(genomic_ids_list) == 0) {
  stop("No genomic datasets contain sample data. Check pathway separation logic.")
}

# Find sample IDs that exist in ALL genomic datasets
common_genomic_ids <- Reduce(intersect, genomic_ids_list)
cat("Sample IDs common to all genomic datasets:", length(common_genomic_ids), "\n")

# If we created unique identifiers for duplicate samples, we need to map them to genomic IDs
if (exists("id_mapping") && has_duplicates) {
  # Get unique original IDs that have matches in genomic data
  original_genomic_matches <- intersect(unique(id_mapping$original_id), common_genomic_ids)
  cat("Original sample IDs matching between clinical and genomic data:", length(original_genomic_matches), "\n")
  
  # Create expanded list of genomic IDs to include duplicates
  expanded_genomic_ids <- character(0)
  
  for (group_name in names(separated_pathways)) {
    group_data <- separated_pathways[[group_name]]
    if (nrow(group_data) > 0 && genomic_id_col %in% colnames(group_data)) {
      # Create a new dataframe to hold the expanded genomic data
      expanded_group_data <- data.frame()
      
      for (id in original_genomic_matches) {
        # Get all new IDs corresponding to this original ID
        new_ids <- id_mapping$new_id[id_mapping$original_id == id]
        # Get the original genomic data row for this ID
        original_row <- group_data[group_data[[genomic_id_col]] == id, , drop = FALSE]
        
        if (nrow(original_row) > 0) {
          # For each new ID, create a duplicate row with the new ID
          for (new_id in new_ids) {
            new_row <- original_row
            new_row[[genomic_id_col]] <- new_id
            expanded_group_data <- rbind(expanded_group_data, new_row)
          }
        }
      }
      
      # Replace the original genomic data with expanded data
      if (nrow(expanded_group_data) > 0) {
        separated_pathways[[group_name]] <- expanded_group_data
        cat(group_name, "genomic data expanded from", nrow(group_data), "to", nrow(expanded_group_data), "samples\n")
      }
    }
  }
  
  # Update individual dataset references
  hallmark <- separated_pathways[["HALLMARK"]]
  kegg <- separated_pathways[["KEGG"]]
  reactome <- separated_pathways[["REACTOME"]]
  biocarta <- separated_pathways[["BIOCARTA"]]
  
  # Now work with the new IDs for further processing
  final_common_ids <- clinical_ids
} else {
  # Standard processing if no duplicates
  # Find sample IDs that exist in both clinical data and all genomic datasets
  final_common_ids <- intersect(clinical_ids, common_genomic_ids)
}

cat("Final sample IDs for analysis:", length(final_common_ids), "\n")

# Initialize variables to avoid "object not found" errors
clinical_only_ids <- character(0)
genomic_only_ids <- character(0)

# Report removed samples (always calculate these variables)
clinical_only_ids <- setdiff(clinical_ids, final_common_ids)
genomic_only_ids <- setdiff(common_genomic_ids, final_common_ids)

# Only print the first report if we're not handling duplicates
if (!exists("id_mapping") || !has_duplicates) {
  cat("Samples in clinical but not in genomic data:", length(clinical_only_ids), "\n")
  cat("Samples in genomic data but not in clinical:", length(genomic_only_ids), "\n")
}

# Always print these lines since the variables are now defined
cat("Samples in clinical but not in genomic data:", length(clinical_only_ids), "\n")
cat("Samples in genomic data but not in clinical:", length(genomic_only_ids), "\n")

# ---- FILTER ALL DATASETS TO COMMON SAMPLE IDs ----
if (exists("id_mapping") && has_duplicates) {
  # If we've already handled duplicates, we've already created the right set of samples
  # No additional filtering needed for the genomic data
  cat("\nUsing expanded genomic data with duplicate samples matching clinical data\n")
} else {
  # Standard filtering if no duplicates
  # Filter clinical data
  clinical <- clinical[clinical[[clinical_id_col]] %in% final_common_ids, , drop = FALSE]
  
  # Filter genomic datasets
  for (group_name in names(separated_pathways)) {
    group_data <- separated_pathways[[group_name]]
    if (nrow(group_data) > 0 && genomic_id_col %in% colnames(group_data)) {
      filtered_data <- group_data[group_data[[genomic_id_col]] %in% final_common_ids, , drop = FALSE]
      separated_pathways[[group_name]] <- filtered_data
    }
  }
}

# Update individual dataset references
hallmark <- separated_pathways[["HALLMARK"]]
kegg <- separated_pathways[["KEGG"]]
reactome <- separated_pathways[["REACTOME"]]
biocarta <- separated_pathways[["BIOCARTA"]]

# ---- SUMMARY AFTER FILTERING ----
cat("\nFinal harmonized data dimensions:\n")
cat("Clinical:", nrow(clinical), "samples\n")
cat("Hallmark:", nrow(hallmark), "samples x", ncol(hallmark) - 1, "pathways\n")
cat("KEGG:", nrow(kegg), "samples x", ncol(kegg) - 1, "pathways\n")
cat("Reactome:", nrow(reactome), "samples x", ncol(reactome) - 1, "pathways\n")
cat("BioCarta:", nrow(biocarta), "samples x", ncol(biocarta) - 1, "pathways\n")

# Verify that dimensions match between clinical and genomic data
if (nrow(clinical) != nrow(hallmark) && nrow(hallmark) > 0) {
  cat("\nWARNING: Dimension mismatch between clinical and Hallmark genomic data!\n")
  cat("Clinical:", nrow(clinical), "vs. Hallmark:", nrow(hallmark), "\n")
}
if (nrow(clinical) != nrow(kegg) && nrow(kegg) > 0) {
  cat("\nWARNING: Dimension mismatch between clinical and KEGG genomic data!\n")
  cat("Clinical:", nrow(clinical), "vs. KEGG:", nrow(kegg), "\n")
}
if (nrow(clinical) != nrow(reactome) && nrow(reactome) > 0) {
  cat("\nWARNING: Dimension mismatch between clinical and Reactome genomic data!\n")
  cat("Clinical:", nrow(clinical), "vs. Reactome:", nrow(reactome), "\n")
}
if (nrow(clinical) != nrow(biocarta) && nrow(biocarta) > 0) {
  cat("\nWARNING: Dimension mismatch between clinical and BioCarta genomic data!\n")
  cat("Clinical:", nrow(clinical), "vs. BioCarta:", nrow(biocarta), "\n")
}

# ---- WRITE OUTPUT ----
# Extract cancer type from clinical file name for output file naming
clinical_basename <- basename(clinical_file)
cancer_type <- strsplit(clinical_basename, "_")[[1]][1]

cat("\nWriting harmonized datasets...\n")

# Write clinical dataset
clinical_output <- file.path(output_dir, paste0(cancer_type, "_harmonized_clinical.csv"))
fwrite(clinical, clinical_output, sep = ",", quote = TRUE, na = "NA")
cat("Clinical data:", clinical_output, "\n")

# Write genomic datasets (only if they contain data)
genomic_outputs <- list()
for (group_name in names(separated_pathways)) {
  group_data <- separated_pathways[[group_name]]
  group_lower <- tolower(group_name)
  
  output_file <- file.path(output_dir, paste0(cancer_type, "_harmonized_", group_lower, "_genomics.csv"))
  
  if (nrow(group_data) > 0) {
    fwrite(group_data, output_file, sep = ",", quote = TRUE, na = "NA")
    cat(paste0(group_name, " genomics: ", output_file, "\n"))
    genomic_outputs[[group_name]] <- output_file
  } else {
    # Create empty file for consistency
    empty_df <- data.frame(SampleID = character(0))
    fwrite(empty_df, output_file, sep = ",", quote = TRUE, na = "NA")
    cat(paste0(group_name, " genomics (empty): ", output_file, "\n"))
    genomic_outputs[[group_name]] <- output_file
  }
}

cat("\nHarmonization completed successfully!\n")
cat("All datasets now have matching sample IDs for downstream analysis.\n")

# Additional information about duplicate handling if applicable
if (exists("id_mapping") && has_duplicates) {
  cat("\nDuplicate sample handling summary:\n")
  cat("- Found", length(unique(duplicate_ids)), "unique sample IDs that were duplicated\n")
  cat("- Created", length(unique(id_mapping$new_id)) - length(unique(id_mapping$original_id)), "additional unique identifiers\n")
  cat("- Duplicated corresponding genomic data to ensure dimension matching\n")
  cat("- Final clinical samples:", nrow(clinical), "\n")
}

cat("Harmonized files written to:", output_dir, "\n")
cat("Files created with prefix:", cancer_type, "\n")
cat("\nSummary of outputs:\n")
cat("- 1 clinical file\n")
cat("- 4 genomics pathway files (HALLMARK, KEGG, REACTOME, BIOCARTA)\n")
cat("- All files in samples × features format\n")