library(data.table)

extract_data <- function(input_dir, metadata_file, output_file) {
  # Read metadata file to map File Name to Case ID
  metadata <- fread(metadata_file, sep = "\t")
  metadata[, `Case ID` := sapply(`Case ID`, function(x) strsplit(x, ",")[[1]][1])]  # Take the first ID
  file_to_case_id <- setNames(metadata$`Case ID`, metadata$`File Name`)

  # Initialize the final dataframe
  final_df <- data.table(gene_name = character())

  # Track processed sample IDs
  processed_samples <- character()

  # Loop through subfolders in the input directory
  folders <- list.dirs(input_dir, recursive = FALSE)
  for (folder in folders) {
    # Ensure folder exists and contains files
    if (!dir.exists(folder)) next

    files <- list.files(folder, pattern = "\\.tsv$", full.names = TRUE)
    for (file_path in files) {
      file <- basename(file_path)
      sample_id <- ifelse(file %in% names(file_to_case_id), file_to_case_id[file], file)

      # Skip duplicate samples
      if (sample_id %in% processed_samples) {
        next
      }
      processed_samples <- c(processed_samples, sample_id)

      # Read the .tsv file, skipping the first row
      data <- tryCatch(
        fread(file_path, skip = 1, select = c("gene_name", "gene_type", "tpm_unstranded")),
        error = function(e) {
          cat(sprintf("Error reading file %s: %s\n", file_path, e$message))
          return(NULL)
        }
      )
      if (is.null(data)) next

      # Ensure required columns exist
      required_columns <- c("gene_name", "gene_type", "tpm_unstranded")
      if (!all(required_columns %in% colnames(data))) {
        cat(sprintf("Skipping file %s: Missing required columns.\n", file_path))
        next
      }

      # Filter rows where gene_type is "protein_coding" and gene_name is not null
      filtered_data <- data[gene_type == "protein_coding" & !is.na(gene_name)]

      # Skip if no valid rows are found
      if (nrow(filtered_data) == 0) {
        cat(sprintf("No valid rows found in file %s.\n", file_path))
        next
      }

      # Extract gene_name and tpm_unstranded columns
      extracted_data <- filtered_data[, .(gene_name, tpm_unstranded)]
      setnames(extracted_data, "tpm_unstranded", sample_id)

      # Remove duplicate gene_names
      extracted_data <- extracted_data[!duplicated(extracted_data$gene_name), ]

      # Merge with the final dataframe
      final_df <- merge(final_df, extracted_data, by = "gene_name", all = TRUE)
    }
  }

  # Remove duplicate gene_names in final_df
  final_df <- final_df[!duplicated(final_df$gene_name), ]

  # Save the final dataframe to a CSV file
  if (nrow(final_df) > 0) {
    fwrite(final_df, output_file)
    cat(sprintf("Output file saved to: %s\n", output_file))
  } else {
    cat("No data was parsed. The output file will be blank.\n")
  }
}

# Example usage
input_dir <- "/Users/jackie-mac/Desktop/VSCode/data/BRCA_raw_RNAseq/"
metadata_file <- "/Users/jackie-mac/Desktop/VSCode/data/BRCA_sample_sheet.2025-06-18.tsv"
output_file <- "/Users/jackie-mac/Desktop/VSCode/outputs/filtered_RNAseq/Protein_coding/BRCA_filtered_RNAseq.csv"

extract_data(input_dir, metadata_file, output_file)