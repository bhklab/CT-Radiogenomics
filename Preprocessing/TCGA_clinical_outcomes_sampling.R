suppressPackageStartupMessages(library(data.table))
library(tools)

# ---- USER INPUTS ----
dataset_files <- c(
  "/Users/jackie-mac/Desktop/VSCode/data/clinical/Linkedomics_clinical/Human__TCGA_BRCA__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
  "/Users/jackie-mac/Desktop/VSCode/data/clinical/Linkedomics_clinical/Human__TCGA_ESCA__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
  "/Users/jackie-mac/Desktop/VSCode/data/clinical/Linkedomics_clinical/Human__TCGA_GBM__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
  "/Users/jackie-mac/Desktop/VSCode/data/clinical/Linkedomics_clinical/Human__TCGA_KIRC__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
  "/Users/jackie-mac/Desktop/VSCode/data/clinical/Linkedomics_clinical/Human__TCGA_LGG__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
  "/Users/jackie-mac/Desktop/VSCode/data/clinical/Linkedomics_clinical/Human__TCGA_LIHC__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi"
)
sample_id_files <- c(
  "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/BRCA_patientids.txt",
  "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/ESCA_patientids.txt",
  "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/GBM_patientids.txt",
  "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/KIRC_patientids.txt",
  "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/LGG_patientids.txt",
  "/Users/jackie-mac/Desktop/VSCode/data/sample_ids/LIHC_patientids.txt"
)
output_file <- "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/2combined_TCGA_clinical_outcomes.csv"

# ---- MATCH DATASETS TO SAMPLE ID FILES ----
get_name_from_patientids <- function(filename) {
  toupper(sub("_patientids.*$", "", file_path_sans_ext(basename(filename))))
}
get_name_from_dataset <- function(filename) {
  # Extract three or four letters after "TCGA_"
  m <- regexpr("TCGA_([A-Z]{3,4})", filename)
  if (m[1] != -1) {
    return(toupper(sub(".*TCGA_([A-Z]{3,4}).*", "\\1", filename)))
  }
  return(NA_character_)
}

sample_id_names <- sapply(sample_id_files, get_name_from_patientids)
dataset_names <- sapply(dataset_files, get_name_from_dataset)
matched_pairs <- list()

for (i in seq_along(dataset_files)) {
  dataset_name <- dataset_names[i]
  match_idx <- which(sample_id_names == dataset_name)
  if (length(match_idx) > 0) {
    matched_pairs[[dataset_files[i]]] <- sample_id_files[match_idx[1]]
  }
}

# ---- EXTRACT DATA ----
output_list <- list()
row_counts <- integer(length(matched_pairs))

i <- 1
for (dataset_path in names(matched_pairs)) {
  sample_id_file <- matched_pairs[[dataset_path]]
  sample_ids_raw <- scan(sample_id_file, what = character(), quiet = TRUE)
  sample_ids <- toupper(trimws(sample_ids_raw))
  dt <- fread(dataset_path, sep = "\t", header = TRUE, na.strings = c("NA", ""))
  attr_col <- names(dt)[1]
  os_row <- dt[get(attr_col) == "overall_survival"]
  if (nrow(os_row) == 0) {
    warning(paste("No 'overall_survival' row found in", dataset_path))
    next
  }
  os_row_samples <- os_row[, -1, with = FALSE]
  all_sample_ids_raw <- colnames(os_row_samples)
  all_sample_ids <- toupper(trimws(all_sample_ids_raw))
  matched_idx <- which(all_sample_ids %in% sample_ids)
  matched_samples <- all_sample_ids_raw[matched_idx]
  missing_samples <- setdiff(sample_ids, all_sample_ids)
  if (length(missing_samples) > 0) {
    cat("Samples not found in dataset", dataset_path, ":\n")
    cat(missing_samples, sep = "\n")
    cat("\n")
  }
  if (length(matched_samples) == 0) {
    cat("No matched samples for", dataset_path, "\n")
    next
  }
  os_values <- as.character(os_row_samples[, ..matched_samples])
  name_prefix <- get_name_from_dataset(dataset_path)
  dt_out <- data.table(
    SampleID = matched_samples,
    overall_survival = as.vector(os_values)
  )
  setnames(dt_out, c(paste0(name_prefix, "_SampleID"), paste0(name_prefix, "_overall_survival")))
  row_counts[i] <- nrow(dt_out)
  output_list[[i]] <- dt_out
  i <- i + 1
}

# ---- PAD AND COMBINE WITH SPACERS ----
max_rows <- max(row_counts, 1)
combined_list <- list()
for (i in seq_along(output_list)) {
  dt <- output_list[[i]]
  n_pad <- max_rows - nrow(dt)
  if (n_pad > 0) {
    pad_dt <- as.data.table(matrix(NA, nrow = n_pad, ncol = ncol(dt)))
    setnames(pad_dt, names(dt))
    dt <- rbind(dt, pad_dt)
  }
  combined_list[[length(combined_list) + 1]] <- dt
  if (i != length(output_list)) {
    combined_list[[length(combined_list) + 1]] <- data.table(" " = rep("", max_rows))
  }
}

if (length(combined_list) == 0) {
  cat("No data to write. No matched samples found in any dataset.\n")
} else {
  final_output <- as.data.table(do.call(cbind, combined_list))
  fwrite(final_output, output_file, na = "NA", quote = TRUE)
  cat("Extraction complete. Output written to:", output_file, "\n")
}