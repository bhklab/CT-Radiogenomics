# =========================
# Feature Filtering Script: Remove Highly Correlated Features (by Pathway Group or Radiomics)
# =========================

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Parse Command Line Arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop("Usage: Rscript feature_filtering.R <kegg_cor_mat> <hallmark_cor_mat> <reactome_cor_mat> <biocarta_cor_mat> <radiomics_cor_mat> <genomics_data> <radiomics_data> <output_dir> <prefix>")
}
cor_mat_files <- setNames(args[1:4], c("KEGG", "HALLMARK", "REACTOME", "BIOCARTA"))
radiomics_cor_mat <- args[5]
genomics_data_file <- args[6]
radiomics_data_file <- args[7]
output_dir <- args[8]
prefix <- args[9]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

genomics_data <- fread(genomics_data_file, data.table = FALSE, check.names = FALSE)
radiomics_data <- fread(radiomics_data_file, data.table = FALSE, check.names = FALSE)

# ---- Genomics Feature Filtering by Group ----
for (group in names(cor_mat_files)) {
  cat("Processing group:", group, "\n")
  cor_mat_file <- cor_mat_files[[group]]
  cat("  Reading correlation matrix from:", cor_mat_file, "\n")
  cor_mat <- read.csv(cor_mat_file, row.names = 1, check.names = FALSE)
  cor_mat <- as.matrix(cor_mat)
  cat("  Correlation matrix dimensions:", dim(cor_mat), "\n")

  # ---- Collect highly correlated feature pairs (|correlation| > 0.9, excluding self-correlation) ----
  results <- data.frame(
    Feature1 = character(),
    Feature2 = character(),
    Correlation = numeric(),
    stringsAsFactors = FALSE
  )

  cat("  Identifying highly correlated feature pairs...\n")
  for (i in 1:(nrow(cor_mat) - 1)) {
    for (j in (i + 1):ncol(cor_mat)) {
      val <- cor_mat[i, j]
      if (!is.na(val) && abs(val) > 0.9 && abs(val) < 1) {
        results <- rbind(
          results,
          data.frame(
            Feature1 = rownames(cor_mat)[i],
            Feature2 = colnames(cor_mat)[j],
            Correlation = val,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  cat("    Found", nrow(results), "highly correlated pairs.\n")

  # Write all highly correlated pairs to CSV
  out_corr_file <- file.path(output_dir, paste0(prefix, "_", group, "_0.9correlation.csv"))
  write.csv(
    results,
    file = out_corr_file,
    row.names = FALSE
  )
  cat("    Written highly correlated pairs to:", out_corr_file, "\n")

  # ---- Create a subset keeping only one feature from each pair ----
  removed <- character(0)
  subset_results <- data.frame(
    Representative = character(),
    Removed = character(),
    Correlation = numeric(),
    stringsAsFactors = FALSE
  )

  cat("  Filtering features...\n")
  for (k in seq_len(nrow(results))) {
    f1 <- results$Feature1[k]
    f2 <- results$Feature2[k]
    val <- results$Correlation[k]
    if (f2 %in% removed && !(f1 %in% removed)) {
      subset_results <- rbind(subset_results, data.frame(
        Representative = f2,
        Removed = f1,
        Correlation = val,
        stringsAsFactors = FALSE
      ))
      removed <- c(removed, f1)
    } else if (!(f1 %in% removed) && !(f2 %in% removed)) {
      subset_results <- rbind(subset_results, data.frame(
        Representative = f1,
        Removed = f2,
        Correlation = val,
        stringsAsFactors = FALSE
      ))
      removed <- c(removed, f2)
    }
  }

  to_remove_unique <- unique(subset_results$Removed)
  cat("    Features to remove for", group, ":", length(to_remove_unique), "\n")

  # Remove rows corresponding to removed features from the original genomics data
  # Assume first column is pathway name, rest are samples
  pathway_col <- 1
  keep_pathways <- setdiff(genomics_data[[pathway_col]], to_remove_unique)

  # Now, filter to only pathways belonging to this group
  group_pattern <- paste0("^", group, "_")
  group_pathways <- grep(group_pattern, keep_pathways, value = TRUE)
  filtered_genomics <- genomics_data[genomics_data[[pathway_col]] %in% group_pathways, , drop = FALSE]

  # Transpose so samples are rows and pathways are columns
  rownames(filtered_genomics) <- filtered_genomics[[1]]
  filtered_genomics <- filtered_genomics[, -1, drop = FALSE]
  filtered_genomics_t <- as.data.frame(t(filtered_genomics))
  filtered_genomics_t <- cbind(SampleID = rownames(filtered_genomics_t), filtered_genomics_t)
  rownames(filtered_genomics_t) <- NULL

  out_filtered_file <- file.path(output_dir, paste0(prefix, "_", group, "_features_filtered.csv"))
  write.csv(
    filtered_genomics_t,
    file = out_filtered_file,
    row.names = FALSE
  )
  cat("    Filtered genomics data written for", group, "to:", out_filtered_file, "\n")
}

# ---- Radiomics Feature Filtering ----
cat("Processing radiomics features...\n")
cat("  Reading radiomics correlation matrix from:", radiomics_cor_mat, "\n")
radiomics_mat <- read.csv(radiomics_cor_mat, row.names = 1, check.names = FALSE)
radiomics_mat <- as.matrix(radiomics_mat)
cat("  Radiomics correlation matrix dimensions:", dim(radiomics_mat), "\n")

# ---- Collect highly correlated radiomics feature pairs (|correlation| > 0.9, excluding self-correlation) ----
results <- data.frame(
  Feature1 = character(),
  Feature2 = character(),
  Correlation = numeric(),
  stringsAsFactors = FALSE
)

cat("  Identifying highly correlated radiomics feature pairs...\n")
for (i in 1:(nrow(radiomics_mat) - 1)) {
  for (j in (i + 1):ncol(radiomics_mat)) {
    val <- radiomics_mat[i, j]
    if (!is.na(val) && abs(val) > 0.9 && abs(val) < 1) {
      results <- rbind(
        results,
        data.frame(
          Feature1 = rownames(radiomics_mat)[i],
          Feature2 = colnames(radiomics_mat)[j],
          Correlation = val,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}
cat("    Found", nrow(results), "highly correlated radiomics pairs.\n")

# Write all highly correlated pairs to CSV
out_corr_file <- file.path(output_dir, paste0(prefix, "_radiomics_0.9correlation.csv"))
write.csv(
  results,
  file = out_corr_file,
  row.names = FALSE
)
cat("    Written highly correlated radiomics pairs to:", out_corr_file, "\n")

# ---- Create a subset keeping only one radiomics feature from each pair ----
removed <- character(0)
subset_results <- data.frame(
  Representative = character(),
  Removed = character(),
  Correlation = numeric(),
  stringsAsFactors = FALSE
)

cat("  Filtering radiomics features...\n")
for (k in seq_len(nrow(results))) {
  f1 <- results$Feature1[k]
  f2 <- results$Feature2[k]
  val <- results$Correlation[k]
  if (f2 %in% removed && !(f1 %in% removed)) {
    subset_results <- rbind(subset_results, data.frame(
      Representative = f2,
      Removed = f1,
      Correlation = val,
      stringsAsFactors = FALSE
    ))
    removed <- c(removed, f1)
  } else if (!(f1 %in% removed) && !(f2 %in% removed)) {
    subset_results <- rbind(subset_results, data.frame(
      Representative = f1,
      Removed = f2,
      Correlation = val,
      stringsAsFactors = FALSE
    ))
    removed <- c(removed, f2)
  }
}

to_remove_unique <- unique(subset_results$Removed)
cat("    Radiomics features to remove:", length(to_remove_unique), "\n")

# ---- Remove columns corresponding to removed features from the original radiomics data ----
# Assume first column is sample ID, rest are features
feature_cols <- 2:ncol(radiomics_data)
keep_feature_names <- setdiff(colnames(radiomics_data)[feature_cols], to_remove_unique)

# Only keep features with allowed prefixes
allowed_prefixes <- c("original_", "wavelet-", "square_", "squareroot_", "logarithm_", "exponential_", "gradient_")
prefix_pattern <- paste0("^(", paste(allowed_prefixes, collapse = "|"), ")")
filtered_feature_names <- grep(prefix_pattern, keep_feature_names, value = TRUE)

# Combine sample ID column with filtered feature columns
filtered_radiomics <- radiomics_data[, c(1, which(colnames(radiomics_data) %in% filtered_feature_names)), drop = FALSE]

out_filtered_file <- file.path(output_dir, paste0(prefix, "_radiomics_features_filtered.csv"))
write.csv(
  filtered_radiomics,
  file = out_filtered_file,
  row.names = FALSE
)
cat("    Filtered radiomics data written to:", out_filtered_file, "\n")

cat("Feature filtering complete.\n")