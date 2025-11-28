# ===============================================================================
# BioCarta Pathway Enrichment Analysis using GSVA
# ===============================================================================
# 
# Purpose: Performs Gene Set Variation Analysis (GSVA) using BioCarta pathway
#          gene sets to calculate enrichment scores for cellular and molecular
#          pathways with focus on signal transduction networks.
#
# Description:
#   This script calculates sample-wise enrichment scores for BioCarta pathways,
#   which represent manually curated pathway maps focusing on molecular 
#   interactions and signal transduction cascades. BioCarta pathways are
#   particularly useful for understanding cellular communication networks.
#
# Input Requirements:
#   1. Gene expression matrix: Genes as rows, samples as columns
#   2. BioCarta pathway gene sets: GMT format file containing pathway definitions
#   3. Gene symbols must match between expression data and pathway annotations
#
# Output:
#   - BioCarta pathway enrichment score matrix: Pathways as rows, samples as columns
#   - Enrichment scores represent pathway activity levels for each sample
#   - Positive scores indicate pathway activation, negative scores indicate inhibition
#   - Results formatted for downstream radiogenomic correlation analysis
#
# Analysis Method:
#   - Uses GSVA algorithm for single-sample pathway enrichment scoring
#   - Employs rank-based statistics for robust score calculation
#   - Normalizes scores to enable cross-sample and cross-pathway comparisons
#   - Focuses on signal transduction and metabolic pathways
#
# Usage:
#   Rscript biocarta_enrichment_GMT.R <input_expression_file> <biocarta_gmt_file> <output_file>
#
# Dependencies: data.table, GSEABase, GSVA
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(GSVA))

perform_biocarta_samplewise_gsva <- function(input_file, gmt_file, output_file) {
  cat("Loading dataset from:", input_file, "\n")
  expression_data <- fread(input_file, data.table = FALSE, check.names = FALSE)
  
  # Set rownames from the first column and remove it
  rownames(expression_data) <- expression_data[[1]]
  expression_data[[1]] <- NULL
  
  # Read BioCarta gene sets from GMT file
  gmt <- getGmt(gmt_file)
  gene_set_list <- geneIds(gmt)
  
  # ---- Filter gene sets: keep only those with >=80% genes present in RNAseq data ----
  rnaseq_genes <- rownames(expression_data)
  filtered_gene_set_list <- lapply(gene_set_list, function(genes) {
    intersect(genes, rnaseq_genes)
  })
  filtered_gene_set_list <- filtered_gene_set_list[
    sapply(seq_along(filtered_gene_set_list), function(i) {
      length(filtered_gene_set_list[[i]]) / length(gene_set_list[[i]]) >= 0.8
    })
  ]

  cat("Number of gene sets after filtering for >=80% gene coverage:", length(filtered_gene_set_list), "\n")
  
  # Prepare expression matrix for GSVA (genes as rows, samples as columns)
  expr_mat <- as.matrix(expression_data)
  gsvaPar <- gsvaParam(expr_mat, filtered_gene_set_list)
  print(gsvaPar)
  cat("Running GSVA using gsvaParam object...\n")
  gsva_results <- gsva(gsvaPar, verbose = TRUE)
  results <- t(gsva_results)
  colnames(results) <- rownames(gsva_results)  # Gene set names become column names
  rownames(results) <- colnames(gsva_results)  # Sample names become row names
  
  # Debug: Print first 5 rows and columns of results matrix
  cat("Debug: First 5x5 of results matrix:\n")
  n_rows <- min(5, nrow(results))
  n_cols <- min(5, ncol(results))
  if (n_rows > 0 && n_cols > 0) {
    print(as.matrix(results[seq_len(n_rows), seq_len(n_cols)]))
  }
  
  write.csv(results, output_file, row.names = TRUE, quote = FALSE)
  cat("Sample-wise GSVA results saved to:", output_file, "\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript biocarta_enrichment_GMT.R <input_file> <gmt_file> [output_file]")
}
input_file <- args[1]
gmt_file <- args[2]
output_file <- ifelse(length(args) >= 3, args[3], "biocarta_gsva_GMT_results.csv")

cat("Starting sample-wise BioCarta GSVA analysis...\n")
perform_biocarta_samplewise_gsva(input_file, gmt_file, output_file)
cat("Process completed.\n")