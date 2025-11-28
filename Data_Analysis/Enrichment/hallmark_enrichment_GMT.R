# ===============================================================================
# MSigDB Hallmark Pathway Enrichment Analysis using GSVA
# ===============================================================================
# 
# Purpose: Performs Gene Set Variation Analysis (GSVA) using MSigDB Hallmark 
#          gene sets to calculate pathway enrichment scores representing key
#          biological processes and cancer-related pathways.
#
# Description:
#   This script calculates sample-wise enrichment scores for MSigDB Hallmark
#   pathways, which represent well-defined biological processes with coherent
#   expression patterns. Hallmark pathways are refined from canonical pathways
#   to reduce redundancy and focus on core biological processes.
#
# Input Requirements:
#   1. Gene expression matrix: Genes as rows, samples as columns
#   2. MSigDB Hallmark gene sets: GMT format file (h.all.v7.4.symbols.gmt)
#   3. Gene symbols must match between expression data and pathway definitions
#
# Output:
#   - Hallmark pathway enrichment score matrix: Pathways as rows, samples as columns
#   - Enrichment scores indicate pathway activity levels for each sample
#   - Positive scores indicate pathway activation, negative scores indicate suppression
#   - Results formatted for downstream radiogenomic correlation analysis
#
# Analysis Method:
#   - Uses GSVA algorithm for robust single-sample enrichment scoring
#   - Employs kernel density estimation for score calculation
#   - Normalizes scores to enable cross-sample comparisons
#   - Focuses on 50 well-characterized Hallmark pathways
#
# Usage:
#   Rscript hallmark_enrichment_GMT.R <input_expression_file> <hallmark_gmt_file> <output_file>
#
# Dependencies: data.table, GSEABase, GSVA
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(GSVA))

perform_hallmark_samplewise_fgsea <- function(input_file, gmt_file, output_file) {
  cat("Loading dataset from:", input_file, "\n")
  expression_data <- fread(input_file, data.table = FALSE, check.names = FALSE)
  
  # Set rownames from the first column and remove it
  rownames(expression_data) <- expression_data[[1]]
  expression_data[[1]] <- NULL
  
  # Read Hallmark gene sets from GMT file
  gmt <- getGmt(gmt_file)
  gene_set_list <- geneIds(gmt)
  
  # ---- Filter gene sets: keep only those with >=80% genes present in RNAseq data ----
  rnaseq_genes <- rownames(expression_data)
  n.cutoff <- 15  # Example minimum sample count
  filtered_gene_set_list <- lapply(gene_set_list, function(genes) {
    intersect(genes, rnaseq_genes)
  })
  filtered_gene_set_list <- filtered_gene_set_list[
    sapply(seq_along(filtered_gene_set_list), function(i) {
      overlap <- length(filtered_gene_set_list[[i]])
      original_size <- length(gene_set_list[[i]])
      proportion <- ifelse(original_size == 0, 0, overlap / original_size)
      proportion >= 0.8 & ncol(expression_data) >= n.cutoff
    })
  ]
  cat("Number of gene sets after filtering for >=80% gene coverage:", length(filtered_gene_set_list), "\n")
  
  # Prepare expression matrix for GSVA (genes as rows, samples as columns)
  expr_mat <- as.matrix(expression_data)
  # Set up GSVA parameters using gsvaParam with expression matrix and gene sets
  gsvaPar <- gsvaParam(expr_mat, filtered_gene_set_list)
  print(gsvaPar)
  cat("Running GSVA using gsvaParam object...\n")
  gsva_results <- gsva(gsvaPar, verbose = TRUE)
  # Transpose to match previous output: samples as rows, pathways as columns
  results <- t(gsva_results)
  # Ensure column and row names are correct
  colnames(results) <- rownames(gsva_results)  # Gene set names become column names
  rownames(results) <- colnames(gsva_results)  # Sample names become row names
  
  # Debug: Print first 5 rows and columns of results matrix
  cat("Debug: First 5x5 of results matrix:\n")
  n_rows <- min(5, nrow(results))
  n_cols <- min(5, ncol(results))
  if (n_rows > 0 && n_cols > 0) {
    print(as.matrix(results[seq_len(n_rows), seq_len(n_cols)]))
  }

  # Save results
  write.csv(results, output_file, row.names = TRUE, quote = FALSE)
  cat("Sample-wise fgsea results saved to:", output_file, "\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript hallmark_enrichmenet_GMT.R <input_file> <gmt_file> [output_file]")
}
input_file <- args[1]
gmt_file <- args[2]
output_file <- ifelse(length(args) >= 3, args[3], "hallmark_fgsea_results.csv")

cat("Starting sample-wise hallmark fgsea analysis...\n")
perform_hallmark_samplewise_fgsea(input_file, gmt_file, output_file)
cat("Process completed.\n")

# Rscript hallmark_enrichment_GMT.R /Users/jackie-mac/Desktop/VSCode/outputs/filtered_RNAseq/Protein_coding/BRCA_filtered.csv /Users/jackie-mac/Desktop/VSCode/data/GMT_files/h.all.v2024.1.Hs.symbols.gmt /Users/jackie-mac/Desktop/VSCode/outputs/enrichment_scores/hallmark_fgsea_results.csv
