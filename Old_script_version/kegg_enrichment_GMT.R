# ===============================================================================
# KEGG Pathway Enrichment Analysis using GSVA
# ===============================================================================
# 
# Purpose: Performs Gene Set Variation Analysis (GSVA) using KEGG pathway 
#          gene sets to calculate pathway enrichment scores for each sample.
#
# Description:
#   This script takes RNA-seq gene expression data and calculates sample-wise
#   enrichment scores for KEGG (Kyoto Encyclopedia of Genes and Genomes) pathways.
#   GSVA transforms gene-level expression data into pathway-level scores, enabling
#   pathway-centric analysis and reducing dimensionality.
#
# Input Requirements:
#   1. Gene expression matrix: Genes as rows, samples as columns
#   2. KEGG pathway gene sets: GMT format file containing pathway definitions
#   3. Gene identifiers must match between expression data and pathway definitions
#
# Output:
#   - Pathway enrichment score matrix: Pathways as rows, samples as columns
#   - Enrichment scores represent pathway activity levels for each sample
#   - Scores can be positive (upregulated) or negative (downregulated)
#   - Results saved as CSV for downstream radiogenomic correlation analysis
#
# Analysis Method:
#   - Uses GSVA algorithm for single-sample pathway enrichment
#   - Employs empirical cumulative distribution functions
#   - Normalizes scores across samples for comparative analysis
#   - Handles missing genes gracefully
#
# Usage:
#   1. Ensure KEGG GMT file is available and properly formatted
#   2. Configure input file paths in the script
#   3. Run: Rscript kegg_enrichment_GMT.R
#
# Dependencies: data.table, GSEABase, GSVA
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(GSVA))

perform_kegg_samplewise_gsva <- function(input_file, gmt_file, output_file) {
  cat("Loading dataset from:", input_file, "\n")
  expression_data <- fread(input_file, data.table = FALSE, check.names = FALSE)
  
  # Set rownames from the first column and remove it
  rownames(expression_data) <- expression_data[[1]]
  expression_data[[1]] <- NULL
  
  # Read KEGG gene sets from GMT file
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
  cat("Sample-wise KEGG GSVA results saved to:", output_file, "\n")
}

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript kegg_enrichment_GMT.R <input_file> <gmt_file> [output_file]")
}
input_file <- args[1]
gmt_file <- args[2]
output_file <- ifelse(length(args) >= 3, args[3], "kegg_gsva_results.csv")

cat("Starting KEGG sample-wise GSVA analysis...\n")
perform_kegg_samplewise_gsva(input_file, gmt_file, output_file)
cat("Process completed.\n")

# Rscript kegg_enrichment_GMT.R /Users/jackie-mac/Desktop/VSCode/outputs/filtered_RNAseq/Protein_coding/BRCA_filtered.csv /Users/jackie-mac/Desktop/VSCode/data/GMT_files/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt /Users/jackie-mac/Desktop/VSCode/outputs/enrichment_scores/BRCA_kegg_fgsea_results.csv
