
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(fgsea))

perform_kegg_samplewise_fgsea <- function(input_file, gmt_file, output_file) {
  cat("Loading dataset from:", input_file, "\n")
  expression_data <- fread(input_file, data.table = FALSE, check.names = FALSE)
  
  # Set rownames from the first column and remove it
  rownames(expression_data) <- expression_data[[1]]
  expression_data[[1]] <- NULL
  
  # Read KEGG gene sets from GMT file
  gmt <- getGmt(gmt_file)
  gene_set_list <- geneIds(gmt)
  pathway_names <- names(gene_set_list)
  
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
  pathway_names <- names(filtered_gene_set_list) # Update pathway names after filtering
  cat("Number of gene sets after filtering for >=80% gene coverage:", length(filtered_gene_set_list), "\n")
  
  # Prepare results data frame: samples as rows, pathways as columns
  results <- data.frame(matrix(nrow = ncol(expression_data), ncol = length(pathway_names)))
  colnames(results) <- pathway_names
  rownames(results) <- colnames(expression_data)
  
  # For each sample, run fgsea using gene expression as ranking
  for (sample_idx in seq_along(colnames(expression_data))) {
    sample_name <- colnames(expression_data)[sample_idx]
    cat("Processing sample:", sample_idx, ":", sample_name, "\n")
    ranks <- expression_data[, sample_idx]
    names(ranks) <- rownames(expression_data)
    ranks <- sort(ranks, decreasing = TRUE)
    fgseaRes <- fgsea(pathways = filtered_gene_set_list, stats = ranks)
    # Align ES to the correct pathway order
    es_vec <- fgseaRes$ES[match(pathway_names, fgseaRes$pathway)]
    results[sample_name, ] <- es_vec
  }
  
  # Save results
  write.csv(results, output_file, row.names = TRUE, quote = FALSE)
  cat("Sample-wise KEGG fgsea results saved to:", output_file, "\n")
}

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript kegg_enrichment_GMT.R <input_file> <gmt_file> [output_file]")
}
input_file <- args[1]
gmt_file <- args[2]
output_file <- ifelse(length(args) >= 3, args[3], "kegg_fgsea_results.csv")

cat("Starting KEGG sample-wise fgsea analysis...\n")
perform_kegg_samplewise_fgsea(input_file, gmt_file, output_file)
cat("Process completed.\n")

# Rscript kegg_enrichment_GMT.R /Users/jackie-mac/Desktop/VSCode/outputs/filtered_RNAseq/Protein_coding/BRCA_filtered.csv /Users/jackie-mac/Desktop/VSCode/data/GMT_files/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt /Users/jackie-mac/Desktop/VSCode/outputs/enrichment_scores/BRCA_kegg_fgsea_results.csv
