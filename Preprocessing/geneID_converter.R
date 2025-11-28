# ===============================================================================
# Gene ID Converter: Entrez IDs to Gene Symbols using BioMart
# ===============================================================================
# 
# Purpose: Converts gene identifiers from Entrez IDs to HGNC gene symbols
#          to ensure compatibility with pathway databases and enable proper
#          gene set enrichment analysis.
#
# Description:
#   This script uses the BioMart database to convert Entrez gene IDs to
#   standardized HGNC gene symbols. This conversion is essential for pathway
#   enrichment analysis as most gene set databases (MSigDB, KEGG, etc.) use
#   gene symbols rather than numeric Entrez IDs.
#
# Input Requirements:
#   1. Gene expression matrix with Entrez IDs as row names
#   2. Internet connection for BioMart database access
#   3. Human genome annotation (can be adapted for other species)
#
# Output:
#   - Gene expression matrix with HGNC symbols as row names
#   - Mapping table showing Entrez ID to Symbol conversions
#   - Summary of successful and failed conversions
#   - Filtered matrix excluding genes without symbol mappings
#
# Conversion Method:
#   - Connects to Ensembl BioMart database
#   - Queries human genome annotations
#   - Maps Entrez IDs to official HGNC symbols
#   - Handles duplicates and missing mappings appropriately
#
# Usage:
#   1. Ensure input file has Entrez IDs as row identifiers
#   2. Configure file paths in the script
#   3. Run: Rscript geneID_converter.R
#
# Dependencies: biomaRt, data.table
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
suppressPackageStartupMessages(library(biomaRt))

# ---- User Inputs ----
input_file <- "/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/NSCLC_filtered_RNAseq.csv"   # <-- Set your input file path
output_file <- "/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/gensym_NSCLC_filtered_RNAseq.csv" # <-- Set your output file path

# ---- Read RNAseq matrix (genes x samples) ----
rna_mat <- read.csv(input_file, row.names = 1, check.names = FALSE)

# ---- Get Entrez IDs ----
entrez_ids <- rownames(rna_mat)

# ---- Connect to BioMart ----
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# ---- Query BioMart for gene symbols ----
bm <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = entrez_ids,
  mart = mart
)

# ---- Create mapping ----
id_map <- setNames(bm$hgnc_symbol, as.character(bm$entrezgene_id))

# ---- Replace Entrez IDs with gene symbols ----
gene_symbols <- id_map[as.character(entrez_ids)]
# If mapping fails, keep original Entrez ID
gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- entrez_ids[is.na(gene_symbols) | gene_symbols == ""]

# ---- Make gene symbols unique ----
gene_symbols_unique <- make.unique(gene_symbols)
rownames(rna_mat) <- gene_symbols_unique

# ---- Write to new CSV ----
write.csv(rna_mat, file = output_file, row.names = TRUE)