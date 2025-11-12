# =========================
# Filter for protein-coding genes using BioMart
# =========================

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
suppressPackageStartupMessages(library(biomaRt))

# ---- User Inputs ----
input_files <- c(
  "/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/Protein_coding/BRCA_filtered_RNAseq.csv",
  "/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/Protein_coding/LGG_filtered_RNAseq.csv",
  "/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/Protein_coding/GBM_filtered_RNAseq.csv",
  "/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/Protein_coding/KIRC_filtered_RNAseq.csv",
  "/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/Protein_coding/CCRCC_filtered_RNAseq.csv"
  #"/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/Protein_coding/HNSCC_filtered_RNAseq.csv",
  #"/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/Protein_coding/PDA_filtered_RNAseq.csv"
  #"/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/PRotein_coding/NSCLC_filtered_RNAseq.csv"
)
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/Protein_coding/MissingRemoved/" # Output directory

# ---- Connect to BioMart ----
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl", 
                host = 'https://www.ensembl.org')

# ---- Query BioMart for protein-coding genes ----
genes <- getBM(
  attributes = c("external_gene_name", "transcript_biotype"),
  filters = "transcript_biotype",
  values = "protein_coding",
  mart = mart
)

# ---- Filter for protein-coding genes ----
protein_coding_genes <- unique(genes$external_gene_name)

# ---- Loop through each input file ----
for (input_file in input_files) {
  # Extract cancer type from the filename (assumes filename format like "BRCA_filtered_RNAseq.csv")
  cancer_type <- strsplit(basename(input_file), "_")[[1]][1]
  
  # Define output file path
  output_file <- file.path(output_dir, paste0(cancer_type, "_filtered_TPM.csv"))
  
  # ---- Read RNAseq matrix (genes x samples, gene symbols as rownames) ----
  rna_mat <- read.csv(input_file, row.names = 1, check.names = FALSE)
  gene_symbols <- rownames(rna_mat)
  
  # ---- Subset the RNAseq matrix ----
  rna_mat_protein_coding <- rna_mat[rownames(rna_mat) %in% protein_coding_genes, , drop = FALSE]
  
  # ---- Filter genes based on expression data ----
  # Calculate for each gene (row) the proportion of samples with missing or zero expression
  missing_or_zero_prop <- apply(rna_mat_protein_coding, 1, function(x) {
    mean(is.na(x) | x == 0)
  })
  
  # Keep only genes with <50% missing/zero values across samples
  rna_mat_filtered <- rna_mat_protein_coding[missing_or_zero_prop < 0.5, , drop = FALSE]
  
  # ---- Write to new CSV ----
  write.csv(rna_mat_filtered, file = output_file, row.names = TRUE)
  cat("Filtered data for ", cancer_type, " saved to: ", output_file, "\n")
}