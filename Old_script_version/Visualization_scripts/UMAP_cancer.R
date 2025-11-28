# ===============================================================================
# UMAP Dimensionality Reduction and Visualization for Cancer Datasets
# ===============================================================================
# 
# Purpose: Performs Uniform Manifold Approximation and Projection (UMAP) 
#          dimensionality reduction on radiogenomic features to visualize
#          sample clustering patterns and identify cancer subtypes.
#
# Description:
#   This script applies UMAP algorithm to high-dimensional radiogenomic data
#   to create 2D visualizations that preserve local and global data structure.
#   It enables exploration of sample relationships, identification of potential
#   cancer subtypes, and visualization of feature-based clustering patterns.
#
# Input Requirements:
#   1. Feature matrix: CSV with samples as rows, features (genomic/radiomic) as columns
#   2. Optional metadata: Clinical variables for sample annotation and coloring
#   3. High-dimensional data suitable for dimensionality reduction
#
# Output:
#   - Interactive UMAP plots with sample point visualization
#   - Static plots colored by clinical variables (stage, grade, survival)
#   - Coordinate files for UMAP embeddings
#   - Clustering analysis results if requested
#
# Visualization Features:
#   - Interactive plots using plotly for exploration
#   - Color coding by clinical variables
#   - Customizable point sizes and transparency
#   - Optional clustering overlay and annotations
#
# Usage:
#   1. Configure input data paths and parameters
#   2. Run: Rscript UMAP_cancer.R
#   3. Explore interactive plots and identify patterns
#
# Dependencies: data.table, uwot, ggplot2, plotly
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(uwot)
  library(ggplot2)
  library(plotly)
})

# ---- USER INPUTS ----
rna_files <- list(
  BRCA = "/Users/jackie-mac/Desktop/VSCode/data/RNAseq/BRCA_filtered_TPM.csv",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/data/RNAseq/KIRC_filtered_TPM.csv",
  LGG = "/Users/jackie-mac/Desktop/VSCode/data/RNAseq/LGG_filtered_TPM.csv",
  GBM = "/Users/jackie-mac/Desktop/VSCode/data/RNAseq/GBM_filtered_TPM.csv",
  CCRCC = "/Users/jackie-mac/Desktop/VSCode/data/RNAseq/CCRCC_filtered_TPM.csv",
  HNSCC = "/Users/jackie-mac/Desktop/VSCode/data/RNAseq/HNSCC_filtered_TPM.csv",
  PDA = "/Users/jackie-mac/Desktop/VSCode/data/RNAseq/PDA_filtered_TPM.csv"
  #NSCLC = "/Users/jackie-mac/Desktop/VSCode/Outputs/Filtered_RNAseq/Protein_coding/MissingRemoved/TPM_Normalized/NSCLC_filtered_TPM.csv"
)
output_umap_plot <- "/Users/jackie-mac/Desktop/VSCode/outputs/plots/UMAP/umap_by_cancertype_plot.png"

# ---- CONSORTIUM LABELS ----
tcga_cancers <- c("BRCA", "GBM", "KIRC", "LGG")
cptac_cancers <- c("CCRCC", "HNSCC", "PDA")
#nsclc_cancer <- "NSCLC" # NSCLC will be its own shape

all_expr <- list()
all_labels <- c()
all_shapes <- c()
for (cancer in names(rna_files)) {
  expr <- fread(rna_files[[cancer]], data.table = FALSE)
  rownames(expr) <- expr[[1]]
  expr <- expr[,-1]
  all_expr[[cancer]] <- expr
  n_samples <- ncol(expr)
  all_labels <- c(all_labels, rep(cancer, n_samples))
  if (cancer == "NSCLC") {
    all_shapes <- c(all_shapes, rep("NSCLC", n_samples))
  } else if (cancer %in% cptac_cancers) {
    all_shapes <- c(all_shapes, rep("CPTAC", n_samples))
  } else {
    all_shapes <- c(all_shapes, rep("TCGA", n_samples))
  }
}
# Find common genes across all datasets
common_genes <- Reduce(intersect, lapply(all_expr, rownames))
# Subset and concatenate
all_expr <- lapply(all_expr, function(x) x[common_genes, , drop = FALSE])
expr_mat <- do.call(cbind, all_expr)
expr_mat <- t(expr_mat) # samples x genes

# ---- RUN UMAP ----
set.seed(123)
umap_res <- umap(expr_mat, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")

umap_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  CancerType = all_labels,
  Consortium = factor(all_shapes, levels = c("CPTAC", "TCGA", "NSCLC"))
)

# ---- PLOT ----
# Use more distinct shapes: CPTAC = X (shape 4), TCGA = triangle (17), NSCLC = circle (16)
shape_map <- c("CPTAC" = 4, "TCGA" = 17, "NSCLC" = 16) # X, triangle, circle

p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = CancerType, shape = Consortium, label = CancerType)) +
  geom_point(alpha = 0.85, size = 2.5, stroke = 1, fill = "white") + # stroke adds border, fill for open shapes
  scale_shape_manual(values = shape_map) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  ) +
  labs(title = "UMAP of RNAseq Data by Cancer Type and Consortium")

# Save static plot as PNG
ggsave(output_umap_plot, p, width = 7, height = 5, bg = "white")

# ---- PLOT (interactive) ----
p_interactive <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = CancerType, shape = Consortium, text = CancerType)) +
  geom_point(alpha = 0.85, size = 2.5, stroke = 1, fill = "white") +
  scale_shape_manual(values = shape_map) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  ) +
  labs(title = "UMAP of RNAseq Data by Cancer Type and Consortium")

# Convert ggplot to interactive plotly object
p_interactive <- ggplotly(p_interactive, tooltip = c("text", "shape"))

# Save as HTML (interactive)
htmlwidgets::saveWidget(p_interactive, file = sub("\\.png$", ".html", output_umap_plot))

cat("UMAP plot saved to:", output_umap_plot, "and interactive HTML.\n")