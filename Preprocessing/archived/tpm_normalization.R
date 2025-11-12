#!/usr/bin/env Rscript

# TPM normalization for RNA-seq raw counts using gene lengths retrieved via biomaRt.
# - Accepts CSV/TSV input with either genes x samples or samples x genes orientation.
# - Detects orientation and ensures genes are rows before normalization.
# - Genes are expected to be HGNC symbols; gene lengths are derived from Ensembl via biomaRt
#   by aggregating transcript lengths per gene (median by default).
# - Outputs a CSV with GeneSymbol as the first column and samples as subsequent columns.

# NOTE: This script queries Ensembl via biomaRt and therefore requires outbound Internet access.
#       Do not run on H4H or other remote clusters without Internet connectivity. Run locally
#       on your laptop/desktop, or supply a cached gene-lengths file via --cache-lengths to
#       avoid live queries on subsequent runs.

# Load lightweight CLI dependency quietly; auto-install locally if missing.
if (!requireNamespace("optparse", quietly = TRUE)) {
	if (interactive()) {
		message("Installing missing package: optparse ...")
		install.packages("optparse")
	} else {
		stop("Package 'optparse' is required. Please install it (install.packages('optparse')).")
	}
}
suppressPackageStartupMessages(library(optparse, quietly = TRUE, warn.conflicts = FALSE))

# Ensure a required package is installed and loadable.
# Args:
#   pkg: Character scalar, the package name to load.
# Raises:
#   error if the package is not installed or cannot be loaded.
safe_require <- function(pkg, auto_install = FALSE) {
	# Install CRAN or Bioconductor packages if missing when auto_install=TRUE or in interactive mode.
	if (!requireNamespace(pkg, quietly = TRUE)) {
		if (auto_install || interactive()) {
			if (pkg %in% c("biomaRt")) {
				if (!requireNamespace("BiocManager", quietly = TRUE)) {
					message("Installing missing package: BiocManager ...")
					install.packages("BiocManager")
				}
				message(sprintf("Installing missing Bioconductor package: %s ...", pkg))
				BiocManager::install(pkg, ask = FALSE, update = FALSE)
			} else {
				message(sprintf("Installing missing package: %s ...", pkg))
				install.packages(pkg)
			}
		} else {
			stop(sprintf("Package '%s' is required. Please install it first (e.g., install.packages('%s') or BiocManager::install('%s')).", pkg, pkg, pkg))
		}
	}
	suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE))
}

# Guess the delimiter for a delimited text file by peeking the first line.
# Args:
#   path: Character path to a CSV/TSV file.
# Returns:
#   A single character: "\t" if tab is detected, otherwise ",".
guess_delim <- function(path) {
	con <- file(path, "r"); on.exit(close(con), add = TRUE)
	first <- readLines(con, n = 1)
	if (grepl("\t", first)) "\t" else ","
}

# Read a delimited counts matrix preserving column names and string values.
# Args:
#   path: Character path to CSV/TSV file.
# Returns:
#   A list with elements:
#     df  - data.frame of the parsed file
#     sep - the delimiter detected ("\t" or ",")
read_counts <- function(path) {
	sep <- guess_delim(path)
	df <- read.table(path, sep = sep, header = TRUE, check.names = FALSE, quote = "\"", comment.char = "", stringsAsFactors = FALSE)
	list(df = df, sep = sep)
}

# Heuristic to determine if strings resemble gene symbols.
# Args:
#   x: character vector of tokens to test.
# Returns:
#   logical vector, TRUE where a token contains at least one letter and isn't blank.
looks_like_gene_symbol <- function(x) {
	# Heuristic: contains letters and not purely numeric; allow '-', '.', numbers.
	grepl("[A-Za-z]", x) & !grepl("^\\s*$", x)
}

# Detect whether the input table is genes x samples (genes in rows) or samples x genes.
# Args:
#   df: data.frame of the input table with an identifier column first.
# Returns:
#   character scalar: "genes_rows" if genes are rows, "genes_cols" if genes are columns.
detect_orientation <- function(df) {
	# Returns "genes_rows" or "genes_cols"; also indicates the name of the identifier header if present.
	cols <- colnames(df)

	# Heuristics: if many column names (excluding the first) look like gene symbols, it's samples x genes
	gene_like_cols <- mean(looks_like_gene_symbol(cols[-1]))

	# For the first few rows, check if entries across columns (excluding first) are largely numeric
	n_check <- min(50, nrow(df))
	part <- df[seq_len(n_check), -1, drop = FALSE]
	numeric_frac_rows <- mean(
		suppressWarnings(
			sapply(part, function(v) mean(!is.na(as.numeric(v)), na.rm = TRUE))
		),
		na.rm = TRUE
	)

	# Check first column values for gene-like strings
	first_col_gene_like <- mean(looks_like_gene_symbol(head(df[[1]], 50)))

	# Simple decisions:
	# - If many gene-like column names, assume samples x genes
	# - Else if first column is gene-like and the rest are numeric -> genes x samples
	if (!is.nan(gene_like_cols) && gene_like_cols > 0.6) {
		return("genes_cols")
	}
	if (!is.nan(first_col_gene_like) && first_col_gene_like > 0.6 && !is.nan(numeric_frac_rows) && numeric_frac_rows > 0.6) {
		return("genes_rows")
	}
	# Fallback: if there are many columns (suggesting samples) treat as genes_rows
	if (ncol(df) >= nrow(df)) return("genes_rows")
	"genes_cols"
}

# Ensure the counts table is oriented as genes (rows) x samples (columns).
# If the table is samples x genes, transpose it.
# Args:
#   df: data.frame with first column being either GeneSymbol (genes x samples) or SampleID (samples x genes).
# Returns:
#   list with elements:
#     counts      - data.frame of numeric-like values with genes as rows and samples as columns
#     genes       - character vector of gene symbols (row order of counts)
#     orientation - "genes_rows" or "genes_cols" indicating original orientation
ensure_genes_by_samples <- function(df) {
	ori <- detect_orientation(df)
	if (ori == "genes_rows") {
		genes <- df[[1]]
		mat <- as.data.frame(df[, -1, drop = FALSE], check.names = FALSE)
		rownames(mat) <- make.names(genes, unique = TRUE)
		list(counts = mat, genes = as.character(genes), orientation = ori)
	} else {
		# samples x genes: first column is sample ID, other columns are genes
		samples <- df[[1]]
		mat <- as.data.frame(df[, -1, drop = FALSE], check.names = FALSE)
		rownames(mat) <- as.character(samples)
		# transpose to genes x samples
		mat_t <- as.data.frame(t(mat), stringsAsFactors = FALSE, check.names = FALSE)
		genes <- rownames(mat_t)
		list(counts = mat_t, genes = as.character(genes), orientation = ori)
	}
}

# Retrieve gene lengths for HGNC symbols from Ensembl via biomaRt.
# Args:
#   genes            - character vector of HGNC symbols
#   species_dataset  - Ensembl dataset (default: "hsapiens_gene_ensembl")
#   host             - optional Ensembl host URL (e.g., archive URL) for reproducibility
#   method           - aggregation of transcript lengths per gene: one of "median", "mean", "max"
# Returns:
#   data.frame with columns: hgnc_symbol, gene_length (in base pairs)
# Notes:
#   Genes without transcript_length entries are omitted.
fetch_gene_lengths_biomart <- function(genes, species_dataset = "hsapiens_gene_ensembl", host = NULL, method = c("median", "mean", "max"), auto_install = FALSE) {
	method <- match.arg(method)
	safe_require("biomaRt", auto_install = auto_install)
	if (is.null(host) || nchar(host) == 0) {
		mart <- biomaRt::useEnsembl(biomart = "genes", dataset = species_dataset)
	} else {
		mart <- biomaRt::useEnsembl(biomart = "genes", dataset = species_dataset, host = host)
	}
	attrs <- c("hgnc_symbol", "ensembl_gene_id", "transcript_length")
	res <- biomaRt::getBM(attributes = attrs, filters = "hgnc_symbol", values = unique(genes), mart = mart)
	res <- res[!is.na(res$hgnc_symbol) & nzchar(res$hgnc_symbol), , drop = FALSE]
	res <- res[!is.na(res$transcript_length), , drop = FALSE]
	if (nrow(res) == 0) {
		warning("biomaRt returned no lengths for provided genes.")
		return(data.frame(hgnc_symbol = character(0), gene_length = numeric(0), stringsAsFactors = FALSE))
	}
	agg_fun <- switch(method,
										median = function(x) stats::median(x, na.rm = TRUE),
										mean = function(x) mean(x, na.rm = TRUE),
										max = function(x) max(x, na.rm = TRUE))
	gene_len <- stats::aggregate(transcript_length ~ hgnc_symbol, data = res, FUN = agg_fun)
	colnames(gene_len) <- c("hgnc_symbol", "gene_length")
	gene_len
}

# Convert a data.frame of values to a numeric matrix, coercing non-finite to 0.
# Args:
#   df: data.frame of values to coerce.
# Returns:
#   numeric matrix with NA/Inf converted to 0, suitable for arithmetic.
to_numeric_matrix <- function(df) {
	m <- as.matrix(df)
	storage.mode(m) <- "double"
	m[!is.finite(m)] <- NA_real_
	m[is.na(m)] <- 0
	m
}

# Compute TPM from raw counts and gene lengths.
# Args:
#   counts_df      - data.frame with genes as rows, samples as columns (raw counts)
#   gene_lengths_bp- named numeric vector of gene lengths in base pairs; names must match rownames(counts_df)
# Returns:
#   data.frame of TPM values (genes x samples), rownames preserved from counts_df.
compute_tpm <- function(counts_df, gene_lengths_bp) {
	# counts_df: data.frame genes x samples
	# gene_lengths_bp: named numeric vector (names are gene symbols)
	mat <- to_numeric_matrix(counts_df)
	genes <- rownames(counts_df)
	# Ensure lengths align to matrix rows
	lens <- gene_lengths_bp[genes]
	keep <- !is.na(lens) & lens > 0
	dropped <- sum(!keep)
	if (dropped > 0) {
		message(sprintf("Dropping %d genes without valid length from TPM computation.", dropped))
	}
	mat <- mat[keep, , drop = FALSE]
	lens <- lens[keep]
	# RPK = counts / (length_kb)
	length_kb <- lens / 1000
	rpk <- sweep(mat, 1, length_kb, FUN = "/")
	# Per-sample scaling to 1e6
	scale <- colSums(rpk, na.rm = TRUE)
	# Avoid divide by zero
	scale[scale == 0] <- NA_real_
	tpm <- sweep(rpk, 2, scale, FUN = "/") * 1e6
	tpm[!is.finite(tpm)] <- 0
	as.data.frame(tpm, check.names = FALSE)
}

main <- function() {
	option_list <- list(
		make_option(c("-i", "--input"), type = "character", help = "Path to raw counts matrix (CSV/TSV)."),
		make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output CSV path for TPM matrix."),
		make_option(c("--species-dataset"), type = "character", default = "hsapiens_gene_ensembl", help = "Ensembl dataset (e.g., hsapiens_gene_ensembl)."),
		make_option(c("--ensembl-host"), type = "character", default = NULL, help = "Optional Ensembl host (e.g., https://dec2021.archive.ensembl.org)."),
		make_option(c("--length-method"), type = "character", default = "median", help = "How to aggregate transcript lengths per gene: median|mean|max [default: %default]."),
		make_option(c("--auto-install"), action = "store_true", default = FALSE, help = "Auto-install missing R packages when running locally/interactive."),
			make_option(c("--cache-lengths"), type = "character", default = NULL, help = "Optional CSV path to read/write cached gene lengths (columns: hgnc_symbol,gene_length)."),
			make_option(c("--missing-lengths-out"), type = "character", default = NULL, help = "Optional CSV to write genes missing lengths from biomaRt (column: hgnc_symbol).")
	)
	parser <- OptionParser(option_list = option_list, description = "TPM-normalize RNA-seq counts using gene lengths via biomaRt. Accepts genes x samples or samples x genes and re-orients as needed.")
	args <- parse_args(parser, positional_arguments = FALSE)

	if (is.null(args$input) || !nzchar(args$input)) {
		if (interactive()) {
			message("Select an input counts file (CSV/TSV)...")
			args$input <- file.choose()
		} else {
			print_help(parser)
			stop("--input is required")
		}
	}

	io <- read_counts(args$input)
	df <- io$df

	oriented <- ensure_genes_by_samples(df)
	counts <- oriented$counts
	genes <- oriented$genes

	message(sprintf("Input orientation: %s; genes: %d; samples: %d", oriented$orientation, nrow(counts), ncol(counts)))

		# Fetch gene lengths (use cache if available)
		gene_len_df <- NULL
		cache_path <- args$`cache-lengths`
		if (!is.null(cache_path) && nzchar(cache_path) && file.exists(cache_path)) {
			message(sprintf("Reading cached gene lengths: %s", cache_path))
			gene_len_df <- tryCatch({
				utils::read.csv(cache_path, stringsAsFactors = FALSE, check.names = FALSE)
			}, error = function(e) {
				message("Failed to read cache, will query biomaRt. Error: ", conditionMessage(e))
				NULL
			})
			if (!is.null(gene_len_df) && !all(c("hgnc_symbol", "gene_length") %in% colnames(gene_len_df))) {
				message("Cache file missing required columns (hgnc_symbol,gene_length). Ignoring cache.")
				gene_len_df <- NULL
			}
		}
			if (is.null(gene_len_df)) {
			gene_len_df <- fetch_gene_lengths_biomart(
				genes,
				species_dataset = args$`species-dataset`,
				host = args$`ensembl-host`,
				method = args$`length-method`,
				auto_install = isTRUE(args$`auto-install`)
			)
			if (!is.null(cache_path) && nzchar(cache_path)) {
				dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
				utils::write.csv(gene_len_df, file = cache_path, row.names = FALSE)
				message(sprintf("Wrote gene lengths cache: %s", cache_path))
			}
		}

			# Track genes without retrieved lengths
			got_syms <- unique(as.character(gene_len_df$hgnc_symbol))
			need_syms <- unique(as.character(genes))
			missing_syms <- setdiff(need_syms, got_syms)
			if (length(missing_syms) > 0) {
				message(sprintf("Note: %d genes had no length from biomaRt and will be dropped in TPM.", length(missing_syms)))
				if (!is.null(args$`missing-lengths-out`) && nzchar(args$`missing-lengths-out`)) {
					out_miss <- data.frame(hgnc_symbol = missing_syms, stringsAsFactors = FALSE)
					dir.create(dirname(args$`missing-lengths-out`), recursive = TRUE, showWarnings = FALSE)
					utils::write.csv(out_miss, file = args$`missing-lengths-out`, row.names = FALSE)
					message(sprintf("Wrote missing-length genes: %s", args$`missing-lengths-out`))
				}
			}
	if (nrow(gene_len_df) == 0) {
		stop("No gene lengths retrieved from biomaRt; cannot compute TPM.")
	}
	# Build named vector for lengths
	# Map by original gene symbols; biomaRt returns hgnc_symbol as provided
	gene_lengths <- setNames(gene_len_df$gene_length, gene_len_df$hgnc_symbol)

	# TPM computation
	tpm_df <- compute_tpm(counts, gene_lengths)

	# Prepare output: GeneSymbol + sample columns
	out <- cbind(data.frame(GeneSymbol = rownames(tpm_df), check.names = FALSE), tpm_df)

	out_path <- args$output
	if (is.null(out_path) || !nzchar(out_path)) {
		stem <- sub("\\.[^.]*$", "", basename(args$input))
		out_path <- file.path(dirname(args$input), paste0(stem, "_TPM.csv"))
	}
	dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
	utils::write.csv(out, file = out_path, row.names = FALSE)
	message(sprintf("Wrote TPM matrix: %s", out_path))
}

if (identical(environment(), globalenv())) {
	main()
}

