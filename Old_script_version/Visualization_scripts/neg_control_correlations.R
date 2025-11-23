#!/usr/bin/env Rscript

# Pearson correlations between original and negative-control FMCIB features.
# Scenarios:
#   - Regions: "full", "roi"
#   - Controls (vs original): "randomized", "sampled", "shuffled"
# Inputs:
#   - Feature matrix (CSV) with columns:
#       SampleID, readii_Permutation, readii_Region, and 4096 feature columns named like pred_0..pred_4095
# Outputs:
#   - Per-region correlation CSV: features (rows) x controls (cols)
#   - Per-region heatmap PNG with colors: purple=-1, green=0, red=1

# ---- Helpers: package loading ----
utils::globalVariables(c("Control", "Feature", "Correlation", "x", "y", "z"))
auto_require <- function(pkg, auto_install = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (auto_install || interactive()) install.packages(pkg)
    else stop(sprintf("Missing package '%s'. Install it with install.packages('%s').", pkg, pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ---- Correlation core ----
find_feature_cols <- function(df, prefix = "pred_") {
  # strict: pred_<digits>, case-insensitive, ignore leading/trailing spaces in names
  nm <- trimws(colnames(df))
  sel <- grepl(paste0("(?i)^", gsub("([\\^$.|?*+(){}])","\\\\\\1", prefix), "[0-9]+$"), nm, perl = TRUE)
  nm[sel]
}

order_features_numeric <- function(features, prefix = "pred_") {
  idx <- suppressWarnings(as.integer(sub(sprintf("(?i)^%s", gsub("([\\^$.|?*+(){}])","\\\\\\1", prefix)), "", features, perl = TRUE)))
  o <- order(idx, na.last = TRUE)
  features[o]
}

compute_cor_for_pair <- function(df_region, id_col, perm_col, control_name, feature_cols) {
	# Split original vs control within region and merge by SampleID
	orig <- df_region[df_region[[perm_col]] == "original", c(id_col, feature_cols), drop = FALSE]
	ctrl <- df_region[df_region[[perm_col]] == control_name, c(id_col, feature_cols), drop = FALSE]
	if (nrow(orig) == 0 || nrow(ctrl) == 0) {
		# Diagnostics: show unique values and counts to aid debugging
		perm_vals <- tryCatch(unique(df_region[[perm_col]]), error = function(e) character())
		perm_tab <- tryCatch(utils::capture.output(print(table(df_region[[perm_col]], useNA = "ifany"))), error = function(e) NULL)
		cols <- tryCatch(colnames(df_region), error = function(e) NULL)
		message(sprintf("compute_cor_for_pair: perm_col='%s', requested control='%s'", perm_col, control_name))
		if (!is.null(cols)) message("df_region columns: ", paste(head(cols, 20), collapse = ", "), if (length(cols) > 20) " ..." else "")
		if (length(perm_vals)) message("unique perms: ", paste(head(perm_vals, 10), collapse = ", "), if (length(perm_vals) > 10) " ..." else "")
		if (!is.null(perm_tab)) message(paste(perm_tab, collapse = "\n"))
		warning(sprintf("No rows for '%s' or 'original' in region subset; correlations set to NA.", control_name))
		return(rep(NA_real_, length(feature_cols)))
	}
	merged <- merge(orig, ctrl, by = id_col, suffixes = c(".orig", ".ctrl"))
	if (nrow(merged) == 0) {
		warning(sprintf("No overlapping %s between original and %s; correlations set to NA.", id_col, control_name))
		return(rep(NA_real_, length(feature_cols)))
	}
	cors <- vapply(seq_along(feature_cols), function(i) {
		f <- feature_cols[i]
		x <- suppressWarnings(as.numeric(merged[[paste0(f, ".orig")]]))
		y <- suppressWarnings(as.numeric(merged[[paste0(f, ".ctrl")]]))
		if (all(!is.finite(x)) || all(!is.finite(y))) return(NA_real_)
		suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs", method = "pearson"))
	}, numeric(1))
	cors
}

# Feature-wise correlations (vector) for reference; original vs same feature in control.
correlation_matrix_by_region <- function(df, region_name, id_col, perm_col, region_col, controls, feature_cols) {
	# Filter region (case-insensitive)
	region_lc <- trimws(tolower(region_name))
	df$.__perm__ <- trimws(tolower(df[[perm_col]]))
	df$.__region__ <- trimws(tolower(df[[region_col]]))
	sub <- df[df$.__region__ == region_lc, , drop = FALSE]
	if (nrow(sub) == 0) {
		warning(sprintf("No rows for region '%s'", region_name))
		mat <- matrix(NA_real_, nrow = length(feature_cols), ncol = length(controls))
		dimnames(mat) <- list(feature_cols, controls)
		return(as.data.frame(mat, check.names = FALSE))
	}
	# Ensure features are numeric-like
	for (f in feature_cols) sub[[f]] <- suppressWarnings(as.numeric(sub[[f]]))
	# Compute correlations for each control vs original
	mat <- sapply(controls, function(ctrl_name) {
		compute_cor_for_pair(sub, id_col, ".__perm__", tolower(ctrl_name), feature_cols)
	})
	mat <- as.data.frame(mat, check.names = FALSE)
	rownames(mat) <- feature_cols
	mat
}

compute_cor_square_for_pair <- function(df_region, id_col, perm_col, control_name, feature_cols) {
	orig <- df_region[df_region[[perm_col]] == "original", c(id_col, feature_cols), drop = FALSE]
	ctrl <- df_region[df_region[[perm_col]] == control_name, c(id_col, feature_cols), drop = FALSE]
	if (nrow(orig) == 0 || nrow(ctrl) == 0) {
		warning(sprintf("No rows for '%s' or 'original' in region subset; returning NA matrix.", control_name))
		mat_na <- matrix(NA_real_, nrow = length(feature_cols), ncol = length(feature_cols))
		dimnames(mat_na) <- list(feature_cols, feature_cols)
		return(mat_na)
	}
	merged <- merge(orig, ctrl, by = id_col, suffixes = c(".orig", ".ctrl"))
	if (nrow(merged) == 0) {
		warning(sprintf("No overlapping %s between original and %s; returning NA matrix.", id_col, control_name))
		mat_na <- matrix(NA_real_, nrow = length(feature_cols), ncol = length(feature_cols))
		dimnames(mat_na) <- list(feature_cols, feature_cols)
		return(mat_na)
	}
	# Build numeric matrices: samples x features
	X <- sapply(feature_cols, function(f) suppressWarnings(as.numeric(merged[[paste0(f, ".orig")]])))
	Y <- sapply(feature_cols, function(f) suppressWarnings(as.numeric(merged[[paste0(f, ".ctrl")]])))
	# Ensure matrix type when a single feature is present
	X <- as.matrix(X); Y <- as.matrix(Y)
	colnames(X) <- feature_cols; colnames(Y) <- feature_cols
	C <- suppressWarnings(stats::cor(X, Y, use = "pairwise.complete.obs", method = "pearson"))
	rownames(C) <- feature_cols; colnames(C) <- feature_cols
	C
}

# Self feature x feature correlation for a single permutation within a region
compute_cor_square_self <- function(df_region, id_col, perm_col, perm_name, feature_cols) {
	subp <- df_region[df_region[[perm_col]] == perm_name, c(id_col, feature_cols), drop = FALSE]
	if (nrow(subp) == 0) {
		warning(sprintf("No rows for permutation '%s' in region subset; returning NA matrix.", perm_name))
		mat_na <- matrix(NA_real_, nrow = length(feature_cols), ncol = length(feature_cols))
		dimnames(mat_na) <- list(feature_cols, feature_cols)
		return(mat_na)
	}
	X <- sapply(feature_cols, function(f) suppressWarnings(as.numeric(subp[[f]])))
	X <- as.matrix(X)
	colnames(X) <- feature_cols
	C <- suppressWarnings(stats::cor(X, X, use = "pairwise.complete.obs", method = "pearson"))
	rownames(C) <- feature_cols; colnames(C) <- feature_cols
	C
}

# General pairwise feature x feature correlation between two permutations within a region
compute_cor_square_between <- function(df_region, id_col, perm_col, perm_a, perm_b, feature_cols) {
	A <- df_region[df_region[[perm_col]] == perm_a, c(id_col, feature_cols), drop = FALSE]
	B <- df_region[df_region[[perm_col]] == perm_b, c(id_col, feature_cols), drop = FALSE]
	if (nrow(A) == 0 || nrow(B) == 0) {
		warning(sprintf("No rows for '%s' or '%s' in region subset; returning NA matrix.", perm_a, perm_b))
		mat_na <- matrix(NA_real_, nrow = length(feature_cols), ncol = length(feature_cols))
		dimnames(mat_na) <- list(feature_cols, feature_cols)
		return(mat_na)
	}
	merged <- merge(A, B, by = id_col, suffixes = c(paste0(".", perm_a), paste0(".", perm_b)))
	if (nrow(merged) == 0) {
		warning(sprintf("No overlapping %s between %s and %s; returning NA matrix.", id_col, perm_a, perm_b))
		mat_na <- matrix(NA_real_, nrow = length(feature_cols), ncol = length(feature_cols))
		dimnames(mat_na) <- list(feature_cols, feature_cols)
		return(mat_na)
	}
	X <- sapply(feature_cols, function(f) suppressWarnings(as.numeric(merged[[paste0(f, ".", perm_a)]])))
	Y <- sapply(feature_cols, function(f) suppressWarnings(as.numeric(merged[[paste0(f, ".", perm_b)]])))
	X <- as.matrix(X); Y <- as.matrix(Y)
	colnames(X) <- feature_cols; colnames(Y) <- feature_cols
	C <- suppressWarnings(stats::cor(X, Y, use = "pairwise.complete.obs", method = "pearson"))
	rownames(C) <- feature_cols; colnames(C) <- feature_cols
	C
}

# ---- Plotting ----
save_heatmap <- function(cor_df, region_name, out_png, width = 14, height_per_1000 = 14) {
  auto_require("ggplot2")
  auto_require("tidyr")
  auto_require("tibble")

  tmp <- tibble::rownames_to_column(cor_df, var = "Feature")
  long <- tidyr::pivot_longer(tmp, cols = -"Feature", names_to = "Control", values_to = "Correlation")

  n_feat <- length(unique(long$Feature))
  height <- max(4, ceiling(n_feat / 1000) * height_per_1000)

	# Smooth rainbow-like palette with more stops, centered at green for 0
	palette_stops <- c(
		"#7300c0",  # -1 deep purple
		"#2f5391",
		"#2580a1",
		"#2a9e96",
		"#00a6a6",
		"#00be00",  #  0 green center
		"#3da03a",
		"#aad367",
		"#ffe600",
		"#ff8800",
		"#ff001e"   # +1 red
	)

		p <- ggplot2::ggplot(long, ggplot2::aes_string(x = "Control", y = "Feature", fill = "Correlation")) +
		ggplot2::geom_tile() +
		ggplot2::scale_fill_gradientn(
			colors = palette_stops,
			limits = c(-1, 1),
			oob = scales::squish,
				breaks = seq(-1, 1, by = 0.25),
				guide = ggplot2::guide_colorbar(
					frame.colour = "black",
					ticks.colour = "black",
					barwidth = grid::unit(0.2, "cm"),
					barheight = grid::unit(8, "cm")
				)
		) +
    ggplot2::labs(
      title = sprintf("Original vs %s (%s)", ifelse(ncol(cor_df) == 1, colnames(cor_df)[1], "Control"), region_name),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  ggplot2::ggsave(out_png, p, width = width, height = height, dpi = 150, limitsize = FALSE, bg = "white")
}

# Save a large square heatmap (features x features) using base graphics (efficient for 4096x4096).
# Colors: purple (-1), green (0), red (1). White background enforced.
save_heatmap_square <- function(cor_mat, title, out_png, width_px = 2400, height_px = 2400, x_axis_name = NULL, y_axis_name = NULL, lower_triangle_only = FALSE, triangle_corner = "bottom-left") {
	# Color palette and device: smooth rainbow-like palette centered at green for 0
	palette_stops <- c(
		"#5e4fa2",  # -1 deep purple
		"#3b66b1",
		"#2c99c2",
		"#35c9c0",
		"#00a6a6",
		"#00a600",  #  0 green center
		"#66bd63",
		"#a6d96a",
		"#fee08b",
		"#f46d43",
		"#d53e4f"   # +1 red
	)
	pal_func <- grDevices::colorRampPalette(palette_stops, interpolate = "linear")
	pal <- pal_func(1024)
	grDevices::png(filename = out_png, width = width_px, height = height_px, res = 150, bg = "white")
	op <- par(no.readonly = TRUE); on.exit({par(op); grDevices::dev.off()}, add = TRUE)
	# Split device: left = heatmap, right = thinner colorbar
	graphics::layout(matrix(c(1,2), nrow = 1), widths = c(0.96, 0.04))
	on.exit(graphics::layout(1), add = TRUE)
	par(bg = "white", mar = c(3.5, 4, 3, 1))

	# Clamp values to [-1,1] and build plot matrix in display orientation
	m <- cor_mat
	m[!is.finite(m)] <- NA_real_
	m[m > 1] <- 1; m[m < -1] <- -1
	# Transform to display orientation (so masking aligns with what's drawn)
	plot_m <- t(m)[, nrow(m):1, drop = FALSE]
	# Optionally mask to show only one triangle in display space
	if (isTRUE(lower_triangle_only) && is.matrix(plot_m) && nrow(plot_m) == ncol(plot_m)) {
		pr <- nrow(plot_m); pc <- ncol(plot_m)
		ridx <- row(plot_m); cidx <- col(plot_m)
		# In display orientation, the main diagonal of the original matrix maps to the anti-diagonal (r + c = pr + 1).
		# Keep one half relative to this anti-diagonal based on requested corner.
		keep <- switch(trimws(tolower(triangle_corner)),
			"bottom-right" = (ridx + cidx >= pr + 1),
			# default: bottom-left
			(ridx + cidx <= pr + 1)
		)
		plot_m[!keep] <- NA_real_
	}
	pr <- nrow(plot_m); pc <- ncol(plot_m)
	graphics::image(x = 1:pc, y = 1:pr, z = plot_m, col = pal, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "")
	title(main = title, line = 1)
	if (!is.null(x_axis_name)) graphics::mtext(x_axis_name, side = 1, line = 2)
	if (!is.null(y_axis_name)) graphics::mtext(y_axis_name, side = 2, line = 2)

	# Draw vertical colorbar matching heatmap height on the right panel (very thin, so keep margins small)
	par(bg = "white", mar = c(3, 0.2, 2, 2))
	zseq <- seq(-1, 1, length.out = 256)
	# Build a single-column gradient: use boundary semantics for x (length 2), 1 row in z
	xbar <- c(0, 1)
	legend_mat <- matrix(zseq, nrow = 1)
	graphics::image(x = xbar, y = zseq, z = legend_mat, col = pal, xlab = "", ylab = "", axes = FALSE)
	legend_breaks <- seq(-1, 1, by = 0.25)
	graphics::axis(side = 4, at = legend_breaks, labels = format(legend_breaks, trim = TRUE), las = 1, cex.axis = 0.7)
	graphics::box()
}

# ---- CLI ----
auto_require("optparse")

option_list <- list(
	optparse::make_option(c("-i", "--input"), type = "character", help = "Path to feature matrix (CSV/TSV)."),
	optparse::make_option(c("-o", "--outdir"), type = "character", default = "data/procdata/radiomic_features/correlations", help = "Output directory for CSVs/PNGs [default: %default]"),
	optparse::make_option(c("--id-col"), type = "character", default = "SampleID", help = "Sample ID column name [default: %default]"),
	optparse::make_option(c("--perm-col"), type = "character", default = "readii_Permutation", help = "Permutation column name [default: %default]"),
	optparse::make_option(c("--region-col"), type = "character", default = "readii_Region", help = "Region column name [default: %default]"),
	optparse::make_option(c("--regions"), type = "character", default = "full,roi,non_roi", help = "Comma-separated list of regions [default: %default]"),
	optparse::make_option(c("--controls"), type = "character", default = "randomized,sampled,shuffled", help = "Comma-separated controls compared to 'original' [default: %default]"),
	optparse::make_option(c("--feature-prefix"), type = "character", default = "pred_", help = "Prefix for feature columns [default: %default]"),
	 optparse::make_option(c("--only-original-csv"), action = "store_true", default = FALSE, help = "Write only the original self-correlation matrix as CSV; do not generate any PNGs [default: %default]"),
	 optparse::make_option(c("--full-matrix"), action = "store_true", default = FALSE, help = "Render full self-correlation matrix (disable lower-triangle-only view) [default: %default]"),
	 optparse::make_option(c("--triangle-corner"), type = "character", default = "bottom-left", help = "Triangle corner to keep when masking self-correlation heatmaps: 'bottom-left' or 'bottom-right' [default: %default]")
)

	parser <- optparse::OptionParser(
		option_list = option_list,
		description = "Compute Pearson correlations between original and negative control features per region and plot heatmaps (input: CSV)."
	)
args <- optparse::parse_args(parser, positional_arguments = FALSE)

only_original_csv <- isTRUE(args$`only-original-csv`)
lower_triangle_only <- !isTRUE(args$`full-matrix`)
triangle_corner <- trimws(tolower(as.character(args$`triangle-corner`)))
if (!triangle_corner %in% c("bottom-left", "bottom-right")) triangle_corner <- "bottom-left"

if (is.null(args$input) || !nzchar(args$input)) {
	if (interactive()) {
		message("Select a feature matrix file (CSV/TSV)...")
		args$input <- file.choose()
	} else {
		optparse::print_help(parser)
		stop("--input is required")
	}
}

# Read CSV
if (!grepl("\\.csv$", tolower(args$input))) {
  warning("Input does not end with .csv; attempting to read as CSV anyway.")
}
df <- utils::read.csv(
  args$input,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
# Normalize column names (strip spaces)
colnames(df) <- trimws(colnames(df))

if (!all(c(args$`id-col`, args$`perm-col`, args$`region-col`) %in% colnames(df))) {
  stop("Input is missing required columns: ", paste(c(args$`id-col`, args$`perm-col`, args$`region-col`), collapse = ", "))
}

feat_cols <- find_feature_cols(df, prefix = args$`feature-prefix`)
if (length(feat_cols) == 0) {
  stop("No feature columns found matching strict pattern '^", args$`feature-prefix`, "[0-9]+$'.")
}

# Optional: drop non-numeric or constant features (defensive)
non_num <- vapply(feat_cols, function(f) !is.numeric(type.convert(df[[f]], as.is = TRUE)), logical(1))
if (any(non_num)) {
  message("Dropping non-numeric feature columns: ", paste(head(feat_cols[non_num], 10), collapse = ", "),
          if (sum(non_num) > 10) " ..." else "")
}
feat_cols <- feat_cols[!non_num]
const_feat <- vapply(feat_cols, function(f) {
  x <- suppressWarnings(as.numeric(df[[f]]))
  length(unique(x[is.finite(x)])) <= 1
}, logical(1))
if (any(const_feat)) {
  message("Dropping constant feature columns (no variance in dataset subset): ",
          paste(head(feat_cols[const_feat], 10), collapse = ", "),
          if (sum(const_feat) > 10) " ..." else "")
}
feat_cols <- feat_cols[!const_feat]

if (length(feat_cols) == 0) stop("After filtering, no usable feature columns remain.")

message(sprintf("Using %d feature columns (prefix '%s'): e.g., %s",
                length(feat_cols), args$`feature-prefix`,
                paste(head(order_features_numeric(feat_cols, args$`feature-prefix`), 5), collapse = ", ")))

# Normalize labels: trim whitespace and lowercase for robust matching
df[[args$`perm-col`]] <- trimws(tolower(as.character(df[[args$`perm-col`]])))
df[[args$`region-col`]] <- trimws(tolower(as.character(df[[args$`region-col`]])))

regions <- strsplit(args$regions, ",", fixed = TRUE)[[1]]
regions <- trimws(tolower(regions))
controls <- strsplit(args$controls, ",", fixed = TRUE)[[1]]
controls <- trimws(tolower(controls))

outdir <- args$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Order features by numeric index if possible
feat_ordered <- order_features_numeric(feat_cols, prefix = args$`feature-prefix`)

## Preflight: report available permutations by region
perm_tab <- tryCatch({
	as.data.frame.matrix(table(df[[args$`perm-col`]], df[[args$`region-col`]]))
}, error = function(e) NULL)
if (!is.null(perm_tab)) {
	message("Available rows by permutation x region (non-zero counts):")
	nonzero <- perm_tab
	nonzero[nonzero == 0] <- NA
	utils::capture.output(print(nonzero)) |> paste(collapse = "\n") |> message()
}


original_self_done <- FALSE
for (reg in regions) {
	# Prepare region subset first (used for both summary and square modes)
	df$.__perm__ <- trimws(tolower(df[[args$`perm-col`]]))
	df$.__region__ <- trimws(tolower(df[[args$`region-col`]]))
	sub <- df[df$.__region__ == trimws(tolower(reg)), , drop = FALSE]

	# If 'original' does not exist in this region, log diagnostics but still proceed with controls
	if (!any(sub$.__perm__ == "original", na.rm = TRUE)) {
		message(sprintf("Region '%s': no 'original' rows after normalization. Available perms:", reg))
		utils::capture.output(print(table(sub$.__perm__, useNA = "ifany"))) |> paste(collapse = "\n") |> message()
	}

	present_ctrls <- intersect(controls, unique(sub$.__perm__))
	missing_ctrls <- setdiff(controls, present_ctrls)
	if (length(missing_ctrls) > 0) {
		message(sprintf("Skipping missing controls for region '%s': %s", reg, paste(missing_ctrls, collapse = ", ")))
		message("Perm counts in subset:")
		utils::capture.output(print(table(sub$.__perm__, useNA = "ifany"))) |> paste(collapse = "\n") |> message()
	}

	# Self-correlation heatmaps for negative controls present in this region (run regardless of 'original' presence)
	if (!only_original_csv) {
		for (perm_nm in present_ctrls) {
			cor_self <- compute_cor_square_self(sub, id_col = args$`id-col`, perm_col = ".__perm__", perm_name = tolower(perm_nm), feature_cols = feat_ordered)
			# Filename: <permutation>_<region>_featurecorr_heatmap_self.png
			out_png_self <- file.path(outdir, sprintf("%s_%s_featurecorr_heatmap_self.png", perm_nm, reg))
			x_axis <- sprintf("%s_%s", perm_nm, reg)
			y_axis <- sprintf("%s_%s", perm_nm, reg)
			save_heatmap_square(cor_self, sprintf("Self feature correlations: %s (%s)", perm_nm, reg), out_png_self, x_axis_name = x_axis, y_axis_name = y_axis, lower_triangle_only = lower_triangle_only, triangle_corner = triangle_corner)
		}
	}

	# Exactly one original self-correlation overall (prefer 'full' region if available)
	produced_original <- FALSE
	if (!original_self_done && any(sub$.__perm__ == "original", na.rm = TRUE)) {
		cor_self_orig <- compute_cor_square_self(sub, id_col = args$`id-col`, perm_col = ".__perm__", perm_name = "original", feature_cols = feat_ordered)
		# Always write CSV for original self-correlation
		out_csv_orig <- file.path(outdir, sprintf("original_%s_featurecorr_self_matrix.csv", reg))
		utils::write.csv(cor_self_orig, file = out_csv_orig, row.names = TRUE)

		# Optionally write PNG (skip if only-original-csv)
		if (!only_original_csv) {
			# Filename: original_<region>_featurecorr_heatmap_self.png
			out_png_self_orig <- file.path(outdir, sprintf("original_%s_featurecorr_heatmap_self.png", reg))
			x_axis <- sprintf("%s_%s", "original", reg)
			y_axis <- sprintf("%s_%s", "original", reg)
			save_heatmap_square(cor_self_orig, sprintf("Self feature correlations: original (%s)", reg), out_png_self_orig, x_axis_name = x_axis, y_axis_name = y_axis, lower_triangle_only = lower_triangle_only, triangle_corner = triangle_corner)
		}
		original_self_done <- TRUE
		produced_original <- TRUE
	}

	if (!only_original_csv) {
		message(sprintf("Region '%s': wrote %d self heatmaps for controls%s",
				reg, length(present_ctrls), if (produced_original) " (+1 original)" else ""))
	} else if (produced_original) {
		message(sprintf("Region '%s': wrote original self-correlation CSV", reg))
	}
}

