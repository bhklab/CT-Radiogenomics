#!/usr/bin/env python3
"""
Feature selection pipeline for FMCIB radiomic features using PCA + SparsePCA.

Steps:
1) Read CSV
2) Filter rows (readii_Permutation == 'original', readii_Region == 'full')
3) Keep SampleID + pred_*
4) Drop low-variance features (percentile threshold)
5) Standardize; run PCA to compute cumulative explained variance
6) Pick K via elbow on cumulative variance, then enforce a minimum cumulative variance
7) Fit SparsePCA with K components
8) Rank original features by variance-weighted absolute loadings; keep top percentage
9) Write filtered matrix, kept/removed lists, PCA cumvar (CSV/PNG), SparsePCA scores/loadings, selected features (+scores)
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Tuple
import time

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA, SparsePCA
try:
	import matplotlib.pyplot as plt  # for cumulative variance plot
	MATPLOTLIB_AVAILABLE = True
except Exception:
	MATPLOTLIB_AVAILABLE = False


def parse_args() -> argparse.Namespace:
	p = argparse.ArgumentParser(description="Filter FMCIB features and run Sparse PCA.")
	p.add_argument("--input", "-i", required=True, help="Path to input CSV file")
	p.add_argument("--outdir", "-o", required=True, help="Output directory")
	# Variance filter: absolute cutoff only
	p.add_argument("--var-threshold", type=float, default=0.15, help="Absolute variance cutoff: keep features with variance >= this value (computed on raw features with population variance ddof=0) [default: %(default)s]")
	p.add_argument("--min-cumvar", type=float, default=0.90, help="Minimum cumulative explained variance required after elbow selection [default: %(default)s]")
	p.add_argument("--select-top-perc", type=float, default=0.20, help="Top percentage of features to keep by variance-weighted absolute loadings (0–1 or 0–100); default 20%%")
	p.add_argument("--verbose", action="store_true", help="Enable verbose debug logging during each step")
	return p.parse_args()


def ensure_outdir(path: str) -> None:
	os.makedirs(path, exist_ok=True)


def stem_from_path(path: str) -> str:
	name = os.path.basename(path)
	if name.lower().endswith(".csv"):
		name = name[:-4]
	return name


def normalize_str_series(s: pd.Series) -> pd.Series:
	return s.astype(str).str.strip().str.lower()


def filter_rows_and_columns(df: pd.DataFrame, id_col: str, perm_col: str, region_col: str, feature_prefix: str) -> Tuple[pd.DataFrame, List[str]]:
	# Normalize labels and filter rows
	df = df.copy()
	df[perm_col] = normalize_str_series(df[perm_col])
	df[region_col] = normalize_str_series(df[region_col])
	sub = df[(df[perm_col] == "original") & (df[region_col] == "full")].copy()

	if sub.empty:
		raise ValueError("No rows after filtering for readii_Permutation=='original' and readii_Region=='full'.")

	# Keep SampleID + pred_* columns
	keep_cols = [c for c in sub.columns if c.startswith(feature_prefix)]
	if id_col in sub.columns:
		keep_cols = [id_col] + keep_cols
	if not any(c.startswith(feature_prefix) for c in keep_cols):
		raise ValueError(f"No feature columns starting with '{feature_prefix}' after filtering.")

	sub = sub.loc[:, keep_cols]

	feature_cols = [c for c in keep_cols if c != id_col]
	return sub, feature_cols


def to_numeric(df: pd.DataFrame, feature_cols: List[str]) -> pd.DataFrame:
	out = df.copy()
	for c in feature_cols:
		out[c] = pd.to_numeric(out[c], errors="coerce")
	return out


def variance_filter(
	X: np.ndarray,
	feature_cols: List[str],
	threshold: float,
	debug: bool = False,
) -> Tuple[np.ndarray, List[str], List[str]]:
	"""
	Filter features by absolute variance cutoff.
	Keeps features with population variance (ddof=0) >= threshold.
	"""
	# Compute per-feature population variance
	vars_ = np.var(X, axis=0, ddof=0)
	cutoff = float(threshold)

	if debug:
		vmin = float(np.min(vars_))
		vmed = float(np.median(vars_))
		vq85 = float(np.quantile(vars_, 0.85))
		vmax = float(np.max(vars_))
		print(f"[DEBUG] variance_filter: mode=abs, threshold={threshold}, cutoff={cutoff:.6g}, var[min/50%/85%/max]=[{vmin:.4g}/{vmed:.4g}/{vq85:.4g}/{vmax:.4g}]")

	mask = vars_ >= cutoff
	kept = [f for f, m in zip(feature_cols, mask) if bool(m)]
	removed = [f for f, m in zip(feature_cols, mask) if not bool(m)]
	if not any(mask):
		# Fallback to keep the single highest-variance feature to avoid empty set
		top_idx = int(np.argmax(vars_))
		mask = np.zeros_like(vars_, dtype=bool)
		mask[top_idx] = True
		kept = [feature_cols[top_idx]]
		removed = [f for i, f in enumerate(feature_cols) if i != top_idx]

	X_sel = X[:, mask]
	return X_sel, kept, removed


def fit_sparse_pca(X: np.ndarray, n_components: int):
	# Hard-coded SparsePCA with fixed seed for reproducibility
	model = SparsePCA(n_components=n_components, alpha=1.0, random_state=42)
	codes = model.fit_transform(X)
	components = model.components_  # shape (n_components, n_features)
	return model, codes, components


def elbow_from_cumvar(cumvar: np.ndarray, debug: bool = False) -> int:
	"""Elbow selection from a given cumulative variance array (length k_max).
	Returns k >= 1.
	"""
	if cumvar.size == 0:
		return 1
	x = np.arange(1, len(cumvar) + 1, dtype=float)
	x_norm = (x - x[0]) / (x[-1] - x[0]) if len(x) > 1 else x * 0.0
	y_norm = (cumvar - cumvar[0]) / (cumvar[-1] - cumvar[0] + 1e-12)
	x0, y0 = 0.0, 0.0
	x1, y1 = 1.0, 1.0
	denom = np.hypot(x1 - x0, y1 - y0)
	if denom <= 0:
		return 1
	distances = np.abs((y1 - y0) * x_norm - (x1 - x0) * y_norm + x1 * y0 - y1 * x0) / denom
	idx = int(np.argmax(distances))
	k = int(idx + 1)
	if debug:
		print(f"[DEBUG] elbow(cumvar): selected_k={k}, cumvar_at_k={float(cumvar[idx]):.4f}")
	return max(1, k)


def compute_feature_importance_weighted(components: np.ndarray, comp_weights: np.ndarray) -> np.ndarray:
	"""Variance-weighted importance: sum_j |loading_{j,feature}| * weight_j
	components: shape (n_components, n_features)
	comp_weights: shape (n_components,) — typically PCA explained_variance_ratio_ for first K comps
	Returns: shape (n_features,)
	"""
	A = np.abs(components)
	w = np.asarray(comp_weights, dtype=float).reshape(-1, 1)
	return np.sum(A * w, axis=0)


def main() -> int:
	args = parse_args()
	ensure_outdir(args.outdir)

	# Simple verbose logger with timestamps
	def vlog(msg: str) -> None:
		if args.verbose:
			print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] [DEBUG] {msg}")

	print(f"Reading: {args.input}")
	read_t0 = time.perf_counter()
	df = pd.read_csv(args.input)
	read_dt = time.perf_counter() - read_t0
	vlog(f"read_csv done in {read_dt:.2f}s; shape={df.shape}")

	# Constants for column names and feature prefix (hardcoded for consistency)
	ID_COL = "SampleID"
	PERM_COL = "readii_Permutation"
	REGION_COL = "readii_Region"
	FEATURE_PREFIX = "pred_"

	# Basic checks
	required_cols = [ID_COL, PERM_COL, REGION_COL]
	missing = [c for c in required_cols if c not in df.columns]
	if missing:
		print(f"ERROR: Missing required columns: {missing}", file=sys.stderr)
		return 2

	# Filter rows and columns
	flt_t0 = time.perf_counter()
	sub, feature_cols = filter_rows_and_columns(df, ID_COL, PERM_COL, REGION_COL, FEATURE_PREFIX)
	flt_dt = time.perf_counter() - flt_t0
	print(f"After row/column filter: {sub.shape[0]} rows, {len(feature_cols)} features")
	vlog(f"row/column filter took {flt_dt:.2f}s; kept feature prefix='{FEATURE_PREFIX}'")

	# Numeric conversion
	num_t0 = time.perf_counter()
	sub_num = to_numeric(sub, feature_cols)
	X = sub_num[feature_cols].to_numpy(dtype=float)
	num_dt = time.perf_counter() - num_t0
	vlog(f"numeric conversion took {num_dt:.2f}s; matrix shape={X.shape}")
	# Validate no missing or non-finite values
	if not np.isfinite(X).all():
		bad_nan = int(np.isnan(X).sum())
		bad_inf = int(np.isinf(X).sum())
		print(
			f"ERROR: Non-finite values detected after numeric conversion (NaN={bad_nan}, Inf={bad_inf}).",
			file=sys.stderr,
		)
		return 3

	# Variance thresholding on raw values
	var_t0 = time.perf_counter()
	if X.shape[1] == 0:
		print("ERROR: No features available before variance filtering.", file=sys.stderr)
		return 3
	X_vf, kept_feats, removed_feats = variance_filter(X, feature_cols, threshold=args.var_threshold, debug=args.verbose)
	disp_thr = f">= {args.var_threshold:.4g} variance"
	print(f"Variance filter: kept {len(kept_feats)}, removed {len(removed_feats)} (threshold {disp_thr})")
	var_dt = time.perf_counter() - var_t0
	vlog(f"variance filtering took {var_dt:.2f}s; new shape={X_vf.shape}")
	if X_vf.shape[1] == 0:
		print("ERROR: All features removed by variance thresholding.", file=sys.stderr)
		return 4

	# Standardize features (post-threshold)
	scale_t0 = time.perf_counter()
	scaler = StandardScaler(with_mean=True, with_std=True)
	X_std = scaler.fit_transform(X_vf)
	scale_dt = time.perf_counter() - scale_t0
	vlog(f"standardization took {scale_dt:.2f}s; matrix shape={X_std.shape}")

	# Determine number of components: run regular PCA, plot cumulative variance, elbow -> k, then use k for SparsePCA
	n_samples, n_features = X_std.shape
	# Hardcoded elbow scan cap for consistency
	k_max = int(max(1, min(1024, n_samples, n_features)))
	pca_t0 = time.perf_counter()
	pca = PCA(n_components=k_max, svd_solver="randomized", random_state=42)
	pca.fit(X_std)
	pca_dt = time.perf_counter() - pca_t0
	vlog(f"PCA (for cumvar) fit took {pca_dt:.2f}s; k_max={k_max}")
	evrs = np.asarray(pca.explained_variance_ratio_, dtype=float)
	cumvar = np.cumsum(evrs)
	# Save cumulative variance plot and CSV
	base = stem_from_path(args.input)
	out_cumvar_csv = os.path.join(args.outdir, f"{base}_pca_cumulative_variance.csv")
	out_cumvar_png = os.path.join(args.outdir, f"{base}_pca_cumulative_variance.png")
	pd.DataFrame({"component": np.arange(1, cumvar.size + 1), "cum_explained_variance": cumvar}).to_csv(out_cumvar_csv, index=False)
	if MATPLOTLIB_AVAILABLE:
		plt.figure(figsize=(6,4), dpi=150)
		plt.plot(np.arange(1, cumvar.size + 1), cumvar, marker='o', linewidth=1)
		plt.xlabel('Number of components')
		plt.ylabel('Cumulative explained variance')
		plt.title('PCA cumulative explained variance')
		plt.grid(True, alpha=0.3)
		plt.tight_layout()
		plt.savefig(out_cumvar_png)
		plt.close()
		vlog(f"saved PCA cumulative variance plot: {out_cumvar_png}")
	else:
		vlog("matplotlib not available; saved CSV but skipped plot image")
	# Elbow from cumulative variance
	elbow_t0 = time.perf_counter()
	n_comp = elbow_from_cumvar(cumvar, debug=args.verbose)
	n_comp = min(n_comp, n_samples, n_features)
	print(f"Elbow-selected n_components (from PCA cumvar) = {n_comp} (shape={X_std.shape})")
	vlog(f"elbow selection from cumvar took {time.perf_counter()-elbow_t0:.2f}s")
	# Ensure minimum cumulative explained variance threshold
	cum_at_k = float(cumvar[n_comp - 1]) if n_comp >= 1 else 0.0
	if cum_at_k < float(args.min_cumvar):
		# Find smallest k achieving the threshold
		meets = np.where(cumvar >= float(args.min_cumvar))[0]
		if meets.size > 0:
			k2 = int(meets[0] + 1)
			old_k = n_comp
			n_comp = min(k2, n_samples, n_features)
			print(f"Adjusting n_components to meet min cumvar {args.min_cumvar:.2f}: {old_k} -> {n_comp} (cumvar={cumvar[n_comp-1]:.4f})")
			vlog(f"min-cumvar adjustment: selected k bumped to {n_comp}")
		else:
			vlog("min-cumvar requested but not achievable within k_max; proceeding with elbow k")
	# no manual n_components override for consistency

	# Fit SparsePCA
	print(f"Fitting SparsePCA: n_components={n_comp}, alpha=1.0")
	spca_t0 = time.perf_counter()
	model, codes, components = fit_sparse_pca(X_std, n_components=n_comp)
	spca_dt = time.perf_counter() - spca_t0
	vlog(f"SparsePCA fit_transform took {spca_dt:.2f}s; codes shape={codes.shape}")

	# Prepare outputs
	base = stem_from_path(args.input)
	out_filtered_csv = os.path.join(args.outdir, f"{base}_original_full_filtered.csv")
	out_lowvar_kept = os.path.join(args.outdir, f"{base}_lowvar_kept_features.txt")
	out_lowvar_removed = os.path.join(args.outdir, f"{base}_lowvar_removed_features.txt")
	out_scores = os.path.join(args.outdir, f"{base}_sparsepca_scores.csv")
	out_loadings = os.path.join(args.outdir, f"{base}_sparsepca_loadings.csv")
	out_selected_list = os.path.join(args.outdir, f"{base}_spca_selected_features.txt")
	out_selected_matrix = os.path.join(args.outdir, f"{base}_spca_selected_original_features.csv")

	# Filtered matrix (rows: samples; cols: SampleID + kept features after variance filter)
	io_t0 = time.perf_counter()
	filtered_df = pd.concat([sub_num[[ID_COL]].reset_index(drop=True), pd.DataFrame(X_vf, columns=kept_feats)], axis=1)
	filtered_df.to_csv(out_filtered_csv, index=False)

	# Feature lists
	with open(out_lowvar_kept, "w") as f:
		for feat in kept_feats:
			f.write(f"{feat}\n")
	with open(out_lowvar_removed, "w") as f:
		for feat in removed_feats:
			f.write(f"{feat}\n")

	# Sample scores (codes): add SampleID
	scores_df = pd.DataFrame(codes, columns=[f"SPC{i+1}" for i in range(codes.shape[1])])
	scores_df.insert(0, ID_COL, sub_num[ID_COL].to_list())
	scores_df.to_csv(out_scores, index=False)

	# Component loadings: features x components
	loadings_df = pd.DataFrame(components.T, index=kept_feats, columns=[f"SPC{i+1}" for i in range(components.shape[0])])
	loadings_df.to_csv(out_loadings)

	# Feature selection from SparsePCA loadings (keep original features)
	# Variance-weighted importance per feature: sum_j |loading_j| * PCA_EV_ratio_j
	comp_weights = np.asarray(evrs[:n_comp], dtype=float)
	feat_importance = compute_feature_importance_weighted(components, comp_weights)  # shape: (n_features,)
	if args.verbose:
		imp_min = float(np.min(feat_importance))
		imp_p25 = float(np.quantile(feat_importance, 0.25))
		imp_med = float(np.median(feat_importance))
		imp_p75 = float(np.quantile(feat_importance, 0.75))
		imp_max = float(np.max(feat_importance))
		print(f"[DEBUG] importance metric=var-weighted-sumabs: min/p25/med/p75/max = {imp_min:.4g}/{imp_p25:.4g}/{imp_med:.4g}/{imp_p75:.4g}/{imp_max:.4g}")
	selection_note = ""
	selected_feats: List[str]
	perc = float(args.select_top_perc)
	if perc > 1.0:
		perc = perc / 100.0
	perc = max(0.0, min(1.0, perc))
	k = max(1, int(round(perc * len(kept_feats))))
	idx_sorted = np.argsort(feat_importance)[::-1][:k]
	selected_feats = [kept_feats[i] for i in idx_sorted]
	selection_note = f"top-perc={perc*100:.1f}% (K={k}) by var-weighted-sumabs"
	# Persist selection
	with open(out_selected_list, "w") as f:
		for feat in selected_feats:
			f.write(f"{feat}\n")
	selected_df = pd.concat([
		sub_num[[ID_COL]].reset_index(drop=True),
		sub_num[selected_feats].reset_index(drop=True),
	], axis=1)
	selected_df.to_csv(out_selected_matrix, index=False)
	# Persist selection scores as CSV (feature, importance)
	out_selected_scores = os.path.join(args.outdir, f"{base}_spca_selected_features_with_scores.csv")
	sel_scores_df = pd.DataFrame({"feature": [kept_feats[i] for i in idx_sorted], "importance": feat_importance[idx_sorted]})
	sel_scores_df.to_csv(out_selected_scores, index=False)
	vlog(f"feature selection by SparsePCA loadings: {selection_note}; wrote {len(selected_feats)} features and scores")
	io_dt = time.perf_counter() - io_t0
	vlog(f"wrote outputs in {io_dt:.2f}s")

	print("Done.")
	print(f"- Filtered matrix:   {out_filtered_csv}")
	print(f"- Kept features:     {out_lowvar_kept}")
	print(f"- Removed features:  {out_lowvar_removed}")
	print(f"- SparsePCA scores:  {out_scores}")
	print(f"- SparsePCA loadings:{out_loadings}")
	print(f"- Selected features:  {out_selected_list}")
	print(f"- Selected matrix:    {out_selected_matrix}")
	print(f"- PCA cumvar CSV:     {out_cumvar_csv}")
	if MATPLOTLIB_AVAILABLE:
		print(f"- PCA cumvar plot:    {out_cumvar_png}")
	print(f"- Selected scores:    {out_selected_scores}")
	return 0


if __name__ == "__main__":
	sys.exit(main())

