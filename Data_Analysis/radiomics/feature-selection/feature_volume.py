#!/usr/bin/env python3
"""
Concise filter for FMCIB radiomic features:
- Keep only original/full rows (readii_Permutation == 'original', readii_Region == 'full', case-insensitive)
- Keep only SampleID and pred_* feature columns
Writes: <outdir>/<base>_original_full_features_only.csv
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA, SparsePCA
try:
	import matplotlib.pyplot as plt
	MATPLOTLIB_AVAILABLE = True
except Exception:
	MATPLOTLIB_AVAILABLE = False

# Hardcoded schema to match feature_PCA.py across all datasets
ID_COL = "SampleID"
PERM_COL = "readii_Permutation"
REGION_COL = "readii_Region"
FEATURE_PREFIX = "pred_"


def parse_args() -> argparse.Namespace:
	p = argparse.ArgumentParser(description="Filter to original/full + pred_* and run volume-aware feature selection.")
	p.add_argument("--input", "-i", required=True, help="Path to input FMCIB features CSV")
	p.add_argument("--outdir", "-o", required=True, help="Output directory")
	p.add_argument("--volume-csv", required=True, help="Path to MIT index-simple CSV (volumes)")
	# Essential thresholds configurable
	p.add_argument("--var-bottom-perc", type=float, default=0.15, help="Drop bottom percentile of variance (0-1 or 0-100) [default: %(default)s]")
	p.add_argument("--abs-corr-threshold", type=float, default=0.30, help="Drop features with |Spearman rho| > threshold vs volume [default: %(default)s]")
	p.add_argument("--min-cumvar", type=float, default=0.95, help="Cumulative variance target for PCA components [default: %(default)s]")
	p.add_argument("--top-perc", type=float, default=0.20, help="Top percentage of features to select by weighted score (0-1 or 0-100) [default: %(default)s]")
	return p.parse_args()


def stem_from_path(path: str) -> str:
	name = os.path.basename(path)
	return name[:-4] if name.lower().endswith(".csv") else name

def infer_cohort_from_input(path: str) -> str | None:
	"""Infer cohort token from filenames like fmcib_features_fmcib_<COHORT>_index.csv"""
	base = stem_from_path(path)
	# Remove trailing _index suffix if present
	if base.endswith("_index"):
		base = base[:-6]
	# Expect pattern with last segment after 'fmcib_'
	if "fmcib_" in base:
		token = base.split("fmcib_")[-1]
		# Strip optional leading prefixes that may remain
		for pref in ("features_", "fmcib_"):
			if token.startswith(pref):
				token = token[len(pref):]
		return token or None
	return None


def main() -> int:
	args = parse_args()
	os.makedirs(args.outdir, exist_ok=True)

	print(f"Reading features: {args.input}")
	try:
		df = pd.read_csv(args.input)
	except Exception as e:
		print(f"ERROR: Failed to read input CSV: {e}", file=sys.stderr)
		return 2

	# Basic schema check
	missing = [c for c in (ID_COL, PERM_COL, REGION_COL) if c not in df.columns]
	if missing:
		print(f"ERROR: Missing required columns: {missing}", file=sys.stderr)
		return 2

	# Row filter (case-insensitive original/full)
	perm = df[PERM_COL].astype(str).str.strip().str.lower()
	region = df[REGION_COL].astype(str).str.strip().str.lower()
	sub = df.loc[(perm == "original") & (region == "full")].copy()
	if sub.empty:
		print("ERROR: No rows after filtering for original/full.", file=sys.stderr)
		return 3

	# Column filter: SampleID + pred_*
	feat_cols = [c for c in sub.columns if c.startswith(FEATURE_PREFIX)]
	if not feat_cols:
		print(f"ERROR: No feature columns starting with '{FEATURE_PREFIX}' after filtering.", file=sys.stderr)
		return 3
	keep_cols = ([ID_COL] if ID_COL in sub.columns else []) + feat_cols
	sub = sub.loc[:, keep_cols]

	# Numeric conversion and validation for feature columns
	sub[feat_cols] = sub[feat_cols].apply(pd.to_numeric, errors="coerce")
	X = sub[feat_cols].to_numpy(dtype=float)
	if not np.isfinite(X).all():
		n_nan = int(np.isnan(X).sum())
		n_inf = int(np.isinf(X).sum())
		print(f"ERROR: Non-finite values detected after numeric conversion (NaN={n_nan}, Inf={n_inf}).", file=sys.stderr)
		return 4

	out_csv = os.path.join(args.outdir, f"{stem_from_path(args.input)}_original_full_features_only.csv")
	sub.to_csv(out_csv, index=False)

	print(f"Filtered matrix: {out_csv}  (rows={sub.shape[0]}, features={len(feat_cols)})")

	# ===== Feature selection pipeline (volume-aware) =====
	# Read volume CSV from explicit flag
	vol_path = args.volume_csv
	print(f"Reading volumes: {vol_path}")
	try:
		vdf = pd.read_csv(vol_path)
	except Exception as e:
		print(f"ERROR: Failed to read volume CSV: {e}", file=sys.stderr)
		return 5

	# Build SampleID -> volume Series
	if 'SampleID' in vdf.columns:
		sid_series = vdf['SampleID'].astype(str).str.strip()
	elif 'Filepath' in vdf.columns:
		# Infer SampleID from first non-empty path component
		def _infer_sid(p: str) -> str:
			parts = str(p).strip().split('/')
			for comp in parts:
				if comp:
					return comp.split('.')[0].split('_')[0]
			return str(p).strip()
		sid_series = vdf['Filepath'].astype(str).map(_infer_sid)
	else:
		print("ERROR: Volume CSV must contain 'SampleID' or 'Filepath' column.", file=sys.stderr)
		return 5

	# Accept either 'sum' or 'original_shape_VoxelVolume' as the volume column
	if 'sum' in vdf.columns:
		vol_col = 'sum'
	elif 'original_shape_VoxelVolume' in vdf.columns:
		vol_col = 'original_shape_VoxelVolume'
	else:
		print("ERROR: Volume CSV must contain either 'sum' or 'original_shape_VoxelVolume' column.", file=sys.stderr)
		return 5

	vdf = vdf.copy()
	vdf['__SID__'] = sid_series
	# Optional modality filter if present
	if 'Modality' in vdf.columns:
		vdf = vdf[vdf['Modality'].isin(['RTSTRUCT', 'SEG'])]
	vol_series = vdf.groupby('__SID__')[vol_col].max()

	# Align volume to filtered features
	merged = sub.merge(vol_series.rename('volume'), left_on=ID_COL, right_index=True, how='left')
	n_before = merged.shape[0]
	merged = merged.dropna(subset=['volume'])
	n_after = merged.shape[0]
	print(f"Volume alignment: kept {n_after}/{n_before} samples with volume values")
	if merged.empty:
		print("ERROR: No samples with matched volume values after alignment.", file=sys.stderr)
		return 6

	# Step 2: Missing value filter (>50% missing -> drop)
	feat_df = merged[feat_cols]
	miss_frac = feat_df.isna().mean(axis=0)
	MISSING_THRESH = 0.50  # hardcoded per spec
	keep_mask = miss_frac <= MISSING_THRESH
	kept_feats = [f for f, k in zip(feat_cols, keep_mask) if bool(k)]
	if not kept_feats:
		print("ERROR: All features removed by missing value filter.", file=sys.stderr)
		return 7
	feat_df = feat_df[kept_feats]
	print(f"Missing filter: kept {len(kept_feats)} / {len(feat_cols)} features (threshold=50% missing)")

	# Impute remaining missing with feature median (robust for correlation/variance)
	feat_df = feat_df.apply(lambda s: s.fillna(s.median()), axis=0)

	# Step 3: Remove bottom 15% by variance
	var_vals = feat_df.var(axis=0, ddof=0).to_numpy(dtype=float)
	q = float(args.var_bottom_perc)
	q = q/100.0 if q > 1.0 else q
	cutoff = float(np.quantile(var_vals, q))
	keep_mask2 = var_vals >= cutoff
	kept_feats2 = [f for f, k in zip(feat_df.columns.tolist(), keep_mask2) if bool(k)]
	if not kept_feats2:
		# Fallback: keep top-variance feature
		top_idx = int(np.argmax(var_vals))
		kept_feats2 = [feat_df.columns.tolist()[top_idx]]
	feat_df = feat_df[kept_feats2]
	print(f"Variance filter: dropped bottom {q*100:.1f}% by variance; kept {feat_df.shape[1]} features (cutoff={cutoff:.4g})")

	# Step 4: Spearman correlation vs volume (remove |rho| > threshold)
	vol = merged['volume']
	spearman_rhos = feat_df.apply(lambda col: col.corr(vol, method='spearman'))
	keep_mask3 = spearman_rhos.abs() <= float(args.abs_corr_threshold)
	kept_feats3 = spearman_rhos.index[keep_mask3].tolist()
	feat_df = feat_df[kept_feats3]
	print(f"Volume correlation filter: kept {len(kept_feats3)} features with |rho| <= {args.abs_corr_threshold}")
	# Save correlation table
	corr_df = pd.DataFrame({
		'feature': spearman_rhos.index,
		'spearman_rho': spearman_rhos.values,
		'abs_spearman_rho': spearman_rhos.abs().values,
		'kept': keep_mask3.values,
	})
	corr_df.sort_values('abs_spearman_rho', ascending=True, inplace=True)
	corr_csv = os.path.join(args.outdir, f"{stem_from_path(args.input)}_volume_correlations.csv")
	corr_df.to_csv(corr_csv, index=False)

	if feat_df.shape[1] == 0:
		print("ERROR: No features remain after volume correlation filtering.", file=sys.stderr)
		return 8

	# Step 5: Standardize
	scaler = StandardScaler(with_mean=True, with_std=True)
	X_std = scaler.fit_transform(feat_df.values)

	# Step 6: PCA to reach min_cumvar
	n_samples, n_features = X_std.shape
	k_max = int(max(1, min(1024, n_samples, n_features)))
	pca = PCA(n_components=k_max, svd_solver="randomized", random_state=42)
	pca.fit(X_std)
	evr = np.asarray(pca.explained_variance_ratio_, dtype=float)
	cumvar = np.cumsum(evr)
	meets = np.where(cumvar >= float(args.min_cumvar))[0]
	if meets.size == 0:
		K = k_max
	else:
		K = int(meets[0] + 1)
	# Persist PCA cumulative variance
	base = stem_from_path(args.input)
	out_cumvar_csv = os.path.join(args.outdir, f"{base}_pca_cumulative_variance.csv")
	out_cumvar_png = os.path.join(args.outdir, f"{base}_pca_cumulative_variance.png")
	pd.DataFrame({"component": np.arange(1, cumvar.size + 1), "cum_explained_variance": cumvar}).to_csv(out_cumvar_csv, index=False)
	if MATPLOTLIB_AVAILABLE:
		plt.figure(figsize=(6,4), dpi=150)
		plt.plot(np.arange(1, cumvar.size + 1), cumvar, marker='o', linewidth=1)
		plt.axhline(float(args.min_cumvar), color='red', linestyle='--', linewidth=1)
		plt.xlabel('Number of components')
		plt.ylabel('Cumulative explained variance')
		plt.title('PCA cumulative explained variance')
		plt.grid(True, alpha=0.3)
		plt.tight_layout()
		plt.savefig(out_cumvar_png)
		plt.close()
	print(f"PCA: K={K} to reach {args.min_cumvar:.2f} cumulative variance (n_features={n_features})")

	# Step 7: SparsePCA with K components
	spca = SparsePCA(n_components=K, alpha=1.0, random_state=42)
	codes = spca.fit_transform(X_std)
	comps = spca.components_  # shape (K, n_features)

	# Save SPCA artifacts
	out_scores = os.path.join(args.outdir, f"{base}_sparsepca_scores.csv")
	out_loadings = os.path.join(args.outdir, f"{base}_sparsepca_loadings.csv")
	scores_df = pd.DataFrame(codes, columns=[f"SPC{i+1}" for i in range(codes.shape[1])])
	scores_df.insert(0, ID_COL, merged[ID_COL].to_list())
	scores_df.to_csv(out_scores, index=False)
	loadings_df = pd.DataFrame(comps.T, index=feat_df.columns.tolist(), columns=[f"SPC{i+1}" for i in range(comps.shape[0])])
	loadings_df.to_csv(out_loadings)

	# Step 8/9: Weighted scores and ranking
	weights = evr[:K].reshape(-1, 1)
	score_vec = np.sum(np.abs(comps) * weights, axis=0)
	score_df = pd.DataFrame({
		'feature': feat_df.columns.tolist(),
		'weighted_score': score_vec,
	}).sort_values('weighted_score', ascending=False, ignore_index=True)
	score_df['rank'] = np.arange(1, len(score_df) + 1)
	out_scores_ranked = os.path.join(args.outdir, f"{base}_feature_scores_ranked.csv")
	score_df.to_csv(out_scores_ranked, index=False)

	# Step 11: Select top 20% by score
	perc = float(args.top_perc)
	perc = perc/100.0 if perc > 1.0 else perc
	top_k = max(1, int(round(perc * score_df.shape[0])))
	selected = score_df.head(top_k).copy()
	out_selected_list = os.path.join(args.outdir, f"{base}_selected_features.txt")
	with open(out_selected_list, 'w') as f:
		for feat in selected['feature'].tolist():
			f.write(f"{feat}\n")
	# Also write selected feature matrix for convenience
	selected_matrix = pd.concat([merged[[ID_COL]].reset_index(drop=True), feat_df[selected['feature']].reset_index(drop=True)], axis=1)
	out_selected_matrix = os.path.join(args.outdir, f"{base}_selected_original_features.csv")
	selected_matrix.to_csv(out_selected_matrix, index=False)

	print("Done.")
	print(f"- Filtered matrix:         {out_csv}")
	print(f"- Volume correlations:     {corr_csv}")
	print(f"- PCA cumvar CSV/plot:     {out_cumvar_csv}{' / ' + out_cumvar_png if MATPLOTLIB_AVAILABLE else ''}")
	print(f"- SparsePCA scores:        {out_scores}")
	print(f"- SparsePCA loadings:      {out_loadings}")
	print(f"- Ranked feature scores:   {out_scores_ranked}")
	print(f"- Selected features (top {perc*100:.1f}%): {out_selected_list}")
	print(f"- Selected feature matrix: {out_selected_matrix}")
	return 0


if __name__ == "__main__":
	sys.exit(main())

