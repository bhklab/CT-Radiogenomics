#!/usr/bin/env python3
"""
Feature selection by hierarchical clustering of FMCIB radiomic features.
After clustering, the medoid feature of each cluster is selected to represent the cluster (feature selection).

Steps:
1) Read CSV
2) Filter rows (readii_Permutation == 'original', readii_Region == 'full')
3) Keep SampleID + pred_*
4) Drop low-variance features by absolute variance cutoff (population variance, ddof=0)
5) Standardize features (z-score)
6) Compute pairwise correlation distances between features (1 - Spearman ρ)
7) Hierarchical clustering over features; cut into K clusters (K chosen automatically)
	- K selection methods:
	  * elbow (default): knee point of within-cluster dispersion vs k
	  * gap: Gap Statistic (Tibshirani et al., 2001)
	  * silhouette: average silhouette over features (with precomputed distances)
8) Visualize the full dendrogram of the hierarchical clustering (saved to disk)
9) Select medoid per cluster (feature minimizing total distance to others in the cluster)
10) Write outputs: filtered matrix, kept/removed lists, dendrogram image(s), cluster labels, medoid list, medoid matrix
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Tuple, Optional
import time

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score

try:
	# Prefer SciPy for hierarchical clustering
	from scipy.spatial.distance import pdist, squareform
	from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
	SCIPY_AVAILABLE = True
except Exception:
	SCIPY_AVAILABLE = False

try:
	import matplotlib.pyplot as plt  # optional dendrogram or diagnostics later
	MATPLOTLIB_AVAILABLE = True
except Exception:
	MATPLOTLIB_AVAILABLE = False


def parse_args() -> argparse.Namespace:
	p = argparse.ArgumentParser(description="Filter FMCIB features and select medoid features via hierarchical clustering.")
	p.add_argument("--input", "-i", required=True, help="Path to input CSV file")
	p.add_argument("--outdir", "-o", required=True, help="Output directory")
	# Variance filter: absolute cutoff only (consistent with feature_PCA.py)
	p.add_argument("--var-threshold", type=float, default=0.15, help="Absolute variance cutoff: keep features with variance >= this value (population variance, ddof=0) [default: %(default)s]")
	# Clustering controls (distance metric fixed to correlation)
	p.add_argument("--linkage", choices=["average", "complete", "single"], default="average", help="Hierarchical linkage method [default: %(default)s]")
	# K selection method
	p.add_argument("--k-method", choices=["elbow", "gap", "silhouette"], default="elbow", help="Method to select number of clusters K [default: %(default)s]")
	p.add_argument("--kmax", type=int, default=None, help="Maximum K to consider for elbow/silhouette/gap (default: min(100, n_features-1))")
	# Gap statistic controls
	p.add_argument("--gap-B", type=int, default=30, help="Number of reference datasets for Gap Statistic (higher is slower, more stable) [default: %(default)s]")
	# Visualization control
	p.add_argument("--no-dendrogram", action="store_true", help="Skip saving the full hierarchical clustering dendrogram plot")
	p.add_argument("--verbose", action="store_true", help="Enable verbose debug logging")
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


def pairwise_feature_distance(X_std: np.ndarray) -> np.ndarray:
		"""
		Compute pairwise Spearman correlation distance (1 - Spearman ρ) between features.
		Steps:
			- Rank-transform each feature across samples (average ranks for ties)
			- Compute Pearson correlation on the rank matrix (equivalent to Spearman)
			- Convert to distance D = 1 - ρ
		Returns an (n_features x n_features) distance matrix (float64).
		"""
		# Rank-transform each column across samples
		ranks = pd.DataFrame(X_std).rank(axis=0, method="average", na_option="keep").to_numpy(dtype=float)
		# Pearson correlation on ranks equals Spearman correlation on original data
		C = np.corrcoef(ranks, rowvar=False)
		C = np.clip(C, -1.0, 1.0)
		D = 1.0 - C
		np.fill_diagonal(D, 0.0)
		D[D < 0] = 0.0
		return D.astype(float, copy=False)


def hierarchical_cluster_labels(D: np.ndarray, n_clusters: int, linkage_method: str, debug: bool = False) -> Tuple[np.ndarray, Optional[np.ndarray]]:
	"""
	Cluster features using hierarchical clustering given a pairwise distance matrix D.
	Returns (labels, Z) where labels are integer cluster labels (0..n_clusters-1) for each feature (len = n_features),
	and Z is the scipy linkage matrix (or None if SciPy unavailable).
	"""
	n_features = D.shape[0]
	if n_clusters < 1:
		n_clusters = 1
	if n_clusters > n_features:
		n_clusters = n_features

	if SCIPY_AVAILABLE:
		condensed = squareform(D, checks=False)
		Z = linkage(condensed, method=linkage_method)
		# Cut into n_clusters
		labels = fcluster(Z, t=n_clusters, criterion="maxclust") - 1  # 0-based
		return labels.astype(int), Z

	# Fallback to a simple greedy clustering if SciPy is unavailable
	if debug:
		print("[DEBUG] SciPy not available; using a simple greedy clustering fallback.")
	labels = -np.ones(n_features, dtype=int)
	assigned = 0
	current_label = 0
	# seed with first unassigned, assign closest until cluster size approx n_features/n_clusters
	target_size = max(1, int(np.ceil(n_features / float(n_clusters))))
	while assigned < n_features and current_label < n_clusters:
		seed = int(np.argmax(labels == -1))  # first unassigned index
		labels[seed] = current_label
		# assign nearest unassigned up to target_size
		d_row = D[seed]
		order = np.argsort(d_row)
		cnt = 1
		for j in order:
			if labels[j] == -1:
				labels[j] = current_label
				cnt += 1
				assigned += 1
				if cnt >= target_size:
					break
		current_label += 1
	# Any remaining unassigned -> last label
	labels[labels == -1] = min(current_label, n_clusters - 1)
	return labels, None


def save_dendrogram(Z: np.ndarray, out_png: str, out_pdf: Optional[str] = None, debug: bool = False) -> None:
	"""Save full hierarchical clustering dendrogram from linkage matrix Z."""
	if not MATPLOTLIB_AVAILABLE:
		if debug:
			print("[DEBUG] Matplotlib not available; skipping dendrogram plot.")
		return
	# Create large-ish figure; no labels for readability and performance
	fig_h = 6.0
	fig_w = 12.0
	plt.figure(figsize=(fig_w, fig_h))
	dendrogram(Z, no_labels=True, count_sort=False, distance_sort=False, show_leaf_counts=False)
	plt.tight_layout()
	plt.savefig(out_png, dpi=200)
	if out_pdf:
		plt.savefig(out_pdf)
	plt.close()


def compute_linkage_matrix(D: np.ndarray, linkage_method: str) -> Optional[np.ndarray]:
	"""Return the SciPy linkage matrix for hierarchical clustering, or None if SciPy unavailable."""
	if not SCIPY_AVAILABLE:
		return None
	condensed = squareform(D, checks=False)
	return linkage(condensed, method=linkage_method)


def _wk_for_k(Dsq: np.ndarray, Z: Optional[np.ndarray], k: int, linkage_method: str) -> Tuple[float, np.ndarray]:
	"""Compute Wk (within dispersion) for k clusters using either provided Z or recomputing from Dsq."""
	if Z is None:
		if not SCIPY_AVAILABLE:
			labels = np.zeros(Dsq.shape[0], dtype=int)
			return float(np.nan), labels
		Z = compute_linkage_matrix(Dsq, linkage_method)
	labels = fcluster(Z, t=int(k), criterion="maxclust") - 1
	Wk = _within_cluster_dispersion(Dsq, labels)
	return Wk, labels


def _save_elbow_plot(out_png: str, out_pdf: Optional[str], ks: np.ndarray, wk: np.ndarray, k_star: int, log_scale: bool = True) -> None:
	if not MATPLOTLIB_AVAILABLE:
		return
	y = np.log(np.maximum(wk, 1e-12)) if log_scale else wk
	plt.figure(figsize=(7.0, 4.5))
	plt.plot(ks, y, marker="o", label="log Wk" if log_scale else "Wk")
	plt.axvline(k_star, color="red", linestyle="--", label=f"K*={k_star}")
	plt.xlabel("k (number of clusters)")
	plt.ylabel("log(Wk)" if log_scale else "Wk")
	plt.legend()
	plt.tight_layout()
	plt.savefig(out_png, dpi=200)
	if out_pdf:
		plt.savefig(out_pdf)
	plt.close()


def _save_silhouette_plot(out_png: str, out_pdf: Optional[str], ks: np.ndarray, sil: np.ndarray, k_star: int) -> None:
	if not MATPLOTLIB_AVAILABLE:
		return
	plt.figure(figsize=(7.0, 4.5))
	plt.plot(ks, sil, marker="o", label="avg silhouette")
	plt.axvline(k_star, color="red", linestyle="--", label=f"K*={k_star}")
	plt.xlabel("k (number of clusters)")
	plt.ylabel("silhouette score")
	plt.legend()
	plt.tight_layout()
	plt.savefig(out_png, dpi=200)
	if out_pdf:
		plt.savefig(out_pdf)
	plt.close()


def choose_k_via_elbow(Dsq: np.ndarray, Z: Optional[np.ndarray], linkage_method: str, Kmax: int, debug: bool = False) -> Tuple[int, np.ndarray, np.ndarray]:
	"""Elbow method: compute Wk for k=1..Kmax, use knee via maximum curvature on log(Wk).
	Returns (k*, ks, Wk_arr)."""
	ks = np.arange(1, Kmax + 1, dtype=int)
	Wk = np.zeros(ks.size, dtype=float)
	for i, k in enumerate(ks):
		Wk[i], _ = _wk_for_k(Dsq, Z, int(k), linkage_method)
	# Use log transform to stabilize
	y = np.log(np.maximum(Wk, 1e-12))
	# Discrete second derivative (curvature)
	if y.size >= 3:
		curv = y[2:] - 2 * y[1:-1] + y[:-2]
		best_mid = int(np.argmax(curv))  # index in 0..len(curv)-1 corresponds to k=ks[best_mid+1]
		k_star = int(ks[best_mid + 1])
	else:
		# Fallback: largest relative drop in Wk
		delta = -np.diff(y, prepend=y[0])
		k_star = int(ks[int(np.argmax(delta))])
	# Guardrails
	k_star = max(1, min(int(k_star), int(Kmax)))
	if debug:
		print(f"[DEBUG] elbow: Kmax={Kmax}, k*={k_star}")
	return k_star, ks, Wk


def choose_k_via_silhouette(Dsq: np.ndarray, Z: Optional[np.ndarray], linkage_method: str, Kmax: int, debug: bool = False) -> Tuple[int, np.ndarray, np.ndarray]:
	"""Average silhouette method: evaluate k=2..Kmax, choose k with max silhouette.
	Returns (k*, ks_eval, sil_scores)."""
	# silhouette is undefined for k=1; start at 2 if possible
	start_k = 2 if Kmax >= 2 else 1
	ks = np.arange(start_k, Kmax + 1, dtype=int)
	if ks.size == 0:
		return 1, np.array([1]), np.array([0.0])
	sil_scores = np.full(ks.size, -1.0, dtype=float)
	for i, k in enumerate(ks):
		_, labels = _wk_for_k(Dsq, Z, int(k), linkage_method)
		# silhouette_score expects a distance matrix with shape (n_samples, n_samples)
		try:
			sil_scores[i] = silhouette_score(Dsq, labels, metric="precomputed")
		except Exception:
			sil_scores[i] = -1.0
	best_idx = int(np.argmax(sil_scores))
	k_star = int(ks[best_idx])
	if debug:
		print(f"[DEBUG] silhouette: Kmax={Kmax}, best_k={k_star}, best_score={sil_scores[best_idx]:.4f}")
	return k_star, ks, sil_scores


def _within_cluster_dispersion(Dsq: np.ndarray, labels: np.ndarray) -> float:
	"""Compute within-cluster dispersion Wk given a full square distance matrix and cluster labels.
	Wk = sum_r ( 1/(2 n_r) * sum_{i,j in C_r} D[i,j] ). For singleton clusters, contribution is 0."""
	W = 0.0
	for lab in np.unique(labels):
		idx = np.where(labels == lab)[0]
		m = idx.size
		if m <= 1:
			continue
		subD = Dsq[np.ix_(idx, idx)]
		W += float(subD.sum()) / (2.0 * float(m))
	return float(W)


def _permute_columns_independently(X: np.ndarray, rng: np.random.Generator) -> np.ndarray:
	"""Generate a reference dataset by independently permuting rows within each column (feature).
	Preserves per-feature marginal distribution while destroying inter-feature dependence."""
	X_ref = np.empty_like(X)
	for j in range(X.shape[1]):
		X_ref[:, j] = rng.permutation(X[:, j])
	return X_ref


def _save_gap_plot(out_png: str, out_pdf: Optional[str], ks: np.ndarray, gaps: np.ndarray, sk: np.ndarray, k_star: int) -> None:
	if not MATPLOTLIB_AVAILABLE:
		return
	plt.figure(figsize=(7.0, 4.5))
	plt.plot(ks, gaps, marker="o", label="Gap(k)")
	plt.fill_between(ks, gaps - sk, gaps + sk, color="gray", alpha=0.2, label="± s(k)")
	plt.axvline(k_star, color="red", linestyle="--", label=f"K*={k_star}")
	plt.xlabel("k (number of clusters)")
	plt.ylabel("Gap")
	plt.legend()
	plt.tight_layout()
	plt.savefig(out_png, dpi=200)
	if out_pdf:
		plt.savefig(out_pdf)
	plt.close()


def choose_k_via_gap_statistic(
	X_std: np.ndarray,
	linkage_method: str,
	B: int = 10,
	Kmax: Optional[int] = None,
	random_state: int = 42,
	debug: bool = False,
) -> Tuple[int, Optional[np.ndarray], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
	"""
	Compute K using the Gap Statistic on hierarchical clustering with Spearman distance.
	Returns: (K*, Z_original, Wk, gap, sk, ks)
	- Z_original: linkage matrix for the original data (None if SciPy unavailable)
	- Wk: within-cluster dispersion values (length = len(ks))
	- gap: Gap(k)
	- sk: std * sqrt(1+1/B)
	- ks: array of k values evaluated
	"""
	if not SCIPY_AVAILABLE:
		# Fallback if SciPy is absent
		n_features = X_std.shape[1]
		k_fallback = max(1, int(round(min(100, n_features - 1))))
		return k_fallback, None, np.array([np.nan]), np.array([np.nan]), np.array([np.nan]), np.array([k_fallback])

	# Distance and linkage on original data
	D = pairwise_feature_distance(X_std)
	Z = compute_linkage_matrix(D, linkage_method=linkage_method)
	if Z is None:
		# Shouldn't happen if SCIPY_AVAILABLE, but guard anyway
		n_features = X_std.shape[1]
		k_fallback = max(1, int(round(min(100, n_features - 1))))
		return k_fallback, None, np.array([np.nan]), np.array([np.nan]), np.array([np.nan]), np.array([k_fallback])

	n_features = X_std.shape[1]
	if Kmax is None:
		Kmax = max(1, min(100, n_features - 1))
	ks = np.arange(1, Kmax + 1, dtype=int)
	Dsq = D  # already square

	# Wk for original data across ks
	Wk = np.zeros(ks.size, dtype=float)
	for i, k in enumerate(ks):
		labels_k = fcluster(Z, t=int(k), criterion="maxclust") - 1
		Wk[i] = _within_cluster_dispersion(Dsq, labels_k)

	# Reference dispersions
	rng = np.random.default_rng(seed=random_state)
	logWk_ref = np.zeros((B, ks.size), dtype=float)
	for b in range(B):
		X_ref = _permute_columns_independently(X_std, rng)
		D_ref = pairwise_feature_distance(X_ref)
		Z_ref = compute_linkage_matrix(D_ref, linkage_method=linkage_method)
		Dsq_ref = D_ref
		for i, k in enumerate(ks):
			labels_ref = fcluster(Z_ref, t=int(k), criterion="maxclust") - 1
			Wkb = _within_cluster_dispersion(Dsq_ref, labels_ref)
			logWk_ref[b, i] = np.log(max(Wkb, 1e-12))
		if debug:
			print(f"[DEBUG] Gap ref {b+1}/{B} done")

	logWk = np.log(np.maximum(Wk, 1e-12))
	gap = np.mean(logWk_ref, axis=0) - logWk
	sk = np.std(logWk_ref, axis=0, ddof=1) * np.sqrt(1.0 + 1.0 / float(B))

	# 1-SE rule: choose the smallest k such that Gap(k) >= Gap(k+1) - s(k+1)
	k_star = int(ks[np.argmax(gap)])  # fallback
	for i in range(ks.size - 1):
		if gap[i] >= gap[i + 1] - sk[i + 1]:
			k_star = int(ks[i])
			break

	return k_star, Z, Wk, gap, sk, ks


def medoids_from_clusters(D: np.ndarray, labels: np.ndarray) -> List[int]:
	"""
	Given a distance matrix D (n_features x n_features) and cluster labels per feature,
	return the medoid index for each cluster: argmin_i sum_j D[i,j] within cluster.
	"""
	medoid_indices: List[int] = []
	for lab in sorted(np.unique(labels)):
		idx = np.where(labels == lab)[0]
		if idx.size == 1:
			medoid_indices.append(int(idx[0]))
			continue
		subD = D[np.ix_(idx, idx)]
		# Sum distances per candidate
		sums = np.sum(subD, axis=1)
		m_local = int(idx[int(np.argmin(sums))])
		medoid_indices.append(m_local)
	return medoid_indices


def main() -> int:
	args = parse_args()
	ensure_outdir(args.outdir)

	def vlog(msg: str) -> None:
		if args.verbose:
			print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] [DEBUG] {msg}")

	print(f"Reading: {args.input}")
	t0 = time.perf_counter()
	df = pd.read_csv(args.input)
	vlog(f"read_csv done in {time.perf_counter()-t0:.2f}s; shape={df.shape}")

	# Constants
	ID_COL = "SampleID"
	PERM_COL = "readii_Permutation"
	REGION_COL = "readii_Region"
	FEATURE_PREFIX = "pred_"

	# Checks
	required = [ID_COL, PERM_COL, REGION_COL]
	missing = [c for c in required if c not in df.columns]
	if missing:
		print(f"ERROR: Missing required columns: {missing}", file=sys.stderr)
		return 2

	# Filter & numeric
	t1 = time.perf_counter()
	sub, feature_cols = filter_rows_and_columns(df, ID_COL, PERM_COL, REGION_COL, FEATURE_PREFIX)
	vlog(f"row/column filter took {time.perf_counter()-t1:.2f}s; kept {len(feature_cols)} features")
	print(f"After row/column filter: {sub.shape[0]} rows, {len(feature_cols)} features")

	t2 = time.perf_counter()
	sub_num = to_numeric(sub, feature_cols)
	X = sub_num[feature_cols].to_numpy(dtype=float)
	vlog(f"numeric conversion took {time.perf_counter()-t2:.2f}s; matrix shape={X.shape}")
	if not np.isfinite(X).all():
		nan_c = int(np.isnan(X).sum())
		inf_c = int(np.isinf(X).sum())
		print(f"ERROR: Non-finite values detected after numeric conversion (NaN={nan_c}, Inf={inf_c}).", file=sys.stderr)
		return 3

	# Variance filter (absolute)
	if X.shape[1] == 0:
		print("ERROR: No features available before variance filtering.", file=sys.stderr)
		return 3
	t3 = time.perf_counter()
	X_vf, kept_feats, removed_feats = variance_filter(X, feature_cols, threshold=args.var_threshold, debug=args.verbose)
	print(f"Variance filter: kept {len(kept_feats)}, removed {len(removed_feats)} (threshold >= {args.var_threshold:.4g} variance)")
	vlog(f"variance filtering took {time.perf_counter()-t3:.2f}s; new shape={X_vf.shape}")
	if X_vf.shape[1] == 0:
		print("ERROR: All features removed by variance thresholding.", file=sys.stderr)
		return 4

	# Standardize
	t4 = time.perf_counter()
	scaler = StandardScaler(with_mean=True, with_std=True)
	X_std = scaler.fit_transform(X_vf)
	vlog(f"standardization took {time.perf_counter()-t4:.2f}s; matrix shape={X_std.shape}")

	# Compute pairwise distances and linkage once
	t5 = time.perf_counter()
	D = pairwise_feature_distance(X_std)
	Z = compute_linkage_matrix(D, linkage_method=args.linkage)
	n_features = X_std.shape[1]
	Kmax_eff = args.kmax if args.kmax is not None else max(1, min(100, n_features - 1))

	# Select K based on requested method
	if args.k_method == "gap":
		k_star, Z_used, Wk_arr, gap_arr, sk_arr, ks_arr = choose_k_via_gap_statistic(
			X_std, linkage_method=args.linkage, B=args.gap_B, Kmax=Kmax_eff, random_state=42, debug=args.verbose
		)
		# Prefer Z from the direct computation if available
		if Z_used is not None:
			Z = Z_used
			
		print(f"Clustering features: n_features={n_features}, k_method=gap, K={k_star}, linkage={args.linkage}")
		n_clusters = int(k_star)
		vlog(f"k selection (gap) took {time.perf_counter()-t5:.2f}s")
	elif args.k_method == "silhouette":
		k_star, ks_sil, sil_scores = choose_k_via_silhouette(D, Z, linkage_method=args.linkage, Kmax=Kmax_eff, debug=args.verbose)
		print(f"Clustering features: n_features={n_features}, k_method=silhouette, K={k_star}, linkage={args.linkage}")
		n_clusters = int(k_star)
		vlog(f"k selection (silhouette) took {time.perf_counter()-t5:.2f}s")
		# Save silhouette diagnostics
		if not args.no_dendrogram:
			base = stem_from_path(args.input)
			out_sil_csv = os.path.join(args.outdir, f"{base}_hc_silhouette_scores.csv")
			pd.DataFrame({"k": ks_sil, "silhouette": sil_scores}).to_csv(out_sil_csv, index=False)
			out_sil_png = os.path.join(args.outdir, f"{base}_hc_silhouette_curve.png")
			out_sil_pdf = os.path.join(args.outdir, f"{base}_hc_silhouette_curve.pdf")
			_save_silhouette_plot(out_sil_png, out_sil_pdf, ks_sil, sil_scores, n_clusters)
	else:
		# elbow (default)
		k_star, ks_elbow, Wk_elbow = choose_k_via_elbow(D, Z, linkage_method=args.linkage, Kmax=Kmax_eff, debug=args.verbose)
		print(f"Clustering features: n_features={n_features}, k_method=elbow, K={k_star}, linkage={args.linkage}")
		n_clusters = int(k_star)
		vlog(f"k selection (elbow) took {time.perf_counter()-t5:.2f}s")
		# Save elbow diagnostics
		if not args.no_dendrogram:
			base = stem_from_path(args.input)
			out_elbow_csv = os.path.join(args.outdir, f"{base}_hc_elbow_values.csv")
			pd.DataFrame({"k": ks_elbow, "Wk": Wk_elbow, "logWk": np.log(np.maximum(Wk_elbow, 1e-12))}).to_csv(out_elbow_csv, index=False)
			out_elbow_png = os.path.join(args.outdir, f"{base}_hc_elbow_curve.png")
			out_elbow_pdf = os.path.join(args.outdir, f"{base}_hc_elbow_curve.pdf")
			_save_elbow_plot(out_elbow_png, out_elbow_pdf, ks_elbow, Wk_elbow, n_clusters, log_scale=True)

	# Labels using chosen K
	t6 = time.perf_counter()
	if Z is None:
		labels, Z = hierarchical_cluster_labels(D, n_clusters=n_clusters, linkage_method=args.linkage, debug=args.verbose)
	else:
		labels = fcluster(Z, t=n_clusters, criterion="maxclust") - 1
	vlog(f"hierarchical clustering produced {len(np.unique(labels))} clusters in {time.perf_counter()-t6:.2f}s")

	# Save dendrogram before medoid selection
	if not args.no_dendrogram and Z is not None:
		base = stem_from_path(args.input)
		out_dendro_png = os.path.join(args.outdir, f"{base}_hc_dendrogram.png")
		out_dendro_pdf = os.path.join(args.outdir, f"{base}_hc_dendrogram.pdf")
		t7p = time.perf_counter()
		save_dendrogram(Z, out_png=out_dendro_png, out_pdf=out_dendro_pdf, debug=args.verbose)
		vlog(f"saved dendrogram in {time.perf_counter()-t7p:.2f}s -> {out_dendro_png}")

		# Save selection diagnostics if gap method was used earlier in this run
		if 'gap_arr' in locals() and isinstance(gap_arr, np.ndarray) and gap_arr.size > 0 and np.isfinite(gap_arr).all():
			out_gap_csv = os.path.join(args.outdir, f"{base}_hc_gap_values.csv")
			gap_df = pd.DataFrame({
				"k": ks_arr,
				"Wk": Wk_arr if np.ndim(Wk_arr) else [Wk_arr],
				"logWk": np.log(np.maximum(Wk_arr, 1e-12)) if np.ndim(Wk_arr) else [float(np.log(max(float(Wk_arr), 1e-12)))],
				"gap": gap_arr,
				"sk": sk_arr,
			})
			gap_df.to_csv(out_gap_csv, index=False)
			out_gap_png = os.path.join(args.outdir, f"{base}_hc_gap_curve.png")
			out_gap_pdf = os.path.join(args.outdir, f"{base}_hc_gap_curve.pdf")
			_save_gap_plot(out_gap_png, out_gap_pdf, ks_arr, gap_arr, sk_arr, n_clusters)

	# Medoid selection per cluster
	t7 = time.perf_counter()
	medoid_idx = medoids_from_clusters(D, labels)
	medoid_idx_sorted = sorted(medoid_idx)  # stable order by feature index
	selected_feats = [kept_feats[i] for i in medoid_idx_sorted]
	vlog(f"selected {len(selected_feats)} medoid features in {time.perf_counter()-t7:.2f}s")

	# Outputs
	base = stem_from_path(args.input)
	out_filtered_csv = os.path.join(args.outdir, f"{base}_original_full_filtered.csv")
	out_lowvar_kept = os.path.join(args.outdir, f"{base}_lowvar_kept_features.txt")
	out_lowvar_removed = os.path.join(args.outdir, f"{base}_lowvar_removed_features.txt")
	out_cluster_labels = os.path.join(args.outdir, f"{base}_feature_clusters.csv")
	out_selected_list = os.path.join(args.outdir, f"{base}_hc_medoid_selected_features.txt")
	out_selected_matrix = os.path.join(args.outdir, f"{base}_hc_medoid_selected_original_features.csv")

	t8 = time.perf_counter()
	# Filtered matrix
	filtered_df = pd.concat([sub_num[[ID_COL]].reset_index(drop=True), pd.DataFrame(X_vf, columns=kept_feats)], axis=1)
	filtered_df.to_csv(out_filtered_csv, index=False)

	# Feature lists
	with open(out_lowvar_kept, "w") as f:
		for feat in kept_feats:
			f.write(f"{feat}\n")
	with open(out_lowvar_removed, "w") as f:
		for feat in removed_feats:
			f.write(f"{feat}\n")

	# Cluster labels per feature
	cl_df = pd.DataFrame({"feature": kept_feats, "cluster": labels.astype(int)})
	cl_df.to_csv(out_cluster_labels, index=False)

	# Selected medoid list and matrix
	with open(out_selected_list, "w") as f:
		for feat in selected_feats:
			f.write(f"{feat}\n")
	selected_df = pd.concat([
		sub_num[[ID_COL]].reset_index(drop=True),
		sub_num[selected_feats].reset_index(drop=True),
	], axis=1)
	selected_df.to_csv(out_selected_matrix, index=False)
	vlog(f"wrote outputs in {time.perf_counter()-t8:.2f}s")

	print("Done.")
	print(f"- Filtered matrix:    {out_filtered_csv}")
	print(f"- Kept features:      {out_lowvar_kept}")
	print(f"- Removed features:   {out_lowvar_removed}")
	print(f"- Cluster labels:     {out_cluster_labels}")
	print(f"- Selected features:  {out_selected_list}")
	print(f"- Selected matrix:    {out_selected_matrix}")
	return 0


if __name__ == "__main__":
	sys.exit(main())

