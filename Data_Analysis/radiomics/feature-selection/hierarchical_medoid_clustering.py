#!/usr/bin/env python3
"""
Hierarchical medoid clustering of radiomic (or other) features using absolute Spearman correlation.

Pipeline:
 1. Read feature matrix (CSV). Assumes rows = samples, columns = features. Optionally an ID column.
 2. Select feature columns (all numeric columns excluding ID by default, or those with a given prefix).
 3. Remove near-zero variance features (variance <= threshold OR <2 unique values).
 4. Compute Spearman correlation matrix between remaining features.
    - We use rank-transform + Pearson on ranks (Spearman equivalent) for speed and NaN robustness.
 5. Convert to distance matrix: D = 1 - |rho| (clusters ignore correlation sign, emphasize strength).
 6. Hierarchical clustering on D with chosen linkage (average by default).
 7. Cut tree using a distance threshold (default: 0.6). fcluster with criterion='distance'.
    - Features whose inter-cluster linkage distance <= threshold are merged.
    - Distance threshold implies grouping features whose |rho| ≥ 1 - threshold (e.g., threshold=0.6 => |rho| ≥ 0.4).
 8. For each resulting cluster, choose a medoid feature: argmin_i mean_j D[i,j] within cluster.
 9. Write outputs:
    * filtered feature matrix (after NZV removal)
    * distance threshold used
    * cluster membership & medoids CSV
    * list of selected medoid features
    * matrix with only medoid features (one representative per cluster)
    * optional distance & correlation matrices (if --save-matrices)

Usage example:
  python hierarchical_medoid_clustering.py \
      --input data/procdata/radiomic_features/fmcib/original/combined/fmcib_original_features_all_datasets.csv \
      --outdir data/procdata/radiomic_features/fmcib/features/hierarch_medoid_clusters \
      --id-col SampleID --feature-prefix pred_ --nzv-threshold 1e-6 --dist-threshold 0.6
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import List, Tuple
import numpy as np
import pandas as pd

try:
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except Exception:
    MATPLOTLIB_AVAILABLE = False


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Hierarchical clustering of features via absolute Spearman distance and medoid selection.")
    p.add_argument("--input", "-i", required=True, help="Input CSV (rows=samples, columns=features).")
    p.add_argument("--outdir", "-o", required=True, help="Output directory.")
    p.add_argument("--id-col", default="SampleID", help="Optional ID column to exclude from feature set (default: SampleID).")
    p.add_argument("--feature-prefix", default=None, help="If provided, only columns starting with this prefix are treated as features.")
    p.add_argument("--nzv-threshold", type=float, default=1e-8, help="Near-zero variance threshold; remove features with variance <= threshold (default: 1e-8).")
    p.add_argument("--dist-threshold", type=float, default=0.6, help="Distance threshold for cutting tree (criterion='distance') (default: 0.6).")
    p.add_argument("--linkage", choices=["average", "complete", "single", "ward"], default="average", help="Linkage method (default: average). Ward only meaningful for Euclidean distances; used here with Spearman distances cautiously.")
    p.add_argument("--save-matrices", action="store_true", help="Save correlation and distance matrices as CSVs.")
    p.add_argument("--no-dendrogram", action="store_true", help="Skip dendrogram plot.")
    p.add_argument("--exclude-singletons", action="store_true", help="Exclude singleton clusters from medoid selection (leave medoid_feature blank for size==1 clusters).")
    p.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return p.parse_args()


def ensure_outdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def vlog(verbose: bool, msg: str) -> None:
    if verbose:
        print(f"[DEBUG] {msg}")


def load_feature_matrix(path: str) -> pd.DataFrame:
    return pd.read_csv(path)


def select_feature_columns(df: pd.DataFrame, id_col: str, prefix: str | None) -> List[str]:
    cols: List[str] = []
    for c in df.columns:
        if c == id_col:
            continue
        if prefix is not None and not c.startswith(prefix):
            continue
        # Require numeric
        try:
            _ = pd.to_numeric(df[c], errors="coerce")
        except Exception:
            continue
        cols.append(c)
    return cols


def remove_near_zero_variance(df: pd.DataFrame, feature_cols: List[str], threshold: float, verbose: bool) -> Tuple[pd.DataFrame, List[str], List[str]]:
    variances = {}
    keep = []
    drop = []
    for c in feature_cols:
        x = pd.to_numeric(df[c], errors="coerce")
        x = x[np.isfinite(x)]
        if x.size == 0:
            drop.append(c)
            continue
        var = float(np.var(x.values, ddof=0))
        variances[c] = var
        unique_ct = int(pd.Series(x).nunique())
        if unique_ct < 2 or var <= threshold:
            drop.append(c)
        else:
            keep.append(c)
    if verbose:
        if variances:
            all_vars = np.array(list(variances.values()))
            print(f"[DEBUG] NZV: threshold={threshold}, kept={len(keep)}, removed={len(drop)}, var[min/median/max]=[{all_vars.min():.3g}/{np.median(all_vars):.3g}/{all_vars.max():.3g}]")
    return df[keep].copy(), keep, drop


def spearman_abs_distance(X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute absolute Spearman correlation matrix and distance matrix D = 1 - |rho|.
    Returns (rho_abs, D)."""
    # Rank transform columns
    ranks = pd.DataFrame(X).rank(axis=0, method="average", na_option="keep").to_numpy(dtype=float)
    C = np.corrcoef(ranks, rowvar=False)
    C = np.clip(C, -1.0, 1.0)
    rho_abs = np.abs(C)
    D = 1.0 - rho_abs
    np.fill_diagonal(D, 0.0)
    D[D < 0] = 0.0
    return rho_abs, D


def hierarchical_clusters(D: np.ndarray, dist_threshold: float, linkage_method: str, verbose: bool) -> np.ndarray:
    if not SCIPY_AVAILABLE:
        raise RuntimeError("SciPy is required for hierarchical clustering in this script.")
    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method=linkage_method)
    labels = fcluster(Z, t=dist_threshold, criterion='distance') - 1  # zero-based
    vlog(verbose, f"Hierarchical clustering produced {len(np.unique(labels))} clusters (distance threshold={dist_threshold}).")
    return labels, Z


def medoids_from_clusters(D: np.ndarray, labels: np.ndarray, exclude_singletons: bool = False) -> dict:
    """Compute medoid indices per cluster.
    Returns a mapping {cluster_label(int): medoid_index(int)}.
    If exclude_singletons is True, clusters with size==1 are omitted from the mapping."""
    medoid_map: dict = {}
    for lab in sorted(np.unique(labels)):
        idx = np.where(labels == lab)[0]
        if idx.size == 1:
            if exclude_singletons:
                continue
            medoid_map[int(lab)] = int(idx[0])
            continue
        subD = D[np.ix_(idx, idx)]
        avg_dist = subD.mean(axis=1)
        medoid_map[int(lab)] = int(idx[int(np.argmin(avg_dist))])
    return medoid_map


def save_dendrogram(Z, out_png: str, out_pdf: str | None, verbose: bool, color_threshold: float | None = None) -> None:
    """Save dendrogram. If color_threshold is provided, SciPy will color clusters
    by that linkage distance (aligns with fcluster distance cutoff)."""
    if not MATPLOTLIB_AVAILABLE:
        vlog(verbose, "Matplotlib not available; skipping dendrogram plot.")
        return
    plt.figure(figsize=(10, 5))
    dendrogram(Z, no_labels=True, color_threshold=color_threshold)
    # Draw a horizontal line at the threshold for visual reference
    if color_threshold is not None:
        try:
            plt.axhline(y=color_threshold, color="red", linestyle="--", linewidth=1.0, label=f"threshold={color_threshold}")
            # Add legend only if axhline successful
            plt.legend(loc="upper right", frameon=True, fontsize=8)
        except Exception:
            pass
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    if out_pdf:
        plt.savefig(out_pdf)
    plt.close()


def main() -> int:
    args = parse_args()
    ensure_outdir(args.outdir)

    print(f"Reading: {args.input}")
    df_raw = load_feature_matrix(args.input)
    if df_raw.empty:
        print("ERROR: Empty input file.", file=sys.stderr)
        return 2

    feature_cols = select_feature_columns(df_raw, args.id_col, args.feature_prefix)
    if not feature_cols:
        print("ERROR: No numeric feature columns detected.", file=sys.stderr)
        return 3
    vlog(args.verbose, f"Detected {len(feature_cols)} candidate feature columns.")

    # NZV removal
    df_feat, kept_cols, removed_cols = remove_near_zero_variance(df_raw, feature_cols, args.nzv_threshold, args.verbose)
    if len(kept_cols) < 2:
        print("Warning: <2 features remain after near-zero variance filtering; medoid selection trivial.")
        # Output kept feature (if any) and exit early
        base = os.path.splitext(os.path.basename(args.input))[0]
        out_list = os.path.join(args.outdir, f"{base}_medoid_features.txt")
        with open(out_list, 'w') as f:
            for c in kept_cols:
                f.write(f"{c}\n")
        # Save filtered matrix
        filtered_path = os.path.join(args.outdir, f"{base}_filtered_matrix.csv")
        df_out = pd.concat([df_raw[[args.id_col]].reset_index(drop=True), df_raw[kept_cols].reset_index(drop=True)], axis=1) if args.id_col in df_raw.columns else df_raw[kept_cols]
        df_out.to_csv(filtered_path, index=False)
        print(f"Done. Saved trivial medoid feature list: {out_list}")
        return 0

    # Numeric matrix
    X = df_feat[kept_cols].apply(pd.to_numeric, errors='coerce').to_numpy(dtype=float)
    if not np.isfinite(X).all():
        print("ERROR: Non-finite values present after numeric conversion.", file=sys.stderr)
        return 4

    # Correlation & distance
    rho_abs, D = spearman_abs_distance(X)

    # Clustering
    labels, Z = hierarchical_clusters(D, args.dist_threshold, args.linkage, args.verbose)

    # Medoids (optionally excluding singleton clusters)
    cluster_medoids = medoids_from_clusters(D, labels, exclude_singletons=args.exclude_singletons)
    medoid_indices_sorted = sorted(cluster_medoids.values())
    medoid_features = [kept_cols[i] for i in medoid_indices_sorted]
    if args.exclude_singletons:
        singleton_total = int(np.sum([np.sum(labels == lab) == 1 for lab in np.unique(labels)]))
        excluded = int(len(np.unique(labels)) - len(cluster_medoids))
        vlog(args.verbose, f"Excluding {excluded} singleton clusters from medoid selection (total singletons detected: {singleton_total}).")

    # Outputs base
    base = os.path.splitext(os.path.basename(args.input))[0]
    filtered_csv = os.path.join(args.outdir, f"{base}_filtered_matrix.csv")
    removed_txt = os.path.join(args.outdir, f"{base}_nzv_removed_features.txt")
    kept_txt = os.path.join(args.outdir, f"{base}_nzv_kept_features.txt")
    clusters_csv = os.path.join(args.outdir, f"{base}_clusters_medoid.csv")
    medoids_txt = os.path.join(args.outdir, f"{base}_medoid_features.txt")
    medoid_matrix_csv = os.path.join(args.outdir, f"{base}_medoid_feature_matrix.csv")
    dendro_png = os.path.join(args.outdir, f"{base}_dendrogram.png")
    dendro_pdf = os.path.join(args.outdir, f"{base}_dendrogram.pdf")

    # Save filtered matrix
    if args.id_col in df_raw.columns:
        filtered_df = pd.concat([df_raw[[args.id_col]].reset_index(drop=True), df_raw[kept_cols].reset_index(drop=True)], axis=1)
    else:
        filtered_df = df_raw[kept_cols]
    filtered_df.to_csv(filtered_csv, index=False)

    # Save NZV lists
    with open(removed_txt, 'w') as f:
        for c in removed_cols:
            f.write(f"{c}\n")
    with open(kept_txt, 'w') as f:
        for c in kept_cols:
            f.write(f"{c}\n")

    # Cluster membership & medoids
    rows = []
    for lab in sorted(np.unique(labels)):
        idx = np.where(labels == lab)[0]
        members = [kept_cols[i] for i in idx]
        medoid = kept_cols[cluster_medoids[lab]] if lab in cluster_medoids else ''
        rows.append({
            'cluster_id': int(lab),
            'size': int(len(idx)),
            'medoid_feature': medoid,
            'members': ';'.join(members)
        })
    pd.DataFrame(rows).to_csv(clusters_csv, index=False)

    # Medoid list & matrix
    with open(medoids_txt, 'w') as f:
        for m in medoid_features:
            f.write(f"{m}\n")
    if medoid_features:
        medoid_df = pd.concat([df_raw[[args.id_col]].reset_index(drop=True), df_raw[medoid_features].reset_index(drop=True)], axis=1) if args.id_col in df_raw.columns else df_raw[medoid_features]
        medoid_df.to_csv(medoid_matrix_csv, index=False)
    else:
        # Create a skeletal file when no medoids were selected (e.g., all singletons excluded)
        if args.id_col in df_raw.columns:
            df_raw[[args.id_col]].to_csv(medoid_matrix_csv, index=False)
        else:
            pd.DataFrame().to_csv(medoid_matrix_csv, index=False)

    # Optional matrices
    if args.save_matrices:
        corr_csv = os.path.join(args.outdir, f"{base}_abs_spearman_correlation.csv")
        dist_csv = os.path.join(args.outdir, f"{base}_abs_spearman_distance.csv")
        pd.DataFrame(rho_abs, index=kept_cols, columns=kept_cols).to_csv(corr_csv)
        pd.DataFrame(D, index=kept_cols, columns=kept_cols).to_csv(dist_csv)

    # Dendrogram
    if not args.no_dendrogram and MATPLOTLIB_AVAILABLE and SCIPY_AVAILABLE:
        save_dendrogram(Z, dendro_png, dendro_pdf, args.verbose, color_threshold=args.dist_threshold)

    print("Done.")
    print(f"- Filtered feature matrix: {filtered_csv}")
    print(f"- NZV removed features:    {removed_txt}")
    print(f"- NZV kept features:       {kept_txt}")
    print(f"- Cluster membership CSV:  {clusters_csv}")
    print(f"- Medoid feature list:     {medoids_txt}")
    print(f"- Medoid feature matrix:   {medoid_matrix_csv}")
    if args.save_matrices:
        print(f"- Correlation matrix:      {corr_csv}")
        print(f"- Distance matrix:         {dist_csv}")
    if not args.no_dendrogram:
        if MATPLOTLIB_AVAILABLE and SCIPY_AVAILABLE:
            print(f"- Dendrogram plot:         {dendro_png}")
        else:
            print("(Dendrogram skipped: matplotlib or scipy missing)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
