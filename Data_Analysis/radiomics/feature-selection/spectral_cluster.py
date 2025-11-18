#!/usr/bin/env python3
"""
Spectral clustering pipeline for FMCIB radiomic features (feature clustering, pancancer).

Goal: Perform pancancer spectral clustering on radiomic features using gap-statistic-based k selection
to find optimal k clusters, visualize clustering results, and select representative features.

High-level Steps:
 1) Read one or more cohort feature matrices; per cohort: keep original/full rows and pred_* columns
 2) Combine all samples across cohorts into one matrix (align by common feature set)
 3) Standardize features across all combined samples
 4) Variance filter (drop bottom X%)
 5) Build similarity (affinity) matrix on features using RBF kernel
 6) Select k via silhouette across candidate k values using spectral embedding + KMeans
 7) Run spectral clustering with optimal k (clustering features)
 8) Visualize 2D feature embedding
 9) Select one medoid (representative feature) per cluster and export medoid feature matrix

Outputs (written to outdir):
 - <base>_combined_original_full_features_only.csv : raw combined SampleID + pred_* (common features)
 - <base>_kept_features_after_variance.txt         : features kept after variance filter
 - <base>_standardized_matrix.csv                  : SampleID + kept features (z-scored)
 - <base>_gap_statistic.csv                        : Gap(k) table
 - <base>_gap_statistic.png                        : Gap(k) plot
 - <base>_feature_clusters.csv                     : feature -> cluster assignment
 - <base>_feature_embedding_2d.csv                 : 2D spectral embedding coordinates per feature
 - <base>_feature_embedding_2d.png                 : 2D scatter colored by cluster
 - <base>_feature_tsne_2d.csv                      : 2D t-SNE embedding coordinates per feature
 - <base>_feature_tsne_2d.png                      : 2D t-SNE scatter colored by cluster
 - <base>_feature_tsne_3d.csv                      : 3D t-SNE embedding coordinates per feature (if --tsne-3d)
 - <base>_feature_tsne_3d.png                      : 3D t-SNE scatter colored by cluster (if --tsne-3d)
 - <base>_cluster_medoids.csv                      : per-cluster representative feature (medoid)
 - <base>_medoid_feature_matrix.csv                : SampleID + standardized medoid features (one per cluster)
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Tuple

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import SpectralClustering, KMeans
from sklearn.manifold import SpectralEmbedding, TSNE
from sklearn.metrics.pairwise import rbf_kernel, pairwise_distances
from scipy.sparse import csr_matrix

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except Exception:
    MATPLOTLIB_AVAILABLE = False

# Schema
ID_COL = "SampleID"
PERM_COL = "readii_Permutation"
REGION_COL = "readii_Region"
FEATURE_PREFIX = "pred_"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Pancancer spectral clustering of FMCIB features (feature-level). Supports RBF, Euclidean, or Spearman affinity. Uses Gap statistic to choose k.")
    p.add_argument("--inputs", "-i", nargs="+", required=True, help="Paths to input FMCIB features CSVs (one or more cohorts)")
    p.add_argument("--outdir", "-o", required=True, help="Output directory")
    p.add_argument("--name", default="pancancer", help="Base name for outputs [default: %(default)s]")
    p.add_argument("--feature-prefix", default=FEATURE_PREFIX, help="Prefix for feature columns [default: %(default)s]")
    p.add_argument("--var-min-abs", type=float, default=0.15, help="Remove features with variance < this absolute threshold [default: %(default)s]")
    p.add_argument("--affinity", choices=["rbf", "euclidean", "spearman"], default="rbf", help="Affinity type for feature similarity graph [default: %(default)s]")
    p.add_argument(
        "--rbf-gamma",
        default="auto",
        help=(
            "Gamma for RBF kernel over features."
            " Accepts a float, 'auto' (1/n_features post-variance), or 'median' (1/(2*median(d)^2))."
        ),
    )
    p.add_argument("--rbf-threshold", type=float, default=0.0, help="Zero out small affinities below this value to sparsify the graph (applies to any affinity) [default: %(default)s]")
    p.add_argument("--min-clusters", type=int, default=2, help="Minimum number of clusters to consider [default: %(default)s]")
    p.add_argument("--max-clusters", type=int, default=None, help="Maximum number of clusters to consider. Default: min(100, n_features-1)")
    p.add_argument("--gap-b", type=int, default=10, help="Number of reference datasets (B) for Gap statistic [default: %(default)s]")
    # t-SNE options
    p.add_argument("--tsne-perplexity", type=float, default=30.0, help="t-SNE perplexity (adjusted if >= n_features) [default: %(default)s]")
    p.add_argument("--tsne-max-iter", type=int, default=1000, help="t-SNE maximum iterations [default: %(default)s]")
    p.add_argument("--tsne-pca-dims", type=int, default=0, help="Optional PCA precompression dimensions (0 = skip) [default: %(default)s]")
    p.add_argument("--no-tsne", action="store_true", help="Disable t-SNE embedding output")
    p.add_argument("--tsne-3d", action="store_true", help="Also compute a 3D t-SNE embedding (CSV + PNG)")
    p.add_argument("--seed", type=int, default=10, help="Random seed for clustering [default: %(default)s]")
    return p.parse_args()


def stem_from_path(path: str) -> str:
    name = os.path.basename(path)
    return name[:-4] if name.lower().endswith(".csv") else name


def read_and_filter(input_csv: str, feature_prefix: str) -> Tuple[pd.DataFrame, List[str]]:
    df = pd.read_csv(input_csv)
    for c in (ID_COL, PERM_COL, REGION_COL):
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")
    perm = df[PERM_COL].astype(str).str.strip().str.lower()
    region = df[REGION_COL].astype(str).str.strip().str.lower()
    sub = df.loc[(perm == "original") & (region == "full")].copy()
    if sub.empty:
        raise ValueError("No rows after filtering for original/full.")
    feat_cols = [c for c in sub.columns if c.startswith(feature_prefix)]
    if not feat_cols:
        raise ValueError(f"No feature columns starting with '{feature_prefix}' after filtering.")
    keep_cols = ([ID_COL] if ID_COL in sub.columns else []) + feat_cols
    sub = sub.loc[:, keep_cols]
    # Numeric coercion and inf->NaN
    sub[feat_cols] = sub[feat_cols].apply(pd.to_numeric, errors="coerce")
    sub[feat_cols] = sub[feat_cols].replace([np.inf, -np.inf], np.nan)
    return sub, feat_cols


def variance_filter(X: pd.DataFrame, min_abs: float) -> Tuple[pd.DataFrame, List[str], float, int]:
    """Filter features by absolute variance threshold.

    Returns (X_filtered, kept_feature_names, cutoff, removed_count).
    """
    var_vals = X.var(axis=0, ddof=0).to_numpy(dtype=float)
    cutoff = float(min_abs)
    keep_mask = var_vals >= cutoff
    cols = X.columns.tolist()
    kept = [c for c, k in zip(cols, keep_mask) if bool(k)]
    if not kept:
        # Fallback: keep top-variance feature
        top_idx = int(np.argmax(var_vals))
        kept = [cols[top_idx]]
    removed = int(len(cols) - len(kept))
    return X[kept].copy(), kept, cutoff, removed


def auto_neighbors(n_items: int) -> int:
    if n_items <= 2:
        return max(1, n_items - 1)
    k = int(round(np.sqrt(n_items)))
    k = max(2, min(k, n_items - 1))
    return k


def main() -> int:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Read and filter all cohorts
    subs: List[pd.DataFrame] = []
    feat_lists: List[List[str]] = []
    for path in args.inputs:
        print(f"Reading: {path}")
        try:
            sub_i, feats_i = read_and_filter(path, args.feature_prefix)
        except Exception as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 2
        subs.append(sub_i)
        feat_lists.append(feats_i)

    # Align by common feature set across cohorts
    common_feats = sorted(set(feat_lists[0]).intersection(*map(set, feat_lists[1:]))) if len(feat_lists) > 1 else feat_lists[0]
    if len(common_feats) < 2:
        print(f"ERROR: Too few common features after alignment across cohorts: {len(common_feats)}", file=sys.stderr)
        return 3
    # Combine samples (stack rows)
    combined = pd.concat([df[[ID_COL] + common_feats].copy() for df in subs], axis=0, ignore_index=True)

    base = args.name.strip() or "pancancer"
    # Save combined raw matrix
    out_raw = os.path.join(args.outdir, f"{base}_combined_original_full_features_only.csv")
    combined.to_csv(out_raw, index=False)
    print(f"Combined matrix: {out_raw} (rows={combined.shape[0]}, features={len(common_feats)})")

    # Impute missing with median per feature before variance and clustering
    X = combined[common_feats].copy()
    X = X.apply(lambda s: s.fillna(s.median()), axis=0)

    # Variance filter (use pre-standardized values across samples)
    X_var, kept_after_var, cutoff, removed_ct = variance_filter(X, args.var_min_abs)
    with open(os.path.join(args.outdir, f"{base}_kept_features_after_variance.txt"), "w") as f:
        for nm in kept_after_var:
            f.write(f"{nm}\n")
    total_feats = len(X.columns)
    print(
        f"Variance filter: removed {removed_ct} of {total_feats} features with variance < {cutoff:.4g}; kept {len(kept_after_var)} features"
    )

    # Standardize (z-score) kept features across all combined samples
    scaler = StandardScaler(with_mean=True, with_std=True)
    Xz = scaler.fit_transform(X_var.values)  # shape: (n_samples_total, n_features_kept)
    out_std = os.path.join(args.outdir, f"{base}_standardized_matrix.csv")
    pd.concat([combined[[ID_COL]].reset_index(drop=True), pd.DataFrame(Xz, columns=kept_after_var)], axis=1).to_csv(out_std, index=False)

    # Build feature affinity (post-variance filtering)
    F = Xz.T  # shape: (n_features, n_samples)
    n_features = F.shape[0]
    if args.affinity == "rbf":
        gamma_val = None
        if isinstance(args.rbf_gamma, str):
            mode = args.rbf_gamma.lower().strip()
            if mode == "auto":
                gamma_val = 1.0 / max(1, n_features)
            elif mode == "median":
                # Median heuristic: gamma = 1 / (2 * median(d)^2), d = pairwise Euclidean distance between features
                D = pairwise_distances(F, metric="euclidean")
                if D.shape[0] >= 2:
                    tri = D[np.triu_indices(D.shape[0], k=1)]
                    med = float(np.median(tri)) if tri.size > 0 else float(np.median(D))
                else:
                    med = float(np.median(D))
                if med <= 1e-12 or not np.isfinite(med):
                    gamma_val = 1.0 / max(1, n_features)
                    print(
                        "Warning: median distance non-positive/non-finite; falling back to auto gamma (1/n_features).",
                        file=sys.stderr,
                    )
                else:
                    gamma_val = 1.0 / (2.0 * (med ** 2))
            else:
                try:
                    gamma_val = float(args.rbf_gamma)
                except Exception:
                    gamma_val = 1.0 / max(1, n_features)
        else:
            try:
                gamma_val = float(args.rbf_gamma)
            except Exception:
                gamma_val = 1.0 / max(1, n_features)
        print(f"Affinity: RBF (gamma={gamma_val:.6g})")
        W_dense = rbf_kernel(F, gamma=gamma_val)
    elif args.affinity == "euclidean":
        # Euclidean distance -> similarity conversion.
        # Compute pairwise distances between feature vectors (each feature across samples).
        D = pairwise_distances(F, metric="euclidean")  # shape (n_features, n_features)
        # Convert distances to similarities. We use S = 1 / (1 + D) (bounded in (0,1]).
        # Alternative scalings can be explored later (e.g., max(D)-D or exp(-D/scale)).
        W_dense = 1.0 / (1.0 + D)
        print("Affinity: Euclidean distance converted via S = 1/(1 + d).")
    else:
        # Spearman correlation between feature vectors (each feature across samples)
        Xz_df = pd.DataFrame(F.T, columns=kept_after_var)  # rows: samples, cols: features
        W_corr = Xz_df.corr(method="spearman")  # (n_features, n_features)
        W_dense = np.abs(W_corr.values)  # non-negative similarity
        # Replace any NaNs resulting from constant vectors
        W_dense = np.nan_to_num(W_dense, nan=0.0, posinf=1.0, neginf=0.0)
        print("Affinity: Spearman |rho| (non-negative).")
    thr = float(args.rbf_threshold)
    if thr > 0.0:
        before = np.count_nonzero(W_dense)
        W_dense[W_dense < thr] = 0.0
        after = np.count_nonzero(W_dense)
        print(f"Thresholded affinities < {thr} -> zero (nonzeros {before} -> {after}).")
    W = csr_matrix(W_dense)

    # Select k via Gap statistic over candidate range
    max_clusters = int(args.max_clusters) if args.max_clusters is not None else int(min(50, n_features - 1))
    min_clusters = int(max(2, min(args.min_clusters, max_clusters)))
    # Gap statistic on spectral embedding, using KMeans inertia as within-dispersion measure
    B = int(max(1, args.gap_b))
    rng = np.random.RandomState(int(args.seed))
    rows = []
    gap_vals = []
    se_vals = []
    for k in range(min_clusters, max_clusters + 1):
        try:
            embed_k = SpectralEmbedding(n_components=k, affinity="precomputed", random_state=int(args.seed))
            E_k = embed_k.fit_transform(W_dense)  # (n_features, k)
            km = KMeans(n_clusters=k, n_init=10, random_state=int(args.seed))
            km.fit(E_k)
            Wk = float(km.inertia_)
            # Reference datasets
            mins = E_k.min(axis=0)
            maxs = E_k.max(axis=0)
            logWkb = []
            for b in range(B):
                U = rng.uniform(low=mins, high=maxs, size=E_k.shape)
                km_b = KMeans(n_clusters=k, n_init=10, random_state=int(args.seed))
                km_b.fit(U)
                logWkb.append(np.log(max(km_b.inertia_, 1e-12)))
            logWkb = np.array(logWkb, dtype=float)
            gap_k = float(np.mean(logWkb) - np.log(max(Wk, 1e-12)))
            sk = float(np.std(logWkb, ddof=1) * np.sqrt(1.0 + 1.0 / B)) if B > 1 else 0.0
        except Exception as e:
            print(f"Warning: k={k} gap statistic failed: {e}", file=sys.stderr)
            Wk, gap_k, sk = np.nan, np.nan, np.nan
        rows.append({"k": k, "gap": gap_k, "sk": sk})
        gap_vals.append(gap_k)
        se_vals.append(sk)
    gap_df = pd.DataFrame(rows)
    gap_csv = os.path.join(args.outdir, f"{base}_gap_statistic.csv")
    gap_df.to_csv(gap_csv, index=False)
    # Choose k via Tibshirani rule: smallest k s.t. Gap(k) >= Gap(k+1) - s_{k+1}
    k_candidates = list(range(min_clusters, max_clusters))  # up to max-1
    chosen = None
    for idx, k in enumerate(k_candidates):
        gk = gap_vals[idx]
        gk1 = gap_vals[idx + 1]
        sk1 = se_vals[idx + 1]
        if np.isfinite(gk) and np.isfinite(gk1) and np.isfinite(sk1) and (gk >= gk1 - sk1):
            chosen = k
            break
    if chosen is None:
        # Fallback: argmax gap
        finite_idx = [i for i, g in enumerate(gap_vals) if np.isfinite(g)]
        chosen = (min_clusters + int(finite_idx[np.argmax([gap_vals[i] for i in finite_idx])])) if finite_idx else min_clusters
    k_est = int(chosen)
    print(f"Selected k via gap statistic: k={k_est}")
    if MATPLOTLIB_AVAILABLE:
        gap_png = os.path.join(args.outdir, f"{base}_gap_statistic.png")
        plt.figure(figsize=(6,4), dpi=150)
        plt.plot(gap_df["k"], gap_df["gap"], marker='o', label='Gap(k)')
        plt.fill_between(gap_df["k"], gap_df["gap"]-gap_df["sk"], gap_df["gap"]+gap_df["sk"], alpha=0.2, label='±s(k)')
        plt.axvline(k_est, color='tab:red', ls='--', lw=1, label=f"k={k_est}")
        plt.xlabel("k")
        plt.ylabel("Gap")
        plt.title("Gap statistic (spectral embedding)")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(gap_png)
        plt.close()

    # Spectral clustering on FEATURES using precomputed RBF affinity
    spec = SpectralClustering(
        n_clusters=int(k_est),
        affinity="precomputed",
        assign_labels="kmeans",
        n_init=10,
        random_state=int(args.seed),
    )
    feature_labels = spec.fit_predict(W_dense)  # length = n_features

    # 2D spectral embedding for FEATURES (precomputed affinity)
    if MATPLOTLIB_AVAILABLE:
        embedder = SpectralEmbedding(
            n_components=2,
            affinity="precomputed",
            random_state=int(args.seed),
        )
        Y = embedder.fit_transform(W_dense)
        emb_csv = os.path.join(args.outdir, f"{base}_feature_embedding_2d.csv")
        pd.DataFrame({"feature": kept_after_var, "x": Y[:, 0], "y": Y[:, 1], "cluster": feature_labels}).to_csv(emb_csv, index=False)
        emb_png = os.path.join(args.outdir, f"{base}_feature_embedding_2d.png")
        plt.figure(figsize=(6, 5), dpi=150)
        plt.scatter(Y[:, 0], Y[:, 1], c=feature_labels, s=10, cmap="tab10", alpha=0.9)
        plt.xlabel("Spectral component 1")
        plt.ylabel("Spectral component 2")
        plt.title(f"Feature spectral embedding ({args.affinity} affinity)")
        plt.grid(alpha=0.2)
        plt.tight_layout()
        plt.savefig(emb_png)
        plt.close()
        # t-SNE embedding (features as observations) unless disabled
        if not args.no_tsne:
            # t-SNE input: original standardized feature vectors F (n_features x n_samples)
            tsne_input = F.copy()
            if args.tsne_pca_dims and args.tsne_pca_dims > 0:
                try:
                    from sklearn.decomposition import PCA
                    pca_dims = int(args.tsne_pca_dims)
                    pca_dims = max(2, min(pca_dims, tsne_input.shape[1]))
                    P = PCA(n_components=pca_dims, random_state=int(args.seed))
                    tsne_input = P.fit_transform(tsne_input)
                    print(f"t-SNE: PCA precompression to {pca_dims} dims")
                except Exception as e:
                    print(f"Warning: PCA precompression failed: {e}", file=sys.stderr)
            perplexity = float(args.tsne_perplexity)
            # Perplexity must be < n_features; adjust if needed.
            if perplexity >= tsne_input.shape[0]:
                new_perp = max(5.0, tsne_input.shape[0] // 3)
                print(f"t-SNE: adjusting perplexity {perplexity} -> {new_perp} (must be < n_features)")
                perplexity = new_perp
            try:
                tsne = TSNE(
                    n_components=2,
                    perplexity=perplexity,
                    init="pca",
                    learning_rate="auto",
                    max_iter=int(args.tsne_max_iter),
                    random_state=int(args.seed),
                )
                Z = tsne.fit_transform(tsne_input)
                tsne_csv = os.path.join(args.outdir, f"{base}_feature_tsne_2d.csv")
                pd.DataFrame({"feature": kept_after_var, "x": Z[:, 0], "y": Z[:, 1], "cluster": feature_labels}).to_csv(tsne_csv, index=False)
                tsne_png = os.path.join(args.outdir, f"{base}_feature_tsne_2d.png")
                plt.figure(figsize=(6,5), dpi=150)
                plt.scatter(Z[:,0], Z[:,1], c=feature_labels, s=10, cmap='tab10', alpha=0.9)
                plt.xlabel("t-SNE 1")
                plt.ylabel("t-SNE 2")
                plt.title("Feature t-SNE embedding")
                plt.grid(alpha=0.2)
                plt.tight_layout()
                plt.savefig(tsne_png)
                plt.close()
                print(f"t-SNE embedding written: {tsne_csv}")
                # Optional 3D t-SNE
                if args.tsne_3d:
                    tsne3 = TSNE(
                        n_components=3,
                        perplexity=perplexity,
                        init="pca",
                        learning_rate="auto",
                        max_iter=int(args.tsne_max_iter),
                        random_state=int(args.seed),
                    )
                    Z3 = tsne3.fit_transform(tsne_input)
                    tsne3_csv = os.path.join(args.outdir, f"{base}_feature_tsne_3d.csv")
                    pd.DataFrame({"feature": kept_after_var, "x": Z3[:, 0], "y": Z3[:, 1], "z": Z3[:, 2], "cluster": feature_labels}).to_csv(tsne3_csv, index=False)
                    tsne3_png = os.path.join(args.outdir, f"{base}_feature_tsne_3d.png")
                    fig = plt.figure(figsize=(6,5), dpi=150)
                    try:
                        ax = fig.add_subplot(111, projection='3d')
                    except Exception:
                        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
                        ax = fig.add_subplot(111, projection='3d')
                    ax.scatter(Z3[:,0], Z3[:,1], Z3[:,2], c=feature_labels, s=10, cmap='tab10', alpha=0.9)
                    ax.set_xlabel("t-SNE 1")
                    ax.set_ylabel("t-SNE 2")
                    ax.set_zlabel("t-SNE 3")
                    ax.set_title("Feature t-SNE embedding (3D)")
                    fig.tight_layout()
                    fig.savefig(tsne3_png)
                    plt.close(fig)
                    print(f"t-SNE 3D embedding written: {tsne3_csv}")
            except Exception as e:
                print(f"Warning: t-SNE failed: {e}", file=sys.stderr)

    # Save feature cluster assignments
    out_clusters = os.path.join(args.outdir, f"{base}_feature_clusters.csv")
    pd.DataFrame({"feature": kept_after_var, "cluster": feature_labels}).to_csv(out_clusters, index=False)

    # Select medoid per cluster: choose feature with highest within-cluster affinity sum
    medoid_rows = []
    labels = np.asarray(feature_labels)
    for c in np.unique(labels):
        idx = np.where(labels == c)[0]
        if idx.size == 0:
            continue
        W_sub = W_dense[np.ix_(idx, idx)]
        degrees = np.sum(W_sub, axis=1)
        best_local = int(np.argmax(degrees))
        feat_name = kept_after_var[idx[best_local]]
        medoid_rows.append({"cluster": int(c), "feature": feat_name, "cluster_size": int(idx.size), "degree_sum": float(degrees[best_local])})
    medoids_csv = os.path.join(args.outdir, f"{base}_cluster_medoids.csv")
    pd.DataFrame(medoid_rows).sort_values(["cluster"], ascending=True).to_csv(medoids_csv, index=False)

    # Build and save medoid feature matrix (standardized values)
    medoid_features = [r["feature"] for r in sorted(medoid_rows, key=lambda d: d["cluster"]) ]
    Xz_df = pd.DataFrame(Xz, columns=kept_after_var)
    medoid_mat = pd.concat([combined[[ID_COL]].reset_index(drop=True), Xz_df[medoid_features]], axis=1)
    medoid_mat_csv = os.path.join(args.outdir, f"{base}_medoid_feature_matrix.csv")
    medoid_mat.to_csv(medoid_mat_csv, index=False)

    print("Done.")
    print(f"- Combined filtered matrix: {out_raw}")
    print(f"- Kept after variance:      {base}_kept_features_after_variance.txt")
    print(f"- Standardized matrix:      {out_std}")
    print(f"- Gap statistic:            {base}_gap_statistic.csv{' / ' + base + '_gap_statistic.png' if MATPLOTLIB_AVAILABLE else ''}")
    print(f"- Feature clusters:         {out_clusters}")
    if MATPLOTLIB_AVAILABLE:
        print(f"- Feature embedding (2D):   {base}_feature_embedding_2d.csv / {base}_feature_embedding_2d.png")
        if not args.no_tsne:
            print(f"- Feature t-SNE (2D):       {base}_feature_tsne_2d.csv / {base}_feature_tsne_2d.png")
            if args.tsne_3d:
                print(f"- Feature t-SNE (3D):       {base}_feature_tsne_3d.csv / {base}_feature_tsne_3d.png")
    print(f"- Cluster medoids:          {medoids_csv}")
    print(f"- Medoid feature matrix:    {medoid_mat_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
