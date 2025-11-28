#!/usr/bin/env python3
"""
Feature selection pipeline (VIF -> Standardize -> PCA) for FMCIB radiomic features.

Workflow:
 0) Input hygiene: filter to readii_Permutation=='original' AND readii_Region=='full'; keep SampleID + pred_* columns; coerce to numeric; replace ±inf with NaN; drop features with > max-missing-frac; impute remaining NaNs with median.
 1) VIF selection via iterative removal using the diagonal of the precision (inverse covariance) matrix of standardized features (Ledoit-Wolf shrinkage). Drop the highest-VIF features in batches until all remaining VIF <= threshold (default 5.0) OR the number of features is reduced to the target (default 500).
 2) Standardize selected features (z-score) and write the standardized matrix.
 3) PCA on standardized features; save explained variance ratios, cumulative variance, elbow estimate, plots, scores, and loadings.

Outputs (written to outdir):
 - <stem>_original_full_features_only.csv        : raw filtered SampleID + pred_*
 - <stem>_kept_features_after_missing.txt       : features kept after missing-value filter
 - <stem>_vif_selected_features.txt             : the selected feature names after VIF selection
 - <stem>_vif_selected_matrix.csv               : SampleID + selected (raw, pre-standardization)
 - <stem>_standardized_matrix.csv               : SampleID + selected (z-scored)
 - <stem>_pca_explained_variance.csv            : component, evr, cumvar
 - <stem>_pca_variance_plot.png                 : EVR line plot
 - <stem>_pca_cumvar_elbow.png                  : cumulative variance with elbow marked
 - <stem>_pca_scores.csv                        : SampleID + PC scores (all components)
 - <stem>_pca_loadings.csv                      : loadings (features x components)
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Tuple

import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import time

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except Exception:
    MATPLOTLIB_AVAILABLE = False

# Hardcoded schema to match other feature-selection scripts
ID_COL = "SampleID"
PERM_COL = "readii_Permutation"
REGION_COL = "readii_Region"
FEATURE_PREFIX = "pred_"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="VIF -> Standardize -> PCA pipeline for FMCIB features")
    p.add_argument("--input", "-i", required=True, help="Path to input FMCIB features CSV")
    p.add_argument("--outdir", "-o", required=True, help="Output directory")
    p.add_argument("--feature-prefix", default=FEATURE_PREFIX, help="Prefix for feature columns [default: %(default)s]")
    p.add_argument("--max-missing-frac", type=float, default=0.50, help="Drop features with missing fraction > this value [default: %(default)s]")
    p.add_argument("--vif-target", type=int, default=500, help="Maximum number of features after VIF selection; stop early if all VIF <= threshold [default: %(default)s]")
    p.add_argument("--vif-threshold", type=float, default=5.0, help="VIF threshold to stop iterative removal when all remaining features have VIF <= this value [default: %(default)s]")
    p.add_argument("--vif-batch-frac", type=float, default=0.05, help="Fraction of features to drop per iteration among those above threshold [default: %(default)s]")
    p.add_argument("--vif-min-batch", type=int, default=10, help="Minimum number of features to drop per iteration [default: %(default)s]")
    p.add_argument("--vif-max-iters", type=int, default=100, help="Maximum number of VIF removal iterations [default: %(default)s]")
    p.add_argument("--vif-log-csv", action="store_true", help="Write per-iteration VIF stats to CSV next to outputs")
    p.add_argument("--pca-max-components", type=int, default=None, help="Max PCA components (defaults to min(n_samples, n_features))")
    p.add_argument("--seed", type=int, default=42, help="Random seed for PCA where applicable")
    return p.parse_args()


def stem_from_path(path: str) -> str:
    name = os.path.basename(path)
    return name[:-4] if name.lower().endswith(".csv") else name


def _elbow_via_max_distance(x_idx: np.ndarray, y_vals: np.ndarray) -> int:
    if y_vals.size <= 2:
        return int(y_vals.size)
    x1, y1 = float(x_idx[0]), float(y_vals[0])
    x2, y2 = float(x_idx[-1]), float(y_vals[-1])
    dx, dy = x2 - x1, y2 - y1
    denom = np.hypot(dx, dy)
    if denom == 0:
        return 1
    distances = np.abs(dy * (x_idx - x1) - dx * (y_vals - y1)) / denom
    elbow_idx0 = int(np.argmax(distances))  # 0-based
    return int(x_idx[elbow_idx0])


def read_and_filter(input_csv: str, feature_prefix: str) -> Tuple[pd.DataFrame, List[str]]:
    df = pd.read_csv(input_csv)
    # Basic schema check
    for c in (ID_COL, PERM_COL, REGION_COL):
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")
    # Row filter original/full (case-insensitive)
    perm = df[PERM_COL].astype(str).str.strip().str.lower()
    region = df[REGION_COL].astype(str).str.strip().str.lower()
    sub = df.loc[(perm == "original") & (region == "full")].copy()
    if sub.empty:
        raise ValueError("No rows after filtering for original/full.")
    # Keep SampleID + pred_*
    feat_cols = [c for c in sub.columns if c.startswith(feature_prefix)]
    if not feat_cols:
        raise ValueError(f"No feature columns starting with '{feature_prefix}' after filtering.")
    keep_cols = ([ID_COL] if ID_COL in sub.columns else []) + feat_cols
    sub = sub.loc[:, keep_cols]
    # Numeric coercion
    sub[feat_cols] = sub[feat_cols].apply(pd.to_numeric, errors="coerce")
    sub[feat_cols] = sub[feat_cols].replace([np.inf, -np.inf], np.nan)
    return sub, feat_cols


def drop_and_impute_missing(df: pd.DataFrame, feat_cols: List[str], max_missing_frac: float) -> Tuple[pd.DataFrame, List[str]]:
    miss_frac = df[feat_cols].isna().mean(axis=0)
    kept = [f for f in feat_cols if float(miss_frac.get(f, 0.0)) <= float(max_missing_frac)]
    if not kept:
        raise ValueError("All features removed by missing-value threshold.")
    X = df[kept].copy()
    # Impute median per feature
    X = X.apply(lambda s: s.fillna(s.median()), axis=0)
    out = pd.concat([df[[ID_COL]].reset_index(drop=True), X.reset_index(drop=True)], axis=1)
    return out, kept


def vif_select_iterative(
    X: np.ndarray,
    feature_names: List[str],
    max_keep: int,
    vif_threshold: float = 5.0,
    batch_frac: float = 0.05,
    min_batch: int = 10,
    max_iters: int = 100,
    log_csv_path: str | None = None,
    random_state: int = 42,
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Iteratively remove highest-VIF features until all VIF <= threshold or feature count <= max_keep.

    Uses Ledoit-Wolf precision on z-scored features; for standardized variables, VIF ≈ diag(precision).
    Returns the reduced matrix, selected feature names, and the final VIF vector for the selected features.
    """
    n, p0 = X.shape
    max_keep = int(max(1, min(max_keep, p0)))
    current_idx = list(range(p0))
    batch_frac = float(max(0.0, min(1.0, batch_frac)))
    min_batch = int(max(1, min(min_batch, p0)))

    iter_logs = []

    def compute_vif(mat: np.ndarray) -> np.ndarray:
        scaler = StandardScaler(with_mean=True, with_std=True)
        Xz = scaler.fit_transform(mat)
        lw = LedoitWolf(store_precision=True)
        lw.fit(Xz)
        precision = lw.precision_
        return np.diag(precision).astype(float)

    # Iterate removal in batches
    for it in range(1, int(max_iters) + 1):
        t0 = time.perf_counter()
        p = len(current_idx)
        vif_vals = compute_vif(X[:, current_idx]) if p > 1 else np.array([1.0])
        max_vif = float(np.nanmax(vif_vals))
        med_vif = float(np.nanmedian(vif_vals))
        # Log
        iter_logs.append({
            "iter": it,
            "n_features": p,
            "max_vif": max_vif,
            "median_vif": med_vif,
            "threshold": float(vif_threshold),
        })
        print(f"[VIF] iter={it} p={p} maxVIF={max_vif:.3f} medVIF={med_vif:.3f} thr={vif_threshold}", flush=True)
        # Stop conditions
        if (p <= max_keep) or (np.all(np.isfinite(vif_vals)) and max_vif <= float(vif_threshold)):
            break
        # Determine batch to drop: among those above threshold, drop top-k by VIF
        above = np.where(vif_vals > float(vif_threshold))[0]
        if above.size == 0:
            break
        # Compute drop count with safeguards
        drop_k = max(int(max(1, round(p * batch_frac))), int(min_batch))
        drop_k = int(min(drop_k, above.size))
        # Do not go below max_keep
        drop_k = int(min(drop_k, p - max_keep)) if p > max_keep else 0
        if drop_k <= 0:
            break
        # Pick local indices with largest VIF among 'above'
        above_sorted = above[np.argsort(vif_vals[above])[::-1]]
        to_drop_local = np.sort(above_sorted[:drop_k])  # ascending for safe deletion
        # Map to original: delete in reverse order to avoid reindexing issues
        removed_names = []
        for idx_local in to_drop_local[::-1]:
            removed_names.append(feature_names[current_idx[idx_local]])
            del current_idx[idx_local]
        t1 = time.perf_counter()
        print(f"[VIF]   dropped={drop_k} -> new_p={len(current_idx)} (dt={(t1-t0):.2f}s)", flush=True)

    sel_idx = np.array(current_idx, dtype=int)
    sel_names = [feature_names[i] for i in sel_idx]
    X_sel = X[:, sel_idx]
    # Optionally write per-iteration log
    if log_csv_path and len(iter_logs) > 0:
        try:
            pd.DataFrame(iter_logs).to_csv(log_csv_path, index=False)
        except Exception:
            pass
    # Final VIFs for selected set
    final_vif = compute_vif(X_sel) if X_sel.shape[1] > 1 else np.array([1.0])
    return X_sel, sel_names, final_vif


def run_pca(
    Xz: np.ndarray,
    sample_ids: List[str],
    out_base: str,
    outdir: str,
    max_comps: int | None,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, int, int]:
    n, p = Xz.shape
    m = min(n, p)
    if max_comps is not None:
        m = max(1, min(m, int(max_comps)))
    pca = PCA(n_components=m, svd_solver="full", random_state=seed)
    scores = pca.fit_transform(Xz)
    evr = np.asarray(pca.explained_variance_ratio_, dtype=float)
    cumvar = np.cumsum(evr)
    idx = np.arange(1, evr.size + 1, dtype=int)
    elbow_k = _elbow_via_max_distance(idx, cumvar)
    # Components to reach 95% cumulative variance
    k95 = int(np.searchsorted(cumvar, 0.95) + 1)
    # Write EVR CSV
    evr_csv = os.path.join(outdir, f"{out_base}_pca_explained_variance.csv")
    pd.DataFrame({
        "component": idx,
        "explained_variance_ratio": evr,
        "cumulative_variance_ratio": cumvar,
    }).to_csv(evr_csv, index=False)
    # Plots
    if MATPLOTLIB_AVAILABLE:
        # EVR plot
        evr_png = os.path.join(outdir, f"{out_base}_pca_variance_plot.png")
        plt.figure(figsize=(6, 4), dpi=150)
        plt.plot(idx, evr, marker='o', lw=1)
        plt.xlabel("Components")
        plt.ylabel("Explained variance ratio")
        plt.title("PCA explained variance")
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(evr_png)
        plt.close()
        # Cumulative + elbow (+ 95% target)
        cum_png = os.path.join(outdir, f"{out_base}_pca_cumvar_elbow.png")
        plt.figure(figsize=(6, 4), dpi=150)
        plt.plot(idx, cumvar, marker='o', lw=1)
        plt.axvline(elbow_k, color='tab:red', ls='--', lw=1, label=f"elbow={elbow_k}")
        plt.axhline(0.95, color='tab:green', ls='--', lw=1, label="target=0.95")
        plt.axvline(k95, color='tab:blue', ls='--', lw=1, label=f"k95={k95}")
        plt.xlabel("Components")
        plt.ylabel("Cumulative explained variance")
        plt.title("PCA cumulative variance (elbow)")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(cum_png)
        plt.close()
    # Console summary
    print(f"PCA: elbow={elbow_k}, k95={k95}, components_fitted={scores.shape[1]}")
    # Scores and loadings
    scores_csv = os.path.join(outdir, f"{out_base}_pca_scores.csv")
    loadings_csv = os.path.join(outdir, f"{out_base}_pca_loadings.csv")
    scores_df = pd.DataFrame(scores, columns=[f"PC{i}" for i in range(1, scores.shape[1] + 1)])
    scores_df.insert(0, ID_COL, sample_ids)
    scores_df.to_csv(scores_csv, index=False)
    # components_: shape (m, p) -> loadings features x components
    loadings = pca.components_.T
    loadings_df = pd.DataFrame(loadings)
    loadings_df.to_csv(loadings_csv, index=False, header=[f"PC{i}" for i in range(1, loadings.shape[1] + 1)])
    # Return for downstream feature ranking
    return evr, loadings, elbow_k, k95


def main() -> int:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Reading: {args.input}")
    try:
        sub, feat_cols = read_and_filter(args.input, args.feature_prefix)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    base = stem_from_path(args.input)
    out_base = base

    # Save raw filtered matrix
    out_raw = os.path.join(args.outdir, f"{out_base}_original_full_features_only.csv")
    sub.to_csv(out_raw, index=False)
    print(f"Filtered matrix: {out_raw} (rows={sub.shape[0]}, features={len(feat_cols)})")

    # Missing-value handling (drop > max_missing_frac, then impute median)
    try:
        sub_imputed, kept_after_missing = drop_and_impute_missing(sub, feat_cols, args.max_missing_frac)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 3
    with open(os.path.join(args.outdir, f"{out_base}_kept_features_after_missing.txt"), "w") as f:
        for nm in kept_after_missing:
            f.write(f"{nm}\n")

    # VIF selection to target
    X_raw = sub_imputed[kept_after_missing].to_numpy(dtype=float)
    sample_ids = sub_imputed[ID_COL].astype(str).to_list()
    # Optional per-iteration log path
    vif_log_path = os.path.join(args.outdir, f"{out_base}_vif_iter_log.csv") if args.vif_log_csv else None
    X_vif, sel_names, vif_final = vif_select_iterative(
        X_raw,
        kept_after_missing,
        max_keep=int(args.vif_target),
        vif_threshold=float(args.vif_threshold),
        batch_frac=float(args.vif_batch_frac),
        min_batch=int(args.vif_min_batch),
        max_iters=int(args.vif_max_iters),
        log_csv_path=vif_log_path,
        random_state=int(args.seed),
    )
    removed = len(kept_after_missing) - len(sel_names)
    print(
        f"VIF selection: resulting shape = (n_samples={X_vif.shape[0]}, n_features={X_vif.shape[1]}). "
        f"Removed {removed} feature(s). Criteria: stop when all VIF <= {args.vif_threshold} or n_features <= {args.vif_target}.",
        flush=True,
    )

    # Write selected feature names and matrices
    out_vif_names = os.path.join(args.outdir, f"{out_base}_vif_selected_features.txt")
    with open(out_vif_names, "w") as f:
        for nm in sel_names:
            f.write(f"{nm}\n")
    out_vif_matrix = os.path.join(args.outdir, f"{out_base}_vif_selected_matrix.csv")
    pd.concat([sub_imputed[[ID_COL]].reset_index(drop=True), pd.DataFrame(X_vif, columns=sel_names)], axis=1).to_csv(out_vif_matrix, index=False)

    # Standardize selected features
    scaler = StandardScaler(with_mean=True, with_std=True)
    Xz = scaler.fit_transform(X_vif)
    out_std = os.path.join(args.outdir, f"{out_base}_standardized_matrix.csv")
    pd.concat([sub_imputed[[ID_COL]].reset_index(drop=True), pd.DataFrame(Xz, columns=sel_names)], axis=1).to_csv(out_std, index=False)

    # PCA and outputs
    evr, loadings, elbow_k, k95 = run_pca(
        Xz, sample_ids, out_base, args.outdir, args.pca_max_components, seed=int(args.seed)
    )

    # Variance-weighted absolute loadings ranking per feature
    # importance_i = sum_j |loading_{i,j}| * EVR_j
    if loadings.size > 0:
        weights = evr.reshape(-1)  # shape (m,)
        importance = np.sum(np.abs(loadings) * weights.reshape(1, -1), axis=1)
        order = np.argsort(-importance)  # descending
        ranked_df = pd.DataFrame({
            "feature": [sel_names[i] for i in order],
            "weighted_score": importance[order],
        })
        ranked_df["rank"] = np.arange(1, len(ranked_df) + 1)
        out_ranked = os.path.join(args.outdir, f"{out_base}_feature_scores_ranked.csv")
        ranked_df.to_csv(out_ranked, index=False)

        # Also write matrices with columns ordered by ranking
        out_std_ranked = os.path.join(args.outdir, f"{out_base}_standardized_matrix_ranked.csv")
        std_ranked_df = pd.concat([
            sub_imputed[[ID_COL]].reset_index(drop=True),
            pd.DataFrame(Xz[:, order], columns=[sel_names[i] for i in order])
        ], axis=1)
        std_ranked_df.to_csv(out_std_ranked, index=False)

        out_raw_ranked = os.path.join(args.outdir, f"{out_base}_vif_selected_matrix_ranked.csv")
        raw_ranked_df = pd.concat([
            sub_imputed[[ID_COL]].reset_index(drop=True),
            pd.DataFrame(X_vif[:, order], columns=[sel_names[i] for i in order])
        ], axis=1)
        raw_ranked_df.to_csv(out_raw_ranked, index=False)

    print("Done.")
    print(f"- Raw filtered matrix:     {out_raw}")
    print(f"- Kept after missing:      {out_base}_kept_features_after_missing.txt")
    print(f"- VIF-selected features:   {out_vif_names}")
    print(f"- VIF-selected matrix:     {out_vif_matrix}")
    print(f"- Standardized matrix:     {out_std}")
    print(f"- PCA EVR CSV/plots:       {out_base}_pca_explained_variance.csv, *_variance_plot.png, *_cumvar_elbow.png")
    print(f"- PCA scores/loadings:     {out_base}_pca_scores.csv / {out_base}_pca_loadings.csv")
    if loadings.size > 0:
        print(f"- Ranked feature scores:   {out_base}_feature_scores_ranked.csv")
        print(f"- Ranked matrices:         {out_base}_vif_selected_matrix_ranked.csv / {out_base}_standardized_matrix_ranked.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
