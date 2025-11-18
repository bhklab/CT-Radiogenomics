#!/usr/bin/env python3
"""
Preprocess FMCIB radiomic features (partial pipeline: steps 1–4 from feature_PCA.py).

Included steps (mirrors feature_PCA.py up to variance filtering):
 1) Read CSV
 2) Filter rows: readii_Permutation == 'original' AND readii_Region == 'full'
 3) Keep SampleID + pred_* feature columns
 4) Drop low-variance features by absolute variance cutoff (population variance, ddof=0)

Outputs:
 - <stem>_original_full_filtered.csv     : SampleID + retained features after variance filter
 - <stem>_lowvar_kept_features.txt       : list of kept feature names
 - <stem>_lowvar_removed_features.txt    : list of removed feature names
 - <stem>_volume_lowcorr_features.csv    : features with |Spearman rho(feature, sum)| < threshold (default 0.3)
 - <stem>_pca_scores.csv                 : PCA scores for samples (only rows with valid volume used for correlation step)
 - <stem>_pca_loadings.csv               : PCA loadings for features included in PCA
 - <stem>_pca_explained_variance.csv     : Per-component variance ratios and cumulative
 - <stem>_pca_cumvar_elbow.png/.pdf      : Cumulative variance plot with elbow and 95% target marked
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from typing import List, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    # Optional dependencies; only needed if PCA is enabled
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
except Exception:  # pragma: no cover - allow import to fail until PCA actually needed
    PCA = None
    StandardScaler = None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Preprocess FMCIB radiomic features (up to variance filtering).")
    p.add_argument("--input", "-i", required=True, help="Path to input CSV file")
    p.add_argument("--outdir", "-o", required=True, help="Output directory for filtered outputs")
    p.add_argument("--var-threshold", type=float, default=0.15,
                   help="Absolute variance cutoff: keep features with population variance >= this value [default: %(default)s]")
    p.add_argument("--feature-prefix", default="pred_", help="Prefix identifying feature columns (default: %(default)s)")
    # Volume correlation inputs
    p.add_argument("--volume-csv", required=True, help="CSV containing volume 'sum' values and metadata (requires columns: filepath, Modality, sum)")
    p.add_argument("--volume-filepath-col", default="filepath", help="Column in volume CSV containing file path with sample id as the segment before first '/'")
    p.add_argument("--volume-modality-col", default="Modality", help="Column in volume CSV specifying modality (must be RTSTRUCT or SEG)")
    p.add_argument("--volume-sum-col", default="sum", help="Column in volume CSV with volume sum values")
    p.add_argument("--allowed-modalities", default="RTSTRUCT,SEG", help="Comma-separated list of allowed modalities for volume matching")
    p.add_argument("--abs-corr-threshold", type=float, default=0.3, help="Keep features with absolute Spearman correlation to volume below this threshold [default: %(default)s]")
    # PCA options
    p.add_argument("--pca", action="store_true", help="After volume filtering, run PCA on the low-correlation feature subset")
    p.add_argument("--pca-variance-target", type=float, default=0.95, help="Minimum cumulative explained variance to retain [default: %(default)s]")
    p.add_argument("--pca-max-components", type=int, default=None, help="Maximum number of PCA components to consider (defaults to min(n_samples, n_features))")
    p.add_argument("--pca-plot", action="store_true", help="Save a cumulative variance plot with elbow and 95% target")
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


def filter_rows_and_columns(df: pd.DataFrame, id_col: str, perm_col: str, region_col: str, feature_prefix: str,
                            debug: bool = False) -> Tuple[pd.DataFrame, List[str]]:
    df = df.copy()
    df[perm_col] = normalize_str_series(df[perm_col])
    df[region_col] = normalize_str_series(df[region_col])
    sub = df[(df[perm_col] == "original") & (df[region_col] == "full")].copy()
    if sub.empty:
        raise ValueError("No rows after filtering for readii_Permutation=='original' and readii_Region=='full'.")
    keep_cols = [c for c in sub.columns if c.startswith(feature_prefix)]
    if id_col in sub.columns:
        keep_cols = [id_col] + keep_cols
    if not any(c.startswith(feature_prefix) for c in keep_cols):
        raise ValueError(f"No feature columns starting with '{feature_prefix}' after filtering.")
    sub = sub.loc[:, keep_cols]
    feature_cols = [c for c in keep_cols if c != id_col]
    if debug:
        print(f"[DEBUG] filter_rows_and_columns: rows={sub.shape[0]}, features={len(feature_cols)}")
    return sub, feature_cols


def to_numeric(df: pd.DataFrame, feature_cols: List[str], debug: bool = False) -> pd.DataFrame:
    out = df.copy()
    for c in feature_cols:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    if debug:
        non_finite = int(np.isnan(out[feature_cols].to_numpy(dtype=float)).sum()) + int(np.isinf(out[feature_cols].to_numpy(dtype=float)).sum())
        print(f"[DEBUG] to_numeric: converted features; non-finite count={non_finite}")
    return out


def variance_filter(X: np.ndarray, feature_cols: List[str], threshold: float, debug: bool = False) -> Tuple[np.ndarray, List[str], List[str]]:
    vars_ = np.var(X, axis=0, ddof=0)
    cutoff = float(threshold)
    if debug:
        vmin = float(np.min(vars_)); vmed = float(np.median(vars_)); vq85 = float(np.quantile(vars_, 0.85)); vmax = float(np.max(vars_))
        print(f"[DEBUG] variance_filter: threshold={threshold}, cutoff={cutoff:.6g}, var[min/50%/85%/max]=[{vmin:.4g}/{vmed:.4g}/{vq85:.4g}/{vmax:.4g}]")
    mask = vars_ >= cutoff
    kept = [f for f, m in zip(feature_cols, mask) if bool(m)]
    removed = [f for f, m in zip(feature_cols, mask) if not bool(m)]
    if not any(mask):
        top_idx = int(np.argmax(vars_))
        mask = np.zeros_like(vars_, dtype=bool); mask[top_idx] = True
        kept = [feature_cols[top_idx]]
        removed = [f for i, f in enumerate(feature_cols) if i != top_idx]
        if debug:
            print("[DEBUG] variance_filter: all features below threshold; kept top-variance feature as fallback")
    X_sel = X[:, mask]
    return X_sel, kept, removed


def extract_sample_id_from_filepath(path: str) -> str:
    if not isinstance(path, str):
        return ""
    # SampleID is everything before the first '/'
    i = path.find('/')
    return path if i == -1 else path[:i]


def build_volume_series(volume_df: pd.DataFrame, filepath_col: str, modality_col: str, sum_col: str, allowed_modalities: List[str], debug: bool = False) -> pd.Series:
    # Normalize modality and filter to allowed
    df = volume_df.copy()
    df[modality_col] = df[modality_col].astype(str).str.strip().str.upper()
    allowed_set = {m.strip().upper() for m in allowed_modalities}
    df = df[df[modality_col].isin(allowed_set)].copy()
    # Extract sample id from filepath
    df['__sample_id__'] = df[filepath_col].astype(str).map(extract_sample_id_from_filepath)
    # Coerce sum to numeric
    df[sum_col] = pd.to_numeric(df[sum_col], errors='coerce')
    df = df[np.isfinite(df[sum_col])]
    # If multiple rows per sample, take the maximum sum (conservative)
    agg = df.groupby('__sample_id__', as_index=True)[sum_col].max()
    if debug:
        print(f"[DEBUG] build_volume_series: allowed modalities={sorted(list(allowed_set))}; unique samples with volume={agg.size}")
    agg.name = sum_col
    return agg


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    # Compute Spearman rho via rank-transform + Pearson
    # Requires len>=3 and some variability
    if x.size < 3 or y.size < 3:
        return float('nan')
    if np.all(x == x[0]) or np.all(y == y[0]):
        return float('nan')
    xr = pd.Series(x).rank(method='average').to_numpy(dtype=float)
    yr = pd.Series(y).rank(method='average').to_numpy(dtype=float)
    C = np.corrcoef(xr, yr)
    rho = float(C[0, 1]) if C.shape == (2, 2) else float('nan')
    return rho


def _elbow_via_max_distance(x_idx: np.ndarray, y_vals: np.ndarray) -> int:
    """Return 1-based index of elbow using max distance to line between first and last point.

    x_idx: 1..m indices (shape (m,))
    y_vals: monotonically increasing sequence (shape (m,))
    """
    if y_vals.size <= 2:
        return int(y_vals.size)
    # Line from first to last point
    x1, y1 = float(x_idx[0]), float(y_vals[0])
    x2, y2 = float(x_idx[-1]), float(y_vals[-1])
    dx, dy = x2 - x1, y2 - y1
    denom = np.hypot(dx, dy)
    if denom == 0:
        return 1
    # Distance of each point to the line
    distances = np.abs(dy * (x_idx - x1) - dx * (y_vals - y1)) / denom
    elbow_idx0 = int(np.argmax(distances))  # 0-based
    return int(x_idx[elbow_idx0])


def main() -> int:
    args = parse_args()
    ensure_outdir(args.outdir)

    def vlog(msg: str) -> None:
        if args.verbose:
            print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] [DEBUG] {msg}")

    print(f"Reading: {args.input}")
    t_read0 = time.perf_counter()
    df = pd.read_csv(args.input)
    vlog(f"read_csv done in {time.perf_counter()-t_read0:.2f}s; shape={df.shape}")

    # Column constants
    ID_COL = "SampleID"; PERM_COL = "readii_Permutation"; REGION_COL = "readii_Region"; PREFIX = args.feature_prefix
    for col in [ID_COL, PERM_COL, REGION_COL]:
        if col not in df.columns:
            print(f"ERROR: Missing required column: {col}", file=sys.stderr)
            return 2

    # Row & column filtering
    t_flt0 = time.perf_counter()
    try:
        sub, feature_cols = filter_rows_and_columns(df, ID_COL, PERM_COL, REGION_COL, PREFIX, debug=args.verbose)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 3
    print(f"After filter: {sub.shape[0]} rows, {len(feature_cols)} features")
    vlog(f"row/column filtering took {time.perf_counter()-t_flt0:.2f}s")

    # Numeric conversion
    t_num0 = time.perf_counter()
    sub_num = to_numeric(sub, feature_cols, debug=args.verbose)
    X = sub_num[feature_cols].to_numpy(dtype=float)
    vlog(f"numeric conversion took {time.perf_counter()-t_num0:.2f}s; matrix shape={X.shape}")
    if not np.isfinite(X).all():
        print("ERROR: Non-finite values detected after numeric conversion.", file=sys.stderr)
        return 4

    # Variance filtering
    if X.shape[1] == 0:
        print("ERROR: No features available before variance filtering.", file=sys.stderr)
        return 5
    t_var0 = time.perf_counter()
    X_vf, kept_feats, removed_feats = variance_filter(X, feature_cols, threshold=args.var_threshold, debug=args.verbose)
    print(f"Variance filter: kept {len(kept_feats)}, removed {len(removed_feats)} (>= {args.var_threshold:.4g} variance)")
    vlog(f"variance filtering took {time.perf_counter()-t_var0:.2f}s; new shape={X_vf.shape}")
    if X_vf.shape[1] == 0:
        print("ERROR: All features removed by variance thresholding.", file=sys.stderr)
        return 6

    # Outputs
    base = stem_from_path(args.input)
    out_filtered_csv = os.path.join(args.outdir, f"{base}_original_full_filtered.csv")
    out_lowvar_kept = os.path.join(args.outdir, f"{base}_lowvar_kept_features.txt")
    out_lowvar_removed = os.path.join(args.outdir, f"{base}_lowvar_removed_features.txt")
    out_vol_lowcorr_csv = os.path.join(args.outdir, f"{base}_volume_lowcorr_features.csv")

    t_io0 = time.perf_counter()
    filtered_df = pd.concat([sub_num[[ID_COL]].reset_index(drop=True), pd.DataFrame(X_vf, columns=kept_feats)], axis=1)
    filtered_df.to_csv(out_filtered_csv, index=False)
    with open(out_lowvar_kept, "w") as f:
        for feat in kept_feats:
            f.write(f"{feat}\n")
    with open(out_lowvar_removed, "w") as f:
        for feat in removed_feats:
            f.write(f"{feat}\n")
    vlog(f"wrote outputs in {time.perf_counter()-t_io0:.2f}s")

    # Volume CSV: read and match to samples, then compute Spearman correlations per feature
    t_vol0 = time.perf_counter()
    vol_df = pd.read_csv(args.volume_csv)
    allowed_modalities = [m for m in args.allowed_modalities.split(',') if m.strip()]
    vol_series = build_volume_series(
        vol_df,
        filepath_col=args.volume_filepath_col,
        modality_col=args.volume_modality_col,
        sum_col=args.volume_sum_col,
        allowed_modalities=allowed_modalities,
        debug=args.verbose,
    )
    # Align volume to filtered samples
    sample_ids = sub_num[ID_COL].astype(str).to_list()
    vol_values = np.array([vol_series.get(str(sid), np.nan) for sid in sample_ids], dtype=float)
    # Keep only rows with finite volume
    row_mask = np.isfinite(vol_values)
    n_align = int(row_mask.sum())
    if n_align < 3:
        print("WARN: Fewer than 3 samples with valid volume; skipping correlation step.")
        # Still write empty CSV for consistency
        pd.DataFrame({"feature": [], "spearman_rho": []}).to_csv(out_vol_lowcorr_csv, index=False)
        return 0
    X_aligned = X_vf[row_mask, :]
    y_aligned = vol_values[row_mask]
    # Compute rho per feature
    rhos = []
    for j, feat in enumerate(kept_feats):
        xr = X_aligned[:, j]
        rho = spearman_corr(xr, y_aligned)
        rhos.append((feat, rho))
    rho_df = pd.DataFrame(rhos, columns=["feature", "spearman_rho"]).dropna()
    rho_df["abs_spearman_rho"] = np.abs(rho_df["spearman_rho"]).astype(float)
    lowcorr_df = rho_df[rho_df["abs_spearman_rho"] < float(args.abs_corr_threshold)].reset_index(drop=True)
    lowcorr_df.to_csv(out_vol_lowcorr_csv, index=False)
    vlog(f"volume correlation step took {time.perf_counter()-t_vol0:.2f}s; matched_samples={n_align}, kept_lowcorr_features={lowcorr_df.shape[0]} (|rho| < {args.abs_corr_threshold})")

    # Optional PCA on low-correlation subset
    if args.pca:
        if lowcorr_df.empty:
            print("WARN: No features meet the low-correlation threshold; skipping PCA.")
            print(f"- Volume low-corr features (< {args.abs_corr_threshold:.3g}): {out_vol_lowcorr_csv}")
            return 0
        if PCA is None or StandardScaler is None:
            print("ERROR: scikit-learn is required for PCA but is not available.", file=sys.stderr)
            return 7
        # Build matrix with only low-corr features and aligned rows (finite volume)
        lowcorr_feats = [f for f in lowcorr_df["feature"].tolist() if f in kept_feats]
        if len(lowcorr_feats) < 2:
            print("WARN: Fewer than 2 low-correlation features available; skipping PCA.")
            print(f"- Volume low-corr features (< {args.abs_corr_threshold:.3g}): {out_vol_lowcorr_csv}")
            return 0
        X_lowcorr = sub_num[lowcorr_feats].to_numpy(dtype=float)[row_mask, :]
        # Standardize features
        scaler = StandardScaler(with_mean=True, with_std=True)
        Xz = scaler.fit_transform(X_lowcorr)
        m = Xz.shape[1]
        n = Xz.shape[0]
        max_comps = min(n, m)
        if args.pca_max_components is not None:
            max_comps = max(1, min(max_comps, int(args.pca_max_components)))
        # Fit PCA with full range to inspect EVR
        pca_full = PCA(n_components=max_comps, svd_solver="full", random_state=0)
        X_scores = pca_full.fit_transform(Xz)
        evr = np.asarray(pca_full.explained_variance_ratio_, dtype=float)
        cumvar = np.cumsum(evr)
        idx = np.arange(1, evr.size + 1, dtype=int)
        # Elbow and 95% target
        elbow_k = _elbow_via_max_distance(idx, cumvar)
        k95 = int(np.searchsorted(cumvar, float(args.pca_variance_target)) + 1)
        k_sel = int(max(elbow_k, min(k95, max_comps)))
        # Truncate scores/loadings to selected components
        X_scores_sel = X_scores[:, :k_sel]
        loadings = pca_full.components_[:k_sel, :].T  # shape (features, k)
        # Outputs
        out_scores_csv = os.path.join(args.outdir, f"{base}_pca_scores.csv")
        out_loadings_csv = os.path.join(args.outdir, f"{base}_pca_loadings.csv")
        out_evr_csv = os.path.join(args.outdir, f"{base}_pca_explained_variance.csv")
        # Scores with SampleID for aligned rows
        scores_df = pd.DataFrame(X_scores_sel, columns=[f"PC{i}" for i in range(1, k_sel + 1)])
        scores_df.insert(0, ID_COL, np.array(sample_ids, dtype=str)[row_mask])
        scores_df.to_csv(out_scores_csv, index=False)
        # Loadings with feature names
        loadings_df = pd.DataFrame(loadings, index=lowcorr_feats, columns=[f"PC{i}" for i in range(1, k_sel + 1)])
        loadings_df.index.name = "feature"
        loadings_df.to_csv(out_loadings_csv)
        # Explained variance table
        evr_df = pd.DataFrame({
            "component": idx,
            "explained_variance_ratio": evr,
            "cumulative_variance_ratio": cumvar,
        })
        evr_df.to_csv(out_evr_csv, index=False)
        # Optional plot
        if args.pca_plot:
            out_plot_png = os.path.join(args.outdir, f"{base}_pca_cumvar_elbow.png")
            out_plot_pdf = os.path.join(args.outdir, f"{base}_pca_cumvar_elbow.pdf")
            plt.figure(figsize=(6, 4))
            plt.plot(idx, cumvar, marker='o', lw=1.5)
            plt.axhline(float(args.pca_variance_target), color='tab:green', ls='--', lw=1, label=f"target={args.pca_variance_target:.2f}")
            plt.axvline(k_sel, color='tab:red', ls='--', lw=1, label=f"k={k_sel}")
            plt.xlabel("Components")
            plt.ylabel("Cumulative explained variance")
            plt.title("PCA cumulative variance (elbow + target)")
            plt.grid(alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.savefig(out_plot_png, dpi=200)
            plt.savefig(out_plot_pdf)
            plt.close()
        print(f"PCA: selected {k_sel} components (elbow={elbow_k}, k95={k95}, max={max_comps}); saved scores/loadings.")
        print(f"- PCA scores:        {out_scores_csv}")
        print(f"- PCA loadings:      {out_loadings_csv}")
        print(f"- PCA variance:      {out_evr_csv}")

    print("Done.")
    print(f"- Filtered matrix:   {out_filtered_csv}")
    print(f"- Kept features:     {out_lowvar_kept}")
    print(f"- Removed features:  {out_lowvar_removed}")
    print(f"- Volume low-corr features (< {args.abs_corr_threshold:.3g}): {out_vol_lowcorr_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
