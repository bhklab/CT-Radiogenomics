#!/usr/bin/env python3
"""
Multi-task Lasso radiogenomic feature selection and validation with logging.

Pipeline:
 1) Load radiomics and genomic signature matrices
 2) Harmonize and align sample IDs
 3) Train/test split
 4) MultiTaskLassoCV (5-fold) with 1-SE rule for alpha
 5) Tighten alpha to get <= max-features
 6) Record selected features and coefficient norms
 7) Bootstrap on selected features to estimate stability (top 25 printed)
 8) Compute vectorized Spearman correlations on test set
 9) Produce heatmap and save outputs
"""
from __future__ import annotations
import os, re, argparse
from typing import Tuple, Dict
from collections import defaultdict, deque

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import MultiTaskLassoCV, MultiTaskLasso
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA, SparsePCA
from sklearn.linear_model import LinearRegression
from scipy.stats import spearmanr, pearsonr

# -----------------------------
# Parsing and utility functions
# -----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Multi-task Lasso radiogenomic analysis.")
    p.add_argument("-r", "--radiomics", required=True)
    p.add_argument("-g", "--genomics", required=True)
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--name", default=None)
    p.add_argument("--train-frac", type=float, default=0.8)
    p.add_argument("--max-features", type=int, default=500)
    p.add_argument("--bootstrap", type=int, default=200)
    p.add_argument("--seed", type=int, default=10)
    p.add_argument("--alpha-step", type=float, default=1.25)
    p.add_argument("--standardize", action="store_true")
    p.add_argument("--cluster", choices=["both", "rows", "cols", "none"], default="both")
    p.add_argument("--id-suffix-regex", default=r"_(\d{4})$")
    p.add_argument("--duplicate-handling", choices=["keep-first","keep-all","mean","median"], default="keep-all")
    p.add_argument("--compute-pvals", action="store_true", help="Compute Spearman correlation p-values and FDR-adjusted q-values on test set")
    p.add_argument("--sparse-pca", action="store_true", help="Use SparsePCA for the final dimensionality reduction (fit_transform). k is still chosen by dense PCA to reach 95% cumulative variance")
    p.add_argument("--sparse-alpha", type=float, default=1.0, help="Alpha (sparsity) parameter for SparsePCA; larger => sparser components")
    return p.parse_args()

def ensure_outdir(path: str):
    os.makedirs(path, exist_ok=True)

# -----------------------------
# Data loading / alignment
# -----------------------------
def load_radiomics(path: str, id_suffix_regex: str, duplicate_handling: str = "keep-first") -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(path)
    if "SampleID" not in df.columns:
        raise ValueError("Radiomics file must have 'SampleID' column.")
    df["full_id"] = df["SampleID"].astype(str)
    try:
        pat = re.compile(id_suffix_regex)
    except re.error:
        pat = re.compile(r"_(\d{4})$")
    df["base_id"] = df["full_id"].apply(lambda s: pat.sub("", str(s).strip()))
    feat_cols = [c for c in df.columns if c not in ("SampleID","full_id","base_id")]
    rad_df = df.set_index("full_id")[feat_cols]
    rad_map = df[["full_id","base_id"]].copy()
    rad_df.index = rad_df.index.astype(str).str.strip()
    rad_map["full_id"] = rad_map["full_id"].astype(str).str.strip()
    rad_map["base_id"] = rad_map["base_id"].astype(str).str.strip()
    # Handle duplicates
    if not rad_df.index.is_unique:
        if duplicate_handling=="keep-first":
            rad_df = rad_df[~rad_df.index.duplicated(keep="first")]
            rad_map = rad_map.drop_duplicates(subset=["full_id"], keep="first")
        elif duplicate_handling in ("mean","median"):
            if duplicate_handling=="mean":
                rad_df = rad_df.groupby(rad_df.index).mean(numeric_only=True)
            else:
                rad_df = rad_df.groupby(rad_df.index).median(numeric_only=True)
            rad_map = rad_map.drop_duplicates(subset=["full_id"], keep="first")
    return rad_df, rad_map

def load_genomics(path: str) -> pd.DataFrame:
    g = pd.read_csv(path, index_col=0)
    g.replace([np.inf,-np.inf],np.nan,inplace=True)
    g.dropna(axis=0,how="all", inplace=True)
    g.dropna(axis=1,how="all", inplace=True)
    g = g.apply(pd.to_numeric, errors="coerce")
    g.index = g.index.astype(str).str.strip().str.replace('.', '-', regex=False)
    return g

def expand_genomics_to_radiomics(geno: pd.DataFrame, rad_map: pd.DataFrame) -> pd.DataFrame:
    # Create groups of genomic rows per base_id, preserving their original order
    groups = {}
    for base, df_grp in geno.groupby(level=0):
        # df_grp may be a DataFrame with one or multiple rows; store list of row arrays
        rows_list = [row.astype(float) for _, row in df_grp.iterrows()]
        groups[str(base)] = rows_list

    rows_out = []
    new_index = []
    # Maintain a pointer per base_id to cycle through multiple genomic rows if present
    pointers = {b: 0 for b in groups}

    for full_id, base_id in rad_map[["full_id","base_id"]].itertuples(index=False):
        b = str(base_id)
        if b in groups:
            lst = groups[b]
            ptr = pointers[b]
            rows_out.append(np.asarray(lst[ptr]))
            pointers[b] = (ptr + 1) % len(lst)
        else:
            # This should not happen if rad_map was pre-filtered to remove missing base_ids.
            rows_out.append(np.full(geno.shape[1], np.nan, dtype=float))
        new_index.append(full_id)

    G = pd.DataFrame(np.vstack(rows_out), index=new_index, columns=geno.columns)
    return G

def align_on_full_ids(rad_df: pd.DataFrame, G_expanded: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    pos_map: Dict[str,deque] = defaultdict(deque)
    for pos,fid in enumerate(G_expanded.index):
        pos_map[str(fid)].append(pos)
    gen_rows, rad_keep_pos = [], []
    rad_keep_labels = []
    for i,fid in enumerate(rad_df.index.astype(str)):
        dq = pos_map.get(fid)
        if dq and len(dq)>0:
            pos = dq.popleft()
            gen_rows.append(G_expanded.iloc[pos].values)
            rad_keep_labels.append(fid)
            rad_keep_pos.append(i)
    if not gen_rows:
        raise RuntimeError("No common full SampleIDs after alignment.")
    rad_a = rad_df.iloc[rad_keep_pos].copy()
    gen_a = pd.DataFrame(np.vstack(gen_rows), index=rad_keep_labels, columns=G_expanded.columns)
    if not rad_a.index.equals(gen_a.index):
        raise RuntimeError("Indices do not match after alignment.")
    return rad_a, gen_a

# -----------------------------
# Vectorized Spearman correlation
# -----------------------------
def spearman_corr_and_pvals(X: pd.DataFrame, Y: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    genes = list(Y.columns)
    feats = list(X.columns)
    corr = np.zeros((len(genes), len(feats)))
    pvals = np.ones((len(genes), len(feats)))
    for gi, g in enumerate(genes):
        yv = Y[g].values
        for fi, f in enumerate(feats):
            xv = X[f].values
            try:
                r, p = spearmanr(xv, yv)
            except Exception:
                r, p = np.nan, np.nan
            corr[gi, fi] = r
            pvals[gi, fi] = p
    return (pd.DataFrame(corr, index=genes, columns=feats),
            pd.DataFrame(pvals, index=genes, columns=feats))


def pearson_corr_and_pvals(X: pd.DataFrame, Y: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compute Pearson correlation and p-values between columns of X and columns of Y.
    Returns DataFrames with index = Y.columns (signatures) and columns = X.columns (features/PCs).
    """
    genes = list(Y.columns)
    feats = list(X.columns)
    corr = np.zeros((len(genes), len(feats)))
    pvals = np.ones((len(genes), len(feats)))
    for gi, g in enumerate(genes):
        yv = Y[g].values
        for fi, f in enumerate(feats):
            xv = X[f].values
            # pairwise remove NaNs
            mask = ~np.isnan(xv) & ~np.isnan(yv)
            if mask.sum() < 2:
                r, p = np.nan, np.nan
            else:
                try:
                    r, p = pearsonr(xv[mask], yv[mask])
                except Exception:
                    r, p = np.nan, np.nan
            corr[gi, fi] = r
            pvals[gi, fi] = p
    return (pd.DataFrame(corr, index=genes, columns=feats),
            pd.DataFrame(pvals, index=genes, columns=feats))


def select_elbow_k(cumvar: np.ndarray) -> int:
    """Select number of components using the elbow (max distance to line) on cumulative variance.

    cumvar: 1D array of cumulative explained variance (length n_components)
    Returns k (int number of components) >= 1
    """
    if len(cumvar) == 0:
        return 1
    # Points (x, y)
    x = np.arange(1, len(cumvar) + 1)
    y = cumvar
    # Line from first to last point
    x1, y1 = x[0], y[0]
    x2, y2 = x[-1], y[-1]
    # Compute distance from each point to the line
    # line vector
    lx = x2 - x1
    ly = y2 - y1
    if lx == 0 and ly == 0:
        return 1
    # For each point, compute area of triangle *2 / base length = distance*2; use vectorized formula
    distances = np.abs(lx * (y1 - y) - (x1 - x) * ly) / np.sqrt(lx * lx + ly * ly)
    # pick index of max distance
    idx = int(np.argmax(distances))
    k = max(1, idx + 1)
    return k

def benjamini_hochberg(pval_df: pd.DataFrame) -> pd.DataFrame:
    flat = pval_df.values.flatten()
    m = len(flat)
    order = np.argsort(flat)
    ranked = flat[order]
    q = np.empty_like(ranked)
    for i, p in enumerate(ranked, start=1):
        q[i-1] = p * m / i
    for i in range(m-2, -1, -1):
        q[i] = min(q[i], q[i+1])
    q_full = np.empty_like(flat)
    q_full[order] = np.clip(q, 0.0, 1.0)
    return pd.DataFrame(q_full.reshape(pval_df.shape), index=pval_df.index, columns=pval_df.columns)

# -----------------------------
# Heatmap
# -----------------------------
def plot_heatmap(corr_df: pd.DataFrame, out_png: str, cluster: str="both"):
    cmap = "Spectral_r"
    vmin, vmax = -1, 1
    dpi = 300
    save_pdf = True
    mask = corr_df.isna()
    data_for_cluster = corr_df.fillna(0.0)
    n_feats, n_sigs = corr_df.shape[1], corr_df.shape[0]
    cell_w, cell_h = 0.38, 0.3
    fig_w, fig_h = max(10,n_feats*cell_w), max(8,n_sigs*cell_h)
    if cluster != "none" and (n_feats < 2 or n_sigs < 2):
        cluster = "none"
    if cluster=="none":
        plt.figure(figsize=(fig_w,fig_h))
        ax = sns.heatmap(corr_df, mask=mask, cmap=cmap, center=0, vmin=vmin, vmax=vmax)
        plt.tight_layout()
        plt.savefig(out_png,dpi=dpi)
        if save_pdf: plt.savefig(out_png.replace(".png",".pdf"))
        plt.close()
    else:
        row_cluster = cluster in ("both","rows")
        col_cluster = cluster in ("both","cols")
        g = sns.clustermap(data_for_cluster, mask=mask, cmap=cmap, center=0, vmin=vmin, vmax=vmax,
                           row_cluster=row_cluster, col_cluster=col_cluster, figsize=(fig_w,fig_h))
        g.fig.suptitle("Genomic vs Radiomic correlations", y=1.02)
        plt.savefig(out_png,dpi=dpi,bbox_inches="tight")
        if save_pdf: plt.savefig(out_png.replace(".png",".pdf"),bbox_inches="tight")
        plt.close()

# -----------------------------
# Main workflow
# -----------------------------
def main():
    """Simplified workflow: read harmonized radiomics and genomics matrices,
    perform variance filtering on radiomics, compute PCA on radiomics to cover
    95% cumulative variance, project radiomics to PCs, then compute Pearson
    and Spearman correlations between PCs and genomic signatures.
    """
    args = parse_args()
    ensure_outdir(args.outdir)
    base = args.name or f"{os.path.splitext(os.path.basename(args.radiomics))[0]}__x__{os.path.splitext(os.path.basename(args.genomics))[0]}"

    # Load harmonized matrices directly (expects samples as index)
    rad_path = args.radiomics
    geno_path = args.genomics
    print(f"Loading radiomics from: {rad_path}")
    X_rad = pd.read_csv(rad_path, index_col=0)
    # If SampleID column present instead, use it
    if 'SampleID' in X_rad.columns:
        X_rad.index = X_rad['SampleID'].astype(str).str.strip()
        X_rad = X_rad.drop(columns=['SampleID'])
    print(f"Loaded radiomics matrix: {X_rad.shape[0]} samples x {X_rad.shape[1]} columns")

    print(f"Loading genomics from: {geno_path}")
    geno = load_genomics(geno_path)
    print(f"Loaded genomics matrix: {geno.shape[0]} samples x {geno.shape[1]} signatures")

    # Detect cancer_type-like column before coercing types so we can residualize
    cancer_col = None
    for c in X_rad.columns:
        if str(c).lower() == 'cancer_type' or re.search(r'cancer[_ ]?type', str(c), re.I):
            cancer_col = c
            break

    # If cancer_type present, residualize features by cancer type and remove the column
    if cancer_col is not None:
        print(f"Residualizing radiomics by cancer type column: {cancer_col}")
        # Build design matrix of cancer type categories
        cancer_series = X_rad[cancer_col].astype(str).astype('category')
        dmat = pd.get_dummies(cancer_series, drop_first=True)
        # Prepare feature matrix (drop cancer column) and coerce to numeric
        feats_df = X_rad.drop(columns=[cancer_col]).apply(pd.to_numeric, errors='coerce')
        # Drop fully-NaN columns before regression
        feats_df = feats_df.loc[:, ~feats_df.isna().all(axis=0)]
        if dmat.shape[1] < 1:
            print('Only one cancer_type level present; skipping residualization')
            X_rad = feats_df
        else:
            # Align rows (should match)
            dmat = dmat.loc[feats_df.index]
            lr = LinearRegression()
            try:
                lr.fit(dmat.values, feats_df.values)
                preds = lr.predict(dmat.values)
                residuals = feats_df.values - preds
                X_rad = pd.DataFrame(residuals, index=feats_df.index, columns=feats_df.columns)
                # Save residualized radiomics for traceability
                res_path = os.path.join(args.outdir, f"{base}_radiomics_residualized_by_cancer_type.csv")
                X_rad.to_csv(res_path)
                print(f"Saved residualized radiomics to {res_path}")
            except Exception as e:
                print(f"[warn] Residualization failed: {e}; falling back to numeric coercion without residualization")
                X_rad = feats_df
    else:
        # Coerce numeric for radiomics, drop fully-NaN columns
        X_rad = X_rad.apply(pd.to_numeric, errors='coerce')
        X_rad = X_rad.loc[:, ~X_rad.isna().all(axis=0)]

    # Align samples between matrices (intersection)
    common = X_rad.index.intersection(geno.index)
    if common.empty:
        raise RuntimeError('No overlapping samples between radiomics and genomics matrices')
    if len(common) < len(X_rad.index):
        print(f"Warning: dropping {len(X_rad.index) - len(common)} radiomics samples not found in genomics")
    X_rad = X_rad.loc[common]
    geno_aln = geno.loc[common]
    print(f"After alignment: {X_rad.shape[0]} samples")

    # Variance filtering (drop features <= 15th percentile variance)
    feat_vars = X_rad.var(axis=0, ddof=0)
    thresh = np.percentile(feat_vars.values, 15)
    low_var_feats = feat_vars[feat_vars <= thresh].index.tolist()
    if low_var_feats:
        print(f"Dropping {len(low_var_feats)} low-variance features (<= {thresh:.6g})")
        X_rad = X_rad.drop(columns=low_var_feats)
    print(f"After variance filtering: {X_rad.shape[1]} features remain")

    # Standardize and run PCA on radiomics
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_rad.values)
    pca_full = PCA()
    pca_full.fit(X_scaled)
    cumvar = np.cumsum(pca_full.explained_variance_ratio_)
    # Determine k using elbow method on cumulative variance
    k_elbow = select_elbow_k(cumvar)
    k = int(k_elbow)
    print(f"Selected k={k} components by elbow method")

    # Fit final dimensionality reduction with k95 components.
    # If requested, use SparsePCA (note: explained variance was computed using dense PCA above).
    if args.sparse_pca:
        print(f"Using SparsePCA with alpha={args.sparse_alpha} to compute {k} components")
        spca = SparsePCA(n_components=k, alpha=float(args.sparse_alpha), random_state=int(args.seed))
        # SparsePCA returns a transformed matrix of shape (n_samples, n_components)
        X_pca = spca.fit_transform(X_scaled)
        pca_model = spca
    else:
        pca_model = PCA(n_components=k, random_state=int(args.seed))
        X_pca = pca_model.fit_transform(X_scaled)
    pc_cols = [f"PC{i+1}" for i in range(k)]
    X_pca_df = pd.DataFrame(X_pca, index=X_rad.index, columns=pc_cols)
    pca_out = os.path.join(args.outdir, f"{base}_rad_pca_k{k}.csv")
    X_pca_df.to_csv(pca_out)
    print(f"Saved PCA-transformed radiomics to {pca_out}")

    # Compute correlations between PCs (columns of X_pca_df) and genomic signatures (columns of geno_aln)
    # Ensure geno_aln columns are numeric and aligned
    geno_aln = geno_aln.apply(pd.to_numeric, errors='coerce')
    geno_aln = geno_aln.loc[:, ~geno_aln.isna().all(axis=0)]

    # Pearson
    pcorr, pp = pearson_corr_and_pvals(X_pca_df, geno_aln)
    pq = benjamini_hochberg(pp)
    pcorr.to_csv(os.path.join(args.outdir, f"pearson_pca_k{k}.csv"))
    pp.to_csv(os.path.join(args.outdir, f"pearson_pca_k{k}_pvals.csv"))
    pq.to_csv(os.path.join(args.outdir, f"pearson_pca_k{k}_q.csv"))
    plot_heatmap(pcorr, os.path.join(args.outdir, f"pearson_pca_k{k}.png"), cluster=args.cluster)
    print(f"Saved Pearson correlations (PCs vs genomic signatures)")

    # Spearman
    scorr, sp = spearman_corr_and_pvals(X_pca_df, geno_aln)
    sq = benjamini_hochberg(sp)
    scorr.to_csv(os.path.join(args.outdir, f"spearman_pca_k{k}.csv"))
    sp.to_csv(os.path.join(args.outdir, f"spearman_pca_k{k}_pvals.csv"))
    sq.to_csv(os.path.join(args.outdir, f"spearman_pca_k{k}_q.csv"))
    plot_heatmap(scorr, os.path.join(args.outdir, f"spearman_pca_k{k}.png"), cluster=args.cluster)
    print(f"Saved Spearman correlations (PCs vs genomic signatures)")

    print("Analysis complete.")


if __name__=="__main__":
    main()
