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
from sklearn.decomposition import PCA
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
    args = parse_args()
    ensure_outdir(args.outdir)
    base = args.name or f"{os.path.splitext(os.path.basename(args.radiomics))[0]}__x__{os.path.splitext(os.path.basename(args.genomics))[0]}"

    # Load / align
    # Keep all radiomics rows (do not drop duplicates) to preserve original patch-level rows
    rad_df, rad_map = load_radiomics(args.radiomics, args.id_suffix_regex, duplicate_handling='keep-all')
    geno = load_genomics(args.genomics)

    # Identify base IDs present in genomics
    geno_bases = set(geno.index.astype(str))
    rad_bases = set(rad_map['base_id'].astype(str))
    missing_bases = sorted(rad_bases - geno_bases)
    if missing_bases:
        # Drop all radiomics rows whose base_id is not present in genomics
        drop_mask = rad_map['base_id'].astype(str).isin(missing_bases)
        drop_full_ids = rad_map.loc[drop_mask, 'full_id'].astype(str).tolist()
        rad_map = rad_map.loc[~drop_mask].copy()
        rad_df = rad_df.loc[rad_map['full_id']].copy()
        # Save list of dropped radiomics full_ids for traceability
        drop_path = os.path.join(args.outdir, f"{base}_dropped_radiomics_missing_genomics_full_ids.csv")
        pd.DataFrame({'full_id': drop_full_ids, 'base_id': [s.split('_')[0] if '_' in s else s for s in drop_full_ids]}).to_csv(drop_path, index=False)
        print(f"[info] Dropped {len(drop_full_ids)} radiomics rows belonging to {len(missing_bases)} missing base_ids; list saved to {drop_path}")

    # Now expand genomics to radiomics full_ids using the filtered rad_map
    G_exp = expand_genomics_to_radiomics(geno, rad_map)
    rad_aln, gen_aln = align_on_full_ids(rad_df, G_exp)
    print(f"Aligned radiomics: {rad_aln.shape} | Genomics: {gen_aln.shape}")

    # Extract numeric radiomics matrix (coerce non-numeric -> NaN, drop fully-NaN cols)
    X_rad = rad_aln.copy().apply(pd.to_numeric, errors='coerce')
    X_rad = X_rad.loc[:, ~X_rad.isna().all(axis=0)]
    print(f"Numeric radiomics features: {X_rad.shape[1]} columns")

    # Note: do not drop `cancer_type` here — keep it available for residualization below

    # Residualize radiomics features by cancer type (if present)
    cancer_col = None
    for c in rad_aln.columns:
        if str(c).lower() == 'cancer_type' or re.search(r'cancer[_ ]?type', str(c), re.I):
            cancer_col = c
            break
    if cancer_col and cancer_col in rad_aln.columns:
        print(f"Residualizing radiomics by cancer type column: {cancer_col}")
        cancer_series = rad_aln[cancer_col].astype(str).astype('category')
        dmat = pd.get_dummies(cancer_series, drop_first=True)
        if dmat.shape[1] < 1:
            print('Only one cancer_type level present; skipping residualization')
        else:
            lr = LinearRegression()
            # align X_rad rows with dmat
            lr.fit(dmat.values, X_rad.values)
            preds = lr.predict(dmat.values)
            residuals = X_rad.values - preds
            X_rad = pd.DataFrame(residuals, index=X_rad.index, columns=X_rad.columns)
            res_path = os.path.join(args.outdir, f"{base}_radiomics_residualized_by_cancer_type.csv")
            X_rad.to_csv(res_path)
            print(f"Saved residualized radiomics to {res_path}")
    else:
        print('No cancer_type column found; skipping residualization')

    # After residualization, remove the cancer_type column from X_rad if it exists
    if cancer_col and cancer_col in X_rad.columns:
        try:
            X_rad = X_rad.drop(columns=[cancer_col])
            print(f"Dropped '{cancer_col}' from feature matrix before variance filtering/PCA.")
        except Exception:
            # ignore if not present
            pass

    # Variance filtering: drop features with variance below 15th percentile
    feat_vars = X_rad.var(axis=0, ddof=0)
    thresh = np.percentile(feat_vars.values, 15)
    low_var_feats = feat_vars[feat_vars <= thresh].index.tolist()
    if low_var_feats:
        print(f"Dropping {len(low_var_feats)} low-variance features (<= {thresh:.6g})")
        X_rad = X_rad.drop(columns=low_var_feats)
    print(f"After variance filtering: {X_rad.shape[1]} features remain")

    # Save the aligned matrices after variance filtering
    aligned_rad_path = os.path.join(args.outdir, f"{base}_radiomics_aligned_filtered.csv")
    aligned_geno_path = os.path.join(args.outdir, f"{base}_genomics_aligned.csv")
    X_rad.to_csv(aligned_rad_path)
    # Save the in-memory aligned genomics (one-to-one with rad_aln) directly
    gen_aln.to_csv(aligned_geno_path)
    print(f"Saved aligned radiomics (filtered) to {aligned_rad_path}")
    print(f"Saved aligned genomics to {aligned_geno_path}")

    # Standardize (z-score) radiomics
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_rad.values)
    # Fit PCA and pick components to cover 90% cumulative variance
    pca_full = PCA()
    pca_full.fit(X_scaled)
    cumvar = np.cumsum(pca_full.explained_variance_ratio_)
    # k for 90% cumvar
    k90 = int(np.searchsorted(cumvar, 0.90) + 1)
    print(f"Selected k={k90} components to reach 0.90 cumulative variance")

    # k via elbow method
    k_elbow = select_elbow_k(cumvar)
    print(f"Selected k={k_elbow} components by elbow method")

    # k for 95% cumvar (user requested marker)
    k95 = int(np.searchsorted(cumvar, 0.95) + 1)
    print(f"Selected k={k95} components to reach 0.95 cumulative variance")

    # Save cumulative variance + elbow plot with dotted 95% line
    try:
        x = np.arange(1, len(cumvar) + 1)
        plt.figure(figsize=(8, 5))
        plt.plot(x, cumvar, marker='o', linestyle='-')
        # mark elbow point
        elbow_x = k_elbow
        elbow_y = cumvar[elbow_x - 1]
        plt.scatter([elbow_x], [elbow_y], color='red', zorder=5)
        plt.vlines(elbow_x, ymin=0, ymax=elbow_y, colors='red', linestyles='--')
        plt.hlines(elbow_y, xmin=1, xmax=elbow_x, colors='red', linestyles='--')
        # dotted horizontal at 95% and vertical at k95
        plt.hlines(0.95, xmin=1, xmax=len(cumvar), colors='gray', linestyles=':')
        plt.vlines(k95, ymin=0, ymax=0.95, colors='gray', linestyles=':')
        plt.xlabel('Number of components')
        plt.ylabel('Cumulative explained variance')
        plt.title('PCA cumulative explained variance with elbow')
        plt.grid(alpha=0.3)
        savepath = os.path.join(args.outdir, 'pca_cumulative_variance_elbow.png')
        plt.tight_layout()
        plt.savefig(savepath, dpi=300)
        plt.close()
        print(f"Saved cumulative variance elbow plot to {savepath}")
    except Exception as e:
        print(f"[warning] Could not save elbow plot: {e}")

    # Helper to run PCA, save transformed, and compute correlations
    def run_pca_and_corr(k, tag):
        pca = PCA(n_components=k)
        X_pca = pca.fit_transform(X_scaled)
        pc_cols = [f"PC{i+1}" for i in range(k)]
        X_pca_df = pd.DataFrame(X_pca, index=X_rad.index, columns=pc_cols)
        pca_out = os.path.join(args.outdir, f"rad_pca_{tag}.csv")
        X_pca_df.to_csv(pca_out)
        print(f"Saved PCA-transformed radiomics ({tag}) to {pca_out}")

        # Use in-memory aligned genomics (gen_aln) — already aligned row-for-row with rad_aln
        if gen_aln.shape[0] != X_pca_df.shape[0]:
            raise RuntimeError(f"Aligned genomics rows ({gen_aln.shape[0]}) do not match PCA radiomics rows ({X_pca_df.shape[0]}). Aborting.")
        geno_for_corr = gen_aln.copy()

        # Pearson
        pcorr, pp = pearson_corr_and_pvals(X_pca_df, geno_for_corr)
        pq = benjamini_hochberg(pp)
        pcorr.to_csv(os.path.join(args.outdir, f"pearson_{tag}.csv"))
        pp.to_csv(os.path.join(args.outdir, f"pearson_{tag}_pvals.csv"))
        pq.to_csv(os.path.join(args.outdir, f"pearson_{tag}_q.csv"))
        plot_heatmap(pcorr, os.path.join(args.outdir, f"pearson_{tag}.png"), cluster=args.cluster)
        print(f"Saved Pearson ({tag})")

        # Spearman
        scorr, sp = spearman_corr_and_pvals(X_pca_df, geno_for_corr)
        sq = benjamini_hochberg(sp)
        scorr.to_csv(os.path.join(args.outdir, f"spearman_{tag}.csv"))
        sp.to_csv(os.path.join(args.outdir, f"spearman_{tag}_pvals.csv"))
        sq.to_csv(os.path.join(args.outdir, f"spearman_{tag}_q.csv"))
        plot_heatmap(scorr, os.path.join(args.outdir, f"spearman_{tag}.png"), cluster=args.cluster)
        print(f"Saved Spearman ({tag})")

    # Run for 90%-based components and elbow-based components
    run_pca_and_corr(k90, '90')
    run_pca_and_corr(k_elbow, 'elbow')

    print("Analysis complete.")

if __name__=="__main__":
    main()
