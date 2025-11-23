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
from sklearn.linear_model import MultiTaskLassoCV, MultiTaskLasso, LinearRegression
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr

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
    rows, new_index, missing_bases = [], [], []
    for full_id, base_id in rad_map[["full_id","base_id"]].itertuples(index=False):
        if base_id in geno.index:
            rows.append(geno.loc[base_id].values.astype(float))
            new_index.append(full_id)
        else:
            missing_bases.append(base_id)
    if not rows:
        raise RuntimeError("No overlapping sample IDs after expansion.")
    G = pd.DataFrame(np.vstack(rows), index=new_index, columns=geno.columns)
    if missing_bases:
        print(f"[warn] {len(set(missing_bases))} base IDs missing in genomics")
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
    rad_df, rad_map = load_radiomics(args.radiomics, args.id_suffix_regex, args.duplicate_handling)
    geno = load_genomics(args.genomics)
    G_exp = expand_genomics_to_radiomics(geno, rad_map)
    rad_aln, gen_aln = align_on_full_ids(rad_df, G_exp)
    print(f"Aligned radiomics: {rad_aln.shape} | Genomics: {gen_aln.shape}")

    # Train/test split
    rng = np.random.RandomState(args.seed)
    indices = np.arange(rad_aln.shape[0])
    rng.shuffle(indices)
    split_point = int(round(args.train_frac*len(indices)))
    train_idx, test_idx = indices[:split_point], indices[split_point:]
    X_train, X_test = rad_aln.iloc[train_idx], rad_aln.iloc[test_idx]
    Y_train, Y_test = gen_aln.iloc[train_idx], gen_aln.iloc[test_idx]

    # ------------------
    # Control for cancer_type
    # ------------------
    # If a `cancer_type` column exists, fit a linear model on the TRAIN set
    # to predict each radiomic feature from one-hot cancer_type, then
    # subtract predicted values to obtain residualized features. Apply the
    # same transformation to the TEST set using the train-fitted coefficients.
    if "cancer_type" in X_train.columns:
        feat_cols = [c for c in X_train.columns if c != "cancer_type"]
        # design matrices
        D_train = pd.get_dummies(X_train["cancer_type"].astype(str), drop_first=False)
        # fit linear model to predict all features at once
        lr = LinearRegression()
        lr.fit(D_train.values, X_train[feat_cols].values)
        pred_train = lr.predict(D_train.values)
        resid_train = X_train[feat_cols].values - pred_train

        # apply to test: ensure same dummy columns exist
        D_test = pd.get_dummies(X_test["cancer_type"].astype(str), drop_first=False)
        D_test = D_test.reindex(columns=D_train.columns, fill_value=0)
        pred_test = lr.predict(D_test.values)
        resid_test = X_test[feat_cols].values - pred_test

        # convert back to dataframes
        X_train_resid = pd.DataFrame(resid_train, index=X_train.index, columns=feat_cols)
        X_test_resid = pd.DataFrame(resid_test, index=X_test.index, columns=feat_cols)
        print(f"Residualized {len(feat_cols)} radiomic features by cancer_type (train-fit applied to test).")
    else:
        X_train_resid = X_train.copy()
        X_test_resid = X_test.copy()

    # ------------------
    # Variance filtering (drop bottom 15% lowest-variance features by TRAIN set)
    # ------------------
    numeric_train = X_train_resid.select_dtypes(include=[np.number])
    if numeric_train.shape[1] == 0:
        raise RuntimeError("No numeric radiomic features available after residualization.")
    variances = numeric_train.var(axis=0, ddof=0)
    cutoff = np.percentile(variances.values, 15)
    keep_feats = variances[variances > cutoff].index.tolist()
    if len(keep_feats) == 0:
        print("[warn] Variance filter removed all features; falling back to keeping all numeric features.")
        keep_feats = list(numeric_train.columns)
    X_train_filtered = X_train_resid[keep_feats].copy()
    X_test_filtered = X_test_resid.reindex(columns=keep_feats).copy()
    print(f"Variance filtering: kept {len(keep_feats)} / {numeric_train.shape[1]} numeric features (>{cutoff:.6g} var cutoff).")

    # Ensure cancer_type is removed before scaling / modeling
    if "cancer_type" in X_train_filtered.columns:
        X_train_filtered = X_train_filtered.drop(columns=["cancer_type"], errors="ignore")
    if "cancer_type" in X_test_filtered.columns:
        X_test_filtered = X_test_filtered.drop(columns=["cancer_type"], errors="ignore")

    # ------------------
    # Z-score standardize using TRAIN fit (required for Lasso)
    # ------------------
    scaler = StandardScaler()
    X_train_vals = scaler.fit_transform(X_train_filtered.values)
    X_test_vals = scaler.transform(X_test_filtered.values)
    # store updated DataFrame references for later mapping
    X_train = pd.DataFrame(X_train_vals, index=X_train_filtered.index, columns=X_train_filtered.columns)
    X_test = pd.DataFrame(X_test_vals, index=X_test_filtered.index, columns=X_test_filtered.columns)

    # MultiTaskLassoCV
    lcv = MultiTaskLassoCV(cv=5, random_state=args.seed)
    lcv.fit(X_train_vals, Y_train.values)
    alphas, mse_path = lcv.alphas_, lcv.mse_path_
    mean_mse, se_mse = mse_path.mean(axis=0), mse_path.std(axis=0)/np.sqrt(mse_path.shape[0])
    threshold = mean_mse.min() + se_mse[np.argmin(mean_mse)]
    alpha_candidates = [a for a,m in zip(alphas,mean_mse) if m<=threshold]
    alpha_start = max(alpha_candidates) if alpha_candidates else float(lcv.alpha_)
    print(f"Alpha CV best: {lcv.alpha_:.6g}, 1-SE start: {alpha_start:.6g}")

    # Feature tightening
    def fit_count(alpha):
        model = MultiTaskLasso(alpha=alpha, random_state=args.seed)
        model.fit(X_train_vals, Y_train.values)
        coefs = model.coef_
        nz_mask = (coefs!=0).any(axis=0)
        return model, int(nz_mask.sum()), nz_mask

    alpha_curr = alpha_start
    step_factor = max(1.05, args.alpha_step)
    model_curr, feat_count, nz_mask = fit_count(alpha_curr)
    print(f"Alpha {alpha_curr:.6g} -> selected features: {feat_count}")
    while feat_count>args.max_features:
        alpha_curr *= step_factor
        model_curr, feat_count, nz_mask = fit_count(alpha_curr)
        print(f"Alpha {alpha_curr:.6g} -> selected features: {feat_count}")
        if alpha_curr>alphas.max()*1000: break
    print(f"Final alpha: {alpha_curr:.6g}, features: {feat_count}")

    selected_features = list(X_train.columns[nz_mask])
    coefs_full = model_curr.coef_.T[nz_mask]
    coef_norm = np.linalg.norm(coefs_full,axis=1)
    coef_df = pd.DataFrame(coefs_full, index=selected_features, columns=Y_train.columns)
    sel_df = pd.DataFrame({"feature":selected_features,"coef_l2_norm":coef_norm}).sort_values("coef_l2_norm",ascending=False)

    # Bootstrap stability
    B = max(1,args.bootstrap)
    stability_counts = np.zeros(len(selected_features),dtype=int)
    X_train_sel = X_train_vals[:,nz_mask]
    Y_train_vals_arr = Y_train.values
    for b in range(B):
        idx = rng.choice(X_train_sel.shape[0], size=X_train_sel.shape[0], replace=True)
        mb = MultiTaskLasso(alpha=alpha_curr, random_state=args.seed)
        mb.fit(X_train_sel[idx], Y_train_vals_arr[idx])
        nz_b = (mb.coef_!=0).any(axis=0)
        stability_counts += nz_b.astype(int)
        if (b+1) % 50 == 0 or b==B-1:
            freq = stability_counts/(b+1)
            top25 = sel_df.copy()
            top25["freq"] = freq
            top25_sorted = top25.sort_values("freq",ascending=False).head(25)
            print(f"[Bootstrap {b+1}/{B}] Top 25 stable features:\n{top25_sorted[['feature','freq']].to_string(index=False)}\n")

    stability_df = pd.DataFrame({"feature":selected_features,"bootstrap_freq":stability_counts/B}).sort_values("bootstrap_freq",ascending=False)

    # Test-set Spearman correlations (vectorized)
    X_test_sel = X_test[selected_features]
    test_corr_df, test_pval_df = spearman_corr_and_pvals(X_test_sel, Y_test)
    if args.compute_pvals:
        test_qval_df = benjamini_hochberg(test_pval_df)

    # Save outputs
    sel_df.to_csv(os.path.join(args.outdir,f"{base}_selected_features.csv"),index=False)
    coef_df.to_csv(os.path.join(args.outdir,f"{base}_coefficients.csv"),index=True)
    stability_df.to_csv(os.path.join(args.outdir,f"{base}_bootstrap_stability.csv"),index=False)
    test_corr_df.to_csv(os.path.join(args.outdir,f"{base}_test_spearman_corr.csv"),index=True)
    if args.compute_pvals:
        test_pval_df.to_csv(os.path.join(args.outdir,f"{base}_test_spearman_pvals.csv"),index=True)
        test_qval_df.to_csv(os.path.join(args.outdir,f"{base}_test_spearman_qvals_fdr.csv"),index=True)
    plot_heatmap(test_corr_df, os.path.join(args.outdir,f"{base}_test_spearman_corr_heatmap.png"), cluster=args.cluster)

    print("Analysis complete.")
    print(f"Selected features: {len(selected_features)}, saved to {args.outdir}")

if __name__=="__main__":
    main()
