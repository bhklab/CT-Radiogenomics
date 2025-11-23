#!/usr/bin/env python3
"""
Single-outcome LASSO feature selection on multiple genomic signatures with variance filtering,
cancer_type covariate, test set correlations with p-values and FDR, heatmaps, and bootstrap stability.
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, deque
from sklearn.linear_model import LassoCV, Lasso
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr, pearsonr
import argparse
from typing import Tuple, Dict

# -----------------------------
# Argument parsing
# -----------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Single-outcome LASSO per genomic signature with bootstrapping")
    parser.add_argument("--radiomics", required=True, help="Radiomics CSV with SampleID and cancer_type")
    parser.add_argument("--genomics", required=True, help="Genomics CSV with signatures (rows = base IDs)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    # NOTE: feature prefix and ID suffix are hard-coded per project conventions
    # Feature prefix: 'pred_'
    # ID suffix regex (to strip patch suffix): r'_(\d{4})$'
    parser.add_argument("--duplicate-handling", choices=["keep-all","keep-first","mean","median"], default="keep-all")
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--bootstrap", type=int, default=100, help="Number of bootstrap resamples per signature")
    parser.add_argument("--seed", type=int, default=10)
    parser.add_argument("--standardize", action="store_true")
    parser.add_argument("--var-filter-frac", type=float, default=0.15, help="Remove bottom % variance features")
    return parser.parse_args()

# Hard-coded project conventions
FEATURE_PREFIX = "pred_"
ID_SUFFIX_REGEX = r"_(\d{4})$"

# -----------------------------
# Data loading
# -----------------------------
def load_radiomics(path: str, id_suffix_regex: str, duplicate_handling: str = "keep-all") -> Tuple[pd.DataFrame, pd.DataFrame]:
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
# Correlations and FDR
# -----------------------------
def compute_correlations_with_pvals(X, Y):
    pearson_corr = pd.DataFrame(index=X.columns, columns=Y.columns, dtype=float)
    pearson_pval = pd.DataFrame(index=X.columns, columns=Y.columns, dtype=float)
    spearman_corr = pd.DataFrame(index=X.columns, columns=Y.columns, dtype=float)
    spearman_pval = pd.DataFrame(index=X.columns, columns=Y.columns, dtype=float)
    for f in X.columns:
        for g in Y.columns:
            x = X[f].astype(float).values
            y = Y[g].astype(float).values
            mask = np.isfinite(x) & np.isfinite(y)
            if mask.sum() < 2:
                r_p, p_p = 0.0, 1.0
                r_s, p_s = 0.0, 1.0
            else:
                try:
                    if np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
                        r_p, p_p = 0.0, 1.0
                    else:
                        r_p, p_p = pearsonr(x[mask], y[mask])
                except Exception:
                    r_p, p_p = 0.0, 1.0
                try:
                    r_s, p_s = spearmanr(x[mask], y[mask])
                    if np.isnan(r_s):
                        r_s, p_s = 0.0, 1.0
                except Exception:
                    r_s, p_s = 0.0, 1.0
            pearson_corr.loc[f, g] = float(r_p)
            pearson_pval.loc[f, g] = float(p_p)
            spearman_corr.loc[f, g] = float(r_s)
            spearman_pval.loc[f, g] = float(p_s)
    return pearson_corr, pearson_pval, spearman_corr, spearman_pval

def benjamini_hochberg(pval_df):
    flat = pval_df.values.flatten()
    mask_nan = np.isnan(flat)
    # replace nan with 1.0 for ranking, will restore to nan later
    flat_masked = np.where(mask_nan, 1.0, flat)
    m = len(flat_masked)
    order = np.argsort(flat_masked)
    ranked = flat_masked[order]
    q = np.empty_like(ranked)
    for i, p in enumerate(ranked, start=1):
        q[i-1] = p * m / i
    for i in range(m-2, -1, -1):
        q[i] = min(q[i], q[i+1])
    q_full = np.empty_like(flat_masked)
    q_full[order] = np.clip(q, 0.0, 1.0)
    # restore nan positions
    q_full[mask_nan] = np.nan
    return pd.DataFrame(q_full.reshape(pval_df.shape), index=pval_df.index, columns=pval_df.columns)

def plot_heatmap(corr_df, out_png, title):
    plt.figure(figsize=(max(10,corr_df.shape[1]*0.3), max(8,corr_df.shape[0]*0.3)))
    sns.heatmap(corr_df, cmap="Spectral_r", center=0, vmin=-1, vmax=1)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

# -----------------------------
# Main workflow
# -----------------------------
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load data
    rad_df, rad_map = load_radiomics(args.radiomics, ID_SUFFIX_REGEX, args.duplicate_handling)
    geno = load_genomics(args.genomics)

    # Expand genomics to patch-level and align on full ids (preserve duplicates)
    print(f"[info] radiomics rows={rad_df.shape[0]}, unique base_ids={rad_map['base_id'].nunique()}, genomics rows={geno.shape[0]}")
    geno_exp = expand_genomics_to_radiomics(geno, rad_map)
    rad_aln, gen_aln = align_on_full_ids(rad_df, geno_exp)
    print(f"After matching: radiomics {rad_aln.shape}, genomics {gen_aln.shape}")

    # Train/test split (80/20)
    np.random.seed(args.seed)
    indices = np.arange(rad_aln.shape[0])
    np.random.shuffle(indices)
    split = int(0.8*len(indices))
    train_idx, test_idx = indices[:split], indices[split:]
    rad_train, rad_test = rad_aln.iloc[train_idx], rad_aln.iloc[test_idx]
    gen_train, gen_test = gen_aln.iloc[train_idx], gen_aln.iloc[test_idx]

    # Variance filter on training set (exclude cancer_type)
    feat_cols = [c for c in rad_train.columns if c not in ('SampleID','cancer_type','base_id')]
    var_thresh = np.percentile(rad_train[feat_cols].var(axis=0), args.var_filter_frac*100)
    keep_feats = rad_train[feat_cols].var(axis=0) > var_thresh
    rad_train_filtered = rad_train[keep_feats.index[keep_feats]].copy()
    rad_test_filtered = rad_test[keep_feats.index[keep_feats]].copy()
    rad_train_filtered['cancer_type'] = rad_train['cancer_type']
    rad_test_filtered['cancer_type'] = rad_test['cancer_type']

    # Standardize features
    feat_cols = [c for c in rad_train_filtered.columns if c != 'cancer_type']
    if args.standardize:
        scaler = StandardScaler()
        X_train = scaler.fit_transform(rad_train_filtered[feat_cols])
        X_test = scaler.transform(rad_test_filtered[feat_cols])
    else:
        X_train = rad_train_filtered[feat_cols].values
        X_test = rad_test_filtered[feat_cols].values
    feature_names = feat_cols + ['cancer_type']  # Add cancer_type for LASSO

    # Append cancer_type as numeric covariate
    X_train = np.hstack([X_train, rad_train_filtered['cancer_type'].values.reshape(-1,1)])
    X_test = np.hstack([X_test, rad_test_filtered['cancer_type'].values.reshape(-1,1)])

    # LASSO + bootstrap per signature
    selected_features_all = set()
    for sig in gen_train.columns:
        y = gen_train[sig].values
        lasso_cv = LassoCV(cv=args.cv_folds, random_state=args.seed, max_iter=10000)
        lasso_cv.fit(X_train, y)
        coef_mask = np.abs(lasso_cv.coef_) > 0
        selected_features = [f for f, m in zip(feature_names, coef_mask) if m]
        selected_features_all.update(selected_features)

        # Bootstrap stability
        if len(selected_features) == 0:
            print(f"[{sig}] No features selected by LassoCV; skipping bootstrap and saving empty file")
            df_out = pd.DataFrame({"feature": [], "coef": [], "bootstrap_freq": []})
            df_out.to_csv(os.path.join(args.outdir, f"{sig}_selected_features_bootstrap.csv"), index=False)
        else:
            bootstrap_counts = pd.Series(0, index=feature_names)
            np.random.seed(args.seed)
            for b in range(args.bootstrap):
                sample_idx = np.random.choice(X_train.shape[0], X_train.shape[0], replace=True)
                lasso_b = Lasso(alpha=lasso_cv.alpha_, max_iter=10000)
                lasso_b.fit(X_train[sample_idx], y[sample_idx])
                bootstrap_counts += (np.abs(lasso_b.coef_) > 0).astype(int)
            denom = max(1, args.bootstrap)
            bootstrap_freq_all = bootstrap_counts / denom
            bootstrap_freq = bootstrap_freq_all.loc[selected_features]

            # Save selected features
            df_out = pd.DataFrame({
                "feature": selected_features,
                "coef": lasso_cv.coef_[coef_mask],
                "bootstrap_freq": bootstrap_freq.values
            })
        df_out.to_csv(os.path.join(args.outdir, f"{sig}_selected_features_bootstrap.csv"), index=False)
        print(f"[{sig}] alpha={lasso_cv.alpha_:.4g}, features={len(selected_features)}, bootstrap freq calculated")

    # Deduplicate features for correlation
    selected_features_all = sorted(selected_features_all)
    if len(selected_features_all) == 0:
        print("[info] No features selected across signatures; skipping correlation analysis.")
        print("Analysis complete.")
        return
    X_test_sel = pd.DataFrame(X_test, columns=feature_names)[selected_features_all]

    # Correlations
    pearson_corr, pearson_pval, spearman_corr, spearman_pval = compute_correlations_with_pvals(X_test_sel, gen_test)
    pearson_qval = benjamini_hochberg(pearson_pval)
    spearman_qval = benjamini_hochberg(spearman_pval)

    # Save correlations and heatmaps
    pearson_corr.to_csv(os.path.join(args.outdir,"test_pearson_correlations.csv"))
    pearson_pval.to_csv(os.path.join(args.outdir,"test_pearson_pvals.csv"))
    pearson_qval.to_csv(os.path.join(args.outdir,"test_pearson_qvals_fdr.csv"))
    spearman_corr.to_csv(os.path.join(args.outdir,"test_spearman_correlations.csv"))
    spearman_pval.to_csv(os.path.join(args.outdir,"test_spearman_pvals.csv"))
    spearman_qval.to_csv(os.path.join(args.outdir,"test_spearman_qvals_fdr.csv"))

    plot_heatmap(pearson_corr, os.path.join(args.outdir,"test_pearson_heatmap.png"), "Pearson correlations (test)")
    plot_heatmap(spearman_corr, os.path.join(args.outdir,"test_spearman_heatmap.png"), "Spearman correlations (test)")

    print("Analysis complete.")

if __name__ == "__main__":
    main()
