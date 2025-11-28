#!/usr/bin/env python3
"""
Correlate LASSO-selected radiomic features with genomic signatures on an 80/20 test split.

This script assumes you already have CSV(s) listing selected radiomic feature names (one column).
For each provided selected-features CSV, the script will:
 - load radiomics and genomics matrices
 - intersect samples, perform an 80/20 train/test split (random_state=10)
 - fit a StandardScaler on the train radiomics numeric features and save it
 - transform the test radiomics features with the train scaler
 - for the selected features in the CSV, compute Pearson and Spearman correlations
   (and p-values + Benjamini-Hochberg FDR q-values) versus all genomic signatures
 - save CSV correlation matrices and heatmap PNGs to the output directory

Usage example:
./.venv/bin/python3 Data_Analysis/Correlations/LASSO_corr.py \
  -r data/procdata/radiomic_features/.../fmcib_original_features_all_datasets_with_cancer_type.csv \
  -g data/procdata/gene_signatures/combined/hallmarks_all_cohorts_combined.csv \
  -f path/to/selected_features_signatureA.csv path/to/selected_features_signatureB.csv \
  -o data/procdata/radiomic_features/fmcib/features/LASSO_corr_out
"""
from __future__ import annotations
import os
import argparse
import re
from typing import Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr, spearmanr
from joblib import dump


def parse_args():
    p = argparse.ArgumentParser(description="Correlate LASSO-selected radiomic features with genomic signatures")
    p.add_argument("-r", "--radiomics", required=True, help="Radiomics CSV with a 'SampleID' column and feature columns")
    p.add_argument("-g", "--genomics", required=True, help="Genomics CSV (index is sample IDs; columns are signatures)")
    p.add_argument("-f", "--features-csv", nargs='+', required=True, help="One or more CSV files listing selected features (one column)")
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--feature-col", default='feature', help="Name of column in features CSV containing feature names (default: 'feature')")
    p.add_argument("--id-col", default='SampleID', help="Column name in radiomics file containing sample IDs (default: SampleID)")
    p.add_argument("--test-size", type=float, default=0.2)
    p.add_argument("--random-state", type=int, default=10)
    p.add_argument("--cluster", choices=['both','rows','cols','none'], default='both')
    return p.parse_args()


def ensure_outdir(path: str):
    os.makedirs(path, exist_ok=True)


def load_radiomics(path: str, id_col: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if id_col not in df.columns:
        raise ValueError(f"Radiomics file must contain column '{id_col}'")
    df[id_col] = df[id_col].astype(str).str.strip()
    df.set_index(id_col, inplace=True)
    return df


def load_genomics(path: str) -> pd.DataFrame:
    g = pd.read_csv(path, index_col=0)
    g.index = g.index.astype(str).str.strip().str.replace('.', '-', regex=False)
    g = g.apply(pd.to_numeric, errors='coerce')
    return g


def pearson_corr_and_pvals(X: pd.DataFrame, Y: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    genes = list(Y.columns)
    feats = list(X.columns)
    corr = np.zeros((len(genes), len(feats)))
    pvals = np.ones((len(genes), len(feats)))
    for gi, g in enumerate(genes):
        yv = Y[g].values
        for fi, f in enumerate(feats):
            xv = X[f].values
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
    return pd.DataFrame(corr, index=genes, columns=feats), pd.DataFrame(pvals, index=genes, columns=feats)


def spearman_corr_and_pvals(X: pd.DataFrame, Y: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    genes = list(Y.columns)
    feats = list(X.columns)
    corr = np.zeros((len(genes), len(feats)))
    pvals = np.ones((len(genes), len(feats)))
    for gi, g in enumerate(genes):
        yv = Y[g].values
        for fi, f in enumerate(feats):
            xv = X[f].values
            # pairwise mask NaNs to mirror Pearson behavior
            mask = ~np.isnan(xv) & ~np.isnan(yv)
            if mask.sum() < 2:
                r, p = np.nan, np.nan
            else:
                try:
                    r, p = spearmanr(xv[mask], yv[mask])
                except Exception:
                    r, p = np.nan, np.nan
            corr[gi, fi] = r
            pvals[gi, fi] = p
    return pd.DataFrame(corr, index=genes, columns=feats), pd.DataFrame(pvals, index=genes, columns=feats)


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


def plot_heatmap(corr_df: pd.DataFrame, out_png: str, cluster: str='both'):
    vmin, vmax = -1, 1
    mask = corr_df.isna()
    # create colormap and set NaN/masked color
    cmap = sns.color_palette('Spectral_r', as_cmap=True)
    try:
        cmap.set_bad('lightgrey')
    except Exception:
        pass
    cbar_kws = {'ticks': np.linspace(vmin, vmax, 9), 'format': '%.2f'}
    fig_w = max(8, corr_df.shape[1] * 0.3)
    fig_h = max(6, corr_df.shape[0] * 0.25)
    if cluster == 'none' or corr_df.shape[0] < 2 or corr_df.shape[1] < 2:
        plt.figure(figsize=(fig_w, fig_h))
        sns.heatmap(corr_df, mask=mask, cmap=cmap, center=0, vmin=vmin, vmax=vmax, cbar_kws=cbar_kws)
        plt.tight_layout()
        plt.savefig(out_png, dpi=300)
        plt.close()
    else:
        row_cluster = cluster in ('both', 'rows')
        col_cluster = cluster in ('both', 'cols')
        g = sns.clustermap(corr_df, mask=mask, cmap=cmap, center=0, vmin=vmin, vmax=vmax,
                           row_cluster=row_cluster, col_cluster=col_cluster, figsize=(fig_w, fig_h),
                           cbar_kws=cbar_kws)
        # ensure colorbar ticks formatting
        try:
            cax = g.cax
            cax.yaxis.set_ticks(np.linspace(vmin, vmax, 9))
        except Exception:
            pass
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.close()


def read_selected_features(path: str, feature_col: str = None) -> List[str]:
    df = pd.read_csv(path)
    # prefer exact column name; fall back to case-insensitive match or first column
    if feature_col is None:
        feature_col = df.columns[0]
    if feature_col not in df.columns:
        # try case-insensitive match
        cols_lower = {c.lower(): c for c in df.columns}
        if feature_col.lower() in cols_lower:
            feature_col = cols_lower[feature_col.lower()]
        else:
            # fallback to first column and warn
            print(f"Warning: feature column '{feature_col}' not found in {path}; using first column '{df.columns[0]}' instead.")
            feature_col = df.columns[0]
    feats = df[feature_col].dropna().astype(str).str.strip().tolist()
    # remove empty strings
    feats = [f for f in feats if f]
    return feats


def main():
    args = parse_args()
    ensure_outdir(args.outdir)

    # Load data
    rad = load_radiomics(args.radiomics, args.id_col)
    # keep the full list of genomic signatures (columns) from the input file so we can
    # always present all signatures in outputs/plots even if some become all-NaN later
    try:
        full_genome_cols = list(pd.read_csv(args.genomics, index_col=0, nrows=0).columns)
    except Exception:
        full_genome_cols = None
    geno = load_genomics(args.genomics)

    # intersect samples
    common = sorted(set(rad.index.astype(str)).intersection(set(geno.index.astype(str))))
    if not common:
        raise RuntimeError('No intersecting sample IDs between radiomics and genomics')
    rad = rad.loc[common].copy()
    geno = geno.loc[common].copy()

    # Select numeric radiomic features (drop any non-numeric columns)
    rad_num = rad.apply(pd.to_numeric, errors='coerce')
    rad_num = rad_num.loc[:, ~rad_num.isna().all(axis=0)]
    print(f"Using {rad_num.shape[0]} samples and {rad_num.shape[1]} numeric radiomic features (after coercion)")

    # Train/test split on sample IDs
    sample_ids = list(rad_num.index.astype(str))
    train_ids, test_ids = train_test_split(sample_ids, test_size=args.test_size, random_state=args.random_state)
    print(f"Train samples: {len(train_ids)} | Test samples: {len(test_ids)}")

    X_train = rad_num.loc[train_ids]
    X_test = rad_num.loc[test_ids]
    Y_test = geno.loc[test_ids]

    # Standardize using train scaler
    scaler = StandardScaler()
    scaler.fit(X_train.values)
    # save scaler
    scaler_path = os.path.join(args.outdir, 'radiomics_train_scaler.joblib')
    dump(scaler, scaler_path)
    print(f"Saved train scaler to {scaler_path}")

    X_test_scaled = pd.DataFrame(scaler.transform(X_test.values), index=X_test.index, columns=X_test.columns)

    # For each selected-features CSV, compute correlations and save
    for feats_csv in args.features_csv:
        tag = os.path.splitext(os.path.basename(feats_csv))[0]
        feats = read_selected_features(feats_csv, args.feature_col)
        # keep only features present in the radiomics matrix
        found = [f for f in feats if f in X_test_scaled.columns]
        missing = [f for f in feats if f not in X_test_scaled.columns]
        if missing:
            print(f"Warning: {len(missing)} selected features not found in radiomics and will be skipped for '{tag}' (examples: {missing[:5]})")
        if not found:
            print(f"No selected features found in radiomics for {tag}; skipping.")
            continue

        X_sub = X_test_scaled.loc[:, found]

        # Pearson
        p_corr, p_p = pearson_corr_and_pvals(X_sub, Y_test)
        # ensure all genomic signatures are present in the result (preserve original order)
        if full_genome_cols is not None:
            p_corr = p_corr.reindex(index=full_genome_cols, columns=p_corr.columns)
            p_p = p_p.reindex(index=full_genome_cols, columns=p_p.columns)
        print(f"Reindexed Pearson result to {p_corr.shape[0]} genomic signatures x {p_corr.shape[1]} features (filled missing with NaN)")
        p_q = benjamini_hochberg(p_p)
        p_corr.to_csv(os.path.join(args.outdir, f"pearson_{tag}.csv"))
        p_p.to_csv(os.path.join(args.outdir, f"pearson_{tag}_pvals.csv"))
        p_q.to_csv(os.path.join(args.outdir, f"pearson_{tag}_q.csv"))
        plot_heatmap(p_corr, os.path.join(args.outdir, f"pearson_{tag}.png"), cluster=args.cluster)
        print(f"Saved Pearson correlation outputs for {tag} ({p_corr.shape[0]} signatures x {p_corr.shape[1]} features)")

        # Spearman
        s_corr, s_p = spearman_corr_and_pvals(X_sub, Y_test)
        if full_genome_cols is not None:
            s_corr = s_corr.reindex(index=full_genome_cols, columns=s_corr.columns)
            s_p = s_p.reindex(index=full_genome_cols, columns=s_p.columns)
        print(f"Reindexed Spearman result to {s_corr.shape[0]} genomic signatures x {s_corr.shape[1]} features (filled missing with NaN)")
        s_q = benjamini_hochberg(s_p)
        s_corr.to_csv(os.path.join(args.outdir, f"spearman_{tag}.csv"))
        s_p.to_csv(os.path.join(args.outdir, f"spearman_{tag}_pvals.csv"))
        s_q.to_csv(os.path.join(args.outdir, f"spearman_{tag}_q.csv"))
        plot_heatmap(s_corr, os.path.join(args.outdir, f"spearman_{tag}.png"), cluster=args.cluster)
        print(f"Saved Spearman correlation outputs for {tag} ({s_corr.shape[0]} signatures x {s_corr.shape[1]} features)")

    print('LASSO correlation analysis complete.')


if __name__ == '__main__':
    main()
