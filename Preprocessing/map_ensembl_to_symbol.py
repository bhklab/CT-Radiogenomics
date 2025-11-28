#!/usr/bin/env python3
"""Map Ensembl IDs in a genomic matrix to gene symbols using an annotation CSV.

Usage:
    python Preprocessing/map_ensembl_to_symbol.py \
        --annot annotation.csv --matrix input_matrix.tsv --out output_matrix.csv

Behavior:
 - Reads the annotation CSV and uses the columns `Ensembl_ID` and `Gene_Symbol` by default
   (configurable via CLI options) to build a mapping from Ensembl -> Gene Symbol.
 - Strips Ensembl version suffixes (e.g., ENSG00000123456.7 -> ENSG00000123456) when matching.
 - Loads the genomic matrix TSV which must contain an ENSG ID column (default name `ENSGID`).
 - Replaces the values in the ENSGID column with the mapped gene symbol (falls back to the
     original Ensembl ID when no symbol is available) and writes a CSV output.
 - Renames sample columns (all columns except the ENSGID column) by extracting the first
   match of the regex `OCT_\d{6}`. If no match is found the original header is kept.
 - Writes the resulting table back to TSV and prints summary statistics.
"""

from __future__ import annotations

import argparse
import re
import sys
from typing import Dict, Optional

import pandas as pd


# Note: Do not strip Ensembl version suffixes — annotation table keeps full IDs.


def build_mapping(annot_csv: str, ensembl_col: str = 'Ensembl_ID', symbol_col: str = 'Gene_Symbol') -> Dict[str, str]:
    df = pd.read_csv(annot_csv, dtype=str)
    if ensembl_col not in df.columns or symbol_col not in df.columns:
        raise KeyError(f"Annotation CSV must contain columns '{ensembl_col}' and '{symbol_col}'. Found: {df.columns.tolist()}")
    # Use Ensembl IDs exactly as provided in the annotation file (do not strip versions)
    df[ensembl_col] = df[ensembl_col].astype(str)
    mapping = {}
    for ensg, sym in zip(df[ensembl_col].fillna(''), df[symbol_col].fillna('')):
        if ensg == '':
            continue
        mapping[ensg] = sym if sym != '' else ensg
    return mapping


def rename_samples(cols, pattern: str = r'(OCT_\d{6})'):
    regex = re.compile(pattern)
    renamed = []
    renamed_count = 0
    for c in cols:
        m = regex.search(str(c))
        if m:
            new = m.group(1)
            if new != c:
                renamed_count += 1
            renamed.append(new)
        else:
            renamed.append(c)
    return renamed, renamed_count


def process(annot_csv: str, matrix_tsv: str, out_tsv: str, ensgid_col: str = 'ENSGID',
            ensembl_col: str = 'Ensembl_ID', symbol_col: str = 'Gene_Symbol', sample_pattern: str = r'(OCT_\d{6})') -> int:
    mapping = build_mapping(annot_csv, ensembl_col=ensembl_col, symbol_col=symbol_col)

    mat = pd.read_csv(matrix_tsv, sep='\t', dtype=str)
    # find ENSGID column case-insensitive
    ens_col = None
    for col in mat.columns:
        if col.lower() == ensgid_col.lower():
            ens_col = col
            break
    if ens_col is None:
        raise KeyError(f"Matrix TSV must contain an Ensembl ID column named '{ensgid_col}' (case-insensitive). Columns: {mat.columns.tolist()}")

    # Map Ensembl -> Gene Symbol
    original_ids = mat[ens_col].fillna('').astype(str)
    # Match Ensembl IDs exactly as they appear (annotation keeps version suffixes)
    mapped = original_ids.map(lambda x: mapping.get(x, ''))

    # Stats
    n_total = len(mat)
    n_mapped = (mapped != '').sum()
    n_unmapped = n_total - n_mapped

    # Replace ENSGID column values: prefer mapped symbol when available, else use original Ensembl ID
    new_ids = [m if m else orig for orig, m in zip(original_ids, mapped)]
    mat[ens_col] = new_ids

    # --- Duplicate-row handling ---
    # Identify exact duplicate IDs (after mapping). For any duplicate rows (identical ID string),
    # compute each row's average expression across sample columns and keep the row with the
    # higher average. If averages are equal, keep the first occurrence.
    # Note: gene symbols or IDs that differ only by version suffix (e.g., .1/.2) are considered
    # different unless the full ID string is identical.
    other_cols = [c for c in mat.columns if c != ens_col]
    if len(other_cols) > 0:
        # compute numeric row means (do not modify original strings)
        numeric_block = mat[other_cols].apply(pd.to_numeric, errors='coerce')
        row_means = numeric_block.mean(axis=1)

        # group rows by exact ID string
        grouped = mat.groupby(ens_col, sort=False)
        indices_to_drop = []
        n_duplicates_groups = 0
        for gid, group in grouped:
            if len(group) <= 1:
                continue
            n_duplicates_groups += 1
            # select the row (index) with maximum mean; idxmax chooses first on ties
            group_idx = group.index
            # pick index label of maximum mean among these rows
            try:
                idx_keep = row_means.loc[group_idx].idxmax()
            except Exception:
                # fallback: keep first
                idx_keep = group_idx[0]
            # drop all other indices in this group
            for idx in group_idx:
                if idx != idx_keep:
                    indices_to_drop.append(idx)

        if indices_to_drop:
            mat = mat.drop(index=indices_to_drop).reset_index(drop=True)
            # also update mapped stats to reflect removed rows
            removed = len(indices_to_drop)
            print(f"Removed {removed} duplicate rows across {n_duplicates_groups} duplicate ID groups; kept higher-average rows.")
        else:
            print('No duplicate ID rows found to collapse.')
    else:
        print('No sample columns detected; skipping duplicate-row handling.')

    # Rename sample columns (exclude ENSGID column)
    other_cols = [c for c in mat.columns if c != ens_col]
    new_other_cols, renamed_count = rename_samples(other_cols, pattern=sample_pattern)
    # apply renaming
    col_map = {old: new for old, new in zip(other_cols, new_other_cols)}
    mat = mat.rename(columns=col_map)

    # Write output as comma-separated CSV
    mat.to_csv(out_tsv, index=False)

    print(f"Wrote: {out_tsv}")
    print(f"Total rows: {n_total}; mapped gene symbols: {n_mapped}; unmapped: {n_unmapped}")
    print(f"Sample column renames applied: {renamed_count}")
    return 0


def _parse_args(argv=None):
    p = argparse.ArgumentParser(description='Replace Ensembl IDs with Gene Symbols in a genomic matrix TSV and normalize sample headers.')
    p.add_argument('--annot', required=True, help='Annotation CSV path containing Ensembl_ID and Gene_Symbol columns')
    p.add_argument('--matrix', required=True, help='Input genomic matrix TSV (must contain ENSGID column)')
    p.add_argument('--out', required=True, help='Output CSV path')
    p.add_argument('--ensgid-col', default='ENSGID', help='Name of Ensembl ID column in matrix (default: ENSGID)')
    p.add_argument('--annot-ensembl-col', default='Ensembl_ID', help='Column name in annotation CSV for Ensembl IDs')
    p.add_argument('--annot-symbol-col', default='Gene_Symbol', help='Column name in annotation CSV for gene symbols')
    p.add_argument('--sample-pattern', default=r'(OCT_\d{6})', help='Regex pattern to extract sample IDs (default: "(OCT_\\d{6})")')
    return p.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)
    try:
        return process(args.annot, args.matrix, args.out, ensgid_col=args.ensgid_col,
                       ensembl_col=args.annot_ensembl_col, symbol_col=args.annot_symbol_col,
                       sample_pattern=args.sample_pattern)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
