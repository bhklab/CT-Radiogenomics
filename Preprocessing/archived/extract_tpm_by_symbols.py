#!/usr/bin/env python3
"""
Extract TPM rows by gene symbols using an annotation map (Symbol/Synonyms -> NCBI Gene ID),
and align/rename TPM sample IDs to match the sample IDs and order in a reference matrix. (Needed for NSCLC-Radiogenomics data)

Inputs (TSV):
  1) raw_genes_tsv: RNA-seq matrix (genes x samples) or a list containing the genes of interest (gene symbols).
      - First column is assumed to be gene symbols unless --genes-col is provided.
      - Sample column names (excluding the first column) define the desired sample IDs and their order in the output.
  2) tpm_tsv: TPM-normalized RNA-seq matrix indexed by NCBI Gene ID (rows) with sample columns.
      - First column is assumed to be gene ID unless --tpm-id-col is provided.
  3) annot_tsv: annotation table that maps GeneID <-> Symbol and contains optional Synonyms (pipe-separated).
    4) sample_map_tsv: mapping of TPM sample IDs (as they appear in #2) to desired sample IDs (as they appear in #1).
         - For this project, TPM (TSV2) columns are GEO IDs and desired IDs are NSCLC; the map should have columns GEO -> NSCLC.

Pre-filtering step:
    - Before any matching, the raw RNA-seq matrix (TSV1) is checked for genes that have no usable values across all samples
        (all entries are NA, 0, or ±Inf after numeric coercion). Such genes are removed from consideration.

Behavior:
  - For each gene symbol present in raw_genes_tsv (in order), find a matching row in annot_tsv where
    either the Symbol column equals the gene symbol or the gene symbol appears in Synonyms (split by '|').
  - Resolve to a single GeneID preferring a primary Symbol match; if multiple GeneIDs remain, use the first
    and record a warning.
  - Fetch the TPM row for that GeneID and add it to the output matrix, with the row labeled by the gene symbol.
  - Symbols that cannot be resolved (no GeneID match, or GeneID missing in TPM) are reported in a side file.

Output:
    - CSV matrix: first column 'GeneSymbol', remaining columns are the desired sample IDs from #1, in the same order.
  - Optional missing report (CSV) with reason codes.

Usage example:
    python Preprocessing/extract_tpm_by_symbols.py \
    --raw-genes data/rawdata/GSE103584_R01_NSCLC_RNAseq.tsv \
    --tpm data/rawdata/GSE103584_norm_counts_TPM_GRCh38.p13_NCBI.tsv \
    --annot data/rawdata/Human.GRCh38.p13.annot.tsv \
        --sample-map data/rawdata/NSCLC_sampleid_map.tsv \
        --map-from-col GEO --map-to-col NSCLC \
        --output data/procdata/NSCLC_TPM_by_symbols.csv \
        --missing-out data/procdata/NSCLC_TPM_missing.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np


def read_raw_genes_and_samples(path: Path, genes_col: Optional[str]) -> Tuple[List[str], List[str]]:
    """Return (gene_symbols, sample_ids_in_order) from a raw genes matrix/list TSV.
    The first column is gene symbols (or --genes-col). Remaining columns are sample IDs.
    """
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if genes_col is None:
        genes_col = df.columns[0]
    if genes_col not in df.columns:
        raise ValueError(f"--genes-col '{genes_col}' not found in {path}")

    # Sample IDs: all columns except the genes column
    sample_ids = [c for c in df.columns if c != genes_col]

    # If we have sample columns, perform pre-filtering of genes with no usable values
    if len(sample_ids) > 0:
        numeric = df[sample_ids].apply(pd.to_numeric, errors="coerce")
        arr = numeric.to_numpy()
        valid = np.isfinite(arr) & (arr != 0)
        good_rows = valid.any(axis=1)
        removed = int((~good_rows).sum())
        if removed > 0:
            print(
                f"Filtered out {removed} genes with all NA/0/Inf across all samples from {path}",
                file=sys.stderr,
            )
        df = df.loc[good_rows].copy()

    genes = df[genes_col].astype(str).str.strip().tolist()
    genes = [g for g in genes if g != ""]
    return genes, sample_ids


def read_tpm(path: Path, id_col: Optional[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if id_col is None:
        id_col = df.columns[0]
    if id_col not in df.columns:
        raise ValueError(f"TPM id column '{id_col}' not found in {path}")
    df = df.set_index(id_col)
    # Keep as string; numeric conversion is optional if needed later
    return df


def build_annotation_map(
    annot_df: pd.DataFrame,
    id_col: str,
    symbol_col: str,
    synonyms_col: str,
    case_insensitive: bool = True,
) -> Dict[str, List[Tuple[str, str]]]:
    """Return map from normalized symbol -> list of (gene_id, match_type) where match_type is 'symbol' or 'synonym'."""
    for c in (id_col, symbol_col):
        if c not in annot_df.columns:
            raise ValueError(f"Annotation column '{c}' not found")
    if synonyms_col not in annot_df.columns:
        # Create an empty synonyms column if missing
        annot_df[synonyms_col] = ""

    def norm(s: str) -> str:
        s = (s or "").strip()
        return s.lower() if case_insensitive else s

    mapping: Dict[str, List[Tuple[str, str]]] = {}
    for _, row in annot_df.iterrows():
        gid = str(row[id_col]).strip()
        sym = str(row[symbol_col]).strip()
        syn = str(row[synonyms_col]).strip()
        if gid == "" and sym == "":
            continue
        if sym != "":
            mapping.setdefault(norm(sym), []).append((gid, "symbol"))
        if syn != "":
            # Split on '|', trim and add
            for s in syn.split("|"):
                s2 = s.strip()
                if s2:
                    mapping.setdefault(norm(s2), []).append((gid, "synonym"))
    return mapping


def resolve_gene_id(
    symbol: str,
    sym2gid: Dict[str, List[Tuple[str, str]]],
) -> Tuple[Optional[str], Optional[str]]:
    """Return (gene_id, reason). If not found, reason explains why; if found with ambiguity, choose best and reason may note 'ambiguous'."""
    key = symbol
    key_l = symbol.lower()
    # Prefer exact key; else lower-case key
    candidates = sym2gid.get(key) or sym2gid.get(key_l)
    if not candidates:
        return None, "not-in-annotations"
    # Prefer primary symbol matches over synonyms
    primary = [gid for gid, mt in candidates if mt == "symbol"]
    if len(primary) == 1:
        return primary[0], None
    if len(primary) > 1:
        return primary[0], f"ambiguous-symbol({len(primary)})"
    # Use first synonym match
    syns = [gid for gid, mt in candidates if mt == "synonym"]
    if len(syns) == 0:
        return None, "no-valid-candidates"
    if len(syns) == 1:
        return syns[0], None
    return syns[0], f"ambiguous-synonym({len(syns)})"


def default_output_path(raw_genes_path: Path) -> Path:
    stem = raw_genes_path.stem
    return raw_genes_path.with_name(f"{stem}_TPM_by_symbols.csv")


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract TPM rows by Gene Symbols using annotation mapping.")
    p.add_argument("--raw-genes", required=True, type=Path, help="TSV containing gene symbols of interest (first column by default)")
    p.add_argument("--tpm", required=True, type=Path, help="TSV of TPM counts indexed by NCBI Gene ID (first column by default)")
    p.add_argument("--annot", required=True, type=Path, help="TSV annotations with columns for GeneID, Symbol, Synonyms")
    p.add_argument("--output", type=Path, default=None, help="Output CSV path (default: <raw_genes>_TPM_by_symbols.csv)")
    p.add_argument("--missing-out", type=Path, default=None, help="Optional CSV report of missing/ambiguous symbols")
    p.add_argument("--genes-col", type=str, default=None, help="Column name in raw genes TSV containing symbols (default: first column)")
    p.add_argument("--tpm-id-col", type=str, default=None, help="Column name in TPM TSV for gene ID (default: first column)")
    p.add_argument("--annot-id-col", type=str, default="GeneID", help="Annotation column name for gene ID [default: GeneID]")
    p.add_argument("--annot-symbol-col", type=str, default="Symbol", help="Annotation column name for gene symbol [default: Symbol]")
    p.add_argument("--annot-synonyms-col", type=str, default="Synonyms", help="Annotation column name for synonyms [default: Synonyms]")
    p.add_argument("--case-insensitive", action="store_true", help="Case-insensitive matching of symbols and synonyms")
    p.add_argument("--sample-map", type=Path, required=True, help="TSV mapping from TPM sample IDs to desired sample IDs")
    p.add_argument(
        "--map-from-col",
        type=str,
        default="GEO",
        help="Mapping column name for TPM (TSV2) sample ID [default: GEO]",
    )
    p.add_argument(
        "--map-to-col",
        type=str,
        default="NSCLC",
        help="Mapping column name for desired sample ID (from raw TSV1) [default: NSCLC]",
    )
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    genes, desired_samples = read_raw_genes_and_samples(args.raw_genes, args.genes_col)
    if len(genes) == 0:
        print("No gene symbols found in raw genes file", file=sys.stderr)
        return 2

    tpm_df = read_tpm(args.tpm, args.tpm_id_col)
    annot_df = pd.read_csv(args.annot, sep="\t", dtype=str, keep_default_na=False)
    sym2gid = build_annotation_map(
        annot_df,
        id_col=args.annot_id_col,
        symbol_col=args.annot_symbol_col,
        synonyms_col=args.annot_synonyms_col,
        case_insensitive=args.case_insensitive,
    )

    # Build and apply sample ID mapping GEO -> NSCLC, then align columns to desired_samples order
    map_df = pd.read_csv(args.sample_map, sep="\t", dtype=str, keep_default_na=False)
    if args.map_from_col not in map_df.columns or args.map_to_col not in map_df.columns:
        raise ValueError("Sample map must contain columns '{}' and '{}'".format(args.map_from_col, args.map_to_col))
    # Normalize and deduplicate mapping
    m_from = map_df[args.map_from_col].astype(str).str.strip()
    m_to = map_df[args.map_to_col].astype(str).str.strip()
    dedup = pd.DataFrame({"from": m_from, "to": m_to}).drop_duplicates(subset=["from"], keep="first")
    # Warn on duplicate target IDs
    dup_to = dedup[dedup.duplicated(subset=["to"], keep=False)]["to"].unique()
    if len(dup_to) > 0:
        print(
            f"Warning: {len(dup_to)} target sample IDs in '{args.map_to_col}' are duplicated in mapping: "
            f"{', '.join(list(map(str, dup_to))[:10])}{' ...' if len(dup_to) > 10 else ''}",
            file=sys.stderr,
        )
    rename_map = dict(zip(dedup["from"], dedup["to"]))
    # Ensure TPM column names are normalized for matching
    tpm_df.columns = tpm_df.columns.astype(str).str.strip()
    before_cols = set(tpm_df.columns)
    tpm_df = tpm_df.rename(columns=rename_map)
    after_cols = set(tpm_df.columns)
    renamed_count = len(after_cols - before_cols)
    if renamed_count == 0:
        print(
            "Warning: No TPM columns were renamed via the sample map. Check that TPM column names match the '" \
            f"{args.map_from_col}' values.",
            file=sys.stderr,
        )
    # Reindex columns to match desired sample order from the raw matrix
    tpm_df = tpm_df.reindex(columns=desired_samples)

    # Report any desired samples missing after mapping
    missing_samples = [s for s in desired_samples if s not in tpm_df.columns]
    if missing_samples:
        print(f"Warning: {len(missing_samples)} desired samples not present in TPM after mapping: "
              f"{', '.join(missing_samples[:10])}{' ...' if len(missing_samples) > 10 else ''}", file=sys.stderr)

    # Resolve and collect rows in the order of genes list
    records: List[pd.Series] = []
    missing_rows: List[Dict[str, str]] = []  # type: ignore[name-defined]
    sample_cols = list(tpm_df.columns)

    for sym in genes:
        gid, reason = resolve_gene_id(sym, sym2gid)
        if gid is None:
            missing_rows.append({"GeneSymbol": sym, "Reason": reason or "not-found", "GeneID": ""})
            continue
        if gid not in tpm_df.index:
            missing_rows.append({"GeneSymbol": sym, "Reason": "id-not-in-tpm", "GeneID": gid})
            continue
        row = tpm_df.loc[gid]
        # Ensure alignment to sample_cols and coerce to numeric where possible (NaN for missing samples)
        row = row.reindex(sample_cols)
        try:
            row_num = pd.to_numeric(row, errors="coerce")
        except Exception:
            row_num = row
        row_num.name = sym  # label with symbol
        records.append(row_num)

    if len(records) == 0:
        print("No genes resolved to TPM entries. Nothing to write.", file=sys.stderr)
        # Still write missing report if requested
        if args.missing_out is not None and len(missing_rows):
            pd.DataFrame(missing_rows).to_csv(args.missing_out, index=False)
        return 1

    out_df = pd.DataFrame(records)
    out_df.insert(0, "GeneSymbol", out_df.index)
    out_df.index = pd.Index(range(len(out_df)))

    out_path = args.output or default_output_path(args.raw_genes)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)

    if args.missing_out is not None:
        pd.DataFrame(missing_rows).to_csv(args.missing_out, index=False)

    # Basic summary
    resolved = len(records)
    missing = len(missing_rows)
    print(f"Wrote {resolved} genes to {out_path} ({missing} missing)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
