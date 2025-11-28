import argparse
import sys
import pandas as pd


def bed9_to_csv(bed_file_path: str, csv_file_path: str, sep: str = '\t', col_names: list | None = None):
    """Convert a BED-like file to CSV.

    The function will read the input with no header and assign column names based on
    the detected number of columns or the explicit `col_names` argument. If the
    number of columns does not match common expectations it will create generic
    column names `col0`, `col1`, ...

    Args:
        bed_file_path: Path to input BED-like file (tab/sep-delimited).
        csv_file_path: Path where CSV will be written.
        sep: Input field separator (default: '\t').
        col_names: Optional list of column names to assign. If provided, its length
                   must match the number of columns in the file.
    """
    try:
        df = pd.read_csv(bed_file_path, sep=sep, header=None, comment='#', dtype=str)
    except FileNotFoundError:
        print(f"Error: File not found at '{bed_file_path}'", file=sys.stderr)
        raise
    except Exception as e:
        print(f"Error reading '{bed_file_path}': {e}", file=sys.stderr)
        raise

    ncols = df.shape[1]
    # If user supplied column names, use them (but validate length)
    if col_names is not None:
        if len(col_names) != ncols:
            raise ValueError(f"Provided column names length ({len(col_names)}) does not match file columns ({ncols})")
        df.columns = col_names
    else:
        # Heuristics for common BED-like formats
        if ncols == 9:
            df.columns = ['Ensembl_ID', 'V2', 'V3','V5' 'V4', 'V6', 'gene_biotype', 'Gene_Symbol', 'Entrez_ID']
        elif ncols == 8:
            # legacy/alternate format used in this project
            df.columns = ['Ensembl_ID', 'V2', 'V3', 'V4', 'V6', 'gene_biotype', 'Gene_Symbol', 'Entrez_ID']
        else:
            # generic names
            df.columns = [f'col{i}' for i in range(ncols)]

    try:
        df.to_csv(csv_file_path, index=False)
    except Exception as e:
        print(f"Error writing CSV to '{csv_file_path}': {e}", file=sys.stderr)
        raise

    print(f"Successfully converted '{bed_file_path}' ({ncols} columns) to '{csv_file_path}'")


def _parse_args(argv=None):
    p = argparse.ArgumentParser(description='Convert a BED-like file to CSV (tab/sep-delimited).')
    p.add_argument('bed', help='Input BED-like file path')
    p.add_argument('csv', help='Output CSV file path')
    p.add_argument('--sep', default='\t', help="Input separator (default: '\t')")
    p.add_argument('--cols', default=None, help='Comma-separated list of column names to set (must match file columns)')
    return p.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)
    cols = None
    if args.cols:
        cols = [c for c in args.cols.split(',')]
    try:
        bed9_to_csv(args.bed, args.csv, sep=args.sep, col_names=cols)
    except Exception as e:
        print(f"Conversion failed: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())