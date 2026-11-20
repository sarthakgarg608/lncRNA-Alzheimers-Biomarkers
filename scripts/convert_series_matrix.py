"""
convert_series_matrix.py
========================
Converts GEO Series Matrix .txt.gz files into easy-to-use formats under data/processed/.

For each dataset produces:
  {GSE}_meta.csv        -- clinical / sample metadata (one row per sample)
  {GSE}_expr.parquet    -- expression matrix (probes x samples), fast binary format
  {GSE}_expr.csv.gz     -- same as CSV gzipped (opens in Excel, R, Python, etc.)

GSE125583 is special: its expression lives in the raw counts TSV, not the matrix.
"""

import gzip
import re
import sys
from io import StringIO
from pathlib import Path

import pandas as pd

# Force UTF-8 output so Unicode arrows/symbols don't crash on Windows
sys.stdout.reconfigure(encoding="utf-8")

RAW = Path("data/raw")
OUT = Path("data/processed")
OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def parse_series_matrix(path):
    """
    Parse a GEO series matrix file.

    Returns
    -------
    meta_df : pd.DataFrame
        Sample metadata (one row per sample, one column per field).
    expr_df : pd.DataFrame or None
        Expression table (probes as index, GSM IDs as columns).
        None if the table section is empty.
    """
    sample_meta = {}
    table_lines = []
    in_table = False

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")

            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                in_table = False
                continue

            if in_table:
                table_lines.append(line)
            elif line.startswith("!Sample_"):
                parts = line.split("\t")
                key = parts[0].lstrip("!")          # e.g. "Sample_title"
                vals = [v.strip('"') for v in parts[1:]]
                if key not in sample_meta:
                    sample_meta[key] = vals
                else:
                    # Repeated key (e.g. Sample_characteristics_ch1) -> numbered variant
                    n = sum(1 for k in sample_meta
                            if k == key or re.match(r"^" + re.escape(key) + r"_\d+$", k))
                    sample_meta[key + "_" + str(n + 1)] = vals

    # Build metadata DataFrame
    meta_df = pd.DataFrame(sample_meta)

    # Rename characteristics columns to human-readable names
    char_cols = [c for c in meta_df.columns if c.startswith("Sample_characteristics_ch1")]
    rename_map = {}
    for col in char_cols:
        for val in meta_df[col].dropna():
            if ":" in val:
                rename_map[col] = val.split(":")[0].strip().lower().replace(" ", "_")
                break
    meta_df = meta_df.rename(columns=rename_map)

    # Expression matrix
    if not table_lines:
        return meta_df, None

    expr_df = pd.read_csv(
        StringIO("\n".join(table_lines)),
        sep="\t",
        index_col=0,
        dtype=str,
    )
    for col in expr_df.columns:
        expr_df[col] = pd.to_numeric(expr_df[col], errors="coerce")

    return meta_df, expr_df


# ---------------------------------------------------------------------------
# datasets
# ---------------------------------------------------------------------------

DATASETS = ["GSE125583", "GSE140829", "GSE63060", "GSE63061"]


def convert_all():
    for gse in DATASETS:
        matrix_path = RAW / (gse + "_series_matrix.txt.gz")
        if not matrix_path.exists():
            print("[SKIP] " + gse + " -- series matrix not found")
            continue

        print("")
        print("=" * 60)
        print("Processing " + gse + " ...")

        meta_df, expr_df = parse_series_matrix(matrix_path)

        # --- Metadata CSV ---
        meta_out = OUT / (gse + "_meta.csv")
        meta_df.to_csv(meta_out, index=False)
        print("  [meta] " + str(meta_df.shape[0]) + " samples x " + str(meta_df.shape[1]) + " fields  ->  " + meta_out.name)

        # --- Expression ---
        if expr_df is not None and not expr_df.empty:
            parquet_out = OUT / (gse + "_expr.parquet")
            csv_out     = OUT / (gse + "_expr.csv.gz")
            expr_df.to_parquet(parquet_out)
            expr_df.to_csv(csv_out, compression="gzip")
            print("  [expr] " + str(expr_df.shape[0]) + " probes x " + str(expr_df.shape[1]) + " samples")
            print("         Parquet  -> " + parquet_out.name)
            print("         CSV.gz   -> " + csv_out.name)
        else:
            print("  [expr] No table in series matrix.")
            if gse == "GSE125583":
                raw_counts = RAW / "GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz"
                if raw_counts.exists():
                    print("         Reading raw counts from " + raw_counts.name + " ...")
                    counts_df = pd.read_csv(raw_counts, sep="\t", index_col=0, compression="gzip")
                    parquet_out = OUT / (gse + "_expr.parquet")
                    csv_out     = OUT / (gse + "_expr.csv.gz")
                    counts_df.to_parquet(parquet_out)
                    counts_df.to_csv(csv_out, compression="gzip")
                    print("         " + str(counts_df.shape[0]) + " genes x " + str(counts_df.shape[1]) + " samples")
                    print("         Parquet  -> " + parquet_out.name)
                    print("         CSV.gz   -> " + csv_out.name)
                else:
                    print("         Raw counts file not found -- skipping.")

    print("")
    print("=" * 60)
    print("All done! Files written to data/processed/:")
    print("")
    for f in sorted(OUT.iterdir()):
        size_mb = f.stat().st_size / 1e6
        print("  " + f.name.ljust(48) + str(round(size_mb, 1)).rjust(6) + " MB")


if __name__ == "__main__":
    convert_all()
