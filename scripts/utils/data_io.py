"""
scripts/utils/data_io.py
========================
Low-memory, format-aware I/O utilities for lncSignal-AD.

All loaders avoid reading large files entirely into memory when possible.
They handle:
  - GEO Series Matrix format (.txt.gz)   — microarray expression + metadata
  - Compressed TSV (.tsv.gz)             — RNA-seq counts, NCBI annotations
  - Plain tab-separated annotation files — GPL platform tables
  - GENCODE GTF (.gtf.gz)                — skipping # comment header lines

Scientific warnings
-------------------
- Series-matrix expression values are PRE-NORMALIZED microarray intensities.
  They must NOT be passed to PyDESeq2 (which expects raw integer counts).
- GSE125583 raw counts ARE integer counts suitable for PyDESeq2.
- Do NOT mix the two scales in a single matrix.
- Feature IDs differ across platforms; never join on position alone.
"""

from __future__ import annotations

import csv
import gzip
import io
from pathlib import Path
from typing import Generator, Iterator

import pandas as pd

from scripts.utils.paths import resolve_path


# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------
PathLike = str | Path


# ---------------------------------------------------------------------------
# GEO Series Matrix loaders
# ---------------------------------------------------------------------------

def load_series_matrix_header(path: PathLike) -> dict[str, list[str]]:
    """
    Parse the ``!``-prefixed header section of a GEO Series Matrix file
    **without** loading the expression table.

    Parameters
    ----------
    path : str or Path
        Path to a ``*_series_matrix.txt.gz`` file.

    Returns
    -------
    dict[str, list[str]]
        Mapping from field name (e.g. ``"Series_geo_accession"``) to a list
        of its values (most fields have one value; ``Sample_characteristics``
        has one row per characteristics key with multi-sample values).

    Notes
    -----
    - Stops at ``!series_matrix_table_begin`` so the expression table is
      never read.
    - Tab-separated quoted values are stripped of surrounding double-quotes.
    """
    path = resolve_path(path)
    header: dict[str, list[str]] = {}
    opener = gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else open(path, encoding="utf-8", errors="replace")

    with opener as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n\r")
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!"):
                continue
            parts = line.split("\t")
            key = parts[0].lstrip("!").strip()
            values = [v.strip().strip('"') for v in parts[1:]]
            if key in header:
                header[key].append(values)
            else:
                header[key] = [values]
    return header


def load_series_matrix_expression(
    path: PathLike,
    *,
    dtype: str = "float32",
    nrows: int | None = None,
) -> pd.DataFrame:
    """
    Load the expression data table from a GEO Series Matrix file.

    The table begins after ``!series_matrix_table_begin`` and ends before
    ``!series_matrix_table_end``.  It has the layout::

        "ID_REF"  "GSM..."  "GSM..."  ...
        "ILMN_1"  12.3      11.7      ...

    Parameters
    ----------
    path : str or Path
        Path to ``*_series_matrix.txt.gz``.
    dtype : str, optional
        NumPy dtype for expression values.  Default ``"float32"`` saves memory.
    nrows : int or None, optional
        Maximum number of data rows to load.  ``None`` loads all.

    Returns
    -------
    pd.DataFrame
        Expression matrix with shape (n_probes, n_samples).
        Index name: ``"feature_id"`` (probe IDs, e.g. ``ILMN_*``).
        Column names: GSM accession strings.

    Raises
    ------
    ValueError
        If the ``!series_matrix_table_begin`` sentinel is not found.
    """
    path = resolve_path(path)
    opener = gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else open(path, encoding="utf-8", errors="replace")

    rows: list[dict] = []
    in_table = False
    n_read = 0

    with opener as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n\r")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                # Next line is the header
                header_line = next(fh).rstrip("\n\r")
                cols = [c.strip().strip('"') for c in header_line.split("\t")]
                # cols[0] = "ID_REF", cols[1:] = GSM IDs
                continue
            if not in_table:
                continue
            if line.startswith("!series_matrix_table_end") or line == "":
                break
            parts = line.split("\t")
            probe_id = parts[0].strip().strip('"')
            values = parts[1:]
            row: dict = {"feature_id": probe_id}
            for col, val in zip(cols[1:], values):
                v = val.strip().strip('"')
                try:
                    row[col] = float(v)
                except ValueError:
                    row[col] = float("nan")
            rows.append(row)
            n_read += 1
            if nrows is not None and n_read >= nrows:
                break

    if not rows and not in_table:
        raise ValueError(
            f"'!series_matrix_table_begin' sentinel not found in {path}.\n"
            "The file may be corrupt or not a GEO Series Matrix."
        )

    df = pd.DataFrame(rows).set_index("feature_id")
    df = df.astype(dtype)
    return df


# ---------------------------------------------------------------------------
# Compressed TSV loader
# ---------------------------------------------------------------------------

def load_tsv_gz(
    path: PathLike,
    *,
    index_col: int | str | None = 0,
    usecols: list[str | int] | None = None,
    dtype: str | dict | None = None,
    comment: str | None = None,
    nrows: int | None = None,
) -> pd.DataFrame:
    """
    Load a gzip-compressed TSV file (e.g. RNA-seq count matrix or NCBI
    annotation table) into a DataFrame.

    Parameters
    ----------
    path : str or Path
    index_col : int or str or None
        Column to use as the row index.
    usecols : list or None
        Subset of column names or indices to load.
    dtype : str or dict or None
        Data type(s) for one or more columns.
    comment : str or None
        Lines starting with this character are skipped.
    nrows : int or None
        Maximum number of rows to read.

    Returns
    -------
    pd.DataFrame
    """
    path = resolve_path(path)
    kwargs: dict = {
        "sep": "\t",
        "index_col": index_col,
        "usecols": usecols,
        "dtype": dtype,
        "comment": comment,
        "nrows": nrows,
    }
    # Remove None values to let pandas use defaults
    kwargs = {k: v for k, v in kwargs.items() if v is not None}

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        return pd.read_csv(fh, **kwargs)


# ---------------------------------------------------------------------------
# Plain (non-gzip) tab-separated annotation loader
# ---------------------------------------------------------------------------

def load_plain_tsv(
    path: PathLike,
    *,
    comment_prefix: str = "#",
    index_col: int | str | None = None,
    usecols: list[str | int] | None = None,
    nrows: int | None = None,
) -> pd.DataFrame:
    """
    Load a plain-text tab-separated file, skipping ``#``-prefixed comment
    lines at the top.  Used for Illumina GPL platform annotation files
    (``GPL10558-50081.txt``, ``GPL6947-13512.txt``, ``GPL15988.txt``).

    Parameters
    ----------
    path : str or Path
    comment_prefix : str
        Lines starting with this prefix are skipped.
    index_col, usecols, nrows : forwarded to pd.read_csv.

    Returns
    -------
    pd.DataFrame
    """
    path = resolve_path(path)
    kwargs: dict = {
        "sep": "\t",
        "index_col": index_col,
        "usecols": usecols,
        "nrows": nrows,
        "comment": comment_prefix if len(comment_prefix) == 1 else None,
    }
    kwargs = {k: v for k, v in kwargs.items() if v is not None}
    df = pd.read_csv(path, encoding="utf-8", encoding_errors="replace", **kwargs)
    return df


# ---------------------------------------------------------------------------
# GTF loader (header peek only — full GTF is huge)
# ---------------------------------------------------------------------------

def load_gtf_gz_header(path: PathLike, n_comment_lines: int = 10) -> list[str]:
    """
    Return the first ``n_comment_lines`` ``##``-prefixed header lines from a
    gzip-compressed GTF file without reading the data records.

    Parameters
    ----------
    path : str or Path
    n_comment_lines : int
        Maximum number of comment lines to collect.

    Returns
    -------
    list[str]
        Comment lines (with ``##`` prefix retained).
    """
    path = resolve_path(path)
    comments: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            stripped = line.rstrip("\n\r")
            if stripped.startswith("#"):
                comments.append(stripped)
                if len(comments) >= n_comment_lines:
                    break
            else:
                break
    return comments


def iter_gtf_gz_records(
    path: PathLike,
    feature_type: str | None = None,
) -> Generator[dict[str, str], None, None]:
    """
    Iterate over records in a gzip-compressed GTF file without loading it
    entirely into memory.

    Parameters
    ----------
    path : str or Path
    feature_type : str or None
        If given (e.g. ``"gene"`` or ``"transcript"``), only yield records
        whose feature column matches.

    Yields
    ------
    dict
        Keys: ``seqname``, ``source``, ``feature``, ``start``, ``end``,
        ``score``, ``strand``, ``frame``, ``attributes`` (raw string).
    """
    path = resolve_path(path)
    cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n\r").split("\t", 8)
            if len(parts) < 9:
                continue
            record = dict(zip(cols, parts))
            if feature_type is not None and record["feature"] != feature_type:
                continue
            yield record
