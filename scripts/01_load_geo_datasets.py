"""
scripts/01_load_geo_datasets.py
================================
Stage 08 — Load all four GEO datasets from local files.

Datasets processed
------------------
  GSE63061  (blood, Illumina HT-12 v4, GPL10558) — primary discovery
  GSE63060  (blood, Illumina HT-12 v3, GPL6947)  — replication
  GSE140829 (blood, Illumina HT-12 v4, GPL15988) — LOCKED external validation
  GSE125583 (brain RNA-seq, GPL16791)             — biological support

Outputs per dataset (under data/interim/<GSEID>/)
-------------------------------------------------
  expression.parquet   — features × samples, float32  (microarray)
  raw_counts.parquet   — features × samples, int32     (RNA-seq only)
  metadata.csv         — per-sample canonical metadata
  technical_metadata.csv — GSE140829 only; NO diagnosis column
  parsing_report.json  — shape, dtypes, NA counts, notes

Scientific integrity rules enforced
------------------------------------
- GSE140829 diagnosis labels are NEVER extracted or stored.
- MCI and Unknown are preserved; never merged into Control or AD.
- Microarray expression (float) and RNA-seq counts (int) are never mixed.
- No normalization, differential expression, or ML in this script.
- Original probe/gene IDs are preserved exactly.

Usage
-----
  .venv\\Scripts\\python.exe scripts/01_load_geo_datasets.py
"""

from __future__ import annotations

import gzip
import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Project root and utility imports
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from scripts.utils.paths import get_project_root, load_dataset_config, resolve_path
from scripts.utils.data_io import (
    load_series_matrix_header,
    load_series_matrix_expression,
    load_tsv_gz,
)
from scripts.utils.validation import (
    check_duplicate_ids,
    check_expression_numeric,
    check_sample_id_alignment,
    ValidationResult,
)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Canonical diagnosis vocabulary
# ---------------------------------------------------------------------------
VALID_DIAGNOSIS_LABELS = frozenset({"AD", "Control", "MCI", "Unknown"})


# ---------------------------------------------------------------------------
# Diagnosis mapping helpers
# ---------------------------------------------------------------------------

def _map_status_label(raw_value: str) -> str:
    """
    Map a raw GEO 'status' or 'diagnosis' string to the canonical vocabulary.

    Canonical labels: "AD", "Control", "MCI", "Unknown"
    MCI is never merged into Control.  Unknown is never discarded.

    Parameters
    ----------
    raw_value : str
        Stripped string value extracted from !Sample_characteristics_ch1.

    Returns
    -------
    str
        One of "AD", "Control", "MCI", "Unknown".
    """
    normed = raw_value.strip().lower()
    if normed in {"alzheimer's disease", "alzheimers disease", "alzheimer disease",
                  "ad", "alzheimer"}:
        return "AD"
    if normed in {"control", "ctl", "hc", "healthy control", "normal", "healthy"}:
        return "Control"
    if normed in {"mci", "mild cognitive impairment"}:
        return "MCI"
    return "Unknown"


# ---------------------------------------------------------------------------
# Series Matrix metadata parser
# ---------------------------------------------------------------------------

def parse_series_matrix_metadata(
    header: dict[str, list],
    *,
    diagnosis_field: str | None = "status",
    age_field: str | None = "age",
    sex_field: str | None = "gender",
    tissue_field: str | None = "tissue",
    ethnicity_field: str | None = "ethnicity",
    batch_field: str | None = "batch",
    locked: bool = False,
) -> pd.DataFrame:
    """
    Convert the parsed Series Matrix header dict into a per-sample metadata DataFrame.

    The !Sample_characteristics_ch1 field in a GEO Series Matrix has this layout:
        Each row is a tab-separated list of "key: value" strings, one per sample.
        Multiple rows = multiple characteristics.

    Parameters
    ----------
    header : dict
        Output of load_series_matrix_header().
    diagnosis_field : str or None
        The characteristic key containing diagnosis labels (e.g. "status").
    age_field, sex_field, etc. : str or None
        Characteristic keys for each metadata variable.
    locked : bool
        If True (GSE140829), do NOT extract or store diagnosis labels.

    Returns
    -------
    pd.DataFrame
        One row per sample with canonical column schema.
    """
    # Sample IDs (geo_accession per sample)
    sample_ids_raw = header.get("Sample_geo_accession", [[]])
    sample_ids: list[str] = []
    for row in sample_ids_raw:
        if isinstance(row, list):
            sample_ids.extend(row)
        elif isinstance(row, str):
            sample_ids.append(row)
    sample_ids = [s for s in sample_ids if s]

    n_samples = len(sample_ids)
    if n_samples == 0:
        raise ValueError("No sample IDs found in header (Sample_geo_accession missing).")

    # Build a dict: characteristic_key -> list of values, one per sample
    char_rows: list[list[str]] = []
    for entry in header.get("Sample_characteristics_ch1", []):
        if isinstance(entry, list):
            char_rows.append(entry)
        elif isinstance(entry, str):
            char_rows.append([entry])

    # Pivot: build a dict keyed by lowercase characteristic name
    char_dict: dict[str, list[str]] = {}
    for row in char_rows:
        if not row:
            continue
        # Each element of `row` is "key: value" for one sample
        # Parse the key from the first non-empty element
        first_nonempty = next((v for v in row if v.strip()), None)
        if not first_nonempty or ":" not in first_nonempty:
            continue
        key_part = first_nonempty.split(":", 1)[0].strip().lower()
        values: list[str] = []
        for v in row:
            if v.strip() and ":" in v:
                values.append(v.split(":", 1)[1].strip())
            else:
                values.append("")
        # Pad or truncate to n_samples
        while len(values) < n_samples:
            values.append("")
        values = values[:n_samples]
        char_dict[key_part] = values

    def _get_col(field_name: str | list[str] | tuple[str, ...] | None) -> list[str | None]:
        if field_name is None:
            return [None] * n_samples
        if isinstance(field_name, (list, tuple)):
            for fn in field_name:
                res = char_dict.get(fn.lower())
                if res and any(v.strip() for v in res if v):
                    return res
            return [None] * n_samples
        return char_dict.get(field_name.lower(), [None] * n_samples)

    # Diagnosis — only for unlocked datasets
    if locked:
        diagnosis_col: list[str | None] = [None] * n_samples
    else:
        raw_diag = _get_col(diagnosis_field)
        diagnosis_col = [
            _map_status_label(v) if v else "Unknown"
            for v in raw_diag
        ]

    # Sample titles
    titles_raw = header.get("Sample_title", [[]])
    titles: list[str] = []
    for row in titles_raw:
        if isinstance(row, list):
            titles.extend(row)
        elif isinstance(row, str):
            titles.append(row)
    while len(titles) < n_samples:
        titles.append("")
    titles = titles[:n_samples]

    # Platform IDs per sample
    platform_raw = header.get("Sample_platform_id", [[]])
    platforms: list[str] = []
    for row in platform_raw:
        if isinstance(row, list):
            platforms.extend(row)
        elif isinstance(row, str):
            platforms.append(row)
    while len(platforms) < n_samples:
        platforms.append("")
    platforms = platforms[:n_samples]

    # Build DataFrame
    rows: list[dict[str, Any]] = []
    age_vals = _get_col(age_field)
    sex_vals = _get_col(sex_field)
    tissue_vals = _get_col(tissue_field)
    ethnicity_vals = _get_col(ethnicity_field)
    batch_vals = _get_col(batch_field)
    if all(v is None or v == "" for v in batch_vals):
        # Fallback: check Sample_description for multi-center cohorts (e.g. AddNeuroMed: LDZ, PRG, KPO...)
        desc_raw = header.get("Sample_description", [[]])
        descriptions: list[str] = []
        for d_entry in desc_raw:
            if isinstance(d_entry, list):
                descriptions.extend(d_entry)
            elif isinstance(d_entry, str):
                descriptions.append(d_entry)
        if len(descriptions) >= n_samples:
            import re
            extracted_centers = []
            for d in descriptions[:n_samples]:
                m = re.match(r"^([A-Z0-9]{3,4})", str(d).strip().strip('"'))
                extracted_centers.append(m.group(1)[:3] if m else None)
            if any(extracted_centers):
                batch_vals = extracted_centers

    for i, sid in enumerate(sample_ids):
        row: dict[str, Any] = {
            "sample_id": sid,
            "title": titles[i] if i < len(titles) else "",
            "platform_id": platforms[i] if i < len(platforms) else "",
            "age": age_vals[i] if i < len(age_vals) else None,
            "sex": sex_vals[i] if i < len(sex_vals) else None,
            "tissue": tissue_vals[i] if i < len(tissue_vals) else None,
            "ethnicity": ethnicity_vals[i] if i < len(ethnicity_vals) else None,
            "batch": batch_vals[i] if i < len(batch_vals) else None,
        }
        if not locked:
            row["diagnosis"] = diagnosis_col[i]
        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# GSE125583 metadata extraction
# ---------------------------------------------------------------------------

def parse_gse125583_metadata(header: dict[str, list]) -> pd.DataFrame:
    """
    Extract metadata for GSE125583 (brain RNA-seq).

    Diagnosis labels come from sample titles ("AD ##-##" or "CON ##-##")
    or from characteristics_ch1 "disease state" field.

    Returns
    -------
    pd.DataFrame with canonical schema columns.
    """
    # Reuse generic parser with brain-specific field names (both 'diagnosis' and 'disease state', 'age')
    df = parse_series_matrix_metadata(
        header,
        diagnosis_field=["diagnosis", "disease state"],
        age_field=["age", "age at death"],
        sex_field=["sex", "gender"],
        tissue_field="tissue",
        ethnicity_field=None,
        batch_field=None,
        locked=False,
    )

    # Map brain-specific diagnosis strings
    if "diagnosis" in df.columns:
        def _map_brain(val: str | None) -> str:
            if val is None:
                return "Unknown"
            v = str(val).strip().lower()
            if v in {"alzheimer's disease", "alzheimers disease", "alzheimer disease", "ad"}:
                return "AD"
            if v in {"neurologically normal", "normal", "control", "ctl"}:
                return "Control"
            if v in {"mci", "mild cognitive impairment"}:
                return "MCI"
            return "Unknown"

        df["diagnosis"] = df["diagnosis"].apply(_map_brain)

    # Fallback: parse title if all Unknown
    if (df["diagnosis"] == "Unknown").all() or "diagnosis" not in df.columns:
        log.warning("GSE125583: No diagnosis from characteristics; parsing from titles.")
        import re
        def _diag_from_title(title: str) -> str:
            t = str(title).strip().upper()
            if re.match(r"^AD\s", t) or t.startswith("AD_") or t.startswith("AD-"):
                return "AD"
            if re.match(r"^(CON|CTRL|CONTROL)\s", t):
                return "Control"
            return "Unknown"
        df["diagnosis"] = df["title"].apply(_diag_from_title)

    return df


# ---------------------------------------------------------------------------
# Save helpers
# ---------------------------------------------------------------------------

def _save_expression(expr: pd.DataFrame, out_dir: Path, fname: str = "expression.parquet") -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    expr.to_parquet(out_dir / fname, engine="pyarrow", compression="snappy")
    log.info("  Saved %s: shape %s", fname, expr.shape)


def _save_metadata(meta: pd.DataFrame, out_dir: Path, fname: str = "metadata.csv") -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    meta.to_csv(out_dir / fname, index=False)
    log.info("  Saved %s: %d samples, columns: %s", fname, len(meta), list(meta.columns))


def _save_report(report: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "parsing_report.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)
    log.info("  Saved parsing_report.json")


# ---------------------------------------------------------------------------
# Individual dataset loaders
# ---------------------------------------------------------------------------

def load_gse63061(cfg: dict, root: Path) -> dict:
    """Load GSE63061 — primary discovery cohort (blood, GPL10558)."""
    gse_id = "GSE63061"
    log.info("=" * 60)
    log.info("Loading %s (primary discovery)", gse_id)
    log.info("=" * 60)

    dcfg = cfg["datasets"][gse_id]
    matrix_path = root / dcfg["series_matrix_path"]
    out_dir = root / "data" / "interim" / gse_id

    # --- Expression matrix ---
    log.info("  Parsing expression matrix...")
    expr = load_series_matrix_expression(matrix_path, dtype="float32")
    expr.index.name = "feature_id"
    log.info("  Expression shape (features × samples): %s", expr.shape)

    # --- Metadata ---
    log.info("  Parsing metadata...")
    header = load_series_matrix_header(matrix_path)
    meta = parse_series_matrix_metadata(
        header,
        diagnosis_field="status",
        age_field="age",
        sex_field="gender",
        tissue_field="tissue",
        ethnicity_field="ethnicity",
        batch_field=None,
        locked=False,
    )
    meta.insert(1, "dataset_id", gse_id)

    # --- Validation ---
    dup_check = check_duplicate_ids(meta["sample_id"].tolist())
    align_check = check_sample_id_alignment(meta["sample_id"].tolist(), list(expr.columns))
    num_check = check_expression_numeric(expr)
    _log_validation(gse_id, [dup_check, align_check, num_check])

    # Diagnosis counts
    diag_counts = meta["diagnosis"].value_counts(dropna=False).to_dict()
    log.info("  Diagnosis counts: %s", diag_counts)

    # Probe IDs
    probe_prefix = "ILMN_"
    ilmn_count = int(expr.index.str.startswith(probe_prefix).sum())
    log.info("  Probes with ILMN_ prefix: %d / %d", ilmn_count, len(expr))

    # NA count
    na_count = int(expr.isna().sum().sum())
    na_pct = 100 * na_count / (expr.shape[0] * expr.shape[1])
    log.info("  Missing expression values: %d (%.3f%%)", na_count, na_pct)

    # --- Save ---
    _save_expression(expr, out_dir)
    _save_metadata(meta, out_dir)

    report = {
        "dataset_id": gse_id,
        "role": dcfg["role"],
        "platform_id": dcfg["platform_id"],
        "tissue": dcfg["tissue"],
        "n_samples": expr.shape[1],
        "n_features": expr.shape[0],
        "feature_id_type": dcfg["feature_id_type"],
        "expression_scale": dcfg["expression_scale"],
        "expression_dtype": str(expr.dtypes.iloc[0]),
        "na_values": na_count,
        "na_pct": round(na_pct, 4),
        "ilmn_prefix_probes": ilmn_count,
        "diagnosis_counts": diag_counts,
        "validation": {
            "duplicate_sample_ids": dup_check.passed,
            "sample_id_alignment": align_check.passed,
            "expression_numeric": num_check.passed,
        },
        "notes": "Primary discovery. Labels unlocked. MCI preserved as separate class.",
    }
    _save_report(report, out_dir)
    return report


def load_gse63060(cfg: dict, root: Path) -> dict:
    """Load GSE63060 — replication cohort (blood, GPL6947)."""
    gse_id = "GSE63060"
    log.info("=" * 60)
    log.info("Loading %s (replication cohort)", gse_id)
    log.info("=" * 60)

    dcfg = cfg["datasets"][gse_id]
    matrix_path = root / dcfg["series_matrix_path"]
    out_dir = root / "data" / "interim" / gse_id

    log.info("  Parsing expression matrix...")
    expr = load_series_matrix_expression(matrix_path, dtype="float32")
    expr.index.name = "feature_id"
    log.info("  Expression shape (features × samples): %s", expr.shape)

    log.info("  Parsing metadata...")
    header = load_series_matrix_header(matrix_path)
    meta = parse_series_matrix_metadata(
        header,
        diagnosis_field="status",
        age_field="age",
        sex_field="gender",
        tissue_field="tissue",
        ethnicity_field="ethnicity",
        batch_field=None,
        locked=False,
    )
    meta.insert(1, "dataset_id", gse_id)

    dup_check = check_duplicate_ids(meta["sample_id"].tolist())
    align_check = check_sample_id_alignment(meta["sample_id"].tolist(), list(expr.columns))
    num_check = check_expression_numeric(expr)
    _log_validation(gse_id, [dup_check, align_check, num_check])

    diag_counts = meta["diagnosis"].value_counts(dropna=False).to_dict()
    log.info("  Diagnosis counts: %s", diag_counts)

    ilmn_count = int(expr.index.str.startswith("ILMN_").sum())
    log.info("  Probes with ILMN_ prefix: %d / %d", ilmn_count, len(expr))

    na_count = int(expr.isna().sum().sum())
    na_pct = 100 * na_count / (expr.shape[0] * expr.shape[1])
    log.info("  Missing expression values: %d (%.3f%%)", na_count, na_pct)

    _save_expression(expr, out_dir)
    _save_metadata(meta, out_dir)

    report = {
        "dataset_id": gse_id,
        "role": dcfg["role"],
        "platform_id": dcfg["platform_id"],
        "tissue": dcfg["tissue"],
        "n_samples": expr.shape[1],
        "n_features": expr.shape[0],
        "feature_id_type": dcfg["feature_id_type"],
        "expression_scale": dcfg["expression_scale"],
        "expression_dtype": str(expr.dtypes.iloc[0]),
        "na_values": na_count,
        "na_pct": round(na_pct, 4),
        "ilmn_prefix_probes": ilmn_count,
        "diagnosis_counts": diag_counts,
        "validation": {
            "duplicate_sample_ids": dup_check.passed,
            "sample_id_alignment": align_check.passed,
            "expression_numeric": num_check.passed,
        },
        "notes": (
            "Replication cohort. Different Illumina HT-12 v3 platform (GPL6947) vs "
            "GSE63061 (GPL10558). Probes are NOT identical across platforms. "
            "Do not merge probe sets without mapping. Labels unlocked. MCI preserved."
        ),
    }
    _save_report(report, out_dir)
    return report


def load_gse140829(cfg: dict, root: Path) -> dict:
    """Load GSE140829 — LOCKED external validation (blood, GPL15988).

    IMPORTANT: Diagnosis labels are NEVER extracted or stored.
    Only technical/structural metadata is saved.
    """
    gse_id = "GSE140829"
    log.info("=" * 60)
    log.info("Loading %s (LOCKED external validation)", gse_id)
    log.info("IMPORTANT: Diagnosis labels will NOT be extracted.", )
    log.info("=" * 60)

    dcfg = cfg["datasets"][gse_id]
    matrix_path = root / dcfg["series_matrix_path"]
    out_dir = root / "data" / "interim" / gse_id

    log.info("  Parsing expression matrix (this may take a minute)...")
    expr = load_series_matrix_expression(matrix_path, dtype="float32")
    expr.index.name = "feature_id"
    log.info("  Expression shape (features × samples): %s", expr.shape)

    log.info("  Parsing structural metadata (NO diagnosis)...")
    header = load_series_matrix_header(matrix_path)

    # Parse ONLY technical/structural metadata — locked=True prevents diagnosis extraction
    tech_meta = parse_series_matrix_metadata(
        header,
        diagnosis_field="diagnosis",   # specified but ignored because locked=True
        age_field="age_at_draw",
        sex_field="Sex",
        tissue_field=None,
        ethnicity_field=None,
        batch_field="batch",
        locked=True,                    # LOCKED — no diagnosis column
    )
    tech_meta.insert(1, "dataset_id", gse_id)
    # Double-check: ensure no diagnosis column slipped through
    assert "diagnosis" not in tech_meta.columns, (
        "SAFETY VIOLATION: 'diagnosis' column found in GSE140829 technical_metadata!"
    )

    dup_check = check_duplicate_ids(tech_meta["sample_id"].tolist())
    num_check = check_expression_numeric(expr)
    _log_validation(gse_id, [dup_check, num_check])

    # Probe ID inspection
    probe_sample = expr.index[:5].tolist()
    log.info("  First 5 probe IDs: %s", probe_sample)
    ilmn_count = int(expr.index.str.startswith("ILMN_").sum())
    log.info("  Probes with ILMN_ prefix: %d / %d (nuID-style expected)", ilmn_count, len(expr))

    na_count = int(expr.isna().sum().sum())
    na_pct = 100 * na_count / (expr.shape[0] * expr.shape[1])
    log.info("  Missing expression values: %d (%.3f%%)", na_count, na_pct)

    _save_expression(expr, out_dir)
    _save_metadata(tech_meta, out_dir, fname="technical_metadata.csv")

    report = {
        "dataset_id": gse_id,
        "role": dcfg["role"],
        "platform_id": dcfg["platform_id"],
        "tissue": dcfg["tissue"],
        "n_samples": expr.shape[1],
        "n_features": expr.shape[0],
        "feature_id_type": dcfg["feature_id_type"],
        "expression_scale": dcfg["expression_scale"],
        "expression_dtype": str(expr.dtypes.iloc[0]),
        "na_values": na_count,
        "na_pct": round(na_pct, 4),
        "ilmn_prefix_probes": ilmn_count,
        "probe_id_sample": probe_sample,
        "validation": {
            "duplicate_sample_ids": dup_check.passed,
            "expression_numeric": num_check.passed,
            "diagnosis_locked": True,
        },
        "notes": (
            "LOCKED external validation hold-out. "
            "Diagnosis labels NOT extracted or stored (Gate E rule). "
            "Only structural metadata saved to technical_metadata.csv. "
            "GPL15988 uses nuID-style probe IDs, not ILMN_* format. "
            "Probe crosswalk to GPL10558 required at Gate B before final evaluation."
        ),
    }
    _save_report(report, out_dir)
    return report


def load_gse125583(cfg: dict, root: Path) -> dict:
    """Load GSE125583 — brain RNA-seq biological support (GPL16791)."""
    gse_id = "GSE125583"
    log.info("=" * 60)
    log.info("Loading %s (brain RNA-seq support)", gse_id)
    log.info("=" * 60)

    dcfg = cfg["datasets"][gse_id]
    matrix_path = root / dcfg["series_matrix_path"]
    counts_path = root / dcfg["counts_path"]
    out_dir = root / "data" / "interim" / gse_id

    # --- Metadata from series matrix ---
    log.info("  Parsing metadata from series matrix...")
    header = load_series_matrix_header(matrix_path)
    meta = parse_gse125583_metadata(header)
    meta.insert(1, "dataset_id", gse_id)

    diag_counts = meta["diagnosis"].value_counts(dropna=False).to_dict()
    log.info("  Diagnosis counts: %s", diag_counts)

    # --- Raw integer counts ---
    log.info("  Loading RNA-seq count matrix from: %s", counts_path.name)
    counts = load_tsv_gz(counts_path, index_col=0)
    log.info("  Counts shape (genes × samples): %s", counts.shape)

    # Verify: index should be GeneID (Entrez integers)
    if counts.index.name:
        log.info("  Count matrix index name: %s", counts.index.name)

    # Sample count verification
    n_genes, n_samples_counts = counts.shape
    n_samples_meta = len(meta)
    log.info("  Metadata samples: %d, Count matrix samples: %d", n_samples_meta, n_samples_counts)

    # Verify values are non-negative integers
    sample_col = counts.iloc[:, 0]
    all_int = bool(np.issubdtype(counts.values.dtype, np.integer))
    min_val = float(counts.values.min())
    max_val = float(counts.values.max())
    any_negative = min_val < 0
    log.info("  Count dtype: %s, min: %.0f, max: %.0f, integer: %s",
             counts.dtypes.iloc[0], min_val, max_val, all_int)
    if any_negative:
        log.warning("  WARNING: Negative values detected in count matrix!")
    if not all_int:
        log.warning("  WARNING: Count matrix is not integer dtype — verify file format!")

    # NA check
    na_count = int(counts.isna().sum().sum())
    log.info("  Missing count values: %d", na_count)

    # Sample ID alignment check
    count_cols = list(counts.columns)
    meta_ids = meta["sample_id"].tolist()
    align_check = check_sample_id_alignment(meta_ids, count_cols)
    log.info("  %s", align_check)

    # Save counts as int32 to save memory (values are non-negative counts)
    counts_save = counts.astype("int32") if not any_negative else counts
    _save_expression(counts_save, out_dir, fname="raw_counts.parquet")
    _save_metadata(meta, out_dir)

    # Validation
    dup_check = check_duplicate_ids(meta_ids)
    _log_validation(gse_id, [dup_check, align_check])

    report = {
        "dataset_id": gse_id,
        "role": dcfg["role"],
        "platform_id": dcfg["platform_id"],
        "tissue": dcfg["tissue"],
        "n_samples_metadata": n_samples_meta,
        "n_samples_counts": n_samples_counts,
        "n_genes": n_genes,
        "count_index_name": str(counts.index.name),
        "count_dtype": str(counts.dtypes.iloc[0]),
        "count_all_integer": all_int,
        "count_min": min_val,
        "count_max": max_val,
        "na_values": na_count,
        "diagnosis_counts": diag_counts,
        "validation": {
            "duplicate_sample_ids": dup_check.passed,
            "sample_id_alignment": align_check.passed,
            "counts_non_negative": not any_negative,
        },
        "notes": (
            "Brain fusiform gyrus RNA-seq (post-mortem). "
            "Raw integer counts: use PyDESeq2, NOT limma or microarray tools. "
            "Gene IDs are NCBI Entrez integers; GENCODE GTF uses Ensembl IDs. "
            "Crosswalk via Human.GRCh38.p13.annot.tsv.gz needed before merging. "
            "CANNOT serve as independent blood validation — different tissue."
        ),
    }
    _save_report(report, out_dir)
    return report


# ---------------------------------------------------------------------------
# Validation logging helper
# ---------------------------------------------------------------------------

def _log_validation(gse_id: str, results: list[ValidationResult]) -> None:
    for r in results:
        level = logging.INFO if r.passed else logging.WARNING
        log.log(level, "  [%s] %s — %s: %s",
                gse_id, "PASS" if r.passed else "FAIL", r.check_name, r.message)


# ---------------------------------------------------------------------------
# Combined summary table
# ---------------------------------------------------------------------------

def save_dataset_summary(reports: list[dict], root: Path) -> None:
    """Save outputs/tables/dataset_summary.csv from parsing reports."""
    rows = []
    for r in reports:
        row = {
            "dataset_id": r.get("dataset_id"),
            "role": r.get("role"),
            "platform_id": r.get("platform_id"),
            "tissue": r.get("tissue"),
            "n_samples": r.get("n_samples") or r.get("n_samples_counts"),
            "n_features": r.get("n_features") or r.get("n_genes"),
            "expression_scale": r.get("expression_scale", "raw_integer_counts"),
            "feature_id_type": r.get("feature_id_type", "entrez_gene_id"),
            "na_values": r.get("na_values", 0),
            "na_pct": r.get("na_pct", 0.0),
            "diagnosis_locked": r.get("validation", {}).get("diagnosis_locked", False),
            "notes": r.get("notes", ""),
        }
        rows.append(row)

    summary = pd.DataFrame(rows)
    out_path = root / "outputs" / "tables" / "dataset_summary.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_path, index=False)
    log.info("Saved dataset summary: %s", out_path)
    print("\n--- DATASET SUMMARY ---")
    print(summary.to_string(index=False))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    root = get_project_root()
    cfg = load_dataset_config()

    log.info("Project root: %s", root)
    log.info("Loaded config with datasets: %s", list(cfg["datasets"].keys()))

    reports: list[dict] = []

    # GSE63061 — primary discovery
    try:
        r = load_gse63061(cfg, root)
        reports.append(r)
    except Exception as exc:
        log.error("GSE63061 FAILED: %s", exc, exc_info=True)

    # GSE63060 — replication
    try:
        r = load_gse63060(cfg, root)
        reports.append(r)
    except Exception as exc:
        log.error("GSE63060 FAILED: %s", exc, exc_info=True)

    # GSE140829 — LOCKED external validation
    try:
        r = load_gse140829(cfg, root)
        reports.append(r)
    except Exception as exc:
        log.error("GSE140829 FAILED: %s", exc, exc_info=True)

    # GSE125583 — brain RNA-seq
    try:
        r = load_gse125583(cfg, root)
        reports.append(r)
    except Exception as exc:
        log.error("GSE125583 FAILED: %s", exc, exc_info=True)

    # Combined summary
    if reports:
        save_dataset_summary(reports, root)

    log.info("Stage 08 loading complete. %d/%d datasets loaded successfully.",
             len(reports), 4)


if __name__ == "__main__":
    main()
