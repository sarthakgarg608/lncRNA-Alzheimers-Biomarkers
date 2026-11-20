"""
scripts/02_extract_metadata.py
================================
Stage 08 — Standardise metadata and produce audit tables + plots.

This script:
  1. Reads the interim metadata files produced by 01_load_geo_datasets.py.
  2. Standardises them to a canonical schema.
  3. Runs validation checks (duplicates, label vocabulary, missing values).
  4. Writes:
       outputs/tables/sample_metadata_summary.csv
       outputs/tables/missing_values_summary.csv
  5. Produces plots in outputs/figures/:
       sample_counts_by_diagnosis.png
       missing_metadata_pct.png
       expression_distributions.png   (one panel per dataset)
       dataset_dimension_comparison.png

Scientific integrity rules
--------------------------
- GSE140829 diagnosis is NEVER displayed or summarised (locked).
- MCI is always a separate category from Control.
- Blood and brain datasets are never merged.
- No normalisation or ML happens here.
- Original sample IDs are preserved.

Usage
-----
  .venv\\Scripts\\python.exe scripts/02_extract_metadata.py
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from scripts.utils.paths import get_project_root, load_dataset_config
from scripts.utils.validation import (
    check_duplicate_ids,
    check_valid_diagnosis_labels,
    ValidationResult,
    VALID_DIAGNOSIS_LABELS,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Canonical metadata schema
# Columns that every unlocked dataset's metadata.csv must contain.
# ---------------------------------------------------------------------------
CANONICAL_COLUMNS = [
    "sample_id",
    "dataset_id",
    "diagnosis",
    "age",
    "sex",
    "tissue",
    "platform_id",
    "batch",
]

# GSE140829 (locked) only gets these columns:
LOCKED_COLUMNS = [
    "sample_id",
    "dataset_id",
    "age",
    "sex",
    "platform_id",
    "batch",
]


# ---------------------------------------------------------------------------
# Metadata standardisation
# ---------------------------------------------------------------------------

def standardise_metadata(df: pd.DataFrame, dataset_id: str, locked: bool = False) -> pd.DataFrame:
    """
    Ensure canonical column schema is present.
    Missing optional columns are filled with None.
    Diagnosis labels are validated against the allowed vocabulary.

    Parameters
    ----------
    df : pd.DataFrame
        Raw metadata loaded from data/interim/<GSEID>/metadata.csv or
        technical_metadata.csv.
    dataset_id : str
        GEO accession (e.g. "GSE63061").
    locked : bool
        True for GSE140829 — no diagnosis column.

    Returns
    -------
    pd.DataFrame
        Standardised DataFrame.
    """
    expected = LOCKED_COLUMNS if locked else CANONICAL_COLUMNS

    # Ensure dataset_id column exists
    if "dataset_id" not in df.columns:
        df = df.copy()
        df.insert(1, "dataset_id", dataset_id)

    # Add missing columns
    for col in expected:
        if col not in df.columns:
            df[col] = None

    # Subset to canonical columns only (preserve any extras after)
    extra_cols = [c for c in df.columns if c not in expected]
    df = df[expected + extra_cols].copy()

    if not locked and "diagnosis" in df.columns:
        # Validate labels
        result = check_valid_diagnosis_labels(df["diagnosis"].tolist())
        if not result.passed:
            log.warning("[%s] Diagnosis label issue: %s", dataset_id, result.message)

    return df


# ---------------------------------------------------------------------------
# Diagnosis counts helper
# ---------------------------------------------------------------------------

def _diagnosis_counts(meta: pd.DataFrame, dataset_id: str, locked: bool) -> dict[str, Any]:
    if locked or "diagnosis" not in meta.columns:
        return {
            "dataset_id": dataset_id,
            "AD": "LOCKED",
            "Control": "LOCKED",
            "MCI": "LOCKED",
            "Unknown": "LOCKED",
            "total_samples": len(meta),
            "diagnosis_locked": True,
        }
    counts = meta["diagnosis"].value_counts(dropna=False).to_dict()
    return {
        "dataset_id": dataset_id,
        "AD": int(counts.get("AD", 0)),
        "Control": int(counts.get("Control", 0)),
        "MCI": int(counts.get("MCI", 0)),
        "Unknown": int(counts.get("Unknown", 0)),
        "total_samples": len(meta),
        "diagnosis_locked": False,
    }


# ---------------------------------------------------------------------------
# Missing values summary
# ---------------------------------------------------------------------------

def _missing_summary(meta: pd.DataFrame, dataset_id: str) -> list[dict]:
    rows = []
    for col in meta.columns:
        n_missing = int(meta[col].isna().sum())
        pct = 100 * n_missing / len(meta) if len(meta) > 0 else 0.0
        rows.append({
            "dataset_id": dataset_id,
            "column": col,
            "n_missing": n_missing,
            "pct_missing": round(pct, 2),
        })
    return rows


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _plot_sample_counts(summary_rows: list[dict], out_path: Path) -> None:
    """Bar chart: sample counts by diagnosis for unlocked datasets."""
    unlocked = [r for r in summary_rows if not r.get("diagnosis_locked", False)]
    if not unlocked:
        log.warning("No unlocked datasets to plot sample counts.")
        return

    datasets = [r["dataset_id"] for r in unlocked]
    ad_counts = [r["AD"] for r in unlocked]
    ctrl_counts = [r["Control"] for r in unlocked]
    mci_counts = [r["MCI"] for r in unlocked]
    unk_counts = [r["Unknown"] for r in unlocked]

    x = np.arange(len(datasets))
    width = 0.2

    fig, ax = plt.subplots(figsize=(9, 5))
    bars_ad = ax.bar(x - 1.5 * width, ad_counts, width, label="AD", color="#E74C3C")
    bars_ctrl = ax.bar(x - 0.5 * width, ctrl_counts, width, label="Control", color="#2ECC71")
    bars_mci = ax.bar(x + 0.5 * width, mci_counts, width, label="MCI", color="#F39C12")
    bars_unk = ax.bar(x + 1.5 * width, unk_counts, width, label="Unknown", color="#95A5A6")

    ax.set_xlabel("Dataset", fontsize=12)
    ax.set_ylabel("Sample count", fontsize=12)
    ax.set_title("Sample counts by diagnosis\n(GSE140829 locked — not shown)", fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.3)

    for bars in [bars_ad, bars_ctrl, bars_mci, bars_unk]:
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.annotate(str(int(h)), xy=(bar.get_x() + bar.get_width() / 2, h),
                            xytext=(0, 3), textcoords="offset points",
                            ha="center", va="bottom", fontsize=8)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Saved: %s", out_path.name)


def _plot_missing_pct(missing_rows: list[dict], out_path: Path) -> None:
    """Bar chart: percentage of missing metadata per column per dataset."""
    df = pd.DataFrame(missing_rows)
    if df.empty:
        return

    # Focus on key metadata columns
    key_cols = ["diagnosis", "age", "sex", "tissue", "batch", "platform_id"]
    df = df[df["column"].isin(key_cols)]

    if df.empty:
        return

    datasets = df["dataset_id"].unique()
    n_ds = len(datasets)
    fig, axes = plt.subplots(1, n_ds, figsize=(4 * n_ds, 5), sharey=True)
    if n_ds == 1:
        axes = [axes]

    for ax, ds in zip(axes, datasets):
        sub = df[df["dataset_id"] == ds].set_index("column")
        cols_present = [c for c in key_cols if c in sub.index]
        vals = [sub.loc[c, "pct_missing"] for c in cols_present]
        colors = ["#E74C3C" if v > 10 else "#2ECC71" if v == 0 else "#F39C12" for v in vals]
        ax.barh(cols_present, vals, color=colors)
        ax.set_title(ds, fontsize=11)
        ax.set_xlabel("% missing", fontsize=9)
        ax.axvline(10, color="red", linestyle="--", alpha=0.5, linewidth=1)
        ax.set_xlim(0, 105)

    fig.suptitle("Missing metadata by column and dataset", fontsize=13, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: %s", out_path.name)


def _plot_expression_distributions(root: Path, datasets_meta: list[tuple[str, str]]) -> None:
    """
    Violin/box plots of expression value distributions per dataset.
    Loads a random sample of rows from each parquet file to avoid memory issues.

    Parameters
    ----------
    root : Path
    datasets_meta : list of (dataset_id, parquet_filename)
        e.g. [("GSE63061", "expression.parquet"), ("GSE125583", "raw_counts.parquet")]
    """
    fig, axes = plt.subplots(1, len(datasets_meta), figsize=(5 * len(datasets_meta), 5))
    if len(datasets_meta) == 1:
        axes = [axes]

    out_path = root / "outputs" / "figures" / "expression_distributions.png"

    for ax, (ds_id, parquet_fname) in zip(axes, datasets_meta):
        parquet_path = root / "data" / "interim" / ds_id / parquet_fname
        if not parquet_path.exists():
            ax.set_title(f"{ds_id}\n(file missing)")
            continue

        try:
            # Load only the first 500 rows (features) × all samples for speed
            table = pq.read_table(parquet_path)
            df_sample = table.to_pandas().iloc[:500, :]
        except Exception as exc:
            log.warning("Could not load %s for plotting: %s", parquet_path.name, exc)
            ax.set_title(f"{ds_id}\n(load error)")
            continue

        # Flatten values
        vals = df_sample.values.flatten()
        vals = vals[~np.isnan(vals)]
        if len(vals) > 20_000:
            rng = np.random.default_rng(42)
            vals = rng.choice(vals, 20_000, replace=False)

        is_counts = "count" in parquet_fname
        if is_counts:
            # Log1p transform for counts
            vals = np.log1p(vals.astype(float))
            xlabel = "log1p(count)"
        else:
            xlabel = "Expression intensity"

        ax.violinplot(vals, positions=[0], showmedians=True)
        ax.set_title(f"{ds_id}", fontsize=11)
        ax.set_ylabel(xlabel, fontsize=9)
        ax.set_xticks([])

    fig.suptitle("Expression value distributions (sample of ≤500 features)", fontsize=12)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Saved: %s", out_path.name)


def _plot_dataset_dimensions(summary_rows: list[dict], out_path: Path) -> None:
    """Bar chart comparing dataset dimensions (samples and features)."""
    datasets = [r["dataset_id"] for r in summary_rows]
    n_samples = [r.get("total_samples", 0) for r in summary_rows]
    # n_features will come from parsing reports
    fig, ax = plt.subplots(figsize=(8, 4))
    bars = ax.bar(datasets, n_samples, color=["#3498DB", "#1ABC9C", "#9B59B6", "#E67E22"])
    ax.set_ylabel("Number of samples", fontsize=11)
    ax.set_title("Sample counts per dataset\n(GSE140829 count is structural, diagnosis locked)",
                 fontsize=12)
    for bar in bars:
        h = bar.get_height()
        ax.annotate(str(int(h)), xy=(bar.get_x() + bar.get_width() / 2, h),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Saved: %s", out_path.name)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    root = get_project_root()
    cfg = load_dataset_config()

    datasets_cfg = {
        "GSE63061":  {"locked": False, "meta_file": "metadata.csv"},
        "GSE63060":  {"locked": False, "meta_file": "metadata.csv"},
        "GSE140829": {"locked": True,  "meta_file": "technical_metadata.csv"},
        "GSE125583": {"locked": False, "meta_file": "metadata.csv"},
    }

    summary_rows: list[dict] = []
    missing_rows: list[dict] = []
    validation_issues: list[str] = []
    all_meta: dict[str, pd.DataFrame] = {}

    for ds_id, ds_opts in datasets_cfg.items():
        interim_dir = root / "data" / "interim" / ds_id
        meta_path = interim_dir / ds_opts["meta_file"]
        locked = ds_opts["locked"]

        if not meta_path.exists():
            log.warning("[%s] metadata file not found: %s — skipping", ds_id, meta_path)
            validation_issues.append(f"{ds_id}: metadata file missing — run 01_load_geo_datasets.py first")
            continue

        log.info("Processing %s metadata...", ds_id)
        raw_meta = pd.read_csv(meta_path)

        # Standardise
        meta = standardise_metadata(raw_meta, ds_id, locked=locked)
        all_meta[ds_id] = meta

        # Duplicate check
        dup_check = check_duplicate_ids(meta["sample_id"].tolist())
        if not dup_check.passed:
            validation_issues.append(f"{ds_id}: {dup_check.message}")
            log.warning("[%s] %s", ds_id, dup_check)
        else:
            log.info("[%s] No duplicate sample IDs — PASS", ds_id)

        # Diagnosis label check (unlocked only)
        if not locked and "diagnosis" in meta.columns:
            diag_check = check_valid_diagnosis_labels(meta["diagnosis"].dropna().tolist())
            if not diag_check.passed:
                validation_issues.append(f"{ds_id}: {diag_check.message}")
                log.warning("[%s] %s", ds_id, diag_check)
            else:
                log.info("[%s] Diagnosis labels valid — PASS", ds_id)

        # Counts
        sc = _diagnosis_counts(meta, ds_id, locked)
        summary_rows.append(sc)
        log.info("[%s] Sample counts: %s", ds_id,
                 {k: v for k, v in sc.items() if k not in ("dataset_id", "diagnosis_locked")})

        # Missing values
        missing_rows.extend(_missing_summary(meta, ds_id))

    # --- Output tables ---
    out_tables = root / "outputs" / "tables"
    out_tables.mkdir(parents=True, exist_ok=True)

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(out_tables / "sample_metadata_summary.csv", index=False)
        log.info("Saved: sample_metadata_summary.csv")
        print("\n--- SAMPLE METADATA SUMMARY ---")
        print(summary_df.to_string(index=False))

    if missing_rows:
        missing_df = pd.DataFrame(missing_rows)
        missing_df.to_csv(out_tables / "missing_values_summary.csv", index=False)
        log.info("Saved: missing_values_summary.csv")

    # --- Plots ---
    out_figs = root / "outputs" / "figures"
    out_figs.mkdir(parents=True, exist_ok=True)

    _plot_sample_counts(summary_rows, out_figs / "sample_counts_by_diagnosis.png")
    _plot_missing_pct(missing_rows, out_figs / "missing_metadata_pct.png")

    # Expression distribution plots — one panel per dataset
    parquet_map = [
        ("GSE63061",  "expression.parquet"),
        ("GSE63060",  "expression.parquet"),
        ("GSE140829", "expression.parquet"),
        ("GSE125583", "raw_counts.parquet"),
    ]
    existing_parquet = [
        (ds, pf) for ds, pf in parquet_map
        if (root / "data" / "interim" / ds / pf).exists()
    ]
    if existing_parquet:
        _plot_expression_distributions(root, existing_parquet)

    _plot_dataset_dimensions(summary_rows, out_figs / "dataset_dimension_comparison.png")

    # --- Validation summary ---
    if validation_issues:
        log.warning("Validation issues found:")
        for issue in validation_issues:
            log.warning("  - %s", issue)
    else:
        log.info("All validation checks passed.")

    log.info("Stage 08 metadata extraction complete.")


if __name__ == "__main__":
    main()
