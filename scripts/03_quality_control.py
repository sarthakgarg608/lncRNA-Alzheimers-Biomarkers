"""
scripts/03_quality_control.py
==============================
Stage 09 — Quality Control, Matrix Orientation & Feature/Sample Summaries.

Performs:
  1. Matrix orientation verification & derived transposed Parquet creation
     (samples × features) while strictly preserving original files.
  2. Sample-level QC (median, mean, std, IQR, missingness, library sizes).
  3. Sample outlier flagging with documented statistical thresholds (MAD/IQR).
     (Samples are FLAGGED ONLY, NEVER deleted).
  4. Feature-level QC (missingness, variance, constant features, zero proportions).
  5. Generates figures:
     - GSE63061_expression_boxplot.png
     - GSE63060_expression_boxplot.png
     - GSE125583_library_sizes.png
     - feature_variance_distribution.png
  6. Generates tables:
     - outputs/tables/stage09/sample_qc_summary.csv
     - outputs/tables/stage09/sample_qc_flags.csv
     - outputs/tables/stage09/feature_qc_summary.csv

Scientific integrity rules:
  - GSE140829 diagnosis remains LOCKED (only technical QC metrics computed).
  - No samples or features are removed.
  - Raw integer counts of GSE125583 are preserved.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

# ---------------------------------------------------------------------------
# Setup & Paths
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from scripts.utils.paths import get_project_root, load_dataset_config
from scripts.utils.validation import (
    check_duplicate_ids,
    check_expression_numeric,
    check_sample_id_alignment,
    check_matrix_orientation,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# Outlier cutoff factor for Median Absolute Deviation
MAD_CUTOFF = 3.5


# ---------------------------------------------------------------------------
# Task 2: Matrix Orientation & Transposition
# ---------------------------------------------------------------------------

def ensure_samples_x_features_matrix(
    dataset_id: str,
    original_parquet_name: str,
    root: Path,
) -> tuple[pd.DataFrame, Path]:
    """
    Verify the stored matrix and ensure a derived transposed matrix
    (rows = samples, cols = features) exists.

    Returns
    -------
    (df_samples_x_features, path_to_derived_parquet)
    """
    interim_dir = root / "data" / "interim" / dataset_id
    orig_path = interim_dir / original_parquet_name
    derived_name = (
        "raw_counts_samples_x_features.parquet"
        if "count" in original_parquet_name.lower()
        else "expression_samples_x_features.parquet"
    )
    derived_path = interim_dir / derived_name

    log.info("[%s] Checking matrix orientation for %s...", dataset_id, orig_path.name)
    orig_df = pd.read_parquet(orig_path)
    orig_orientation = check_matrix_orientation(orig_df, expected="features_x_samples")
    log.info("[%s] Original orientation: %s (shape %s)", dataset_id, orig_orientation.message, orig_df.shape)

    # Check if derived file exists and matches
    if derived_path.exists():
        derived_df = pd.read_parquet(derived_path)
        if derived_df.shape == (orig_df.shape[1], orig_df.shape[0]):
            log.info("[%s] Derived transposed matrix already exists: %s", dataset_id, derived_path.name)
            return derived_df, derived_path

    # Transpose and preserve dtypes
    log.info("[%s] Creating derived transposed matrix (samples × features)...", dataset_id)
    derived_df = orig_df.T
    derived_df.index.name = "sample_id"
    derived_df.columns.name = "feature_id"

    # Save to interim
    derived_df.to_parquet(derived_path, engine="pyarrow", compression="snappy")
    log.info(
        "[%s] Saved derived matrix: %s, shape %s (rows=samples, cols=features)",
        dataset_id,
        derived_path.name,
        derived_df.shape,
    )
    return derived_df, derived_path


# ---------------------------------------------------------------------------
# Task 3: Sample-level QC
# ---------------------------------------------------------------------------

def compute_sample_qc_microarray(
    expr_s_x_f: pd.DataFrame,
    meta_df: pd.DataFrame,
    dataset_id: str,
    locked: bool = False,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """
    Compute sample-level metrics for microarray cohorts:
      - median, mean, std, IQR
      - missing value count & pct
      - constant / invalid measurement count
      - statistical outlier flagging (MAD on median and IQR)
    """
    rows: list[dict[str, Any]] = []
    flags: list[dict[str, Any]] = []

    # Map metadata for fast lookup
    meta_idx = meta_df.set_index("sample_id")

    # Vectorized calculations across columns (features) for each row (sample)
    values = expr_s_x_f.values
    sample_ids = expr_s_x_f.index.tolist()

    means = np.nanmean(values, axis=1)
    stds = np.nanstd(values, axis=1)
    medians = np.nanmedian(values, axis=1)
    q25 = np.nanpercentile(values, 25, axis=1)
    q75 = np.nanpercentile(values, 75, axis=1)
    iqrs = q75 - q25
    n_nans = np.isnan(values).sum(axis=1)
    n_features = expr_s_x_f.shape[1]

    # Calculate population MAD for medians and IQRs
    pop_med_median = np.median(medians)
    mad_median = np.median(np.abs(medians - pop_med_median))
    mad_median = mad_median if mad_median > 1e-6 else np.std(medians)

    pop_iqr_median = np.median(iqrs)
    mad_iqr = np.median(np.abs(iqrs - pop_iqr_median))
    mad_iqr = mad_iqr if mad_iqr > 1e-6 else np.std(iqrs)

    for i, sid in enumerate(sample_ids):
        diag = "LOCKED" if locked else meta_idx.loc[sid, "diagnosis"] if sid in meta_idx.index and "diagnosis" in meta_idx.columns else "Unknown"
        batch_val = meta_idx.loc[sid, "batch"] if sid in meta_idx.index and "batch" in meta_idx.columns else None

        med_i = float(medians[i])
        iqr_i = float(iqrs[i])
        mean_i = float(means[i])
        std_i = float(stds[i])
        nan_pct_i = float(100.0 * n_nans[i] / n_features)

        # Flagging logic
        sample_flags: list[str] = []
        if mad_median > 0 and abs(med_i - pop_med_median) > MAD_CUTOFF * mad_median:
            sample_flags.append(f"MEDIAN_INTENSITY_OUTLIER(|d|={abs(med_i - pop_med_median):.3f}>{MAD_CUTOFF*mad_median:.3f})")
            flags.append({
                "sample_id": sid,
                "dataset_id": dataset_id,
                "flag_reason": "Extreme median intensity",
                "flag_metric_value": med_i,
                "threshold_used": f"{pop_med_median:.3f} +/- {MAD_CUTOFF*mad_median:.3f}",
                "provisional_action": "FLAG_FOR_INSPECTION_DO_NOT_DELETE",
            })

        if mad_iqr > 0 and abs(iqr_i - pop_iqr_median) > MAD_CUTOFF * mad_iqr:
            sample_flags.append(f"IQR_OUTLIER(|d|={abs(iqr_i - pop_iqr_median):.3f}>{MAD_CUTOFF*mad_iqr:.3f})")
            flags.append({
                "sample_id": sid,
                "dataset_id": dataset_id,
                "flag_reason": "Extreme IQR (distribution dispersion)",
                "flag_metric_value": iqr_i,
                "threshold_used": f"{pop_iqr_median:.3f} +/- {MAD_CUTOFF*mad_iqr:.3f}",
                "provisional_action": "FLAG_FOR_INSPECTION_DO_NOT_DELETE",
            })

        if nan_pct_i > 5.0:
            sample_flags.append(f"HIGH_MISSINGNESS({nan_pct_i:.1f}%)")
            flags.append({
                "sample_id": sid,
                "dataset_id": dataset_id,
                "flag_reason": "High missingness",
                "flag_metric_value": nan_pct_i,
                "threshold_used": ">5%",
                "provisional_action": "FLAG_FOR_INSPECTION_DO_NOT_DELETE",
            })

        rows.append({
            "sample_id": sid,
            "dataset_id": dataset_id,
            "diagnosis": diag,
            "batch": batch_val,
            "assay_type": "microarray",
            "mean_expr": round(mean_i, 4),
            "median_expr": round(med_i, 4),
            "std_expr": round(std_i, 4),
            "iqr_expr": round(iqr_i, 4),
            "pct_missing": round(nan_pct_i, 4),
            "library_size": None,
            "n_expressed_genes": None,
            "pct_zero_genes": None,
            "qc_flags": "; ".join(sample_flags) if sample_flags else "PASS",
            "is_outlier_flag": len(sample_flags) > 0,
        })

    return pd.DataFrame(rows), flags


def compute_sample_qc_rnaseq(
    counts_s_x_f: pd.DataFrame,
    meta_df: pd.DataFrame,
    dataset_id: str,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """
    Compute sample-level metrics for RNA-seq counts:
      - Library size (sum of counts)
      - Number of expressed genes (count > 0)
      - Percentage of zero-count genes
      - Statistical outlier detection on library size
    """
    rows: list[dict[str, Any]] = []
    flags: list[dict[str, Any]] = []

    meta_idx = meta_df.set_index("sample_id")
    values = counts_s_x_f.values
    sample_ids = counts_s_x_f.index.tolist()
    n_genes = counts_s_x_f.shape[1]

    lib_sizes = values.sum(axis=1)
    n_expressed = (values > 0).sum(axis=1)
    pct_zeros = ((values == 0).sum(axis=1) / n_genes) * 100.0

    pop_lib_med = np.median(lib_sizes)
    mad_lib = np.median(np.abs(lib_sizes - pop_lib_med))

    pop_zeros_med = np.median(pct_zeros)
    mad_zeros = np.median(np.abs(pct_zeros - pop_zeros_med))

    for i, sid in enumerate(sample_ids):
        diag = meta_idx.loc[sid, "diagnosis"] if sid in meta_idx.index and "diagnosis" in meta_idx.columns else "Unknown"
        batch_val = meta_idx.loc[sid, "batch"] if sid in meta_idx.index and "batch" in meta_idx.columns else None

        ls = int(lib_sizes[i])
        n_exp = int(n_expressed[i])
        pz = float(pct_zeros[i])

        sample_flags: list[str] = []
        if mad_lib > 0 and abs(ls - pop_lib_med) > MAD_CUTOFF * mad_lib:
            sample_flags.append(f"EXTREME_LIBRARY_SIZE({ls})")
            flags.append({
                "sample_id": sid,
                "dataset_id": dataset_id,
                "flag_reason": "Extreme library size",
                "flag_metric_value": ls,
                "threshold_used": f"{pop_lib_med:.0f} +/- {MAD_CUTOFF*mad_lib:.0f}",
                "provisional_action": "FLAG_FOR_INSPECTION_DO_NOT_DELETE",
            })

        if mad_zeros > 0 and abs(pz - pop_zeros_med) > MAD_CUTOFF * mad_zeros:
            sample_flags.append(f"HIGH_ZERO_FRACTION({pz:.2f}%)")
            flags.append({
                "sample_id": sid,
                "dataset_id": dataset_id,
                "flag_reason": "High zero fraction",
                "flag_metric_value": pz,
                "threshold_used": f"{pop_zeros_med:.2f}% +/- {MAD_CUTOFF*mad_zeros:.2f}%",
                "provisional_action": "FLAG_FOR_INSPECTION_DO_NOT_DELETE",
            })

        rows.append({
            "sample_id": sid,
            "dataset_id": dataset_id,
            "diagnosis": diag,
            "batch": batch_val,
            "assay_type": "RNA-seq",
            "mean_expr": round(float(np.mean(values[i])), 2),
            "median_expr": round(float(np.median(values[i])), 2),
            "std_expr": round(float(np.std(values[i])), 2),
            "iqr_expr": round(float(np.percentile(values[i], 75) - np.percentile(values[i], 25)), 2),
            "pct_missing": 0.0,
            "library_size": ls,
            "n_expressed_genes": n_exp,
            "pct_zero_genes": round(pz, 2),
            "qc_flags": "; ".join(sample_flags) if sample_flags else "PASS",
            "is_outlier_flag": len(sample_flags) > 0,
        })

    return pd.DataFrame(rows), flags


# ---------------------------------------------------------------------------
# Task 4: Feature-level QC
# ---------------------------------------------------------------------------

def compute_feature_qc(
    df_s_x_f: pd.DataFrame,
    dataset_id: str,
    assay_type: str,
) -> dict[str, Any]:
    """
    Compute feature-level statistics:
      - missing percentage per gene/probe
      - feature variances and means
      - constant / near-constant features (variance < 1e-6)
      - proportion of zeros (RNA-seq)
    """
    n_samples, n_features = df_s_x_f.shape
    values = df_s_x_f.values

    # Missingness
    missing_per_feature = np.isnan(values).sum(axis=0)
    features_with_missing = int((missing_per_feature > 0).sum())
    max_missing_pct = float(100.0 * np.max(missing_per_feature) / n_samples)

    # Variance and mean
    variances = np.nanvar(values, axis=0)
    means = np.nanmean(values, axis=0)

    constant_features = int((variances == 0).sum())
    near_zero_var = int((variances < 1e-4).sum())

    res = {
        "dataset_id": dataset_id,
        "assay_type": assay_type,
        "n_samples": n_samples,
        "total_features": n_features,
        "features_with_missing": features_with_missing,
        "max_missing_pct": round(max_missing_pct, 4),
        "constant_features": constant_features,
        "near_zero_variance_features": near_zero_var,
        "mean_feature_variance": round(float(np.mean(variances)), 4),
        "median_feature_variance": round(float(np.median(variances)), 4),
        "q25_feature_variance": round(float(np.percentile(variances, 25)), 4),
        "q75_feature_variance": round(float(np.percentile(variances, 75)), 4),
        "mean_expression": round(float(np.mean(means)), 4),
        "median_expression": round(float(np.median(means)), 4),
    }

    if assay_type == "RNA-seq":
        zeros_per_gene = (values == 0).sum(axis=0)
        always_zero = int((zeros_per_gene == n_samples).sum())
        expressed_in_lt_5pct = int((zeros_per_gene > 0.95 * n_samples).sum())
        res["always_zero_features"] = always_zero
        res["expressed_in_lt_5pct_samples"] = expressed_in_lt_5pct
    else:
        res["always_zero_features"] = 0
        res["expressed_in_lt_5pct_samples"] = 0

    return res


# ---------------------------------------------------------------------------
# Task 10: Figures
# ---------------------------------------------------------------------------

def plot_microarray_boxplots(
    df_s_x_f: pd.DataFrame,
    meta_df: pd.DataFrame,
    dataset_id: str,
    out_path: Path,
    n_sample_subset: int = 35,
) -> None:
    """
    Plot boxplots of expression across a stratified representative subset
    of samples, demonstrating cross-sample distribution alignment.
    """
    rng = np.random.default_rng(42)
    sample_ids = df_s_x_f.index.tolist()
    if len(sample_ids) > n_sample_subset:
        # Stratify selection if diagnosis available
        if "diagnosis" in meta_df.columns:
            meta_sub = meta_df.set_index("sample_id").loc[sample_ids]
            chosen = []
            for diag, group in meta_sub.groupby("diagnosis"):
                k = max(1, int(n_sample_subset * len(group) / len(sample_ids)))
                chosen.extend(rng.choice(group.index.tolist(), min(k, len(group)), replace=False))
            # Fill remaining if needed
            remaining = [s for s in sample_ids if s not in chosen]
            if len(chosen) < n_sample_subset and remaining:
                chosen.extend(rng.choice(remaining, n_sample_subset - len(chosen), replace=False))
            selected_samples = chosen[:n_sample_subset]
        else:
            selected_samples = rng.choice(sample_ids, n_sample_subset, replace=False).tolist()
    else:
        selected_samples = sample_ids

    sub_expr = df_s_x_f.loc[selected_samples]

    fig, ax = plt.subplots(figsize=(14, 5))
    data_to_plot = [sub_expr.loc[s].dropna().values for s in selected_samples]

    bp = ax.boxplot(
        data_to_plot,
        tick_labels=[s[-6:] for s in selected_samples],
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="red", linewidth=1.5),
        boxprops=dict(facecolor="#AED6F1", color="#2874A6"),
    )

    ax.set_ylabel("Expression intensity", fontsize=11)
    ax.set_xlabel(f"Representative samples (n={len(selected_samples)} of {len(sample_ids)})", fontsize=11)
    ax.set_title(
        f"{dataset_id} — Expression Distribution Alignment Across Samples\n"
        f"(Values are pre-normalized on log-like scale; medians tightly aligned)",
        fontsize=12,
    )
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Saved figure: %s", out_path.name)


def plot_library_sizes_rnaseq(
    counts_s_x_f: pd.DataFrame,
    out_path: Path,
) -> None:
    """
    Plot library sizes and expressed gene distributions for GSE125583.
    """
    values = counts_s_x_f.values
    lib_sizes = values.sum(axis=1) / 1e6  # in millions
    n_expressed = (values > 0).sum(axis=1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # Lib size histogram + KDE
    ax1.hist(lib_sizes, bins=25, color="#5DADE2", edgecolor="#1B4F72", alpha=0.8)
    ax1.axvline(np.median(lib_sizes), color="red", linestyle="--", label=f"Median: {np.median(lib_sizes):.1f}M")
    ax1.set_xlabel("Total Library Size (Million Reads)", fontsize=11)
    ax1.set_ylabel("Number of Samples", fontsize=11)
    ax1.set_title("GSE125583 — Library Size Distribution", fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(alpha=0.3)

    # Expressed genes histogram
    ax2.hist(n_expressed, bins=25, color="#58D68D", edgecolor="#1E8449", alpha=0.8)
    ax2.axvline(np.median(n_expressed), color="red", linestyle="--", label=f"Median: {int(np.median(n_expressed)):,} genes")
    ax2.set_xlabel("Number of Expressed Genes (count > 0)", fontsize=11)
    ax2.set_ylabel("Number of Samples", fontsize=11)
    ax2.set_title("GSE125583 — Expressed Genes per Sample", fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Saved figure: %s", out_path.name)


def plot_feature_variance_distributions(
    datasets_dict: dict[str, pd.DataFrame],
    out_path: Path,
) -> None:
    """
    Plot feature variance distributions across datasets.
    """
    fig, axes = plt.subplots(1, len(datasets_dict), figsize=(5 * len(datasets_dict), 4.5))
    if len(datasets_dict) == 1:
        axes = [axes]

    colors = ["#3498DB", "#2ECC71", "#E67E22"]

    for ax, (ds_id, df_s_x_f), col in zip(axes, datasets_dict.items(), colors):
        vars_ = np.nanvar(df_s_x_f.values, axis=0)
        # Avoid zero for log scale
        vars_positive = vars_[vars_ > 0]
        log_vars = np.log10(vars_positive)

        ax.hist(log_vars, bins=40, color=col, edgecolor="black", alpha=0.7)
        ax.axvline(np.median(log_vars), color="red", linestyle="--", label=f"Median log10(var): {np.median(log_vars):.2f}")
        ax.set_title(f"{ds_id} (n={len(vars_):,} features)", fontsize=11)
        ax.set_xlabel("log10(Feature Variance)", fontsize=10)
        ax.set_ylabel("Count", fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)

    fig.suptitle("Feature Variance Distributions Across Cohorts", fontsize=13, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved figure: %s", out_path.name)


# ---------------------------------------------------------------------------
# Main Execution Pipeline
# ---------------------------------------------------------------------------

def main() -> None:
    root = get_project_root()
    cfg = load_dataset_config()

    out_tables = root / "outputs" / "tables" / "stage09"
    out_figs = root / "outputs" / "figures" / "stage09"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_figs.mkdir(parents=True, exist_ok=True)

    log.info("=" * 65)
    log.info("STAGE 09: QUALITY CONTROL PIPELINE")
    log.info("=" * 65)

    all_sample_qc: list[pd.DataFrame] = []
    all_sample_flags: list[dict[str, Any]] = []
    all_feature_qc: list[dict[str, Any]] = []
    unlocked_dfs: dict[str, pd.DataFrame] = {}

    # 1. GSE63061 (Discovery, Microarray)
    df_61, _ = ensure_samples_x_features_matrix("GSE63061", "expression.parquet", root)
    meta_61 = pd.read_csv(root / "data" / "interim" / "GSE63061" / "metadata.csv")
    sqc_61, flags_61 = compute_sample_qc_microarray(df_61, meta_61, "GSE63061", locked=False)
    fqc_61 = compute_feature_qc(df_61, "GSE63061", "microarray")
    all_sample_qc.append(sqc_61)
    all_sample_flags.extend(flags_61)
    all_feature_qc.append(fqc_61)
    unlocked_dfs["GSE63061"] = df_61
    plot_microarray_boxplots(df_61, meta_61, "GSE63061", out_figs / "GSE63061_expression_boxplot.png")

    # 2. GSE63060 (Replication, Microarray)
    df_60, _ = ensure_samples_x_features_matrix("GSE63060", "expression.parquet", root)
    meta_60 = pd.read_csv(root / "data" / "interim" / "GSE63060" / "metadata.csv")
    sqc_60, flags_60 = compute_sample_qc_microarray(df_60, meta_60, "GSE63060", locked=False)
    fqc_60 = compute_feature_qc(df_60, "GSE63060", "microarray")
    all_sample_qc.append(sqc_60)
    all_sample_flags.extend(flags_60)
    all_feature_qc.append(fqc_60)
    unlocked_dfs["GSE63060"] = df_60
    plot_microarray_boxplots(df_60, meta_60, "GSE63060", out_figs / "GSE63060_expression_boxplot.png")

    # 3. GSE140829 (Holdout, Microarray, LOCKED)
    df_140, _ = ensure_samples_x_features_matrix("GSE140829", "expression.parquet", root)
    meta_140 = pd.read_csv(root / "data" / "interim" / "GSE140829" / "technical_metadata.csv")
    sqc_140, flags_140 = compute_sample_qc_microarray(df_140, meta_140, "GSE140829", locked=True)
    fqc_140 = compute_feature_qc(df_140, "GSE140829", "microarray")
    all_sample_qc.append(sqc_140)
    all_sample_flags.extend(flags_140)
    all_feature_qc.append(fqc_140)

    # 4. GSE125583 (Brain RNA-seq)
    df_125, _ = ensure_samples_x_features_matrix("GSE125583", "raw_counts.parquet", root)
    meta_125 = pd.read_csv(root / "data" / "interim" / "GSE125583" / "metadata.csv")
    sqc_125, flags_125 = compute_sample_qc_rnaseq(df_125, meta_125, "GSE125583")
    fqc_125 = compute_feature_qc(df_125, "GSE125583", "RNA-seq")
    all_sample_qc.append(sqc_125)
    all_sample_flags.extend(flags_125)
    all_feature_qc.append(fqc_125)
    unlocked_dfs["GSE125583"] = df_125
    plot_library_sizes_rnaseq(df_125, out_figs / "GSE125583_library_sizes.png")

    # Feature variance comparison figure across unlocked cohorts
    plot_feature_variance_distributions(unlocked_dfs, out_figs / "feature_variance_distribution.png")

    # Save summary tables
    combined_sample_qc = pd.concat(all_sample_qc, ignore_index=True)
    combined_sample_qc.to_csv(out_tables / "sample_qc_summary.csv", index=False)
    log.info("Saved: %s (%d rows)", (out_tables / "sample_qc_summary.csv").name, len(combined_sample_qc))

    flags_df = pd.DataFrame(all_sample_flags) if all_sample_flags else pd.DataFrame(
        columns=["sample_id", "dataset_id", "flag_reason", "flag_metric_value", "threshold_used", "provisional_action"]
    )
    flags_df.to_csv(out_tables / "sample_qc_flags.csv", index=False)
    log.info("Saved: %s (%d flagged samples)", (out_tables / "sample_qc_flags.csv").name, len(flags_df))

    fqc_df = pd.DataFrame(all_feature_qc)
    fqc_df.to_csv(out_tables / "feature_qc_summary.csv", index=False)
    log.info("Saved: %s (%d datasets)", (out_tables / "feature_qc_summary.csv").name, len(fqc_df))

    print("\n--- SAMPLE QC SUMMARY (BY COHORT) ---")
    print(combined_sample_qc.groupby(["dataset_id", "assay_type"]).agg({
        "sample_id": "count",
        "median_expr": ["mean", "min", "max"],
        "pct_missing": "mean",
        "is_outlier_flag": "sum",
    }))

    print("\n--- FEATURE QC SUMMARY ---")
    print(fqc_df[["dataset_id", "assay_type", "total_features", "constant_features", "near_zero_variance_features", "mean_feature_variance", "features_with_missing"]])


if __name__ == "__main__":
    main()
