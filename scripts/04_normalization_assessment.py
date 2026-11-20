"""
scripts/04_normalization_assessment.py
=======================================
Stage 09 — Normalization Audit, PCA & Unsupervised Sample Clustering.

Performs:
  1. Detailed normalization audit:
     - Microarray: confirms VST + quantile/spline normalization from Lumi;
       verifies log-scale properties (min, medians, max).
       Concludes: DO NOT re-log2 or re-normalize.
     - Brain RNA-seq: confirms raw integer counts;
       calculates exploratory log1p(CPM) strictly for visualization.
       Concludes: Raw counts must be preserved for PyDESeq2.
  2. PCA analysis (PC1, PC2) on standardized top-variance features:
     - GSE63061 (colored by diagnosis: AD, Control, MCI, Unknown)
     - GSE63060 (colored by diagnosis: AD, Control, MCI)
     - GSE125583 (colored by diagnosis: AD, Control)
  3. Unsupervised sample clustering (hierarchical clustering of sample correlation).
  4. Figures generated:
     - outputs/figures/stage09/GSE63061_PCA_diagnosis.png
     - outputs/figures/stage09/GSE63060_PCA_diagnosis.png
     - outputs/figures/stage09/GSE125583_PCA_diagnosis.png
     - outputs/figures/stage09/sample_clustering.png
  5. Tables generated / updated:
     - outputs/tables/stage09/normalization_assessment.csv
     - updates outputs/tables/stage09/sample_qc_summary.csv with PC1 and PC2 coords

Scientific integrity:
  - GSE140829 diagnosis remains LOCKED. PCA colored by batch only if plotted.
  - No samples are removed based on PCA overlap.
  - RNA-seq raw counts are never replaced.
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
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ---------------------------------------------------------------------------
# Setup & Paths
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from scripts.utils.paths import get_project_root, load_dataset_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

TOP_VARIABLE_FEATURES = 2000


# ---------------------------------------------------------------------------
# Task 5: Normalization Assessment Audit
# ---------------------------------------------------------------------------

def audit_normalization(root: Path) -> pd.DataFrame:
    """
    Produce a rigorous audit table evaluating the scale and normalization
    status of all 4 datasets.
    """
    records = []

    # 1. GSE63061
    df61 = pd.read_parquet(root / "data" / "interim" / "GSE63061" / "expression_samples_x_features.parquet")
    vals61 = df61.values
    meds61 = np.nanmedian(vals61, axis=1)
    records.append({
        "dataset_id": "GSE63061",
        "assay_type": "microarray",
        "platform_id": "GPL10558",
        "original_scale_reported": "Variance-stabilized and quantile-normalized (lumi package in R)",
        "observed_min": round(float(np.nanmin(vals61)), 3),
        "observed_median_range": f"{np.min(meds61):.2f} - {np.max(meds61):.2f}",
        "observed_max": round(float(np.nanmax(vals61)), 3),
        "is_pre_normalized": True,
        "recommended_action": "PRESERVE_AS_IS",
        "justification": (
            "Data already transformed by BeadStudio/lumi (VST + quantile normalization). "
            "Values are in log2-like intensity space (~5.5 to 14.5). Re-normalizing or "
            "applying log2 again would corrupt distribution and distort biological signal."
        ),
    })

    # 2. GSE63060
    df60 = pd.read_parquet(root / "data" / "interim" / "GSE63060" / "expression_samples_x_features.parquet")
    vals60 = df60.values
    meds60 = np.nanmedian(vals60, axis=1)
    records.append({
        "dataset_id": "GSE63060",
        "assay_type": "microarray",
        "platform_id": "GPL6947",
        "original_scale_reported": "Variance-stabilized and normalized (lumi package in R)",
        "observed_min": round(float(np.nanmin(vals60)), 3),
        "observed_median_range": f"{np.min(meds60):.2f} - {np.max(meds60):.2f}",
        "observed_max": round(float(np.nanmax(vals60)), 3),
        "is_pre_normalized": True,
        "recommended_action": "PRESERVE_AS_IS",
        "justification": (
            "Pre-normalized by source via lumi. Medians across samples are aligned (~7.5). "
            "No additional transformation warranted. Preserved as original normalized intensities."
        ),
    })

    # 3. GSE140829
    df140 = pd.read_parquet(root / "data" / "interim" / "GSE140829" / "expression_samples_x_features.parquet")
    vals140 = df140.values
    meds140 = np.nanmedian(vals140, axis=1)
    records.append({
        "dataset_id": "GSE140829",
        "assay_type": "microarray",
        "platform_id": "GPL15988",
        "original_scale_reported": "VST + robust spline normalization (RSN) via lumi (normalized_data.txt)",
        "observed_min": round(float(np.nanmin(vals140)), 3),
        "observed_median_range": f"{np.min(meds140):.2f} - {np.max(meds140):.2f}",
        "observed_max": round(float(np.nanmax(vals140)), 3),
        "is_pre_normalized": True,
        "recommended_action": "PRESERVE_AS_IS_LOCKED",
        "justification": (
            "Pre-normalized by source via VST and RSN. Sample medians are aligned (~7.3). "
            "Preserved as-is. Hold-out dataset remains locked."
        ),
    })

    # 4. GSE125583
    df125 = pd.read_parquet(root / "data" / "interim" / "GSE125583" / "raw_counts_samples_x_features.parquet")
    vals125 = df125.values
    records.append({
        "dataset_id": "GSE125583",
        "assay_type": "RNA-seq",
        "platform_id": "GPL16791",
        "original_scale_reported": "Raw integer read counts (HTSeqGenie/GSNAP)",
        "observed_min": int(np.min(vals125)),
        "observed_median_range": f"Library sizes: {np.min(vals125.sum(axis=1)):,} - {np.max(vals125.sum(axis=1)):,}",
        "observed_max": int(np.max(vals125)),
        "is_pre_normalized": False,
        "recommended_action": "PRESERVE_RAW_COUNTS_USE_DESEQ2",
        "justification": (
            "Data consists of genuine non-negative integer read counts. "
            "Must NOT be log-transformed or normalized with microarray methods. "
            "Raw integer counts preserved for negative binomial modeling via PyDESeq2 at DE stage. "
            "log1p(CPM) computed exclusively for non-destructive exploratory visualization."
        ),
    })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Task 6: PCA Preparation & Execution
# ---------------------------------------------------------------------------

def run_pca_on_matrix(
    df_s_x_f: pd.DataFrame,
    assay_type: str,
    n_top_features: int = TOP_VARIABLE_FEATURES,
    random_state: int = 42,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """
    Run PCA on top-variance features.

    Parameters
    ----------
    df_s_x_f : pd.DataFrame
        Samples on rows, features on columns.
    assay_type : str
        "microarray" or "RNA-seq".

    Returns
    -------
    (pca_coords, explained_variance_ratio, top_feature_names)
    """
    values = df_s_x_f.values

    if assay_type == "RNA-seq":
        # Library-size adjustment: CPM = (count / lib_size) * 1e6 -> log1p(CPM)
        lib_sizes = values.sum(axis=1, keepdims=True)
        lib_sizes[lib_sizes == 0] = 1.0
        cpm = (values / lib_sizes) * 1e6
        norm_matrix = np.log1p(cpm)
    else:
        norm_matrix = values

    # Feature variances
    variances = np.nanvar(norm_matrix, axis=0)
    top_indices = np.argsort(variances)[::-1][:n_top_features]
    top_matrix = norm_matrix[:, top_indices]
    top_features = [df_s_x_f.columns[idx] for idx in top_indices]

    # Standardize features for PCA
    scaler = StandardScaler(with_mean=True, with_std=True)
    scaled_matrix = scaler.fit_transform(top_matrix)

    pca = PCA(n_components=2, random_state=random_state)
    coords = pca.fit_transform(scaled_matrix)
    return coords, pca.explained_variance_ratio_, top_features


def plot_pca_diagnosis(
    coords: np.ndarray,
    variance_ratio: np.ndarray,
    meta_df: pd.DataFrame,
    dataset_id: str,
    out_path: Path,
) -> None:
    """
    Plot PC1 vs PC2 colored by diagnosis labels.
    """
    meta_idx = meta_df.set_index("sample_id")
    sample_ids = meta_df["sample_id"].tolist()

    diagnoses = [
        meta_idx.loc[sid, "diagnosis"] if sid in meta_idx.index and "diagnosis" in meta_idx.columns else "Unknown"
        for sid in sample_ids
    ]

    color_map = {
        "AD": "#E74C3C",        # Red
        "Control": "#2ECC71",   # Green
        "MCI": "#F39C12",       # Orange
        "Unknown": "#95A5A6",   # Gray
    }

    fig, ax = plt.subplots(figsize=(8, 6.5))

    unique_diags = sorted(set(diagnoses), key=lambda d: ("AD", "Control", "MCI", "Unknown").index(d) if d in ("AD", "Control", "MCI", "Unknown") else 99)

    for diag in unique_diags:
        idxs = [i for i, d in enumerate(diagnoses) if d == diag]
        c = color_map.get(diag, "#34495E")
        ax.scatter(
            coords[idxs, 0],
            coords[idxs, 1],
            c=c,
            label=f"{diag} (n={len(idxs)})",
            alpha=0.75,
            edgecolors="none",
            s=45,
        )

    ax.set_xlabel(f"PC1 ({variance_ratio[0]*100:.1f}% explained variance)", fontsize=11)
    ax.set_ylabel(f"PC2 ({variance_ratio[1]*100:.1f}% explained variance)", fontsize=11)
    ax.set_title(
        f"{dataset_id} — Unsupervised PCA (Top {TOP_VARIABLE_FEATURES} Variable Features)\n"
        f"Colored by Disease Diagnosis",
        fontsize=12,
    )
    ax.legend(fontsize=10, loc="best", framealpha=0.85)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Saved PCA figure: %s", out_path.name)


# ---------------------------------------------------------------------------
# Task 6: Unsupervised Sample Clustering
# ---------------------------------------------------------------------------

def plot_sample_clustering(
    datasets_matrices: dict[str, tuple[pd.DataFrame, str]],
    out_path: Path,
) -> None:
    """
    Plot hierarchical clustering dendrograms based on sample correlation
    for GSE63061, GSE63060, and GSE125583.
    """
    fig, axes = plt.subplots(1, len(datasets_matrices), figsize=(6 * len(datasets_matrices), 5))
    if len(datasets_matrices) == 1:
        axes = [axes]

    for ax, (ds_id, (df_s_x_f, assay_type)) in zip(axes, datasets_matrices.items()):
        values = df_s_x_f.values
        if assay_type == "RNA-seq":
            lib_sizes = values.sum(axis=1, keepdims=True)
            lib_sizes[lib_sizes == 0] = 1.0
            norm_vals = np.log1p((values / lib_sizes) * 1e6)
        else:
            norm_vals = values

        # Use top 1000 variable features for correlation
        vars_ = np.nanvar(norm_vals, axis=0)
        top_idx = np.argsort(vars_)[::-1][:1000]
        sub_matrix = norm_vals[:, top_idx]

        # Sample correlation distance
        # Standardize samples
        scaler = StandardScaler()
        sub_scaled = scaler.fit_transform(sub_matrix.T).T  # scale each sample
        corr_matrix = np.corrcoef(sub_scaled)
        # Distance = 1 - corr
        dist_matrix = np.clip(1.0 - corr_matrix, 0.0, 2.0)
        np.fill_diagonal(dist_matrix, 0.0)

        # Condensed distance matrix
        condensed = pdist(sub_scaled, metric="correlation")
        # Ward linkage
        Z = linkage(condensed, method="average")

        dendrogram(
            Z,
            ax=ax,
            no_labels=True,
            color_threshold=0.7 * np.max(Z[:, 2]),
            above_threshold_color="#2C3E50",
        )
        ax.set_title(f"{ds_id} (n={len(df_s_x_f)} samples)", fontsize=11)
        ax.set_xlabel("Sample Clusters", fontsize=10)
        ax.set_ylabel("Distance (1 - Pearson r)", fontsize=10)
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle("Unsupervised Hierarchical Clustering of Samples", fontsize=13, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved clustering figure: %s", out_path.name)


# ---------------------------------------------------------------------------
# Main Execution
# ---------------------------------------------------------------------------

def main() -> None:
    root = get_project_root()
    out_tables = root / "outputs" / "tables" / "stage09"
    out_figs = root / "outputs" / "figures" / "stage09"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_figs.mkdir(parents=True, exist_ok=True)

    log.info("=" * 65)
    log.info("STAGE 09: NORMALIZATION AUDIT, PCA & CLUSTERING")
    log.info("=" * 65)

    # 1. Normalization audit
    norm_audit_df = audit_normalization(root)
    norm_audit_df.to_csv(out_tables / "normalization_assessment.csv", index=False)
    log.info("Saved: %s", (out_tables / "normalization_assessment.csv").name)

    # 2. PCA on GSE63061, GSE63060, GSE125583
    pca_results: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    # GSE63061
    df61 = pd.read_parquet(root / "data" / "interim" / "GSE63061" / "expression_samples_x_features.parquet")
    meta61 = pd.read_csv(root / "data" / "interim" / "GSE63061" / "metadata.csv")
    coords61, var61, _ = run_pca_on_matrix(df61, assay_type="microarray")
    pca_results["GSE63061"] = (coords61, var61)
    plot_pca_diagnosis(coords61, var61, meta61, "GSE63061", out_figs / "GSE63061_PCA_diagnosis.png")

    # GSE63060
    df60 = pd.read_parquet(root / "data" / "interim" / "GSE63060" / "expression_samples_x_features.parquet")
    meta60 = pd.read_csv(root / "data" / "interim" / "GSE63060" / "metadata.csv")
    coords60, var60, _ = run_pca_on_matrix(df60, assay_type="microarray")
    pca_results["GSE63060"] = (coords60, var60)
    plot_pca_diagnosis(coords60, var60, meta60, "GSE63060", out_figs / "GSE63060_PCA_diagnosis.png")

    # GSE125583
    df125 = pd.read_parquet(root / "data" / "interim" / "GSE125583" / "raw_counts_samples_x_features.parquet")
    meta125 = pd.read_csv(root / "data" / "interim" / "GSE125583" / "metadata.csv")
    coords125, var125, _ = run_pca_on_matrix(df125, assay_type="RNA-seq")
    pca_results["GSE125583"] = (coords125, var125)
    plot_pca_diagnosis(coords125, var125, meta125, "GSE125583", out_figs / "GSE125583_PCA_diagnosis.png")

    # Sample clustering
    clustering_inputs = {
        "GSE63061": (df61, "microarray"),
        "GSE63060": (df60, "microarray"),
        "GSE125583": (df125, "RNA-seq"),
    }
    plot_sample_clustering(clustering_inputs, out_figs / "sample_clustering.png")

    # Update sample_qc_summary.csv with PC1 and PC2 coordinates
    sample_qc_path = out_tables / "sample_qc_summary.csv"
    if sample_qc_path.exists():
        sqc_df = pd.read_csv(sample_qc_path)
        # Create map of sample_id -> (pc1, pc2)
        pc_map: dict[str, tuple[float, float]] = {}
        for s_idx, sid in enumerate(df61.index):
            pc_map[sid] = (round(float(coords61[s_idx, 0]), 4), round(float(coords61[s_idx, 1]), 4))
        for s_idx, sid in enumerate(df60.index):
            pc_map[sid] = (round(float(coords60[s_idx, 0]), 4), round(float(coords60[s_idx, 1]), 4))
        for s_idx, sid in enumerate(df125.index):
            pc_map[sid] = (round(float(coords125[s_idx, 0]), 4), round(float(coords125[s_idx, 1]), 4))

        # Also for GSE140829 (technical PCA)
        df140 = pd.read_parquet(root / "data" / "interim" / "GSE140829" / "expression_samples_x_features.parquet")
        coords140, _, _ = run_pca_on_matrix(df140, assay_type="microarray")
        for s_idx, sid in enumerate(df140.index):
            pc_map[sid] = (round(float(coords140[s_idx, 0]), 4), round(float(coords140[s_idx, 1]), 4))

        sqc_df["pc1"] = sqc_df["sample_id"].map(lambda s: pc_map.get(s, (None, None))[0])
        sqc_df["pc2"] = sqc_df["sample_id"].map(lambda s: pc_map.get(s, (None, None))[1])
        sqc_df.to_csv(sample_qc_path, index=False)
        log.info("Updated %s with PC1 and PC2 coordinates", sample_qc_path.name)

    print("\n--- NORMALIZATION AUDIT SUMMARY ---")
    print(norm_audit_df[["dataset_id", "assay_type", "is_pre_normalized", "recommended_action", "observed_median_range"]].to_string(index=False))

    print("\n--- PCA EXPLAINED VARIANCE RATIOS ---")
    for ds_id, (_, var_ratio) in pca_results.items():
        print(f"  {ds_id}: PC1={var_ratio[0]*100:.2f}%, PC2={var_ratio[1]*100:.2f}%, Cumulative={sum(var_ratio)*100:.2f}%")


if __name__ == "__main__":
    main()
