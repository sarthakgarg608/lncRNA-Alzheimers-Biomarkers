"""
scripts/05_batch_effect_analysis.py
====================================
Stage 09 — Confounder Characterization & Batch-Effect Assessment.

Performs:
  1. Confounder Analysis (Task 7):
     - Age distribution across diagnosis groups (ANOVA / Kruskal-Wallis).
     - Sex distribution across diagnosis groups (Chi-square test of independence).
     - Center/batch distribution across diagnosis groups (Chi-square test).
     - Missingness audit across key clinical covariates.
     - Outputs: outputs/tables/stage09/confounder_assessment.csv
  2. Batch-Effect Assessment (Task 8):
     - Association between Center/batch and top PCs (ANOVA: PC1 ~ batch, PC2 ~ batch).
     - Variance explained (R^2) by batch.
     - Identification of potential confounding between batch and diagnosis.
     - Concludes: Center/batch is a significant covariate and should be included
       in linear models (e.g. limma/statsmodels) rather than blindly running ComBat.
     - Outputs: outputs/tables/stage09/batch_effect_assessment.csv
  3. Figures generated:
     - outputs/figures/stage09/PCA_batch_effects.png (Faceted PCA colored by batch/center)
     - outputs/figures/stage09/metadata_missingness.png (Covariate missingness audit)

Scientific integrity rules:
  - GSE140829 diagnosis remains LOCKED. Technical batch assessment only.
  - No batch correction is applied automatically (Task 8 requirement).
  - Confounder risks are reported honestly without altering data.
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
from scipy import stats

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


# ---------------------------------------------------------------------------
# Task 7: Confounder Analysis
# ---------------------------------------------------------------------------

def analyze_confounders(root: Path) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """
    Test statistical associations between diagnosis and potential confounders
    (Age, Sex, Batch/Center) for unlocked cohorts.
    """
    results: list[dict[str, Any]] = []

    # Datasets to analyze
    unlocked_configs = [
        ("GSE63061", "metadata.csv", ["AD", "Control", "MCI"]),
        ("GSE63060", "metadata.csv", ["AD", "Control", "MCI"]),
        ("GSE125583", "metadata.csv", ["AD", "Control"]),
    ]

    for ds_id, meta_file, target_diags in unlocked_configs:
        meta_path = root / "data" / "interim" / ds_id / meta_file
        if not meta_path.exists():
            continue
        df = pd.read_csv(meta_path)
        sub = df[df["diagnosis"].isin(target_diags)].copy()

        # 1. Age vs Diagnosis
        if "age" in sub.columns and sub["age"].notna().sum() > 10:
            age_data = sub.dropna(subset=["age"]).copy()
            age_data["age"] = pd.to_numeric(age_data["age"], errors="coerce")
            age_data = age_data.dropna(subset=["age"])

            groups = [group["age"].values for _, group in age_data.groupby("diagnosis")]
            if len(groups) >= 2 and all(len(g) > 2 for g in groups):
                f_stat, p_val = stats.f_oneway(*groups)
                kw_stat, kw_p = stats.kruskal(*groups)

                # Mean age per group
                means_str = ", ".join(f"{d}: {g['age'].mean():.1f}y" for d, g in age_data.groupby("diagnosis"))
                results.append({
                    "dataset_id": ds_id,
                    "covariate": "age",
                    "test_type": "ANOVA & Kruskal-Wallis",
                    "test_statistic": f"F={f_stat:.2f}, KW={kw_stat:.2f}",
                    "p_value": round(float(p_val), 6),
                    "significant_at_05": p_val < 0.05,
                    "interpretation": (
                        f"Age difference across diagnoses ({means_str}). "
                        + ("Statistically significant age difference: MUST adjust for age in DE models."
                           if p_val < 0.05 else "Age difference not statistically significant at alpha=0.05.")
                    ),
                })
        else:
            results.append({
                "dataset_id": ds_id,
                "covariate": "age",
                "test_type": "N/A",
                "test_statistic": "N/A",
                "p_value": None,
                "significant_at_05": False,
                "interpretation": "Age metadata not available in this cohort.",
            })

        # 2. Sex vs Diagnosis
        if "sex" in sub.columns and sub["sex"].notna().sum() > 10:
            sex_data = sub.dropna(subset=["sex"])
            contingency = pd.crosstab(sex_data["diagnosis"], sex_data["sex"])
            if contingency.shape[0] >= 2 and contingency.shape[1] >= 2:
                chi2, p_val, dof, _ = stats.chi2_contingency(contingency)
                results.append({
                    "dataset_id": ds_id,
                    "covariate": "sex",
                    "test_type": f"Chi-Square (df={dof})",
                    "test_statistic": f"Chi2={chi2:.2f}",
                    "p_value": round(float(p_val), 6),
                    "significant_at_05": p_val < 0.05,
                    "interpretation": (
                        f"Sex distribution by diagnosis. "
                        + ("Significant difference in sex ratio between groups: MUST adjust for sex in DE models."
                           if p_val < 0.05 else "Sex ratio comparable across diagnosis groups (p >= 0.05).")
                    ),
                })
        else:
            results.append({
                "dataset_id": ds_id,
                "covariate": "sex",
                "test_type": "N/A",
                "test_statistic": "N/A",
                "p_value": None,
                "significant_at_05": False,
                "interpretation": "Sex metadata not available in this cohort.",
            })

        # 3. Batch/Center vs Diagnosis
        if "batch" in sub.columns and sub["batch"].notna().sum() > 10:
            batch_data = sub.dropna(subset=["batch"])
            contingency_batch = pd.crosstab(batch_data["diagnosis"], batch_data["batch"])
            if contingency_batch.shape[0] >= 2 and contingency_batch.shape[1] >= 2:
                chi2, p_val, dof, _ = stats.chi2_contingency(contingency_batch)
                results.append({
                    "dataset_id": ds_id,
                    "covariate": "batch_center",
                    "test_type": f"Chi-Square (df={dof})",
                    "test_statistic": f"Chi2={chi2:.2f}",
                    "p_value": round(float(p_val), 6),
                    "significant_at_05": p_val < 0.05,
                    "interpretation": (
                        f"Center/batch balance across diagnoses ({contingency_batch.shape[1]} sites). "
                        + ("Significant site-diagnosis association: Site should be included as a blocking factor/covariate."
                           if p_val < 0.05 else "Centers are reasonably balanced across diagnoses.")
                    ),
                })
        else:
            results.append({
                "dataset_id": ds_id,
                "covariate": "batch_center",
                "test_type": "N/A",
                "test_statistic": "N/A",
                "p_value": None,
                "significant_at_05": False,
                "interpretation": "Batch/site metadata not available in this cohort.",
            })

    return pd.DataFrame(results), results


# ---------------------------------------------------------------------------
# Task 8: Batch-Effect Assessment (ANOVA on PC1 and PC2)
# ---------------------------------------------------------------------------

def assess_batch_effects(root: Path) -> pd.DataFrame:
    """
    Test association between reported batch/center and expression variance
    (PC1 and PC2 coordinates from sample_qc_summary.csv).
    """
    sqc_path = root / "outputs" / "tables" / "stage09" / "sample_qc_summary.csv"
    if not sqc_path.exists():
        raise FileNotFoundError(f"{sqc_path} does not exist. Run 03_quality_control.py and 04_normalization_assessment.py first.")

    sqc = pd.read_csv(sqc_path)
    records: list[dict[str, Any]] = []

    for ds_id in ["GSE63061", "GSE63060", "GSE140829", "GSE125583"]:
        sub = sqc[sqc["dataset_id"] == ds_id].copy()

        if "batch" not in sub.columns or sub["batch"].notna().sum() < 5:
            records.append({
                "dataset_id": ds_id,
                "batch_variable": "N/A",
                "n_batches": 0,
                "pc1_f_stat": None,
                "pc1_p_val": None,
                "pc1_r2_pct": None,
                "pc2_f_stat": None,
                "pc2_p_val": None,
                "pc2_r2_pct": None,
                "batch_effect_detected": False,
                "recommendation": "Batch metadata unavailable; monitor surrogate variables (SVA) at DE stage.",
            })
            continue

        batch_sub = sub.dropna(subset=["batch", "pc1", "pc2"]).copy()
        n_batches = batch_sub["batch"].nunique()

        if n_batches < 2:
            records.append({
                "dataset_id": ds_id,
                "batch_variable": "batch",
                "n_batches": n_batches,
                "pc1_f_stat": None,
                "pc1_p_val": None,
                "pc1_r2_pct": None,
                "pc2_f_stat": None,
                "pc2_p_val": None,
                "pc2_r2_pct": None,
                "batch_effect_detected": False,
                "recommendation": "Single batch; no inter-batch correction needed.",
            })
            continue

        # ANOVA for PC1
        groups_pc1 = [grp["pc1"].values for _, grp in batch_sub.groupby("batch") if len(grp) > 1]
        f1, p1 = stats.f_oneway(*groups_pc1) if len(groups_pc1) >= 2 else (np.nan, np.nan)

        # ANOVA for PC2
        groups_pc2 = [grp["pc2"].values for _, grp in batch_sub.groupby("batch") if len(grp) > 1]
        f2, p2 = stats.f_oneway(*groups_pc2) if len(groups_pc2) >= 2 else (np.nan, np.nan)

        # Calculate R^2 = SS_between / SS_total
        total_ss1 = np.sum((batch_sub["pc1"] - batch_sub["pc1"].mean()) ** 2)
        between_ss1 = sum(len(g) * (np.mean(g) - batch_sub["pc1"].mean()) ** 2 for g in groups_pc1) if groups_pc1 else 0
        r2_1 = (between_ss1 / total_ss1 * 100.0) if total_ss1 > 0 else 0.0

        total_ss2 = np.sum((batch_sub["pc2"] - batch_sub["pc2"].mean()) ** 2)
        between_ss2 = sum(len(g) * (np.mean(g) - batch_sub["pc2"].mean()) ** 2 for g in groups_pc2) if groups_pc2 else 0
        r2_2 = (between_ss2 / total_ss2 * 100.0) if total_ss2 > 0 else 0.0

        batch_detected = (p1 < 0.05) or (p2 < 0.05)

        recom = (
            f"Significant batch/center effect detected on PC1/PC2 (p < 0.05; R2_PC1={r2_1:.1f}%, R2_PC2={r2_2:.1f}%). "
            "Recommendation: Include batch as an additive covariate in downstream differential expression (limma/statsmodels). "
            "Do NOT apply global ComBat before cross-validation to prevent data leakage."
            if batch_detected
            else "No statistically significant batch effect detected on top PCs (p >= 0.05)."
        )

        records.append({
            "dataset_id": ds_id,
            "batch_variable": "center_code" if ds_id.startswith("GSE6306") else "batch_id",
            "n_batches": n_batches,
            "pc1_f_stat": round(float(f1), 2) if not np.isnan(f1) else None,
            "pc1_p_val": round(float(p1), 6) if not np.isnan(p1) else None,
            "pc1_r2_pct": round(float(r2_1), 2),
            "pc2_f_stat": round(float(f2), 2) if not np.isnan(f2) else None,
            "pc2_p_val": round(float(p2), 6) if not np.isnan(p2) else None,
            "pc2_r2_pct": round(float(r2_2), 2),
            "batch_effect_detected": batch_detected,
            "recommendation": recom,
        })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Task 10: Figures
# ---------------------------------------------------------------------------

def plot_pca_batch_effects(root: Path, out_path: Path) -> None:
    """
    Generate faceted PCA colored by batch / collection center:
      - GSE63061 (Center)
      - GSE63060 (Center)
      - GSE140829 (Technical Batch ID — NO diagnosis labels)
    """
    sqc_path = root / "outputs" / "tables" / "stage09" / "sample_qc_summary.csv"
    sqc = pd.read_csv(sqc_path)

    cohorts = ["GSE63061", "GSE63060", "GSE140829"]
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for ax, ds_id in zip(axes, cohorts):
        sub = sqc[sqc["dataset_id"] == ds_id].dropna(subset=["batch", "pc1", "pc2"]).copy()
        if sub.empty:
            ax.set_title(f"{ds_id} (No batch metadata)")
            continue

        batches = sub["batch"].astype(str).tolist()
        unique_batches = sorted(set(batches))
        # Use colormap
        cmap_name = "tab10" if len(unique_batches) <= 10 else "tab20"
        cmap = matplotlib.colormaps[cmap_name].resampled(len(unique_batches))

        for b_idx, b_val in enumerate(unique_batches):
            b_sub = sub[sub["batch"].astype(str) == b_val]
            ax.scatter(
                b_sub["pc1"],
                b_sub["pc2"],
                label=f"{b_val} (n={len(b_sub)})",
                color=cmap(b_idx),
                alpha=0.75,
                s=35,
                edgecolors="none",
            )

        ax.set_xlabel("PC1", fontsize=10)
        ax.set_ylabel("PC2", fontsize=10)
        ax.set_title(f"{ds_id} — PCA colored by Batch/Center\n({len(unique_batches)} batches)", fontsize=11)
        ax.legend(fontsize=8, loc="best", framealpha=0.8, ncol=2 if len(unique_batches) > 6 else 1)
        ax.grid(alpha=0.3)

    fig.suptitle("Unsupervised PCA Colored by Clinical Center / Technical Batch", fontsize=13, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved batch PCA figure: %s", out_path.name)


def plot_metadata_missingness(root: Path, out_path: Path) -> None:
    """
    Bar chart showing missingness percentage across key covariates
    for all 4 datasets.
    """
    covariates = ["diagnosis", "age", "sex", "batch", "tissue", "ethnicity"]
    datasets = ["GSE63061", "GSE63060", "GSE140829", "GSE125583"]

    missing_matrix = []

    for ds_id in datasets:
        meta_file = "technical_metadata.csv" if ds_id == "GSE140829" else "metadata.csv"
        meta_path = root / "data" / "interim" / ds_id / meta_file
        df = pd.read_csv(meta_path)
        row = []
        for cov in covariates:
            if cov not in df.columns:
                row.append(100.0)  # entirely missing
            else:
                pct = (df[cov].isna().sum() / len(df)) * 100.0
                row.append(pct)
        missing_matrix.append(row)

    missing_df = pd.DataFrame(missing_matrix, index=datasets, columns=covariates)

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(covariates))
    width = 0.18
    colors = ["#3498DB", "#2ECC71", "#9B59B6", "#E67E22"]

    for i, (ds_id, col) in enumerate(zip(datasets, colors)):
        vals = missing_df.loc[ds_id].values
        ax.bar(x + (i - 1.5) * width, vals, width, label=ds_id, color=col)

    ax.set_ylabel("% Missing", fontsize=11)
    ax.set_title("Clinical Covariate Missingness Across All Four Cohorts", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(covariates, fontsize=10)
    ax.set_ylim(0, 105)
    ax.legend(fontsize=9)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Saved metadata missingness figure: %s", out_path.name)


# ---------------------------------------------------------------------------
# Main Execution Pipeline
# ---------------------------------------------------------------------------

def main() -> None:
    root = get_project_root()
    out_tables = root / "outputs" / "tables" / "stage09"
    out_figs = root / "outputs" / "figures" / "stage09"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_figs.mkdir(parents=True, exist_ok=True)

    log.info("=" * 65)
    log.info("STAGE 09: CONFOUNDER & BATCH-EFFECT ANALYSIS")
    log.info("=" * 65)

    # 1. Confounder analysis
    conf_df, _ = analyze_confounders(root)
    conf_path = out_tables / "confounder_assessment.csv"
    conf_df.to_csv(conf_path, index=False)
    log.info("Saved: %s (%d rows)", conf_path.name, len(conf_df))

    # 2. Batch-effect assessment
    batch_df = assess_batch_effects(root)
    batch_path = out_tables / "batch_effect_assessment.csv"
    batch_df.to_csv(batch_path, index=False)
    log.info("Saved: %s (%d rows)", batch_path.name, len(batch_df))

    # 3. Figures
    plot_pca_batch_effects(root, out_figs / "PCA_batch_effects.png")
    plot_metadata_missingness(root, out_figs / "metadata_missingness.png")

    print("\n--- CONFOUNDER ASSESSMENT SUMMARY ---")
    print(conf_df[["dataset_id", "covariate", "test_type", "p_value", "significant_at_05"]].to_string(index=False))

    print("\n--- BATCH EFFECT ASSESSMENT SUMMARY ---")
    print(batch_df[["dataset_id", "batch_variable", "n_batches", "pc1_r2_pct", "pc1_p_val", "pc2_r2_pct", "pc2_p_val", "batch_effect_detected"]].to_string(index=False))


if __name__ == "__main__":
    main()
