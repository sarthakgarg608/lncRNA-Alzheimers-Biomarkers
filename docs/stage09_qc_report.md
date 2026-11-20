# Stage 09 — Quality Control, Normalization Assessment & Batch-Effect Characterization

**Project:** lncSignal-AD — Establishing the Role of lncRNAs as Biomarkers in Alzheimer's Disease  
**Stage:** 09 — Quality Control, Normalization Assessment & Batch-Effect Characterization  
**Date:** 2026-11-21  
**Python Version:** 3.12.10  
**Status:** ✅ CONDITIONAL PASS

---

> Stage 10 corrections to flag counts, ages and batch interpretation are documented in [the preflight audit](stage10_preflight_audit.md).

## 1. Executive Summary

Stage 09 performs rigorous quality control, normalization auditing, PCA-based exploration,
unsupervised clustering, confounder characterization, and batch-effect assessment across all
four GEO datasets. No samples or features have been removed. All 178 flagged samples (202 events) are
provisionally marked for inspection only.

**Verdict: CONDITIONAL PASS — Datasets are scientifically usable for downstream analysis.
Significant confounders (age, center) must be incorporated as covariates in all linear models.**

> [!IMPORTANT]
> GSE140829 diagnosis labels remain **LOCKED** throughout Stage 09. Only technical QC metrics
> and batch identifiers were accessed for this dataset.

---

## 2. Datasets Analysed

| Dataset | Samples | Features | Assay Type | Role |
|---|---:|---:|---|---|
| GSE63061 | 388 | 32,049 | Blood microarray (Illumina HT-12 v3) | Primary discovery |
| GSE63060 | 329 | 38,323 | Blood microarray (Illumina HT-12 v4) | Replication |
| GSE140829 | 587 | 47,008 | Blood microarray (GPL15988) | External validation — **LOCKED** |
| GSE125583 | 289 | 39,376 | Brain RNA-seq (NCBI Entrez IDs) | Biological support |

---

## 3. Matrix Orientation (Task 2)

### Finding

All four Parquet files from Stage 08 were stored in GEO's original orientation:
**features × samples** (probes/genes on rows, samples on columns).

| Dataset | Original Shape (Stage 08) | Derived Transposed Shape |
|---|---|---|
| GSE63061 | (32,049 × 388) | (388 × 32,049) ✅ |
| GSE63060 | (38,323 × 329) | (329 × 38,323) ✅ |
| GSE140829 | (47,008 × 587) | (587 × 47,008) ✅ |
| GSE125583 | (39,376 × 289) | (289 × 39,376) ✅ |

**Action taken:** Derived transposed matrices (`*_samples_x_features.parquet`) created in
`data/interim/<dataset>/` for all downstream analysis. **Original Parquet files are preserved
unchanged.**

---

## 4. Sample-Level QC (Task 3)

### 4.1 Summary Statistics by Dataset

**GSE63061 (Blood microarray, 388 samples)**

- Expression range: ~5.5–14.5 (log-scale, lumi VST + quantile normalized)
- Sample medians: 6.14–6.18 (tightly aligned — confirming normalization quality)
- Mean ± SD per sample: ~6.6 ± 1.6
- Missing values: **0** across all 388 samples
- Constant features: **0**
- Flagged samples: **58** (outlier flags only — NOT deleted)

**GSE63060 (Blood microarray, 329 samples)**

- Expression range: ~5.5–13.5
- Sample medians: 7.54–7.54 (perfectly aligned)
- Mean ± SD per sample: ~7.9 ± 1.5
- Missing values: **0** across all 329 samples
- Constant features: **0**
- Flagged samples: **20**

**GSE140829 (Blood microarray, 587 samples — LOCKED)**

- Expression range: ~4.4–14.5
- Sample medians: 7.32–7.38 (tightly aligned, VST + RSN normalization)
- Technical batch: 19 unique batches (batch explains 36–43% of PC variance)
- Missing values: **0**
- Flagged samples: **47** (technical QC only — no diagnosis label accessed)

**GSE125583 (Brain RNA-seq, 289 samples)**

- Library sizes: 8,057,646 – 93,556,289 (mean ~40M, expected variation for RNA-seq)
- Zero-count genes: 1,181 genes with 0 counts in all samples
- Genes expressed in <5% of samples: 5,524
- Flagged samples: **53** (library size or zero fraction)

### 4.2 Total Flagged Samples

| Dataset | Flagged (of total) | Flag Reasons |
|---|---|---|
| GSE63061 | 58 / 388 (14.9%) | Extreme median intensity, extreme IQR |
| GSE63060 | 20 / 329 (6.1%) | Extreme median intensity |
| GSE140829 | 47 / 587 (8.0%) | Extreme median intensity, extreme IQR (technical only) |
| GSE125583 | 53 / 289 (18.3%) | Extreme library size, high zero fraction |
| **Total** | **178 / 1,593 (11.2%)** | |

> [!NOTE]
> Flagged samples are retained in all interim files. Removal decisions must be made
> explicitly after biological review in Stage 10. The `sample_qc_flags.csv` table documents
> each flag with its metric value and threshold.

---

## 5. Feature-Level QC (Task 4)

| Dataset | Total Features | Features with Missing | Constant Features | Mean Feature Variance |
|---|---:|---:|---:|---|
| GSE63061 | 32,049 | 0 | 0 | 0.0263 |
| GSE63060 | 38,323 | 0 | 0 | 0.0289 |
| GSE140829 | 47,008 | 0 | 0 | 0.0360 |
| GSE125583 | 39,376 | 0 | 1,181 | 1.05 × 10⁸ |

**Key finding:** Microarray datasets have 0 missing and 0 constant features — extremely clean data.
RNA-seq (GSE125583) has 1,181 always-zero genes (3.0% of features), which is expected for
brain RNA-seq and will be filtered at the differential expression stage.

---

## 6. Normalization Assessment (Task 5)

### 6.1 Verdict by Dataset

| Dataset | Assay Type | Normalization Method | `is_pre_normalized` | Action |
|---|---|---|---|---|
| GSE63061 | Microarray | VST + quantile normalization (lumi R package) | ✅ True | **DO NOT re-normalize** |
| GSE63060 | Microarray | VST + quantile normalization (lumi R package) | ✅ True | **DO NOT re-normalize** |
| GSE140829 | Microarray | VST + robust spline normalization (lumi RSN) | ✅ True | **PRESERVE AS IS (LOCKED)** |
| GSE125583 | RNA-seq | Raw integer counts (HTSeqGenie/GSNAP) | ❌ False | **PRESERVE RAW COUNTS → PyDESeq2** |

### 6.2 Evidence for Normalization Status

**GSE63061/GSE63060:**
Values range ~5.5–14.5, consistent with variance-stabilized log-scale output from BeadStudio/lumi.
Sample medians are tightly aligned across all samples (GSE63061: 6.14–6.18; GSE63060: 7.54–7.54),
a hallmark of quantile normalization. **Re-applying log2 would corrupt the signal** (median would
drop from ~7 to ~2.8).

**GSE140829:**
VST + RSN normalization is confirmed by tight median alignment (7.32–7.38). The sample metadata
contains a `normalized_data.txt` source file, confirming pre-processing by lumi RSN.

**GSE125583:**
Values are non-negative integers (0 to 500,000+), library sizes in the millions, and 1,181
always-zero genes — all consistent with raw RNA-seq read counts. These **must not** be
log-transformed or normalized using microarray methods. Raw integer counts will be fed
directly to PyDESeq2 for differential expression analysis.

> [!CAUTION]
> Never pass processed microarray intensities into PyDESeq2 (which expects non-negative integers).
> Never re-log2 or re-quantile-normalize already-normalized microarray data.

---

## 7. PCA Analysis (Task 6)

### 7.1 Explained Variance

| Dataset | PC1 (%) | PC2 (%) | PC1+PC2 Total | Dominant Structure |
|---|---:|---:|---:|---|
| GSE63061 | 33.5 | 11.0 | 44.5% | AD vs. Control partially separated |
| GSE63060 | 26.1 | 13.1 | 39.2% | Moderate AD/Control separation |
| GSE125583 | 46.4 | 24.5 | 70.9% | Strong AD/Control separation in brain |

(GSE140829 PCA colored by batch only — no diagnosis plotted.)

### 7.2 Observations

- **GSE63061:** PC1 captures 33.5% variance. Partial AD/Control separation visible. MCI samples
  overlap with both groups (expected). No extreme outlier clusters.

- **GSE63060:** PC1 explains 26.1%. Moderate separation between AD and Control.
  GSE63060 is a replication cohort; similar structure to GSE63061 supports cross-cohort consistency.

- **GSE125583 (brain):** Highest explained variance (70.9% in 2 PCs), consistent with strong
  biological signal in brain tissue. AD samples cluster distinctly from Controls.

- **GSE140829:** PCA colored by technical batch. PC1 and PC2 heavily structured by batch
  (36% and 43% variance respectively) — showing batch-associated PC structure. **Diagnosis labels not accessed.**

### 7.3 PCA Input Preparation

- **Microarray:** Top 2,000 highest-variance probes selected → StandardScaler applied → PCA(n=2).
- **RNA-seq:** log1p(CPM) transformation applied to raw counts (non-destructive) → top-variance
  genes selected → StandardScaler → PCA(n=2). Raw counts preserved unchanged.

---

## 8. Unsupervised Sample Clustering (Task 6)

Hierarchical clustering of pairwise Pearson correlations among a random subsample of 50 samples
per dataset was performed (Ward linkage on 1-correlation distance).

**Key observations:**
- GSE63061/GSE63060: Samples cluster broadly by diagnosis, with some intermixing (expected for
  blood, which has modest signal compared to brain).
- GSE125583: Stronger cluster structure by diagnosis, consistent with PCA.
- GSE140829: Clustering was examined by technical batch. Diagnosis association was not tested.

---

## 9. Potential Outliers (Task 3, provisional)

202 sample flag events recorded (a sample can receive multiple flags):

| Dataset | Flag Events | Unique Samples Flagged |
|---|---:|---:|
| GSE63061 | 64 flag events | 58 unique samples |
| GSE63060 | 20 flag events | 20 unique samples |
| GSE140829 | 55 flag events | 47 unique samples |
| GSE125583 | 63 flag events | 53 unique samples |

**Flagging thresholds used:**
- Median intensity outlier: |median_i − population_median| > 3.5 × MAD
- IQR outlier: |IQR_i − population_IQR_median| > 3.5 × MAD
- RNA-seq: library size > 3.5 × MAD from population library size median

> [!WARNING]
> These are provisional flags for human review. No samples are automatically excluded.
> Downstream analysis should be performed with and without flagged samples as a sensitivity check.

---

## 10. Confounder Assessment (Task 7)

### 10.1 Age

| Dataset | Test | p-value | Significant | Interpretation |
|---|---|---|---|---|
| GSE63061 | ANOVA + KW | 0.000788 | ✅ Yes | AD older than Control (mean: AD 77.9y, Ctrl 75.3y) |
| GSE63060 | ANOVA + KW | 0.001152 | ✅ Yes | Significant age difference across groups |
| GSE125583 | ANOVA + KW | 0.000737 | ✅ Yes | AD mean 83.4 years, Control mean 87.0 years (288 complete ages) |

**Decision: Age MUST be included as a covariate in all differential expression models.**

### 10.2 Sex

| Dataset | Test | p-value | Significant | Interpretation |
|---|---|---|---|---|
| GSE63061 | Chi-square | 0.970957 | ❌ No | Sex well balanced across diagnosis groups |
| GSE63060 | Chi-square | 0.015495 | ✅ Yes | Sex imbalance — more females in AD group |
| GSE125583 | Chi-square | 0.803187 | ❌ No | Sex balanced across AD/Control |

**Decision: Sex should be included as a covariate in GSE63060 DE models. Optional for GSE63061 and GSE125583.**

### 10.3 Center / Batch

| Dataset | Test | p-value | Significant | Interpretation |
|---|---|---|---|---|
| GSE63061 | Chi-square | 0.001727 | ✅ Yes | 8 centers significantly imbalanced across diagnoses |
| GSE63060 | Chi-square | 0.000248 | ✅ Yes | 6 centers significantly imbalanced across diagnoses |
| GSE125583 | N/A | NaN | ❌ No | No batch metadata available |

**Decision: Center must be included as a blocking factor in GSE63061 and GSE63060 limma models.**

---

## 11. Batch-Effect Assessment (Task 8)

| Dataset | Batch Variable | n Batches | PC1 R² | PC2 R² | Batch Detected | Recommendation |
|---|---|---:|---:|---:|---|---|
| GSE63061 | Center code | 8 | 5.4% | 9.9% | ✅ Yes | Include center as covariate in limma |
| GSE63060 | Center code | 6 | 7.0% | 2.7% | ✅ Yes | Include center as covariate in limma |
| GSE140829 | Batch ID | 19 | 36.0% | 42.7% | ✅ Yes | Batch-associated PC structure — 19 batches |
| GSE125583 | N/A | 0 | — | — | ❌ No | Monitor SVA at DE stage |

> [!IMPORTANT]
> **GSE140829 has batch-associated PC scores** (36–43% PC variance explained by batch).
> This must be explicitly accounted for before any external validation.
> **No ComBat correction is applied at this stage** — that decision belongs at Stage 10 with
> proper cross-validation design.

**Strategy confirmed:**
- Blood cohorts (GSE63061, GSE63060): Include center as additive covariate in limma design formula.
  Do NOT apply global ComBat before CV (data leakage risk).
- GSE140829: Assess surrogate variable analysis (SVA) or design-matrix inclusion of batch.
- GSE125583: No batch metadata; apply SVA or monitor for latent technical effects in PyDESeq2.

---

## 12. QC Figures Generated

All figures saved to `outputs/figures/stage09/`:

| Figure | Description |
|---|---|
| `GSE63061_expression_boxplot.png` | Per-sample expression distribution (box plots, 388 samples). Medians tightly aligned. |
| `GSE63060_expression_boxplot.png` | Per-sample expression distribution (box plots, 329 samples). Medians aligned at ~7.5. |
| `GSE125583_library_sizes.png` | RNA-seq library size distribution. High variability (8M–93M), 41 extreme library sizes. |
| `feature_variance_distribution.png` | Feature (probe/gene) variance distributions per dataset. RNA-seq variance ~10⁸× microarray. |
| `GSE63061_PCA_diagnosis.png` | PCA by diagnosis label (AD/Control/MCI/Unknown). PC1=33.5%, PC2=11.0%. |
| `GSE63060_PCA_diagnosis.png` | PCA by diagnosis label. PC1=26.1%, PC2=13.1%. |
| `GSE125583_PCA_diagnosis.png` | PCA by diagnosis (AD/Control) in brain RNA-seq. PC1=46.4%, PC2=24.5%. |
| `sample_clustering.png` | Hierarchical clustering (Ward linkage, 1 − Pearson correlation). |
| `PCA_batch_effects.png` | Faceted PCA colored by center/batch for GSE63061, GSE63060, GSE140829. |
| `metadata_missingness.png` | Heatmap of missing clinical covariate proportions per dataset. |

---

## 13. Files Created / Modified

### New scripts

| File | Purpose |
|---|---|
| `scripts/03_quality_control.py` | Sample QC, feature QC, matrix orientation, flagging |
| `scripts/04_normalization_assessment.py` | Normalization audit, PCA, clustering |
| `scripts/05_batch_effect_analysis.py` | Confounder ANOVA/chi-square, batch ANOVA, figures |

### New tables (`outputs/tables/stage09/`)

| File | Rows | Description |
|---|---:|---|
| `sample_qc_summary.csv` | 1,593 | Sample-level QC metrics + PC1/PC2 coordinates |
| `sample_qc_flags.csv` | 202 | Provisional outlier flag records |
| `feature_qc_summary.csv` | 4 | Per-dataset feature-level statistics |
| `normalization_assessment.csv` | 4 | Normalization audit verdict per dataset |
| `confounder_assessment.csv` | 9 | Age/sex/batch confounder test results |
| `batch_effect_assessment.csv` | 4 | Batch ANOVA R² and significance per dataset |

### New tests

| File | Tests | Status |
|---|---:|---|
| `tests/test_stage09_qc.py` | 40 | ✅ All pass |
| `tests/test_stage09_normalization.py` | 24 | ✅ All pass |
| `tests/test_stage09_batch.py` | 22 | ✅ All pass |

### Derived interim files (preserving originals)

| File | Shape | Note |
|---|---|---|
| `data/interim/GSE63061/expression_samples_x_features.parquet` | (388, 32,049) | Derived transposed |
| `data/interim/GSE63060/expression_samples_x_features.parquet` | (329, 38,323) | Derived transposed |
| `data/interim/GSE140829/expression_samples_x_features.parquet` | (587, 47,008) | Derived transposed |
| `data/interim/GSE125583/raw_counts_samples_x_features.parquet` | (289, 39,376) | Derived transposed |

---

## 14. Test Results

```
Stage 09 test files:  86 tests → 86 passed, 0 failed
Full test suite:     310 tests → 310 passed, 0 failed

compileall scripts tests -q  →  Exit code 0 (no syntax errors)
pip check                    →  "No broken requirements found."
```

---

## 15. Scientific Limitations and Warnings

| Item | Risk | Mitigation |
|---|---|---|
| Blood microarray lncRNA coverage | Gene-symbol matching alone is **provisional**. ILMN_ probe IDs may not uniquely map to lncRNAs. | Explicit lncRNA annotation cross-walk required (Stage 10 Gate B). |
| Age confounding in all cohorts | Age differs by diagnosis in all three unlocked cohorts. Brain controls are older on average. | Age must be included in all DE model formulas. |
| Center confounding in blood cohorts | 8 centers (GSE63061) and 6 centers (GSE63060) significantly co-vary with diagnosis. | Center as additive covariate in limma design. |
| GSE140829 severe batch effect | 19 technical batches explain 36–43% of top PCA variance. | Plan batch covariate strategy before Gate E evaluation. |
| GSE125583 ≠ blood validation | Brain RNA-seq cannot serve as independent blood biomarker validation. | Used only as biological support, not validation. |
| 178 flagged samples retained | Some flagged samples may be low-quality. | Sensitivity analysis with/without flagged samples at DE stage. |
| GSE125583 identifier namespace | NCBI Entrez integer IDs; GENCODE uses Ensembl IDs. | Cross-walk audit required (Stage 10, Gate B). |

---

## 16. Recommended Preprocessing for Stage 10

| Task | Decision |
|---|---|
| **Microarray normalization** | None. Preserve pre-normalized lumi values. |
| **RNA-seq normalization** | None. Pass raw integer counts to PyDESeq2. |
| **DE model for GSE63061/GSE63060** | limma-voom with formula: `~diagnosis + age + sex + center` |
| **DE model for GSE125583** | PyDESeq2 with formula: `~diagnosis + age + sex` (Braak score optional) |
| **Flagged sample handling** | Retain in primary analysis; mark as sensitivity-check covariate |
| **Feature filtering** | At DE stage: remove always-zero genes (GSE125583: 1,181 genes) |
| **Probe → gene annotation** | Required at Stage 10 (GPL10558/GPL6947 → gene symbol) |
| **lncRNA annotation** | Gate B: probe-to-biotype cross-walk using GENCODE GTF |
| **GSE140829** | Remains locked. No diagnosis access until Gate E (final validation). |

---

## 17. Gate Assessment

| Gate | Status | Evidence |
|---|---|---|
| **Gate A** (Dataset Audit) | ✅ PASS | Real sample/feature counts verified; label dictionaries documented |
| **Gate B** (lncRNA Annotation) | ⏳ NOT TESTED | Deferred to Stage 10 — requires probe-to-biotype cross-walk |
| **Gate C** (QC/DE) | ✅ CONDITIONAL PASS | QC complete; confounders identified; DE models not yet run |
| **Gate D** (ML) | ⏳ NOT TESTED | Deferred |
| **Gate E** (External Validation) | 🔒 LOCKED | GSE140829 outcome labels not accessed |

---

## 18. Stage 09 Verdict

**CONDITIONAL PASS**

All four datasets are scientifically usable for downstream lncRNA differential expression analysis.
No data fabrication, no samples deleted, no batch correction applied, and no holdout labels accessed.

The following conditions must be satisfied before Stage 10 begins:
1. All DE models must include **age** as a covariate.
2. Blood cohort DE models must include **center** as a covariate.
3. Raw integer counts for GSE125583 must be fed to **PyDESeq2** (never to limma or DESeq2 as
   if they were microarray intensities).
4. GSE140829 diagnosis labels must remain sealed until Gate E.
5. lncRNA annotation cross-walk (Gate B) must be completed before any lncRNA-only ML.

---

*Report generated by Stage 09 scripts. All results are from real patient data. 
No values were fabricated. Seed = 42 used throughout.*
