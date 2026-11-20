# Stage 08 Data Audit Report
## lncSignal-AD — GEO Dataset Loading, Metadata Extraction and Data Audit

**Stage:** 08  
**Date:** 2026-11-21  
**Python:** 3.12.10 (.venv)  
**Status:** ✅ PASS — All 4 datasets loaded; 224/224 tests passed

> [!IMPORTANT]
> GSE140829 is the **locked external validation hold-out**. Its outcome labels were NOT extracted, stored, or summarised anywhere in this stage. Only structural/technical metadata was saved. This lock remains until Gate E (after model freeze and feature/probe compatibility confirmation).

---

## 1. Dataset Summary

| Dataset | Role | Tissue | Platform | Samples | Features | Scale | Feature ID | Missing | Diagnosis Locked |
|---|---|---|---|---|---|---|---|---|---|
| GSE63061 | Discovery | Whole blood | GPL10558 | 388 | 32,049 probes | Normalised microarray | `ILMN_*` | 0 (0%) | No |
| GSE63060 | Replication | Whole blood | GPL6947 | 329 | 38,323 probes | Normalised microarray | `ILMN_*` | 0 (0%) | No |
| GSE140829 | Ext. validation | Whole blood | GPL15988 | 587 | 47,008 probes | Normalised microarray | nuID (non-ILMN_) | 0 (0%) | **YES** |
| GSE125583 | Brain support | Fusiform gyrus | GPL16791 | 289 | 39,376 genes | Raw integer counts | Entrez GeneID | 0 (0%) | No |

---

## 2. Diagnosis Label Distribution (Unlocked Datasets Only)

| Dataset | AD | Control | MCI | Unknown | Total |
|---|---|---|---|---|---|
| GSE63061 | 139 | 134 | 109 | 6 | 388 |
| GSE63060 | 145 | 104 | 80 | 0 | 329 |
| GSE140829 | LOCKED | LOCKED | LOCKED | LOCKED | 587 |
| GSE125583 | 219 | 70 | 0 | 0 | 289 |

> [!NOTE]
> GSE63061: 6 "Unknown" samples — these were present in the data with an unrecognised label (possibly "borderline" or empty). They are preserved, never merged into Control.
>
> GSE125583: Diagnosis labels were parsed from sample titles (AD ##-## / CON ##-##) since the `disease state` characteristic field was not found in the metadata. This is documented in the parsing report.

---

## 3. Probe / Gene Identifier Audit

| Dataset | Platform | Probe ID format | Example probe IDs | Notes |
|---|---|---|---|---|
| GSE63061 | GPL10558 | `ILMN_*` (Illumina nuID-based) | `ILMN_1343048`, `ILMN_1651209` | 100% ILMN_ prefix confirmed |
| GSE63060 | GPL6947 | `ILMN_*` (Illumina HT-12 v3) | `ILMN_1343048`, `ILMN_1651209` | 100% ILMN_ prefix confirmed |
| GSE140829 | GPL15988 | nuID alphanumeric (non-ILMN_) | `0..gjSdbqnr.rUSTBI`, `0.a.oQZKQBR330TnQI` | 0 ILMN_ probes — requires gene-symbol/sequence crosswalk to GPL10558 |
| GSE125583 | GPL16791 | NCBI Entrez integer | `100287102`, `2778`, … | Index name: `GeneID`; range: 0 – 9,773,565 counts |

> [!WARNING]
> **GPL15988 uses a non-ILMN_ probe ID namespace.** Direct probe ID merging with GPL10558 (GSE63061) is impossible without a crosswalk. Gate B (probe annotation) must resolve this before any external validation analysis.

---

## 4. Expression Value Validation

| Dataset | dtype | Min | Max | NaN count | Integer? | Microarray safe? |
|---|---|---|---|---|---|---|
| GSE63061 | float32 | >0 | ~15 (log-scale) | 0 | No | ✅ Yes |
| GSE63060 | float32 | >0 | ~15 (log-scale) | 0 | No | ✅ Yes |
| GSE140829 | float32 | >0 | ~15 (log-scale) | 0 | No | ✅ Yes (but locked) |
| GSE125583 | int32 | 0 | 9,773,565 | 0 | ✅ Yes | N/A — RNA-seq counts; use PyDESeq2 |

> [!CAUTION]
> Never pass GSE63061/GSE63060/GSE140829 expression values into PyDESeq2. They are pre-normalised microarray intensities (float), not raw integer counts. Only GSE125583 is suitable for PyDESeq2.

---

## 5. Metadata Quality Summary

| Dataset | Duplicate IDs | Diagnosis labels valid | Age available | Sex available | Batch available |
|---|---|---|---|---|---|
| GSE63061 | PASS (0 dupes) | PASS | ✅ | ✅ | ❌ (not in this dataset) |
| GSE63060 | PASS (0 dupes) | PASS | ✅ | ✅ | ❌ |
| GSE140829 | PASS (0 dupes) | Locked | ✅ | ✅ | ✅ |
| GSE125583 | PASS (0 dupes) | PASS | ❌ (age at death not in series matrix) | ❌ | ❌ |

> [!NOTE]
> **Stage 10 correction.** The original age and sex availability entries above and the historical interpretation below were superseded by the [sample-level preflight audit](stage10_preflight_audit.md). The series matrix contains both characteristics.

> GSE125583 metadata age and sex fields were absent from the series matrix characteristics. The raw counts file contains gene-level counts only. Additional metadata (age at death, sex, RIN) may be available in the supplementary materials on GEO, but were not downloaded in this stage.

---

## 6. Scientific Integrity Checks

| Check | Result |
|---|---|
| GSE140829 diagnosis column absent from technical_metadata.csv | ✅ CONFIRMED |
| GSE140829 parsing_report.json has `diagnosis_locked: true` | ✅ CONFIRMED |
| MCI never merged into Control (all datasets) | ✅ CONFIRMED |
| Blood and brain data saved to separate directories | ✅ CONFIRMED |
| RNA-seq counts non-negative integers | ✅ CONFIRMED (min=0, max=9,773,565, dtype=int32) |
| Microarray expression float dtype | ✅ CONFIRMED (float32) |
| Expression matrix orientation: features × samples | ✅ CONFIRMED |
| Sample ID alignment: metadata ↔ expression columns | ✅ PASS (all 4 datasets) |
| No probe ID modification | ✅ CONFIRMED (ILMN_ prefixes intact, nuID preserved) |

---

## 7. Interim Files Generated

```
data/interim/
├── GSE63061/
│   ├── expression.parquet       # (32049, 388) float32 — features × samples
│   ├── metadata.csv             # 388 rows; columns: sample_id, dataset_id, diagnosis, age, sex, tissue, ethnicity, batch, platform_id
│   └── parsing_report.json
├── GSE63060/
│   ├── expression.parquet       # (38323, 329) float32 — features × samples
│   ├── metadata.csv             # 329 rows
│   └── parsing_report.json
├── GSE140829/
│   ├── expression.parquet       # (47008, 587) float32 — features × samples
│   ├── technical_metadata.csv   # 587 rows; NO diagnosis column
│   └── parsing_report.json      # diagnosis_locked: true
└── GSE125583/
    ├── raw_counts.parquet        # (39376, 289) int32 — genes × samples
    ├── metadata.csv             # 289 rows
    └── parsing_report.json
```

**Output tables:**
- [`outputs/tables/dataset_summary.csv`](../outputs/tables/dataset_summary.csv)
- [`outputs/tables/sample_metadata_summary.csv`](../outputs/tables/sample_metadata_summary.csv)
- [`outputs/tables/missing_values_summary.csv`](../outputs/tables/missing_values_summary.csv)

**Output figures:**
- `outputs/figures/sample_counts_by_diagnosis.png`
- `outputs/figures/missing_metadata_pct.png`
- `outputs/figures/expression_distributions.png`
- `outputs/figures/dataset_dimension_comparison.png`

---

## 8. Scripts Created

| Script | Purpose |
|---|---|
| `scripts/01_load_geo_datasets.py` | Parse all 4 GEO datasets; save interim parquet + JSON reports |
| `scripts/02_extract_metadata.py` | Standardise metadata; produce audit tables and plots |

**Tests:**
| Test file | Tests | Scope |
|---|---|---|
| `tests/test_geo_loading.py` | 22 tests | Sentinel detection, expression shape/dtype, probe ID preservation, metadata extraction, MCI preservation |
| `tests/test_metadata_extraction.py` | 46 tests | Diagnosis label mapping (all variants), schema validation, locked dataset handling |
| `tests/test_stage08_integrity.py` | 64 tests | Locked diagnosis, count integrity, float dtype, no duplicates, blood/brain separation, matrix orientation, parsing report structure |

**Total after Stage 08: 224/224 tests passed (0 failures, 0 skipped)**

---

## 9. Platform Identifier Namespace Issues (Gate B Prerequisite)

| Issue | Impact | Required before |
|---|---|---|
| GPL15988 uses nuID-style alphanumeric probe IDs (not `ILMN_*`) | Cannot directly crosswalk GSE140829 probes to GPL10558 probes by ID | Gate B (lncRNA annotation) and Gate E (external validation) |
| GPL10558 and GPL6947 use different probe sets despite both being `ILMN_*` | Probes cannot be directly compared across the two blood cohorts without mapping | Stage 09 (QC/normalisation) |
| GSE125583 uses NCBI Entrez IDs; GENCODE GTF uses Ensembl IDs | Joining requires crosswalk via `Human.GRCh38.p13.annot.tsv.gz` | Stage 09 (brain annotation) |

---

## 10. Gate C Status

| Gate | Status | Evidence |
|---|---|---|
| Gate A — Dataset audit | ✅ PASS | Real sample and feature counts verified; label dictionary confirmed |
| Gate B — lncRNA annotation | 🔵 NOT TESTED | Requires probe-to-biotype mapping (Stage 09+) |
| Gate C — QC/DE | 🔵 NOT TESTED | Requires normalisation and DE analysis (Stage 10+) |
| Gate D — ML | 🔵 NOT TESTED | Requires feature selection and CV (Stage 11+) |
| Gate E — External validation | 🔵 NOT TESTED | GSE140829 labels remain locked |

---

## 11. What You Can Show Your Professor at This Stage

1. **Real sample counts per dataset**: GSE63061 (388), GSE63060 (329), GSE140829 (587), GSE125583 (289).
2. **Clearly separated diagnosis labels**: AD/Control/MCI/Unknown — MCI never merged into Control.
3. **Expression matrix shapes and dtypes**: microarray as float32, RNA-seq as int32.
4. **Missing value audit**: 0 missing expression values across all datasets.
5. **Scientific integrity**: GSE140829 diagnosis locked; blood and brain data strictly separated.
6. **Platform identifier namespace warning**: GPL15988 probe IDs require crosswalk before Gate E.
7. **Statement**: *"We have loaded all four GEO datasets from local files, verified expression matrix integrity, extracted standardised metadata, and confirmed that the external validation hold-out (GSE140829) remains locked. No normalisation, differential expression, or machine learning has been performed."*

---

## 12. Next Stage

**Stage 09** — Quality control, normalisation, and confounders.
- Inspect expression value distributions per dataset.
- Apply appropriate normalisation (quantile/RMA for microarray; VST/rlog for RNA-seq).
- Identify and characterise batch effects.
- Remove low-quality samples using expression QC metrics.
- Document covariates for subsequent regression/correction.

> [!IMPORTANT]
> Do NOT begin Stage 09 until explicitly instructed. GSE140829 outcome labels must remain locked until Gate E.
