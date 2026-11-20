# lncSignal-AD

An exploratory study of long non-coding RNAs (lncRNAs) in Alzheimer’s disease. We ask whether blood lncRNA measurements differ between people with Alzheimer’s disease and healthy controls, whether those differences recur in a second blood cohort, and which expression patterns appear in postmortem brain tissue.

## What we found

- **Blood.** Three discovery probes met the 5% false discovery rate threshold. Two lacked an equivalent probe in the replication cohort. The third, measuring MIR29B2CHG, showed higher expression in both cohorts but missed the prespecified replication threshold (directional FDR 0.109). The study has yet to establish a replicated blood biomarker. [Differential expression report](docs/stage11_differential_expression_report.md)
- **Brain.** In fusiform gyrus, 1,067 of 2,157 valid lncRNA gene tests met the 5% FDR threshold after adjustment for age and sex. These results describe postmortem tissue. [Brain and blood comparisons](outputs/tables/stage11/brain_blood_comparison.csv)
- **Co-expression.** In the same brain samples, ZFAS1 expression tracked ribosomal protein genes, including RPL34 (adjusted correlation 0.921). Its strongly co-expressed partners were enriched for cytoplasmic translation. Cell-type composition could contribute to these patterns. [Co-expression report](docs/stage12_coexpression_report.md) · [Coding-gene partners](outputs/tables/stage12/edge_summary.csv)

## Data and methods

| Cohort | Measurement | Role and samples |
|---|---|---|
| [GSE63061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63061) | Whole blood, Illumina HT-12 v4 | Discovery, 139 AD and 134 controls in the comparison |
| [GSE63060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63060) | Whole blood, Illumina HT-12 v3 | Related AddNeuroMed replication cohort, 145 AD and 104 controls |
| [GSE140829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140829) | Whole blood, Illumina HT-12 v4 | Reserved external cohort, 587 samples. Outcome labels remain sealed |
| [GSE125583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125583) | Postmortem fusiform gyrus, RNA-seq | 289 samples. The adjusted analyses used 218 AD and 70 controls with complete covariates |

Blood probes were screened against GENCODE v50 transcript sequences. We retained individual probe IDs, modeled processed intensities with covariate-adjusted linear regression, and compared identical probe sequences across blood cohorts. Brain analyses used Entrez-indexed raw counts and PyDESeq2. Co-expression used partial correlations adjusted for diagnosis, age and sex, followed by enrichment against locally mapped Gene Ontology Biological Process terms. [Annotation methods](docs/lncrna_reference_methodology.md) · [Statistical analysis plan](docs/stage11_statistical_analysis_plan.md)

## Getting started

To explore the findings, start with the [reports](docs/), [result tables](outputs/tables/) or [figures](outputs/figures/). To work on the code, clone the repository and run these commands in PowerShell:

```powershell
git clone https://github.com/stackSentinel-32/lncRNA-Alzheimers-Biomarkers.git
cd lncRNA-Alzheimers-Biomarkers
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pytest -q tests/test_validation.py
```

The quick test checks validation logic and dataset configuration using the clone alone. Reproducing the research also requires the local inputs:

1. Obtain the GEO series matrices, array platform tables, brain count matrix, NCBI annotation and GENCODE v50 references listed in the [dataset configuration](config/datasets.yaml). Keep originals under the configured paths. The transcript FASTA belongs in `data/optional_raw_data/`.
2. Follow the [dataset audit](docs/stage08_data_audit.md), [quality control](docs/stage09_qc_report.md), [annotation](docs/stage10_annotation_report.md), [differential expression](docs/stage11_differential_expression_report.md) and [co-expression](docs/stage12_coexpression_report.md) reports in order. They record inputs, commands and expected outputs. Probe sequence verification requires a C++11-compatible `g++` on PATH.
3. After generating the local inputs and outputs, run the full checks:

```powershell
.\.venv\Scripts\python.exe -m pytest -v
.\.venv\Scripts\python.exe -m compileall scripts tests -q
.\.venv\Scripts\python.exe -m pip check
```

## Contributing

Create a branch for one change, add or update tests for the behavior you changed, and open a pull request. In the PR, describe the research question or bug, the command you ran, and any effect on the results. If you use a fork, open the PR against this repository’s `main` branch.

Keep raw downloads, sample-level working files, virtual environments and large generated matrices local. The repository shares source code, configuration, tests, reports, figures and compact research results. Some report links point to artifacts generated during a full local run, including the 1.74 GB pairwise correlation table. Follow the [research conventions](AGENTS.md) when working with annotations, cohort labels or the external holdout.

Computational Genomics and Proteomics (OE3E40).
