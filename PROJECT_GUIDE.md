# Establishing the Role of lncRNAs as Biomarkers in Alzheimer's Disease
## Complete end-to-end computational genomics, transcriptomics, statistics, and machine-learning project guide

**Course:** Computational Genomics and Proteomics (OE3E40)  
**Institute:** PDPM Indian Institute of Information Technology, Design and Manufacturing, Jabalpur  
**Project type:** Academic research project + reproducible data-science portfolio  
**Prepared:** 21 November 2026  
**Status:** Research plan and executable starter framework; **not** a report of completed experiments.


> **Scientific integrity statement.** Nothing in this guide should be presented as a new experimental finding until it has been run and checked. Previously discussed exploratory figures and counts are **provisional**, not final validated biomarker results. Public datasets can be reprocessed, but this work **cannot establish clinical diagnostic effectiveness** without appropriate external and experimental validation. **No synthetic AUC, p-value, or gene hit is reported.**

---

## Contents

1. [The problem in simple English and Hinglish](#1-the-problem-in-simple-english-and-hinglish)
2. [Precise problem statement, scope, and hypotheses](#2-precise-problem-statement-scope-and-hypotheses)
3. [Research design and decision gates](#3-research-design-and-decision-gates)
4. [Public datasets: verified facts and links](#4-public-datasets-verified-facts-and-links)
5. [Exactly what to download and why](#5-exactly-what-to-download-and-why)
6. [Tech stack and installation](#6-tech-stack-and-installation)
7. [Repository layout and data dictionary](#7-repository-layout-and-data-dictionary)
8. [Step 1: Acquire and audit GEO datasets](#8-step-1-acquire-and-audit-geo-datasets)
9. [Step 2: Identify credible lncRNA measurements](#9-step-2-identify-credible-lncrna-measurements)
10. [Step 3: Quality control, normalization, and confounders](#10-step-3-quality-control-normalization-and-confounders)
11. [Step 4: Differential expression analysis](#11-step-4-differential-expression-analysis)
12. [Step 5: Co-expression and pathway analysis](#12-step-5-co-expression-and-pathway-analysis)
13. [Step 6: Machine-learning biomarker panels](#13-step-6-machine-learning-biomarker-panels)
14. [Step 7: External validation, error analysis, explainability](#14-step-7-external-validation-error-analysis-explainability)
15. [Step 8: RNA-seq brain support and optional R analysis](#15-step-8-rna-seq-brain-support-and-optional-r-analysis)
16. [Step 9: Optional MATLAB Bioinformatics Toolbox component](#16-step-9-optional-matlab-bioinformatics-toolbox-component)
17. [What to do if the first approach fails](#17-what-to-do-if-the-first-approach-fails)
18. [Plots, report, presentation, and research claims](#18-plots-report-presentation-and-research-claims)
19. [GitHub, demo application, and portfolio](#19-github-demo-application-and-portfolio)
20. [Suggested four-week timeline and completion checklist](#20-suggested-four-week-timeline-and-completion-checklist)
21. [Viva questions with answers](#21-viva-questions-with-answers)
22. [Two-minute pitch](#22-two-minute-pitch)
23. [References, dataset citations, and verification log](#23-references-dataset-citations-and-verification-log)
24. [Step-by-step Codex and Antigravity prompts](#24-step-by-step-codex-and-antigravity-prompts)

> **Implementation note (uploaded-files version):** Your local `raw/` folder already appears to contain `GSE63061_series_matrix.txt.gz`, `GSE63060_series_matrix.txt.gz`, `GSE140829_series_matrix.txt.gz`, `GSE125583_series_matrix.txt.gz`, `GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz`, three GPL tables, GENCODE v50 lncRNA GTF and transcript FASTA, an NCBI annotation TSV, and `GPL15988_files/`. **Do not download them again blindly.** Their exact content and formatting must still be inspected. For the GSE125583 brain pipeline, **prefer checking your existing `GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz` first**; the older guide’s `GSE125583_RAW.tar` parser is an alternate source and should not be treated as mandatory. An additional comprehensive GENCODE annotation GTF may be necessary for rigorous coding-vs-lncRNA off-target auditing.

> **Agent workflow:** Use the project-level `AGENTS.md` provided with this guide; open the folder containing your data in **Codex** or **Google Antigravity**; copy **one stage prompt at a time** from [Section 24](#24-step-by-step-codex-and-antigravity-prompts). The code already in Sections 8–16 is **starter code**. Agents should inspect local file schemas, adapt it, test it, and stop when a scientific gate fails—not force it to run by fabricating results.

---

## 1. The problem in simple English and Hinglish

### 1.1 What is Alzheimer's disease (AD)?

**English:** Alzheimer's is a disease in which the brain progressively loses some of its ability to remember, think, and function. Two well-studied pathological changes are amyloid-beta accumulation and abnormal tau protein. Inflammation and loss of connections between nerve cells also matter. Diagnosis and research already use established biomarkers, including amyloid/tau-related imaging, cerebrospinal-fluid tests, and increasingly validated blood-based protein tests. Our proposed lncRNA panel is **experimental**; it is **not** a replacement for those established tests.

**Hinglish:** Alzheimer's mein brain ke neurons aur unke connections affect hote hain. Memory aur thinking dheere-dheere kharab ho sakti hai. Brain mein amyloid aur tau se related changes hote hain. Hum ek **naya research signal** dhoondh rahe hain, koi approved medical test replace nahi kar rahe.

### 1.2 DNA, RNA, mRNA, lncRNA

```text
DNA (instruction book)
       |
       | transcription
       v
RNA molecules
   |                         |
   v                         v
mRNA: can help make      lncRNA: usually does not
proteins                  encode a protein, but can
                          regulate gene activity
```

- **DNA:** genetic information.
- **mRNA:** carries protein-making instructions.
- **lncRNA:** a class of mostly non-protein-coding RNA transcripts, conventionally longer than 200 nucleotides. Some influence transcription, RNA stability, chromatin, and cellular signaling. **Not every lncRNA has a proven function.**
- **Gene expression:** a measurement of how much RNA associated with a gene is present in a biological sample.

**Hinglish:** DNA instruction book hai. mRNA usse protein banane ke instructions le jaata hai. lncRNA generally protein nahi banata, lekin genes ke regulation mein role ho sakta hai. Hum measure karenge ki kis patient mein kaunsa lncRNA kitna express ho raha hai.

### 1.3 What is a biomarker?

A **biomarker** is a measurable biological feature associated with a state or process. A useful candidate should be measurable reliably, discriminate between groups in independent data, and add information beyond obvious confounders. For a clinical diagnostic biomarker, further analytical and prospective clinical validation is necessary.

**Hinglish:** Biomarker body ka measurable indicator hota hai. Agar AD aur control patients mein kisi RNA ka level reliably different ho, toh wo **candidate biomarker** ho sakta hai. Sirf ek plot mein difference milna final proof nahi hai.

### 1.4 Why lncRNAs?

The attached 2024 review by **Canoy et al.** consolidated **362 unique lncRNAs** reported in the AD literature, alongside 741 miRNAs. It discusses links to amyloid-beta, tau, apoptosis, and inflammation. This supports **biological motivation**, **not** a claim that 362 lncRNAs are clinically validated biomarkers. Literature-associated examples include **BACE1-AS, NEAT1, and XIST**; their specific direction and mechanism can depend on tissue, assay, and study context. [R1]

### 1.5 One-sentence project explanation

**English:** We will examine public gene-expression data to identify reliably measured lncRNAs associated with Alzheimer's disease, investigate their relationships with other genes, and test whether a small lncRNA panel can distinguish AD from controls in independent samples.

**Hinglish:** Public patient gene-expression data lekar hum check karenge ki kaunse lncRNAs Alzheimer's se associated hain, woh kin genes se connected hain, aur kya unka chhota panel ML se AD aur healthy controls mein difference bata sakta hai.

> **Implementation prompt for this section:** [Agent Stage 01 — Codex and Antigravity](#agent-stage-01)  


---

## 2. Precise problem statement, scope, and hypotheses

### Problem statement

**Can reliably annotated long non-coding RNA measurements show reproducible association with Alzheimer's disease and provide incremental out-of-sample discrimination beyond demographic covariates?**

### Objectives

1. Inspect AD and control sample metadata, expression formats, and platform coverage.
2. Define a reproducible lncRNA annotation process with probe-mapping evidence.
3. Perform quality control and account for available confounders.
4. Identify differentially expressed candidate lncRNAs using valid statistics and false-discovery-rate (FDR) correction.
5. Study **associations** between candidate lncRNAs and coding genes; explore pathway enrichment.
6. Train and compare interpretable biomarker models with leakage-resistant cross-validation.
7. Test a frozen candidate panel on a separate cohort **only if comparable lncRNA measurements exist**.
8. Publish reproducible code, honest results, limitations, and an explanatory dashboard if time permits.

### Testable hypotheses

- **H1:** At least some *reliably measured* lncRNAs are associated with AD after adjustment for available covariates and multiple testing.
- **H2:** A compact lncRNA expression model discriminates AD from controls better than a simple demographic baseline on out-of-sample data.
- **H3:** Some candidate lncRNAs correlate with protein-coding genes in biologically relevant pathways; these correlations **suggest**, but do not prove, biological function.

Null results are valid outcomes. Set performance thresholds (for example, exploratory external AUC target 0.75) **in advance**, label them as *targets*, and report the measured value even if it falls below the target.

### What is in scope / out of scope?

| In scope | Not demonstrated by this project |
|---|---|
| Reanalysis of public human expression profiles | A clinically deployable Alzheimer's diagnostic test |
| Associations, hypothesis generation | Direct lncRNA-mediated causation |
| Internal cross-validation and possible external validation | Proof of *early* diagnosis unless genuinely early-stage longitudinal samples exist |
| Statistical and ML baselines | Equivalence of brain and blood expression |
| Literature-backed biological interpretation | qRT-PCR or independent lab validation (unless separately conducted) |

**Key distinction:** AD vs control classification in already diagnosed subjects does **not** by itself establish a test for pre-symptomatic or early disease.

> **Implementation prompt for this section:** [Agent Stage 02 — Codex and Antigravity](#agent-stage-02)  


---

## 3. Research design and decision gates

### Recommended design: two complementary analyses, not an invalid tissue merge

```text
                           PUBLIC GEO DATA
                                 |
              +------------------+------------------+
              |                                     |
              v                                     v
    BLOOD microarrays                          BRAIN RNA-seq
    GSE63061 (discovery)                       GSE125583 (support)
    GSE63060 (related batch)                    gene counts + metadata
    GSE140829 (external cohort)                 |
              |                                v
    Audit platform probes                 Verify lncRNA gene-ID coverage
              |                                |
        +-----+------+                         v
        |            |                    DE + co-expression
     Adequate      Inadequate                  |
     specific      specific                    v
     lncRNA        lncRNA                 Biological context
     mapping       mapping                (NOT blood validation)
        |            |
        v            v
    Blood DE,     Reframe claim;
    panel ML      use better RNA-seq,
        |         or brain-only study
        v
   Nested CV -> frozen panel -> GSE140829 test (only if comparable)
```

### Critical decision gates

| Gate | Question | Pass action | Fail action |
|---|---|---|---|
| G0 | Are the GEO datasets and actual file types verified? | Record exact sources and hashes | Inspect GEO metadata and supplementary files again |
| G1 | Can probes be **uniquely** attributed to current lncRNA genes? | Continue lncRNA-specific blood analysis | Do not label generic microarray genes as lncRNAs |
| G2 | Are AD/control labels clean and sufficiently numerous? | Analyze those groups; keep MCI separate | Correct metadata; exclude ambiguous samples |
| G3 | Do both discovery and holdout cohorts measure the **same** validated features? | Pre-register shared-feature panel | No lncRNA-specific external-validation claim |
| G4 | Are covariates and batch effects understood? | Run adjusted models and sensitivity checks | State residual confounding; prefer same-cohort analyses |
| G5 | Does an independent result replicate? | Report validated candidate panel | Report negative result and investigate shift without relabeling it a success |

**Pre-declared feasibility criterion (a project decision, not a universal biological rule):** Audit whether there are at least ~20 uniquely mapped, reasonably variable lncRNA genes shared across the intended blood cohorts. If not, prioritise a scientifically honest brain-first study or a more appropriate blood RNA-seq cohort. The exact threshold is not a statistical law and should not be tuned to obtain a better AUC.

### Four principles that must never be violated

- **No leakage:** all label-informed feature selection and hyperparameter tuning must happen **inside** training folds; the external cohort is a final, once-only test.
- **No platform shortcut:** an Illumina probe's older gene symbol or a GENCODE symbol match alone does **not** establish unique lncRNA measurement; inspect the probe sequence and possible cross-hybridization.
- **No tissue shortcut:** findings from postmortem brain tissue cannot be called a validated blood diagnostic panel.
- **No fake clinical claim:** expression-associated markers are *candidate biomarkers* until appropriately validated.

> **Implementation prompt for this section:** [Agent Stage 03 — Codex and Antigravity](#agent-stage-03)  


---

## 4. Public datasets: verified facts and links

Counts below represent **published AD/control subsets** or the GEO series size as specified. They are not necessarily the count after your QC and label exclusions. Verify the actual metadata when you run the scripts.

| GEO accession | Tissue / type | Platform | AD | Controls | Other / total | Planned role | Main caution |
|---|---|---|---:|---:|---|---|---|
| [GSE63061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63061) | Peripheral whole blood; processed **microarray**, not RNA-seq | [GPL10558](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558), Illumina HumanHT-12 v4 | 139 | 134 | 109 MCI plus six borderline/ambiguous records; **388 total GEO samples** | Discovery **after probe-mapping audit** | Older general expression chip; only some lncRNAs detectable |
| [GSE63060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63060) | Whole blood; microarray | [GPL6947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6947), HumanHT-12 v3 | 145 | 104 | MCI and possibly other records | Replication / cohort batch sensitivity | Same AddNeuroMed overall project; different chip version; **not fully independent clinical validation** |
| [GSE140829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140829) | Whole blood; microarray | [GPL15988](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15988), HumanHT-12 v4 **nuID** | 204 | 249 | 134 MCI; **587 total** | Proposed separate-cohort final test | Shared *chip family* is not enough: nuID/probe-ID mapping and normalization must be checked |
| [GSE125583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125583) | Postmortem fusiform gyrus **brain**; RNA-seq | Illumina HiSeq 2500 | 219 | 70 | **289 AD + control** in published analysis | Independent *biological* support in another tissue | Available gene counts use **Entrez IDs**; lncRNA annotation coverage needs auditing |
| [GSE161199](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161199) | Peripheral blood mononuclear cells (PBMC); RNA-seq | Check exact GEO metadata | **Verify** | **Verify** | Also contains ALS and Parkinson's disease | **Contingency dataset** focused on coding and lncRNA expression | Smaller/mixed-condition cohort; do not assume adequate AD sample size |

**Important data corrections:**

1. **GSE63061:** 388 is the total series size; 139 AD + 134 control = **273 binary-classification cases** before further exclusions. Do **not** label its MCI or ambiguous records as healthy controls. The earlier exploratory “136 controls” count was an inconsistent label-handling result, not a valid final figure.
2. **GSE63060 + GSE63061** are batches within an AddNeuroMed context. A model replicated across them is useful but not the same as cross-study validation.
3. **GSE140829:** verify probe equivalence and cohort/sample preprocessing rather than assuming v4 and v4 nuID are immediately interchangeable.
4. **GSE125583:** published counts are Entrez-coded, so blindly merging with GENCODE's Ensembl gene IDs can yield **zero matches for a technical reason**. Conversely, successful matching does not guarantee that old quantification detected all lncRNAs.
5. The dataset counts above come from GEO-indexed records and published reuse articles, cited in [Section 23](#23-references-dataset-citations-and-verification-log). Recheck each series for sample exclusions, paired/duplicate participants, and changes in metadata.

### How to choose the final study design

- **For a clean short college demonstration:** begin with **GSE63061** for dataset loading, class counts, PCA, and differential-expression *exploration*, explicitly stating that lncRNA assignment is still under validation.
- **For a strong lncRNA scientific claim:** **first prove annotation adequacy**. If the blood arrays fail, shift the primary lncRNA investigation to a cohort with demonstrable RNA-seq lncRNA coverage. Prefer scientifically correct scope over a good-looking but incorrect AUC.
- **For a respectable ML portfolio:** use simple interpretable models, a strict unseen cohort if compatible, a README with caveats, and reproducible charts.

> **Implementation prompt for this section:** [Agent Stage 04 — Codex and Antigravity](#agent-stage-04)  


---

## 5. Exactly what to download and why

### The dataset files (only these first)

| Download | Link | Which item on the page? | Why |
|---|---|---|---|
| Blood discovery | [GSE63061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63061) | **Series Matrix File(s)** (or download via GEOparse) | Expression values + GEO sample metadata |
| Related blood batch | [GSE63060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63060) | **Series Matrix File(s)** | Replication/robustness, **not** independent proof |
| External blood | [GSE140829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140829) | **Series Matrix File(s)** | Final test after annotation/platform compatibility gate |
| Brain support | [GSE125583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125583) | **Supplementary** `GSE125583_RAW.tar` (per-sample TSV, ~87.4 MB compressed, verify size on GEO) | Processed per-sample gene-count/statistics files; **not** raw sequencing FASTQ |
| Contingency PBMC RNA-seq | [GSE161199](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161199) | Check supplementary counts and metadata **before downloading** | Alternative with explicit lncRNA profiling intent |

### Annotation files

| Download | Link | Why |
|---|---|---|
| Blood chip annotation 1 | [GPL10558](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558) | Map GSE63061 probe IDs to vendor gene IDs / symbols / sequences |
| Blood chip annotation 2 | [GPL6947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6947) | Map GSE63060 probes |
| Blood chip annotation 3 | [GPL15988](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15988) | Interpret GSE140829's **nuID** identifier scheme |
| Human gene annotation | [GENCODE human release page](https://www.gencodegenes.org/human/) | Download the **comprehensive gene annotation GTF** and/or long-noncoding-RNA annotation GTF; record release + genome build |
| Probe sequence resource | In platform files above or published manufacturer's annotation | Verify whether probes uniquely target the proposed lncRNA rather than also binding coding transcripts |

**Do not download now:** massive raw SRA FASTQ libraries, whole-genome FASTA, gigabytes of unrelated data, or paid software. **Do not download the whole transcriptome FASTA unless sequence-level remapping is actually needed** and you have chosen a reliable alignment/mapping procedure.

**Very important:** GENCODE gene types may have evolved since the original arrays were designed. Record your exact annotation release; do not silently treat modern gene symbols as equivalent to old probe annotations.

> **Implementation prompt for this section:** [Agent Stage 05 — Codex and Antigravity](#agent-stage-05)  


---

## 6. Tech stack and installation

### Recommended stack

| Layer | Technology | Reason |
|---|---|---|
| Main language | **Python 3.11** (recommended, consistently supported) | One language for data analysis, ML, plots, and app |
| Interactive work | Jupyter Notebook / VS Code / Google Colab | Inspect metadata, tables, plots |
| GEO import | `GEOparse` | Download and parse GEO series and platforms |
| Core processing | `pandas`, `numpy`, `scipy` | Tables, arrays, calculations |
| Statistics | `statsmodels` | Gene-wise linear models, FDR adjustment |
| Plots | `matplotlib`, `seaborn` | QC, PCA, boxplots, heatmaps, ROC |
| ML | `scikit-learn` | Imputation, scaling, LASSO, SVM, RF, nested CV |
| Pathway exploration | `gseapy` | Online GO/Reactome/Enrichr exploratory analysis; network access required |
| RNA-seq DE (optional) | R + `DESeq2` **or** compatible `pydeseq2` release | Correct raw-count negative-binomial modeling |
| R optional for microarrays | `limma` | Preferred microarray moderated differential-expression inference |
| MATLAB for course | MATLAB Bioinformatics Toolbox | GEO exploration and genomic interpretation; some modeling uses Statistics Toolbox |
| Version control | Git + GitHub | Transparent reproducible portfolio |
| Dashboard (optional) | `streamlit` | Research-only lncRNA explorer |

### Terminal commands (Windows PowerShell)

```powershell
mkdir ad-lncrna-biomarkers
cd ad-lncrna-biomarkers
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install GEOparse pandas numpy scipy matplotlib seaborn statsmodels scikit-learn gseapy jupyter joblib streamlit
python -m pip freeze > requirements-lock.txt
```

If PowerShell blocks activation, launch the virtual environment's Python directly, or adjust your user-level execution policy after confirming its security implications. In Colab, omit the venv commands and run the `pip install` line prefixed with `!`.

**Optional R**, only if you need preferred microarray or RNA-seq differential-expression statistics:

```r
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "DESeq2", "edgeR"))
```

**Avoid immediately installing every bioinformatics package:** some are large and version-sensitive. Install `rtracklayer`, `sva`, `clusterProfiler`, `pydeseq2`, SHAP, or WGCNA **only if that specific analysis is justified and feasible**.

> **Implementation prompt for this section:** [Agent Stage 06 — Codex and Antigravity](#agent-stage-06)  


---

## 7. Repository layout and data dictionary

### Create these folders

```text
ad-lncrna-biomarkers/
├── README.md
├── PROJECT_GUIDE.md            # this document
├── requirements-lock.txt
├── .gitignore
├── data/
│   ├── raw/                    # GEO original downloads (usually Git-ignored)
│   ├── annotations/            # GPL + GENCODE + mapping audit
│   ├── interim/                # parsed GEO matrices / metadata
│   └── processed/              # curated matrices (generated, reproducible)
├── scripts/
│   ├── 01_fetch_geo.py
│   ├── 02_propose_annotation.py
│   └── 03_analyze_blood.py
├── notebooks/                   # optional experiments and explanations
├── outputs/
│   ├── figures/
│   ├── tables/
│   └── models/
├── docs/
│   ├── methods.md
│   ├── dataset_audit.md
│   ├── probe_mapping_audit.md
│   ├── interpretation.md
│   └── references.md
└── app.py                       # optional research explorer
```

### Data dictionary (verify field names for *each* dataset)

| Logical field | Example of acceptable format | Why needed | Important warning |
|---|---|---|---|
| `sample_id` | `GSM...` | Match expression columns to people/samples | Detect duplicated technical/biological samples |
| `diagnosis` | `AD`, `Control`, `MCI`, `Unknown` | Outcome label | **Never** convert MCI/Unknown into controls |
| `age` | numeric years | Confounder | Missingness must be documented |
| `sex` | female/male/other/unknown as actually reported | Confounder | XIST can reflect sex rather than AD |
| `brain_region` | fusiform gyrus | Tissue-specific inference | Never merge different brain regions indiscriminately |
| `tissue` | whole blood, PBMC, brain | Interpret sample matrix | Different tissues have different cell composition |
| `rin` | RNA integrity number, if reported | Technical RNA-quality covariate | Missing in many GEO records |
| `batch` | site, processing run, platform | Batch correction / stratification | Batch fully confounded with disease cannot be repaired mathematically |
| `probe_id` | `ILMN_...` or nuID string | Identify array feature | Probe not necessarily a unique gene |
| `gene_symbol` | `NEAT1`, `XIST` | Human-readable biology | May be outdated/ambiguous |
| `gene_id` | Ensembl `ENSG...` **or** Entrez numeric ID | Stable identifier with namespace | Never merge Entrez and Ensembl IDs as if identical |
| `gene_type` | `lncRNA`, `protein_coding`, etc. | Separate RNA classes | Use a recorded annotation release |
| `expression` | numeric normalized intensity **or** raw count | Input to analysis | Methods depend fundamentally on measurement type |

**Matrix orientation:** the retrieval stage saves **features × samples**. ML stage transposes to **samples × genes/features**. Every merge must be verified by `sample_id` and feature ID, never by position alone.

### Minimum audit log

Create `docs/dataset_audit.md` and record series ID, retrieval date, sample count by exact original label, number excluded by rule, platform, annotation release, duplicates, covariate availability, expression scale, and file SHA256 checksums. This log is as valuable as a pretty ROC curve.

> **Implementation prompt for this section:** [Agent Stage 07 — Codex and Antigravity](#agent-stage-07)  


---

## 8. Step 1: Acquire and audit GEO datasets

**Goal:** Save sample-level metadata and processed expression matrices before you do any statistical testing.

**Important:** The following program uses a real, well-established Python GEO parser, but it has **not been executed against live GEO downloads in this document-generation environment**. GEO download availability and returned field names may vary, so the script intentionally prints audit outputs and refuses to guess missing labels.

Save as `scripts/01_fetch_geo.py`:

```python
"""Download a GEO series, preserve original metadata, and flag labels safely.
Usage: python scripts/01_fetch_geo.py GSE63061
"""
from pathlib import Path
import re
import sys
import json
import GEOparse
import pandas as pd

GSE = sys.argv[1] if len(sys.argv) > 1 else "GSE63061"
raw = Path("data/raw")
out = Path("data/interim") / GSE
raw.mkdir(parents=True, exist_ok=True)
out.mkdir(parents=True, exist_ok=True)

series = GEOparse.get_GEO(geo=GSE, destdir=str(raw), silent=False)

# Preserves original normalized/preprocessed values from GEO.
expression = series.pivot_samples("VALUE")  # rows: probes; columns: GSMs
expression = expression.apply(pd.to_numeric, errors="coerce")
expression.index.name = "feature_id"
expression.to_csv(out / "expression.csv")


def label_from_metadata(md):
    """Strict labels: refuse 'CTL to AD', borderline, unknown, mixed labels."""
    characteristics = md.get("characteristics_ch1", [])
    values = []
    for record in characteristics:
        if isinstance(record, str) and ":" in record:
            key, value = record.split(":", 1)
            if key.strip().lower() in {"status", "diagnosis", "disease state"}:
                values.append(value.strip())
    if len(values) != 1:
        return "Unknown"
    s = re.sub(r"\s+", " ", values[0]).strip().lower()
    if s in {"ad", "alzheimer", "alzheimer's disease", "alzheimers disease",
             "alzheimer disease"}:
        return "AD"
    if s in {"control", "ctl", "hc", "healthy control", "normal", "healthy"}:
        return "Control"
    if s in {"mci", "mild cognitive impairment"}:
        return "MCI"
    return "Unknown"


def get_characteristic(md, keys):
    for record in md.get("characteristics_ch1", []):
        if isinstance(record, str) and ":" in record:
            key, value = record.split(":", 1)
            if key.strip().lower() in keys:
                return value.strip()
    return None

rows = []
for gsm_id, gsm in series.gsms.items():
    md = gsm.metadata
    rows.append({
        "sample_id": gsm_id,
        "diagnosis": label_from_metadata(md),
        "original_characteristics": " | ".join(md.get("characteristics_ch1", [])),
        "age": get_characteristic(md, {"age", "age (years)", "age at death"}),
        "sex": get_characteristic(md, {"sex", "gender"}),
        "batch": get_characteristic(md, {"batch", "study batch", "site"}),
        "title": " | ".join(md.get("title", [])),
        "platform_id": " | ".join(md.get("platform_id", [])),
    })
metadata = pd.DataFrame(rows)
metadata.to_csv(out / "metadata.csv", index=False)

# Keep raw sample metadata JSON for manual inspection and re-labeling.
with open(out / "raw_sample_metadata.json", "w", encoding="utf-8") as f:
    json.dump({sid: gsm.metadata for sid, gsm in series.gsms.items()},
              f, indent=2, default=str)

for gpl_id, gpl in series.gpls.items():
    gpl.table.to_csv(out / f"{gpl_id}_platform.csv", index=False)

print("Series:", GSE)
print("Expression shape (features x samples):", expression.shape)
print("All diagnosis labels:\n", metadata["diagnosis"].value_counts(dropna=False))
print("Platform IDs:\n", metadata["platform_id"].value_counts(dropna=False))
print("Missing expression values:", int(expression.isna().sum().sum()))
print("Saved to:", out)
```

Run:

```bash
python scripts/01_fetch_geo.py GSE63061
python scripts/01_fetch_geo.py GSE63060
python scripts/01_fetch_geo.py GSE140829
```

**Do not proceed if the printed class counts disagree with the labels you can inspect on GEO.** For GSE140829 the disease field may use a different key or abbreviation, so first check `raw_sample_metadata.json`; update the parser with **documented exact values**, not arbitrary substring matches. Preserve unknown cases instead of silently reclassifying them.

### First sanity-check commands

```python
from pathlib import Path
import pandas as pd

p = Path("data/interim/GSE63061")
m = pd.read_csv(p / "metadata.csv")
x = pd.read_csv(p / "expression.csv", index_col=0)
print("Samples:", len(m), "Feature matrix:", x.shape)
print(m["diagnosis"].value_counts(dropna=False))
assert x.columns.is_unique
assert m["sample_id"].is_unique
assert set(x.columns).issubset(set(m["sample_id"]))
print(m.head(5).to_string(index=False))
```

### What you can legitimately show your professor at this stage

1. Number of downloaded GEO samples and clearly separated AD/control/MCI/ambiguous labels.
2. A small table showing accession, tissue, platform, gene/probe columns, and available metadata.
3. A screenshot of the matrix shape and missing-value check.
4. A statement: *"We have acquired the expression datasets and are auditing lncRNA-specific measurement coverage before making biomarker claims."*

> **Implementation prompt for this section:** [Agent Stage 08 — Codex and Antigravity](#agent-stage-08)  


---

## 9. Step 2: Identify credible lncRNA measurements

### This is the highest-risk part of the entire project

Most GEO blood datasets above are **general-purpose gene-expression microarrays**. A probe's historical symbol may refer to a protein-coding transcript, an ambiguous locus, or an RNA feature that has since been reannotated. Therefore:

```text
GEO microarray feature ID
      |
      v
Retrieve vendor GPL annotation + probe sequence
      |
      v
Map probe to a current gene / transcript annotation
      |
      v
Check gene biotype AND uniqueness / cross-hybridization
      |
  +---+----+
  |        |
  v        v
Reliable  Ambiguous
lncRNA    or coding hit
  |        |
  v        v
Keep      Exclude from
          lncRNA-specific claim
```

**Do not equate** a simple `Gene Symbol` match to *uniquely validated* probe specificity. That match is only a **candidate** until its provenance and sequence are reviewed. The previously mentioned count of ~152 possible blood lncRNAs was a preliminary exploratory annotation result, **not** 152 confirmed lncRNA measurements.

### Use GENCODE gene types carefully

Download the **comprehensive human GTF** from [GENCODE](https://www.gencodegenes.org/human/) and save it as something like `data/annotations/gencode.annotation.gtf.gz`. Note the release and genome build in `docs/probe_mapping_audit.md`.

Gene-type categories such as `lncRNA`, `lincRNA`, `antisense`, or `sense_intronic` occur across annotation releases. Avoid auto-adding all `processed_transcript` genes to the lncRNA set without checking their definition. If you download only GENCODE's lncRNA-only GTF, you will still need a comprehensive annotation to detect alternative protein-coding targets.

Save as `scripts/02_propose_annotation.py`:

```python
"""Propose lncRNA annotations by gene symbol.
IMPORTANT: this is a CANDIDATE MAP, not proof of probe specificity.
Usage: python scripts/02_propose_annotation.py data/annotations/gencode.annotation.gtf.gz GSE63061 GPL10558
"""
from pathlib import Path
import gzip
import re
import sys
import pandas as pd

gtf_path = Path(sys.argv[1])
gse_id = sys.argv[2]
gpl_id = sys.argv[3]
platform_path = Path("data/interim") / gse_id / f"{gpl_id}_platform.csv"
output_dir = Path("data/annotations")
output_dir.mkdir(parents=True, exist_ok=True)

lnc_biotypes = {
    "lncRNA", "lincRNA", "antisense", "sense_intronic", "sense_overlapping",
    "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "macro_lncRNA",
}


def attr(s, key):
    match = re.search(r'\b' + re.escape(key) + r' "([^"]+)"', str(s))
    return match.group(1) if match else None

# Read in chunks: GTF can be large. Keep only the 'gene' feature.
genes = []
for chunk in pd.read_csv(
    gtf_path, sep="\t", comment="#", header=None, compression="infer",
    usecols=[2, 8], names=["feature", "attributes"], chunksize=100000,
    dtype=str, low_memory=False
):
    rows = chunk[chunk["feature"] == "gene"]
    for attrs in rows["attributes"]:
        gene_id = attr(attrs, "gene_id")
        symbol = attr(attrs, "gene_name")
        biotype = attr(attrs, "gene_type") or attr(attrs, "gene_biotype")
        if gene_id and symbol and biotype:
            genes.append((gene_id.split(".")[0], symbol, biotype))

annotation = pd.DataFrame(genes, columns=["gene_id", "gene_symbol", "gene_type"])
annotation["is_lnc_candidate"] = annotation["gene_type"].isin(lnc_biotypes)

# Remove ambiguous symbols mapped to more than one current locus.
unique_symbols = annotation.groupby("gene_symbol")["gene_id"].nunique()
annotation = annotation[annotation["gene_symbol"].map(unique_symbols) == 1]

platform = pd.read_csv(platform_path, dtype=str)
print("PLATFORM COLUMNS:", platform.columns.tolist())
probe_col = next((c for c in ["ID", "ID_REF", "ProbeID"] if c in platform.columns), None)
symbol_col = next((c for c in ["Gene Symbol", "GeneSymbol", "SYMBOL", "Symbol",
                               "gene_symbol", "ILMN_Gene"] if c in platform.columns), None)
if not probe_col or not symbol_col:
    raise ValueError("Inspect the platform columns and choose the correct ID/symbol columns manually.")

p = platform[[probe_col, symbol_col]].copy()
p.columns = ["feature_id", "gene_symbol"]
p["gene_symbol"] = p["gene_symbol"].astype(str).str.strip()
# Exclude multiple/unclear assignments rather than choosing one arbitrarily.
p = p[p["gene_symbol"].str.fullmatch(r"[A-Za-z0-9_.-]+", na=False)]
proposal = p.merge(annotation, on="gene_symbol", how="inner")
proposal = proposal[proposal["is_lnc_candidate"]].copy()
proposal["mapping_status"] = "PROVISIONAL_SYMBOL_ONLY"
proposal["mapping_evidence"] = "Gene-symbol join; sequence specificity NOT verified"
proposal = proposal.drop_duplicates(subset=["feature_id", "gene_id"])
path = output_dir / f"{gse_id}_{gpl_id}_candidate_lncrna_probes.csv"
proposal.to_csv(path, index=False)
print("Possible probes (NOT validated):", len(proposal))
print("Possible unique lncRNA genes (NOT validated):", proposal["gene_id"].nunique())
print("Saved:", path)
```

**Required manual verification BEFORE ML:**

1. Inspect the platform table's actual `gene_symbol` column. `ILMN_Gene` may be a **vendor ID**, not a modern symbol; if so, the automatic join is invalid and must be replaced with proper cross-references.
2. Inspect the actual probe sequence and mapping from GPL10558/GPL15988 (nuIDs may need decoding or a published crosswalk).
3. Check potential matches to **protein-coding as well as lncRNA transcripts**. Exclude probes that map ambiguously or cross-hybridize. A real sequence alignment method (e.g., BLAT/Bowtie) with transparent parameters is preferable to naive string matching; probe mismatches and strand orientation matter.
4. Record stable gene/transcript IDs, transcriptome annotation release, genome build, match specificity, and any deprecated symbols.
5. Produce `data/processed/GSE63061_validated_probe_map.csv` **only for probes that pass the audit**. Do this separately for GSE140829, then join by a proven stable gene identifier (and document probe equivalence).
6. **Holdout label access:** annotation may use public platform sequences; do not use external AD/control outcomes to pick features.

A validated mapping CSV must look like:

```csv
feature_id,gene_id,gene_symbol,mapping_evidence
<actual_GPL_probe_ID>,<verified_gene_ID>,<verified_symbol>,<sequence_and_annotation_evidence>
```

The angle-bracket examples above are **placeholders**, not real rows. Never invent mapping results. A probe table with no validated rows is a valid reason to change datasets.

**What to report:** number of starting probes; candidate lncRNA mappings; number excluded for missing/ambiguous symbols; number excluded for cross-hybridization; number of uniquely mapped genes; number overlapping with the external platform.

> **Implementation prompt for this section:** [Agent Stage 09 — Codex and Antigravity](#agent-stage-09)  


---

## 10. Step 3: Quality control, normalization, and confounders

### 10.1 QC questions

- Are there samples with unusually low/high overall expression or many missing values?
- Do boxplots show consistent distributions **for data already normalized by the source**?
- Does PCA reveal extreme outliers or platform/batch differences?
- Do the AD and control groups differ in age, sex, batch, RIN, brain region, or cell composition?
- Are biological/technical replicate samples being double-counted?

### 10.2 Normalization: different data need different methods

| Data | Correct approach | Incorrect shortcut |
|---|---|---|
| Already variance-stabilized + quantile-normalized microarray | Inspect the provided processing notes and distributions; generally **do not normalize or log-transform blindly again** | Feeding normalized intensities into DESeq2 as raw counts |
| Raw RNA-seq integer counts | DESeq2 size-factor normalization and negative-binomial model, or an appropriate count-based alternative | Running a normal t-test on unscaled raw counts |
| ML from processed array data | Impute missing features and standardize **inside each training fold** | Standardizing on all training + external test samples together |
| Multi-study arrays | Harmonize probe identity and independently documented preprocessing; assess study shift | Combining unrelated platforms and forcing batch correction with test labels |

For example, GSE63061's GEO sample processing describes **variance-stabilization and quantile normalization**. Do not assume its matrix needs an additional `log2(x+1)` step. The numeric values in sample records have already been transformed. [D1]

### 10.3 Confounders

A **confounder** is a factor related both to diagnosis and measured expression, making a gene appear disease-associated when it may be responding to something else.

- **Age:** AD subjects may be older.
- **Sex:** XIST expression often strongly reflects sex; it is an excellent **confounding example**, not automatic proof of AD effect.
- **RNA integrity (RIN):** degraded RNA can change measurements.
- **Batch/site:** laboratory differences can look like biological differences.
- **Blood-cell composition:** gene-expression change can reflect more/fewer neutrophils or lymphocytes, not changes within a cell.
- **Brain region and cell-type composition:** brain measurements depend on which cells and regions were sampled.

**Hinglish:** Agar AD group mein females zyada hain aur XIST females mein naturally high hota hai, model XIST ko “AD gene” samajh sakta hai. Isliye age/sex aur other confounders check karna zaroori hai.

### 10.4 When to batch-correct

Use `batch` as a covariate in differential expression when possible. For multi-cohort validation, keep the external cohort untouched by a correction trained jointly with its disease labels. A batch fully aligned with diagnosis (`all AD in batch A; all controls in batch B`) makes diagnosis and batch effects unidentifiable. Do not claim ComBat solves that.

### 10.5 Plots you must make

1. Sample counts and missing values chart.
2. Sample expression boxplot (small representative subset if hundreds of samples).
3. PCA plot colored by diagnosis.
4. PCA plot colored by batch/site, if available.
5. Age distribution and sex-by-diagnosis chart.
6. Clustering/outlier plot; document actual removal rules.

**Interpretation examples (not actual observations):**

- “The median intensity distributions are comparable, *consistent with*, but not conclusive proof of, successful preprocessing.”
- “PCA does not show complete AD/control separation. That is not surprising; disease effects can be subtle and confounders may be large.”
- “An extreme sample was flagged based on a pre-specified technical QC rule, not because it made the model score worse.”

> **Implementation prompt for this section:** [Agent Stage 10 — Codex and Antigravity](#agent-stage-10)  


---

## 11. Step 4: Differential expression analysis

### Simple explanation

**Question:** For a given validated lncRNA, is its expression systematically different between AD and controls, after accounting for available covariates?

- **Difference in expression:** how much group levels differ.
- **p-value:** strength of statistical evidence against a null hypothesis, conditional on model assumptions; not the probability that the hypothesis is true.
- **FDR / adjusted p-value:** controls the expected fraction of false discoveries in the selected list under the procedure's assumptions.
- **Effect size:** how large the difference is (often more informative than significance alone).

**Hinglish:** Har valid lncRNA ka AD vs control expression compare karo. Lekin hazaar tests karne par kuch randomly significant milenge. Isliye **Benjamini-Hochberg FDR correction** use karte hain.

### Use the right inference method

- **Normalized blood microarray:** the preferred established method is **limma** in R, because of its empirical-Bayes moderated variance estimates. A Python `statsmodels` per-feature linear model is a valid *transparent baseline* but **not exactly equivalent to limma**.
- **Brain RNA-seq integer counts:** **DESeq2** (or a compatible `pydeseq2` implementation) is appropriate. Do **not** use DESeq2 on array intensities or normalized TPM/RPKM.
- **Fold-change label:** for a VST-like processed array, a linear-model coefficient is an **expression-scale difference**. Call it `log2FoldChange` **only if** you have verified the actual scale is log2 intensity and that interpretation is justified.

### Python statistical example (for already prepared lncRNA matrix)

The integrated script in Section 13 produces this analysis automatically, after you supply validated feature mappings. Conceptually:

```python
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

# expr: samples x validated lncRNA genes; labels: index-aligned Series (0/1)
# Optional covariates should be included only after correct encoding and audit.

def differential_expression(expr: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    labels = labels.loc[expr.index]
    design = sm.add_constant(pd.DataFrame({"AD": labels.astype(float)}, index=expr.index))
    results = []
    for gene in expr.columns:
        y = pd.to_numeric(expr[gene], errors="coerce")
        valid = y.notna() & np.isfinite(y)
        if valid.sum() < 10 or y[valid].nunique() < 3:
            continue
        fit = sm.OLS(y[valid], design.loc[valid]).fit()
        results.append({
            "gene": gene,
            "AD_effect_on_supplied_scale": fit.params["AD"],
            "p_value": fit.pvalues["AD"],
            "n": int(valid.sum()),
        })
    table = pd.DataFrame(results)
    if table.empty:
        raise ValueError("No genes passed the minimum quality requirements")
    table["fdr_bh"] = multipletests(table["p_value"], method="fdr_bh")[1]
    return table.sort_values(["fdr_bh", "p_value"])
```

The example above is **unadjusted** for age, sex, and batch; use the adjusted implementation or `limma` when reliable covariate data exist. A useful initial candidate criterion may be `FDR < 0.05` **plus meaningful effect size** and consistency in replication. It is not valid to call a gene “confirmed” solely from its FDR.

### Report outputs

- `outputs/tables/differential_expression.csv`: feature/gene, group effect, p-value, BH FDR, sample size, covariates.
- Volcano plot: significance vs effect size; label y-axis and x-axis honestly.
- Heatmap of **selected exploratory** genes (for descriptive analysis; not a validation estimate).
- Top-hit table with literature references and mapping-confidence score.

**ML anti-leakage rule:** differential-expression screening across all training samples may be used for **exploratory biology**, but if it determines ML features, re-run that selection **inside each training fold**. Never select the “top 20” using the full discovery dataset and then advertise cross-validation as unbiased.

> **Implementation prompt for this section:** [Agent Stage 11 — Codex and Antigravity](#agent-stage-11)  


---

## 12. Step 5: Co-expression and pathway analysis

### Why do this?

Differential expression answers “what changes?” Co-expression asks “which genes tend to change together?” Enrichment asks “are co-expressed protein-coding genes over-represented in known biological pathways?”

**Hinglish:** Agar lncRNA `L` aur protein-coding genes `A`, `B`, `C` patients mein saath change hote hain, toh ho sakta hai related biology ho. But **correlation causation nahi hota**.

### Recommended compact analysis for laptop/Colab

1. Keep a **same-tissue** expression matrix of lncRNAs and mRNAs from the **same patients**.
2. Pre-filter low-information genes without consulting the outcomes where possible.
3. For a small list of biologically plausible candidate lncRNAs, compute Pearson or Spearman correlation against expressed mRNAs.
4. Prefer residual or partial correlations adjusting for available age, sex, batch, and major cell-composition variables.
5. Control FDR across tested lncRNA-mRNA pairs, not just over the chosen significant subset.
6. For each lncRNA, obtain a list of strongly co-expressed protein-coding genes. Feed these into a GO/Reactome enrichment tool **with a suitable background gene universe**.
7. Interpret terms as **hypotheses**, not directly proven targets or interactions.

### Minimal correlation example

```python
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import pandas as pd

# Same-sample, already QC'd, appropriately normalized expression matrices.
# Rows = samples; columns = gene IDs. Do not use different tissues/cohorts here.

def correlate_one_lncrna(lnc_values, mrna_matrix):
    rows = []
    for coding_gene in mrna_matrix.columns:
        pair = pd.concat([lnc_values, mrna_matrix[coding_gene]], axis=1).dropna()
        if len(pair) < 20 or pair.iloc[:, 0].nunique() < 3 or pair.iloc[:, 1].nunique() < 3:
            continue
        rho, p = spearmanr(pair.iloc[:, 0], pair.iloc[:, 1])
        rows.append((coding_gene, rho, p, len(pair)))
    res = pd.DataFrame(rows, columns=["coding_gene", "rho", "p_value", "n"])
    if res.empty:
        return res
    res["fdr_bh"] = multipletests(res.p_value, method="fdr_bh")[1]
    return res.sort_values("fdr_bh")
```

**Caution:** this snippet is an unadjusted *illustration*. In the formal analysis, residualize relevant confounders or use an appropriate partial-correlation model. For enrichment, first ensure IDs correspond to protein-coding genes, not arbitrary probe names; background universe should be all genes that were reasonably measurable/tested in that experiment.

### Pathway enrichment (optional online step)

```python
# Requires an internet connection and accepted gene-symbol mapping.
import gseapy as gp

# significant_mrna_symbols = [...]  # replace with actual validated results
# universe_symbols = [...]          # measured protein-coding background
# enr = gp.enrichr(gene_list=significant_mrna_symbols,
#                  gene_sets=["GO_Biological_Process_2023", "Reactome_2022"],
#                  organism="Human", outdir="outputs/enrichment")
```

Before running, verify which gene-set collection years are currently offered by the service. Include the assay-specific background where the tool permits. **Pathway names appearing in literature do not mean your analysis found those pathways.**

### WGCNA (optional extension)

WGCNA finds modules of co-expressed genes. It can be informative but adds many tuning decisions (sample clustering, soft power, module size, preservation, confounding). On limited time and hardware, **transparent lncRNA-mRNA correlation + appropriately adjusted enrichment is preferable** to a superficial WGCNA figure. Only use WGCNA if you can explain module preservation and have enough well-QC'd samples.

> **Implementation prompt for this section:** [Agent Stage 12 — Codex and Antigravity](#agent-stage-12)  


---

## 13. Step 6: Machine-learning biomarker panels

### 13.1 The ML question

**Can a small number of confidently measured lncRNAs predict AD vs Control on subjects whose labels were not used to train the model?**

We compare:

| Model | Idea in plain English | Why include it? | Risk |
|---|---|---|---|
| **Age + sex logistic baseline** | Use basic demographic information only | Shows whether RNA adds any information beyond simple confounders | Available ages/sex may be incomplete or uneven |
| **L1/LASSO logistic regression** | Penalize coefficients, naturally shrink some to zero | Main interpretable **sparse panel** | Correlated genes may lead to unstable selected features |
| **Random forest with in-fold selection** | Combine decisions from multiple trees | Checks non-linear structure | Overfits easily; importance can be misleading |
| **RBF SVM with in-fold selection** | Find non-linear separation between groups | Strong small/medium dataset benchmark | Harder to interpret and calibrate |

**Why not use deep learning as the default?** A few hundred subjects and limited validated lncRNA features are poor conditions for a complex network. Transparent, regularized models are easier to validate scientifically.

### 13.2 What do the hyperparameters mean?

- **LASSO `C`:** inverse of regularization strength. Smaller `C` = stronger penalty = usually fewer selected lncRNAs. Try a short grid, e.g. `[0.01, 0.1, 1, 10]`, not 1,000 values.
- **Random forest `n_estimators`:** number of trees.
- **`max_depth`:** maximum tree complexity.
- **`min_samples_leaf`:** minimum samples per leaf; larger values constrain overfitting.
- **SVM `C`:** trade-off between classification mistakes and margin simplicity.
- **SVM `gamma`:** influence radius of each training sample with the RBF kernel.

**Hinglish:** Hyperparameter model ka setting hai jo training se **pehle** choose/tune hota hai. Grid Search alag settings try karta hai. Nested CV mein andar wala loop tuning ke liye, aur bahar wala loop performance check karne ke liye hota hai. Bahar ke test fold ko feature selection ke time touch nahi karna.

### 13.3 Nested CV schematic

```text
FULL DISCOVERY COHORT (GSE63061)
          |
          v
  Outer 5-fold CV  <------ estimates out-of-sample performance
          |
          +---- Outer training fold (80%)
          |          |
          |          v
          |      Inner CV (3-fold)
          |      -> imputation
          |      -> scaling
          |      -> feature selection
          |      -> hyperparameter search
          |      -> fit best model
          |
          +---- Outer unseen fold (20%)
                     |
                     v
                   predict
          |
          v
    Aggregate all outer-fold predictions
          |
          v
    Choose model family, refit on *all discovery*
          |
          v
   FINAL EXTERNAL COHORT (GSE140829), ONCE
```

For formal performance reporting, the inner and outer splits must be **patient/group aware** if the same person has multiple samples; use `StratifiedGroupKFold` or group-preserving alternatives as appropriate. The simple code below assumes one independent sample per subject and requires you to audit that assumption.

### 13.4 Complete Python starter for QC, DE baseline, nested ML, and optional external test

**Inputs:**

- `data/interim/GSE63061/expression.csv` and `metadata.csv` from Step 1.
- `data/processed/GSE63061_validated_probe_map.csv` with **real, reviewed** mappings from Step 2.
- For external testing: matching GSE140829 files + validated mapping. You must **explicitly authorize** the holdout only after checking source preprocessing and probe equivalence.

Save as `scripts/03_analyze_blood.py`. This is a **general working script**, not a promise that a particular GEO file's metadata and probe identifiers will match without the documented audit. With no validated map, **it intentionally stops instead of manufacturing lncRNA results**.

```python
"""AD lncRNA discovery: QC, exploratory DE, nested CV, optional external holdout.

Prerequisites:
  python scripts/01_fetch_geo.py GSE63061
  Manually curate: data/processed/GSE63061_validated_probe_map.csv
Run:
  python scripts/03_analyze_blood.py

Data rules:
  - Input matrices have PROBES x SAMPLES, as exported by 01_fetch_geo.py.
  - A validated probe map must be supported by sequence/annotation review.
  - No synthetic results, no diagnosis imputation, no holdout-based selection.
"""
from pathlib import Path
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import (
    roc_auc_score, roc_curve, confusion_matrix, f1_score,
    recall_score, precision_score, accuracy_score
)
import joblib

SEED = 42
DISCOVERY = "GSE63061"
EXTERNAL = "GSE140829"
# Change to True ONLY after a documented probe/scale/comparability audit.
APPROVE_EXTERNAL_TEST = False
MIN_SHARED_FEATURES_FOR_EXTERNAL = 20  # project feasibility gate, not universal rule
OUTPUT = Path("outputs")
FIG = OUTPUT / "figures"
TAB = OUTPUT / "tables"
MODEL = OUTPUT / "models"
for folder in [FIG, TAB, MODEL]:
    folder.mkdir(parents=True, exist_ok=True)


def read_cohort(accession):
    base = Path("data/interim") / accession
    map_path = Path("data/processed") / f"{accession}_validated_probe_map.csv"
    if not map_path.exists():
        raise FileNotFoundError(
            f"{map_path} is missing. Review probe specificity; do NOT use "
            "the provisional symbol-only annotation as verified lncRNAs."
        )

    exp = pd.read_csv(base / "expression.csv", index_col=0)
    exp.index = exp.index.astype(str)
    exp.columns = exp.columns.astype(str)
    exp = exp.apply(pd.to_numeric, errors="coerce")
    meta = pd.read_csv(base / "metadata.csv", dtype={"sample_id": str})
    meta = meta.drop_duplicates(subset=["sample_id"]).set_index("sample_id")

    mp = pd.read_csv(map_path, dtype=str)
    need = {"feature_id", "gene_id", "gene_symbol", "mapping_evidence"}
    if not need.issubset(mp.columns):
        raise ValueError(f"Mapping file must include {need}")
    mp = mp.dropna(subset=list(need))
    mp = mp[mp["mapping_evidence"].str.len() > 5].copy()
    mp = mp[mp["feature_id"].isin(exp.index)]
    if mp.empty:
        raise ValueError("No validated probe IDs match the expression matrix")
    if mp.groupby("feature_id")["gene_id"].nunique().max() != 1:
        raise ValueError("A probe maps to multiple gene IDs; resolve ambiguity")
    if mp["mapping_evidence"].str.contains(
        "PROVISIONAL|NOT verified|PLACEHOLDER", case=False, regex=True
    ).any():
        raise ValueError("Provisional/placeholder mappings must not enter ML")

    # Pre-declared median aggregation of probes measuring the same gene.
    # Review isoform/probe consistency; different isoforms are a limitation.
    mp = mp.drop_duplicates(subset=["feature_id"])
    selected = exp.loc[mp["feature_id"]].copy()
    selected.index = mp.set_index("feature_id").loc[selected.index, "gene_id"].values
    expression_genes = selected.groupby(level=0).median().T

    # Parse labels strictly; no MCI/Unknown in binary analysis.
    meta = meta.loc[meta.index.intersection(expression_genes.index)].copy()
    meta = meta[meta["diagnosis"].isin(["AD", "Control"])].copy()
    expression_genes = expression_genes.loc[meta.index]
    labels = (meta["diagnosis"] == "AD").astype(int)
    if labels.nunique() < 2:
        raise ValueError("Only one class remains; audit metadata parsing")
    if len(meta) < 40:
        warnings.warn("Small binary cohort: nested-CV estimates may be unstable")
    print(accession, "samples:", len(meta), "gene features:", expression_genes.shape[1])
    print("Class counts:", meta["diagnosis"].value_counts().to_dict())
    return expression_genes, labels, meta, mp


def cleanup_features(X):
    # Unsupervised per-cohort eligibility; outer-fold-specific variation
    # screening is handled in the model pipeline, not on external labels.
    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.loc[:, X.notna().mean(axis=0) >= 0.80]
    X = X.loc[:, X.nunique(dropna=True) > 1]
    if X.shape[1] < 3:
        raise ValueError("Too few usable validated lncRNA features")
    return X


def qc_plots(X, y, accession):
    # QC/PCA are descriptive; this is not supervised feature selection.
    rng = np.random.default_rng(SEED)
    subset_ids = rng.choice(X.index.to_numpy(), size=min(32, len(X)), replace=False)
    fig, ax = plt.subplots(figsize=(13, 4))
    sns.boxplot(data=X.loc[subset_ids].T, ax=ax, color="lightgray", showfliers=False)
    ax.tick_params(axis="x", labelrotation=90, labelsize=6)
    ax.set_title(f"{accession}: representative sample distributions (source units)")
    ax.set_ylabel("Processed expression, original supplied scale")
    fig.tight_layout()
    fig.savefig(FIG / f"{accession}_boxplot.png", dpi=180)
    plt.close(fig)

    # PCA on most variable features is for *visual QC*, not for feature selection.
    variable = X.var().nlargest(min(300, X.shape[1])).index
    a = SimpleImputer(strategy="median").fit_transform(X[variable])
    a = StandardScaler().fit_transform(a)
    coordinates = PCA(n_components=2, random_state=SEED).fit_transform(a)
    chart = pd.DataFrame({"PC1": coordinates[:, 0], "PC2": coordinates[:, 1],
                          "diagnosis": y.map({0: "Control", 1: "AD"}).values})
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.scatterplot(data=chart, x="PC1", y="PC2", hue="diagnosis", ax=ax)
    ax.set_title(f"{accession}: exploratory PCA (not a classifier)")
    fig.tight_layout()
    fig.savefig(FIG / f"{accession}_pca.png", dpi=180)
    plt.close(fig)


def exploratory_de(X, y, meta):
    # This table is DESCRIPTIVE. It must NOT drive features of an
    # already-reported cross-validation unless run again INSIDE every fold.
    d = pd.DataFrame({"AD": y.astype(float)}, index=X.index)
    covariates = []
    age = pd.to_numeric(meta.get("age", pd.Series(index=X.index, dtype=float)),
                        errors="coerce")
    if age.notna().mean() >= 0.8 and age.nunique() > 3:
        d["age"] = age
        covariates.append("age")
    sex = meta.get("sex", pd.Series(index=X.index, dtype=str))
    sex = sex.astype(str).str.lower().str.strip().map(
        {"female": 1.0, "f": 1.0, "male": 0.0, "m": 0.0}
    )
    if sex.notna().mean() >= 0.8 and sex.nunique() == 2:
        d["female"] = sex
        covariates.append("female")
    d = d.dropna()
    if d["AD"].nunique() != 2:
        raise ValueError("Covariate complete-case subset loses a class")
    design = sm.add_constant(d, has_constant="add")
    rows = []
    for gene in X.columns:
        measurement = pd.to_numeric(X.loc[d.index, gene], errors="coerce")
        keep = measurement.notna() & np.isfinite(measurement)
        if keep.sum() < 20 or measurement[keep].nunique() < 3:
            continue
        if d.loc[keep, "AD"].nunique() != 2:
            continue
        try:
            fit = sm.OLS(measurement[keep], design.loc[keep]).fit()
            rows.append({
                "gene_id": gene,
                "effect_AD_on_source_scale": float(fit.params["AD"]),
                "p_value": float(fit.pvalues["AD"]),
                "n": int(keep.sum()),
                "covariates": ",".join(covariates) if covariates else "none",
            })
        except (ValueError, np.linalg.LinAlgError):
            continue
    result = pd.DataFrame(rows)
    if result.empty:
        print("No genes were suitable for this DE model")
        return result
    result["fdr_bh"] = multipletests(result["p_value"], method="fdr_bh")[1]
    result = result.sort_values("fdr_bh")
    result.to_csv(TAB / "discovery_exploratory_de.csv", index=False)
    return result


def make_models(n_features):
    common = [
        ("fill", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
    ]
    # The supervised univariate screening below is INSIDE the CV pipeline;
    # L1 logistic regression itself performs sparse feature selection.
    k_grid = sorted(set([min(5, n_features), min(10, n_features)]))
    return {
        "lasso_logistic": (
            Pipeline(common + [("model", LogisticRegression(
                penalty="l1", solver="liblinear",
                random_state=SEED, max_iter=4000
            ))]),
            {"model__C": [0.01, 0.1, 1.0, 10.0]}
        ),
        "random_forest": (
            Pipeline(common + [
                ("select", SelectKBest(score_func=f_classif, k=k_grid[0])),
                ("model", RandomForestClassifier(random_state=SEED, n_jobs=1))
            ]),
            {"select__k": k_grid,
             "model__n_estimators": [200, 400],
             "model__max_depth": [3, 6],
             "model__min_samples_leaf": [3, 7]}
        ),
        "rbf_svm": (
            Pipeline(common + [
                ("select", SelectKBest(score_func=f_classif, k=k_grid[0])),
                ("model", SVC(kernel="rbf", probability=True,
                              random_state=SEED))
            ]),
            {"select__k": k_grid,
             "model__C": [0.1, 1.0, 10.0],
             "model__gamma": ["scale", 0.01]}
        ),
    }


def evaluate(y_true, probability):
    decision = (np.asarray(probability) >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, decision, labels=[0, 1]).ravel()
    return {
        "auc": float(roc_auc_score(y_true, probability)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else None,
        "specificity": float(tn / (tn + fp)) if tn + fp else None,
        "precision": float(precision_score(y_true, decision, zero_division=0)),
        "f1": float(f1_score(y_true, decision, zero_division=0)),
        "accuracy": float(accuracy_score(y_true, decision)),
        "TN": int(tn), "FP": int(fp), "FN": int(fn), "TP": int(tp)
    }


def nested_cv(X, y, models):
    # A class must have enough subjects for the chosen CV design.
    if y.value_counts().min() < 10:
        raise ValueError("Insufficient subjects per class for this nested CV setup")
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    inner = StratifiedKFold(n_splits=3, shuffle=True, random_state=SEED + 1)
    summaries = []
    oof_all = {}
    for name, (pipeline, grid) in models.items():
        oof = pd.Series(index=X.index, dtype=float)
        best_per_fold = []
        for fold, (train_idx, valid_idx) in enumerate(outer.split(X, y), start=1):
            search = GridSearchCV(
                pipeline, grid, scoring="roc_auc", cv=inner,
                n_jobs=1, refit=True, error_score=np.nan
            )
            search.fit(X.iloc[train_idx], y.iloc[train_idx])
            prob = search.predict_proba(X.iloc[valid_idx])[:, 1]
            oof.iloc[valid_idx] = prob
            best_per_fold.append({"fold": fold, "settings": search.best_params_})
        if oof.isna().any():
            raise RuntimeError(f"Missing out-of-fold predictions for {name}")
        metrics = evaluate(y, oof)
        metrics["model"] = name
        summaries.append(metrics)
        oof_all[name] = oof
        pd.DataFrame({"sample_id": X.index, "AD": y, "oof_probability": oof}).to_csv(
            TAB / f"{name}_nested_cv_predictions.csv", index=False
        )
        (TAB / f"{name}_fold_tuning.json").write_text(
            json.dumps(best_per_fold, indent=2), encoding="utf-8"
        )
        print(name, "out-of-fold metrics:", metrics)

    summary = pd.DataFrame(summaries).sort_values("auc", ascending=False)
    summary.to_csv(TAB / "nested_cv_model_comparison.csv", index=False)
    fig, ax = plt.subplots(figsize=(7, 6))
    for name, probabilities in oof_all.items():
        fpr, tpr, _ = roc_curve(y, probabilities)
        ax.plot(fpr, tpr, label=f"{name} (OOF AUC {roc_auc_score(y, probabilities):.3f})")
    ax.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Chance")
    ax.set(xlabel="False positive rate", ylabel="True positive rate",
           title="Discovery cohort nested-CV ROC (OOF estimates)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "discovery_nested_cv_roc.png", dpi=180)
    plt.close(fig)
    return summary, oof_all


def stratified_auc_interval(labels, probabilities, repetitions=1000):
    # Exploratory patient-level non-parametric stratified bootstrap interval.
    rng = np.random.default_rng(SEED)
    labels = np.asarray(labels)
    probabilities = np.asarray(probabilities)
    ids0 = np.where(labels == 0)[0]
    ids1 = np.where(labels == 1)[0]
    vals = []
    for _ in range(repetitions):
        sample = np.r_[rng.choice(ids0, len(ids0), replace=True),
                       rng.choice(ids1, len(ids1), replace=True)]
        vals.append(roc_auc_score(labels[sample], probabilities[sample]))
    return np.quantile(vals, [0.025, 0.975]).tolist()


def run():
    X, y, meta, probe_map = read_cohort(DISCOVERY)
    X = cleanup_features(X)
    qc_plots(X, y, DISCOVERY)
    de = exploratory_de(X, y, meta)
    if not de.empty:
        print("Exploratory FDR<0.05 genes:", int((de.fdr_bh < 0.05).sum()))
    if X.shape[1] < 5:
        raise ValueError("Insufficient feature count for a meaningful panel study")

    # External study audit can identify available feature IDs WITHOUT using its labels.
    external_files = (
        (Path("data/interim") / EXTERNAL / "expression.csv").exists()
        and (Path("data/processed") / f"{EXTERNAL}_validated_probe_map.csv").exists()
    )
    X_test = y_test = meta_test = None
    if external_files:
        # Read metadata here to check count/coverage; do not tune on outcomes.
        X_test, y_test, meta_test, _ = read_cohort(EXTERNAL)
        common_features = sorted(set(X.columns) & set(X_test.columns))
        print("Validated shared lncRNA gene IDs:", len(common_features))
        if len(common_features) < MIN_SHARED_FEATURES_FOR_EXTERNAL:
            print("External gate FAIL: too few shared validated genes.")
            X_test = y_test = None
        elif (X_test[common_features].notna().mean(axis=0) < 0.80).any():
            print("External gate FAIL: insufficient measurable data in shared genes.")
            X_test = y_test = None
        else:
            X = X[common_features]
            X_test = X_test[common_features]
    else:
        print("External files/verified map unavailable: discovery analysis only.")

    models = make_models(X.shape[1])
    summary, oof = nested_cv(X, y, models)
    selected_name = summary.iloc[0]["model"]
    chosen_pipeline, chosen_grid = models[selected_name]
    final_search = GridSearchCV(
        chosen_pipeline, chosen_grid, scoring="roc_auc",
        cv=StratifiedKFold(5, shuffle=True, random_state=SEED + 2),
        n_jobs=1, refit=True
    )
    final_search.fit(X, y)
    joblib.dump(final_search.best_estimator_, MODEL / "discovery_fitted_model.joblib")
    (MODEL / "frozen_model_manifest.json").write_text(json.dumps({
        "discovery": DISCOVERY, "external": EXTERNAL,
        "model_family": selected_name,
        "parameters": final_search.best_params_,
        "feature_gene_ids_in_order": X.columns.tolist(),
        "random_seed": SEED,
        "external_approved": APPROVE_EXTERNAL_TEST,
        "note": "Do not interpret this file as clinical software"
    }, indent=2), encoding="utf-8")
    print("Chosen family by discovery CV:", selected_name)

    if selected_name == "lasso_logistic":
        weights = final_search.best_estimator_.named_steps["model"].coef_[0]
        selected = pd.DataFrame({"gene_id": X.columns, "coefficient": weights})
        selected = selected[selected.coefficient != 0]
        selected = selected.reindex(selected.coefficient.abs().sort_values(ascending=False).index)
        selected.to_csv(TAB / "lasso_selected_genes.csv", index=False)
        print("Non-zero LASSO features:", len(selected))
        # A panel with 0 or >10 features may require a tighter in-CV penalty grid;
        # do not retune after reading external outcomes.

    if X_test is not None:
        if not APPROVE_EXTERNAL_TEST:
            print("HOLDOUT LOCKED. Audit comparable probes and preprocessing first.")
        else:
            probabilities = final_search.best_estimator_.predict_proba(X_test)[:, 1]
            results = evaluate(y_test, probabilities)
            results["auc_bootstrap_95pct"] = stratified_auc_interval(y_test, probabilities)
            print("FINAL EXTERNAL RESULTS:", results)
            (TAB / "external_results.json").write_text(
                json.dumps(results, indent=2), encoding="utf-8")
            pd.DataFrame({
                "sample_id": X_test.index,
                "true_AD": y_test,
                "predicted_probability": probabilities,
                "predicted_AD_at_0_5": (probabilities >= 0.5).astype(int)
            }).to_csv(TAB / "external_predictions.csv", index=False)
            fpr, tpr, _ = roc_curve(y_test, probabilities)
            fig, ax = plt.subplots(figsize=(6, 5))
            ax.plot(fpr, tpr, label=f"External AUC={results['auc']:.3f}")
            ax.plot([0, 1], [0, 1], linestyle="--", color="gray")
            ax.set(xlabel="False positive rate", ylabel="True positive rate",
                   title="Locked independent-cohort evaluation")
            ax.legend()
            fig.tight_layout()
            fig.savefig(FIG / "external_roc.png", dpi=180)
            plt.close(fig)

    print("Saved figures, tables and trained research model in", OUTPUT)


if __name__ == "__main__":
    run()
```

**Practical cautions before running this analysis script:**

1. The program deliberately **refuses to run with no validated lncRNA mapping**. This is desirable: supplying fabricated probe rows would invalidate the science.
2. `random_forest` and `rbf_svm` use a supervised `SelectKBest` **inside each fold**; the principal sparse panel is LASSO logistic regression. If you alter selection rules or tuning grids, re-run the **entire nested CV** and record the reason. Never use external outcomes to select a more favorable configuration.
3. The feature-exclusion function is an exploratory eligibility filter. For strict internal CV, any sample-dependent variance/missingness filter whose outcome could affect scoring should also be moved into the pipeline (in-fitted training folds). **The final inference should use the strict in-fold version.** The code already handles imputation/scaling/label-informed selection in folds, but the cohort-wide eligibility screen is a simplification to document.
4. **The `read_cohort(EXTERNAL)` call sees metadata labels before external approval** because it audits class counts. For an especially strict locked-test protocol, audit external feature IDs without loading outcome labels until final evaluation. Keep the final holdout labels inaccessible to everyone tuning the method.
5. Choosing the top family by an outer-CV leaderboard is a practical model-development step; internal comparison can still be optimistic due to model-family selection. Only the untouched external test offers a stronger final estimate.
6. If normalization differs greatly across cohorts, the external score may reflect platform shift more than biology. **Never** recalibrate or reselect features on final-test labels and then reuse the same test score as confirmation.
7. For accurate “five-to-ten lncRNAs” results, enforce a panel-size constraint **within the inner CV** (or use a pre-specified family/grid that yields that size). Do not post-hoc choose genes based on the final holdout outcome.
8. If two samples belong to the same donor or repeated visit, replace the fold splitter with a group-aware one.
9. The output file `discovery_fitted_model.joblib` is a research artifact, **not** clinical software. Only load joblib files you created and trust.

### 13.5 Baseline comparison: age + sex

A biomarker panel is more credible if it improves over an age/sex-only model on the **same patients**. Do this only when metadata completeness supports a fair comparison. A minimal implementation:

```python
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

# Use only verified age/sex measurements and the SAME subject subset used
# for the lncRNA comparison. Cross-validate with the same outer folds.
# Sex encoding shown for datasets using the explicit male/female strings.
demo = pd.DataFrame(index=meta.index)
demo["age"] = pd.to_numeric(meta["age"], errors="coerce")
demo["female"] = meta["sex"].astype(str).str.lower().map(
    {"female": 1, "f": 1, "male": 0, "m": 0}
)
if demo.notna().mean().min() >= 0.8:
    baseline = make_pipeline(SimpleImputer(strategy="median"),
                             StandardScaler(),
                             LogisticRegression(max_iter=2000, random_state=42))
    # Train and test the baseline with the SAME nested/outer folds as RNA model.
else:
    print("Age/sex baseline not comparable: too much missing metadata")
```

**Do not compare** an age/sex score from one set of patients with an RNA score from a different patient subset. Compare RNA-only, demographics-only, and (optionally) combined models on identical splits and feature-selection procedures.

### 13.6 Optional stability selection (useful novelty)

Repeat the **inner-training-only** feature-selection step under many bootstrap or repeated-fold samples. Report the proportion of resamples selecting each lncRNA, not just one LASSO coefficient list. A stable gene chosen in 80% of resamples is usually easier to defend than one chosen only once—**but the 80% is an example, not a measured value**. Stability does not remove biological confounding or prove causality.

> **Implementation prompt for this section:** [Agent Stage 13 — Codex and Antigravity](#agent-stage-13)  


---

## 14. Step 7: External validation, error analysis, explainability

### 14.1 Why accuracy alone is insufficient

If most samples are controls, a model can appear accurate simply by predicting "control" often. Report:

| Metric | Formula / meaning | Why |
|---|---|---|
| Sensitivity / recall | `TP / (TP + FN)` | How many AD subjects were detected? |
| Specificity | `TN / (TN + FP)` | How often were controls correctly recognized? |
| Precision / PPV | `TP / (TP + FP)` | Among predicted AD, what proportion were AD in this sample? Prevalence-dependent |
| F1 | `2 × precision × recall / (precision + recall)` | Balances recall and precision |
| ROC-AUC | Probability a randomly chosen AD case gets a higher score than a control | Threshold-free discrimination summary |
| Confusion matrix | TN, FP, FN, TP | Reveals types of error |
| Calibration (optional) | Predicted probability vs observed frequency | Does 0.8 predicted risk behave like 80% event frequency in this cohort? |

### 14.2 Correct external validation protocol

1. Choose the candidate lncRNA **gene IDs and probe mappings** without examining external AD/control outcomes.
2. Align the same proven comparable features and document data preprocessing differences.
3. Tune preprocessing and hyperparameters **on discovery data only**.
4. Freeze model family, coefficients, selected panel, decision rule, and analysis code.
5. Unlock the external labels, run exactly one pre-specified evaluation, report every metric and a confidence interval.
6. Document **all** failures and subgroup disparities, even if the external AUC is disappointing.

**Why GSE140829 is not automatically valid:** GPL15988 uses **nuIDs**, while GSE63061 is associated with GPL10558. Even if both are Illumina HumanHT-12 v4 family, a correct gene/probe mapping must be proven. Identical *gene symbols* can conceal different probes. Independently normalized data may also be on different scales.

### 14.3 Error analysis

After final testing, inspect:

- False negatives: AD classified as control; were they early/atypical cases, or lower RNA expression for technical reasons?
- False positives: controls classified as AD; could other illnesses, age, or inflammation explain them?
- Calibration and distribution shift.
- AUC/sensitivity by age or sex **only if each subgroup is large enough**. Small subgroup scores are unstable; do not overinterpret them.
- Confidence intervals: use patient-level bootstrap rather than reporting AUC to four decimals without uncertainty.

### 14.4 Explainability

- **LASSO:** coefficients show direction and magnitude **after the model's scaling**, not causal effects. Correlated features make coefficients unstable.
- **Random forest:** permutation importance on held-out fold data is more directly interpretable than naive impurity importance; still correlational.
- **SHAP:** an optional visualization for the final candidate model if installed and compatible; it is a model explanation, not a molecular causal mechanism.

### 14.5 Independent validation vs technical replication

- **GSE63060** checks a related AddNeuroMed batch (useful replication, weaker independence).
- **GSE140829** is a separate cohort and stronger independent validation **if actual feature/assay compatibility is documented**.
- **GSE125583** is brain tissue, so it is **biological context**, not direct blood diagnostic external validation.
- **qRT-PCR** on new samples would be a logical laboratory validation step but is beyond this computational-only scope.

### 14.6 Clinically responsible conclusion

Do **not** say: “We developed a test that detects Alzheimer's early.”

You may say, **if supported by actual results**: “We identified a preliminary panel of correctly annotated lncRNA expression features associated with AD and evaluated its discrimination under nested cross-validation and an independent cohort. Prospective and qRT-PCR validation are still required.”

> **Implementation prompt for this section:** [Agent Stage 14 — Codex and Antigravity](#agent-stage-14)  


---

## 15. Step 8: RNA-seq brain support and optional R analysis

### 15.1 Why we need the brain cohort

The blood arrays may provide only partial lncRNA coverage. GSE125583 is RNA-seq from postmortem **fusiform gyrus** tissue, with a published AD/control analysis of **219 AD and 70 controls**. It provides valuable biological evidence *in brain tissue*, but it is **not** a patient-matched blood and brain study. [D4]

### 15.2 Critical GSE125583 file-format difference

GEO lists a supplementary archive `GSE125583_RAW.tar` (around 87.4 MB compressed) containing sample-level TSV files. This is **processed count/normalized-statistics data**, not SRA FASTQ reads. GEO sample processing says genes were counted and supplementary IDs are **Entrez Gene IDs**, with counts and nRPKM statistics. Therefore:

- Inspect the first file to determine **exact column names** before coding a parser.
- Select the actual **raw integer count** field for DESeq2, not nRPKM.
- Use an authoritative **Entrez-to-Ensembl** annotation crosswalk (gene IDs have different namespaces).
- Audit what fraction of modern lncRNAs are present before assuming full lncRNA coverage.
- Check all 289 samples for clinical/technical covariates and excluded samples.

**Local-file-first correction:** Your downloaded `GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz` is not the same artifact as the older `GSE125583_RAW.tar` sample archive. Before running the TAR-specific code below, inspect your existing compressed TSV with `pandas.read_csv(..., sep="\t", compression="gzip", nrows=5)`; determine actual header rows, whether the first column is Entrez gene IDs, and whether the remaining columns represent raw integer sample counts. Also compare column/sample identifiers against the GSE125583 Series Matrix. **Use the TSV directly if it contains suitable verified counts; there is no need to download the TAR solely to satisfy the example script.** Do not infer count integrality from the filename.

### 15.3 Python: inspect supplementary archive schema (runnable once downloaded)

```python
from pathlib import Path
import tarfile
import gzip
import io
import pandas as pd

archive = Path("data/raw/GSE125583_RAW.tar")
assert archive.exists(), "Download the supplementary GSE125583_RAW.tar first"

with tarfile.open(archive, "r:*") as tar:
    files = [f for f in tar.getmembers() if f.isfile() and
             (f.name.endswith(".tsv") or f.name.endswith(".tsv.gz"))]
    print("Number of sample TSV files:", len(files))
    print("First five names:", [f.name for f in files[:5]])
    if not files:
        raise RuntimeError("No TSV files in the archive; inspect its contents")
    member = files[0]
    raw = tar.extractfile(member)
    if member.name.endswith(".gz"):
        raw = gzip.GzipFile(fileobj=raw)
    preview = pd.read_csv(raw, sep="\t", nrows=8)
    print("Columns:", preview.columns.tolist())
    print(preview.to_string(index=False))
```

**Next:** write a parsing adapter using the **real** gene-ID and integer-count columns that appear. This guide does not invent those column names. If the available data only contain normalized nRPKM rather than integer counts, do not falsely run DESeq2 on them.

### 15.4 Recommended brain analysis

1. Parse actual counts and GEO sample labels.
2. Align sample identifiers and verify no duplicate individuals.
3. Check library sizes, missingness, detection rates, age/sex/RIN/batch.
4. Map Entrez Gene IDs using a documented release-matched crosswalk.
5. Annotate gene biotypes and report how many lncRNAs remain measurable.
6. Filter rare genes with an objective count threshold consistent with design (e.g., retain genes with at least 10 counts in a defensible number of samples, adjusted for class size).
7. Run `DESeq2` on integer counts with an adjusted design (if covariates support it).
8. For normalized visualization, use VST/rlog-like transformation as appropriate; **do not use transformed values as counts**.
9. Compare direction/interpretation with credible literature and, only where truly comparable, the blood findings.
10. Use same-brain-cohort coding genes for co-expression and pathway exploration.

### 15.5 Optional R: correct RNA-seq differential-expression model

Assume you have created **real** `brain_counts.csv` with columns = exact sample IDs, rows = genes, and **non-negative integer counts**, plus `brain_metadata.csv` with `sample_id` and `diagnosis` (`AD`/`Control`) and any supported covariates.

```r
library(DESeq2)
counts <- read.csv("data/processed/brain_counts.csv",
                   row.names = 1, check.names = FALSE)
meta <- read.csv("data/processed/brain_metadata.csv",
                 stringsAsFactors = FALSE)
meta <- meta[meta$diagnosis %in% c("AD", "Control"), ]
rownames(meta) <- meta$sample_id
stopifnot(all(colnames(counts) %in% rownames(meta)))
meta <- meta[colnames(counts), , drop = FALSE]
stopifnot(identical(colnames(counts), rownames(meta)))
stopifnot(all(is.finite(as.matrix(counts))))
stopifnot(all(as.matrix(counts) >= 0))
stopifnot(all(as.matrix(counts) == round(as.matrix(counts))))

meta$diagnosis <- relevel(factor(meta$diagnosis), ref = "Control")
# Add age/sex/RIN/batch here ONLY if correctly encoded and estimable.
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = meta,
  design = ~ diagnosis
)
keep <- rowSums(counts(dds) >= 10) >= 10  # illustrative threshold; justify it
# Better: choose minimum samples based on group sizes and detectability.
dds <- dds[keep, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("diagnosis", "AD", "Control"),
               alpha = 0.05)
out <- as.data.frame(res)
out$gene_id <- rownames(out)
write.csv(out, "outputs/tables/brain_deseq2_all_genes.csv", row.names = FALSE)
```

For high-quality final inference, add verified covariates to the `design` formula **before** `DESeq(dds)` and ensure the model matrix has full rank. Apply the lncRNA biotype filter using a **verified matching ID namespace**, documenting the number of noncoding genes in the analyzed count matrix. If proper integer counts are unavailable, select a method for the actual data scale instead.

### 15.6 Optional R: preferred limma analysis for blood arrays

After preparing a **verified** log-like normalized feature matrix and metadata:

```r
library(limma)
expr <- read.csv("data/processed/blood_lncrna_expression.csv",
                 row.names = 1, check.names = FALSE)
meta <- read.csv("data/processed/blood_metadata.csv")
rownames(meta) <- meta$sample_id
meta <- meta[colnames(expr), , drop = FALSE]
stopifnot(identical(colnames(expr), rownames(meta)))
meta$diagnosis <- factor(meta$diagnosis, levels = c("Control", "AD"))
stopifnot(!anyNA(meta$diagnosis))

design <- model.matrix(~ diagnosis, data = meta)
fit <- lmFit(as.matrix(expr), design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "diagnosisAD", number = Inf,
                    adjust.method = "BH")
write.csv(results, "outputs/tables/blood_limma_results.csv")
```

Include age, sex, and batch terms **only when their metadata and design justify them**. If expression is not log-like normalized intensity, verify the appropriate limma preprocessing before fitting. **Do not mix `DESeq2` and `limma` p-values indiscriminately across different tissues and platforms.**

> **Implementation prompt for this section:** [Agent Stage 15 — Codex and Antigravity](#agent-stage-15)  


---

## 16. Step 9: Optional MATLAB Bioinformatics Toolbox component

MATLAB can showcase course-specific skills without replacing more appropriate statistical methods. Import the **already curated** expression matrix and run a demonstration of sample QC and a single-gene comparison. This example assumes the expression CSV has **genes in rows** and **GSM samples in columns**. Certain operations below use the **Statistics and Machine Learning Toolbox**, not only the Bioinformatics Toolbox; verify installed products.

```matlab
% MATLAB demonstration on already validated, curated expression CSV.
T = readtable('data/processed/blood_lncrna_expression.csv', ...
              'ReadRowNames', true, 'VariableNamingRule', 'preserve');
X = table2array(T);      % genes x samples
names = T.Properties.RowNames;

figure;
boxplot(X(:, 1:min(25, size(X,2))));
title('Sample expression distributions (subset, supplied scale)');
ylabel('Normalized expression');
xtickangle(90);

% PCA expects observations x features; transpose X.
X_samples = X';
% Handle missing values with a documented method BEFORE PCA.
if any(isnan(X_samples), 'all')
    error('Missing values: apply a justified imputation/QC policy first');
end
[coeff, score, latent, ~, explained] = pca(X_samples);
figure;
scatter(score(:,1), score(:,2), 25, 'filled');
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('PCA - gene-expression sample QC');

% For a named, validated lncRNA, find its ACTUAL gene ID/row first.
% Do not assume 'NEAT1' is present simply because it is well known.
idx = find(strcmp(names, 'YOUR_VERIFIED_GENE_ID'), 1);
if ~isempty(idx)
    fprintf('Expression of selected gene loaded for %d samples.\n', size(X,2));
else
    disp('Selected validated gene ID not found; check the annotation mapping.');
end
```

**If your professor asks “why MATLAB plus Python?”:** MATLAB provides familiar genomic signal processing, matrix operations, and course-specific demonstrations. Python provides an accessible end-to-end data and ML stack; R/limma or DESeq2 are used where they offer more established genomics-specific statistical inference.

> **Implementation prompt for this section:** [Agent Stage 16 — Codex and Antigravity](#agent-stage-16)  


---

## 17. What to do if the first approach fails

**Failure is an experimental result, not a reason to hide or manipulate data.** Establish contingency choices **before** inspecting external performance.

| Failure scenario | Scientific reason | Appropriate alternative | What *not* to do |
|---|---|---|---|
| Too few verified blood lncRNA probes | Old chips were not designed for comprehensive lncRNA detection | Primary brain RNA-seq study; or audit GSE161199 / another suitable PBMC RNA-seq dataset | Rebrand coding genes as lncRNAs |
| No lncRNAs survive FDR | Effects weak, many tests, possible noise | Report negative findings; use literature-guided **pre-specified** candidate tests as clearly labeled exploratory work | Drop FDR correction to force significance |
| Poor internal AUC | Limited signal, sample size, confounding, overfitting | Simplify panel, improve covariates/QC, repeat properly nested procedure on discovery data | Tune hyperparameters on final holdout |
| High internal but poor external AUC | Cohort/platform/domain shift or overfit | Show both scores, examine platform comparability, report failure and scope | Hide external results |
| Same genes not measurable externally | Platform/probe mismatch | Restrict validation to genuinely common verified genes if adequate; otherwise state no direct validation | Treat same gene names as identical probe measurements |
| Sample age/sex missing | Confounding cannot be fully adjusted | Use available covariates, state limitations, investigate alternative datasets | Fabricate missing covariate values |
| Brain results disagree with blood | Tissue and cell composition differ | Frame as tissue specificity, not contradiction requiring deletion | Pool brain and blood samples as identical features |
| Batch matches diagnosis exactly | Non-identifiable confounding | Restrict to within-batch studies or seek better cohort | Claim ComBat proves biological separation |
| No independent validated dataset | Inadequate external test design | Strong nested CV and explicit limitation; plan future validation | Call within-study cross-validation independent external validation |

### A scientifically credible fallback project

**Title:** *An audit of lncRNA biomarker measurability and Alzheimer's-associated expression signatures across legacy blood microarrays and brain RNA-seq.*

A negative feasibility result can be surprisingly strong if you quantify why legacy arrays cannot provide trustworthy broad lncRNA discovery. It demonstrates **assay literacy, reproducibility, domain shift analysis, and honest biomedical ML**. A recruiter or professor may value this more than an unsupported “AUC 0.99” claim.

### A research extension if time allows

Instead of blindly adding more models, evaluate whether **annotation confidence + fold-stable feature selection + covariate adjustment** improve replicability. Report a small multi-objective comparison:

- Robustness: feature-selection frequency across discovery resamples.
- Biological support: evidence from published studies in the *same tissue* where possible.
- Predictive utility: external AUC + confidence interval, when true compatible external measurement exists.
- Incremental utility: added performance above age/sex baseline.

**Do not optimize these criteria against the final external labels.**

> **Implementation prompt for this section:** [Agent Stage 17 — Codex and Antigravity](#agent-stage-17)  


---

## 18. Plots, report, presentation, and research claims

### 18.1 Minimum publication-quality figures

| Figure | What it communicates | Honest caption reminder |
|---|---|---|
| 1. Study/design flowchart | Which samples and datasets enter each analysis | Show exclusions, cohort independence, platform differences |
| 2. Sample composition table/chart | AD, control, MCI, ambiguous labels | Distinguish total GEO size vs binary subset |
| 3. Probe-mapping audit funnel | Starting probes → candidate genes → unique validated lncRNAs | Show how many are **provisional** vs **verified** |
| 4. QC sample intensity boxplots | Compare distributions and possible outliers | Similar boxes do **not** prove all confounding removed |
| 5. PCA by diagnosis and batch | Structure, possible technical effects | PCA is exploratory, not a classifier |
| 6. Differential-expression volcano | Direction/evidence for lncRNAs | Exact expression scale and adjusted p-values |
| 7. Candidate heatmap | Patterns of exploratory features | If chosen using labels, not evidence of out-of-sample prediction |
| 8. Co-expression graph | lncRNA–mRNA correlations | Network edges mean **association**, not direct molecular interaction |
| 9. Nested-CV ROC | Performance without each outer fold's training leakage | Label as **out-of-fold**, not independent clinical validation |
| 10. External ROC/confusion matrix | Genuine holdout performance | Only create if external mapping/scale gate passes |
| 11. Feature-selection stability | Reproducibility of panel across resamples | Selection frequency is not biological causality |
| 12. Covariate comparison | Age/sex baseline vs RNA-only vs combined | Same patients and folds required |

### 18.2 Suggested report structure (15–25 pages, adaptable)

1. **Title page:** project title, OE3E40, five names/roll numbers, institution, submission date.
2. **Abstract (200–300 words):** question, dataset(s), methods, **actual** findings if available, limitations.
3. **Introduction:** Alzheimer's context, established biomarkers, gene expression and lncRNAs.
4. **Literature review:** Canoy et al. and a few primary studies; what was truly established.
5. **Research gap and objectives:** focus on verified lncRNA identification and leakage-resistant validation.
6. **Data:** platform, tissue, sample counts, metadata, exclusions, accession citations.
7. **Methodology:** audit, normalization, mapping, confounders, differential expression, correlation, ML, validation.
8. **Results:** actual exported tables/plots only. Include negative results and missing-data statistics.
9. **Discussion:** biological hypotheses; gene/assay uncertainties; relation to prior literature.
10. **Limitations:** cross-sectional labels, microarray coverage, tissue mismatch, age/sex, RIN, cell mixture, cohort shift, unavailable qRT-PCR.
11. **Conclusion:** carefully calibrated to evidence.
12. **Future work:** PBMC RNA-seq replication, larger cohorts, target assays, qRT-PCR, prospective clinical studies.
13. **References:** DOI, GEO accession, package versions.
14. **Appendix:** parameter settings, code-link, data dictionary, audit and reproducibility records.

### 18.3 Abstract template — fill only after running experiments

> **Background:** Long non-coding RNAs (lncRNAs) have been associated with Alzheimer's disease (AD), but their potential as reproducible biomarkers depends on assay specificity and independent validation. **Methods:** We analyzed publicly available [datasets], documented [number] eligible AD cases and [number] controls, and audited lncRNA measurements against [annotation release]. We performed [QC/DE method], explored [co-expression analysis], and evaluated [model families] using [nested CV design] and [external study, if validated]. **Results:** [Insert true measured numbers, confidence intervals, and whether independent replication succeeded.] **Conclusions:** [Clearly state what the data support and what they do not.] Further independent assay and laboratory validation are required.

### 18.4 Suggested presentation (10 slides, 7–10 minutes)

1. Title, course and group names.
2. Problem: what is AD and why study lncRNAs?
3. DNA→RNA→lncRNA diagram and an illustrative biomarker mechanism.
4. Research gap, precise objectives, **clinical-scope warning**.
5. Selected public datasets and critical chip/coverage limitations.
6. Complete workflow flowchart and quality/decision gates.
7. Actual QC and annotation-audit results.
8. Actual DE/co-expression results **or** feasibility findings.
9. ML panel, nested CV, hyperparameter tuning, external validation outcomes (if available).
10. Honest limitations, fallback path, conclusion, references.

**Never decorate the results slide with made-up volcano points or an idealized ROC curve without marking it prominently “Illustration, not results.”**

> **Implementation prompt for this section:** [Agent Stage 18 — Codex and Antigravity](#agent-stage-18)  


---

## 19. GitHub, demo application, and portfolio

### 19.1 README outline

```markdown
# Alzheimer's lncRNA Biomarker Discovery

## Research question
Can carefully annotated lncRNA expression features discriminate AD
from controls in independent public cohorts?

## Status
Research prototype. Not a diagnostic or medical device.

## Data sources
GSE63061, GSE63060, GSE140829, GSE125583 [links and descriptions].

## Method
GEO metadata audit -> verified lncRNA mapping -> QC -> DE ->
co-expression -> nested CV -> optional independent validation.

## Reproduce
1. Create Python environment.
2. Run scripts/01_fetch_geo.py for approved series.
3. Review and curate platform probe mappings (document evidence).
4. Run scripts/03_analyze_blood.py.
5. Review output tables and audit log.

## Results
Only verified measured findings with confidence intervals.
[No numerical claim until code has been run.]

## Limitations
Legacy array lncRNA coverage; potential confounding; tissue-specific
expression; small samples; no wet-lab or clinical validation.

## License and citations
Code license, GEO dataset acknowledgements, paper citations.
```

### 19.2 `.gitignore`

```gitignore
.venv/
__pycache__/
.ipynb_checkpoints/
.env
*.pyc
*.joblib
*.pkl
# Downloadable public data should be fetched from original source;
# avoid bloating Git and check any dataset reuse restrictions.
data/raw/
data/interim/
data/processed/
# Keep small selected aggregate figures/tables if their reuse is allowed.
```

**Reproducibility:** never publish a result table without its generation code, annotation version, sample inclusion rules, random seed, and model parameters. Keep a small provenance/checksum file. For public datasets, cite original submitters and reuse papers. Confirm any project-specific data reuse constraints before redistribution.

### 19.3 Optional Streamlit demo (research-only)

A lightweight app can let a visitor select a **validated** gene and view distribution plots. This is preferable to a medical “Enter patient details to diagnose AD” tool, which could mislead people about clinical validity.

Save as `app.py`:

```python
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

st.title("Alzheimer's lncRNA Expression Explorer")
st.caption("Educational / research visualization only. Not a medical test.")

@st.cache_data
def load_data():
    # Prepare this CSV from your audited discovery cohort:
    # columns = sample_id, diagnosis, gene_id, expression
    return pd.read_csv("outputs/tables/explorer_long_data.csv")

try:
    df = load_data()
    required = {"sample_id", "diagnosis", "gene_id", "expression"}
    if not required.issubset(df.columns):
        st.error("The explorer data file is missing required columns.")
        st.stop()
    gene = st.selectbox("Select verified lncRNA gene ID", sorted(df["gene_id"].unique()))
    d = df[df["gene_id"] == gene]
    fig, ax = plt.subplots()
    d.boxplot(column="expression", by="diagnosis", ax=ax)
    fig.suptitle("")
    ax.set_title(f"{gene}: descriptive group expression")
    ax.set_ylabel("Processed expression (source scale)")
    st.pyplot(fig)
    st.caption("The plot shows an association, not an individual diagnosis.")
except FileNotFoundError:
    st.info("Generate outputs/tables/explorer_long_data.csv from verified data first.")
```

Run `streamlit run app.py`. The `explorer_long_data.csv` is an **output contract** to build from your validated feature matrix, not a file claimed to exist already.

### 19.4 Resume bullet templates — use only after actual completion

- “Developed a reproducible Python pipeline to audit Alzheimer's gene-expression cohorts, validate lncRNA probe mapping, and perform statistical differential-expression analysis.”
- “Compared sparse logistic regression, random forest, and SVM using nested cross-validation; reported ROC-AUC, sensitivity, specificity, and model-selection stability with leakage controls.”
- “Investigated lncRNA–mRNA co-expression and biological pathway hypotheses, documenting annotation limitations and external cohort shift.”

**After actual execution**, add true verified cohort counts, gene counts, and metrics. Never insert an invented AUC or pretend that a research proposal was implemented.

### 19.5 Optional LinkedIn post draft (publish only when you have results)

> I worked on a computational genomics project investigating whether long non-coding RNA expression can serve as a candidate biomarker signal in Alzheimer's disease. The project emphasized a challenge that biomedical ML often overlooks: **the expression platform must actually measure the biological features we claim to model**. Our workflow combined public GEO data auditing, lncRNA annotation, differential expression, co-expression exploration, and leakage-aware machine learning. We documented both promising results and limitations rather than treating an internal AUC as clinical validation. [Add genuine result figure and GitHub link when available.] This remains an exploratory research project, not a diagnostic tool.

> **Implementation prompt for this section:** [Agent Stage 19 — Codex and Antigravity](#agent-stage-19)  


---

## 20. Suggested four-week timeline and completion checklist

**Assumption only:** deadline and daily availability were not supplied. This is an adaptable four-week schedule for a small academic project; it contains **no assignments to individual group members**.

| Time | Focus | Required output | Go/no-go rule |
|---|---|---|---|
| Days 1–3 | Understand AD/lncRNA; audit source papers and GEO metadata | Dataset audit, exact labels, code environment | Never proceed with unknown label meaning |
| Days 4–7 | Download data; GENCODE/GPL annotation feasibility | Probe mapping audit + lncRNA measurement feasibility | If no unique lncRNA coverage, activate fallback |
| Days 8–10 | QC, age/sex/batch, exploratory PCA | QC figures and recorded exclusions | Exclusions must be technical/predefined |
| Days 11–14 | DE and descriptive biomarker shortlist | FDR table, effect plots, preliminary literature table | Keep interpretation limited to available tissue |
| Days 15–17 | Co-expression, GO/Reactome analysis | Validated mRNA list, correlation table, enrichment | Correlation cannot be called direct regulation |
| Days 18–21 | Nested-CV ML, tuning and stability analysis | ROC, metrics, frozen panel, selection frequency | All selection/tuning within discovery folds |
| Days 22–24 | External validation if mapping passes; error analysis | External results **or** documented infeasibility | No second attempt using same external labels as a fresh test |
| Days 25–27 | Report, references, PPT, reproduction test | Finished methods and discussion, real exported plots | All figures reproducible |
| Days 28–30 | Viva, README, optional Streamlit app | Demo, GitHub, question bank | No unsupported clinical claims |

### End-to-end checklist

- [ ] Dataset accession, tissue, platform, original publication, and download date recorded.
- [ ] All original labels audited; MCI/ambiguous excluded from binary task.
- [ ] Repeated patients/related samples checked.
- [ ] GEO normalization notes captured and respected.
- [ ] Platform annotation **and probe specificity** reviewed.
- [ ] GENCODE version and gene-ID namespace recorded.
- [ ] Number of uniquely validated lncRNA features counted.
- [ ] Age, sex, RIN, batch, tissue and cell-composition limitations discussed.
- [ ] QC plots and explicit outlier policy exported.
- [ ] Appropriate statistical model + BH FDR run.
- [ ] Shortlist includes stable IDs and mapping evidence.
- [ ] Co-expression generated from **same-tissue, same-sample** RNA data.
- [ ] Any enrichment uses coding genes and an appropriate background universe.
- [ ] Imputation, scaling and supervised feature selection performed inside CV.
- [ ] Model hyperparameters tuned with inner folds.
- [ ] External validation held out and only run when comparability verified.
- [ ] AUC, sensitivity, specificity, F1 and confidence interval reported honestly.
- [ ] Demographic-only baseline compared on same patient subset when available.
- [ ] Report distinguishes exploratory discovery, replication and independent testing.
- [ ] GitHub README, requirements lock, reproducibility records, caveats completed.
- [ ] Final presentation uses only genuine outcomes or clearly marked illustrations.

> **Implementation prompt for this section:** [Agent Stage 20 — Codex and Antigravity](#agent-stage-20)  


---

## 21. Viva questions with answers

**Q1. What does lncRNA stand for?**  
Long non-coding RNA: generally an RNA transcript longer than 200 nucleotides without established protein-coding function. It can influence gene regulation.  
*Hinglish:* Lamba RNA jo usually protein nahi banata lekin gene regulation mein role ho sakta hai.

**Q2. Why are lncRNAs worth studying in Alzheimer's?**  
Prior literature links some lncRNAs to amyloid, tau, neuronal stress, inflammation and other AD-related pathways. Association alone is not clinical validation.

**Q3. What is a biomarker?**  
A measurable indicator associated with a biological condition; a candidate must show reliable measurement and independent predictive utility before clinical use.

**Q4. What is the difference between mRNA and lncRNA?**  
mRNA normally encodes protein; lncRNAs usually do not but can regulate expression or RNA/protein interactions.

**Q5. What is GEO?**  
NCBI Gene Expression Omnibus: a public repository for gene-expression and related functional genomics experiments.

**Q6. What is a microarray probe?**  
A short sequence designed to bind a complementary RNA-associated sequence. Probe identity and binding specificity are necessary before claiming the probe measures a particular lncRNA.

**Q7. Why isn't GSE63061 an RNA-seq dataset?**  
It uses the Illumina HumanHT-12 v4 **expression microarray**; it stores processed fluorescence-derived feature intensities, not raw sequencing read counts.

**Q8. Why can we not use every feature as an lncRNA?**  
A general microarray also includes protein-coding or ambiguously mapped probes; some may cross-hybridize. We need reliable biotype and probe mapping.

**Q9. Why use GENCODE?**  
It provides curated gene/transcript annotations, including biotypes. We must document the release, ID namespace and relationship to the array's platform annotation.

**Q10. What is the difference between GSE63060 and GSE63061?**  
They are related AddNeuroMed batches with different platform versions, useful for robustness/replication but not fully independent external clinical validation.

**Q11. What is MCI, and why not group it with controls?**  
Mild cognitive impairment: a distinct clinical category; silently labeling it “healthy” corrupts the binary reference labels. It can be considered as a separate exploratory classification question.

**Q12. What is PCA?**  
Principal Component Analysis reduces many expression features to a few axes that capture large variation. In this project it is mainly a QC visualization, not proof of diagnosis.

**Q13. Why not normalize the GEO microarray again automatically?**  
The deposited data may already be variance-stabilized and quantile-normalized; an unnecessary second transformation can distort the data.

**Q14. What is a confounder?**  
A variable related to the studied condition and measurements that can create or distort an apparent association, such as age, sex, batch or blood-cell composition.

**Q15. Why is XIST useful in QC, and why is it risky for classification?**  
XIST expression often reflects sex. A strong sex difference may confirm the biology/probe behavior but can masquerade as an AD effect if case/control sex distributions differ.

**Q16. What is differential expression?**  
Testing whether gene expression differs between groups using a statistical model with appropriate covariate adjustment and multiple-testing control.

**Q17. Why FDR instead of just p < 0.05?**  
Thousands of tests produce false positives by chance. BH-FDR controls the expected fraction of false discoveries among rejected hypotheses under its assumptions.

**Q18. Why use limma on microarrays and DESeq2 on counts?**  
limma's empirical-Bayes models suit normalized microarray intensities; DESeq2 explicitly models integer RNA-seq counts with a negative-binomial framework and size factors.

**Q19. What is co-expression?**  
A statistical tendency for RNA abundances to vary together in the same samples. It does not by itself establish direct interaction or causation.

**Q20. Why do GO/Reactome enrichment?**  
To test if co-expressed coding genes are over-represented in annotated biological functions relative to a defined background; the result is hypothesis-generating.

**Q21. What is LASSO?**  
L1-penalized regression; it shrinks some coefficients to zero, producing a sparse and interpretable set of candidate genes.

**Q22. What is hyperparameter tuning?**  
Searching over model settings (e.g., LASSO `C`, forest depth, SVM `gamma`) using training-only cross-validation to improve generalization rather than fit the test data.

**Q23. What is nested cross-validation?**  
An outer loop estimates held-out performance while an inner loop chooses hyperparameters and performs supervised feature selection. It reduces optimistic performance estimates.

**Q24. What is data leakage?**  
Using information from a validation/test sample during training, scaling, gene ranking or parameter selection. It makes scores look better than genuine unseen performance.

**Q25. Why compare against age and sex?**  
To show whether expression contributes predictive information beyond demographic differences that might otherwise drive a seemingly high AUC.

**Q26. What is sensitivity vs specificity?**  
Sensitivity = fraction of true AD cases detected. Specificity = fraction of controls correctly recognized.

**Q27. What does ROC-AUC = 0.5 mean?**  
It indicates approximately chance-level ranking for the binary task; 1.0 is perfect ranking on the evaluated sample. AUC does not tell us whether predicted probabilities are calibrated.

**Q28. Why isn't GSE125583 an external blood validation dataset?**  
It measures postmortem brain tissue with RNA-seq, not the same kind of blood specimen or assay. It can support biological plausibility, not straightforwardly validate a blood test.

**Q29. What if you obtain no statistically significant lncRNAs?**  
Report it. Assess measurement coverage and study power, describe limitations, and propose better-suited RNA-seq or experimental follow-up without relaxing standards to force hits.

**Q30. Can this detect Alzheimer's early or diagnose patients now?**  
No such claim follows from retrospective AD/control gene-expression reanalysis. A genuine diagnostic or early-detection claim requires appropriate independent, prospective, clinical, and analytical validation.

> **Implementation prompt for this section:** [Agent Stage 21 — Codex and Antigravity](#agent-stage-21)  


---

## 22. Two-minute pitch

> Good morning. Our project is called **Establishing the Role of Long Non-coding RNAs as Biomarkers in Alzheimer's Disease**. Alzheimer's is a progressive neurological disease affecting memory and cognition. Biological changes involve amyloid-beta, tau, inflammation, and neuronal damage. We are investigating whether some long non-coding RNA molecules—which can regulate genes even though they usually do not make proteins—show reproducible expression patterns associated with Alzheimer's.
>
> We will use **real public gene-expression datasets** from NCBI GEO. The blood datasets GSE63061, GSE63060, and GSE140829 are based on microarray technologies. We also examine brain RNA-seq dataset GSE125583 for biological context. A major strength of our approach is that we will not automatically call every microarray signal a lncRNA: we first verify each feature using platform annotation, GENCODE biotypes, and probe specificity.
>
> After inspecting patient labels and quality, we will compare expression between Alzheimer's cases and controls using methods appropriate for the actual data type. We will then investigate correlations between candidate lncRNAs and protein-coding genes and explore biological pathways. For machine learning, we will compare a sparse LASSO logistic model against random forest and SVM with **nested cross-validation** and hyperparameter tuning. If the same validated features are available in a truly separate cohort, we will perform one locked independent evaluation and report AUC, sensitivity, specificity, and confidence intervals.
>
> Our expected contribution is not a clinically approved test. It is a **reproducible, honest computational evaluation** of whether lncRNA expression provides a promising Alzheimer's-associated signature, where current datasets are strong, and where better RNA-seq and laboratory qRT-PCR validation are still needed. Thank you.

### Ultra-short Hinglish version

*“Hum AD vs control ke public RNA data ko analyze kar rahe hain. Sabse pehle check karenge ki microarray chip sach mein lncRNA measure kar rahi hai ya nahi. Phir QC, differential expression, gene correlation aur ML using LASSO, Random Forest, SVM karenge. Nested CV se leakage/overfitting control karenge, aur compatible independent cohort pe final test karenge. Isko candidate-biomarker research bolenge, clinical diagnosis nahi.”*

> **Implementation prompt for this section:** [Agent Stage 22 — Codex and Antigravity](#agent-stage-22)  


---

## 23. References, dataset citations, and verification log

### Core biological literature

**[R1]** Canoy, R. J., et al. (2024). *Non-coding RNAs involved in the molecular pathology of Alzheimer's disease: a systematic review*. **Frontiers in Neuroscience**, 18, 1421675. DOI: <https://doi.org/10.3389/fnins.2024.1421675>. **This is the supplied 12-page review.** Its Table 1 on page 6 reports 362 lncRNAs in reviewed literature; the review does not establish their independent clinical validity.

**[R2]** A current research background to interpret alongside historical reviews: National Institute on Aging / Alzheimer's Association biomarker framework and clinical diagnostic literature. **Check current guidance** before writing statements about which blood-protein biomarkers are in clinical use. lncRNA research remains a distinct, exploratory question. An accessible background starting point is <https://pmc.ncbi.nlm.nih.gov/articles/PMC11775455/> (verify publication/context).

### Data source records and peer-reviewed descriptions

**[D1]** GSE63061 — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63061>. Platform GPL10558 — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558>. GEO reports 388 series samples; binary AD/control subsets (139/134) are reported in published analyses. See <https://pmc.ncbi.nlm.nih.gov/articles/PMC11380308/> and <https://pmc.ncbi.nlm.nih.gov/articles/PMC12964574/>. The GEO sample processing notes indicate variance-stabilization and quantile normalization. **Verify labels in the actual matrix**.

**[D2]** GSE63060 — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63060>. GPL6947 — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6947>. Published AD/control subset 145/104; related to AddNeuroMed batch 1; some additional MCI records. Published comparison: <https://pmc.ncbi.nlm.nih.gov/articles/PMC11380308/>. Note differences between total GEO series size and binary subset.

**[D3]** GSE140829 — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140829>. GPL15988 — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15988>. A reported 587 samples comprise 204 AD, 249 controls and 134 MCI. Published summary: <https://pmc.ncbi.nlm.nih.gov/articles/PMC11633037/> and <https://pmc.ncbi.nlm.nih.gov/articles/PMC9249435/>. **Prove nuID-to-validated gene mapping before transferring models.**

**[D4]** GSE125583 — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125583>. Postmortem fusiform-gyrus RNA-seq. Study reporting 219 AD and 70 controls: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10225579/>. GEO lists `GSE125583_RAW.tar` (~87.4 MB) and sample-level Entrez-gene-ID counts/nRPKM. Examine each TSV's actual schema before selecting integer-count columns.

**[D5]** GSE161199 — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161199>. GEO describes coding and long non-coding RNA profiles in PBMC samples from AD, PD, ALS and controls; **the AD/control sample breakdown was not verified for this guide**, so no count is quoted.

**[D6]** GENCODE — <https://www.gencodegenes.org/human/>. Current landing page provides GTF annotations and gene/transcript resources. Record exact release and reference assembly, plus any ID-crosswalk release used. See GENCODE consortium description: <https://doi.org/10.1093/nar/gky955>.

### Statistical and computational methods

**[M1]** Love, M. I., Huber, W., & Anders, S. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. **Genome Biology**, 15, 550. <https://doi.org/10.1186/s13059-014-0550-8>.

**[M2]** Ritchie, M. E., et al. (2015). *limma powers differential expression analyses for RNA-sequencing and microarray studies*. **Nucleic Acids Research**, 43(7), e47. <https://doi.org/10.1093/nar/gkv007>.

**[M3]** Benjamini, Y., & Hochberg, Y. (1995). *Controlling the false discovery rate: a practical and powerful approach to multiple testing*. **Journal of the Royal Statistical Society B**, 57(1), 289–300. <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>.

**[M4]** Langfelder, P., & Horvath, S. (2008). *WGCNA: an R package for weighted correlation network analysis*. **BMC Bioinformatics**, 9, 559. <https://doi.org/10.1186/1471-2105-9-559>.

**[M5]** Pedregosa, F., et al. (2011). *Scikit-learn: Machine Learning in Python*. **JMLR**, 12, 2825–2830. <https://jmlr.org/papers/v12/pedregosa11a.html>.

**[M6]** Johnson, W. E., Li, C., & Rabinovic, A. (2007). *Adjusting batch effects in microarray expression data using empirical Bayes methods*. **Biostatistics**, 8(1), 118–127. <https://doi.org/10.1093/biostatistics/kxj037>.

**[M7]** Official library documentation: GEOparse <https://geoparse.readthedocs.io/>, Bioconductor <https://bioconductor.org/>, scikit-learn <https://scikit-learn.org/stable/>, statsmodels <https://www.statsmodels.org/stable/>. **Check installed library versions and API changes before running code.**

### Evidence status and how to verify it

| Claim type | This document's status | How to verify before final submission |
|---|---|---|
| Dataset accession/tissue/platform | Grounded in GEO-indexed listings and published reports | Open the accession and platform pages; save date/screenshots |
| GSE63061 / GSE140829 AD/control counts | Reported by published analyses | Recompute from exact original GEO metadata, excluding ambiguities |
| GSE125583 219/70 | Reported by published reanalysis | Inspect the actual sample manifest and participant IDs |
| 362 lncRNAs linked to AD literature | Reported in the **user-supplied review** | Check page 6, Table 1 of Canoy et al. (2024) |
| Specific lncRNAs measured uniquely on blood chips | **Unverified until your probe audit** | Check probe sequences / unique transcript mapping |
| Number of DE lncRNAs, pathways, correlations | **Not calculated here** | Run verified pipeline, retain tables and adjusted p-values |
| CV AUC, external AUC, sensitivity | **Not calculated here** | Run discovery with true labels; test external once if compatible |
| Clinical efficacy / early detection | **Not established** | Requires stronger prospective / assay / clinical validation |

### Final reproducibility declaration to paste into your report

> All datasets were obtained from publicly accessible repositories. We recorded accession IDs, processing notes, feature-annotation versions, and inclusion/exclusion criteria. The project distinguishes preliminary gene-symbol assignment from verified lncRNA probe specificity, exploratory biological associations from validated prediction, and related-cohort replication from independent external testing. No performance or statistical findings are claimed without a reproducible corresponding analysis output. All model development, preprocessing estimation, label-informed feature selection, and hyperparameter selection were restricted to training data under cross-validation. This is a computational hypothesis-generation study and not a clinical diagnostic device.

> **Implementation prompt for this section:** [Agent Stage 23 — Codex and Antigravity](#agent-stage-23)  


---

## Final one-page summary

**Research goal:** Can real, correctly measured lncRNAs help distinguish AD vs controls?

**Data:** GSE63061, GSE63060, GSE140829 blood microarrays; GSE125583 brain RNA-seq; GSE161199 contingency.

**Order that matters:** Download → audit labels and source-scale → validate lncRNA measurement → QC/confounders → DE/FDR → co-expression/enrichment → nested-CV LASSO/RF/SVM → locked independent test (if comparable) → report honest results.

**Most important risk:** incorrect lncRNA probe mapping and cross-cohort comparability; **never optimize around it**.

**Most important strengths:** verifiable sources, explicit decision gates, statistically valid tests, leakage-resistant ML, independent-validation logic, and transparent limitations.

**Deliverables:** reproducible code, dataset audit, candidate mapping table, QC figures, differential-expression table, optional co-expression/pathway results, ML comparison and any actual external results, report, slides, README, and optional educational app.

**Conclusion to remember:** *A robust negative result with a correct experimental design is more scientifically valuable than an impressive but invalid classifier score.*


---

## 24. Step-by-step Codex and Antigravity prompts

### How to use this prompt pack (read once)

1. Download this enhanced guide and place it in the repository root as **`PROJECT_GUIDE.md`**. Place the downloadable companion **`AGENTS.md`** beside it. Both agents can use the guide as a project specification; Codex also uses repo-level `AGENTS.md` as persistent instructions.
2. Open the actual project folder in Codex or in Google Antigravity. The user’s original files may be in `raw/`, not `data/raw/`; the first agent task must discover both locations and **not move/delete data without permission**.
3. Start with Stage **05** (file inventory) if you already have a project folder and downloads; otherwise start with Stages **01–07** for the research background and setup. Then use Stages **08–15** in scientific order. Optional MATLAB, dashboard, report and viva stages follow.
4. Copy the **entire Codex prompt** into Codex or the **entire Antigravity prompt** into the Antigravity agent chat. Do **not** paste all 23 prompts as one mega-task. Ask for one verified stage at a time. Both prompts aim at the same outputs but have different execution emphases.
5. After each stage, require the agent to show **file edits, executed test commands, actual result paths, logs, blockers and a PASS/BLOCKED gate**. Review data-label decisions and probe mapping yourself. If there is a blocker, run the recovery prompt below, rather than asking the agent to guess.
6. Do **not** paste sensitive credentials or patient-identifiable information into agents. Public GEO data should be kept in Git-ignored folders. Ask approval before live downloads, installations, destructive operations, pushing to GitHub and unlocking external validation.

**Tip:** In Codex you can keep repository rules in `AGENTS.md`; Google's Antigravity also works with project-context/instruction files when provided in the workspace. The exact UI may change—always attach or reference `PROJECT_GUIDE.md` and `AGENTS.md` in the first agent request.

**Official agent documentation:** [Codex AGENTS.md guidance](https://learn.chatgpt.com/docs/agent-configuration/agents-md), [Google Antigravity getting started](https://antigravity.google/docs/getting-started), [Google Antigravity project guidance](https://antigravity.google/docs/projects).

### First prompt to send to either agent (once per new workspace)

```text
You are implementing a scientifically rigorous academic research repository about lncRNAs as exploratory Alzheimer's biomarkers. Read PROJECT_GUIDE.md and AGENTS.md fully before making changes. You have access to already downloaded public GEO data in raw/ or data/raw/. Do not download, move, delete or overwrite those files without checking and receiving approval. Start by listing the repository, checking the available inputs, and making a short execution plan with scientific gates. Work on only the stage I specify. Every change needs a runnable command, a test/validation step, an explicit PASS/BLOCKED status, and a list of unverified scientific assumptions. Never fabricate outcomes or use GSE140829 labels until the frozen-model external validation stage is explicitly approved.
```

<a id="agent-stage-01"></a>
### Agent Stage 01 — Understand the problem, biological scope and glossary

**Guide section:** [1](#1-the-problem-in-simple-english-and-hinglish)

**Prompt for Codex (copy all):**

```text
Read Section 1 of PROJECT_GUIDE.md and AGENTS.md. Create `docs/biology_basics.md`: explain Alzheimer’s, DNA->RNA->lncRNA, gene expression, blood vs postmortem brain, and experimental biomarkers in beginner English plus simple Hinglish. Include a small Mermaid diagram and a glossary with 15 terms. Use Canoy et al. 2024 as background; explain that association is not diagnosis or causation. Do not invent study results or describe a proposed panel as clinically proven. Include 3 viva Q&As. Check Markdown headings/links and list the files changed.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Open the project and read Section 1 and AGENTS.md. Before writing, show a compact plan for a student-friendly biological introduction. Create `docs/biology_basics.md` with a clear flow diagram, glossary, English and Hinglish narration, and 3 viva questions. Use only supported background claims, with placeholders rather than fabricated new results. Review the file visually in the editor and show the modified files and any claims needing a citation.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-02"></a>
### Agent Stage 02 — Lock the research question and hypotheses

**Guide section:** [2](#2-precise-problem-statement-scope-and-hypotheses)

**Prompt for Codex (copy all):**

```text
Read Section 2 and AGENTS.md. Create `docs/research_protocol.md` with one primary question, H1/H2/H3, null hypotheses, outcome definition AD vs Control, exclusion of MCI/Unknown, primary/secondary endpoints, exploratory vs confirmatory labels, predefined FDR rule, and why an AD-versus-control study cannot claim presymptomatic detection. Include a small decision table describing what null or negative findings mean. Do not mention measured AUC or discovered lncRNAs. Add `tests/test_protocol_files.py` to check required headings only, and run it.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Review Section 2 and write `docs/research_protocol.md` as a pre-analysis protocol. Show the intended claim hierarchy and the exact evidence required for each claim. Record what the project does NOT establish (causality, early-stage diagnosis, clinical performance). Add 3 viva questions and an explicit “not yet analyzed” banner. Preview the document and report any ambiguities rather than guessing.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-03"></a>
### Agent Stage 03 — Define architecture and scientific decision gates

**Guide section:** [3](#3-research-design-and-decision-gates)

**Prompt for Codex (copy all):**

```text
Read Section 3. Make `docs/analysis_plan.md` containing a Mermaid architecture with separate blood and brain branches; G0–G5 gates; criteria for continuing or stopping blood-specific lncRNA ML; risk register for ambiguous probes, platform shift, subject duplicates, confounding and overfit. Add machine-readable `config/analysis_policy.yaml` holding policies, not observed numbers. Include no pipeline that merges blood with brain as a single training table. Validate YAML loading and mark unverified assumptions as pending.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Inspect the guide’s decision gates and propose an implementation plan before any code changes. Draw a simple blood/brain branching workflow in `docs/analysis_plan.md`, with an approval gate after each key scientific stage. Record what Antigravity may automate versus what requires human evidence (especially probe specificity and opening the external test). Deliver a status checklist and ask for approval before any external-label evaluation.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-04"></a>
### Agent Stage 04 — Verify public GEO sources and published cohort facts

**Guide section:** [4](#4-public-datasets-verified-facts-and-links)

**Prompt for Codex (copy all):**

```text
Read Section 4 and AGENTS.md. Implement `scripts/00_verify_sources.py`: maintain a `config/datasets.yaml` with GSE63061, GSE63060, GSE140829, GSE125583, associated GPL IDs, sample type, intended role and URLs. Do not assume published subset counts equal actual local labels. Prefer inspection of local files; network metadata checks are optional with approval. Export `docs/dataset_sources.md` and a CSV of source provenance plus verification status. Clearly distinguish 'from guide', 'checked locally', and 'confirmed externally'. Test YAML/config parsing.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Use Section 4 to create a dataset evidence ledger (`docs/dataset_sources.md`) with the four GEO links, tissue/platform/role and confidence in each sample-count claim. First look at local file presence and sizes; do not redownload. Explicitly flag GEO totals versus AD/control subsets and the related nature of GSE63060/61. Present sources separately from any unverified local observations. Supply the list of points needing manual GEO confirmation.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-05"></a>
### Agent Stage 05 — Audit exactly what is downloaded, and what is missing

**Guide section:** [5](#5-exactly-what-to-download-and-why)

**Prompt for Codex (copy all):**

```text
Read Section 5. Inspect `raw/` and `data/raw/` for the exact filenames in AGENTS.md. Make `scripts/00_inventory_local.py` that prints presence, byte size, gzip validity, SHA256 streaming hash, and 3–5 safe header lines per text/TSV file without extracting entire archives. Inspect `GPL15988_files/` contents recursively without deleting. Write `outputs/tables/local_inventory.csv` and `docs/local_file_inventory.md` grouping 'needed now', 'later', and 'potentially missing'. The existing `GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz` must be considered before asking for `GSE125583_RAW.tar`. Flag full GENCODE comprehensive GTF as conditionally needed to test coding off-targets. Run the inventory on available files; if files are absent, show exact commands and mark NOT RUN.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Scan the project explorer for `raw/` versus `data/raw/` and the 12 downloaded files/folders described in AGENTS.md. Create a visual keep/optional/missing table and an inventory script; inspect first rows, archive structure and types without reading all data into RAM. Do not delete, move or redownload anything without approval. Specifically inspect `GPL15988_files/` and the compressed GSE125583 NCBI count matrix. Report any schema/filename mismatch and the smallest next action needed.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-06"></a>
### Agent Stage 06 — Install dependencies in a clean, reproducible environment

**Guide section:** [6](#6-tech-stack-and-installation)

**Prompt for Codex (copy all):**

```text
Read Section 6. Detect OS, Python version and any existing venv before changing anything. Create `pyproject.toml` or `requirements.txt` with only immediately required Python packages and compatible Python 3.11 versions; optional extras for R, MATLAB, Streamlit, enrichment and SHAP. Add `scripts/check_environment.py` to report versions/import successes without triggering downloads. Provide Windows PowerShell and Colab setup commands in `docs/setup.md`. Avoid installing heavyweight optional libraries prematurely. Run `python -m compileall` and lightweight import checks where possible; describe failures without pretending to resolve them.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Review the existing Python interpreter and package state through Antigravity’s terminal before writing install instructions. Prepare `docs/setup.md` for Windows Python 3.11 and Colab, separating required libraries from optional limma/DESeq2, MATLAB, and dashboard packages. Add a single environment smoke test. Ask approval before installing packages, changing PowerShell policies or internet access. Show exactly which imports work in the current environment.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-07"></a>
### Agent Stage 07 — Build repo layout, schemas and data dictionary

**Guide section:** [7](#7-repository-layout-and-data-dictionary)

**Prompt for Codex (copy all):**

```text
Read Section 7, AGENTS.md, and local filenames. Create missing project directories under `data/`, `scripts/`, `tests/`, `docs/`, and `outputs/` without touching originals. Create `docs/data_dictionary.md` listing ID namespace, source expression scale, feature orientation, diagnosis categories, age/sex/batch/RIN and permitted missing values. Add `src/ad_lncrna/schema.py` with reusable checks: unique sample IDs; numeric expression; genes/features x samples in interim; sample alignment by ID; explicit type mismatch errors. Test these validators on synthetic fixtures ONLY. Ensure `.gitignore` excludes raw data and large intermediate artifacts.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Inspect the existing folder tree; propose minimal additions instead of recreating or overwriting working files. Implement `docs/data_dictionary.md` and simple schema validations with sample ID joins, units and ID types recorded. Preview the tree and run a small smoke test. Keep the user's original `raw/` files immutable and list any move/copy needed for later stage compatibility without performing it silently.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-08"></a>
### Agent Stage 08 — Load GEO series, inspect metadata and prove class labels

**Guide section:** [8](#8-step-1-acquire-and-audit-geo-datasets)

**Prompt for Codex (copy all):**

```text
Read Section 8. First inspect the first several lines of each local `GSE*_series_matrix.txt.gz` instead of re-downloading via GEOparse. Implement `scripts/01_audit_geo.py` capable of parsing local GEO Series Matrix headers (`!Sample_*` fields), `!series_matrix_table_begin/end` and numeric expression into `data/interim/<GSE>/expression` and `metadata` files, with streaming/chunking where practical. Preserve ALL original metadata fields and the original diagnosis text for discovery/replication cohorts only. Build a mapping dictionary only after printing unique labels; exclude MCI, borderline and unknown from binary analysis. For external GSE140829, inspect source/probe/sample-ID schema only; seal diagnosis labels until the approved final external test. Assert metadata/column IDs match and no duplicate IDs; show original vs eligible counts. Include tests for parser edge cases and missing phenotype keys. Do not fabricate the published sample counts. Export `docs/dataset_audit.md`, tables and checksums, then stop at Gate A.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Open local GSE63061, GSE63060 and GSE125583 Series Matrix files to audit diagnosis metadata; for GSE140829 inspect only platform, feature and sample-ID schema while keeping diagnosis labels sealed until external testing. Begin with a read-only preview and display every unique diagnosis value for the discovery/replication cohorts and their expression table shapes. After proposing the mapping rules for review, implement a local-file-first parser with assertions, sample exclusions and provenance output. Generate the honest class-count table and stop at Gate A if metadata cannot be unambiguously labeled. Show terminal command and sample outputs only from real execution.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-09"></a>
### Agent Stage 09 — Establish credible lncRNA feature annotation and probe specificity

**Guide section:** [9](#9-step-2-identify-credible-lncrna-measurements)

**Prompt for Codex (copy all):**

```text
Read Section 9 and the exact GPL/GENCODE files. Write `scripts/02_annotation_audit.py`: inspect actual GPL column names, probe ID conventions, available probe sequences and GENCODE v50 gene/transcript attributes. Parse candidate lncRNA gene IDs but label symbol-only matches PROVISIONAL. Detect ambiguous/multiple assignments and report absent sequences/nuID crosswalks. For GSE125583, inspect Entrez IDs and a trustworthy NCBI-to-Ensembl crosswalk; never match numerical Entrez IDs to ENSG strings. A validated blood map is only written when documented sequence uniqueness is actually tested against an appropriate comprehensive human transcript reference, including protein-coding transcripts. If alignment tooling/data is missing, output `BLOCKED` and a manual evidence checklist, not fake 'validated' CSV. Write `docs/probe_mapping_audit.md` with counts and confidence labels; add tests for duplicates/collisions.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Inspect GPL10558, GPL6947, GPL15988 (including `GPL15988_files/`), the GENCODE lncRNA GTF and transcript FASTA. Show real available sequence/ID fields first. Prepare an annotation evidence funnel with provisional vs verified categories, exact reason codes, and an approval gate. Never infer specific lncRNA probes solely from symbols or a list of famous genes. If comprehensive annotation or a nuID crosswalk is missing, stop with a precise file/step request. Do not proceed to lncRNA-only blood ML until Gate B passes.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-10"></a>
### Agent Stage 10 — Perform QC, scale checks and confounder analysis

**Guide section:** [10](#10-step-3-quality-control-normalization-and-confounders)

**Prompt for Codex (copy all):**

```text
Read Section 10. Once Gate A passes and the relevant expression matrix exists, add `scripts/03_qc.py`. Confirm whether values are raw counts, already normalized VST/quantile-array intensities, or another scale BEFORE transformations. Produce missingness summary, per-sample distributions, PCA on adequately filtered/appropriately transformed features, age/sex/batch balance plots when metadata exist, possible duplicate flags, and an explicit outlier review table. No label-based outlier deletion; no blind log2 or batch correction. Export figures to `outputs/figures/qc/`, tables to `outputs/tables/qc/` and `docs/qc_report.md` with interpretations that describe actual observed plots. Add shape/type tests and stop at Gate C if unresolved scale or confounding prevents valid analysis.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Review the parsed real matrices and their GEO processing notes. Build a QC plan with separate branches for normalized arrays and integer brain counts. Implement and preview sample boxplots, PCA colored by diagnosis/site where available and covariate-balance charts; save true figures and a QC decision log. Do not normalize twice, use holdout labels to choose outliers, or call PCA separation a predictive result. Request review of flagged samples before exclusion.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-11"></a>
### Agent Stage 11 — Run differential expression with correct data-model assumptions

**Guide section:** [11](#11-step-4-differential-expression-analysis)

**Prompt for Codex (copy all):**

```text
Read Section 11. Implement a unit-tested DE interface with distinct methods: normalized blood microarray => limma in R if available, otherwise statsmodels adjusted linear model clearly labeled as a baseline; genuine integer RNA-seq counts => DESeq2/PyDESeq2, never OLS on raw counts as final inference. Handle available age/sex/batch covariates with rank/full-case checks. Include log/measurement-scale labeling, BH-FDR across the tested universe, prespecified gene filters, effect sizes, tested N and model diagnostics. Output `outputs/tables/blood_de.csv` or brain file as applicable, a truthful volcano plot, and `docs/de_methods.md`. Exploratory full-data DE must NOT be reused as pre-CV feature selection. If no verified blood lncRNA probes, stop or run gene-level exploratory analysis under the correct non-lncRNA label.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Use the QC and annotation status to select scientifically appropriate statistical methods. Before executing, show a one-page method decision: input type, groups, adjustment terms, test family, FDR scope, exclusion criteria. Write/execute DE code only on real validated features, with explicit outputs and zero-result handling. Label coefficients with correct source units; do not call unverified processed-array differences 'log2 fold changes'. Preview exported tables and volcano plot, including an interpretation of null results.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-12"></a>
### Agent Stage 12 — Investigate same-tissue co-expression and pathway enrichment

**Guide section:** [12](#12-step-5-co-expression-and-pathway-analysis)

**Prompt for Codex (copy all):**

```text
Read Section 12. Build `scripts/05_coexpression.py` using same-cohort, same-subject lncRNA and protein-coding matrices only. Exclude unverified IDs; prefilter low-information genes; compute Spearman or justified adjusted correlations; apply BH correction over *all* tested pairs. Export lncRNA-mRNA edges plus sample counts and parameters. Enrichment is optional: if allowed internet and gene IDs are valid, query current GO/Reactome resources using an assay-specific measured-gene universe; otherwise export an input gene list and a documented manual step. Produce one network plot with association (not causal) edge labels, and `docs/pathways.md` noting possible cell-composition confounding. Tests must include missing-data and constant-gene cases.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Plan a laptop-friendly co-expression workflow before considering WGCNA. Check that every lncRNA and mRNA column comes from the same measured samples and tissue. Create edge table, FDR-adjusted statistics, a clear network diagram and pathway enrichment with an honest background universe. Ask for permission before calling external enrichment APIs. Mark all edges as correlations, not interactions or mechanisms, and provide a failure report if the cohort lacks viable coding/lncRNA pairs.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-13"></a>
### Agent Stage 13 — Train leakage-resistant ML and compare baselines

**Guide section:** [13](#13-step-6-machine-learning-biomarker-panels)

**Prompt for Codex (copy all):**

```text
Read Section 13. Implement `scripts/06_train_ml.py` and `tests/test_ml_leakage.py`. Inputs MUST be a real, verified lncRNA expression matrix plus audited AD/Control labels; halt if validated map missing. StratifiedGroupKFold if subjects repeat, otherwise justified StratifiedKFold. Outer 5-fold nested CV with inner folds for imputation, scaling, feature selection AND tuning; fixed seed 42. Compare age/sex-only on the exact same subject set and splits, L1 logistic, RF, and RBF SVM where feasible. Save out-of-fold scores, AUC, sensitivity, specificity, F1, confusion matrices, fold sizes, parameter grids and fold-specific selected gene IDs; record CI caveats. Train a final frozen model on discovery only, with serialized mapping+feature order+preprocessing. Do not access GSE140829 diagnosis labels. Label scores as internal CV and return Gate D status.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
First show a diagram of outer versus inner CV and which transformations are fitted in each training fold. Build a reproducible training run on the approved discovery cohort with baselines and constrained models, honest seed and patient-group handling. Run tests ensuring no test-fold samples enter fit/feature selection and the external cohort is not read. Produce model comparison charts and an explanation of how every selected feature earned validated lncRNA status. Ask approval before changing panel size rules or running external evaluation.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-14"></a>
### Agent Stage 14 — Validate on an independent blood cohort without peeking

**Guide section:** [14](#14-step-7-external-validation-error-analysis-explainability)

**Prompt for Codex (copy all):**

```text
Read Section 14. Prepare `scripts/07_external_validate.py` with a `--dry-run` mode that verifies frozen model, input feature IDs, GPL10558↔GPL15988/nuID equivalence evidence, source normalization, subject independence and prechosen decision threshold, WITHOUT reading external labels. Keep an explicit `--approve-final-test` switch and require user approval for the actual evaluation. On approval, run the one-time test with the *frozen* gene ordering and training-fitted transforms, export ROC-AUC with patient-level CI, sensitivity/specificity/F1, confusion matrix, calibration caveats and documented shift analysis. Never reselect features or tune on the holdout; if the mapping fails, write `docs/external_validation_blocked.md` and report no valid external assessment.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Construct a validation readiness dashboard/checklist rather than automatically scoring GSE140829. Display shared verified gene IDs, provenance of probe equivalence, preprocessing comparability, model freeze timestamp/hash and data independence. In the first pass inspect only features and non-outcome schema; do not open holdout diagnosis values. Request explicit approval for the final test, run it once if approved, and record all successes and disappointments without retraining to improve the test score.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-15"></a>
### Agent Stage 15 — Analyze GSE125583 brain RNA-seq count data correctly

**Guide section:** [15](#15-step-8-rna-seq-brain-support-and-optional-r-analysis)

**Prompt for Codex (copy all):**

```text
Read Section 15 AND inspect the existing `GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz` BEFORE assuming `GSE125583_RAW.tar` is needed. Create `scripts/08_brain_inspect.py` that outputs the actual columns, first five rows, feature namespace, count-integrality check and GSM/sample name matching with `GSE125583_series_matrix.txt.gz`. Use `Human.GRCh38.p13.annot.tsv.gz` and a justified NCBI↔Ensembl mapping to annotate gene types; quantify coverage and uncertain IDs. If true nonnegative integer raw counts and clean metadata exist, implement DESeq2 or compatible count-based analysis with age/sex/RIN if available; export true brain lncRNA DE, QC and matched same-tissue co-expression. Otherwise stop and explain why count assumptions failed. Do not call brain support independent blood validation.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Use the file explorer to inspect the predownloaded GSE125583 `.tsv.gz` matrix and GEO Series Matrix header; confirm orientation, whether numeric columns are raw integer counts, sample IDs and identifier namespace. Create an inspection report FIRST. Only after schema checks, plan brain lncRNA mapping, count filtering, QC and DESeq2. Make a separate branch for brain results, never concatenate blood and brain samples. Show count-gene coverage before claiming brain lncRNA findings.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-16"></a>
### Agent Stage 16 — Add an optional course-specific MATLAB demonstration

**Guide section:** [16](#16-step-9-optional-matlab-bioinformatics-toolbox-component)

**Prompt for Codex (copy all):**

```text
Read Section 16. Create `matlab/brain_blood_qc_demo.m` or clearly named script that reads the verified small curated CSV and its sample metadata. Verify which steps need Statistics and Machine Learning Toolbox (e.g. `pca`, `boxplot`) versus Bioinformatics Toolbox. Show PCA, a QC plot and (if verified gene exists) a single lncRNA expression plot; use actual gene IDs, sample-aligned labels and honest axis units. Detect absent toolboxes/files gracefully. Include `matlab/README.md` with usage and manual checks. Do not claim a MATLAB run unless MATLAB is installed and execution actually happened.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Prepare an optional MATLAB mini-demo using curated outputs from the earlier stages. Before editing, inventory the imported CSV schema and available MATLAB toolbox functions. Generate a readable `.m` script and short course-oriented explanation; do not silently replace validated Python/R analyses with inappropriate MATLAB statistics. Clearly mark whether execution was actually tested in MATLAB or only code-reviewed.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-17"></a>
### Agent Stage 17 — Handle scientific blockers and fallback designs

**Guide section:** [17](#17-what-to-do-if-the-first-approach-fails)

**Prompt for Codex (copy all):**

```text
Read Section 17. Implement `scripts/09_gate_summary.py` that summarizes decision gates from previous audit artifacts without fabricating missing ones. If too few uniquely mapped blood lncRNAs or absent shared features, generate `docs/fallback_decision.md` with a brain-first/assay-audit study design and exact implications for model claims. If no FDR hits, summarize negative finding without changing thresholds. If external AUC disappoints, preserve the result and describe shift. Add tests for missing artifacts and failed gates. Do not work around blockers by rebranding coding genes as lncRNAs.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Inspect earlier gate reports and produce a decision memo: continue original blood ML, shift to brain-first lncRNA analysis, or publish an annotation/feasibility audit. For each option show what data really support, deliverables and remaining unknowns. Treat negative results as useful. Do not generate simulated 'good results' or re-run the final test to chase a score. Ask the user to select the scientifically justified branch when choices remain.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-18"></a>
### Agent Stage 18 — Assemble research report, figures and presentation

**Guide section:** [18](#18-plots-report-presentation-and-research-claims)

**Prompt for Codex (copy all):**

```text
Read Section 18. Build `docs/final_report.md` following abstract, introduction, literature, data, methods, actual results, discussion, limitations, conclusion and references; build `docs/presentation_outline.md` for 10 slides with 30–60 second speaker notes each. Link only to existing verified figure/table paths; if a result is missing, write `NOT YET RUN` with a clear explanation rather than a made-up chart. Add English and Hinglish caption explanations for central diagrams. Generate a script to check internal links/figure paths and an evidence-to-claim table. Do not create fake citations or threshold claims.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Inspect completed outputs, not just the plan. Assemble a human-friendly report outline and 10-slide speaker script with simple diagrams (DNA→lncRNA, blood/brain workflow, QC→DE→ML→validation) and short explanations. For every graph, display its file source and whether it is real or illustrative. Ask for a college format only if supplied; otherwise use a neutral academic format. Preview the result and flag every placeholder before declaring it presentation-ready.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-19"></a>
### Agent Stage 19 — Prepare GitHub README and optional educational Streamlit app

**Guide section:** [19](#19-github-demo-application-and-portfolio)

**Prompt for Codex (copy all):**

```text
Read Section 19. Create `README.md` (aim, data links, gate warnings, setup, reproducible run commands, real results only), `.gitignore`, `LICENSE` suggestion without selecting a license without user approval, and optional `app.py` that reads `outputs/tables/explorer_long_data.csv` if it exists. The app must display research-only group expression distributions for *verified* lncRNAs and a visible 'not medical diagnosis' notice, without a patient prediction form. Add an app smoke test and check that private/large data, credentials and joblib artifacts are not tracked. Do not push to GitHub without explicit authorization.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Prepare portfolio packaging with a clear README, reproducibility checklist and educational, non-diagnostic Streamlit explorer if real curated data exist. Preview the interface and verify missing-file behavior; if no validated output exists, show a friendly empty-state instead of fake charts. Review Git status for raw data or secrets and request approval before any remote push. Draft resume bullet placeholders but never populate them with invented AUC/feature numbers.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-20"></a>
### Agent Stage 20 — Run the end-to-end reproducibility and delivery checklist

**Guide section:** [20](#20-suggested-four-week-timeline-and-completion-checklist)

**Prompt for Codex (copy all):**

```text
Read Section 20. Implement `scripts/10_verify_project.py` that checks expected inputs/outputs for each stage, validates hashes or manifests where available, verifies CLI help and imports, and generates `docs/reproduction_report.md` with PASS/FAIL/BLOCKED/NOT RUN. Add a `Makefile` or `scripts/run_steps.ps1` with SAFE stage-specific commands; no automatic external unlock. Ensure no hardcoded user paths, no silent network downloads, and no cross-platform destructive commands. Generate a simple four-week/milestone checklist without assigning individual people. Execute smoke tests, report failures exactly, and give next commands.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Perform a repository-wide review: can another student reproduce dataset audit, annotation, QC, DE, ML and figures in the stated order? Produce a stage-by-stage execution checklist and evidence report, showing each artifact’s creation command, status and blocker. Do not certify unexecuted analyses as completed. Make a minimal runbook for Windows PowerShell and Colab; do not run the final holdout test again.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-21"></a>
### Agent Stage 21 — Prepare biology and ML viva questions with real project status

**Guide section:** [21](#21-viva-questions-with-answers)

**Prompt for Codex (copy all):**

```text
Read Section 21 and the stage logs. Write `docs/viva_qa.md` with the original guide’s question coverage plus practical answers on gene-expression matrices, lncRNA vs coding, probe specificity, FDR, DESeq2 vs limma, PCA, nested CV, data leakage, external validation, AUC and negative results. Keep answers simple (3–6 sentences) in English with a one-line Hinglish explanation. Add 3 'show me your actual script/plot' viva tasks using actual filenames but do not claim outputs exist if they do not. Include a short glossary and a clear 'what we completed' block sourced from execution logs.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Turn the scientific sections and real execution logs into a 20–30 question viva rehearsal sheet. Have answers in simple English with Hinglish reminders, including why blood probes are not automatically lncRNAs and why internal CV is not clinical validation. Add likely professor follow-up questions and identify which findings can be shown live. Do not present planned experiments as already completed.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-22"></a>
### Agent Stage 22 — Create a truthful two-minute pitch and demo narration

**Guide section:** [22](#22-two-minute-pitch)

**Prompt for Codex (copy all):**

```text
Read Section 22 and completed artifacts. Prepare `docs/pitch_2min.md` with a 2-minute English script, a 45-second Hinglish recap, and a minute-by-minute demo order: research question -> one real GEO table -> actual QC -> mapping audit -> results/limitations. Insert references to real output filenames conditionally and say 'in progress' where necessary. Add two alternatives for positive and negative scientific outcomes without plugging in fabricated numbers. Verify spoken length by word count (roughly 230–300 English words for two minutes depending on delivery).
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Draft a clear, human presentation for a student audience. Use existing evidence only and create a 2-minute pitch with a simple visual cue per section and a Hinglish explanation. Build two truthful endings: candidate results if validated, or a rigorous feasibility audit if not. Rehearse the sequence against actual files in the workspace and call out any missing demo screenshots before presentation.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

<a id="agent-stage-23"></a>
### Agent Stage 23 — Audit references, citations, final claims and release readiness

**Guide section:** [23](#23-references-dataset-citations-and-verification-log)

**Prompt for Codex (copy all):**

```text
Read Section 23, the entire guide and all generated reports. Create `docs/evidence_registry.csv` with claim, source URL/DOI/local artifact, evidence type, verification date, and status (VERIFIED/PROVISIONAL/UNSUPPORTED). Check each dataset accession URL, paper DOI formatting and internal document links where network access is permitted; otherwise mark external verification pending. Write `docs/final_claim_audit.md` identifying every clinical, causal, performance or sample-count claim that exceeds evidence. Confirm results/figures derive from real output artifacts, not illustrative data; do not manufacture citations. Produce a final release checklist and a reproducibility declaration accurately scoped to what was actually executed.
Before finishing, list actual commands run, test outcomes, files changed, input/output schemas, scientific caveats and whether this stage is PASS, BLOCKED or NOT RUN. Do not advance to another stage unless I explicitly ask.
```

**Prompt for Google Antigravity (copy all):**

```text
Run a final evidence review in the editor. Inspect the source table, citation list, results section and all figures; distinguish GEO facts, published findings, experimental outputs and unverified plans. Show an evidence table and request corrections for any unsupported biomarker, AUC, early-diagnosis or mechanistic claim. Mark the project READY only if scripts, reproducibility proof, mappings and output labels satisfy their gates; otherwise return a precise remaining-issues list.
End with a concise evidence/artifact review: created paths, how to execute, tests truly run, assumptions, scientific decision gate and the next stage; do not launch another stage automatically.
```

### Recovery prompt — when an agent hits an error (use instead of vague “fix it”)

```text
Stop changing unrelated files. The current stage has failed. First reproduce the exact error using the command shown in the log. Print the full traceback, input file paths and existence checks, first 5 rows or header metadata, data shapes, column dtypes and relevant package versions; avoid printing entire large files. Distinguish a file-path issue, sample-ID mismatch, platform ID mismatch, data-format/count-scale issue, insufficient unique lncRNA annotation, library API/version problem or statistical design problem. Propose the smallest scientifically sound fix and explain why it is correct. Add a regression test reproducing the failure. Do not replace missing biological information with guesses; do not bypass probe-specificity, validation or leakage gates. After the fix, rerun only the failed command and its tests and report PASS/BLOCKED with evidence.
```

### Handoff prompt — switching between Codex and Antigravity

```text
You are taking over an existing scientific repository. Read PROJECT_GUIDE.md, AGENTS.md, docs/dataset_audit.md, docs/probe_mapping_audit.md and the latest stage report (where present). Run `git status` and inspect recent files without editing. List stages completed, tests actually run, PASS/BLOCKED gate statuses, unresolved scientific questions, input paths and the exact safest next task. Do not assume prior agent output is correct; verify manifests and tests. Do not overwrite existing work, redownload GEO files or open final GSE140829 test labels. Ask me for the next stage number after presenting your handoff summary.
```

### Required final deliverables from agent execution

- `docs/dataset_audit.md` and `outputs/tables/local_inventory.csv` based on actual local files.
- Probe mapping audit with **verified vs provisional** status and real feature counts; a genuine validated map only if specificity was demonstrated.
- Clear, runnable stage scripts + automated parser/schema/leakage tests.
- QC plots, valid differential-expression and pathway/co-expression results **where supported**; block or explain otherwise.
- Internal nested-CV report and, **only with explicit approval and compatible mapped features**, an external-validation report.
- Final `README.md`, honest report/presentation, optional MATLAB demo, optional research-only Streamlit explorer, and a reproducibility report.

**Final rule:** Strong prompts reduce misunderstandings, but no coding agent can guarantee error-free bioinformatics. Real file inspection, gate approvals, independent validation and human review remain essential.
