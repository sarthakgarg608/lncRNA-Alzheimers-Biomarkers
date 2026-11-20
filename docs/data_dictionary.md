# Data Dictionary & Bioinformatics Schemas

**Project:** Establishing the Role of lncRNAs as Biomarkers in Alzheimer's Disease  
**Course:** Computational Genomics and Proteomics (OE3E40)  
**Target Audience:** Undergraduate / B.Tech researchers in Computational Genomics  
**Version:** 1.0 (Stage 07 Deliverable)  

---

## 1. Introduction: What Are We Measuring?

To understand biomarker discovery, we must first be clear about biological entities and measurements:

- **Sample:** Biological material drawn from a specific human subject at a point in time (e.g., peripheral whole blood from patient `GSM1539409`). Each sample is assigned an accession ID by NCBI Gene Expression Omnibus (GEO).
- **Gene:** A functional physical unit of heredity. In transcriptomics, genes are transcribed into RNA. Some RNAs code for proteins (mRNA), while others do not (non-coding RNAs).
- **lncRNA (Long Non-Coding RNA):** An RNA molecule longer than 200 nucleotides that is not translated into a functional protein. Many lncRNAs act as epigenetic and transcriptional regulators. In Alzheimer's disease (AD), specific lncRNAs (such as *BACE1-AS*, *NEAT1*, and *MALAT1*) are known to regulate amyloid precursor processing, tau phosphorylation, or neuroinflammatory responses.
- **Microarray Probe:** A short, synthetic single-stranded oligonucleotide (usually 50 nucleotides on Illumina HumanHT-12 arrays) fixed onto a beadchip surface. If an RNA molecule in the patient sample contains a complementary sequence, it hybridizes to the probe, producing a fluorescent signal.

---

## 2. Standard Sample Metadata Schema

Every sample across our discovery, replication, and validation cohorts is standardized into the following schema:

| Column Name | Data Type | Permitted Values / Format | Description & Scientific Role |
|---|---|---|---|
| `sample_id` | `string` | Unique `GSM...` string | NCBI GEO sample accession. Links expression columns to donor phenotypes. |
| `dataset_id` | `string` | `GSE63061`, `GSE63060`, `GSE140829`, `GSE125583` | Cohort identifier. Essential for accounting for batch effects and study design. |
| `diagnosis` | `string` | `AD`, `Control`, `MCI`, `Unknown` | Standardized clinical diagnosis classification. |
| `original_label` | `string` | As reported in GEO | Exact raw label from GEO (e.g., `"Alzheimer's disease"`, `"neurologically normal"`). Kept for auditability. |
| `age` | `float` / `null` | Numeric years (or `null`) | Patient age at blood draw or death. A critical confounder: AD risk increases exponentially with age. |
| `sex` | `string` / `null` | `Female`, `Male`, `Unknown`, `null` | Biological sex. Critical confounder: lncRNA *XIST* is expressed only in females, so unadjusted analysis could rediscover sex rather than AD. |
| `tissue` | `string` | `whole_blood`, `fusiform_gyrus_brain` | Anatomical source of RNA. Blood and brain have completely different baseline cell compositions. |
| `brain_region` | `string` / `null` | `fusiform_gyrus` (or `null` for blood) | Specific brain region sampled for post-mortem tissue in `GSE125583`. |
| `batch` | `string` / `null` | Laboratory plate / site identifier | Technical processing run or hybridisation batch. Needed for batch-effect modeling. |
| `platform_id` | `string` | `GPL10558`, `GPL6947`, `GPL15988`, `GPL16791` | GEO platform ID denoting array technology or sequencer model. |

---

## 3. Why AD, MCI, and Controls Must Remain Strictly Separate

One of the most dangerous errors in biomedical data science is **label leakage and forced binary classification**:

1. **MCI (Mild Cognitive Impairment) is a transitional, heterogeneous state:**
   - Some patients with MCI convert to Alzheimer's dementia within 2–3 years (amnestic MCI converters).
   - Other MCI patients remain stable for a decade, or their cognitive issues stem from vascular disease, depression, or metabolic disorders.
   - **Rule:** Never re-label MCI patients as "Control" or "AD". Doing so pollutes the baseline and destroys the biological specificity of the biomarker signal.

2. **Unknown / Borderline Cases:**
   - Any sample where the original record indicates doubtful or mixed pathology must remain classified as `Unknown` and excluded from supervised training sets.

3. **External Hold-Out Locking (GSE140829):**
   - The diagnosis labels of `GSE140829` must remain **sealed** while building and selecting the lncRNA panel. If external labels are used to select features, the resulting Area Under the ROC Curve (AUC) is artificially inflated (optimism bias) and cannot replicate in actual clinical practice.

---

## 4. What Is an Expression Matrix?

An expression matrix is a numerical two-dimensional table representing RNA abundances:

### Interim / Storage Convention: Features × Samples
```text
                  GSM1539409   GSM1539410   GSM1539411   GSM1539412
ILMN_1343048         12.45        12.31        12.60        12.18
ILMN_1343049          8.12         8.05         7.98         8.22
...
```
- **Rows:** Genomic features (microarray probe IDs or gene identifiers).
- **Columns:** Individual biological samples (`GSM...` accession IDs).
- **Values:** Measured abundance.

### Machine-Learning Convention: Samples × Features
When training models (LASSO, Random Forests, SVM), the matrix is transposed so that:
- **Rows:** Observations (samples).
- **Columns:** Features (candidate lncRNA biomarkers).
- **Target Vector $y$:** Diagnosis (`AD` vs. `Control`).

---

## 5. Microarray Intensities vs. RNA-Seq Counts: Fundamental Differences

A common catastrophic mistake in computational biology is passing microarray data into RNA-seq software like DESeq2. The physics and mathematics are completely different:

| Property | Illumina Microarray (`GSE63061`, `GSE63060`, `GSE140829`) | RNA-seq (`GSE125583`) |
|---|---|---|
| **Measurement Principle** | Fluorescent hybridization intensity on bead arrays | Direct high-throughput cDNA fragment sequencing |
| **Data Nature** | Continuous, log-transformed intensity values | Discrete non-negative integer read counts ($0, 1, 2, ...$) |
| **Underlying Distribution** | Approximately log-normal / Gaussian after normalization | Negative Binomial (overdispersed Poisson) |
| **Dynamic Range** | Limited by saturation and background fluorescence ($10^2$ to $10^5$) | High dynamic range, limited only by sequencing depth |
| **Appropriate Differential Expression Tools** | `limma` (moderated $t$-statistic) or `statsmodels` OLS | `DESeq2` / `pydeseq2` (Negative Binomial Generalized Linear Models) |
| **Fatal Misuse** | **Passing microarray intensities into DESeq2 fails mathematical assumptions** | Treating raw counts as normalized continuous values without dispersion modeling |

---

## 6. Why Gene Identifier Namespaces and Crosswalks Matter

Biological databases use different identifier namespaces that cannot be merged naively:

1. **Illumina Probe IDs (`ILMN_...` or `nuID`):**
   - Microarray chips measure *probes*, not genes directly.
   - A single gene may have 3 different probes measuring different exons or isoforms.
   - An older probe designed in 2008 may cross-hybridize to multiple transcripts or to an unannotated pseudogene.
   - **Crucial Rule:** Matching a probe's reported symbol (e.g., "NEAT1") in an annotation table is only a *provisional* candidate. We must verify probe sequence uniqueness against GENCODE transcript databases before asserting it is an lncRNA biomarker.

2. **NCBI Entrez IDs vs. Ensembl Gene IDs:**
   - `GSE125583_raw_counts` indexes genes by **NCBI Entrez numeric IDs** (e.g., `100287102`).
   - `gencode.v50.long_noncoding_RNAs.gtf` uses **Ensembl Gene IDs** (e.g., `ENSG00000290825.2`).
   - `Human.GRCh38.p13.annot.tsv.gz` provides the official crosswalk bridge between Entrez and Ensembl.
   - Merging tables by position or assuming numerical IDs equal Ensembl IDs leads to complete data corruption.

---

## 7. What Information Is Missing from Public Datasets?

When interpreting our research findings, we must be honest about gaps in public repositories:

1. **APOE $\varepsilon4$ Genotype:**
   - Apolipoprotein E ($\varepsilon4$ allele) is the strongest genetic risk factor for sporadic AD.
   - While `GSE140829` records `apoe_status` (`E3_E3`, `E3_E4`, etc.), `GSE63061` and `GSE63060` have incomplete public APOE records.
2. **Detailed Cognitive Scores (MMSE, CDR):**
   - Exact mini-mental state examination scores are not uniformly provided across all public sample records.
3. **Blood Cell Subtype Proportions:**
   - Whole blood RNA is a mixture of neutrophils, lymphocytes, monocytes, and platelets. An apparent increase in an lncRNA could reflect an altered neutrophil-to-lymphocyte ratio rather than neurodegeneration. Computational deconvolution or cell-type marker auditing must be considered.
4. **RNA Integrity Number (RIN):**
   - Degraded RNA in post-mortem brain samples or archived blood tubes can distort transcript quantification if not controlled.
