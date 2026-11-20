# Stage 12 · lncRNA co-expression in GSE125583 fusiform gyrus

**PASS for the prespecified within-cohort analysis.** After adjustment for diagnosis, age and sex, ZFAS1 expression tracked multiple ribosomal protein genes, including RPL34 (partial r = 0.921). LY86-AS1 tracked SYN2 (r = 0.926), and PCED1B-AS1 tracked CD68 (r = 0.872). Their strong positive partner lists were enriched for cytoplasmic translation, chemical synaptic transmission and immune response, respectively. These associations persisted after the severe-QC exclusions. The samples are bulk postmortem tissue: differences in neuronal, glial or immune-cell proportions could contribute to all three patterns.

## Inputs and eligibility

The raw NCBI GRCh38.p13 count matrix contained 39,376 Entrez genes and 289 GSMs. The Stage 10 matrices had the same 289 GSMs in the same order, with distinct case codes, unique and disjoint Entrez axes, no missing counts, and exact agreement with every corresponding nonnegative integer in the original TSV. Eligibility used observed counts, independent of Stage 11 significance. [Input SHA256 manifest](../outputs/tables/stage12/sources.csv) · [Prespecified plan](stage12_analysis_plan.md)

| Stage 10 category | Input genes | Below 10 counts in ≥10 samples¹ | Constant raw genes² | Tested genes |
|---|---:|---:|---:|---:|
| lncRNA | 4,811 | 2,653 | 119 | 2,158 |
| Protein-coding | 19,346 | 2,547 | 185 | 16,799 |

¹ The exclusion count is the number of genes with fewer than 10 primary samples reaching 10 counts. ² Constant genes overlap the low-count exclusions. All retained genes had estimable residual variation in both fits. The coding background for GO comprised all 16,799 tested coding genes. [Per-gene audit](../outputs/tables/stage12/feature_audit.csv)

The primary complete cases matched Stage 11 exactly: 218 AD and 70 controls. GSM3577605 lacked age and sex. Primary analysis retained all 52 QC-flagged complete cases. The sensitivity fit excluded six AD samples and one control with severe concurrent low-depth and elevated-zero-fraction flags, leaving 212 AD and 69 controls. All sample decisions are recorded. Raw GEO characteristics identified tissue, case ID, diagnosis, age, sex and Braak score. Neither the source characteristics nor the saved model metadata supplied measured cell-type fractions. Batch, RIN and postmortem interval were also unavailable for adjustment. [Sample ledger](../outputs/tables/stage12/sample_audit.csv)

## Normalization check

For each fit, whole-transcriptome positive-count median-ratio factors were estimated from the original counts of all measured genes in that fit. Positive entries contributed log ratios against the gene reference, and sample factors were centered to geometric mean one. Reapplying this implementation to the Stage 11 lncRNA-only gene set reproduced the saved Stage 11 factors in both fits. Correlations used log2(1 + raw count / whole-transcriptome factor).

| Fit | Whole-transcriptome factor range | Spearman correlation with Stage 11 factors | Median log2(whole / lnc-only), AD / Control | Median log2 ratio, flagged / unflagged |
|---|---|---:|---:|---:|
| Primary | 0.209–3.451 | 0.982 | +0.0079 / −0.0161 | −0.0212 / +0.0074 |
| Severe-QC sensitivity | 0.352–3.364 | 0.981 | +0.0085 / −0.0185 | −0.0138 / +0.0058 |

Individual primary log2 ratios ranged from −0.341 to +0.268. In the joint descriptive HC3 regression, the adjusted AD coefficient was +0.012 log2 units (95% CI −0.010 to +0.034) and the any-QC-flag coefficient was −0.031 (−0.065 to +0.002). After severe-QC exclusions the corresponding estimates were +0.014 and −0.016, both with intervals spanning zero. There was little evidence of a systematic diagnosis-related shift. Flagged samples tended toward smaller whole-transcriptome factors relative to the lncRNA-only reference, with uncertainty about the size of the shift. This comparison evaluates scaling and leaves the original Stage 11 DE results intact. [Factor coefficients and intervals](../outputs/tables/stage12/factor_comparison.csv) · [Normalization figure](../outputs/figures/stage12/normalization_comparison.png)

## Pairwise associations and QC stability

A four-column full-rank design projected out the intercept, diagnosis, centered age and sex. Pearson correlation of the two residual vectors estimates the association conditional on those covariates. The two-sided t test used n − rank(design) − 1 degrees of freedom: 283 primary and 276 sensitivity. Blocks of 32 lncRNAs controlled memory use. One global BH correction per fit included all valid pairwise p-values across all blocks and both signs.

| Pairs at global BH q < 0.05 | Positive | Negative | Total | Distinct lncRNAs / coding genes |
|---|---:|---:|---:|---:|
| Primary | 11,337,570 | 11,372,608 | 22,710,178 | 2,158 / 16,799 |
| Sensitivity | 11,189,330 | 11,311,579 | 22,500,909 | 2,158 / 16,799 |
| Significant in both, same sign | 10,829,305 | 10,903,195 | 21,732,500 | 2,158 / 16,799 |

The denominator was 36,252,242 pairs per fit. Of primary discoveries, 95.69% retained significance and direction. A further 977,675 lost significance, 768,406 gained significance, and three passed both thresholds with opposite signs. Including sensitivity pairs below threshold, 470 primary significant pairs changed sign. The sensitivity samples are a subset of the primary cohort, so this measures internal stability. The dense result (62.64% of primary pairs significant) is consistent with broad shared variation. The number of discoveries alone offers little specificity for an individual biological mechanism. [Per-lncRNA counts](../outputs/tables/stage12/edges_per_lncrna.csv)

| lncRNA | Coding partner | Primary partial r | Sensitivity partial r |
|---|---|---:|---:|
| ZFAS1 | RPL34 | +0.921 | +0.912 |
| LY86-AS1 | SYN2 | +0.926 | +0.927 |
| PCED1B-AS1 | CD68 | +0.872 | +0.878 |
| LOC100652999 | GADD45GIP1 | −0.929 | −0.924 |
| NPTN-IT1 | KRTCAP2 | −0.927 | −0.924 |
| SLX1B-SULT1A4 | SULT1A4 | +0.992 | +0.992 |

The last pair illustrates an important limitation of the network figure. NCBI describes SLX1B-SULT1A4 as a readthrough and places it at 29,455,091–29,464,980 on NC_000016.10. SULT1A4 spans 29,459,886–29,464,980 on the same reference. Shared sequence or ambiguous read assignment can contribute to their near-perfect association. Their individual contributions cannot be distinguished with the gene-count matrix. The figure displays eight selected associations; the readable summary contains up to three positive and three negative partners per lncRNA, 12,948 rows in total. [Network figure](../outputs/figures/stage12/coexpression_network.png) · [Readable edge summary](../outputs/tables/stage12/edge_summary.csv)

## Local GO Biological Process enrichment

The GO source was the project’s NCBI `Human.GRCh38.p13.annot.tsv.gz` snapshot, SHA256 `a0b6ff49b0851334c709edbd52ae9c80ae998d66d9a30185de364ffb4a59efea`. Its supplied GOProcessID and GOProcess fields aligned directly to Entrez GeneID and agreed with the Stage 10 annotation. The checksum identifies the precise mapping snapshot; the local file gives no separate GO ontology release date. Only supplied Biological Process memberships were used. No external calls or downloads were needed.

Strong foregrounds required q < 0.05, consistent direction and |r| ≥ 0.75 in **both** fits. Across 925 nonempty lncRNA/direction lists there were 88,089 memberships from 7,206 distinct coding genes. Of 16,799 eligible background genes, 15,074 had supplied BP annotations and 1,725 lacked them. All 16,799 remained in the denominator. Of 12,239 observed terms, 1,826 had 15–500 measured members. The 560 lists containing at least ten genes had 86,984 memberships: 78,571 mapped to any BP term and 74,197 to at least one size-eligible term. [Measured background](../outputs/tables/stage12/enrichment_background.csv) · [Per-list mapping coverage](../outputs/tables/stage12/foreground_coverage.csv)

One-sided hypergeometric tests covered all 560 × 1,826 list–term combinations, including zero overlaps. Global BH across **1,022,560 tests** yielded **3,258 significant list–term associations** covering 357 distinct terms and 278 lncRNAs. Positive lists contributed 2,833 associations and negative lists 425. Multiple terms can share genes, and different lncRNAs can share partners.

| lncRNA partner list | GO Biological Process | Overlap / foreground | Background term size | Fold enrichment | Global q |
|---|---|---:|---:|---:|---:|
| ZFAS1 · positive | Cytoplasmic translation (GO:0002181) | 77 / 176 | 87 | 84.48 | 4.38 × 10⁻¹⁴³ |
| LY86-AS1 · positive | Chemical synaptic transmission (GO:0007268) | 70 / 725 | 221 | 7.34 | 4.07 × 10⁻³⁷ |
| GIHCG · positive | Mitochondrial translation (GO:0032543) | 28 / 283 | 95 | 17.50 | 3.11 × 10⁻²³ |
| SMARCA5-AS1 · negative | Proton motive force-driven mitochondrial ATP synthesis (GO:0042776) | 17 / 113 | 63 | 40.12 | 3.30 × 10⁻¹⁹ |
| PCED1B-AS1 · positive | Immune response (GO:0006955) | 24 / 141 | 251 | 11.39 | 7.28 × 10⁻¹⁵ |

These tests ask whether a process is overrepresented among a particular lncRNA’s co-expressed coding partners. They do not test pathway activity differences between AD and controls. Dependence among genes, expression-dependent detectability, uneven GO coverage and missing cell fractions limit interpretation of enrichment magnitudes and nominal probabilities. [Significant GO results](../outputs/tables/stage12/GO_BP_significant.csv.gz)

## Artifacts and reproduction

- **Implementation and documentation** — `scripts/13_brain_coexpression.py`, `scripts/utils/stage12_go.py`, `tests/test_stage12_coexpression.py`, `docs/stage12_analysis_plan.md` and this report.
- **Inputs and normalization** — `sources.csv`, `sample_audit.csv`, `feature_audit.csv`, `size_factors.csv`, `factor_comparison.csv`.
- **Associations** — `all_pairs.parquet`, `edges_per_lncrna.csv`, `edge_summary.csv`, `network_edges.csv`. The 1.74 GB Zstandard Parquet has one row per Entrez pair, with primary and sensitivity r, two-sided p and global BH q. Records are grouped by lncRNA so selected columns and row groups can be read without loading the entire file.
- **Enrichment and verification** — `enrichment_foreground.csv.gz`, `enrichment_background.csv`, `foreground_coverage.csv`, `GO_BP_all_tests.parquet`, `GO_BP_significant.csv.gz`, `run_summary.json` and `pytest.log`. Tables are under `outputs/tables/stage12/`. Both PNGs are under `outputs/figures/stage12/`.

Run from the repository root in PowerShell with the existing Python 3.12.10 environment:

```powershell
.venv\Scripts\python.exe scripts\13_brain_coexpression.py
.venv\Scripts\python.exe -m pytest -v
.venv\Scripts\python.exe -m compileall scripts tests -q
.venv\Scripts\python.exe -m pip check
```

**367 tests passed**, including six Stage 12 tests. Compilation succeeded and pip found no broken requirements. The independent GO check reproduced the ZFAS1 overlap and hypergeometric p-value from the raw annotation after explicitly reading the CSV Entrez keys as strings.

All 90 Parquet row groups were read back. All saved correlations and probabilities were finite and within bounds, and the positive and negative significant counts reconciled to the per-lncRNA summaries in both fits. The two figures were visually inspected. The tests cover sample alignment, missing and altered counts, constant genes, covariate adjustment, reference regression, global BH across block boundaries, real size factors, real pair statistics and the GO enrichment denominator. [Test log](../outputs/tables/stage12/pytest.log)

**Verdict: PASS.** The computational objectives were met with local reference data. Interpretation remains conditional on bulk tissue composition and gene-level read assignment. A useful follow-up would examine the strongest distinct-locus associations in cell-type-resolved or independent brain data. GSE140829 outcomes remained sealed. Work stopped at Stage 12.

