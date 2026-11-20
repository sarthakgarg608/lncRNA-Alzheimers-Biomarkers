# Stage 10 scientific report · lncRNA annotation and probe verification

**The three blood cohorts share 349 identical, transcriptome-screened probe sequences measuring 342 lncRNA genes.** This exceeds the project's approximate 20-gene feasibility heuristic. The separate brain RNA-seq analysis yields **4,811 lncRNA genes**. Both results support further exploratory research. They establish feature identity within the stated references. Association with Alzheimer's disease, quantitative assay transfer and diagnostic performance remain unevaluated.

## Sources and preflight corrections

All six annotation sources were present locally. No reference download or package installation was required.

| Source | Records | Use |
|---|---:|---|
| GPL10558-50081.txt | 48,107 | Discovery ILMN probe ID, sequence, accession and legacy gene annotations |
| GPL6947-13512.txt | 49,576 | Replication ILMN probe ID, sequence, nuID and legacy annotations |
| GPL15988.txt | 47,322 | External nuID, GenBank accession, gene symbol and description |
| gencode.v50.long_noncoding_RNAs.gtf.gz | 965,508 | GRCh38, GENCODE v50 / Ensembl 116. Gene records yielded 36,476 lncRNA and 1,084 TEC |
| gencode.v50.transcripts.fa.gz | 670,670 | Version-matched sequences for 86,411 genes, including coding, pseudogene and other competing transcripts |
| Human.GRCh38.p13.annot.tsv.gz | 39,376 | Entrez GeneID, GeneType, status, description and explicit Ensembl cross-references |

The platform tables have unique primary IDs. The NCBI table has 27,077 populated Ensembl links and 203 rows sharing an Ensembl target with another row. Version and assembly details are limited to what the local files and registry support. All six files have recorded SHA256 digests. [Source manifest with actual schemas and missingness](../outputs/tables/stage10/annotation_sources.json)

Stage 09's 202 flag records refer to 178 distinct samples. The original GSE125583 series matrix supplies 288 numeric ages in 289 samples. Its age values and missing entry matched the Stage 08 metadata by GSM. Controls averaged 87.0 years and AD samples 83.4 years. The Stage 09 ANOVA p-value of 0.000737 reproduced. The previously stated direction was reversed. In GSE140829, batch explains 35.96% and 42.72% of PC1 and PC2 score variance respectively. Batch–diagnosis confounding remains untested. The historical reports carry correction notices. [Preflight counts, source fields and test implementation](stage10_preflight_audit.md)

## Blood mapping and confidence

Every measured probe joined its platform table. Each distinct measured 50-mer was scanned on both strands against the local transcript FASTA. A validated assignment requires an exact forward match to a GENCODE lncRNA gene, no other gene within two substitutions across the full probe and no identified annotation conflict. This operational specificity rule was selected for this analysis. The hit ledger contains 81,593 sequence–gene records with supporting transcript, orientation and distance. Positive evidence for coding and other RNA targets was checked independently of the lncRNA-only GTF. [Reference and mapping methodology](lncrna_reference_methodology.md)

| Cohort | Measured probes | Any lncRNA candidate | Validated probes | Validated genes | Provisional | Ambiguous | Non-lncRNA | Unmapped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| GSE63061 discovery | 32,049 | 5,221 | 605 | 572 | 1,598 | 8,508 | 18,054 | 3,284 |
| GSE63060 replication | 38,323 | 5,990 | 1,125 | 1,087 | 2,683 | 7,960 | 20,020 | 6,535 |
| GSE140829 external | 47,008 | 8,088 | 940 | 845 | 2,166 | 12,562 | 25,436 | 5,904 |

The five probe statuses partition each measured axis. Candidate counts span those statuses and count each probe once. Among discovery candidates, 374 were provisional and 4,242 ambiguous. Corresponding counts were 1,146 and 3,719 in replication, and 507 and 6,641 in the external cohort. Across all probes, provisional includes unresolved unique historical assignments or incomplete sequence evidence. Ambiguous includes plausible alternative genes and conflicts between modern and legacy annotations. Non-lncRNA requires a unique exact assignment with positive evidence for another class. Unmapped means that the evidence did not identify a modern target. The detailed mapping rows retain original identifiers, sequences, linked genes, accessions, confidence and conflict indicators.

All accepted array probes have complete, nonconstant Stage 09 intensities. Multiple probes targeting the same gene remain separate. Median Pearson correlations across within-gene probe pairs were 0.095 (37 pairs) in discovery, 0.017 (39) in replication and 0.048 (110) externally. These values give little support for default averaging. The retained measurement unit is the original probe. [Probe concordance](../outputs/tables/stage10/within_gene_probe_concordance.csv)

## Compatibility across cohorts

All 47,008 measured GPL15988 nuIDs decoded with valid checksums to 50-base sequences. An independent check on GPL6947 reproduced the supplied sequence in all 48,802 platform rows containing both fields. Of the external measured probes, 31,868 have an identical sequence among all discovery probes. Comparing candidate annotations alone gives 5,607 potentially shared genes. The evidence contracts substantially when both measurements must qualify as unique lncRNA probes. [Bioconductor lumi manual](https://www.bioconductor.org/packages/release/bioc/manuals/lumi/man/lumi.pdf) · [External probe-level crosswalk](../outputs/tables/stage10/external_feature_compatibility.csv)

| Comparison | Shared validated genes | Genes with identical validated probes | Distinct identical sequence–gene pairs |
|---|---:|---:|---:|
| Discovery and replication | 357 | 348 | 355 |
| Discovery and external | 561 | 560 | 592 |
| Replication and external | 465 | 454 | 466 |
| All three blood cohorts | 351 | 342 | 349 |

For external planning, 592 discovery–external pairs have comparable probe identity. Actual external testing requires the later frozen-model and transfer conditions at Gate E. No external diagnoses were accessed, and neither assay calibration nor model performance was tested. [All single- and multi-cohort overlaps](../outputs/tables/stage10/lncrna_overlap_summary.csv)

## Brain RNA-seq

All 39,376 count-matrix GeneIDs joined the NCBI table. Linking used the supplied Entrez-to-Ensembl field. Original gene names and NCBI annotations remain visible for release comparisons. The lncRNA count matrix retains Entrez IDs, including genes without a current GENCODE cross-reference.

| Measured gene category | Genes | Evidence |
|---|---:|---|
| lncRNA | 4,811 | 4,046 active ncRNA records linked to GENCODE lncRNA and 765 with explicit NCBI long RNA descriptions |
| Protein-coding | 19,346 | Active coding assignment after conflict and crosswalk checks |
| Other | 4,086 | 1,913 specific short RNA descriptions and 2,173 other specific NCBI classes |
| Unresolved | 11,133 | Includes 8,598 generic ncRNAs, 1,610 inactive records and 203 rows with shared crosswalk targets |

Of the 4,811 brain lncRNA genes, 4,693 had nonconstant counts. All 4,811 are retained. The brain lncRNA set includes 4,093 distinct supplied Ensembl links, of which 4,046 have concordant v50 lncRNA annotations. Discovery shares 262 validated gene IDs with the brain subset, replication 394 and external 382. The four-way intersection contains 136 genes. These figures describe tissue coverage only. [Discrepancy ledger, including unresolved records and changed symbols](../outputs/tables/stage10/annotation_discrepancies.csv)

## Derived data and figures

The three blood matrices contain 388 × 605, 329 × 1,125 and 587 × 940 original processed intensities respectively. The two brain matrices contain 289 × 4,811 and 289 × 19,346 original integer counts. Rows retain sample IDs and source order throughout. No clinical columns enter the matrices. The inventory contains 7,481 cohort-specific lncRNA features. [Validated feature inventory](../outputs/tables/stage10/validated_feature_inventory.csv)

Six figures show observed data. The GTF figure identifies TEC records excluded from the reference. The confidence, ambiguity and candidate figures show the difference between a candidate and an attributable measurement. The overlap figure shows the smaller subset surviving identical-sequence matching. The brain figure shows the size of the unresolved category alongside usable count subsets. [Stage 10 figures](../outputs/figures/stage10)

## Reproduction and verification

From the repository root in PowerShell, using the existing Python 3.12.10 environment and g++ 6.3.0 on PATH:

```powershell
.venv\Scripts\python.exe scripts\06_build_lncrna_reference.py
.venv\Scripts\python.exe scripts\07_annotate_blood_probes.py
.venv\Scripts\python.exe scripts\08_annotate_brain_genes.py
.venv\Scripts\python.exe -m scripts.stage10_summary
.venv\Scripts\python.exe -m pytest -v
.venv\Scripts\python.exe -m compileall scripts tests -q
.venv\Scripts\python.exe -m pip check
```

The four analysis entry points ran against real local inputs. The full suite passed **343 tests**, including 33 new Stage 10 cases. Tests cover GTF parsing, host genes and pseudogenes, duplicate probes, namespaces, checksum failure, exact and near matches, reverse complements and confidence assignment. Integration checks compared every blood output cell and sample ID against Stages 08 and 09 and every brain output count against the raw TSV. The raw input hashes match the existing inventory. The external accession smoke test now stops at the series accession before sample characteristics. Compilation succeeded and pip reported no broken requirements. [Full pytest log](../outputs/tables/stage10/pytest.log)

## Files created or modified

- **Implementation** · `scripts/06_build_lncrna_reference.py`, `scripts/07_annotate_blood_probes.py`, `scripts/08_annotate_brain_genes.py`, `scripts/stage10_summary.py`, `scripts/stage10_match.cpp`, `scripts/utils/stage10.py`.
- **Tests** · `tests/test_stage10_gtf.py`, `tests/test_stage10_probe_mapping.py`, `tests/test_stage10_brain_annotation.py`, `tests/test_stage10_crosswalk.py`, `tests/test_stage10_integrity.py`. Updated `tests/test_data_io.py` to protect holdout outcomes during full-suite verification.
- **Documents** · This report, `docs/stage10_preflight_audit.md`, `docs/lncrna_reference_methodology.md`, `docs/stage10_gate_b_assessment.md`. Updated `docs/stage08_data_audit.md` and `docs/stage09_qc_report.md` with traceable corrections.
- **data/annotations/** · `lncrna_reference.parquet`, `transcript_gene_catalog.parquet`, `GSE63061_probe_mapping.csv`, `GSE63060_probe_mapping.csv`, `GSE140829_probe_mapping.csv`, `GSE125583_gene_biotypes.csv`.
- **data/processed/stage10/** · `GSE63061_validated_lncrna.parquet`, `GSE63060_validated_lncrna.parquet`, `GSE140829_validated_lncrna.parquet`, `GSE125583_lncrna_raw_counts.parquet`, `GSE125583_protein_coding_raw_counts.parquet`.
- **outputs/tables/stage10/** · `annotation_sources.json`, `gencode_biotype_summary.csv`, `sequence_hits.tsv`, `GSE63061_mapping_summary.csv`, `GSE63060_mapping_summary.csv`, `GSE140829_mapping_summary.csv`, `GSE125583_lncrna_summary.csv`, `blood_coverage.csv`, `annotation_discrepancies.csv`, `external_feature_compatibility.csv`, `lncrna_overlap_summary.csv`, `validated_feature_inventory.csv`, `within_gene_probe_concordance.csv`, `pytest.log`.
- **outputs/figures/stage10/** · `gencode_biotypes.png`, `mapping_confidence.png`, `candidate_vs_validated.png`, `shared_blood_genes.png`, `brain_biotypes.png`, `ambiguity.png`.

## Gate decision and next research stage

B1 and B2 pass for discovery annotation and replicated probe identity. B3 passes for the documented sequence-matched subset and remains conditional for quantitative external testing. B4 passes for the supported brain subset. The two-substitution transcriptome screen does not measure gapped, partial or genomic off-target hybridization, assay background, or reproducibility between laboratories. Detection above background and reliable quantitative transfer require separate evidence. Weak within-gene concordance reinforces the choice to retain probe IDs. [Gate B assessment](stage10_gate_b_assessment.md)

The recommended next stage is a separately authorized Gate C review of assay detectability, sample independence and covariate design in discovery and replication, with brain counts analyzed separately. Stage 10 ended after annotation and verification. External outcomes remain sealed.

**CONDITIONAL PASS:** The shared blood probe subset and the brain lncRNA count subset support exploratory analysis. Unresolved annotations and quantitative external transfer remain open.

