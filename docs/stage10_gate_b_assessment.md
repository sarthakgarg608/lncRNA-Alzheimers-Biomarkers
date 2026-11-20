# Stage 10 Gate B assessment

| Gate | Decision | Evidence and scope |
|---|---|---|
| B1 · Blood discovery | PASS for computational probe annotation | GSE63061 has 605 sequence-specific, nonconstant probes for 572 lncRNA genes. All 388 samples remain. |
| B2 · Related-cohort replication | PASS for matched-probe coverage | GPL6947 measures 1,087 validated genes. Discovery and replication share 357 genes, with identical validated 50-mers for 348 genes (355 distinct sequence–gene pairs). |
| B3 · External technical compatibility | PASS for matched probe identity, conditional for future testing | GPL15988 measures 845 validated genes. Of 561 shared with discovery, 560 have an identical validated sequence (592 pairs). The nuID checksum and independent GPL6947 sequence check support decoding. Cross-cohort intensity calibration and a frozen analysis remain outstanding. |
| B4 · Brain RNA-seq | PASS for the annotated count subset | Of 39,376 Entrez features, 4,811 qualify as lncRNA, including 4,046 with concordant NCBI and GENCODE evidence. There are 19,346 protein-coding features. Unresolved genes are excluded from both derived matrices. |

Across the three blood cohorts, 351 validated genes overlap and 342 have at least one identical validated 50-mer across all three. Those 349 distinct sequence–gene pairs exceed the predeclared planning heuristic of approximately 20 shared genes. All validated array probes have complete, nonconstant values in their respective cohorts. Above-background detection and analytical reproducibility still require independent evidence. For genes with multiple probes, the median within-gene pairwise correlation is 0.095 in discovery, 0.017 in replication and 0.048 in the external cohort. Retain probe-level measurements.

These gates address computational annotation within the local GENCODE v50 transcriptome and the documented two-substitution screen. External outcomes remain sealed. The 136 Ensembl genes present in all three blood inventories and the brain lncRNA subset provide biological context across tissues. A reasonable next research stage is a separately authorized Gate C design and assay-quality review of the discovery and replication data, with an independent brain analysis. The external cohort remains reserved for Gate E.

**CONDITIONAL PASS** for Stage 10. The matched blood subset and brain count subset are usable. Quantitative cross-platform transfer and the unresolved features remain open.

