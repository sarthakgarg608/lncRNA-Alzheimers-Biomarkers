# Stage 11 · Differential expression of lncRNAs in Alzheimer’s disease

**Three discovery blood probes showed adjusted AD–Control differences at FDR < 0.05. None met the prespecified replication criterion.** Two lacked an identical sequence-supported probe in the related replication cohort. The third, a MIR29B2CHG probe, had a similar positive standardized effect in both cohorts but a directional replication FDR of 0.109. In a separate analysis of fusiform gyrus, 1,067 of 2,157 valid brain lncRNA gene tests met FDR < 0.05. The brain findings describe bulk postmortem tissue. External GSE140829 outcomes remained sealed.

## Question, samples and annotation

For each eligible probe or gene, the two-sided null was a zero adjusted diagnosis coefficient, with Control as reference. Positive coefficients indicate higher expression in AD. The plan was written before feature-level significance testing. [Statistical analysis plan](stage11_statistical_analysis_plan.md)

| Dataset | Original samples | Primary AD / Control | Stage 10 features | Valid tests | FDR < 0.05 | Exploratory FDR < 0.10 |
|---|---:|---:|---:|---:|---:|---:|
| GSE63061 · discovery blood | 388 | 139 / 134 | 605 probes | 605 | 3 probes, 3 genes represented | 4 probes, 4 genes represented |
| GSE63060 · related replication blood | 329 | 145 / 104 | 1,125 probes | 1,125 | 6 probes, 6 genes represented | 13 probes, 13 genes represented |
| GSE125583 · brain | 289 | 218 / 70 | 4,811 Entrez genes | 2,157 | 1,067 genes | 1,196 genes |

Discovery excluded 109 MCI and six Unknown samples. Replication excluded 80 MCI. The brain series contained 219 AD and 70 controls. Checking the original age characteristic against the GSM-indexed metadata reproduced all 289 entries. GSM3577605, an AD sample, lacked both age and sex and was excluded by the complete-case rule. All other eligible samples, including flagged samples, entered the primary analyses. The three sample ledgers identify every sample, decision and QC reason.

Only Stage 10 `VALIDATED_UNIQUE_LNCRNA` blood probes were eligible. Their sequences had a forward exact match and a single GENCODE gene within two substitutions in both orientations. Their source tables contained no accepted conflicts. Each probe retained its ID. The 605 and 1,125 probes represent 572 and 1,087 distinct genes respectively. No provisional, ambiguous or unmapped probe entered a model. In brain, 4,046 input genes had concordant NCBI/GENCODE lncRNA evidence and 765 qualified by explicit NCBI long-RNA description. Counts use Entrez IDs throughout. Stage 10 evidence classes remain in the DE tables. Before testing, every matrix value, ordered sample and feature ID was checked against its Stage 09 projection and eligible annotation. Matrix and mapping SHA256 digests are in the fit logs.

## Models and design checks

**Blood.** Rscript was unavailable in the accessible environment. The prespecified fallback fit separate statsmodels OLS models for processed probe intensity against diagnosis, centered age in years, sex and collection-center proxy. The three-character center code came from the saved sample description. No additional normalization was applied. Standard errors used the finite-sample-corrected chip-clustered sandwich estimator, with chip inferred from the numeric array-title barcode. Discovery had 34 chip groups and 33 t-reference degrees of freedom. Replication had 29 and 28, respectively. No cross-probe empirical-Bayes variance moderation was used. The coefficients are adjusted differences in processed intensity units. Standardization for matched-probe comparison divides each coefficient by its own model’s residual SD.

The primary designs had full rank, 11 columns in discovery and nine in replication, with condition numbers 109.24 and 38.39 in the specified units. The discovery center proxy had eight levels. OA0 contained one control and no AD samples, and TLS had seven AD and two controls. These small strata contribute little within-center information about diagnosis. Replication had six centers with both groups at each center. Sex counts were 166 female / 107 male in discovery and 161 / 88 in replication. Clustered errors accommodate within-chip correlation; systematic chip-level differences associated with diagnosis can still bias a coefficient. The documented blood sample-description codes were unique within each cohort and had zero literal overlap across cohorts. Their formats differ and do not establish patient independence across AddNeuroMed batches.

**Brain.** PyDESeq2 0.6.0 fit a negative-binomial model to nonnegative integer counts. Unique raw case codes and GSM order were checked. The covariates were centered age and sex, with 130 female and 158 male samples. The four-column design had full rank and condition number 25.27. Batch was unavailable. Braak pathology was excluded from the primary adjustment because conditioning on disease pathology changes the biological question. The saved model metadata lacked RIN, postmortem interval and cell-composition estimates.

A fixed filter required at least 10 counts in at least 10 primary samples. It retained 2,158 genes and left 2,653 labeled low-count rows in the output. Positive-count size factors on these genes ranged from 0.260 to 3.615. The actual fitted dispersion trend was parametric in both primary and sensitivity fits. Primary final dispersions ranged from 0.0141 to 9.291. The unshrunk adjusted log2 fold changes, log2-scale SE and two-sided Wald statistics are recorded. Independent filtering and count replacement were disabled. Gene-wise dispersion, MAP dispersion and LFC convergence were checked individually. Entrez 101928673 failed gene-wise dispersion convergence and was excluded from primary inference. In sensitivity analysis, Entrez 100507178 also failed MAP dispersion convergence. Each failure remains in the diagnostics.

Whole-transcriptome library sizes ranged from 8.06 to 93.56 million counts in primary samples. One gene had a maximum Cook distance exceeding the nominal 99th-percentile F cutoff in each fit. The primary maximum was 6.559 for TTTY14 in GSM3577761. The package’s replicate-aware Cook filter masked zero genes. The diagnostics retain per-gene maxima and the corresponding GSM; this observation warrants caution when interpreting TTTY14. No flagged count was replaced.

**Multiplicity.** Each cohort and each sensitivity fit has its own BH correction over finite, converged, eligible two-sided tests. The directional replication family was specified separately. Blood gene counts describe the distinct stable IDs represented by significant probes. No minimum-p probe selection, gene-level omnibus test or gene-level FDR was performed.

## Blood effects and replication

| Discovery probe | Gene | AD − Control | Clustered SE | Two-sided FDR | Matched replication |
|---|---|---:|---:|---:|---|
| ILMN_3200414 | ENSG00000304921 | −0.25288 | 0.04471 | 0.00161 | No eligible identical sequence |
| ILMN_3297577 | ENSG00000303192 | −0.23211 | 0.04692 | 0.00652 | No eligible identical sequence |
| ILMN_1907042 | MIR29B2CHG | +0.08122 | 0.01748 | 0.01049 | +0.06334, directional FDR 0.10927 |

The GENCODE symbol field contains the stable identifier for the first two genes. In replication’s own two-sided analysis, six probes reached FDR < 0.05: PSMD6-AS2 (+0.12817, FDR 0.000837), ANXA2R-AS1 (+0.11641), ATXN7L3-AS1 (+0.10488), LINC01331 (+0.02689), LINC02762 (+0.06437) and SNHG5 (−0.50195). The latter five had FDR 0.04268. Each value is in its dataset’s processed intensity units. [Discovery results](../outputs/tables/stage11/GSE63061_differential_expression.csv) · [Replication-cohort results](../outputs/tables/stage11/GSE63060_differential_expression.csv)

Across 355 identical-sequence and identical-gene pairs, 198 had consistent signs and 157 had opposite signs. Another 250 discovery probes lacked an eligible match. Standardized effects correlated at Pearson r = 0.350 across matched pairs. Ten pairs passed the replication-directional BH threshold on the replication side, but none also met discovery FDR < 0.05. The complete prespecified replication rule therefore yielded zero supported pairs. For MIR29B2CHG, the matched standardized coefficients were +0.474 and +0.397 residual SD. Its one-sided replication p-value was 0.00616 before correction over the 355-pair family. Another discovery probe for that same gene, ILMN_1747531, had coefficient −0.00368 and FDR 0.995. Across matched pairs, 14 genes had opposing same-gene probe signs in at least one cohort. A single gene-level direction would obscure these observations. [Probe-level replication evidence](../outputs/tables/stage11/blood_replication_summary.csv)

## Brain and exploratory cross-tissue context

Of 1,067 significant brain genes, 539 had higher and 528 had lower adjusted expression in AD. Among concordant NCBI/GENCODE records, 929 of 1,899 valid tests were significant. Among explicit NCBI long-RNA records, 138 of 258 were significant.

| Brain gene | Entrez ID | Adjusted log2 fold change | SE | FDR | Sensitivity log2 fold change |
|---|---:|---:|---:|---:|---:|
| TAGAP-AS1 | 105378083 | −0.545 | 0.059 | 5.26 × 10⁻¹⁷ | −0.548 |
| LINC02937 | 100507103 | +0.709 | 0.078 | 1.12 × 10⁻¹⁶ | +0.700 |
| KDM7A-DT | 100134229 | −0.536 | 0.062 | 3.01 × 10⁻¹⁵ | −0.538 |

These are the three smallest primary brain p-values among valid tests. The number of associations depends on the lncRNA-specific normalization reference, available covariates and composition of postmortem tissue. [Full brain results and annotation evidence](../outputs/tables/stage11/GSE125583_differential_expression.csv)

| Distribution of valid adjusted effects | Minimum | 25th percentile | Median | 75th percentile | Maximum |
|---|---:|---:|---:|---:|---:|
| Discovery · processed intensity | −0.2529 | −0.0043 | +0.0001 | +0.0046 | +0.1325 |
| Replication · processed intensity | −0.5019 | −0.0058 | +0.0007 | +0.0073 | +0.1780 |
| Brain · log2 fold change | −1.7543 | −0.2176 | +0.0220 | +0.2269 | +2.1789 |

An explicit Ensembl crosswalk linked 281 discovery and 411 replication probes to a brain lncRNA. Both tissues had valid tests for 188 discovery probes (175 genes) and 206 replication probes (193 genes). Among discovery–brain pairs, 85 signs agreed and 103 opposed. None were significant in both tissues; 109 were significant only in brain and 79 in neither. Among replication–brain pairs, 102 signs agreed and 104 opposed. Three were significant in both tissues, two only in blood, 110 only in brain and 91 in neither. PSMD6-AS2 and ANXA2R-AS1 were higher in AD in both replication blood and brain. LINC01331 was higher in replication blood and lower in brain. Tissue composition, assay scale and participants differ. All unmatched or untestable rows remain marked as such. [Cross-tissue comparisons](../outputs/tables/stage11/brain_blood_comparison.csv)

## QC sensitivity and figures

All primary eligible samples with Stage 09 flags were retained. Discovery included 47 flagged samples (20 AD and 27 controls); replication included 14 (five AD and nine controls); brain included 52 (39 AD and 13 controls). The sensitivity rule excluded six discovery samples with concurrent extreme median and IQR flags, zero replication samples and seven brain samples with extreme low library size and elevated zero fraction. Stage 09 used a two-sided zero-fraction cutoff despite naming its flag “High zero fraction”; the Stage 11 rule explicitly verified the high direction. The flags identify technical extremes and can also reflect biology. They do not establish assay failure.

| Cohort | Sensitivity AD / Control | Same sign among jointly valid tests | Primary hits retained at FDR < 0.05 and same sign | Lost / gained |
|---|---:|---:|---:|---:|
| Discovery | 135 / 132 | 568 / 605 | 3 / 3 | 0 / 1 |
| Replication | 145 / 104 | 1,125 / 1,125 | 6 / 6 | 0 / 0 |
| Brain | 212 / 69 | 2,126 / 2,156 | 1,046 / 1,067 | 21 / 21 |

The replication sample set remained identical under the sensitivity criterion, so its coefficients are identical by construction. In brain, equal total discovery counts in both fits conceal 42 changes in FDR-threshold membership. Primary and sensitivity results remain side by side in every DE table. [QC sensitivity summary](../outputs/tables/stage11/qc_sensitivity_summary.csv)

The three volcano plots show adjusted effect against −log10 BH FDR. The p-value histograms show a much larger concentration of small p-values in brain. The matched blood effect scatter shows modest agreement, while the cross-tissue plot retains both tissue-specific units. The primary–sensitivity scatter illustrates how far individual coefficients moved. The boxplots display the lowest-FDR feature in each cohort and show unadjusted processed array intensities or visualization-only log2(1 + count / size factor) in brain. Statistical conclusions come from the covariate-adjusted models. [Eight Stage 11 figures](../outputs/figures/stage11)

## Deliverables and verification

The association table contains 1,076 feature–cohort rows (3 + 6 + 1,067). Each retains its evidence, direction, coefficient, raw p-value, FDR, replication context and sensitivity result. These full-cohort, label-informed results cannot be used as a fixed feature set in future cross-validation. Any later supervised selection must be learned within training folds. [Candidate association table](../outputs/tables/stage11/candidate_lncrna_summary.csv)

Files created under the repository root:

- `docs/` — `stage11_statistical_analysis_plan.md`, `stage11_differential_expression_report.md`.
- `scripts/` — `09_differential_expression_blood.py`, `10_differential_expression_brain.py`, `11_replication_analysis.py`, `12_de_visualization.py`, `utils/stage11.py`.
- `tests/` — `test_stage11_blood_de.py`, `test_stage11_brain_de.py`, `test_stage11_replication.py`, `test_stage11_statistics.py`, `test_stage11_integrity.py`.
- `outputs/tables/stage11/` — for each of GSE63061, GSE63060 and GSE125583, `<dataset>_differential_expression.csv` and `<dataset>_samples.csv`. For brain, `GSE125583_primary_diagnostics.csv`, `GSE125583_sensitivity_diagnostics.csv`, `GSE125583_primary_size_factors.csv`, `GSE125583_sensitivity_size_factors.csv`. Shared files are `blood_replication_summary.csv`, `brain_blood_comparison.csv`, `candidate_lncrna_summary.csv`, `qc_sensitivity_summary.csv`, `blood_fit.log`, `brain_fit.log` and `pytest.log`.
- `outputs/figures/stage11/` — `GSE63061_volcano.png`, `GSE63060_volcano.png`, `GSE125583_volcano.png`, `blood_effect_replication.png`, `top_lncrna_expression_boxplots.png`, `brain_blood_comparison.png`, `pvalue_distribution.png` and `qc_sensitivity_comparison.png`.

PowerShell reproduction from the repository root, with the existing Python 3.12.10 environment:

```powershell
.venv\Scripts\python.exe scripts\09_differential_expression_blood.py
.venv\Scripts\python.exe scripts\10_differential_expression_brain.py
.venv\Scripts\python.exe scripts\11_replication_analysis.py
.venv\Scripts\python.exe scripts\12_de_visualization.py
.venv\Scripts\python.exe -m pytest -v
.venv\Scripts\python.exe -m compileall scripts tests -q
.venv\Scripts\python.exe -m pip check
```

All four analyses ran locally. **361 tests passed**, including 18 new Stage 11 cases. Real-data checks covered sample and feature order, matrix provenance, original raw-file checksums, nonnegative integer counts, valid-fit convergence, Wald and clustered-t p-values, BH correction and the external lock. An independent least-squares and cluster-sandwich calculation reproduced one actual coefficient and SE in each blood cohort. The candidate list reconciled exactly to primary FDR-qualified rows. The eight PNGs were rendered and inspected. Compilation passed and pip found no broken requirements. [Full test log](../outputs/tables/stage11/pytest.log)

## Scientific decision

**CONDITIONAL PASS.** Gate C supports reporting adjusted associations for the valid tests and a negative result for the prespecified blood replication criterion. The small blood discovery set, incomplete replication-platform coverage, uncertain cross-cohort patient linkage, unmeasured biological and technical covariates, and brain normalization assumptions limit generalization. Computational probe specificity also leaves assay detectability and quantitative transfer unresolved. No external outcomes or classifier performance were evaluated.

The recommended next stage is a separately authorized review of the robust brain associations and assay quality, with priority given to whole-transcriptome normalization sensitivity, tissue composition and independent brain replication. A blood follow-up would require a better-matched assay and independently documented participants. Stage 11 ends here.

