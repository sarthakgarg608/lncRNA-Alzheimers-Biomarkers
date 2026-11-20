# AGENTS.md — Alzheimer's lncRNA Biomarker Research

These instructions apply to any AI coding agent working on this repository, including OpenAI Codex and Google Antigravity. Read `PROJECT_GUIDE.md` (or the supplied updated Markdown guide) before changes and follow the requested stage **only**.

## Scientific correctness
- This is an exploratory research project, **not** a clinical diagnostic tool. Do not fabricate patient metadata, p-values, biomarkers, citations, plots, or AUC scores. Illustrations must be labeled as simulated, never passed off as results.
- Preserve original GEO data and original sample labels. Never relabel MCI, borderline, mixed, or unknown samples as AD or Control. Do not silently drop records.
- Blood microarrays (GSE63061, GSE63060, GSE140829) have uncertain lncRNA coverage. Gene-symbol matching alone makes **provisional** candidates, not validated lncRNA probes. Confirm unique probe/sequence specificity before calling a feature a lncRNA biomarker or using it in lncRNA-only ML.
- GSE125583 brain RNA-seq cannot be treated as independent **blood** validation. Do not merge tissue-level observations as if comparable patients.
- GENCODE uses Ensembl IDs; the available GSE125583 count matrix may use NCBI/Entrez IDs. Audit identifier namespaces and document crosswalk evidence before joining.
- Never pass processed microarray intensities into DESeq2 as if they were raw integer RNA-seq counts. Use DESeq2/PyDESeq2 for genuine integer count data; consider limma for normalized microarrays.
- No leakage: data-dependent feature selection, imputation, standardization, model tuning and thresholds must be learned on training folds, not held-out folds. Keep GSE140829 outcomes sealed until the final model/panel is frozen and compatibility is proved. Verify subject independence for group-aware CV.
- Report negative or empty results faithfully. Do not lower significance thresholds post hoc or bypass quality gates to produce a favorable result.

## Existing local files (look for `raw/` or `data/raw/` first)
- `GSE63061_series_matrix.txt.gz`
- `GSE63060_series_matrix.txt.gz`
- `GSE140829_series_matrix.txt.gz`
- `GSE125583_series_matrix.txt.gz`
- `GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- `GPL10558-50081.txt`, `GPL6947-13512.txt`, `GPL15988.txt`
- `gencode.v50.long_noncoding_RNAs.gtf.gz`, `gencode.v50.transcripts.fa.gz`
- `Human.GRCh38.p13.annot.tsv.gz`, `GPL15988_files/`
These are reported filenames, not proof of file contents. Inspect size, first rows, columns, gene IDs, checksums, and GEO metadata. Do not re-download if a usable copy exists; do not permanently delete or overwrite originals. A full/comprehensive GENCODE GTF may be required for reliable coding-vs-lncRNA auditing; request it only when necessary.

## Implementation discipline
1. First inspect the repo, existing tests, data paths, installed versions and the exact files relevant to the stage. Make a compact plan before editing.
2. Work on **one numbered stage at a time**, leaving subsequent stages untouched unless strictly needed to unblock an interface.
3. Prefer Python 3.11, pathlib, pandas, scipy, statsmodels, scikit-learn, pytest and clear type/shape assertions. Avoid unnecessary heavyweight packages. Add R/DESeq2/limma or MATLAB only when the relevant stage calls for them.
4. Keep file formats, ID namespaces, units, source scale, sample order, exclusion rules, random seed (42), and provenance explicit in logs and docstrings.
5. Create reproducible scripts, no hardcoded user's absolute paths, no guessed GEO field names, no network download without first checking local files and asking if necessary.
6. Add small synthetic **unit tests only for parser/algorithm correctness**. Never include synthetic observations in actual research outputs. Keep test fixtures under `tests/`.
7. External network calls, destructive operations, package installs, git push, changing permissions and actual external holdout evaluation require approval unless the user explicitly authorized them.
8. Show meaningful error messages if required files or mapping evidence are missing; stop instead of generating fake results.
9. For each stage deliver: (a) files changed, (b) command(s) to run on Windows PowerShell and optionally Colab, (c) assertions/tests actually run and their pass/fail status, (d) real outputs or a precise blocker, and (e) recommended next stage.
10. Never claim execution or data validation that did not occur. Treat the guide as a plan until actually run on the user's machine.

## Approval gates
- **Gate A:** dataset audit: real per-series sample/feature counts and label dictionary verified.
- **Gate B:** lncRNA annotation: evidence-backed, unique probe-to-gene/biotype mappings; otherwise halt lncRNA-only blood ML.
- **Gate C:** QC/DE: valid expression scale, covariates, sample identity and false-discovery-rate correction.
- **Gate D:** ML: leakage-resistant nested CV and patient-aware splitting where needed.
- **Gate E:** external validation: proven identical or meaningfully comparable validated features, frozen model, no prior use of external outcome labels. Only one final evaluation.

## Status report template
**Stage:** ...  
**Inputs found:** ...  
**Files created/modified:** ...  
**Commands and test results:** ...  
**Scientific / data warnings:** ...  
**Gate:** PASS / BLOCKED / NOT TESTED, with evidence  
**Next action:** ...
