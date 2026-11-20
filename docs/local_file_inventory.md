# Local File Inventory — lncSignal-AD

> Generated: 2026-11-21 02:47:39  
> Stage: 05 — Local Dataset Inventory and Verification

---

## Summary

| Metric | Count |
|--------|-------|
| Total catalogued | 12 |
| ✅ PASS | 11 |
| ❌ FAIL | 0 |
| ⚠️  WARN | 1 |
| 🔴 MISSING | 0 |

---

## Full Inventory Table

| Filename | Location | File Type | Data Purpose | Req? | Exists | Size | SHA256 (full) | Gzip | Rows | Status | Problems |
|----------|----------|-----------|--------------|------|--------|------|---------------|------|------|--------|----------|
| `GSE63061_series_matrix.txt.gz` | `data\raw` | GEO Series Matrix (gzip) | Discovery cohort metadata + expression matrix – AddNeuroMed blood, Illumina HT-1… | Required | YES | 60.8 MB | `c6de1af5bd36bca3…` | OK | 32049 | ✅ PASS | None |
| `GSE63060_series_matrix.txt.gz` | `data\raw` | GEO Series Matrix (gzip) | Replication cohort metadata + expression matrix – AddNeuroMed blood, Illumina HT… | Required | YES | 59.9 MB | `5233e79ab7eece13…` | OK | 38323 | ✅ PASS | None |
| `GSE140829_series_matrix.txt.gz` | `data\raw` | GEO Series Matrix (gzip) | External validation cohort metadata + expression – UCLA/UCSF blood, Illumina HT-… | Required | YES | 215.8 MB | `be4f08a2cb01fe91…` | OK | 47008 | ✅ PASS | None |
| `GSE125583_series_matrix.txt.gz` | `data\raw` | GEO Series Matrix (gzip) | Brain support cohort metadata – fusiform gyrus RNA-seq, 219 AD / 70 ctrl (Genent… | Required | YES | 14.1 KB | `6817baa614825e99…` | OK | 1 | ✅ PASS | None |
| `GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz` | `data\raw` | RNA-seq counts TSV (gzip) | Raw integer counts for GSE125583 (289 samples × 39,376 genes, GRCh38.p13 Entrez … | Required | YES | 11.5 MB | `07f2763c4f7d3909…` | OK | 39376 | ✅ PASS | None |
| `Human.GRCh38.p13.annot.tsv.gz` | `data\raw` | NCBI gene annotation TSV (gzip) | NCBI gene annotation table for GRCh38.p13: GeneID, Symbol, GeneType, Ensembl ID,… | Required | YES | 4.2 MB | `a0b6ff49b0851334…` | OK | 39376 | ✅ PASS | None |
| `GPL10558-50081.txt` | `data\raw` | Illumina probe annotation (plain TSV) | Probe-level annotation for HT-12 v4 (GPL10558) – 48,107 probes, 30 columns incl.… | Required | YES | 68.8 MB | `0f5d64833f263f83…` | N/A (plain) | 48107 | ✅ PASS | None |
| `GPL6947-13512.txt` | `data\raw` | Illumina probe annotation (plain TSV) | Probe-level annotation for HT-12 v3 (GPL6947) – 49,576 probes, 30 columns incl. … | Required | YES | 24.2 MB | `7b5406a5d26a608d…` | N/A (plain) | 49576 | ✅ PASS | None |
| `GPL15988.txt` | `data\raw` | HTML page (INCORRECT FORMAT) | Intended: probe annotation for HT-12 v4 (GPL15988) – ACTUAL: HTML page saved fro… | Required | YES | 2.8 MB | `7f810ba37308029d…` | N/A (plain) | 47322 | ✅ PASS | None |
| `gencode.v50.long_noncoding_RNAs.gtf.gz` | `data\raw` | GENCODE GTF annotation (gzip) | GENCODE v50 lncRNA gene/transcript annotation for GRCh38 – 965,508 GTF records (… | Required | YES | 12.8 MB | `88ef3d23a2d6cf34…` | OK | 965508 | ✅ PASS | None |
| `gencode.v50.transcripts.fa.gz` | `data\optional_raw_data` | GENCODE transcript FASTA (gzip) | GENCODE v50 all transcript sequences – used for probe-to-lncRNA BLAST/alignment … | Optional | YES | 175.1 MB | `5a320f524d73b579…` | OK |  | ✅ PASS | None |
| `GPL15988_files/` | `data\optional_raw_data` | Browser-saved web asset folder | JS, CSS, GIF assets saved alongside GPL15988.txt HTML page – NOT biologically us… | Optional | YES | 9 files | `` | N/A (directory) | 9 | ⚠️ WARN – web assets only, no probe TSV | GPL15988_files/ contains only browser-saved web assets (JS, CSS, GIF) – NOT usable probe annotation. GPL15988.txt is an  |

---

## File-by-File Notes

### `GSE63061_series_matrix.txt.gz`

- **Location:** `data\raw`
- **Type:** GEO Series Matrix (gzip)
- **Purpose:** Discovery cohort metadata + expression matrix – AddNeuroMed blood, Illumina HT-12 v4 (GPL10558), 139 AD / 134 ctrl
- **Required:** Required
- **Exists:** YES
- **Size:** 60.8 MB (63713724 bytes)
- **SHA-256:** `c6de1af5bd36bca3de0921d53051bdefd05666f6e94a6906a26c8cb461dd55f5`
- **Gzip integrity:** OK
- **Data rows:** 32049
- **Status:** PASS
- **Problems:** None
- **First line preview:** `!Series_title	"Alzheimer, Mild Cognitive impairment and control samples from AddneuroMed Cohort (batch 2)"`

### `GSE63060_series_matrix.txt.gz`

- **Location:** `data\raw`
- **Type:** GEO Series Matrix (gzip)
- **Purpose:** Replication cohort metadata + expression matrix – AddNeuroMed blood, Illumina HT-12 v3 (GPL6947), 145 AD / 104 ctrl
- **Required:** Required
- **Exists:** YES
- **Size:** 59.9 MB (62803866 bytes)
- **SHA-256:** `5233e79ab7eece137269153f37fa5ce8b39f000ce72fb3d93240a0d24b63779f`
- **Gzip integrity:** OK
- **Data rows:** 38323
- **Status:** PASS
- **Problems:** None
- **First line preview:** `!Series_title	"Alzheimer, MCI and control samples from AddneuroMed Cohort (batch 1)"`

### `GSE140829_series_matrix.txt.gz`

- **Location:** `data\raw`
- **Type:** GEO Series Matrix (gzip)
- **Purpose:** External validation cohort metadata + expression – UCLA/UCSF blood, Illumina HT-12 v4 (GPL15988), 204 AD / 249 ctrl
- **Required:** Required
- **Exists:** YES
- **Size:** 215.8 MB (226272058 bytes)
- **SHA-256:** `be4f08a2cb01fe912a764df66fd654dec70d9e1ecd36bded00b675b27e7f83f3`
- **Gzip integrity:** OK
- **Data rows:** 47008
- **Status:** PASS
- **Problems:** None
- **First line preview:** `!Series_title	"Systems-level analysis of peripheral blood gene expression in dementia patients reveals an innate immune response shared across multiple disorders [ad_mci]"`

### `GSE125583_series_matrix.txt.gz`

- **Location:** `data\raw`
- **Type:** GEO Series Matrix (gzip)
- **Purpose:** Brain support cohort metadata – fusiform gyrus RNA-seq, 219 AD / 70 ctrl (Genentech/UCSF)
- **Required:** Required
- **Exists:** YES
- **Size:** 14.1 KB (14446 bytes)
- **SHA-256:** `6817baa614825e99e25c29ed501ae2da6ac74e7f89aeb43c902190dadb182667`
- **Gzip integrity:** OK
- **Data rows:** 1
- **Status:** PASS
- **Problems:** None
- **First line preview:** `!Series_title	"RNA Expression in fusiform gyrus of  AD and Control subjects"`

### `GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz`

- **Location:** `data\raw`
- **Type:** RNA-seq counts TSV (gzip)
- **Purpose:** Raw integer counts for GSE125583 (289 samples × 39,376 genes, GRCh38.p13 Entrez IDs)
- **Required:** Required
- **Exists:** YES
- **Size:** 11.5 MB (12084520 bytes)
- **SHA-256:** `07f2763c4f7d3909eb8d74bf1749b90ac5897ca9e189e1749c7122bfc31cd37e`
- **Gzip integrity:** OK
- **Data rows:** 39376
- **Status:** PASS
- **Problems:** None
- **First line preview:** `GeneID	GSM3577568	GSM3577569	GSM3577570	GSM3577571	GSM3577572	GSM3577573	GSM3577574	GSM3577575	GSM3577576	GSM3577577	GSM3577578	GSM3577579	GSM3577580	GSM3577581	GSM3577582	GSM3577583	GSM3577584	GSM357`

### `Human.GRCh38.p13.annot.tsv.gz`

- **Location:** `data\raw`
- **Type:** NCBI gene annotation TSV (gzip)
- **Purpose:** NCBI gene annotation table for GRCh38.p13: GeneID, Symbol, GeneType, Ensembl ID, GO terms (39,376 genes, 18 columns)
- **Required:** Required
- **Exists:** YES
- **Size:** 4.2 MB (4417321 bytes)
- **SHA-256:** `a0b6ff49b0851334c709edbd52ae9c80ae998d66d9a30185de364ffb4a59efea`
- **Gzip integrity:** OK
- **Data rows:** 39376
- **Status:** PASS
- **Problems:** None
- **First line preview:** `GeneID	Symbol	Description	Synonyms	GeneType	EnsemblGeneID	Status	ChrAcc	ChrStart	ChrStop	Orientation	Length	GOFunctionID	GOProcessID	GOComponentID	GOFunction	GOProcess	GOComponent`

### `GPL10558-50081.txt`

- **Location:** `data\raw`
- **Type:** Illumina probe annotation (plain TSV)
- **Purpose:** Probe-level annotation for HT-12 v4 (GPL10558) – 48,107 probes, 30 columns incl. Entrez, Symbol, SEQUENCE
- **Required:** Required
- **Exists:** YES
- **Size:** 68.8 MB (72095087 bytes)
- **SHA-256:** `0f5d64833f263f8304e66b5279664a69b203977de2277d426584267dd8c40c9f`
- **Gzip integrity:** N/A (plain)
- **Data rows:** 48107
- **Status:** PASS
- **Problems:** None
- **First line preview:** `#ID = Unique identifier for the probe (across all products and species)`

### `GPL6947-13512.txt`

- **Location:** `data\raw`
- **Type:** Illumina probe annotation (plain TSV)
- **Purpose:** Probe-level annotation for HT-12 v3 (GPL6947) – 49,576 probes, 30 columns incl. Entrez, Symbol, SEQUENCE
- **Required:** Required
- **Exists:** YES
- **Size:** 24.2 MB (25399113 bytes)
- **SHA-256:** `7b5406a5d26a608d1d2eb75ea7a7b0f573b0c49075ddb5884f1122cefd077cc0`
- **Gzip integrity:** N/A (plain)
- **Data rows:** 49576
- **Status:** PASS
- **Problems:** None
- **First line preview:** `#ID = Unique identifier for the probe (across all products and species)`

### `GPL15988.txt`

- **Location:** `data\raw`
- **Type:** HTML page (INCORRECT FORMAT)
- **Purpose:** Intended: probe annotation for HT-12 v4 (GPL15988) – ACTUAL: HTML page saved from NCBI GEO website, NOT usable as TSV
- **Required:** Required
- **Exists:** YES
- **Size:** 2.8 MB (2921503 bytes)
- **SHA-256:** `7f810ba37308029d87162dc1e81be590315688ad5cf2aaeb09fcc9a48850fde1`
- **Gzip integrity:** N/A (plain)
- **Data rows:** 47322
- **Status:** PASS
- **Problems:** None
- **First line preview:** `#ID =`

### `gencode.v50.long_noncoding_RNAs.gtf.gz`

- **Location:** `data\raw`
- **Type:** GENCODE GTF annotation (gzip)
- **Purpose:** GENCODE v50 lncRNA gene/transcript annotation for GRCh38 – 965,508 GTF records (Ensembl 116)
- **Required:** Required
- **Exists:** YES
- **Size:** 12.8 MB (13377098 bytes)
- **SHA-256:** `88ef3d23a2d6cf3422e188a4d5506ba5386b2c05283a99c916794b5cf9f2b6bf`
- **Gzip integrity:** OK
- **Data rows:** 965508
- **Status:** PASS
- **Problems:** None
- **First line preview:** `##description: evidence-based annotation of the human genome (GRCh38), version 50 (Ensembl 116) - long non-coding RNAs`

### `gencode.v50.transcripts.fa.gz`

- **Location:** `data\optional_raw_data`
- **Type:** GENCODE transcript FASTA (gzip)
- **Purpose:** GENCODE v50 all transcript sequences – used for probe-to-lncRNA BLAST/alignment (optional in early stages)
- **Required:** Optional
- **Exists:** YES
- **Size:** 175.1 MB (183554921 bytes)
- **SHA-256:** `5a320f524d73b5793518eb19b118829033713443d0f42af20a67bb31cc06cf56`
- **Gzip integrity:** OK
- **Data rows:** 
- **Status:** PASS
- **Problems:** None
- **First line preview:** `>ENST00000832824.1|ENSG00000290825.2|-|-|DDX11L16-260|DDX11L16|1379|lncRNA|`

### `GPL15988_files/`

- **Location:** `data\optional_raw_data`
- **Type:** Browser-saved web asset folder
- **Purpose:** JS, CSS, GIF assets saved alongside GPL15988.txt HTML page – NOT biologically useful annotation data
- **Required:** Optional
- **Exists:** YES
- **Size:** 9 files ( bytes)
- **SHA-256:** ``
- **Gzip integrity:** N/A (directory)
- **Data rows:** 9
- **Status:** WARN – web assets only, no probe TSV
- **Problems:** GPL15988_files/ contains only browser-saved web assets (JS, CSS, GIF) – NOT usable probe annotation. GPL15988.txt is an HTML page, not a TSV. ACTION REQUIRED: download GPL15988 probe annotation as TSV.
- **First line preview:** `bg_main3.gif; CoreModule.js.download; gtm.js.download; InstrumentPageStarterJS.js.download; js`

---

## Action Items

### 🔴 Critical (blocks downstream analysis)

1. **GPL15988.txt is an HTML page, not a TSV.**  
   The file `data/raw/GPL15988.txt` was saved as a browser HTML page.  
   The `GPL15988_files/` folder contains only browser JS/CSS assets — no probe data.  
   **Fix:** Download the proper GPL15988 probe annotation TSV from:
   ```
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15988
   ```
   Click **Download full table** → plain text (not the SOFT .gz).  
   Save as `data/raw/GPL15988.txt` (overwriting the HTML).  
   This file is required for the GSE140829 external validation cohort.

### ✅ No action needed

- All 10 other files pass integrity checks.
- `gencode.v50.transcripts.fa.gz` is present in `data/optional_raw_data/` (183 MB, gzip OK).
- `GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz` is present (289 samples × 39,376 genes).
- Python environment has all required packages.

---

## Environment Check

| Package | Version | Status |
|---------|---------|--------|
| Python  | 3.12.10 | ✅ |
| pandas  | 3.1.0   | ✅ |
| numpy   | 2.5.4   | ✅ |
| scipy   | 1.18.1  | ✅ |
| statsmodels | 0.15.0 | ✅ |
| biopython | 1.88 | ✅ |
| scikit-learn | 1.9.1 | ✅ |
| gseapy | 1.3.2 | ✅ |
| pydeseq2 | 0.6.0 | ✅ |
| matplotlib | 3.11.2 | ✅ |
| seaborn | 0.13.2 | ✅ |
| pyarrow | 26.0.0 | ✅ |
