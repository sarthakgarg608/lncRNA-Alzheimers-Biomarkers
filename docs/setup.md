# Environment Setup and Verification Guide — lncSignal-AD

> **Project:** Establishing the Role of lncRNAs as Biomarkers in Alzheimer's Disease  
> **Stage:** 06 — Environment Verification  
> **Last Verified:** 2026-11-21  

---

## 1. Overview and Python Version Policy

This project utilizes a standardized scientific computing environment for processing microarray intensity data (Illumina HumanHT-12 v3/v4), fusiform gyrus RNA-seq counts, GENCODE transcript annotations, differential expression modeling, and machine-learning evaluation.

### Python Version Consistency Note

- **Documentation Baseline (`AGENTS.md`, `PROJECT_GUIDE.md`):** Recommends Python 3.11 as a conservative baseline.
- **Active Environment:** **Python 3.12.10** (`AMD64`, Windows 11 Build 10.0.26300).
- **Compatibility Status:** All 18 packages in `requirements.txt` (including `pydeseq2`, `gseapy`, `scikit-learn`, `biopython`, and `pyarrow`) have binary wheels and full compatibility under Python 3.12. `python -m pip check` confirms **no broken requirements**.
- **Recommendation:** Retain the active Python 3.12.10 virtual environment (`.venv`) as it satisfies all required version constraints without modification.

---

## 2. Windows PowerShell Setup (Local Workstation)

### Step 1: Initialize Virtual Environment

From the project root (`C:\Projects\lncSignal-AD`):

```powershell
# Using Python 3.12 (or 3.11 if preferred)
py -3.12 -m venv .venv

# Activate environment
.\.venv\Scripts\Activate.ps1
```

*Note on PowerShell Execution Policy:* If script execution is restricted, activate using:
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
.\.venv\Scripts\Activate.ps1
```
Or execute Python directly via `.\.venv\Scripts\python.exe`.

### Step 2: Install Required Dependencies

```powershell
# Upgrade core packaging tools
.\.venv\Scripts\python.exe -m pip install --upgrade pip setuptools wheel

# Install all dependencies from requirements.txt
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

### Step 3: Run Environment Smoke Test

Run the automated verification script:

```powershell
.\.venv\Scripts\python.exe scripts/check_environment.py
```

Run test suite:

```powershell
.\.venv\Scripts\pytest.exe tests
```

---

## 3. Google Colab Setup (Cloud GPU / High-RAM)

To replicate the analysis or run compute-intensive steps in Google Colab:

```python
# 1. Clone repository and navigate to root
!git clone https://github.com/<user>/lncSignal-AD.git
%cd lncSignal-AD

# 2. Install requirements in Colab runtime
!pip install --upgrade pip
!pip install -r requirements.txt

# 3. Execute environment verification smoke test
!python scripts/check_environment.py

# 4. Run automated test suite
!pytest tests
```

---

## 4. Dependency Tiers & Optional Packages

To keep the computing environment clean and avoid dependency conflicts, packages are grouped into distinct tiers:

### Tier 1: Essential Python Packages (Installed & Verified)

These packages are mandatory for pipeline stages 05 through 12:

| Package | Minimum Required | Installed Version | Purpose |
|---|---|---|---|
| `pandas` | `>=2.0` | 3.1.0 | Microarray & metadata tabular data manipulation |
| `numpy` | `>=1.24` | 2.5.4 | Numerical arrays, matrix operations |
| `scipy` | `>=1.10` | 1.18.1 | Statistical distributions and hypothesis testing |
| `statsmodels` | `>=0.14` | 0.15.0 | Linear models, Benjamini-Hochberg FDR correction |
| `GEOparse` | `>=2.0` | 2.0.4 | Parsing GEO series matrices and platform records |
| `pydeseq2` | `>=0.4` | 0.6.0 | Negative binomial GLM for GSE125583 RNA-seq counts |
| `biopython` | `>=1.81` | 1.88 | FASTA parsing, sequence matching, BLAST coordination |
| `gseapy` | `>=1.1` | 1.3.2 | Gene ontology (GO) and pathway enrichment analysis |
| `scikit-learn` | `>=1.3` | 1.9.1 | LASSO, Random Forest, SVM, nested CV, ROC/AUC |
| `joblib` | `>=1.3` | 1.6.0 | Model serialization and parallelized cross-validation |
| `matplotlib` | `>=3.7` | 3.11.2 | Visualization (boxplots, PCA, volcano plots, ROC curves) |
| `seaborn` | `>=0.12` | 0.13.2 | Statistical heatmaps, cluster maps |
| `networkx` | `>=3.0` | 3.7 | Biological interaction networks and graphs |
| `pyarrow` | `>=14.0` | 26.0.0 | High-performance Parquet format storage for expression |
| `jupyter` | `>=1.0` | 1.1.1 | Interactive notebooks (`notebooks.ipynb`) |
| `pytest` | `>=8.0` | 9.2.0 | Unit and regression testing (`tests/`) |

### Tier 2: Optional Python Packages (Installed & Verified)

| Package | Minimum Required | Installed Version | Purpose |
|---|---|---|---|
| `streamlit` | `>=1.30` | 1.67.0 | Interactive research dashboard (`app.py`) |
| `shap` | `>=0.45` | 0.54.1 | SHAP explainability analysis for ML biomarkers |

### Tier 3: External Non-Python Environments (Optional / As Needed)

Do **not** install these until specifically called for in downstream stages:

1. **R / Bioconductor (Optional for Limma / DESeq2 comparison):**
   ```r
   install.packages("BiocManager")
   BiocManager::install(c("GEOquery", "limma", "DESeq2", "edgeR"))
   ```
   *Role:* Used only if benchmarking Python linear models against standard Bioconductor `limma` for normalized microarrays.

2. **MATLAB Bioinformatics Toolbox (Course-specific optional):**
   *Role:* Available if exploring GEO datasets using MATLAB routines.

---

## 5. Verification Checklist

Before advancing beyond Stage 06, ensure:

- [x] Stage 05 (Local Dataset Inventory) is confirmed **PASS** with all 11 required files verified.
- [x] Python version is 3.11 or 3.12 (active: 3.12.10).
- [x] `python -m pip check` reports `No broken requirements found.`
- [x] `python scripts/check_environment.py` returns exit code 0 and writes `outputs/tables/environment_report.csv`.
- [x] `python -m compileall scripts` compiles without syntax errors.
- [x] `pytest tests` passes all unit tests (17 passed).
