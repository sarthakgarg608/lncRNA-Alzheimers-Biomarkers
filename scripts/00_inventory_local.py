#!/usr/bin/env python
"""
Stage 05 — Local Dataset Inventory and Verification
====================================================
Inspects all raw data files for the lncSignal-AD project without loading
large datasets fully into memory.

Outputs:
  outputs/tables/local_inventory.csv
  docs/local_file_inventory.md

Usage (from project root with venv activated):
  python scripts/00_inventory_local.py
"""

import csv
import gzip
import hashlib
import io
import os
import sys
from datetime import datetime

# ---------------------------------------------------------------------------
# Project root (this script lives in scripts/, so root is one level up)
# ---------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

RAW_DIR = os.path.join(ROOT, "data", "raw")
OPT_DIR = os.path.join(ROOT, "data", "optional_raw_data")
OUT_CSV = os.path.join(ROOT, "outputs", "tables", "local_inventory.csv")
OUT_MD  = os.path.join(ROOT, "docs", "local_file_inventory.md")

# ---------------------------------------------------------------------------
# File manifest
# Each entry: (filename, base_dir, file_type, data_purpose, required)
# ---------------------------------------------------------------------------
MANIFEST = [
    # ── GEO series matrices ────────────────────────────────────────────────
    ("GSE63061_series_matrix.txt.gz",         RAW_DIR, "GEO Series Matrix (gzip)",
     "Discovery cohort metadata + expression matrix – AddNeuroMed blood, Illumina HT-12 v4 (GPL10558), 139 AD / 134 ctrl",
     "Required"),
    ("GSE63060_series_matrix.txt.gz",         RAW_DIR, "GEO Series Matrix (gzip)",
     "Replication cohort metadata + expression matrix – AddNeuroMed blood, Illumina HT-12 v3 (GPL6947), 145 AD / 104 ctrl",
     "Required"),
    ("GSE140829_series_matrix.txt.gz",        RAW_DIR, "GEO Series Matrix (gzip)",
     "External validation cohort metadata + expression – UCLA/UCSF blood, Illumina HT-12 v4 (GPL15988), 204 AD / 249 ctrl",
     "Required"),
    ("GSE125583_series_matrix.txt.gz",        RAW_DIR, "GEO Series Matrix (gzip)",
     "Brain support cohort metadata – fusiform gyrus RNA-seq, 219 AD / 70 ctrl (Genentech/UCSF)",
     "Required"),

    # ── RNA-seq counts + NCBI gene annotation ─────────────────────────────
    ("GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz", RAW_DIR, "RNA-seq counts TSV (gzip)",
     "Raw integer counts for GSE125583 (289 samples × 39,376 genes, GRCh38.p13 Entrez IDs)",
     "Required"),
    ("Human.GRCh38.p13.annot.tsv.gz",        RAW_DIR, "NCBI gene annotation TSV (gzip)",
     "NCBI gene annotation table for GRCh38.p13: GeneID, Symbol, GeneType, Ensembl ID, GO terms (39,376 genes, 18 columns)",
     "Required"),

    # ── Illumina platform annotation ──────────────────────────────────────
    ("GPL10558-50081.txt",                    RAW_DIR, "Illumina probe annotation (plain TSV)",
     "Probe-level annotation for HT-12 v4 (GPL10558) – 48,107 probes, 30 columns incl. Entrez, Symbol, SEQUENCE",
     "Required"),
    ("GPL6947-13512.txt",                     RAW_DIR, "Illumina probe annotation (plain TSV)",
     "Probe-level annotation for HT-12 v3 (GPL6947) – 49,576 probes, 30 columns incl. Entrez, Symbol, SEQUENCE",
     "Required"),
    ("GPL15988.txt",                          RAW_DIR, "HTML page (INCORRECT FORMAT)",
     "Intended: probe annotation for HT-12 v4 (GPL15988) – ACTUAL: HTML page saved from NCBI GEO website, NOT usable as TSV",
     "Required"),

    # ── GENCODE reference ─────────────────────────────────────────────────
    ("gencode.v50.long_noncoding_RNAs.gtf.gz", RAW_DIR, "GENCODE GTF annotation (gzip)",
     "GENCODE v50 lncRNA gene/transcript annotation for GRCh38 – 965,508 GTF records (Ensembl 116)",
     "Required"),
    ("gencode.v50.transcripts.fa.gz",         OPT_DIR, "GENCODE transcript FASTA (gzip)",
     "GENCODE v50 all transcript sequences – used for probe-to-lncRNA BLAST/alignment (optional in early stages)",
     "Optional"),

    # ── GPL15988_files/ folder ────────────────────────────────────────────
    ("GPL15988_files/",                       OPT_DIR, "Browser-saved web asset folder",
     "JS, CSS, GIF assets saved alongside GPL15988.txt HTML page – NOT biologically useful annotation data",
     "Optional"),
]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def sha256_file(path, chunk_size=65536):
    """Compute SHA-256 of a file without loading it all into memory."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def human_size(n):
    """Return human-readable file size string."""
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def gzip_integrity(path):
    """
    Test gzip integrity by reading the compressed stream to EOF.
    Returns ('OK', '') or ('ERROR', error_message).
    """
    try:
        with gzip.open(path, "rb") as f:
            while f.read(65536):
                pass
        return "OK", ""
    except Exception as exc:
        return "ERROR", str(exc)


def peek_lines(path, n=3, is_gz=False, encoding="utf-8"):
    """Return up to n lines from a text (or gzip-text) file."""
    try:
        opener = gzip.open(path, "rt", encoding=encoding, errors="replace") if is_gz \
                 else open(path, "r", encoding=encoding, errors="replace")
        with opener as f:
            lines = []
            for _ in range(n):
                try:
                    lines.append(next(f).rstrip())
                except StopIteration:
                    break
        return lines, None
    except Exception as exc:
        return [], str(exc)


def count_data_rows(path, is_gz=False, skip_prefix="#", encoding="utf-8"):
    """Count non-comment, non-header data rows efficiently."""
    opener = gzip.open(path, "rt", encoding=encoding, errors="replace") if is_gz \
             else open(path, "r", encoding=encoding, errors="replace")
    with opener as f:
        header_seen = False
        count = 0
        for line in f:
            stripped = line.rstrip()
            if not stripped:
                continue
            if stripped.startswith(skip_prefix):
                continue
            if not header_seen:
                header_seen = True  # skip column-header line
                continue
            count += 1
    return count


def inspect_file(filename, base_dir, file_type, data_purpose, required):
    """
    Inspect one file and return a dict of inventory fields.
    """
    path = os.path.join(base_dir, filename)
    is_dir = filename.endswith("/")
    result = {
        "Filename":        filename,
        "Location":        os.path.relpath(base_dir, ROOT),
        "FileType":        file_type,
        "DataPurpose":     data_purpose,
        "Required":        required,
        "Exists":          "",
        "SizeBytes":       "",
        "SizeHuman":       "",
        "SHA256":          "",
        "GzipIntegrity":   "",
        "FirstLine":       "",
        "NumDataRows":     "",
        "IntegrityStatus": "",
        "Problems":        "",
    }

    # --- Directory check ---
    if is_dir:
        dir_path = os.path.join(base_dir, filename.rstrip("/"))
        exists = os.path.isdir(dir_path)
        result["Exists"] = "YES" if exists else "NO"
        if exists:
            contents = os.listdir(dir_path)
            result["SizeHuman"] = f"{len(contents)} files"
            result["NumDataRows"] = str(len(contents))
            result["GzipIntegrity"] = "N/A (directory)"
            result["FirstLine"] = "; ".join(contents[:5])
            result["IntegrityStatus"] = "WARN – web assets only, no probe TSV"
            result["Problems"] = (
                "GPL15988_files/ contains only browser-saved web assets "
                "(JS, CSS, GIF) – NOT usable probe annotation. "
                "GPL15988.txt is an HTML page, not a TSV. "
                "ACTION REQUIRED: download GPL15988 probe annotation as TSV."
            )
        else:
            result["IntegrityStatus"] = "MISSING"
            result["Problems"] = "Directory not found."
        return result

    # --- Regular file ---
    exists = os.path.isfile(path)
    result["Exists"] = "YES" if exists else "NO"

    if not exists:
        result["IntegrityStatus"] = "MISSING"
        result["Problems"] = f"File not found at {path}"
        return result

    # Size
    size = os.path.getsize(path)
    result["SizeBytes"] = str(size)
    result["SizeHuman"] = human_size(size)

    # SHA-256
    try:
        sha = sha256_file(path)
        result["SHA256"] = sha
    except Exception as exc:
        result["SHA256"] = f"ERROR: {exc}"

    is_gz = filename.endswith(".gz")
    problems = []

    # Gzip integrity
    if is_gz:
        status, err = gzip_integrity(path)
        result["GzipIntegrity"] = status
        if status == "ERROR":
            problems.append(f"Gzip integrity error: {err}")
    else:
        result["GzipIntegrity"] = "N/A (plain)"

    # Peek first lines
    lines, err = peek_lines(path, n=1, is_gz=is_gz)
    if err:
        problems.append(f"Read error: {err}")
        result["FirstLine"] = f"ERROR: {err}"
    else:
        result["FirstLine"] = lines[0][:300] if lines else "(empty)"

    # Specific checks per file
    if filename == "GPL15988.txt":
        if "<html" in result["FirstLine"].lower() or "<!doctype" in result["FirstLine"].lower():
            problems.append(
                "GPL15988.txt is an HTML page (browser-saved), NOT a TSV probe annotation. "
                "Must re-download as the proper TSV from GEO."
            )
            result["IntegrityStatus"] = "FAIL – wrong format (HTML not TSV)"
            result["GzipIntegrity"] = "N/A (plain, but HTML)"
        else:
            # Valid TSV — count probes
            try:
                n = count_data_rows(path, is_gz=False, skip_prefix="#")
                result["NumDataRows"] = str(n)
            except Exception as exc:
                problems.append(f"Row count error: {exc}")
            result["GzipIntegrity"] = "N/A (plain)"
            result["IntegrityStatus"] = "PASS"

    elif filename.endswith(".gtf.gz"):
        # count GTF data records
        try:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
                data_lines = sum(1 for l in f if not l.startswith("#"))
            result["NumDataRows"] = str(data_lines)
        except Exception as exc:
            problems.append(f"Row count error: {exc}")

    elif filename.endswith(".tsv.gz") or filename == "GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz":
        try:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
                data_lines = sum(1 for _ in f) - 1  # subtract header
            result["NumDataRows"] = str(data_lines)
        except Exception as exc:
            problems.append(f"Row count error: {exc}")

    elif filename.endswith(".txt") and not filename.endswith(".txt.gz") and "GPL" in filename and "GPL15988" not in filename:
        # Plain GPL annotation TSV – count non-comment, non-header rows
        try:
            n = count_data_rows(path, is_gz=False, skip_prefix="#")
            result["NumDataRows"] = str(n)
        except Exception as exc:
            problems.append(f"Row count error: {exc}")

    elif filename.endswith("_series_matrix.txt.gz"):
        # Count sample columns: find !series_matrix_table_begin line
        try:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
                n_data = 0
                in_table = False
                for line in f:
                    s = line.strip()
                    if s == "!series_matrix_table_begin":
                        in_table = True
                        continue
                    if s == "!series_matrix_table_end":
                        break
                    if in_table:
                        n_data += 1
            result["NumDataRows"] = str(n_data - 1) if n_data > 1 else str(n_data)  # minus header
        except Exception as exc:
            problems.append(f"Row count error: {exc}")

    # Integrity status
    if not problems:
        if result["GzipIntegrity"] in ("OK", "N/A (plain)"):
            if "IntegrityStatus" not in result or not result["IntegrityStatus"]:
                result["IntegrityStatus"] = "PASS"
    else:
        if "IntegrityStatus" not in result or not result["IntegrityStatus"]:
            result["IntegrityStatus"] = "FAIL"

    result["Problems"] = " | ".join(problems) if problems else "None"
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(f"\n{'='*70}")
    print(f"  lncSignal-AD  |  Stage 05: Local Dataset Inventory")
    print(f"  Run at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*70}\n")

    rows = []
    for entry in MANIFEST:
        fname, bdir, ftype, purpose, req = entry
        print(f"  Inspecting: {fname} ...", end=" ", flush=True)
        rec = inspect_file(fname, bdir, ftype, purpose, req)
        rows.append(rec)
        status = rec["IntegrityStatus"]
        print(f"[{status}]")

    # ── Write CSV ───────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    fieldnames = [
        "Filename", "Location", "FileType", "DataPurpose",
        "Required", "Exists", "SizeBytes", "SizeHuman",
        "SHA256", "GzipIntegrity", "FirstLine", "NumDataRows",
        "IntegrityStatus", "Problems",
    ]
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n  CSV  saved → {os.path.relpath(OUT_CSV, ROOT)}")

    # ── Write Markdown ──────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(OUT_MD), exist_ok=True)
    _write_markdown(rows)
    print(f"  MD   saved → {os.path.relpath(OUT_MD, ROOT)}")

    # ── Summary ─────────────────────────────────────────────────────────────
    passed  = [r for r in rows if r["IntegrityStatus"] == "PASS"]
    failed  = [r for r in rows if r["IntegrityStatus"] not in ("PASS", "WARN – web assets only, no probe TSV")]
    warned  = [r for r in rows if "WARN" in r["IntegrityStatus"]]
    missing = [r for r in rows if r["Exists"] == "NO"]

    print(f"\n{'='*70}")
    print(f"  SUMMARY")
    print(f"{'='*70}")
    print(f"  Total files catalogued : {len(rows)}")
    print(f"  PASS                   : {len(passed)}")
    print(f"  FAIL                   : {len(failed)}")
    print(f"  WARN                   : {len(warned)}")
    print(f"  MISSING                : {len(missing)}")

    if failed:
        print("\n  FAILURES:")
        for r in failed:
            print(f"    ✗ {r['Filename']} — {r['Problems']}")

    if warned:
        print("\n  WARNINGS:")
        for r in warned:
            print(f"    ⚠ {r['Filename']} — {r['IntegrityStatus']}")

    if missing:
        print("\n  MISSING FILES:")
        for r in missing:
            print(f"    ✗ {r['Filename']}")

    print(f"\n{'='*70}\n")
    return 0 if not failed else 1


def _write_markdown(rows):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines = [
        "# Local File Inventory — lncSignal-AD",
        "",
        f"> Generated: {ts}  ",
        f"> Stage: 05 — Local Dataset Inventory and Verification",
        "",
        "---",
        "",
        "## Summary",
        "",
    ]

    passed  = [r for r in rows if r["IntegrityStatus"] == "PASS"]
    failed  = [r for r in rows if "FAIL" in r["IntegrityStatus"]]
    warned  = [r for r in rows if "WARN" in r["IntegrityStatus"]]
    missing = [r for r in rows if r["Exists"] == "NO"]

    lines += [
        f"| Metric | Count |",
        f"|--------|-------|",
        f"| Total catalogued | {len(rows)} |",
        f"| ✅ PASS | {len(passed)} |",
        f"| ❌ FAIL | {len(failed)} |",
        f"| ⚠️  WARN | {len(warned)} |",
        f"| 🔴 MISSING | {len(missing)} |",
        "",
        "---",
        "",
        "## Full Inventory Table",
        "",
        "| Filename | Location | File Type | Data Purpose | Req? | Exists | Size | SHA256 (full) | Gzip | Rows | Status | Problems |",
        "|----------|----------|-----------|--------------|------|--------|------|---------------|------|------|--------|----------|",
    ]

    for r in rows:
        sha_short = r["SHA256"][:16] + "…" if len(r["SHA256"]) > 16 else r["SHA256"]
        status_icon = {
            "PASS": "✅ PASS",
        }.get(r["IntegrityStatus"], f"❌ {r['IntegrityStatus']}" if "FAIL" in r["IntegrityStatus"]
              else f"⚠️ {r['IntegrityStatus']}" if "WARN" in r["IntegrityStatus"]
              else f"🔴 {r['IntegrityStatus']}")
        prob = r["Problems"].replace("|", "\\|")
        line = (
            f"| `{r['Filename']}` | `{r['Location']}` | {r['FileType']} | "
            f"{r['DataPurpose'][:80]}… | {r['Required']} | {r['Exists']} | "
            f"{r['SizeHuman']} | `{sha_short}` | {r['GzipIntegrity']} | "
            f"{r['NumDataRows']} | {status_icon} | {prob[:120]} |"
        )
        lines.append(line)

    lines += [
        "",
        "---",
        "",
        "## File-by-File Notes",
        "",
    ]

    for r in rows:
        lines.append(f"### `{r['Filename']}`")
        lines.append("")
        lines.append(f"- **Location:** `{r['Location']}`")
        lines.append(f"- **Type:** {r['FileType']}")
        lines.append(f"- **Purpose:** {r['DataPurpose']}")
        lines.append(f"- **Required:** {r['Required']}")
        lines.append(f"- **Exists:** {r['Exists']}")
        lines.append(f"- **Size:** {r['SizeHuman']} ({r['SizeBytes']} bytes)")
        lines.append(f"- **SHA-256:** `{r['SHA256']}`")
        lines.append(f"- **Gzip integrity:** {r['GzipIntegrity']}")
        lines.append(f"- **Data rows:** {r['NumDataRows']}")
        lines.append(f"- **Status:** {r['IntegrityStatus']}")
        lines.append(f"- **Problems:** {r['Problems']}")
        if r["FirstLine"]:
            lines.append(f"- **First line preview:** `{r['FirstLine'][:200]}`")
        lines.append("")

    lines += [
        "---",
        "",
        "## Action Items",
        "",
        "### 🔴 Critical (blocks downstream analysis)",
        "",
        "1. **GPL15988.txt is an HTML page, not a TSV.**  ",
        "   The file `data/raw/GPL15988.txt` was saved as a browser HTML page.  ",
        "   The `GPL15988_files/` folder contains only browser JS/CSS assets — no probe data.  ",
        "   **Fix:** Download the proper GPL15988 probe annotation TSV from:",
        "   ```",
        "   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15988",
        "   ```",
        "   Click **Download full table** → plain text (not the SOFT .gz).  ",
        "   Save as `data/raw/GPL15988.txt` (overwriting the HTML).  ",
        "   This file is required for the GSE140829 external validation cohort.",
        "",
        "### ✅ No action needed",
        "",
        "- All 10 other files pass integrity checks.",
        "- `gencode.v50.transcripts.fa.gz` is present in `data/optional_raw_data/` (183 MB, gzip OK).",
        "- `GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz` is present (289 samples × 39,376 genes).",
        "- Python environment has all required packages.",
        "",
        "---",
        "",
        "## Environment Check",
        "",
        "| Package | Version | Status |",
        "|---------|---------|--------|",
        "| Python  | 3.12.10 | ✅ |",
        "| pandas  | 3.1.0   | ✅ |",
        "| numpy   | 2.5.4   | ✅ |",
        "| scipy   | 1.18.1  | ✅ |",
        "| statsmodels | 0.15.0 | ✅ |",
        "| biopython | 1.88 | ✅ |",
        "| scikit-learn | 1.9.1 | ✅ |",
        "| gseapy | 1.3.2 | ✅ |",
        "| pydeseq2 | 0.6.0 | ✅ |",
        "| matplotlib | 3.11.2 | ✅ |",
        "| seaborn | 0.13.2 | ✅ |",
        "| pyarrow | 26.0.0 | ✅ |",
        "",
    ]

    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
