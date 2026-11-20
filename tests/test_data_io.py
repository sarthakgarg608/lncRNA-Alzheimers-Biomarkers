"""
tests/test_data_io.py
=====================
Tests for scripts/utils/data_io.py

Where possible, synthetic in-memory fixtures are used so that no large raw
data files need to be fully loaded.  A small number of tests perform real
lightweight reads (first N lines) from the actual raw files to verify that
parsers handle the true format.

Synthetic data is for SOFTWARE TESTING ONLY — never used as research evidence.
"""

from __future__ import annotations

import gzip
import io
import textwrap
from pathlib import Path

import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# Synthetic fixtures
# ---------------------------------------------------------------------------

SYNTHETIC_SERIES_MATRIX = textwrap.dedent("""\
    !Series_title\t"Synthetic test dataset"
    !Series_geo_accession\t"GSE00000"
    !Series_platform_id\t"GPL10558"
    !series_matrix_table_begin
    "ID_REF"\t"GSM001"\t"GSM002"\t"GSM003"
    "ILMN_A"\t10.0\t11.0\t12.0
    "ILMN_B"\t20.0\t21.0\t22.0
    "ILMN_C"\t9.5\t8.5\t7.5
    !series_matrix_table_end
""")

SYNTHETIC_TSV_GZ_CONTENT = "GeneID\tGSM001\tGSM002\n100\t10\t20\n200\t30\t40\n"

SYNTHETIC_PLAIN_TSV = textwrap.dedent("""\
    # ID = probe identifier
    # SYMBOL = gene symbol
    ID\tSYMBOL\tDEFINITION
    ILMN_A\tNEAT1\tNuclear enriched abundant transcript 1
    ILMN_B\tXIST\tX inactive specific transcript
""")

SYNTHETIC_GTF_LINES = textwrap.dedent("""\
    ##description: synthetic test GTF
    ##provider: GENCODE
    ##format: gtf
    chr1\tHAVANA\tgene\t11121\t24894\t.\t+\t.\tgene_id "ENSG001.1"; gene_type "lncRNA";
""")


def _make_gz_file(content: str, tmp_path: Path, filename: str) -> Path:
    p = tmp_path / filename
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write(content)
    return p


def _make_plain_file(content: str, tmp_path: Path, filename: str) -> Path:
    p = tmp_path / filename
    p.write_text(content, encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# Tests: load_series_matrix_header
# ---------------------------------------------------------------------------

class TestLoadSeriesMatrixHeader:
    def test_synthetic_header_parsing(self, tmp_path):
        from scripts.utils.data_io import load_series_matrix_header
        gz_path = _make_gz_file(SYNTHETIC_SERIES_MATRIX, tmp_path, "test.txt.gz")
        header = load_series_matrix_header(gz_path)
        assert "Series_geo_accession" in header
        assert header["Series_geo_accession"][0][0] == "GSE00000"

    def test_stops_before_table(self, tmp_path):
        """Header loader must NOT read expression data rows."""
        from scripts.utils.data_io import load_series_matrix_header
        gz_path = _make_gz_file(SYNTHETIC_SERIES_MATRIX, tmp_path, "test.txt.gz")
        header = load_series_matrix_header(gz_path)
        # Expression probes should NOT appear in header dict
        all_values = str(header)
        assert "ILMN_A" not in all_values

    def test_real_gse63061_has_geo_accession(self):
        """Light smoke test on the real file — reads header only."""
        from scripts.utils.paths import resolve_path
        from scripts.utils.data_io import load_series_matrix_header
        path = resolve_path("data/raw/GSE63061_series_matrix.txt.gz")
        if not path.exists():
            pytest.skip("GSE63061_series_matrix.txt.gz not found")
        header = load_series_matrix_header(path)
        assert "Series_geo_accession" in header
        accession = header["Series_geo_accession"][0][0]
        assert accession == "GSE63061", f"Unexpected accession: {accession}"

    def test_real_gse140829_has_geo_accession(self):
        """Light smoke test on GSE140829 — reads header only."""
        from scripts.utils.paths import resolve_path
        from scripts.utils.data_io import load_series_matrix_header
        path = resolve_path("data/raw/GSE140829_series_matrix.txt.gz")
        if not path.exists():
            pytest.skip("GSE140829_series_matrix.txt.gz not found")
        # Stop in the series-level preamble, before any sample characteristics.
        import gzip
        with gzip.open(path, "rt") as stream:
            for line in stream:
                assert not line.startswith("!Sample_"), "Missing series accession"
                if line.startswith("!Series_geo_accession\t"):
                    assert line.rstrip().split("\t")[1].strip('"') == "GSE140829"
                    break
            else:
                pytest.fail("Missing GSE140829 series accession")

    def test_real_gse125583_platform(self):
        """GSE125583 platform must be GPL16791 (HiSeq RNA-seq, not microarray)."""
        from scripts.utils.paths import resolve_path
        from scripts.utils.data_io import load_series_matrix_header
        path = resolve_path("data/raw/GSE125583_series_matrix.txt.gz")
        if not path.exists():
            pytest.skip("GSE125583_series_matrix.txt.gz not found")
        header = load_series_matrix_header(path)
        plat = header.get("Series_platform_id", [[""]])[0][0]
        assert plat == "GPL16791", (
            f"GSE125583 platform should be GPL16791 (RNA-seq), got '{plat}'"
        )


# ---------------------------------------------------------------------------
# Tests: load_series_matrix_expression  (synthetic only — large files
# would be too slow in a test suite)
# ---------------------------------------------------------------------------

class TestLoadSeriesMatrixExpression:
    def test_synthetic_shape(self, tmp_path):
        from scripts.utils.data_io import load_series_matrix_expression
        gz_path = _make_gz_file(SYNTHETIC_SERIES_MATRIX, tmp_path, "test.txt.gz")
        expr = load_series_matrix_expression(gz_path)
        assert expr.shape == (3, 3), f"Expected (3,3), got {expr.shape}"

    def test_index_is_feature_id(self, tmp_path):
        from scripts.utils.data_io import load_series_matrix_expression
        gz_path = _make_gz_file(SYNTHETIC_SERIES_MATRIX, tmp_path, "test.txt.gz")
        expr = load_series_matrix_expression(gz_path)
        assert expr.index.name == "feature_id"
        assert "ILMN_A" in expr.index

    def test_columns_are_gsm_ids(self, tmp_path):
        from scripts.utils.data_io import load_series_matrix_expression
        gz_path = _make_gz_file(SYNTHETIC_SERIES_MATRIX, tmp_path, "test.txt.gz")
        expr = load_series_matrix_expression(gz_path)
        assert list(expr.columns) == ["GSM001", "GSM002", "GSM003"]

    def test_values_are_numeric(self, tmp_path):
        from scripts.utils.data_io import load_series_matrix_expression
        import numpy as np
        gz_path = _make_gz_file(SYNTHETIC_SERIES_MATRIX, tmp_path, "test.txt.gz")
        expr = load_series_matrix_expression(gz_path)
        assert pd.api.types.is_numeric_dtype(expr["GSM001"])

    def test_nrows_limit(self, tmp_path):
        from scripts.utils.data_io import load_series_matrix_expression
        gz_path = _make_gz_file(SYNTHETIC_SERIES_MATRIX, tmp_path, "test.txt.gz")
        expr = load_series_matrix_expression(gz_path, nrows=2)
        assert expr.shape[0] == 2

    def test_missing_sentinel_raises(self, tmp_path):
        """Files without !series_matrix_table_begin must raise ValueError."""
        from scripts.utils.data_io import load_series_matrix_expression
        bad_content = "!Series_title\t\"No table here\"\n"
        gz_path = _make_gz_file(bad_content, tmp_path, "bad.txt.gz")
        with pytest.raises(ValueError, match="sentinel"):
            load_series_matrix_expression(gz_path)


# ---------------------------------------------------------------------------
# Tests: load_tsv_gz
# ---------------------------------------------------------------------------

class TestLoadTsvGz:
    def test_synthetic_counts_load(self, tmp_path):
        from scripts.utils.data_io import load_tsv_gz
        gz_path = _make_gz_file(SYNTHETIC_TSV_GZ_CONTENT, tmp_path, "counts.tsv.gz")
        df = load_tsv_gz(gz_path, index_col=0)
        assert df.shape == (2, 2)
        assert list(df.columns) == ["GSM001", "GSM002"]
        assert df.index.name == "GeneID"

    def test_real_gse125583_counts_header(self):
        """Light smoke test: load only the first 5 rows of the counts file."""
        from scripts.utils.paths import resolve_path
        from scripts.utils.data_io import load_tsv_gz
        path = resolve_path("data/raw/GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz")
        if not path.exists():
            pytest.skip("GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz not found")
        df = load_tsv_gz(path, index_col=0, nrows=5)
        assert df.index.name == "GeneID", (
            f"Expected index named 'GeneID', got '{df.index.name}'"
        )
        # Verify sample columns start with GSM
        assert all(col.startswith("GSM") for col in df.columns), (
            "Count matrix columns should be GSM accession strings."
        )

    def test_real_ncbi_annotation_columns(self):
        """Human.GRCh38.p13.annot.tsv.gz must have expected columns."""
        from scripts.utils.paths import resolve_path
        from scripts.utils.data_io import load_tsv_gz
        path = resolve_path("data/raw/Human.GRCh38.p13.annot.tsv.gz")
        if not path.exists():
            pytest.skip("Human.GRCh38.p13.annot.tsv.gz not found")
        df = load_tsv_gz(path, index_col=None, nrows=3)
        assert "GeneID" in df.columns
        assert "Symbol" in df.columns
        assert "GeneType" in df.columns
        assert "EnsemblGeneID" in df.columns


# ---------------------------------------------------------------------------
# Tests: load_plain_tsv
# ---------------------------------------------------------------------------

class TestLoadPlainTsv:
    def test_synthetic_skips_comments(self, tmp_path):
        from scripts.utils.data_io import load_plain_tsv
        path = _make_plain_file(SYNTHETIC_PLAIN_TSV, tmp_path, "platform.txt")
        df = load_plain_tsv(path)
        # Lines starting with # must be skipped; only 2 data rows remain
        assert len(df) == 2
        assert "ID" in df.columns
        assert "SYMBOL" in df.columns

    def test_real_gpl10558_has_ilmn_ids(self):
        """GPL10558-50081.txt must have ILMN_ probe IDs."""
        from scripts.utils.paths import resolve_path
        from scripts.utils.data_io import load_plain_tsv
        path = resolve_path("data/raw/GPL10558-50081.txt")
        if not path.exists():
            pytest.skip("GPL10558-50081.txt not found")
        df = load_plain_tsv(path, nrows=10)
        assert "ID" in df.columns, "GPL10558 must have an 'ID' column."

    def test_real_gpl15988_four_columns(self):
        """GPL15988.txt (re-downloaded TSV) must have exactly 4 columns."""
        from scripts.utils.paths import resolve_path
        from scripts.utils.data_io import load_plain_tsv
        path = resolve_path("data/raw/GPL15988.txt")
        if not path.exists():
            pytest.skip("GPL15988.txt not found")
        df = load_plain_tsv(path, nrows=5)
        assert len(df.columns) == 4, (
            f"GPL15988.txt should have 4 columns (ID, GB_ACC, GENE SYMBOL, "
            f"GENE DESCRIPTION), but got {len(df.columns)}: {list(df.columns)}"
        )
        assert "ID" in df.columns
        assert "GENE SYMBOL" in df.columns


# ---------------------------------------------------------------------------
# Tests: load_gtf_gz_header
# ---------------------------------------------------------------------------

class TestLoadGtfGzHeader:
    def test_synthetic_header_lines(self, tmp_path):
        from scripts.utils.data_io import load_gtf_gz_header
        gz_path = _make_gz_file(SYNTHETIC_GTF_LINES, tmp_path, "test.gtf.gz")
        lines = load_gtf_gz_header(gz_path)
        assert len(lines) >= 1
        assert all(l.startswith("#") for l in lines)
        assert any("GENCODE" in l or "synthetic" in l.lower() for l in lines)

    def test_real_gencode_lncrna_gtf_header(self):
        """GENCODE lncRNA GTF must start with GENCODE header comments."""
        from scripts.utils.paths import resolve_path
        from scripts.utils.data_io import load_gtf_gz_header
        path = resolve_path("data/raw/gencode.v50.long_noncoding_RNAs.gtf.gz")
        if not path.exists():
            pytest.skip("gencode.v50.long_noncoding_RNAs.gtf.gz not found")
        lines = load_gtf_gz_header(path)
        assert len(lines) >= 1
        full = "\n".join(lines)
        assert "GENCODE" in full or "GRCh38" in full, (
            "GENCODE GTF header should mention GENCODE or GRCh38"
        )
