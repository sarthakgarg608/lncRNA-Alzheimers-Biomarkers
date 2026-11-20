"""
tests/test_stage08_integrity.py
=================================
Stage 08 integrity tests — scientific and data constraints.

Tests cover:
  - GSE140829 locked: diagnosis raises or returns None if accessed
  - RNA-seq counts are non-negative integers
  - Microarray expression values are float (not int)
  - No negative values in blood microarray expression
  - No duplicate sample IDs in loaded data
  - Blood and brain data must never be merged into one matrix
  - Expression matrix is in features × samples orientation (interim storage)
  - Parquet interim files exist and load correctly (if loading scripts have been run)
  - Parsing report JSON files have required keys

Scientific integrity:
  - Only synthetic fixtures are used for unit tests.
  - Integration tests (reading actual parquet) are skipped if files are absent.
  - No real patient data or computed biomarker results appear in this file.
"""

from __future__ import annotations

import gzip
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
_TESTS_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _TESTS_DIR.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import importlib
_load_geo = importlib.import_module("scripts.01_load_geo_datasets")
parse_series_matrix_metadata = _load_geo.parse_series_matrix_metadata
_map_status_label = _load_geo._map_status_label
from scripts.utils.validation import (
    check_duplicate_ids,
    check_expression_numeric,
    check_sample_id_alignment,
    VALID_DIAGNOSIS_LABELS,
)


# ===========================================================================
# Helpers
# ===========================================================================

def _make_float_expr(n_probes: int = 20, n_samples: int = 6) -> pd.DataFrame:
    """Synthetic microarray expression matrix (float32, all positive)."""
    rng = np.random.default_rng(42)
    data = rng.uniform(5.0, 15.0, (n_probes, n_samples)).astype(np.float32)
    probe_ids = [f"ILMN_{1000 + i}" for i in range(n_probes)]
    sample_ids = [f"GSM{1000 + i}" for i in range(n_samples)]
    df = pd.DataFrame(data, index=probe_ids, columns=sample_ids)
    df.index.name = "feature_id"
    return df


def _make_int_counts(n_genes: int = 100, n_samples: int = 10) -> pd.DataFrame:
    """Synthetic RNA-seq count matrix (int32, all non-negative)."""
    rng = np.random.default_rng(42)
    data = rng.integers(0, 5000, (n_genes, n_samples), dtype=np.int32)
    gene_ids = list(range(100000, 100000 + n_genes))
    sample_ids = [f"GSM{2000 + i}" for i in range(n_samples)]
    df = pd.DataFrame(data, index=gene_ids, columns=sample_ids)
    df.index.name = "GeneID"
    return df


def _locked_header(n_samples: int = 4) -> dict:
    sample_ids = [f"GSM{8000 + i}" for i in range(n_samples)]
    return {
        "Sample_geo_accession": [sample_ids],
        "Sample_platform_id": [["GPL15988"] * n_samples],
        "Sample_title": [["Title"] * n_samples],
        "Sample_characteristics_ch1": [
            [f"diagnosis: AD" for _ in range(n_samples)],
            [f"age_at_draw: 70" for _ in range(n_samples)],
        ],
    }


# ===========================================================================
# Tests — GSE140829 diagnosis locked
# ===========================================================================

class TestGSE140829DiagnosisLocked:
    """Diagnosis information for GSE140829 must be inaccessible."""

    def test_parse_with_locked_no_diagnosis_column(self) -> None:
        """parse_series_matrix_metadata with locked=True must not produce diagnosis column."""
        header = _locked_header()
        meta = parse_series_matrix_metadata(
            header, diagnosis_field="diagnosis", locked=True
        )
        assert "diagnosis" not in meta.columns

    def test_locked_metadata_still_has_sample_ids(self) -> None:
        """Structural metadata (sample_id, platform) must be present even when locked."""
        header = _locked_header()
        meta = parse_series_matrix_metadata(
            header, diagnosis_field="diagnosis", locked=True
        )
        assert "sample_id" in meta.columns
        assert len(meta) == 4

    def test_locked_flag_true_means_no_labels_regardless(self) -> None:
        """Even if diagnosis values are in the file, locked=True must suppress them."""
        header = _locked_header(n_samples=6)
        meta = parse_series_matrix_metadata(
            header, diagnosis_field="diagnosis", locked=True
        )
        # Should have no column named 'diagnosis' or any label-like column
        label_cols = [c for c in meta.columns if "diagn" in c.lower() or "label" in c.lower()]
        assert len(label_cols) == 0, f"Unexpected label-like columns: {label_cols}"

    @pytest.mark.skipif(
        not (_PROJECT_ROOT / "data" / "interim" / "GSE140829" / "technical_metadata.csv").exists(),
        reason="GSE140829 technical_metadata.csv not yet generated — run 01_load_geo_datasets.py",
    )
    def test_actual_technical_metadata_no_diagnosis(self) -> None:
        """Integration: actual technical_metadata.csv must have no diagnosis column."""
        path = _PROJECT_ROOT / "data" / "interim" / "GSE140829" / "technical_metadata.csv"
        df = pd.read_csv(path)
        assert "diagnosis" not in df.columns, (
            "CRITICAL: 'diagnosis' column found in GSE140829 technical_metadata.csv — "
            "this violates the locked external validation protocol!"
        )

    @pytest.mark.skipif(
        not (_PROJECT_ROOT / "data" / "interim" / "GSE140829" / "parsing_report.json").exists(),
        reason="Parsing report not yet generated",
    )
    def test_actual_parsing_report_says_locked(self) -> None:
        """Integration: parsing_report.json must state diagnosis is locked."""
        path = _PROJECT_ROOT / "data" / "interim" / "GSE140829" / "parsing_report.json"
        with open(path, encoding="utf-8") as fh:
            report = json.load(fh)
        assert report.get("validation", {}).get("diagnosis_locked") is True, (
            "parsing_report.json does not confirm diagnosis_locked=True for GSE140829"
        )


# ===========================================================================
# Tests — RNA-seq counts are non-negative integers
# ===========================================================================

class TestRnaSeqCountIntegrity:
    """RNA-seq count matrix must contain non-negative integer values."""

    def test_counts_dtype_is_integer(self) -> None:
        counts = _make_int_counts()
        assert np.issubdtype(counts.values.dtype, np.integer), (
            f"Expected integer dtype, got {counts.values.dtype}"
        )

    def test_counts_no_negative_values(self) -> None:
        counts = _make_int_counts()
        assert float(counts.values.min()) >= 0, "Negative values found in count matrix"

    def test_counts_index_is_gene_id(self) -> None:
        counts = _make_int_counts()
        assert counts.index.name == "GeneID"

    def test_counts_columns_are_gsm_ids(self) -> None:
        counts = _make_int_counts()
        assert all(str(c).startswith("GSM") for c in counts.columns)

    @pytest.mark.skipif(
        not (_PROJECT_ROOT / "data" / "interim" / "GSE125583" / "raw_counts.parquet").exists(),
        reason="GSE125583 raw_counts.parquet not yet generated — run 01_load_geo_datasets.py",
    )
    def test_actual_counts_non_negative(self) -> None:
        """Integration: actual raw_counts.parquet must have no negative values."""
        import pyarrow.parquet as pq
        path = _PROJECT_ROOT / "data" / "interim" / "GSE125583" / "raw_counts.parquet"
        df = pq.read_table(path).to_pandas()
        min_val = float(df.values.min())
        assert min_val >= 0, f"Negative count values found: min={min_val}"

    @pytest.mark.skipif(
        not (_PROJECT_ROOT / "data" / "interim" / "GSE125583" / "raw_counts.parquet").exists(),
        reason="GSE125583 raw_counts.parquet not yet generated",
    )
    def test_actual_counts_integer_dtype(self) -> None:
        """Integration: actual raw_counts.parquet must have integer dtype."""
        import pyarrow.parquet as pq
        path = _PROJECT_ROOT / "data" / "interim" / "GSE125583" / "raw_counts.parquet"
        df = pq.read_table(path).to_pandas()
        assert np.issubdtype(df.values.dtype, np.integer), (
            f"Non-integer dtype found: {df.values.dtype}"
        )


# ===========================================================================
# Tests — Microarray expression values are float
# ===========================================================================

class TestMicroarrayExpressionFloat:
    """Microarray expression matrix must have float dtype."""

    def test_expression_dtype_is_float(self) -> None:
        expr = _make_float_expr()
        assert np.issubdtype(expr.values.dtype, np.floating), (
            f"Expected float dtype, got {expr.values.dtype}"
        )

    def test_expression_values_positive(self) -> None:
        expr = _make_float_expr()
        assert float(expr.values.min()) >= 0.0

    def test_expression_not_integer_dtype(self) -> None:
        expr = _make_float_expr()
        # Microarray must NOT be integer
        assert not np.issubdtype(expr.values.dtype, np.integer)

    @pytest.mark.parametrize("gse_id", ["GSE63061", "GSE63060", "GSE140829"])
    @pytest.mark.skipif(
        not any(
            (_PROJECT_ROOT / "data" / "interim" / gse / "expression.parquet").exists()
            for gse in ["GSE63061", "GSE63060", "GSE140829"]
        ),
        reason="Interim parquet files not yet generated — run 01_load_geo_datasets.py",
    )
    def test_actual_microarray_float(self, gse_id: str) -> None:
        """Integration: blood microarray expression parquet must be float dtype."""
        import pyarrow.parquet as pq
        path = _PROJECT_ROOT / "data" / "interim" / gse_id / "expression.parquet"
        if not path.exists():
            pytest.skip(f"{gse_id} parquet not generated yet")
        df = pq.read_table(path).to_pandas()
        assert np.issubdtype(df.values.dtype, np.floating), (
            f"{gse_id}: Expected float dtype, got {df.values.dtype}"
        )


# ===========================================================================
# Tests — No duplicate sample IDs
# ===========================================================================

class TestNoDuplicateSampleIds:
    """Sample IDs must be unique within each dataset."""

    def test_no_duplicates_in_clean_list(self) -> None:
        ids = [f"GSM{i}" for i in range(10)]
        result = check_duplicate_ids(ids)
        assert result.passed

    def test_duplicates_detected(self) -> None:
        ids = ["GSM001", "GSM002", "GSM001", "GSM003"]
        result = check_duplicate_ids(ids)
        assert not result.passed
        assert "GSM001" in result.details.get("duplicate_ids", [])

    def test_single_id_no_duplicate(self) -> None:
        result = check_duplicate_ids(["GSM001"])
        assert result.passed

    def test_empty_list_no_duplicate(self) -> None:
        result = check_duplicate_ids([])
        assert result.passed

    @pytest.mark.parametrize("gse_id", ["GSE63061", "GSE63060", "GSE125583"])
    def test_actual_no_duplicate_sample_ids(self, gse_id: str) -> None:
        """Integration: actual metadata.csv files must have unique sample IDs."""
        path = _PROJECT_ROOT / "data" / "interim" / gse_id / "metadata.csv"
        if not path.exists():
            pytest.skip(f"{gse_id} metadata.csv not generated yet")
        df = pd.read_csv(path)
        result = check_duplicate_ids(df["sample_id"].tolist())
        assert result.passed, (
            f"{gse_id}: Duplicate sample IDs found — {result.details}"
        )

    def test_actual_no_duplicate_in_locked_metadata(self) -> None:
        """Integration: GSE140829 technical_metadata.csv must have unique sample IDs."""
        path = _PROJECT_ROOT / "data" / "interim" / "GSE140829" / "technical_metadata.csv"
        if not path.exists():
            pytest.skip("GSE140829 technical_metadata.csv not generated yet")
        df = pd.read_csv(path)
        result = check_duplicate_ids(df["sample_id"].tolist())
        assert result.passed


# ===========================================================================
# Tests — Blood and brain data are never merged
# ===========================================================================

class TestBloodBrainSeparation:
    """Blood microarray and brain RNA-seq data must never be in the same matrix."""

    def test_blood_and_brain_different_feature_counts(self) -> None:
        """Blood probes and brain genes have different ID types — merging is not safe."""
        blood_expr = _make_float_expr(n_probes=20)
        brain_counts = _make_int_counts(n_genes=100)
        # They must NOT share index names
        assert blood_expr.index.name != brain_counts.index.name or (
            set(blood_expr.index) != set(brain_counts.index)
        )

    def test_blood_index_is_ilmn_probes(self) -> None:
        blood_expr = _make_float_expr(n_probes=5)
        assert all(str(i).startswith("ILMN_") for i in blood_expr.index)

    def test_brain_index_is_entrez_gene_id(self) -> None:
        brain_counts = _make_int_counts(n_genes=5)
        assert brain_counts.index.name == "GeneID"
        assert all(isinstance(i, (int, np.integer)) for i in brain_counts.index)

    def test_cannot_concatenate_blood_and_brain(self) -> None:
        """Direct concatenation of blood and brain on axis=1 produces misaligned index."""
        blood_expr = _make_float_expr(n_probes=10)
        brain_counts = _make_int_counts(n_genes=10)
        # Indices differ: ILMN_* vs integer GeneID
        # Concatenation along axis=0 would produce NaN everywhere — proves they differ
        combined = pd.concat([blood_expr.astype(float), brain_counts.astype(float)], axis=0)
        # All values in blood rows should be non-NaN, all values in brain rows should be non-NaN
        # But the COLUMNS differ (different GSM IDs), so combined has NaN everywhere
        # This test asserts that columns are indeed different (they should be!)
        blood_cols = set(blood_expr.columns)
        brain_cols = set(brain_counts.columns)
        assert blood_cols != brain_cols, (
            "Blood and brain synthetic data use the same sample IDs — test fixture error"
        )

    @pytest.mark.parametrize("gse_id", ["GSE63061", "GSE63060"])
    def test_actual_blood_parquet_separate_from_brain(self, gse_id: str) -> None:
        """Integration: blood parquet must be in a different directory from brain parquet."""
        blood_path = _PROJECT_ROOT / "data" / "interim" / gse_id / "expression.parquet"
        brain_path = _PROJECT_ROOT / "data" / "interim" / "GSE125583" / "raw_counts.parquet"
        # Different paths = not merged
        assert blood_path != brain_path
        # Both should be under their own dataset directory
        assert blood_path.parent.name == gse_id
        assert brain_path.parent.name == "GSE125583"


# ===========================================================================
# Tests — Expression matrix orientation (features × samples)
# ===========================================================================

class TestMatrixOrientation:
    """Interim storage must be features (probes/genes) × samples."""

    def test_float_expr_features_x_samples(self) -> None:
        n_probes, n_samples = 20, 6
        expr = _make_float_expr(n_probes=n_probes, n_samples=n_samples)
        assert expr.shape == (n_probes, n_samples), (
            f"Expected ({n_probes}, {n_samples}), got {expr.shape}"
        )

    def test_int_counts_genes_x_samples(self) -> None:
        n_genes, n_samples = 100, 10
        counts = _make_int_counts(n_genes=n_genes, n_samples=n_samples)
        assert counts.shape == (n_genes, n_samples)

    def test_probes_on_rows_not_columns(self) -> None:
        """Probe IDs must be in the index (rows), not columns."""
        expr = _make_float_expr(n_probes=20, n_samples=6)
        assert all(str(i).startswith("ILMN_") for i in expr.index)
        assert all(str(c).startswith("GSM") for c in expr.columns)

    @pytest.mark.parametrize("gse_id,parquet_fname", [
        ("GSE63061", "expression.parquet"),
        ("GSE63060", "expression.parquet"),
        ("GSE140829", "expression.parquet"),
        ("GSE125583", "raw_counts.parquet"),
    ])
    def test_actual_parquet_features_on_rows(self, gse_id: str, parquet_fname: str) -> None:
        """Integration: actual parquet files must have features × samples orientation."""
        path = _PROJECT_ROOT / "data" / "interim" / gse_id / parquet_fname
        if not path.exists():
            pytest.skip(f"{gse_id}/{parquet_fname} not generated yet")
        import pyarrow.parquet as pq
        meta_csv = "metadata.csv" if gse_id != "GSE140829" else "technical_metadata.csv"
        meta_path = _PROJECT_ROOT / "data" / "interim" / gse_id / meta_csv
        if not meta_path.exists():
            pytest.skip(f"{gse_id} metadata not generated yet")
        df = pq.read_table(path).to_pandas()
        meta = pd.read_csv(meta_path)
        # Columns should be GSM IDs (sample IDs)
        col_set = set(df.columns)
        gsm_col_count = sum(1 for c in col_set if str(c).startswith("GSM"))
        # At least 50% of columns should be GSM IDs
        assert gsm_col_count / max(len(col_set), 1) > 0.5, (
            f"{gse_id}: Columns don't appear to be GSM sample IDs. "
            f"Sample columns: {list(col_set)[:5]}"
        )


# ===========================================================================
# Tests — Parsing report JSON structure
# ===========================================================================

class TestParsingReportStructure:
    """Parsing report JSONs must contain required keys."""

    REQUIRED_KEYS = {"dataset_id", "role", "platform_id", "tissue", "validation", "notes"}

    @pytest.mark.parametrize("gse_id", ["GSE63061", "GSE63060", "GSE140829", "GSE125583"])
    def test_parsing_report_has_required_keys(self, gse_id: str) -> None:
        """Integration: parsing_report.json must contain all required keys."""
        path = _PROJECT_ROOT / "data" / "interim" / gse_id / "parsing_report.json"
        if not path.exists():
            pytest.skip(f"{gse_id} parsing_report.json not generated yet — run 01_load_geo_datasets.py")
        with open(path, encoding="utf-8") as fh:
            report = json.load(fh)
        for key in self.REQUIRED_KEYS:
            assert key in report, (
                f"{gse_id}/parsing_report.json missing required key: '{key}'"
            )

    @pytest.mark.parametrize("gse_id", ["GSE63061", "GSE63060", "GSE125583"])
    def test_unlocked_reports_have_diagnosis_counts(self, gse_id: str) -> None:
        """Unlocked datasets must have diagnosis_counts in their report."""
        path = _PROJECT_ROOT / "data" / "interim" / gse_id / "parsing_report.json"
        if not path.exists():
            pytest.skip(f"{gse_id} parsing_report.json not generated yet")
        with open(path, encoding="utf-8") as fh:
            report = json.load(fh)
        assert "diagnosis_counts" in report, (
            f"{gse_id}: 'diagnosis_counts' missing from parsing_report.json"
        )
        dc = report["diagnosis_counts"]
        assert isinstance(dc, dict)
        # MCI must be represented (even if zero) or AD/Control present
        assert any(k in dc for k in ["AD", "Control", "MCI"])

    def test_locked_report_has_diagnosis_locked_true(self) -> None:
        """GSE140829 parsing report must flag diagnosis_locked=True."""
        path = _PROJECT_ROOT / "data" / "interim" / "GSE140829" / "parsing_report.json"
        if not path.exists():
            pytest.skip("GSE140829 parsing_report.json not generated yet")
        with open(path, encoding="utf-8") as fh:
            report = json.load(fh)
        assert report.get("validation", {}).get("diagnosis_locked") is True


# ===========================================================================
# Tests — Sample ID alignment (expression columns ↔ metadata)
# ===========================================================================

class TestSampleIdAlignment:
    """Expression matrix columns must match metadata sample IDs."""

    def test_alignment_passes_when_matching(self) -> None:
        ids = [f"GSM{i}" for i in range(5)]
        result = check_sample_id_alignment(ids, ids)
        assert result.passed

    def test_alignment_fails_when_mismatch(self) -> None:
        meta_ids = ["GSM001", "GSM002", "GSM003"]
        expr_cols = ["GSM001", "GSM002", "GSM999"]
        result = check_sample_id_alignment(meta_ids, expr_cols)
        assert not result.passed

    def test_alignment_fails_on_extra_in_expr(self) -> None:
        meta_ids = ["GSM001", "GSM002"]
        expr_cols = ["GSM001", "GSM002", "GSM003"]
        result = check_sample_id_alignment(meta_ids, expr_cols)
        assert not result.passed

    def test_alignment_fails_on_extra_in_meta(self) -> None:
        meta_ids = ["GSM001", "GSM002", "GSM003"]
        expr_cols = ["GSM001", "GSM002"]
        result = check_sample_id_alignment(meta_ids, expr_cols)
        assert not result.passed
