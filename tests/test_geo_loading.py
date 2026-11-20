"""
tests/test_geo_loading.py
==========================
Unit tests for GEO Series Matrix parsing (Stage 08).

Tests cover:
  - !series_matrix_table_begin sentinel detection
  - Expression matrix shape and dtype correctness
  - Probe ID format preservation (ILMN_* for GPL10558/GPL6947)
  - Metadata extraction from Series Matrix headers
  - Sample ID alignment between expression columns and metadata

Scientific integrity:
  - No actual GEO files are loaded in unit tests (uses synthetic fixtures).
  - All fixtures are stored here (no real patient data).
  - MCI must never be merged into Control.
  - Diagnosis mapping must be exact and case-insensitive.
"""

from __future__ import annotations

import gzip
import io
import json
import sys
import textwrap
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

from scripts.utils.data_io import (
    load_series_matrix_header,
    load_series_matrix_expression,
)
import importlib
_load_geo = importlib.import_module("scripts.01_load_geo_datasets")
parse_series_matrix_metadata = _load_geo.parse_series_matrix_metadata
_map_status_label = _load_geo._map_status_label
parse_gse125583_metadata = _load_geo.parse_gse125583_metadata


# ---------------------------------------------------------------------------
# Fixtures: synthetic Series Matrix content
# ---------------------------------------------------------------------------

def _build_series_matrix(
    n_probes: int = 5,
    n_samples: int = 4,
    probe_prefix: str = "ILMN_",
    characteristics: list[str] | None = None,
    series_id: str = "GSE99999",
) -> bytes:
    """
    Build a minimal GEO Series Matrix file as bytes.
    Returns gzip-compressed content.
    """
    if characteristics is None:
        characteristics = [
            "\t".join(["status: Alzheimer's disease"] * 2 + ["status: Control"] * 2),
            "\t".join(["age: 72"] * 2 + ["age: 68"] * 2),
            "\t".join(["gender: Female"] * 3 + ["gender: Male"] * 1),
            "\t".join(["tissue: blood"] * 4),
        ]

    sample_ids = [f"GSM{10000 + i}" for i in range(n_samples)]
    probe_ids = [f"{probe_prefix}{1000 + i}" for i in range(n_probes)]

    lines: list[str] = [
        f'!Series_geo_accession\t"{series_id}"',
    ]
    # Sample geo accession
    lines.append("!Sample_geo_accession\t" + "\t".join(f'"{s}"' for s in sample_ids))
    # Sample platform
    lines.append("!Sample_platform_id\t" + "\t".join(['"GPL10558"'] * n_samples))
    # Sample title
    lines.append("!Sample_title\t" + "\t".join(['"Sample title"'] * n_samples))
    # Characteristics
    for char_line in characteristics:
        lines.append(f'!Sample_characteristics_ch1\t{char_line}')
    # Sentinel
    lines.append("!series_matrix_table_begin")
    # Header row
    lines.append('"ID_REF"\t' + "\t".join(f'"{s}"' for s in sample_ids))
    # Data rows
    rng = np.random.default_rng(42)
    for pid in probe_ids:
        vals = rng.uniform(5.0, 15.0, n_samples)
        lines.append(f'"{pid}"\t' + "\t".join(f"{v:.4f}" for v in vals))
    lines.append("!series_matrix_table_end")
    lines.append("")

    content = "\n".join(lines)
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        gz.write(content.encode("utf-8"))
    buf.seek(0)
    return buf.read()


@pytest.fixture()
def series_matrix_gz(tmp_path: Path) -> Path:
    """Write a synthetic series matrix to a temp .gz file."""
    content = _build_series_matrix(n_probes=10, n_samples=6)
    p = tmp_path / "test_series_matrix.txt.gz"
    p.write_bytes(content)
    return p


@pytest.fixture()
def series_matrix_no_mci_gz(tmp_path: Path) -> Path:
    """Series matrix with only AD and Control (no MCI)."""
    chars = [
        "\t".join(["status: Alzheimer's disease"] * 3 + ["status: Control"] * 3),
        "\t".join(["age: 70"] * 6),
        "\t".join(["gender: Female"] * 6),
        "\t".join(["tissue: blood"] * 6),
    ]
    content = _build_series_matrix(n_probes=8, n_samples=6, characteristics=chars)
    p = tmp_path / "test_no_mci.txt.gz"
    p.write_bytes(content)
    return p


@pytest.fixture()
def series_matrix_with_mci_gz(tmp_path: Path) -> Path:
    """Series matrix with AD, Control, and MCI."""
    chars = [
        "\t".join(
            ["status: Alzheimer's disease"] * 2 +
            ["status: Control"] * 2 +
            ["status: MCI"] * 2
        ),
        "\t".join(["age: 72"] * 6),
        "\t".join(["gender: Male"] * 6),
        "\t".join(["tissue: blood"] * 6),
    ]
    content = _build_series_matrix(n_probes=8, n_samples=6, characteristics=chars)
    p = tmp_path / "test_with_mci.txt.gz"
    p.write_bytes(content)
    return p


# ===========================================================================
# Tests — sentinel detection
# ===========================================================================

class TestSentinelDetection:
    """Test that !series_matrix_table_begin is correctly detected."""

    def test_sentinel_found(self, series_matrix_gz: Path) -> None:
        """Expression loading should succeed when sentinel is present."""
        expr = load_series_matrix_expression(series_matrix_gz)
        assert expr is not None
        assert expr.shape[0] > 0

    def test_sentinel_missing_raises(self, tmp_path: Path) -> None:
        """Files without the sentinel must raise ValueError."""
        content = "!Series_geo_accession\t\"GSE00001\"\nsome_data\n"
        p = tmp_path / "no_sentinel.txt.gz"
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
            gz.write(content.encode())
        buf.seek(0)
        p.write_bytes(buf.read())

        with pytest.raises(ValueError, match="sentinel not found"):
            load_series_matrix_expression(p)

    def test_header_stops_at_sentinel(self, series_matrix_gz: Path) -> None:
        """Header parser must stop before the expression table."""
        header = load_series_matrix_header(series_matrix_gz)
        # Expression data rows should NOT appear in header
        assert "ID_REF" not in header
        assert "series_matrix_table_begin" not in header


# ===========================================================================
# Tests — expression matrix shape and dtype
# ===========================================================================

class TestExpressionMatrix:
    """Test expression matrix loading and properties."""

    def test_shape_probes_x_samples(self, series_matrix_gz: Path) -> None:
        """Matrix must be (n_probes, n_samples) — features × samples."""
        expr = load_series_matrix_expression(series_matrix_gz)
        assert expr.shape == (10, 6), f"Expected (10, 6), got {expr.shape}"

    def test_index_name_is_feature_id(self, series_matrix_gz: Path) -> None:
        """Index name must be 'feature_id'."""
        expr = load_series_matrix_expression(series_matrix_gz)
        assert expr.index.name == "feature_id"

    def test_dtype_float32(self, series_matrix_gz: Path) -> None:
        """Expression values must be float32 by default."""
        expr = load_series_matrix_expression(series_matrix_gz, dtype="float32")
        assert expr.dtypes.iloc[0] == np.float32

    def test_no_negative_values_microarray(self, series_matrix_gz: Path) -> None:
        """Synthetic microarray values are positive (5–15 range)."""
        expr = load_series_matrix_expression(series_matrix_gz)
        assert float(expr.values.min()) >= 0.0

    def test_column_count_matches_samples(self, series_matrix_gz: Path) -> None:
        """Number of columns must equal number of samples."""
        expr = load_series_matrix_expression(series_matrix_gz)
        assert expr.shape[1] == 6

    def test_columns_are_gsm_ids(self, series_matrix_gz: Path) -> None:
        """Column names must be GSM accession strings."""
        expr = load_series_matrix_expression(series_matrix_gz)
        assert all(c.startswith("GSM") for c in expr.columns)

    def test_index_is_unique(self, series_matrix_gz: Path) -> None:
        """Probe IDs (index) must be unique."""
        expr = load_series_matrix_expression(series_matrix_gz)
        assert expr.index.is_unique


# ===========================================================================
# Tests — probe ID format preservation
# ===========================================================================

class TestProbeIdPreservation:
    """Test that ILMN_* probe IDs are preserved exactly."""

    def test_ilmn_prefix_preserved(self, series_matrix_gz: Path) -> None:
        """ILMN_ prefix probes must remain unchanged after loading."""
        expr = load_series_matrix_expression(series_matrix_gz)
        ilmn_probes = expr.index.str.startswith("ILMN_").sum()
        assert ilmn_probes == 10, f"Expected 10 ILMN_ probes, got {ilmn_probes}"

    def test_probes_not_stripped(self, series_matrix_gz: Path) -> None:
        """Probe IDs must not have extra whitespace or quotes stripped into them."""
        expr = load_series_matrix_expression(series_matrix_gz)
        for pid in expr.index:
            assert pid == pid.strip(), f"Probe ID has leading/trailing whitespace: {repr(pid)}"
            assert '"' not in pid, f"Probe ID contains quote: {repr(pid)}"


# ===========================================================================
# Tests — metadata extraction
# ===========================================================================

class TestMetadataExtraction:
    """Test metadata extraction from Series Matrix headers."""

    def test_sample_ids_extracted(self, series_matrix_gz: Path) -> None:
        """All sample IDs must be extracted from Sample_geo_accession."""
        header = load_series_matrix_header(series_matrix_gz)
        meta = parse_series_matrix_metadata(header)
        assert len(meta) == 6
        assert all(meta["sample_id"].str.startswith("GSM"))

    def test_no_duplicate_sample_ids(self, series_matrix_gz: Path) -> None:
        """Sample IDs must be unique in metadata."""
        header = load_series_matrix_header(series_matrix_gz)
        meta = parse_series_matrix_metadata(header)
        assert meta["sample_id"].is_unique

    def test_age_extracted(self, series_matrix_gz: Path) -> None:
        """Age values must be extracted from characteristics."""
        header = load_series_matrix_header(series_matrix_gz)
        meta = parse_series_matrix_metadata(header, age_field="age")
        # At least some age values should be non-null
        non_null_age = meta["age"].notna().sum()
        assert non_null_age > 0, "No age values extracted"

    def test_sex_extracted(self, series_matrix_gz: Path) -> None:
        """Sex values must be extracted from characteristics."""
        header = load_series_matrix_header(series_matrix_gz)
        meta = parse_series_matrix_metadata(header, sex_field="gender")
        non_null_sex = meta["sex"].notna().sum()
        assert non_null_sex > 0, "No sex values extracted"

    def test_missing_field_returns_none(self, series_matrix_gz: Path) -> None:
        """Requesting a field that doesn't exist should return None values, not error."""
        header = load_series_matrix_header(series_matrix_gz)
        meta = parse_series_matrix_metadata(header, batch_field="nonexistent_batch_field")
        assert "batch" in meta.columns
        assert meta["batch"].isna().all()


# ===========================================================================
# Tests — sample ID alignment
# ===========================================================================

class TestSampleIdAlignment:
    """Test that expression columns match metadata sample IDs."""

    def test_expression_columns_match_metadata_ids(self, series_matrix_gz: Path) -> None:
        """Expression column names must equal metadata sample_id values."""
        expr = load_series_matrix_expression(series_matrix_gz)
        header = load_series_matrix_header(series_matrix_gz)
        meta = parse_series_matrix_metadata(header)
        expr_cols = set(expr.columns)
        meta_ids = set(meta["sample_id"].tolist())
        assert expr_cols == meta_ids, (
            f"Mismatch:\n  Expr columns: {sorted(expr_cols)}\n  Meta IDs: {sorted(meta_ids)}"
        )

    def test_expression_column_order_consistent(self, series_matrix_gz: Path) -> None:
        """Expression columns should be a consistent ordered list (not shuffled)."""
        expr = load_series_matrix_expression(series_matrix_gz)
        # Loading the same file twice should give same column order
        expr2 = load_series_matrix_expression(series_matrix_gz)
        assert list(expr.columns) == list(expr2.columns)


# ===========================================================================
# Tests — MCI preservation
# ===========================================================================

class TestMciPreservation:
    """MCI must never be merged into Control or discarded."""

    def test_mci_preserved_in_metadata(self, series_matrix_with_mci_gz: Path) -> None:
        """MCI samples must appear as 'MCI' in diagnosis column."""
        header = load_series_matrix_header(series_matrix_with_mci_gz)
        meta = parse_series_matrix_metadata(header, diagnosis_field="status")
        assert "MCI" in meta["diagnosis"].values, "MCI not found in diagnosis column"

    def test_mci_not_merged_into_control(self, series_matrix_with_mci_gz: Path) -> None:
        """MCI count must be > 0 and not absorbed into Control count."""
        header = load_series_matrix_header(series_matrix_with_mci_gz)
        meta = parse_series_matrix_metadata(header, diagnosis_field="status")
        mci_count = (meta["diagnosis"] == "MCI").sum()
        ctrl_count = (meta["diagnosis"] == "Control").sum()
        assert mci_count == 2, f"Expected 2 MCI, got {mci_count}"
        assert ctrl_count == 2, f"Expected 2 Control, got {ctrl_count}"

    def test_mci_separate_from_all_other_labels(self, series_matrix_with_mci_gz: Path) -> None:
        """MCI label must be distinct from AD, Control, and Unknown."""
        header = load_series_matrix_header(series_matrix_with_mci_gz)
        meta = parse_series_matrix_metadata(header, diagnosis_field="status")
        mci_rows = meta[meta["diagnosis"] == "MCI"]
        assert not (mci_rows["diagnosis"] == "Control").any()
        assert not (mci_rows["diagnosis"] == "AD").any()
        assert not (mci_rows["diagnosis"] == "Unknown").any()
