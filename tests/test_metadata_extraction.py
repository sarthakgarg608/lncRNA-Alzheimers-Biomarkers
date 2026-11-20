"""
tests/test_metadata_extraction.py
===================================
Unit tests for diagnosis label mapping and metadata standardisation (Stage 08).

Tests cover:
  - Exact diagnosis label mapping (raw GEO strings → canonical vocabulary)
  - MCI preservation (never merged into Control or AD)
  - Case handling of label mapping
  - Unknown assignment for unrecognised labels
  - Canonical schema column presence
  - GSE140829 locked metadata has no diagnosis column
  - Metadata standardisation utility (02_extract_metadata.py)

Scientific integrity:
  - All fixtures are synthetic.
  - No real GEO data is loaded in unit tests.
  - Canonical diagnosis vocabulary: {"AD", "Control", "MCI", "Unknown"} only.
"""

from __future__ import annotations

import sys
from pathlib import Path

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
_map_status_label = _load_geo._map_status_label
parse_series_matrix_metadata = _load_geo.parse_series_matrix_metadata

_extract_meta = importlib.import_module("scripts.02_extract_metadata")
standardise_metadata = _extract_meta.standardise_metadata
CANONICAL_COLUMNS = _extract_meta.CANONICAL_COLUMNS
LOCKED_COLUMNS = _extract_meta.LOCKED_COLUMNS
from scripts.utils.validation import VALID_DIAGNOSIS_LABELS, check_valid_diagnosis_labels


# ===========================================================================
# Tests — exact diagnosis label mapping
# ===========================================================================

class TestDiagnosisLabelMapping:
    """Test _map_status_label for all expected raw GEO string values."""

    # --- AD variants ---
    @pytest.mark.parametrize("raw,expected", [
        ("Alzheimer's disease", "AD"),
        ("alzheimer's disease", "AD"),
        ("ALZHEIMER'S DISEASE", "AD"),
        ("Alzheimers disease", "AD"),
        ("Alzheimer Disease", "AD"),
        ("AD", "AD"),
        ("ad", "AD"),
        ("Alzheimer", "AD"),
    ])
    def test_ad_variants(self, raw: str, expected: str) -> None:
        assert _map_status_label(raw) == expected, f"Input: {repr(raw)}"

    # --- Control variants ---
    @pytest.mark.parametrize("raw,expected", [
        ("Control", "Control"),
        ("control", "Control"),
        ("CONTROL", "Control"),
        ("CTL", "Control"),
        ("ctl", "Control"),
        ("HC", "Control"),
        ("hc", "Control"),
        ("Healthy Control", "Control"),
        ("healthy control", "Control"),
        ("Normal", "Control"),
        ("normal", "Control"),
        ("Healthy", "Control"),
        ("healthy", "Control"),
    ])
    def test_control_variants(self, raw: str, expected: str) -> None:
        assert _map_status_label(raw) == expected, f"Input: {repr(raw)}"

    # --- MCI variants ---
    @pytest.mark.parametrize("raw,expected", [
        ("MCI", "MCI"),
        ("mci", "MCI"),
        ("Mild Cognitive Impairment", "MCI"),
        ("mild cognitive impairment", "MCI"),
    ])
    def test_mci_variants(self, raw: str, expected: str) -> None:
        assert _map_status_label(raw) == expected, f"Input: {repr(raw)}"

    # --- Unknown / unrecognised ---
    @pytest.mark.parametrize("raw", [
        "Unknown",
        "unknown",
        "borderline",
        "mixed",
        "unclear",
        "other",
        "n/a",
        "",
        "   ",
        "Lewy body dementia",
        "Vascular dementia",
        "FTD",
    ])
    def test_unknown_for_unrecognised(self, raw: str) -> None:
        result = _map_status_label(raw)
        assert result == "Unknown", f"Expected 'Unknown' for {repr(raw)}, got {repr(result)}"

    def test_mci_is_not_control(self) -> None:
        """MCI must map to 'MCI', not 'Control'."""
        assert _map_status_label("MCI") != "Control"
        assert _map_status_label("Mild Cognitive Impairment") != "Control"

    def test_ad_is_not_control(self) -> None:
        """AD must not map to Control."""
        assert _map_status_label("Alzheimer's disease") != "Control"

    def test_all_outputs_in_canonical_vocab(self) -> None:
        """All possible outputs of _map_status_label are in canonical vocabulary."""
        test_inputs = [
            "Alzheimer's disease", "Control", "MCI", "borderline", "unknown", "healthy",
        ]
        for raw in test_inputs:
            result = _map_status_label(raw)
            assert result in VALID_DIAGNOSIS_LABELS, (
                f"Output {repr(result)} for input {repr(raw)} is not in canonical vocabulary"
            )

    def test_exact_geo_value_alzheimers_disease(self) -> None:
        """The exact string seen in GSE63061/63060 must map to 'AD'."""
        assert _map_status_label("Alzheimer's disease") == "AD"

    def test_exact_geo_value_control(self) -> None:
        """The exact string 'Control' seen in GSE63061/63060 must map to 'Control'."""
        assert _map_status_label("Control") == "Control"

    def test_exact_geo_value_mci(self) -> None:
        """The exact string 'MCI' seen in GSE63061/63060 must map to 'MCI'."""
        assert _map_status_label("MCI") == "MCI"

    def test_whitespace_stripped(self) -> None:
        """Leading/trailing whitespace in raw value must be ignored."""
        assert _map_status_label("  Alzheimer's disease  ") == "AD"
        assert _map_status_label("  Control  ") == "Control"
        assert _map_status_label("  MCI  ") == "MCI"


# ===========================================================================
# Tests — parse_series_matrix_metadata output schema
# ===========================================================================

class TestMetadataSchema:
    """Test that parse_series_matrix_metadata produces the expected schema."""

    def _make_header(
        self,
        n_samples: int = 4,
        statuses: list[str] | None = None,
    ) -> dict:
        """Build a minimal synthetic header dict."""
        if statuses is None:
            statuses = ["Alzheimer's disease"] * 2 + ["Control"] * 2
        sample_ids = [f"GSM{9000 + i}" for i in range(n_samples)]
        char_lines: list[list[str]] = []
        # status row
        char_lines.append([f"status: {s}" for s in statuses])
        # age row
        char_lines.append([f"age: {60 + i}" for i in range(n_samples)])
        # gender row
        char_lines.append(["gender: Female"] * n_samples)
        # tissue row
        char_lines.append(["tissue: blood"] * n_samples)

        return {
            "Sample_geo_accession": [sample_ids],
            "Sample_platform_id": [["GPL10558"] * n_samples],
            "Sample_title": [["Title"] * n_samples],
            "Sample_characteristics_ch1": char_lines,
        }

    def test_sample_id_column_present(self) -> None:
        header = self._make_header()
        meta = parse_series_matrix_metadata(header)
        assert "sample_id" in meta.columns

    def test_diagnosis_column_present_when_unlocked(self) -> None:
        header = self._make_header()
        meta = parse_series_matrix_metadata(header, locked=False)
        assert "diagnosis" in meta.columns

    def test_diagnosis_column_absent_when_locked(self) -> None:
        """When locked=True, no 'diagnosis' column must appear."""
        header = self._make_header()
        meta = parse_series_matrix_metadata(header, locked=True)
        assert "diagnosis" not in meta.columns

    def test_correct_number_of_rows(self) -> None:
        header = self._make_header(n_samples=6)
        meta = parse_series_matrix_metadata(header)
        assert len(meta) == 6

    def test_diagnosis_values_correct(self) -> None:
        """Diagnosis values in the DataFrame must match the mapped canonical labels."""
        header = self._make_header(
            n_samples=4,
            statuses=["Alzheimer's disease", "Control", "MCI", "borderline"],
        )
        meta = parse_series_matrix_metadata(header, diagnosis_field="status")
        expected = ["AD", "Control", "MCI", "Unknown"]
        assert list(meta["diagnosis"]) == expected, (
            f"Expected {expected}, got {list(meta['diagnosis'])}"
        )

    def test_mci_preserved_in_schema(self) -> None:
        header = self._make_header(statuses=["MCI"] * 2 + ["Control"] * 2)
        meta = parse_series_matrix_metadata(header, diagnosis_field="status")
        mci_count = (meta["diagnosis"] == "MCI").sum()
        assert mci_count == 2

    def test_sample_ids_start_with_gsm(self) -> None:
        header = self._make_header(n_samples=4)
        meta = parse_series_matrix_metadata(header)
        assert all(meta["sample_id"].str.startswith("GSM"))


# ===========================================================================
# Tests — standardise_metadata (02_extract_metadata)
# ===========================================================================

class TestStandardiseMetadata:
    """Test metadata standardisation to canonical column schema."""

    def _make_raw_meta(
        self,
        n: int = 5,
        with_diagnosis: bool = True,
        locked: bool = False,
    ) -> pd.DataFrame:
        df = pd.DataFrame({
            "sample_id": [f"GSM{i}" for i in range(n)],
            "dataset_id": ["GSE99999"] * n,
            "age": [str(60 + i) for i in range(n)],
            "sex": ["Female"] * n,
            "tissue": ["blood"] * n,
            "platform_id": ["GPL10558"] * n,
            "batch": [None] * n,
        })
        if with_diagnosis and not locked:
            df["diagnosis"] = ["AD"] * 2 + ["Control"] * 2 + ["MCI"] * 1
        return df

    def test_canonical_columns_present_unlocked(self) -> None:
        raw = self._make_raw_meta(with_diagnosis=True)
        std = standardise_metadata(raw, "GSE99999", locked=False)
        for col in CANONICAL_COLUMNS:
            assert col in std.columns, f"Missing canonical column: {col}"

    def test_no_diagnosis_column_for_locked(self) -> None:
        raw = self._make_raw_meta(with_diagnosis=False, locked=True)
        std = standardise_metadata(raw, "GSE140829", locked=True)
        assert "diagnosis" not in std.columns

    def test_missing_optional_column_filled_with_none(self) -> None:
        """If 'ethnicity' is absent, it should be absent (not in schema)."""
        raw = self._make_raw_meta(with_diagnosis=True)
        # Remove 'batch' to test graceful handling
        raw = raw.drop(columns=["batch"])
        std = standardise_metadata(raw, "GSE99999", locked=False)
        # batch should be re-added as None column
        assert "batch" in std.columns
        assert std["batch"].isna().all()

    def test_dataset_id_added_if_missing(self) -> None:
        raw = self._make_raw_meta()
        raw = raw.drop(columns=["dataset_id"])
        std = standardise_metadata(raw, "GSE77777", locked=False)
        assert "dataset_id" in std.columns
        assert (std["dataset_id"] == "GSE77777").all()

    def test_valid_diagnosis_labels_pass(self) -> None:
        """Labels from canonical vocabulary should not trigger warnings."""
        raw = self._make_raw_meta(with_diagnosis=True)
        # This should not raise
        std = standardise_metadata(raw, "GSE99999", locked=False)
        result = check_valid_diagnosis_labels(std["diagnosis"].dropna().tolist())
        assert result.passed

    def test_invalid_diagnosis_label_reported(self) -> None:
        """Labels outside canonical vocabulary should fail the validation check."""
        raw = self._make_raw_meta(with_diagnosis=True)
        raw["diagnosis"] = ["AD", "Control", "borderline", "MCI", "Unknown"]
        std = standardise_metadata(raw, "GSE99999", locked=False)
        result = check_valid_diagnosis_labels(std["diagnosis"].dropna().tolist())
        assert not result.passed


# ===========================================================================
# Tests — GSE140829 locked diagnosis
# ===========================================================================

class TestGSE140829Locked:
    """Specific tests ensuring GSE140829 diagnosis is never accessible."""

    def _make_locked_header(self, n_samples: int = 4) -> dict:
        sample_ids = [f"GSM{8000 + i}" for i in range(n_samples)]
        return {
            "Sample_geo_accession": [sample_ids],
            "Sample_platform_id": [["GPL15988"] * n_samples],
            "Sample_title": [["Title"] * n_samples],
            "Sample_characteristics_ch1": [
                [f"diagnosis: {'AD' if i < 2 else 'Control'}" for i in range(n_samples)],
                [f"age_at_draw: {65 + i}" for i in range(n_samples)],
                [f"Sex: {'Male' if i % 2 == 0 else 'Female'}" for i in range(n_samples)],
                [f"batch: {i // 2 + 1}" for i in range(n_samples)],
            ],
        }

    def test_locked_metadata_has_no_diagnosis_column(self) -> None:
        header = self._make_locked_header()
        meta = parse_series_matrix_metadata(
            header,
            diagnosis_field="diagnosis",
            age_field="age_at_draw",
            sex_field="Sex",
            batch_field="batch",
            locked=True,
        )
        assert "diagnosis" not in meta.columns, (
            "SAFETY VIOLATION: 'diagnosis' column found in locked GSE140829 metadata!"
        )

    def test_locked_metadata_has_structural_fields(self) -> None:
        header = self._make_locked_header()
        meta = parse_series_matrix_metadata(
            header,
            diagnosis_field="diagnosis",
            age_field="age_at_draw",
            sex_field="Sex",
            batch_field="batch",
            locked=True,
        )
        assert "sample_id" in meta.columns
        assert "age" in meta.columns
        assert "sex" in meta.columns
        assert "batch" in meta.columns

    def test_locked_standardise_no_diagnosis(self) -> None:
        """standardise_metadata with locked=True must not produce diagnosis column."""
        raw = pd.DataFrame({
            "sample_id": [f"GSM{i}" for i in range(4)],
            "dataset_id": ["GSE140829"] * 4,
            "age": ["70", "72", "68", "75"],
            "sex": ["Male", "Female", "Male", "Female"],
            "platform_id": ["GPL15988"] * 4,
            "batch": ["1", "1", "2", "2"],
        })
        std = standardise_metadata(raw, "GSE140829", locked=True)
        assert "diagnosis" not in std.columns
