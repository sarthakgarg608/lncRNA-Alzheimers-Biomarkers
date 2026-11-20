"""
tests/test_validation.py
========================
Tests for scripts/utils/validation.py

All test data is SYNTHETIC.  These tests verify software correctness only.
No results from these tests should ever be treated as biological evidence.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scripts.utils.validation import (
    ValidationResult,
    check_duplicate_ids,
    check_expression_dimensions,
    check_expression_numeric,
    check_expression_scale_flag,
    check_matrix_orientation,
    check_missing_columns,
    check_sample_id_alignment,
    check_valid_diagnosis_labels,
    REQUIRED_METADATA_COLUMNS,
    VALID_DIAGNOSIS_LABELS,
)


# ---------------------------------------------------------------------------
# ValidationResult dataclass
# ---------------------------------------------------------------------------

class TestValidationResult:
    def test_str_pass(self):
        r = ValidationResult("test", True, "all good")
        assert "[PASS]" in str(r)

    def test_str_fail(self):
        r = ValidationResult("test", False, "something wrong")
        assert "[FAIL]" in str(r)


# ---------------------------------------------------------------------------
# check_duplicate_ids
# ---------------------------------------------------------------------------

class TestCheckDuplicateIds:
    def test_no_duplicates_passes(self):
        result = check_duplicate_ids(["GSM001", "GSM002", "GSM003"])
        assert result.passed

    def test_duplicates_fails(self):
        result = check_duplicate_ids(["GSM001", "GSM002", "GSM001"])
        assert not result.passed
        assert "GSM001" in str(result.details)

    def test_empty_list_passes(self):
        result = check_duplicate_ids([])
        assert result.passed

    def test_all_unique_large_list(self):
        ids = [f"GSM{i:06d}" for i in range(1000)]
        result = check_duplicate_ids(ids)
        assert result.passed

    def test_reports_all_duplicates(self):
        ids = ["A", "B", "A", "C", "B"]
        result = check_duplicate_ids(ids)
        assert not result.passed
        dups = result.details["duplicate_ids"]
        assert "A" in dups
        assert "B" in dups


# ---------------------------------------------------------------------------
# check_sample_id_alignment
# ---------------------------------------------------------------------------

class TestCheckSampleIdAlignment:
    def test_exact_match_passes(self):
        meta = ["GSM001", "GSM002", "GSM003"]
        expr = ["GSM001", "GSM002", "GSM003"]
        result = check_sample_id_alignment(meta, expr)
        assert result.passed

    def test_mismatch_fails(self):
        meta = ["GSM001", "GSM002", "GSM004"]
        expr = ["GSM001", "GSM002", "GSM003"]
        result = check_sample_id_alignment(meta, expr)
        assert not result.passed

    def test_subset_allowed(self):
        meta = ["GSM001", "GSM002"]
        expr = ["GSM001", "GSM002", "GSM003"]
        result = check_sample_id_alignment(meta, expr, allow_subset=True)
        assert result.passed

    def test_extra_in_meta_fails_even_subset_allowed(self):
        meta = ["GSM001", "GSM002", "GSM999"]
        expr = ["GSM001", "GSM002"]
        result = check_sample_id_alignment(meta, expr, allow_subset=True)
        assert not result.passed

    def test_order_does_not_matter(self):
        meta = ["GSM003", "GSM001", "GSM002"]
        expr = ["GSM001", "GSM002", "GSM003"]
        result = check_sample_id_alignment(meta, expr)
        assert result.passed


# ---------------------------------------------------------------------------
# check_missing_columns
# ---------------------------------------------------------------------------

class TestCheckMissingColumns:
    def test_all_present_passes(self):
        df = pd.DataFrame(columns=REQUIRED_METADATA_COLUMNS)
        result = check_missing_columns(df)
        assert result.passed

    def test_missing_column_fails(self):
        df = pd.DataFrame(columns=["sample_id", "diagnosis"])
        result = check_missing_columns(df, required=["sample_id", "diagnosis", "tissue"])
        assert not result.passed
        assert "tissue" in result.details["missing_columns"]

    def test_extra_columns_ok(self):
        cols = REQUIRED_METADATA_COLUMNS + ["extra_column", "another_col"]
        df = pd.DataFrame(columns=cols)
        result = check_missing_columns(df)
        assert result.passed

    def test_empty_required_always_passes(self):
        df = pd.DataFrame()
        result = check_missing_columns(df, required=[])
        assert result.passed


# ---------------------------------------------------------------------------
# check_valid_diagnosis_labels
# ---------------------------------------------------------------------------

class TestCheckValidDiagnosisLabels:
    def test_all_valid_passes(self):
        labels = ["AD", "Control", "MCI", "Unknown", "AD", "Control"]
        result = check_valid_diagnosis_labels(labels)
        assert result.passed

    def test_invalid_label_fails(self):
        labels = ["AD", "Control", "Healthy"]  # "Healthy" not in vocabulary
        result = check_valid_diagnosis_labels(labels)
        assert not result.passed
        assert "Healthy" in str(result.details)

    def test_mci_does_not_fail(self):
        """MCI is a valid label — it must NOT be reclassified."""
        labels = ["MCI", "MCI", "Control"]
        result = check_valid_diagnosis_labels(labels)
        assert result.passed

    def test_null_values_flagged_but_no_fail(self):
        """Null/NaN labels should be flagged in message but not cause FAIL."""
        labels = ["AD", None, "Control"]
        result = check_valid_diagnosis_labels(labels)
        # None values are not from illegal vocabulary, so should PASS
        assert result.passed
        assert "null" in result.message.lower() or "nan" in result.message.lower() or "unknown" in result.message.lower()

    def test_empty_labels_passes(self):
        result = check_valid_diagnosis_labels([])
        assert result.passed

    def test_case_sensitive(self):
        """Labels are case-sensitive: 'control' (lower) is NOT valid."""
        labels = ["ad", "control"]
        result = check_valid_diagnosis_labels(labels)
        assert not result.passed


# ---------------------------------------------------------------------------
# check_expression_numeric
# ---------------------------------------------------------------------------

class TestCheckExpressionNumeric:
    def _make_expr(self, data) -> pd.DataFrame:
        return pd.DataFrame(data, index=["ILMN_A", "ILMN_B"], columns=["GSM001", "GSM002"])

    def test_all_numeric_passes(self):
        df = self._make_expr([[10.0, 11.0], [20.0, 21.0]])
        result = check_expression_numeric(df)
        assert result.passed

    def test_non_numeric_column_fails(self):
        df = pd.DataFrame(
            {"GSM001": [10.0, 20.0], "GSM002": ["bad", "data"]},
            index=["ILMN_A", "ILMN_B"]
        )
        result = check_expression_numeric(df)
        assert not result.passed
        assert "GSM002" in str(result.details)

    def test_high_nan_fraction_fails(self):
        data = [[float("nan")] * 2] * 2
        df = self._make_expr(data)
        result = check_expression_numeric(df, max_nan_fraction=0.1)
        assert not result.passed

    def test_low_nan_fraction_passes(self):
        df = self._make_expr([[10.0, float("nan")], [20.0, 21.0]])
        result = check_expression_numeric(df, max_nan_fraction=0.5)
        assert result.passed

    def test_empty_dataframe_fails(self):
        df = pd.DataFrame()
        result = check_expression_numeric(df)
        assert not result.passed


# ---------------------------------------------------------------------------
# check_expression_dimensions
# ---------------------------------------------------------------------------

class TestCheckExpressionDimensions:
    def test_ok_dimensions_pass(self):
        df = pd.DataFrame(
            np.random.rand(500, 50),
            index=[f"ILMN_{i}" for i in range(500)],
            columns=[f"GSM{i:04d}" for i in range(50)]
        )
        result = check_expression_dimensions(df, min_features=100, min_samples=10)
        assert result.passed

    def test_too_few_features_fails(self):
        df = pd.DataFrame(np.random.rand(5, 20))
        result = check_expression_dimensions(df, min_features=100, min_samples=10)
        assert not result.passed

    def test_too_few_samples_fails(self):
        df = pd.DataFrame(np.random.rand(500, 3))
        result = check_expression_dimensions(df, min_features=100, min_samples=10)
        assert not result.passed


# ---------------------------------------------------------------------------
# check_matrix_orientation
# ---------------------------------------------------------------------------

class TestCheckMatrixOrientation:
    def test_features_x_samples_correct(self):
        df = pd.DataFrame(
            {"GSM001": [1.0], "GSM002": [2.0]},
            index=pd.Index(["ILMN_A"], name="feature_id")
        )
        result = check_matrix_orientation(df, expected="features_x_samples")
        assert result.passed

    def test_gsm_in_index_flagged_as_transposed(self):
        """If index contains GSM IDs, the matrix looks transposed."""
        df = pd.DataFrame(
            {"ILMN_A": [1.0], "ILMN_B": [2.0]},
            index=pd.Index(["GSM001"], name="sample_id")
        )
        result = check_matrix_orientation(df, expected="features_x_samples")
        assert not result.passed

    def test_samples_x_features_correct(self):
        df = pd.DataFrame(
            {"ILMN_A": [1.0, 2.0], "ILMN_B": [3.0, 4.0]},
            index=pd.Index(["GSM001", "GSM002"], name="sample_id")
        )
        result = check_matrix_orientation(df, expected="samples_x_features")
        assert result.passed


# ---------------------------------------------------------------------------
# check_expression_scale_flag
# ---------------------------------------------------------------------------

class TestCheckExpressionScaleFlag:
    def test_normalized_microarray_passes(self):
        cfg = {"expression_scale": "normalized_microarray"}
        result = check_expression_scale_flag(cfg)
        assert result.passed
        assert "limma" in result.message or "statsmodels" in result.message

    def test_raw_integer_counts_passes(self):
        cfg = {"expression_scale": "raw_integer_counts"}
        result = check_expression_scale_flag(cfg)
        assert result.passed
        assert "DESeq2" in result.message or "PyDESeq2" in result.message

    def test_missing_field_fails(self):
        cfg = {}
        result = check_expression_scale_flag(cfg)
        assert not result.passed

    def test_unknown_scale_fails(self):
        cfg = {"expression_scale": "unknown_magic_scale"}
        result = check_expression_scale_flag(cfg)
        assert not result.passed

    def test_dataset_config_scales_all_valid(self):
        """All four datasets in config must have valid expression_scale."""
        from scripts.utils.paths import load_dataset_config
        cfg = load_dataset_config()
        for ds_id, ds_cfg in cfg["datasets"].items():
            result = check_expression_scale_flag(ds_cfg)
            assert result.passed, (
                f"Dataset '{ds_id}' has invalid expression_scale. "
                f"Validation message: {result.message}"
            )
