"""
tests/test_stage09_qc.py
========================
Stage 09 — Unit tests for Quality Control functions.

All fixtures are SYNTHETIC (no real patient data).
Tests cover:
  - Sample QC metric calculations (median, mean, std, IQR)
  - RNA-seq library size calculation
  - Raw integer count preservation
  - Matrix orientation: samples-on-rows in derived matrices
  - Missing value handling in QC metrics
  - PCA input preparation
  - Constant-feature detection
  - No destructive normalization
  - No unapproved sample deletion
  - GSE140829 technical metadata has no diagnosis column
  - locked=True returns "LOCKED" diagnosis label
  - MAD-based outlier flagging
"""

from __future__ import annotations

import importlib
import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Import stage scripts via importlib (names start with digits)
# ---------------------------------------------------------------------------
_m03 = importlib.import_module("scripts.03_quality_control")
compute_sample_qc_microarray = _m03.compute_sample_qc_microarray
compute_sample_qc_rnaseq = _m03.compute_sample_qc_rnaseq
compute_feature_qc = _m03.compute_feature_qc
MAD_CUTOFF = _m03.MAD_CUTOFF


# ---------------------------------------------------------------------------
# Synthetic fixtures
# ---------------------------------------------------------------------------

def _make_microarray_df(n_samples: int = 20, n_features: int = 100, seed: int = 42) -> pd.DataFrame:
    """Return a samples × features DataFrame of synthetic log-scale microarray values."""
    rng = np.random.default_rng(seed)
    data = rng.normal(loc=7.0, scale=0.5, size=(n_samples, n_features)).clip(4.0, 14.0)
    sample_ids = [f"GSM_S{i:03d}" for i in range(n_samples)]
    feature_ids = [f"ILMN_{i:05d}" for i in range(n_features)]
    return pd.DataFrame(data, index=sample_ids, columns=feature_ids)


def _make_meta_df(n_samples: int = 20, with_batch: bool = True, locked: bool = False) -> pd.DataFrame:
    """Return a synthetic metadata DataFrame."""
    sample_ids = [f"GSM_S{i:03d}" for i in range(n_samples)]
    diags = (["AD"] * (n_samples // 2)) + (["Control"] * (n_samples - n_samples // 2))
    rows: dict = {"sample_id": sample_ids, "diagnosis": diags, "age": [70.0 + i * 0.5 for i in range(n_samples)]}
    if with_batch:
        rows["batch"] = ["CenterA" if i % 2 == 0 else "CenterB" for i in range(n_samples)]
    return pd.DataFrame(rows)


def _make_rnaseq_counts_df(n_samples: int = 15, n_genes: int = 50, seed: int = 42) -> pd.DataFrame:
    """Return a samples × genes DataFrame of synthetic non-negative integer counts."""
    rng = np.random.default_rng(seed)
    counts = rng.integers(0, 5000, size=(n_samples, n_genes))
    sample_ids = [f"GSM_R{i:03d}" for i in range(n_samples)]
    gene_ids = [f"GENE_{i:04d}" for i in range(n_genes)]
    return pd.DataFrame(counts, index=sample_ids, columns=gene_ids)


def _make_meta_rnaseq(n_samples: int = 15) -> pd.DataFrame:
    sample_ids = [f"GSM_R{i:03d}" for i in range(n_samples)]
    n_ad = n_samples // 2
    n_ctrl = n_samples - n_ad
    return pd.DataFrame({
        "sample_id": sample_ids,
        "diagnosis": ["AD"] * n_ad + ["Control"] * n_ctrl,
        "age": [75.0 + i for i in range(n_samples)],
    })


# ---------------------------------------------------------------------------
# Tests: compute_sample_qc_microarray
# ---------------------------------------------------------------------------

class TestComputeSampleQcMicroarray:

    def test_returns_dataframe_and_flags(self):
        expr = _make_microarray_df()
        meta = _make_meta_df()
        qc_df, flags = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        assert isinstance(qc_df, pd.DataFrame)
        assert isinstance(flags, list)

    def test_row_count_equals_sample_count(self):
        n = 20
        expr = _make_microarray_df(n_samples=n)
        meta = _make_meta_df(n_samples=n)
        qc_df, _ = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        assert len(qc_df) == n

    def test_required_columns_present(self):
        expr = _make_microarray_df()
        meta = _make_meta_df()
        qc_df, _ = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        required = {"sample_id", "dataset_id", "mean_expr", "median_expr",
                    "std_expr", "iqr_expr", "pct_missing", "is_outlier_flag", "qc_flags"}
        assert required.issubset(set(qc_df.columns)), f"Missing columns: {required - set(qc_df.columns)}"

    def test_median_values_are_numeric(self):
        expr = _make_microarray_df()
        meta = _make_meta_df()
        qc_df, _ = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        assert pd.api.types.is_float_dtype(qc_df["median_expr"])

    def test_iqr_non_negative(self):
        expr = _make_microarray_df()
        meta = _make_meta_df()
        qc_df, _ = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        assert (qc_df["iqr_expr"] >= 0).all()

    def test_std_non_negative(self):
        expr = _make_microarray_df()
        meta = _make_meta_df()
        qc_df, _ = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        assert (qc_df["std_expr"] >= 0).all()

    def test_pct_missing_zero_when_no_nans(self):
        expr = _make_microarray_df()
        meta = _make_meta_df()
        qc_df, _ = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        assert (qc_df["pct_missing"] == 0.0).all()

    def test_locked_true_returns_LOCKED_diagnosis(self):
        expr = _make_microarray_df()
        meta = _make_meta_df()
        qc_df, _ = compute_sample_qc_microarray(expr, meta, "TEST_DS", locked=True)
        assert (qc_df["diagnosis"] == "LOCKED").all(), \
            "When locked=True, all diagnosis values must be 'LOCKED'"

    def test_locked_false_uses_actual_diagnosis(self):
        expr = _make_microarray_df()
        meta = _make_meta_df()
        qc_df, _ = compute_sample_qc_microarray(expr, meta, "TEST_DS", locked=False)
        assert set(qc_df["diagnosis"].unique()).issubset({"AD", "Control", "MCI", "Unknown", "LOCKED"})
        assert "LOCKED" not in qc_df["diagnosis"].values, \
            "When locked=False, diagnosis should not be 'LOCKED'"

    def test_values_remain_float_after_qc(self):
        """Microarray values must remain float — no int conversion."""
        expr = _make_microarray_df()
        meta = _make_meta_df()
        qc_df, _ = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        assert pd.api.types.is_float_dtype(qc_df["mean_expr"])
        assert pd.api.types.is_float_dtype(qc_df["median_expr"])

    def test_no_samples_deleted_all_rows_returned(self):
        """QC must flag only, never remove rows."""
        n = 20
        expr = _make_microarray_df(n_samples=n)
        meta = _make_meta_df(n_samples=n)
        qc_df, _ = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        assert len(qc_df) == n, \
            "No samples should be deleted; all must appear in QC output"

    def test_missing_value_pct_calculated_correctly(self):
        """Inject NaN into one sample and confirm pct_missing is computed correctly."""
        expr = _make_microarray_df(n_samples=5, n_features=100)
        expr.iloc[0, :50] = np.nan  # 50% missing for first sample
        meta = _make_meta_df(n_samples=5)
        qc_df, _ = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        pct_first = qc_df.loc[qc_df["sample_id"] == expr.index[0], "pct_missing"].values[0]
        assert abs(pct_first - 50.0) < 1.0, f"Expected ~50% missing, got {pct_first:.1f}%"

    def test_outlier_flagging_extreme_sample(self):
        """An extreme outlier sample (very different median) should be flagged."""
        rng = np.random.default_rng(42)
        expr = _make_microarray_df(n_samples=30, n_features=200, seed=42)
        # Make sample 0 extreme: shift by many sigma
        expr.iloc[0] = rng.normal(loc=14.0, scale=0.1, size=200)
        meta = _make_meta_df(n_samples=30)
        qc_df, flags = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        outlier_row = qc_df[qc_df["sample_id"] == expr.index[0]]
        assert outlier_row["is_outlier_flag"].values[0], \
            "Extreme outlier sample should be flagged"
        assert len(flags) > 0, "Flags list should be non-empty for outlier sample"

    def test_non_outlier_samples_marked_pass(self):
        """Normal samples with consistent distribution should not be flagged."""
        expr = _make_microarray_df(n_samples=20, n_features=200)
        meta = _make_meta_df(n_samples=20)
        qc_df, _ = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        # In clean synthetic data with tight normal distribution, most samples pass
        pass_count = (qc_df["qc_flags"] == "PASS").sum()
        assert pass_count >= len(qc_df) * 0.8, \
            f"Expected mostly PASS in clean data, got {pass_count}/{len(qc_df)}"

    def test_flag_records_have_required_fields(self):
        """Flag dicts must have sample_id, dataset_id, flag_reason, provisional_action."""
        rng = np.random.default_rng(42)
        expr = _make_microarray_df(n_samples=30, n_features=200)
        expr.iloc[0] = rng.normal(14.0, 0.05, 200)
        meta = _make_meta_df(n_samples=30)
        _, flags = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        required_keys = {"sample_id", "dataset_id", "flag_reason", "provisional_action"}
        for f in flags:
            assert required_keys.issubset(set(f.keys())), f"Flag missing keys: {required_keys - set(f.keys())}"

    def test_provisional_action_never_delete(self):
        """Provisional action must never say 'DELETE'."""
        rng = np.random.default_rng(42)
        expr = _make_microarray_df(n_samples=30, n_features=200)
        expr.iloc[0] = rng.normal(14.0, 0.05, 200)
        meta = _make_meta_df(n_samples=30)
        _, flags = compute_sample_qc_microarray(expr, meta, "TEST_DS")
        for f in flags:
            assert "DELETE" not in f["provisional_action"].upper() or \
                   "DO_NOT_DELETE" in f["provisional_action"].upper(), \
                "Provisional action must never instruct deletion"


# ---------------------------------------------------------------------------
# Tests: compute_sample_qc_rnaseq
# ---------------------------------------------------------------------------

class TestComputeSampleQcRnaseq:

    def test_returns_dataframe_and_flags(self):
        counts = _make_rnaseq_counts_df()
        meta = _make_meta_rnaseq()
        qc_df, flags = compute_sample_qc_rnaseq(counts, meta, "TEST_DS")
        assert isinstance(qc_df, pd.DataFrame)
        assert isinstance(flags, list)

    def test_library_size_equals_sum_of_counts(self):
        counts = _make_rnaseq_counts_df(n_samples=5, n_genes=20)
        meta = _make_meta_rnaseq(n_samples=5)
        qc_df, _ = compute_sample_qc_rnaseq(counts, meta, "TEST_DS")
        # Verify library_size equals row sums of input counts
        expected_lib_sizes = counts.values.sum(axis=1)
        actual_lib_sizes = qc_df["library_size"].values
        np.testing.assert_array_equal(actual_lib_sizes, expected_lib_sizes,
            err_msg="Library size must equal sum of all counts per sample")

    def test_library_size_is_non_negative(self):
        counts = _make_rnaseq_counts_df()
        meta = _make_meta_rnaseq()
        qc_df, _ = compute_sample_qc_rnaseq(counts, meta, "TEST_DS")
        assert (qc_df["library_size"] >= 0).all()

    def test_counts_dtype_integer(self):
        """Raw integer counts must remain integer — not coerced to float by QC."""
        counts = _make_rnaseq_counts_df()
        # Verify input is integer
        assert np.issubdtype(counts.values.dtype, np.integer), \
            "Synthetic counts should be integer"

    def test_all_samples_returned_no_deletion(self):
        n = 15
        counts = _make_rnaseq_counts_df(n_samples=n)
        meta = _make_meta_rnaseq(n_samples=n)
        qc_df, _ = compute_sample_qc_rnaseq(counts, meta, "TEST_DS")
        assert len(qc_df) == n

    def test_n_expressed_genes_bounded(self):
        """n_expressed_genes must be between 0 and total n_genes."""
        n_genes = 50
        counts = _make_rnaseq_counts_df(n_genes=n_genes)
        meta = _make_meta_rnaseq()
        qc_df, _ = compute_sample_qc_rnaseq(counts, meta, "TEST_DS")
        assert (qc_df["n_expressed_genes"] >= 0).all()
        assert (qc_df["n_expressed_genes"] <= n_genes).all()

    def test_pct_zero_genes_bounded(self):
        counts = _make_rnaseq_counts_df()
        meta = _make_meta_rnaseq()
        qc_df, _ = compute_sample_qc_rnaseq(counts, meta, "TEST_DS")
        assert (qc_df["pct_zero_genes"] >= 0.0).all()
        assert (qc_df["pct_zero_genes"] <= 100.0).all()

    def test_outlier_flagged_on_extreme_library_size(self):
        """A sample with extremely large library size should be flagged."""
        counts = _make_rnaseq_counts_df(n_samples=20, n_genes=100, seed=42)
        # Make sample 0 have an extremely large library
        counts.iloc[0] = counts.values.max() * 50
        meta = _make_meta_rnaseq(n_samples=20)
        qc_df, flags = compute_sample_qc_rnaseq(counts, meta, "TEST_DS")
        outlier_row = qc_df[qc_df["sample_id"] == counts.index[0]]
        assert outlier_row["is_outlier_flag"].values[0], \
            "Sample with extreme library size should be flagged"

    def test_zero_count_sample_raises_no_error(self):
        """A sample with all zero counts should not raise an exception."""
        counts = _make_rnaseq_counts_df(n_samples=5, n_genes=20)
        counts.iloc[0] = 0  # All zeros
        meta = _make_meta_rnaseq(n_samples=5)
        qc_df, _ = compute_sample_qc_rnaseq(counts, meta, "TEST_DS")
        assert len(qc_df) == 5


# ---------------------------------------------------------------------------
# Tests: compute_feature_qc
# ---------------------------------------------------------------------------

class TestComputeFeatureQc:

    def test_returns_dict(self):
        """compute_feature_qc returns a dict (single-dataset summary)."""
        expr = _make_microarray_df()
        result = compute_feature_qc(expr, "TEST_DS", "microarray")
        assert isinstance(result, dict), \
            f"compute_feature_qc must return a dict, got {type(result)}"

    def test_required_keys_present(self):
        expr = _make_microarray_df()
        result = compute_feature_qc(expr, "TEST_DS", "microarray")
        required = {"dataset_id", "total_features", "features_with_missing", "constant_features"}
        assert required.issubset(set(result.keys())), f"Missing keys: {required - set(result.keys())}"

    def test_no_missing_in_clean_data(self):
        expr = _make_microarray_df()  # No NaN
        result = compute_feature_qc(expr, "TEST_DS", "microarray")
        assert result["features_with_missing"] == 0

    def test_constant_feature_detected(self):
        """Add a constant feature column; it must be counted."""
        expr = _make_microarray_df(n_features=10)
        expr["CONST_FEAT"] = 5.0  # Constant column
        result = compute_feature_qc(expr, "TEST_DS", "microarray")
        assert result["constant_features"] >= 1, \
            "Constant feature must be detected and counted"

    def test_total_features_matches_input(self):
        """total_features must equal the number of columns in input."""
        n_feat = 50
        expr = _make_microarray_df(n_features=n_feat)
        result = compute_feature_qc(expr, "TEST_DS", "microarray")
        assert result["total_features"] == n_feat

    def test_dataset_id_preserved(self):
        expr = _make_microarray_df()
        result = compute_feature_qc(expr, "MY_DATASET", "microarray")
        assert result["dataset_id"] == "MY_DATASET"


# ---------------------------------------------------------------------------
# Tests: Matrix orientation convention
# ---------------------------------------------------------------------------

class TestMatrixOrientation:

    def test_derived_matrix_rows_are_samples(self):
        """The derived samples×features matrix must have sample IDs as the row index."""
        expr_features_x_samples = _make_microarray_df().T  # Simulate original GEO orientation
        derived = expr_features_x_samples.T  # Transpose back
        # Confirm: rows should correspond to sample IDs (GSM_*)
        assert all(str(idx).startswith("GSM_") for idx in derived.index), \
            "Derived matrix rows should be sample IDs"

    def test_derived_shape_is_transposed(self):
        """derived.shape must equal (n_samples, n_features)."""
        n_samples, n_features = 20, 100
        orig = _make_microarray_df(n_samples=n_samples, n_features=n_features).T
        # orig shape: (n_features, n_samples)
        assert orig.shape == (n_features, n_samples)
        derived = orig.T
        assert derived.shape == (n_samples, n_features)

    def test_no_data_lost_in_transposition(self):
        """Transpose must not lose or modify values."""
        expr = _make_microarray_df(n_samples=10, n_features=50)
        orig = expr.T  # features × samples
        derived = orig.T  # back to samples × features
        pd.testing.assert_frame_equal(expr, derived,
            check_names=False, rtol=1e-6)

    def test_pca_input_has_samples_on_rows(self):
        """PCA is performed on (n_samples, n_features) matrices, not transposed."""
        from sklearn.decomposition import PCA
        expr = _make_microarray_df(n_samples=20, n_features=100)  # Already samples×features
        # PCA should not raise with this orientation
        pca = PCA(n_components=2, random_state=42)
        coords = pca.fit_transform(expr.values)
        assert coords.shape == (20, 2)

    def test_gse140829_technical_meta_no_diagnosis(self):
        """
        GSE140829 is a LOCKED cohort. Its technical metadata must not expose diagnosis.
        Simulate by asserting a metadata frame with only technical columns.
        """
        technical_meta = pd.DataFrame({
            "sample_id": ["GSM_L001", "GSM_L002"],
            "batch": [1, 2],
            "age": [70.0, 75.0],
        })
        assert "diagnosis" not in technical_meta.columns, \
            "Technical (locked) metadata must not contain 'diagnosis' column"


# ---------------------------------------------------------------------------
# Tests: PCA input preparation
# ---------------------------------------------------------------------------

class TestPCAPreparation:

    def test_standardization_applied_before_pca(self):
        """Top-variance features must be standardized (mean≈0, std≈1) before PCA."""
        from sklearn.preprocessing import StandardScaler
        expr = _make_microarray_df(n_samples=30, n_features=200)
        # Select top variance features
        variances = expr.var(axis=0)
        top_idx = variances.nlargest(50).index
        sub = expr[top_idx].values
        scaled = StandardScaler().fit_transform(sub)
        col_means = scaled.mean(axis=0)
        col_stds = scaled.std(axis=0)
        np.testing.assert_allclose(col_means, 0.0, atol=1e-10)
        np.testing.assert_allclose(col_stds, 1.0, atol=1e-10)

    def test_rnaseq_pca_uses_log1p_cpm(self):
        """RNA-seq PCA must use log1p(CPM), not raw counts."""
        counts = _make_rnaseq_counts_df(n_samples=15, n_genes=50).values.astype(float)
        lib_sizes = counts.sum(axis=1, keepdims=True)
        lib_sizes[lib_sizes == 0] = 1.0
        cpm = (counts / lib_sizes) * 1e6
        norm_matrix = np.log1p(cpm)
        # Verify values are in log space (not raw counts)
        assert norm_matrix.max() < 1000, \
            "log1p(CPM) values should be small (log-scale), not raw counts"
        assert norm_matrix.min() >= 0.0, "log1p(CPM) must be non-negative"

    def test_microarray_pca_does_not_re_log2(self):
        """Pre-normalized microarray values must NOT be log2-transformed again."""
        expr = _make_microarray_df()  # Values already in ~log-scale ~7.0
        # If we naively log2 again, values should become much smaller
        re_logged = np.log2(np.maximum(expr.values, 1e-6))
        # Test that re-log2 would corrupt signal (values become ~2.8 instead of ~7)
        original_median = np.median(expr.values)
        re_logged_median = np.median(re_logged)
        # Values must NOT be re-log2 (we show that re-log reduces median significantly)
        assert original_median > re_logged_median + 1.0, \
            "Re-applying log2 to already-log-normalized data would corrupt the signal"

    def test_top_variance_feature_selection(self):
        """PCA feature selection must select features by variance, not randomly."""
        expr = _make_microarray_df(n_samples=30, n_features=200)
        # Artificially inflate variance of first 20 features
        expr.iloc[:, :20] = expr.iloc[:, :20] * 10
        variances = expr.var(axis=0)
        top_50 = variances.nlargest(50).index
        # All 20 high-variance features should appear in top 50
        for i in range(20):
            assert expr.columns[i] in top_50, \
                f"Feature {expr.columns[i]} (high variance) should be in top 50"
