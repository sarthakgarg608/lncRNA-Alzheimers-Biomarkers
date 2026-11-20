"""
scripts/utils/validation.py
============================
Data-quality validation helpers for the lncSignal-AD pipeline.

None of the functions here modify the data — they only inspect and report.
Every check returns a ``ValidationResult`` with a boolean ``passed`` flag
and a human-readable ``message`` so that callers can decide how to handle
failures (log, raise, or accumulate).

Scientific rules enforced
--------------------------
- Diagnosis labels are restricted to: AD, Control, MCI, Unknown.
  Samples labelled MCI or Unknown must NEVER be silently merged with Control.
- Expression values must be numeric; non-numeric entries indicate a parser
  or format issue, not a valid missing value.
- Raw integer counts (RNA-seq) and normalized microarray intensities must NOT
  be mixed in the same matrix or confused.
- Matrix orientation conventions: features (probes/genes) on rows, samples on
  columns — this is the convention for raw/interim storage.  ML stages
  transpose to samples × features.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class ValidationResult:
    """Outcome of a single validation check."""
    check_name: str
    passed: bool
    message: str
    details: dict = field(default_factory=dict)

    def __str__(self) -> str:
        status = "PASS" if self.passed else "FAIL"
        return f"[{status}] {self.check_name}: {self.message}"


# ---------------------------------------------------------------------------
# Allowed diagnosis label vocabulary
# ---------------------------------------------------------------------------

VALID_DIAGNOSIS_LABELS: frozenset[str] = frozenset({"AD", "Control", "MCI", "Unknown"})

# Labels that must NOT be silently merged into controls
MUST_NOT_MERGE_INTO_CONTROL: frozenset[str] = frozenset({"MCI", "Unknown"})


# ---------------------------------------------------------------------------
# Sample / ID checks
# ---------------------------------------------------------------------------

def check_duplicate_ids(
    ids: Sequence[str],
    *,
    id_label: str = "sample_id",
) -> ValidationResult:
    """
    Detect duplicate values in a sequence of identifiers.

    Parameters
    ----------
    ids : sequence of str
    id_label : str
        Human-readable name for the type of ID (used in messages).

    Returns
    -------
    ValidationResult
    """
    seen: set[str] = set()
    duplicates: list[str] = []
    for i in ids:
        if i in seen:
            duplicates.append(i)
        seen.add(i)

    if duplicates:
        return ValidationResult(
            check_name="check_duplicate_ids",
            passed=False,
            message=(
                f"{len(duplicates)} duplicate {id_label}(s) detected.  "
                f"Duplicates may indicate technical replicates or a parser error."
            ),
            details={"duplicate_ids": duplicates[:20]},
        )
    return ValidationResult(
        check_name="check_duplicate_ids",
        passed=True,
        message=f"No duplicate {id_label}s found (n={len(ids)}).",
    )


def check_sample_id_alignment(
    metadata_ids: Sequence[str],
    expression_ids: Sequence[str],
    *,
    allow_subset: bool = False,
) -> ValidationResult:
    """
    Verify that sample IDs in the metadata and expression matrix match.

    Parameters
    ----------
    metadata_ids : sequence of str
        Sample IDs from the metadata table (e.g. GSM accessions).
    expression_ids : sequence of str
        Sample IDs from the expression matrix columns.
    allow_subset : bool
        If True, metadata can be a strict subset of expression columns.

    Returns
    -------
    ValidationResult
    """
    meta_set = set(metadata_ids)
    expr_set = set(expression_ids)

    only_meta = meta_set - expr_set
    only_expr = expr_set - meta_set

    if allow_subset:
        if meta_set <= expr_set:
            return ValidationResult(
                check_name="check_sample_id_alignment",
                passed=True,
                message=(
                    f"All {len(meta_set)} metadata IDs are present in expression "
                    f"matrix ({len(expr_set)} total columns)."
                ),
            )
    else:
        if not only_meta and not only_expr:
            return ValidationResult(
                check_name="check_sample_id_alignment",
                passed=True,
                message=f"Metadata and expression columns match exactly (n={len(meta_set)}).",
            )

    return ValidationResult(
        check_name="check_sample_id_alignment",
        passed=False,
        message=(
            f"{len(only_meta)} IDs only in metadata, "
            f"{len(only_expr)} only in expression matrix."
        ),
        details={
            "only_in_metadata": sorted(only_meta)[:10],
            "only_in_expression": sorted(only_expr)[:10],
        },
    )


# ---------------------------------------------------------------------------
# Metadata column checks
# ---------------------------------------------------------------------------

REQUIRED_METADATA_COLUMNS: list[str] = [
    "sample_id",
    "dataset_id",
    "diagnosis",
    "tissue",
    "platform_id",
]

OPTIONAL_METADATA_COLUMNS: list[str] = [
    "age",
    "sex",
    "batch",
    "brain_region",
    "original_label",
]


def check_missing_columns(
    df: pd.DataFrame,
    required: Sequence[str] = REQUIRED_METADATA_COLUMNS,
    *,
    df_label: str = "metadata",
) -> ValidationResult:
    """
    Verify that all required columns are present in ``df``.

    Parameters
    ----------
    df : pd.DataFrame
    required : sequence of str
    df_label : str

    Returns
    -------
    ValidationResult
    """
    missing = [c for c in required if c not in df.columns]
    if missing:
        return ValidationResult(
            check_name="check_missing_columns",
            passed=False,
            message=f"{len(missing)} required column(s) absent from {df_label}: {missing}",
            details={"missing_columns": missing},
        )
    return ValidationResult(
        check_name="check_missing_columns",
        passed=True,
        message=f"All {len(required)} required columns present in {df_label}.",
    )


def check_valid_diagnosis_labels(
    labels: Sequence[str | None],
    *,
    allowed: frozenset[str] = VALID_DIAGNOSIS_LABELS,
) -> ValidationResult:
    """
    Confirm that every diagnosis label is one of the allowed values.

    Raises a FAIL if any label outside the allowed vocabulary is found.
    ``None`` / NaN values are reported but do not fail (they should be
    mapped to ``"Unknown"`` by the parsing layer).

    IMPORTANT: MCI and Unknown must never be silently reclassified as
    Control.  This check does NOT enforce that — it only checks vocabulary.

    Parameters
    ----------
    labels : sequence of str or None
    allowed : frozenset[str]

    Returns
    -------
    ValidationResult
    """
    illegal: list[str] = []
    n_null = 0
    for lab in labels:
        if lab is None or (isinstance(lab, float) and np.isnan(lab)):
            n_null += 1
            continue
        if lab not in allowed:
            illegal.append(lab)

    issues: list[str] = []
    if n_null:
        issues.append(f"{n_null} null/NaN label(s) (should be mapped to 'Unknown')")
    if illegal:
        issues.append(f"{len(illegal)} label(s) outside vocabulary {allowed}: {list(set(illegal))[:5]}")

    if illegal:
        return ValidationResult(
            check_name="check_valid_diagnosis_labels",
            passed=False,
            message="; ".join(issues),
            details={"illegal_labels": list(set(illegal))},
        )
    return ValidationResult(
        check_name="check_valid_diagnosis_labels",
        passed=True,
        message=(
            "All diagnosis labels are valid"
            + (f" ({n_null} mapped to Unknown)" if n_null else "")
            + "."
        ),
    )


# ---------------------------------------------------------------------------
# Expression matrix checks
# ---------------------------------------------------------------------------

def check_expression_numeric(
    expr: pd.DataFrame,
    *,
    max_nan_fraction: float = 0.3,
) -> ValidationResult:
    """
    Verify that expression values are numeric and that the NaN fraction is
    below a threshold.

    Parameters
    ----------
    expr : pd.DataFrame
        Expression matrix (features × samples or samples × features).
    max_nan_fraction : float
        Maximum acceptable fraction of NaN values.

    Returns
    -------
    ValidationResult
    """
    n_total = expr.size
    if n_total == 0:
        return ValidationResult(
            check_name="check_expression_numeric",
            passed=False,
            message="Expression matrix is empty.",
        )

    # Check for non-numeric dtype columns
    non_numeric_cols = [c for c in expr.columns if not pd.api.types.is_numeric_dtype(expr[c])]
    if non_numeric_cols:
        return ValidationResult(
            check_name="check_expression_numeric",
            passed=False,
            message=f"{len(non_numeric_cols)} non-numeric column(s): {non_numeric_cols[:5]}",
            details={"non_numeric_columns": non_numeric_cols[:20]},
        )

    n_nan = int(expr.isna().sum().sum())
    nan_frac = n_nan / n_total

    if nan_frac > max_nan_fraction:
        return ValidationResult(
            check_name="check_expression_numeric",
            passed=False,
            message=(
                f"NaN fraction {nan_frac:.2%} exceeds threshold {max_nan_fraction:.2%}."
            ),
            details={"n_nan": n_nan, "nan_fraction": nan_frac},
        )
    return ValidationResult(
        check_name="check_expression_numeric",
        passed=True,
        message=f"Expression is numeric. NaN fraction: {nan_frac:.4%} ({n_nan}/{n_total}).",
        details={"n_nan": n_nan, "nan_fraction": nan_frac},
    )


def check_expression_dimensions(
    expr: pd.DataFrame,
    *,
    min_features: int = 100,
    min_samples: int = 10,
    orientation: str = "features_x_samples",
) -> ValidationResult:
    """
    Check that the expression matrix has plausible dimensions.

    Parameters
    ----------
    expr : pd.DataFrame
    min_features : int
        Minimum expected number of rows (features).
    min_samples : int
        Minimum expected number of columns (samples).
    orientation : str
        Expected orientation: ``"features_x_samples"`` or ``"samples_x_features"``.

    Returns
    -------
    ValidationResult
    """
    n_rows, n_cols = expr.shape
    issues: list[str] = []

    if orientation == "features_x_samples":
        if n_rows < min_features:
            issues.append(f"Only {n_rows} rows (features); expected ≥{min_features}.")
        if n_cols < min_samples:
            issues.append(f"Only {n_cols} columns (samples); expected ≥{min_samples}.")
    elif orientation == "samples_x_features":
        if n_cols < min_features:
            issues.append(f"Only {n_cols} columns (features); expected ≥{min_features}.")
        if n_rows < min_samples:
            issues.append(f"Only {n_rows} rows (samples); expected ≥{min_samples}.")

    if issues:
        return ValidationResult(
            check_name="check_expression_dimensions",
            passed=False,
            message="; ".join(issues),
            details={"shape": (n_rows, n_cols)},
        )
    return ValidationResult(
        check_name="check_expression_dimensions",
        passed=True,
        message=f"Dimensions OK: {n_rows} rows × {n_cols} columns ({orientation}).",
        details={"shape": (n_rows, n_cols)},
    )


def check_matrix_orientation(
    expr: pd.DataFrame,
    *,
    expected: str = "features_x_samples",
) -> ValidationResult:
    """
    Heuristic check that the matrix is oriented as expected.

    Uses the index name and column name patterns as heuristics:
    - ``"ID_REF"`` or ``"feature_id"`` index → features on rows.
    - Index values starting with ``"GSM"`` → samples on rows (transposed).

    Parameters
    ----------
    expr : pd.DataFrame
    expected : str
        ``"features_x_samples"`` or ``"samples_x_features"``.

    Returns
    -------
    ValidationResult
    """
    index_name = str(expr.index.name or "").lower()
    sample_prefix_in_index = any(
        str(v).startswith("GSM") for v in expr.index[:5]
    )
    sample_prefix_in_cols = any(
        str(c).startswith("GSM") for c in expr.columns[:5]
    )

    # Features-on-rows: index looks like probe IDs, columns look like GSM
    if expected == "features_x_samples":
        if sample_prefix_in_cols and not sample_prefix_in_index:
            return ValidationResult(
                check_name="check_matrix_orientation",
                passed=True,
                message="Matrix orientation appears correct: features × samples.",
            )
        elif sample_prefix_in_index:
            return ValidationResult(
                check_name="check_matrix_orientation",
                passed=False,
                message=(
                    "Matrix index contains GSM-like values — matrix may be transposed. "
                    "Expected features × samples orientation."
                ),
            )

    elif expected == "samples_x_features":
        if sample_prefix_in_index and not sample_prefix_in_cols:
            return ValidationResult(
                check_name="check_matrix_orientation",
                passed=True,
                message="Matrix orientation appears correct: samples × features.",
            )
        elif sample_prefix_in_cols:
            return ValidationResult(
                check_name="check_matrix_orientation",
                passed=False,
                message=(
                    "Matrix columns contain GSM-like values — matrix may need transposition. "
                    "Expected samples × features orientation."
                ),
            )

    return ValidationResult(
        check_name="check_matrix_orientation",
        passed=True,
        message=f"Orientation check inconclusive (expected: {expected}). "
                "Verify manually using index/column names.",
    )


def check_expression_scale_flag(
    dataset_cfg: dict,
) -> ValidationResult:
    """
    Verify that the ``expression_scale`` field in the dataset configuration
    is set to a known value and is consistent with the downstream method.

    Valid values and appropriate tools:
    - ``"normalized_microarray"``  →  use limma / statsmodels (not DESeq2)
    - ``"raw_integer_counts"``     →  use PyDESeq2 / DESeq2

    Parameters
    ----------
    dataset_cfg : dict
        One dataset block from ``config/datasets.yaml``.

    Returns
    -------
    ValidationResult
    """
    known_scales = {
        "normalized_microarray": "limma or statsmodels (not DESeq2 / PyDESeq2)",
        "raw_integer_counts": "PyDESeq2 or DESeq2 (not limma)",
    }
    scale = dataset_cfg.get("expression_scale", None)

    if scale is None:
        return ValidationResult(
            check_name="check_expression_scale_flag",
            passed=False,
            message="'expression_scale' field is missing from dataset configuration.",
        )
    if scale not in known_scales:
        return ValidationResult(
            check_name="check_expression_scale_flag",
            passed=False,
            message=f"Unknown expression_scale '{scale}'. "
                    f"Valid values: {list(known_scales.keys())}",
        )
    return ValidationResult(
        check_name="check_expression_scale_flag",
        passed=True,
        message=f"expression_scale='{scale}'. Appropriate DE tool: {known_scales[scale]}.",
    )
