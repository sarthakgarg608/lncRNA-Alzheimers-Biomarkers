"""
scripts/utils/__init__.py
=========================
Makes scripts/utils a Python package so utilities can be imported as:

    from scripts.utils.paths import get_project_root, resolve_dataset_paths
    from scripts.utils.data_io import load_series_matrix_header, load_counts_tsv
    from scripts.utils.validation import check_duplicate_ids, validate_metadata

This file intentionally re-exports the most-used public symbols.
"""

from scripts.utils.paths import (
    get_project_root,
    get_config_path,
    load_dataset_config,
    resolve_path,
)

from scripts.utils.data_io import (
    load_series_matrix_header,
    load_series_matrix_expression,
    load_tsv_gz,
    load_plain_tsv,
    load_gtf_gz_header,
)

from scripts.utils.validation import (
    check_duplicate_ids,
    check_missing_columns,
    check_valid_diagnosis_labels,
    check_expression_numeric,
    check_expression_dimensions,
    check_sample_id_alignment,
    check_matrix_orientation,
    check_expression_scale_flag,
    ValidationResult,
)

__all__ = [
    # paths
    "get_project_root",
    "get_config_path",
    "load_dataset_config",
    "resolve_path",
    # data_io
    "load_series_matrix_header",
    "load_series_matrix_expression",
    "load_tsv_gz",
    "load_plain_tsv",
    "load_gtf_gz_header",
    # validation
    "check_duplicate_ids",
    "check_missing_columns",
    "check_valid_diagnosis_labels",
    "check_expression_numeric",
    "check_expression_dimensions",
    "check_sample_id_alignment",
    "check_matrix_orientation",
    "check_expression_scale_flag",
    "ValidationResult",
]
