"""
tests/test_paths.py
===================
Tests for scripts/utils/paths.py

Uses only file-system operations and synthetic YAML fixtures.
No real GEO datasets are read in these tests.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
import yaml


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _project_root() -> Path:
    """Return the project root (parent of tests/)."""
    return Path(__file__).parent.parent


# ---------------------------------------------------------------------------
# Tests: get_project_root
# ---------------------------------------------------------------------------

class TestGetProjectRoot:
    def test_returns_path_object(self):
        from scripts.utils.paths import get_project_root
        root = get_project_root()
        assert isinstance(root, Path)

    def test_root_contains_agents_md(self):
        from scripts.utils.paths import get_project_root
        root = get_project_root()
        assert (root / "AGENTS.md").exists(), (
            f"AGENTS.md not found at {root}.  "
            "Project root detection may be wrong."
        )

    def test_root_contains_config_dir(self):
        from scripts.utils.paths import get_project_root
        root = get_project_root()
        assert (root / "config").is_dir(), "config/ directory not found in project root."

    def test_env_var_override(self, tmp_path, monkeypatch):
        """LNCSIGNAL_ROOT env var should override auto-detection."""
        from scripts.utils.paths import get_project_root
        monkeypatch.setenv("LNCSIGNAL_ROOT", str(tmp_path))
        # tmp_path exists, so it should be returned directly
        result = get_project_root()
        assert result == tmp_path.resolve()


# ---------------------------------------------------------------------------
# Tests: get_config_path
# ---------------------------------------------------------------------------

class TestGetConfigPath:
    def test_returns_path_to_yaml(self):
        from scripts.utils.paths import get_config_path
        p = get_config_path()
        assert isinstance(p, Path)
        assert p.suffix == ".yaml"
        assert p.name == "datasets.yaml"

    def test_yaml_file_exists(self):
        from scripts.utils.paths import get_config_path
        p = get_config_path()
        assert p.exists(), f"datasets.yaml not found at {p}"


# ---------------------------------------------------------------------------
# Tests: load_dataset_config
# ---------------------------------------------------------------------------

class TestLoadDatasetConfig:
    def test_loads_yaml_successfully(self):
        from scripts.utils.paths import load_dataset_config
        cfg = load_dataset_config()
        assert isinstance(cfg, dict)

    def test_has_datasets_key(self):
        from scripts.utils.paths import load_dataset_config
        cfg = load_dataset_config()
        assert "datasets" in cfg, "config/datasets.yaml must have a 'datasets' top-level key."

    def test_has_all_four_datasets(self):
        from scripts.utils.paths import load_dataset_config
        cfg = load_dataset_config()
        datasets = cfg["datasets"]
        for expected_id in ["GSE63061", "GSE63060", "GSE140829", "GSE125583"]:
            assert expected_id in datasets, f"{expected_id} missing from datasets.yaml"

    def test_each_dataset_has_required_keys(self):
        from scripts.utils.paths import load_dataset_config
        required = {"geo_accession", "role", "tissue", "technology",
                    "platform_id", "series_matrix_path", "expression_scale",
                    "feature_id_type", "locked"}
        cfg = load_dataset_config()
        for ds_id, ds in cfg["datasets"].items():
            missing = required - set(ds.keys())
            assert not missing, (
                f"Dataset '{ds_id}' in datasets.yaml missing keys: {missing}"
            )

    def test_gse140829_is_locked(self):
        """External validation labels must be locked."""
        from scripts.utils.paths import load_dataset_config
        cfg = load_dataset_config()
        ds = cfg["datasets"]["GSE140829"]
        assert ds.get("locked") is True, (
            "GSE140829 (external validation) must have locked=true in datasets.yaml."
        )

    def test_gse125583_has_counts_path(self):
        """Brain RNA-seq dataset must specify a counts_path."""
        from scripts.utils.paths import load_dataset_config
        cfg = load_dataset_config()
        ds = cfg["datasets"]["GSE125583"]
        assert "counts_path" in ds, (
            "GSE125583 config must include 'counts_path' pointing to the RNA-seq counts file."
        )

    def test_expression_scale_values(self):
        """All datasets must declare a known expression_scale."""
        from scripts.utils.paths import load_dataset_config
        valid_scales = {"normalized_microarray", "raw_integer_counts"}
        cfg = load_dataset_config()
        for ds_id, ds in cfg["datasets"].items():
            scale = ds.get("expression_scale")
            assert scale in valid_scales, (
                f"Dataset '{ds_id}' has invalid expression_scale='{scale}'. "
                f"Must be one of: {valid_scales}"
            )

    def test_missing_config_raises(self, tmp_path, monkeypatch):
        """load_dataset_config should raise FileNotFoundError if YAML is absent."""
        from scripts.utils import paths as paths_mod
        monkeypatch.setattr(paths_mod, "_YAML_AVAILABLE", True)
        monkeypatch.setenv("LNCSIGNAL_ROOT", str(tmp_path))
        with pytest.raises(FileNotFoundError):
            paths_mod.load_dataset_config()


# ---------------------------------------------------------------------------
# Tests: resolve_path
# ---------------------------------------------------------------------------

class TestResolvePath:
    def test_relative_path_resolved_to_absolute(self):
        from scripts.utils.paths import resolve_path
        p = resolve_path("data/raw")
        assert p.is_absolute()

    def test_absolute_path_returned_unchanged(self):
        from scripts.utils.paths import resolve_path
        abs_path = Path(__file__).resolve()
        result = resolve_path(abs_path)
        assert result == abs_path

    def test_series_matrix_paths_exist(self):
        """All series_matrix_path values in config must point to existing files."""
        from scripts.utils.paths import load_dataset_config, resolve_path
        cfg = load_dataset_config()
        for ds_id, ds in cfg["datasets"].items():
            sm_path = resolve_path(ds["series_matrix_path"])
            assert sm_path.exists(), (
                f"series_matrix_path for {ds_id} does not exist: {sm_path}"
            )

    def test_annotation_paths_exist(self):
        """All annotation_path values must point to existing files."""
        from scripts.utils.paths import load_dataset_config, resolve_path
        cfg = load_dataset_config()
        for ds_id, ds in cfg["datasets"].items():
            ann = ds.get("annotation_path")
            if ann:
                p = resolve_path(ann)
                assert p.exists(), (
                    f"annotation_path for {ds_id} does not exist: {p}"
                )
