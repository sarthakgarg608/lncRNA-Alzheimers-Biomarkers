"""
scripts/utils/paths.py
======================
Path resolution utilities for the lncSignal-AD project.

All paths are built relative to the project root so that the code is portable
across machines and operating systems.  No absolute paths are hard-coded.

Usage example
-------------
>>> from scripts.utils.paths import get_project_root, load_dataset_config
>>> root = get_project_root()
>>> cfg  = load_dataset_config()
>>> raw  = resolve_path(cfg["datasets"]["GSE63061"]["series_matrix_path"])
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Lazy import for yaml — avoids hard dependency at module load time for
# callers that only need get_project_root().
# ---------------------------------------------------------------------------
try:
    import yaml as _yaml
    _YAML_AVAILABLE = True
except ImportError:
    _yaml = None          # type: ignore[assignment]
    _YAML_AVAILABLE = False


# ---------------------------------------------------------------------------
# Project root detection
# ---------------------------------------------------------------------------

def get_project_root() -> Path:
    """
    Return the absolute path of the project root directory.

    Strategy (in order of priority):
    1. ``LNCSIGNAL_ROOT`` environment variable (useful for CI / Colab).
    2. Walk up from this file's location until a directory containing
       ``config/datasets.yaml`` **or** ``AGENTS.md`` is found.
    3. Fall back to the current working directory (works when the script is
       run from the project root, which is the documented usage).

    Returns
    -------
    Path
        Absolute project root path.

    Raises
    ------
    FileNotFoundError
        If none of the heuristics locate a recognisable project root.
    """
    # 1. Environment variable override
    env_root = os.environ.get("LNCSIGNAL_ROOT")
    if env_root:
        p = Path(env_root).resolve()
        if p.is_dir():
            return p

    # 2. Walk up from this file looking for a marker
    markers = {"AGENTS.md", "config"}
    current = Path(__file__).resolve()
    for parent in [current, *current.parents]:
        if any((parent / m).exists() for m in markers):
            return parent

    # 3. Fall back to cwd
    cwd = Path.cwd().resolve()
    if (cwd / "config").exists() or (cwd / "AGENTS.md").exists():
        return cwd

    raise FileNotFoundError(
        "Could not detect the project root.  "
        "Set the LNCSIGNAL_ROOT environment variable to the project directory, "
        "or run scripts from the project root."
    )


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def get_config_path() -> Path:
    """Return the path to ``config/datasets.yaml``."""
    return get_project_root() / "config" / "datasets.yaml"


def load_dataset_config() -> dict[str, Any]:
    """
    Load and return the contents of ``config/datasets.yaml`` as a dictionary.

    Returns
    -------
    dict
        Parsed YAML document.

    Raises
    ------
    ImportError
        If PyYAML is not installed.
    FileNotFoundError
        If ``config/datasets.yaml`` does not exist.
    """
    if not _YAML_AVAILABLE:
        raise ImportError(
            "PyYAML is required to load the dataset configuration.  "
            "Install it with: pip install pyyaml"
        )
    cfg_path = get_config_path()
    if not cfg_path.is_file():
        raise FileNotFoundError(
            f"Dataset configuration not found at: {cfg_path}\n"
            "Create config/datasets.yaml from the project template."
        )
    with cfg_path.open("r", encoding="utf-8") as fh:
        return _yaml.safe_load(fh)


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

def resolve_path(relative_path: str | Path) -> Path:
    """
    Resolve a path that is relative to the project root.

    Parameters
    ----------
    relative_path : str or Path
        Path as stored in ``config/datasets.yaml`` (e.g.
        ``"data/raw/GSE63061_series_matrix.txt.gz"``).

    Returns
    -------
    Path
        Absolute path.
    """
    p = Path(relative_path)
    if p.is_absolute():
        return p
    return get_project_root() / p


def resolve_dataset_paths(dataset_cfg: dict[str, Any]) -> dict[str, Path]:
    """
    Given a single dataset's config block, return a dict of *existing* file
    paths keyed by their YAML field name (e.g. ``"series_matrix_path"``).
    Missing optional files are skipped; missing required files raise an error.

    Parameters
    ----------
    dataset_cfg : dict
        One dataset's block from ``config/datasets.yaml``.

    Returns
    -------
    dict[str, Path]
        Only those path fields that resolve to existing files.
    """
    path_keys = [k for k in dataset_cfg if k.endswith("_path")]
    resolved: dict[str, Path] = {}
    for key in path_keys:
        raw_val = dataset_cfg[key]
        if raw_val is None:
            continue
        p = resolve_path(raw_val)
        if p.exists():
            resolved[key] = p
        # If not found, skip silently — callers can check completeness.
    return resolved
