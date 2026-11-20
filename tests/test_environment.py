"""
Environment and dependency unit tests for lncSignal-AD.
"""

import importlib
import importlib.metadata
import platform
import subprocess
import sys
import pytest
from packaging.version import parse as parse_version

ESSENTIAL_PACKAGES = [
    ("pandas", "pandas", "2.0"),
    ("numpy", "numpy", "1.24"),
    ("scipy", "scipy", "1.10"),
    ("statsmodels", "statsmodels", "0.14"),
    ("GEOparse", "GEOparse", "2.0"),
    ("pydeseq2", "pydeseq2", "0.4"),
    ("biopython", "Bio", "1.81"),
    ("gseapy", "gseapy", "1.1"),
    ("scikit-learn", "sklearn", "1.3"),
    ("joblib", "joblib", "1.3"),
    ("matplotlib", "matplotlib", "3.7"),
    ("seaborn", "seaborn", "0.12"),
    ("networkx", "networkx", "3.0"),
    ("pyarrow", "pyarrow", "14.0"),
    ("pytest", "pytest", "8.0"),
]


def test_python_version_compatibility():
    """Verify running on supported Python 3.11 or 3.12."""
    assert sys.version_info.major == 3
    assert sys.version_info.minor in (11, 12), f"Unsupported Python minor version: {sys.version_info.minor}"


@pytest.mark.parametrize("pip_name,import_name,min_ver", ESSENTIAL_PACKAGES)
def test_essential_packages_importable_and_version(pip_name, import_name, min_ver):
    """Verify each essential package is importable and meets minimum version requirement."""
    mod = importlib.import_module(import_name)
    assert mod is not None, f"Failed to import {import_name}"

    try:
        ver = importlib.metadata.version(pip_name)
    except Exception:
        ver = getattr(mod, "__version__", None)

    assert ver is not None, f"Could not determine version for {pip_name}"
    assert parse_version(ver) >= parse_version(min_ver), (
        f"{pip_name} version {ver} is below required {min_ver}"
    )


def test_pip_check_clean():
    """Verify pip check reports no broken requirements."""
    res = subprocess.run(
        [sys.executable, "-m", "pip", "check"],
        capture_output=True,
        text=True,
        check=False
    )
    assert res.returncode == 0, f"pip check failed: {res.stdout} {res.stderr}"
