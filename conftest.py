# conftest.py — project root conftest for pytest
# Ensures the project root is on sys.path so that
# `from scripts.utils.paths import ...` works in all tests.

import sys
from pathlib import Path

# Insert project root at the front of sys.path
PROJECT_ROOT = Path(__file__).parent.resolve()
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
