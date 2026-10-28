"""!
@file conftest.py
@brief Build the C++ models before collecting their tests; skip them without cmake or pybind11.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_model import build_model  # noqa: E402

collect_ignore_glob = []

_build = build_model()
if _build is None:
    collect_ignore_glob.append("test_*.py")
else:
    sys.path.insert(0, str(_build))
