"""!
@file conftest.py
@brief The RTL tests need the SystemC model (kessc) and Verilator; skip them otherwise.
"""

import shutil
import sys
from pathlib import Path

MODEL = Path(__file__).resolve().parents[1] / "model"
sys.path.insert(0, str(MODEL))
from build_model import build_model  # noqa: E402

collect_ignore_glob = []

_build = None if shutil.which("verilator") is None else build_model()
if _build is None or not any(_build.glob("kessc*")):
    collect_ignore_glob.append("test_*.py")
else:
    sys.path.insert(0, str(_build))
