"""!
@file build_model.py
@brief Build the kesmodel and kessc extensions with CMake into model/build.

Used by the pytest configuration of the model tests and of the RTL testbench.
The build directory is reused, so only changed sources are compiled again.
"""

import shutil
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Optional

MODEL = Path(__file__).resolve().parent
BUILD = MODEL / "build"


def build_model() -> Optional[Path]:
    """!
    @brief Configure and build the extensions.
    @return The build directory, or None when cmake or pybind11 are missing.
    @throws RuntimeError If configuring or building fails.
    """
    try:
        import pybind11
    except ImportError:
        warnings.warn("pybind11 is not installed: the C++ models are not built")
        return None
    if shutil.which("cmake") is None:
        warnings.warn("cmake was not found: the C++ models are not built")
        return None
    configure = ["cmake", "-S", str(MODEL), "-B", str(BUILD), "-DCMAKE_BUILD_TYPE=Release",
                 f"-DPython_EXECUTABLE={sys.executable}", f"-Dpybind11_DIR={pybind11.get_cmake_dir()}"]
    if shutil.which("ninja") and not (BUILD / "CMakeCache.txt").exists():
        configure += ["-G", "Ninja"]
    for cmd in (configure, ["cmake", "--build", str(BUILD)]):
        done = subprocess.run(cmd, capture_output=True, text=True)
        if done.returncode:
            raise RuntimeError(f"{' '.join(cmd)} failed:\n{done.stdout}\n{done.stderr}")
    return BUILD
