#!/usr/bin/env bash
# Check the container: tool versions, SystemC called from Python, Verilator with cocotb and
# pyuvm (including SAIF output), and the repository tests when they are present.
# Usage, from the host: docker/run.sh docker/smoke/smoke_test.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
build="${TMPDIR:-/tmp}/fec-smoke"
rm -rf "$build"
mkdir -p "$build"
pytest_opts=(-q -p no:cacheprovider)

echo "== versions"
verilator --version
cmake --version | head -n 1
g++ --version | head -n 1
python3 - << 'PY'
from importlib.metadata import version
import sys
print("python", sys.version.split()[0], "|",
      " | ".join(f"{p} {version(p)}" for p in ("cocotb", "pyuvm", "cocotb-coverage", "pybind11", "pytest")))
PY

echo "== SystemC + pybind11"
cmake -S "$here/sc_bind" -B "$build/sc_bind" -G Ninja -DCMAKE_BUILD_TYPE=Release > "$build/sc_bind.log"
cmake --build "$build/sc_bind" >> "$build/sc_bind.log"
PYTHONPATH="$build/sc_bind" python3 -m pytest "${pytest_opts[@]}" "$here/sc_bind"

echo "== Verilator + cocotb + pyuvm"
FEC_SMOKE_BUILD="$build/sim" python3 -m pytest "${pytest_opts[@]}" "$here/sim"

if [ -d "$repo/algebra" ]; then
    echo "== algebra"
    (cd "$repo/algebra" && python3 -m pytest "${pytest_opts[@]}" test_galois_unit.py)
fi
if [ -d "$repo/unidad7" ]; then
    echo "== unidad7"
    (cd "$repo/unidad7" && python3 -m pytest "${pytest_opts[@]}")
fi

echo "== all good"
