#!/usr/bin/env bash
# Build Verilator from source (Ubuntu packages are too old for cocotb 2.x).
# Usage: verilator.sh <tag> <prefix> [jobs]
set -euo pipefail

tag="${1:?Verilator tag, e.g. v5.052}"
prefix="${2:?install prefix}"
jobs="${3:-$(nproc)}"

src="$(mktemp -d)"
git clone --quiet --depth 1 --branch "$tag" https://github.com/verilator/verilator.git "$src"
cd "$src"
autoconf
./configure --prefix="$prefix"
make -j "$jobs"
make install
cd /
rm -rf "$src"
