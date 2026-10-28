#!/usr/bin/env bash
# Build SystemC from the Accellera repository as a shared library (needed by pybind11).
# Usage: systemc.sh <tag> <prefix> [jobs]
set -euo pipefail

tag="${1:?SystemC tag, e.g. 3.0.2}"
prefix="${2:?install prefix}"
jobs="${3:-$(nproc)}"

src="$(mktemp -d)"
git clone --quiet --depth 1 --branch "$tag" https://github.com/accellera-official/systemc.git "$src"
cmake -S "$src" -B "$src/build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$prefix" \
    -DCMAKE_CXX_STANDARD=17 \
    -DBUILD_SHARED_LIBS=ON \
    -DDISABLE_COPYRIGHT_MESSAGE=ON \
    -DENABLE_EXAMPLES=OFF
cmake --build "$src/build" -j "$jobs"
cmake --install "$src/build"
rm -rf "$src"
