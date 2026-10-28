#!/usr/bin/env bash
# System packages: compilers, CMake, Python and the build dependencies of Verilator.
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build pkg-config git ca-certificates curl \
    python3 python3-dev python3-venv \
    autoconf flex bison help2man perl libfl2 libfl-dev zlib1g zlib1g-dev ccache \
    less
rm -rf /var/lib/apt/lists/*
