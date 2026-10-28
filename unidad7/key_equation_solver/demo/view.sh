#!/usr/bin/env bash
# Open the demo waveform with its view. GTKWave reads the alpha filter relative to the
# current directory, so this script runs it from the demo directory.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
if [ ! -f sim_build/dump.vcd ]; then
    echo "sim_build/dump.vcd does not exist: run run_demo.py first" >&2
    exit 1
fi
exec gtkwave sim_build/dump.vcd zhang_demo.gtkw
