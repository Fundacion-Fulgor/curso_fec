#!/usr/bin/env python3
"""!
@file run_demo.py
@brief Build and run the demo of Example 31 of Zhang on the RTL.

Prints the riBM registers of every iteration next to Table 4.3 of the book and
leaves the waveform in demo/sim_build/dump.vcd, to open with the view
zhang_demo.gtkw. From unidad7:

    ../docker/run.sh python3 key_equation_solver/demo/run_demo.py
"""

import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from cocotb_tools.runner import get_runner

HERE = Path(__file__).resolve().parent
KES = HERE.parent
UNIDAD7 = KES.parent
SOURCES = [KES / "rtl" / f for f in ("kes_gf_mul.v", "kes_pe0.v", "kes_pe1.v", "kes_ctrl.v", "key_equation_solver.v")] + [HERE / "zhang_demo.v"]


def run(build_dir: Path = HERE / "sim_build") -> bool:
    """!
    @brief Build the demo with Verilator and run it with waveforms.
    @param build_dir Where the simulator and dump.vcd are written.
    @return True if the RTL matched the book.
    """
    pythonpath = os.pathsep.join(str(p) for p in (HERE, UNIDAD7, UNIDAD7.parent / "algebra"))
    runner = get_runner("verilator")
    runner.build(sources=SOURCES, hdl_toplevel="zhang_demo", build_dir=build_dir, waves=True)
    results = runner.test(hdl_toplevel="zhang_demo", test_module="zhang_demo", build_dir=build_dir, test_dir=build_dir, waves=True, extra_env={"PYTHONPATH": pythonpath})
    cases = list(ET.parse(results).getroot().iter("testcase"))
    return len(cases) == 1 and cases[0].find("failure") is None and cases[0].find("error") is None


def main() -> int:
    ok = run()
    print()
    print("RTL and Table 4.3 match" if ok else "the RTL does NOT match Table 4.3")
    print(f"waveform: {HERE / 'sim_build' / 'dump.vcd'}")
    print(f"view it from {HERE} with:  gtkwave sim_build/dump.vcd zhang_demo.gtkw")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
