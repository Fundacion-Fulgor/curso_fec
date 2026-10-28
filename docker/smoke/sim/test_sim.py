"""!
@file test_sim.py
@brief Build the smoke DUT with Verilator and run the pyuvm test through the cocotb runner.

Also checks that Verilator writes a SAIF activity file, used later for power estimation.
"""

import os
import xml.etree.ElementTree as ET
from pathlib import Path

from cocotb_tools.runner import get_runner

HERE = Path(__file__).resolve().parent


def test_verilator_cocotb_pyuvm(tmp_path):
    build = Path(os.environ.get("FEC_SMOKE_BUILD", tmp_path))
    runner = get_runner("verilator")
    runner.build(sources=[HERE / "smoke_dut.sv"], hdl_toplevel="smoke_dut",
                 build_dir=build, build_args=["--trace-saif"], waves=True)
    results = runner.test(hdl_toplevel="smoke_dut", test_module="smoke_tb",
                          build_dir=build, test_dir=build, waves=True,
                          extra_env={"PYTHONPATH": str(HERE)})

    cases = list(ET.parse(results).getroot().iter("testcase"))
    bad = [c.get("name") for c in cases
           if c.find("failure") is not None or c.find("error") is not None
           or c.find("skipped") is not None]
    assert len(cases) == 1 and not bad, f"test failed or did not run: {bad}"
    assert (build / "dump.saif").stat().st_size > 0
