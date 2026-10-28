"""!
@file test_kes_rtl.py
@brief Build key_equation_solver with Verilator and run the pyuvm tests of kes_tb.py.

One build per configuration: "qr" (M = 8, T = 8, the TP3 module) and "tp3"
(M = 4, T = 3, the reduced code). Set KES_WAVES=1 to dump dump.vcd and
KES_SAIF=1 to dump dump.saif into sim_build/<code>.
"""

import os
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from cocotb_tools.runner import get_runner

HERE = Path(__file__).resolve().parent
KES = HERE.parent
UNIDAD7 = KES.parent
RTL = [KES / "rtl" / f for f in ("kes_gf_mul.v", "kes_pe0.v", "kes_pe1.v", "kes_ctrl.v",
                                 "key_equation_solver.v")]
CODES = {"qr": {"M": 8, "POLY": 0x1D, "T": 8}, "tp3": {"M": 4, "POLY": 0x3, "T": 3}}
TESTS = ["KesReferenceTest", "KesSuiteTest", "KesRandomTest"]


@pytest.mark.parametrize("code", sorted(CODES))
def test_key_equation_solver(code):
    saif = os.environ.get("KES_SAIF") == "1"
    waves = saif or os.environ.get("KES_WAVES") == "1"
    build = HERE / "sim_build" / code
    pythonpath = [HERE, UNIDAD7, UNIDAD7.parent / "algebra", KES / "model" / "build"]

    runner = get_runner("verilator")
    runner.build(sources=RTL, hdl_toplevel="key_equation_solver", parameters=CODES[code],
                 build_dir=build, waves=waves, build_args=["--trace-saif"] if saif else [])
    results = runner.test(hdl_toplevel="key_equation_solver", test_module="kes_tb",
                          build_dir=build, test_dir=build, waves=waves,
                          extra_env={"KES_CODE": code,
                                     "PYTHONPATH": os.pathsep.join(str(p) for p in pythonpath)})

    cases = list(ET.parse(results).getroot().iter("testcase"))
    failed = [c.get("name") for c in cases if c.find("failure") is not None or c.find("error") is not None]
    assert sorted(c.get("name") for c in cases) == sorted(TESTS)
    assert not failed, f"failed: {failed}"
