"""!
@file smoke_tb.py
@brief Minimal pyuvm test on Verilator: drive the registered XOR and check it.
"""

import cocotb
import pyuvm
from cocotb.clock import Clock
from cocotb.triggers import FallingEdge, ReadOnly, RisingEdge
from pyuvm import uvm_test

PAIRS = [(0x01, 0x02), (0xFF, 0x0F), (0x55, 0xAA), (0x00, 0x00)]


@pyuvm.test()
class SmokeTest(uvm_test):
    """!
    @brief Reset, then one input pair per clock with the result checked after the edge.
    """

    async def run_phase(self):
        self.raise_objection()
        dut = cocotb.top
        Clock(dut.i_clk, 10, unit="ns").start()
        dut.i_rst.value = 1
        dut.i_a.value = 0
        dut.i_b.value = 0
        await RisingEdge(dut.i_clk)
        for a, b in PAIRS:
            await FallingEdge(dut.i_clk)
            dut.i_rst.value = 0
            dut.i_a.value = a
            dut.i_b.value = b
            await RisingEdge(dut.i_clk)
            await ReadOnly()
            assert int(dut.o_q.value) == a ^ b
        self.logger.info("checked %d pairs", len(PAIRS))
        self.drop_objection()
