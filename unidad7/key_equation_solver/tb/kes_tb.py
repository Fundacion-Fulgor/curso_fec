"""!
@file kes_tb.py
@brief pyuvm testbench of key_equation_solver (riBM).

The scoreboard runs the SystemC model (kessc) in lockstep with the RTL: on every
clock edge both see the same inputs, and the outputs and the internal registers
must match. Every result is also decoded with the Chien search and the Forney
formula and compared with the errors injected by the sequence.

KES_CODE selects the configuration ("qr": M = 8, T = 8; "tp3": M = 4, T = 3) and
KES_SEED the seed of the random sequences.
"""

import os
import random
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import cocotb
import pyuvm
from cocotb.clock import Clock
from cocotb.triggers import ReadOnly, RisingEdge
from pyuvm import (uvm_analysis_port, uvm_component, uvm_driver, uvm_env, uvm_sequence, uvm_sequence_item, uvm_sequencer, uvm_subscriber, uvm_test)

import kessc
from decode_ref import chien_forney
from rs_code import RSCode
from stimuli import Case, StimulusGenerator

CODE_NAME = os.environ.get("KES_CODE", "qr")
CODE = RSCode.qr_2m() if CODE_NAME == "qr" else RSCode.tp3_reduced()
SEED = int(os.environ.get("KES_SEED", "1"))
T = CODE.t
ITERATIONS = 2 * T

# syndromes -> stimuli.Case, filled by the sequences so the scoreboard knows the injected errors
KNOWN_CASES: Dict[tuple, Case] = {}


class KesItem(uvm_sequence_item):
    """!
    @brief gap idle clock edges followed by one action: "load", "reset" or "idle".
    """

    def __init__(self, name="kes_item", kind="load", syndromes=None, gap=0):
        super().__init__(name)
        self.kind = kind
        self.syndromes = syndromes
        self.gap = gap


class KesSeqBase(uvm_sequence):
    """!
    @brief Helpers shared by the sequences. Every sequence starts with a reset
           and ends with enough idle edges to drain the last word.
    """

    def __init__(self, name, seed=SEED):
        super().__init__(name)
        self.gen = StimulusGenerator(CODE, seed=seed)
        self.rng = random.Random(seed)

    async def send(self, kind="load", syndromes=None, gap=0):
        item = KesItem("item", kind, syndromes, gap)
        await self.start_item(item)
        await self.finish_item(item)

    async def load_case(self, case, gap):
        syndromes = case.syndromes[: 2 * T]
        KNOWN_CASES[tuple(syndromes)] = case
        await self.send("load", syndromes, gap)

    async def body(self):
        await self.send("reset")
        await self.words()
        await self.send("idle", gap=ITERATIONS + 2)

    async def words(self):
        raise NotImplementedError


class KesReferenceSeq(KesSeqBase):
    """!
    @brief Known words: the TP3 exercise 6 vector with two errors (qr), or the
           words r1 and r2 of exercise 3 (tp3), followed by an all-zero syndrome.
    """

    async def words(self):
        if CODE_NAME == "qr":
            cases = [self.gen.case({40: 0x55, 3: 0xA7}, list(range(28)), "reference")]
        else:
            msg = CODE.unpack(0x52D09F17C, 9)
            cw = CODE.encode(msg)
            cases = []
            for word in (0x52D0AF17CB1EE64, 0x56D09F171B1E06F):
                rx = CODE.unpack(word, 15)
                errors = {CODE.n - 1 - i: a ^ b for i, (a, b) in enumerate(zip(rx, cw)) if a != b}
                cases.append(self.gen.case(errors, msg, "exercise_3"))
        cases.append(self.gen.clean("zeros"))
        for case in cases:
            await self.load_case(case, gap=ITERATIONS + 3)


class KesSuiteSeq(KesSeqBase):
    """!
    @brief The deterministic suite of stimuli.py, back to back (one word every 2T clocks).
    """

    async def words(self):
        for case in self.gen.suite():
            await self.load_case(case, gap=ITERATIONS - 1)


class KesRandomSeq(KesSeqBase):
    """!
    @brief Random words at random distances: back to back, with gaps, too early
           (the word in progress is dropped) and with resets in the middle.
    """

    def __init__(self, name, words=300, seed=SEED):
        super().__init__(name, seed)
        self.n_words = words

    async def words(self):
        for _ in range(self.n_words):
            kind = self.rng.random()
            if kind < 0.1:
                await self.send("reset", gap=self.rng.randint(0, ITERATIONS))
            pick = self.rng.random()
            if pick < 0.55:
                case = self.gen.random_errors(self.rng.randint(0, T))
            elif pick < 0.85:
                case = self.gen.random_errors(self.rng.randint(T + 1, min(2 * T, CODE.n)))
            else:
                case = None
            gap = self.rng.choice([ITERATIONS - 1, ITERATIONS - 1,
                                   ITERATIONS + self.rng.randint(0, 5),
                                   self.rng.randint(0, ITERATIONS - 2)])
            if case is None:
                await self.send("load", self.gen.random_syndromes(), gap)
            else:
                await self.load_case(case, gap)


class KesDriver(uvm_driver):
    """!
    @brief Drives one item per call. Inputs are written right after a rising edge
           and sampled by the next one; i_valid and i_rst return to 0 after it.
    """

    def build_phase(self):
        self.dut = cocotb.top

    async def run_phase(self):
        while True:
            item = await self.seq_item_port.get_next_item()
            for _ in range(item.gap):
                await self.edge()
            if item.kind == "load":
                await self.edge(valid=True, syndromes=item.syndromes)
            elif item.kind == "reset":
                await self.edge(rst=True)
            else:
                await self.edge()
            self.seq_item_port.item_done()

    async def edge(self, rst=False, valid=False, syndromes=None):
        dut = self.dut
        dut.i_rst.value = int(rst)
        dut.i_valid.value = int(valid)
        if syndromes is not None:
            dut.i_syn.value = CODE.pack(syndromes, ascending=True)
        await RisingEdge(dut.i_clk)
        dut.i_rst.value = 0
        dut.i_valid.value = 0


@dataclass
class KesCycle:
    """!
    @brief What the monitor sees on one edge: the sampled inputs, the outputs
           after the edge and the internal registers.
    """

    rst: bool
    valid: bool
    syndromes: Optional[List[int]]
    o_valid: bool
    o_lambda: int
    o_omega: int
    busy: bool
    state: int
    regs: Dict[str, object] = field(default_factory=dict)


class KesMonitor(uvm_component):
    """!
    @brief Samples the inputs on each rising edge and, once the edge has settled,
           the outputs and the registers lam_q, b_q, delta_q, theta_q, gamma and k.
    """

    def build_phase(self):
        self.ap = uvm_analysis_port("ap", self)
        self.dut = cocotb.top

    async def run_phase(self):
        dut = self.dut
        ctrl = dut.u_ctrl
        while True:
            await RisingEdge(dut.i_clk)
            rst = bool(dut.i_rst.value)
            valid = bool(dut.i_valid.value)
            syndromes = CODE.unpack(int(dut.i_syn.value), 2 * T, ascending=True) if valid else None
            await ReadOnly()
            regs = {
                "lam": CODE.unpack(int(dut.lam_q.value), T + 1, ascending=True),
                "b": CODE.unpack(int(dut.b_q.value), T + 1, ascending=True),
                "delta": CODE.unpack(int(dut.delta_q.value), 2 * T, ascending=True),
                "theta": CODE.unpack(int(dut.theta_q.value), 2 * T, ascending=True),
                "gamma": int(ctrl.gamma_q.value),
                "k": ctrl.k_q.value.to_signed(),
                "mc": bool(dut.mc.value),
            }
            self.ap.write(KesCycle(rst, valid, syndromes, bool(dut.o_valid.value),
                                   int(dut.o_lambda.value), int(dut.o_omega.value),
                                   bool(ctrl.busy_q.value), int(ctrl.count_q.value), regs))


class KesScoreboard(uvm_subscriber):
    """!
    @brief Lockstep comparison with kessc and decoding of every result.
    """

    def build_phase(self):
        self.ref = kessc.harness()
        self.errors = 0
        self.cycles = 0
        self.results = []
        self.decoded = 0
        self.flagged = 0
        self.beyond_t = 0
        self.current = None
        self.prev_busy = False
        self.prev_regs = None

    def write(self, c):
        self.cycles += 1
        p = self.ref.clock(CODE_NAME, "ribm", rst=c.rst, valid=c.valid, syndromes=c.syndromes)
        diffs = []
        if c.o_valid != p["o_valid"]:
            diffs.append(f"o_valid {c.o_valid} != {p['o_valid']}")
        if c.o_lambda != p["lambda_bus"]:
            diffs.append(f"o_lambda {c.o_lambda:x} != {p['lambda_bus']:x}")
        if c.o_omega != p["omega_bus"]:
            diffs.append(f"o_omega {c.o_omega:x} != {p['omega_bus']:x}")
        if c.busy != p["busy"]:
            diffs.append(f"busy {c.busy} != {p['busy']}")
        if p["busy"]:
            s = self.ref.snapshot(CODE_NAME, "ribm")
            if c.state != s["r"]:
                diffs.append(f"state {c.state} != {s['r']}")
            for name, value in c.regs.items():
                if value != s[name]:
                    diffs.append(f"{name} {value} != {s[name]}")
        # design intent, not protocol: the datapath holds while idle so it does not toggle
        if not self.prev_busy and not c.valid and self.prev_regs is not None:
            changed = [n for n in ("lam", "b", "delta", "theta", "gamma", "k")
                       if c.regs[n] != self.prev_regs[n]]
            if changed:
                diffs.append(f"registers changed while idle: {changed}")
        self.prev_busy = c.busy
        self.prev_regs = c.regs
        if diffs:
            self.errors += 1
            if self.errors <= 10:
                self.logger.error(f"cycle {self.cycles}: " + "; ".join(diffs))

        if c.o_valid:
            self.check_result(c)
        if c.rst:
            self.current = None
        elif c.valid:
            self.current = c.syndromes

    def check_result(self, c):
        """!
        @brief Decode the result of the word that was in progress.
        """
        self.results.append((c.o_lambda, c.o_omega))
        case = KNOWN_CASES.get(tuple(self.current)) if self.current else None
        if case is None:
            return
        lam = CODE.unpack(c.o_lambda, T + 1, ascending=True)
        omega = CODE.unpack(c.o_omega, T, ascending=True)
        cf = chien_forney(CODE.field, CODE.n, lam, omega, "high")
        if case.correctable:
            if cf.fail or cf.errors != CODE.error_word(case.errors):
                self.errors += 1
                self.logger.error(f"wrong decoding of {case.kind}: {case.errors} -> "
                                  f"{dict(zip(cf.positions, cf.values))}, fail {cf.fail}")
            else:
                self.decoded += 1
        else:
            self.beyond_t += 1
            self.flagged += cf.fail

    def report_phase(self):
        self.logger.info(f"{self.cycles} cycles, {len(self.results)} results, {self.decoded} decoded, "
                         f"{self.flagged}/{self.beyond_t} beyond t flagged, {self.errors} errors")


class KesCoverage(uvm_subscriber):
    """!
    @brief Functional coverage: degree of each result, protocol situations and MC decisions.
    """

    def build_phase(self):
        self.hits = Counter()
        self.prev_busy = False

    def write(self, c):
        if c.o_valid:
            lam = CODE.unpack(c.o_lambda, T + 1, ascending=True)
            degree = max((i for i, v in enumerate(lam) if v), default=-1)
            self.hits[f"degree_{degree}"] += 1
            if c.valid:
                self.hits["back_to_back"] += 1
        if c.valid and not c.rst and self.prev_busy and not c.o_valid:
            self.hits["restart"] += 1
        if c.rst and self.prev_busy:
            self.hits["reset_while_busy"] += 1
        if c.busy and c.regs["mc"]:
            self.hits["mc"] += 1
        self.prev_busy = c.busy and not c.rst

    def report_phase(self):
        self.logger.info("coverage: " + ", ".join(f"{k} {v}" for k, v in sorted(self.hits.items())))


class KesEnv(uvm_env):
    def build_phase(self):
        self.seqr = uvm_sequencer("seqr", self)
        self.driver = KesDriver("driver", self)
        self.monitor = KesMonitor("monitor", self)
        self.scoreboard = KesScoreboard("scoreboard", self)
        self.coverage = KesCoverage("coverage", self)

    def connect_phase(self):
        self.driver.seq_item_port.connect(self.seqr.seq_item_export)
        self.monitor.ap.connect(self.scoreboard.analysis_export)
        self.monitor.ap.connect(self.coverage.analysis_export)


class KesTestBase(uvm_test):
    """!
    @brief Runs one sequence and checks the scoreboard and the required coverage.
    """

    required: set = set()

    def build_phase(self):
        self.env = KesEnv("env", self)

    def make_sequence(self):
        raise NotImplementedError

    async def run_phase(self):
        self.raise_objection()
        Clock(cocotb.top.i_clk, 10, unit="ns").start()
        await self.make_sequence().start(self.env.seqr)
        await RisingEdge(cocotb.top.i_clk)
        self.drop_objection()

    def check_phase(self):
        sb, hits = self.env.scoreboard, self.env.coverage.hits
        assert sb.errors == 0, f"{sb.errors} scoreboard errors"
        assert sb.results, "no result was produced"
        missing = {name for name in self.required if not hits[name]}
        assert not missing, f"coverage holes: {sorted(missing)}"


@pyuvm.test()
class KesReferenceTest(KesTestBase):
    """!
    @brief Known vectors, with the expected buses of the QR reference word.
    """

    def make_sequence(self):
        return KesReferenceSeq("reference")

    def check_phase(self):
        super().check_phase()
        sb = self.env.scoreboard
        if CODE_NAME == "qr":
            assert sb.results[0] == (0x6674BB, 0xB837), [hex(v) for v in sb.results[0]]
        assert sb.decoded == len(sb.results), "every reference word must decode"


@pyuvm.test()
class KesSuiteTest(KesTestBase):
    """!
    @brief Every category of stimuli.py, back to back.
    """

    required = {f"degree_{d}" for d in range(T + 1)} | {"back_to_back", "mc"}

    def make_sequence(self):
        return KesSuiteSeq("suite")


@pyuvm.test()
class KesRandomTest(KesTestBase):
    """!
    @brief Random protocol stream with words beyond t and random syndromes.
    """

    required = {f"degree_{d}" for d in range(T + 1)} | {"back_to_back", "restart", "reset_while_busy", "mc"}

    def make_sequence(self):
        return KesRandomSeq("random")
