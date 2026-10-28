"""!
@file test_systemc.py
@brief The SystemC model (kessc) against kes_timing.py, cycle by cycle.

Random protocol streams mix back-to-back words, gaps, words dropped by an
early load and resets. On every edge the ports, the buses and the internal
state of both models must match.
"""

import random

import pytest

from kes_timing import KesTimed
from rs_code import RSCode
from stimuli import StimulusGenerator

kessc = pytest.importorskip("kessc")

CODES = {"qr": RSCode.qr_2m(), "tp3": RSCode.tp3_reduced()}
INSTANCES = [(c, a) for c in CODES for a in ("ibm", "ribm", "ribm_systolic")]
IDLE = (False, False, None)


@pytest.fixture(scope="module")
def harness():
    return kessc.harness()


def random_stream(code, iterations, seed, words=40):
    """!
    @brief Edges with loads at random distances: back-to-back, with gaps,
           too early (the previous word is dropped), and occasional resets.
    """
    rng = random.Random(seed)
    gen = StimulusGenerator(code, seed=seed)
    stream = [(True, False, None)]
    for _ in range(words):
        if rng.random() < 0.5:
            syn = gen.random_errors(rng.randint(0, min(2 * code.t, code.n))).syndromes
        else:
            syn = gen.random_syndromes()
        stream.append((False, True, syn[: 2 * code.t]))
        kind = rng.random()
        if kind < 0.4:
            gap = iterations - 1
        elif kind < 0.7:
            gap = iterations + rng.randint(0, 5)
        elif kind < 0.9:
            gap = rng.randint(0, iterations - 2)
        else:
            gap = rng.randint(0, iterations)
            stream += [IDLE] * gap + [(True, False, None)]
            continue
        stream += [IDLE] * gap
    return stream + [IDLE] * (iterations + 2)


def run_both(harness, code_name, arch, stream):
    """!
    @brief Drive both models and check every edge. Returns the number of results.
    """
    code = CODES[code_name]
    ref = KesTimed(arch, code.field, code.t)
    results = 0
    for edge, (rst, valid, syn) in enumerate(stream):
        p = harness.clock(code_name, arch, rst=rst, valid=valid, syndromes=syn)
        r = ref.clock(rst, valid, syn)
        where = f"{code_name}/{arch} edge {edge}"
        assert p["o_valid"] == r.o_valid, where
        assert (p["o_lambda"], p["o_omega"]) == (r.o_lambda, r.o_omega), where
        assert (p["busy"], p["state"]) == (r.busy, r.state), where
        assert p["lambda_bus"] == code.pack(r.o_lambda, ascending=True), where
        assert p["omega_bus"] == code.pack(r.o_omega, ascending=True), where
        assert harness.snapshot(code_name, arch) == ref.kes.snapshot(), where
        results += r.o_valid
    return results


@pytest.mark.parametrize("code_name, arch", INSTANCES)
@pytest.mark.parametrize("seed", [1, 2, 3])
def test_random_streams(harness, code_name, arch, seed):
    iterations = harness.config(code_name, arch)["iterations"]
    stream = random_stream(CODES[code_name], iterations, seed)
    assert run_both(harness, code_name, arch, stream) >= 10


@pytest.mark.parametrize("code_name, arch", INSTANCES)
def test_back_to_back(harness, code_name, arch):
    code = CODES[code_name]
    n = harness.config(code_name, arch)["iterations"]
    gen = StimulusGenerator(code, seed=9)
    stream = [(True, False, None)]
    for e in range(code.t + 1):
        stream += [(False, True, gen.random_errors(e).syndromes[: 2 * code.t])] + [IDLE] * (n - 1)
    stream += [IDLE] * 2
    assert run_both(harness, code_name, arch, stream) == code.t + 1


def test_other_instances_stay_idle(harness):
    gen = StimulusGenerator(CODES["qr"], seed=10)
    harness.clock("qr", "ibm", rst=True)
    before = harness.ports("qr", "ibm")
    harness.clock("qr", "ribm", valid=True, syndromes=gen.random_syndromes())
    for _ in range(20):
        harness.clock("qr", "ribm")
    assert harness.ports("qr", "ibm") == before
    assert not before["busy"]


def test_time_and_cycles(harness):
    t0, c0 = harness.time_ns, harness.cycles
    for _ in range(7):
        harness.clock("tp3", "ribm")
    assert harness.time_ns == t0 + 70
    assert harness.cycles == c0 + 7


def test_configuration(harness):
    assert kessc.harness() is harness
    assert sorted(harness.instances) == sorted(INSTANCES)
    qr = harness.config("qr", "ribm_systolic")
    assert (qr["m"], qr["poly"], qr["t"], qr["iterations"], qr["widths"]) == (8, 0x1D, 8, 16, (128, 72, 64))
    tp3 = harness.config("tp3", "ibm")
    assert (tp3["m"], tp3["t"], tp3["iterations"], tp3["widths"]) == (4, 3, 9, (24, 16, 12))
    assert kessc.release().startswith("3.0")


def test_errors_do_not_advance_time(harness):
    t0 = harness.time_ns
    with pytest.raises(ValueError, match="no instance"):
        harness.clock("qr", "euclid")
    with pytest.raises(ValueError, match="no syndromes"):
        harness.clock("qr", "ribm", valid=True)
    with pytest.raises(ValueError):
        harness.clock("qr", "ribm", valid=True, syndromes=[0] * 15)
    with pytest.raises(ValueError):
        harness.clock("tp3", "ribm", valid=True, syndromes=[0] * 5 + [16])
    with pytest.raises(ValueError):
        harness.ports("tp3", "euclid")
    assert harness.time_ns == t0
