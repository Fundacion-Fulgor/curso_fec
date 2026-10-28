"""!
@file test_kes_timing.py
@brief Unit tests for kes_timing.py: the port protocol of key_equation_solver.
"""

import pytest

from kes_ref import ARCHITECTURES, make_kes
from kes_timing import KesTimed, simulate
from rs_code import RSCode
from stimuli import StimulusGenerator

NAMES = sorted(ARCHITECTURES)
IDLE = (False, False, None)


@pytest.fixture(scope="module")
def qr():
    return RSCode.qr_2m()


@pytest.fixture(scope="module")
def words(qr):
    gen = StimulusGenerator(qr, seed=41)
    return [gen.random_errors(e).syndromes for e in (0, 3, 8, 9, 1)]


def expected(name, qr, syn):
    res = make_kes(name, qr.field, qr.t).run(syn)
    return res.lam, res.omega


@pytest.mark.parametrize("name", NAMES)
class TestTiming:
    """!
    @brief Protocol checks, run for every architecture.
    """

    def test_latency(self, name, qr, words):
        n = make_kes(name, qr.field, qr.t).iterations
        out = simulate(name, qr.field, qr.t, [(False, True, words[2])] + [IDLE] * (n + 3))
        assert [i for i, p in enumerate(out) if p.o_valid] == [n]
        assert (out[n].o_lambda, out[n].o_omega) == expected(name, qr, words[2])
        dut = KesTimed(name, qr.field, qr.t)
        assert dut.latency == dut.interval == n

    def test_back_to_back(self, name, qr, words):
        n = make_kes(name, qr.field, qr.t).iterations
        stimulus = []
        for syn in words:
            stimulus += [(False, True, syn)] + [IDLE] * (n - 1)
        stimulus += [IDLE] * 2
        out = simulate(name, qr.field, qr.t, stimulus)
        valid = [i for i, p in enumerate(out) if p.o_valid]
        assert valid == [n * (j + 1) for j in range(len(words))]
        for i, syn in zip(valid, words):
            assert (out[i].o_lambda, out[i].o_omega) == expected(name, qr, syn)

    def test_new_word_restarts(self, name, qr, words):
        n = make_kes(name, qr.field, qr.t).iterations
        stimulus = [(False, True, words[1])] + [IDLE] * 4 + [(False, True, words[2])] + [IDLE] * n
        out = simulate(name, qr.field, qr.t, stimulus)
        assert [i for i, p in enumerate(out) if p.o_valid] == [5 + n]
        assert (out[5 + n].o_lambda, out[5 + n].o_omega) == expected(name, qr, words[2])

    def test_outputs_hold(self, name, qr, words):
        n = make_kes(name, qr.field, qr.t).iterations
        out = simulate(name, qr.field, qr.t, [(False, True, words[1])] + [IDLE] * (2 * n))
        result = (out[n].o_lambda, out[n].o_omega)
        for p in out[n + 1:]:
            assert not p.o_valid and (p.o_lambda, p.o_omega) == result and not p.busy

    def test_reset(self, name, qr, words):
        n = make_kes(name, qr.field, qr.t).iterations
        stimulus = [(False, True, words[1])] + [IDLE] * n + [(False, True, words[2])] + [IDLE] * 3
        stimulus += [(True, False, None)] + [IDLE] * (n + 2)
        out = simulate(name, qr.field, qr.t, stimulus)
        assert [i for i, p in enumerate(out) if p.o_valid] == [n]
        after = out[n + 5]
        assert not after.busy and after.state is None
        assert after.o_lambda == [0] * (qr.t + 1) and after.o_omega == [0] * qr.t

    def test_reset_wins_over_valid(self, name, qr, words):
        out = simulate(name, qr.field, qr.t, [(True, True, words[1])] + [IDLE] * 30)
        assert not any(p.o_valid or p.busy for p in out)

    def test_internal_state_follows_the_trace(self, name, qr, words):
        trace = make_kes(name, qr.field, qr.t).run(words[3], trace=True).trace
        dut = KesTimed(name, qr.field, qr.t)
        dut.clock(valid=True, syndromes=words[3])
        for r in range(dut.latency):
            assert dut.ports().state == r
            assert dut.kes.snapshot() == trace[r]
            dut.clock()

    def test_buses(self, name, qr, words):
        dut = KesTimed(name, qr.field, qr.t)
        dut.clock(valid=True, syndromes=words[1])
        while not dut.o_valid:
            dut.clock()
        lam, omega = expected(name, qr, words[1])
        assert dut.lambda_bus() == qr.pack(lam, ascending=True)
        assert dut.omega_bus() == qr.pack(omega, ascending=True)


class TestArguments:
    """!
    @brief Invalid inputs.
    """

    def test_valid_needs_syndromes(self, qr):
        with pytest.raises(ValueError):
            KesTimed("ribm_systolic", qr.field, qr.t).clock(valid=True)

    def test_bad_syndromes(self, qr):
        with pytest.raises(ValueError):
            KesTimed("ribm", qr.field, qr.t).clock(valid=True, syndromes=[0] * 15)

    def test_unknown_architecture(self, qr):
        with pytest.raises(ValueError):
            KesTimed("euclid", qr.field, qr.t)

    def test_idle(self, qr):
        out = simulate("ibm", qr.field, qr.t, [IDLE] * 40)
        assert not any(p.o_valid or p.busy for p in out)
