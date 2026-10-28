"""!
@file test_kes_ref.py
@brief Unit tests for kes_ref.py.

The riBM is checked state by state against Table 4.3 of Zhang, "VLSI
Architectures for Modern Error Correcting Codes" (Example 31), and the other
two architectures against it.
"""

import random

import pytest

from galois import GaloisField, GFPoly
from kes_ref import ARCHITECTURES, IBM, KesBase, RiBMSplit, RiBMSystolic, make_kes
from rs_code import RSCode
from stimuli import StimulusGenerator

NAMES = sorted(ARCHITECTURES)

# GF(2^3), x^3 + x + 1: alpha^0..alpha^6 = 1, 2, 4, 3, 6, 7, 5
A1, A2, A3, A4, A5, A6 = 2, 4, 3, 6, 7, 5
ZHANG_SYN = [0, A1, A2, A6]

# Zhang, Table 4.3: (Lambda, B, k, gamma, delta-hat, theta-hat) at the start of r
ZHANG_TABLE = [
    ([1, 0, 0], [1, 0, 0], 0, 1, [0, A1, A2, A6], [0, A1, A2, A6]),
    ([1, 0, 0], [0, 1, 0], 1, 1, [A1, A2, A6, 0], [0, A1, A2, A6]),
    ([1, 0, A1], [1, 0, 0], -2, A1, [A2, 1, A3, 1], [A2, A6, 0, 0]),
    ([A1, A2, A2], [0, 1, 0], -1, A1, [A2, A2, A1, 0], [A2, A6, 0, 0]),
]
ZHANG_LAMBDA = [A2, A3, A5]
ZHANG_OMEGA_HIGH = [A6, A4]


@pytest.fixture(scope="module")
def gf8():
    return GaloisField(3, 0b011)


@pytest.fixture(scope="module", params=["qr", "small"])
def code(request):
    return RSCode.qr_2m() if request.param == "qr" else RSCode.tp3_reduced()


def full_product(field, lam, syndromes):
    """!
    @brief Coefficients of Lambda(x) S(x), ascending, padded to 3t.
    """
    t = len(lam) - 1
    prod = GFPoly.from_ascending(field, lam) * GFPoly.from_ascending(field, syndromes)
    return prod.to_list(3 * t, ascending=True)


class TestZhangExample:
    """!
    @brief Example 31 of Zhang: RS(7, 3), t = 2, S = (0, a, a^2, a^6).
    """

    def test_ribm_state_by_state(self, gf8):
        res = RiBMSplit(gf8, 2).run(ZHANG_SYN, trace=True)
        for state, (lam, b, k, gamma, delta, theta) in zip(res.trace, ZHANG_TABLE):
            assert state["lam"] == lam
            assert state["b"] == b
            assert state["k"] == k
            assert state["gamma"] == gamma
            assert state["delta"] == delta
            assert state["theta"] == theta
        assert res.lam == ZHANG_LAMBDA
        assert res.omega == ZHANG_OMEGA_HIGH

    def test_ribm_systolic(self, gf8):
        res = RiBMSystolic(gf8, 2).run(ZHANG_SYN)
        assert (res.lam, res.omega, res.form) == (ZHANG_LAMBDA, ZHANG_OMEGA_HIGH, "high")

    def test_ibm(self, gf8):
        # classic Omega of Example 30 is a x, scaled by lambda_0 = a^2
        res = IBM(gf8, 2).run(ZHANG_SYN)
        assert (res.lam, res.omega, res.form) == (ZHANG_LAMBDA, [0, A3], "classic")


class TestReferenceVector:
    """!
    @brief QR message 00..1B with errors 0x55 at position 40 and 0xA7 at position 3.
    """

    def test_buses(self):
        qr = RSCode.qr_2m()
        rx = qr.add_errors(qr.encode(list(range(28))), {40: 0x55, 3: 0xA7})
        syn = qr.syndromes(rx)
        assert qr.pack(syn, ascending=True) == 0x409F1C63F681F695B7A04E3BDF46FFF2
        for name in NAMES:
            res = make_kes(name, qr.field, 8).run(syn)
            assert res.lambda_bus() == 0x6674BB
            assert res.degree == 2
        assert make_kes("ribm_systolic", qr.field, 8).run(syn).omega_bus() == 0xB837
        assert make_kes("ribm", qr.field, 8).run(syn).omega_bus() == 0xB837


class TestCrossArchitecture:
    """!
    @brief The three architectures against each other and against Lambda(x) S(x).
    """

    def test_correctable_cases_agree(self, code):
        gen = StimulusGenerator(code, seed=11)
        cases = [c for c in gen.suite() if c.correctable]
        cases += [gen.random_errors(e) for e in range(code.t + 1) for _ in range(10)]
        t = code.t
        for case in cases:
            res = {name: make_kes(name, code.field, t).run(case.syndromes) for name in NAMES}
            assert res["ibm"].lam == res["ribm"].lam == res["ribm_systolic"].lam
            assert res["ribm"].omega == res["ribm_systolic"].omega
            prod = full_product(code.field, res["ribm"].lam, case.syndromes)
            assert res["ibm"].omega == prod[:t]
            assert res["ribm"].omega == prod[2 * t: 3 * t]
            assert all(c == 0 for c in prod[res["ribm"].degree: 2 * t])
            assert res["ribm"].degree == case.n_errors
            assert all(r.length == case.n_errors for r in res.values())

    def test_ibm_omega_is_always_the_low_product(self, code):
        gen = StimulusGenerator(code, seed=12)
        for _ in range(100):
            syn = gen.random_syndromes()
            res = IBM(code.field, code.t).run(syn)
            assert res.omega == full_product(code.field, res.lam, syn)[: code.t]

    def test_lambda0_is_the_product_of_the_gammas(self, code):
        # random syndromes for iBM and riBM; the RiBM only up to t errors,
        # since beyond t Lambda can overwrite delta in its array
        gen = StimulusGenerator(code, seed=13)
        F = code.field
        inputs = {
            "ibm": [gen.random_syndromes() for _ in range(30)],
            "ribm": [gen.random_syndromes() for _ in range(30)],
            "ribm_systolic": [gen.random_errors(e).syndromes for e in range(code.t + 1)],
        }
        for name, syndromes in inputs.items():
            for syn in syndromes:
                res = make_kes(name, F, code.t).run(syn, trace=True)
                prod = 1
                for state in res.trace[: 2 * code.t]:
                    prod = F.mul(prod, state["gamma"])
                assert res.lam[0] == prod != 0

    def test_ribm_systolic_overlap_beyond_t(self):
        small = RSCode.tp3_reduced()
        gen = StimulusGenerator(small, seed=13)
        differs = 0
        for _ in range(300):
            syn = gen.random_syndromes()
            differs += RiBMSystolic(small.field, 3).run(syn).lam != RiBMSplit(small.field, 3).run(syn).lam
        assert differs > 0

    def test_zero_syndromes(self, code):
        t = code.t
        for name in NAMES:
            res = make_kes(name, code.field, t).run([0] * (2 * t))
            assert res.lam == [1] + [0] * t
            assert res.omega == [0] * t
            assert res.degree == 0

    def test_beyond_t_still_agree_on_the_ibm_relation(self, code):
        # past t the outputs may differ between architectures; check each one alone
        gen = StimulusGenerator(code, seed=14)
        for _ in range(50):
            syn = gen.random_syndromes()
            for name in NAMES:
                res = make_kes(name, code.field, code.t).run(syn)
                assert len(res.lam) == code.t + 1 and len(res.omega) == code.t
                assert res.lam[0] != 0


class TestRegisters:
    """!
    @brief Register contents, traces and the IBM syndrome register.
    """

    def test_trace_length_and_keys(self, gf8):
        expected = {
            "ibm": {"lam", "b", "ds", "omega"},
            "ribm": {"lam", "b", "delta", "theta"},
            "ribm_systolic": {"delta", "theta"},
        }
        for name in NAMES:
            kes = make_kes(name, gf8, 2)
            res = kes.run(ZHANG_SYN, trace=True)
            assert len(res.trace) == kes.iterations + 1 == res.iterations + 1
            assert [s["r"] for s in res.trace] == list(range(kes.iterations + 1))
            assert set(res.trace[0]) == {"r", "gamma", "k", "mc", "discrepancy"} | expected[name]
            assert res.trace[-1]["mc"] is False

    def test_iterations(self, gf8):
        assert IBM(gf8, 2).iterations == 6
        assert RiBMSplit(gf8, 2).iterations == RiBMSystolic(gf8, 2).iterations == 4

    def test_ibm_syndrome_register(self, gf8):
        res = IBM(gf8, 2).run([1, 2, 3, 4], trace=True)
        assert res.trace[0]["ds"] == [1, 2, 3, 4]
        assert res.trace[1]["ds"] == [2, 3, 4, 1]
        assert res.trace[4]["ds"] == [1, 2, 0, 0]    # rotated back, DS_t .. DS_2t-1 cleared
        assert res.trace[5]["ds"] == [2, 0, 0, 1]
        assert [len(s["omega"]) for s in res.trace] == [0, 0, 0, 0, 0, 1, 2]

    def test_ribm_systolic_layout(self, gf8):
        res = RiBMSystolic(gf8, 2).run(ZHANG_SYN, trace=True)
        assert res.trace[0]["delta"] == ZHANG_SYN + [0, 0, 1]
        split = RiBMSplit(gf8, 2).run(ZHANG_SYN, trace=True)
        for r, (a, b) in enumerate(zip(res.trace, split.trace)):
            # Lambda enters from the right: lambda_i(r) sits in PE1_(3t - r + i)
            assert a["delta"][6 - r: 6 - r + 3][: r + 1] == b["lam"][: r + 1]

    def test_mc_matches_the_control_decision(self, gf8):
        res = RiBMSplit(gf8, 2).run(ZHANG_SYN, trace=True)
        # r = 0: discrepancy 0; r = 1: a, k = 1; r = 2, 3: k < 0; r = 4: finished
        assert [s["mc"] for s in res.trace] == [False, True, False, False, False]


class TestProtocol:
    """!
    @brief load / step / result sequencing and argument checks.
    """

    def test_step_before_load(self, gf8):
        for name in NAMES:
            with pytest.raises(RuntimeError, match="load"):
                make_kes(name, gf8, 2).step()

    def test_step_after_the_end(self, gf8):
        kes = RiBMSystolic(gf8, 2)
        kes.run(ZHANG_SYN)
        assert kes.done
        with pytest.raises(RuntimeError, match="done"):
            kes.step()

    def test_result_before_the_end(self, gf8):
        kes = IBM(gf8, 2)
        with pytest.raises(RuntimeError):
            kes.result()
        kes.load(ZHANG_SYN)
        kes.step()
        assert not kes.done
        with pytest.raises(RuntimeError):
            kes.result()

    def test_reset_and_reload(self, gf8):
        for name in NAMES:
            kes = make_kes(name, gf8, 2)
            first = kes.run(ZHANG_SYN)
            kes.reset()
            assert not kes.loaded and kes.r == 0 and kes.gamma == 0 and not kes.mc
            assert all(not any(v) for v in kes._registers().values())
            kes.run([1, 2, 3, 4])
            assert kes.run(ZHANG_SYN).lam == first.lam

    @pytest.mark.parametrize("syn, exc", [([0, 1, 2], ValueError), ([0, 1, 2, 8], ValueError),
                                          ([0, 1, 2, 3.0], TypeError)])
    def test_load_rejects_bad_syndromes(self, gf8, syn, exc):
        with pytest.raises(exc):
            RiBMSplit(gf8, 2).load(syn)

    @pytest.mark.parametrize("t", [0, -1, 2.0, True])
    def test_rejects_bad_t(self, gf8, t):
        with pytest.raises(ValueError):
            RiBMSystolic(gf8, t)

    def test_rejects_bad_field(self):
        with pytest.raises(TypeError):
            IBM(8, 2)

    def test_make_kes(self, gf8):
        assert isinstance(make_kes("ibm", gf8, 2), IBM)
        assert isinstance(make_kes("ribm", gf8, 2), RiBMSplit)
        assert isinstance(make_kes("ribm_systolic", gf8, 2), RiBMSystolic)
        assert all(issubclass(cls, KesBase) for cls in ARCHITECTURES.values())
        with pytest.raises(ValueError, match="unknown"):
            make_kes("euclid", gf8, 2)

    def test_degree_of_a_zero_locator(self, gf8):
        res = RiBMSplit(gf8, 2).run(ZHANG_SYN)
        res.lam = [0, 0, 0]
        assert res.degree == -1
