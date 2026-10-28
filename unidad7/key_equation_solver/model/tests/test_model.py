"""!
@file test_model.py
@brief The C++ model (kesmodel) against galois.py and kes_ref.py.

The traces are compared state by state, including every register and the
control signals, for correctable cases and for random syndromes beyond t.
"""

import pytest

import kesmodel
from galois import GaloisField
from kes_ref import ARCHITECTURES, make_kes
from rs_code import RSCode
from stimuli import StimulusGenerator

NAMES = sorted(ARCHITECTURES)
FIELDS = [(1, 0x1), (3, 0b011), (4, 0x3), (8, 0x1D), (8, 0x2B)]


def cpp_field(field):
    return kesmodel.GaloisField(field.m, field.poly)


@pytest.fixture(scope="module", params=["qr", "small", "zhang"])
def code(request):
    if request.param == "qr":
        return RSCode.qr_2m()
    if request.param == "small":
        return RSCode.tp3_reduced()
    return RSCode(GaloisField(3, 0b011), 7, 3, first_root=1)


class TestField:
    """!
    @brief GF(2^m) arithmetic against galois.py.
    """

    @pytest.mark.parametrize("m, poly", FIELDS)
    def test_matches_galois(self, m, poly):
        ref, cpp = GaloisField(m, poly), kesmodel.GaloisField(m, poly)
        assert (cpp.m, cpp.poly, cpp.order, cpp.n) == (m, poly, ref.order, ref.n)
        elems = range(ref.order)
        assert all(cpp.mul(a, b) == ref.mul(a, b) for a in elems for b in elems)
        assert all(cpp.inverse(a) == ref.inverse(a) for a in range(1, ref.order))
        assert all(cpp.log(a) == ref.log(a) for a in range(1, ref.order))
        assert all(cpp.exp(k) == ref.exp(k) for k in range(-2 * ref.n, 2 * ref.n + 1))

    def test_larger_field(self):
        ref, cpp = GaloisField(16, 0x100B), kesmodel.GaloisField(16, 0x100B)
        assert all(cpp.exp(k) == ref.exp(k) for k in range(0, ref.n, 97))
        assert cpp.mul(0x1234, 0xBEEF) == ref.mul(0x1234, 0xBEEF)

    @pytest.mark.parametrize("m, poly", [(0, 1), (17, 1), (4, 16), (4, 0b0010), (4, 0b1111), (4, 0b0101)])
    def test_rejects_bad_fields(self, m, poly):
        with pytest.raises(ValueError):
            kesmodel.GaloisField(m, poly)

    def test_errors(self):
        F = kesmodel.GaloisField(8, 0x1D)
        with pytest.raises(ZeroDivisionError):
            F.inverse(0)
        with pytest.raises(ZeroDivisionError):
            F.log(0)
        for call in (lambda: F.mul(256, 1), lambda: F.inverse(256), lambda: F.log(300)):
            with pytest.raises(ValueError):
                call()
        assert "m=8" in repr(F)


@pytest.mark.parametrize("name", NAMES)
class TestTraces:
    """!
    @brief Every state of every architecture, compared with kes_ref.
    """

    def check(self, name, code, syndromes):
        ref = make_kes(name, code.field, code.t).run(syndromes, trace=True)
        cpp = kesmodel.make_kes(name, cpp_field(code.field), code.t).run(syndromes, trace=True)
        assert cpp.trace == ref.trace
        assert (cpp.lam, cpp.omega, cpp.form, cpp.iterations, cpp.length, cpp.degree) == \
               (ref.lam, ref.omega, ref.form, ref.iterations, ref.length, ref.degree)

    def test_suite(self, name, code):
        gen = StimulusGenerator(code, seed=51)
        for case in gen.suite():
            self.check(name, code, case.syndromes[: 2 * code.t])

    def test_random_syndromes(self, name, code):
        gen = StimulusGenerator(code, seed=52)
        for _ in range(300):
            self.check(name, code, gen.random_syndromes())

    def test_zero_syndromes(self, name, code):
        self.check(name, code, [0] * (2 * code.t))

    def test_step_by_step(self, name, code):
        syn = StimulusGenerator(code, seed=53).random_errors(code.t).syndromes[: 2 * code.t]
        ref = make_kes(name, code.field, code.t)
        cpp = kesmodel.make_kes(name, cpp_field(code.field), code.t)
        assert cpp.snapshot() == ref.snapshot()
        ref.load(syn)
        cpp.load(syn)
        while not ref.done:
            assert cpp.snapshot() == ref.snapshot()
            assert (cpp.r, cpp.gamma, cpp.k, cpp.mc, cpp.discrepancy) == \
                   (ref.r, ref.gamma, ref.k, ref.mc, ref.discrepancy)
            ref.step()
            cpp.step()
        assert cpp.done and cpp.snapshot() == ref.snapshot()
        cpp.reset()
        ref.reset()
        assert cpp.snapshot() == ref.snapshot() and not cpp.loaded

    def test_untraced_run(self, name, code):
        res = kesmodel.make_kes(name, cpp_field(code.field), code.t).run([1] * (2 * code.t))
        assert res.trace is None


class TestInterface:
    """!
    @brief Properties and error handling of the solvers.
    """

    def test_properties(self):
        F = kesmodel.GaloisField(8, 0x1D)
        for name, form, n in (("ibm", "classic", 24), ("ribm", "high", 16), ("ribm_systolic", "high", 16)):
            kes = kesmodel.make_kes(name, F, 8)
            assert (kes.name, kes.form, kes.iterations, kes.t) == (name, form, n, 8)
            assert not kes.loaded and not kes.done and kes.r == 0
        assert kesmodel.ARCHITECTURES == ("ibm", "ribm", "ribm_systolic")

    def test_errors(self):
        F = kesmodel.GaloisField(8, 0x1D)
        kes = kesmodel.make_kes("ribm_systolic", F, 8)
        with pytest.raises(RuntimeError, match="load"):
            kes.step()
        with pytest.raises(RuntimeError):
            kes.result()
        with pytest.raises(ValueError):
            kes.load([0] * 15)
        with pytest.raises(ValueError):
            kes.load([0] * 15 + [256])
        with pytest.raises(TypeError):
            kes.load([0] * 15 + [-1])
        kes.run([0] * 16)
        with pytest.raises(RuntimeError, match="done"):
            kes.step()
        with pytest.raises(ValueError, match="unknown"):
            kesmodel.make_kes("euclid", F, 8)
        with pytest.raises(ValueError):
            kesmodel.make_kes("ibm", F, 0)

    def test_qr_reference_vector(self):
        qr = RSCode.qr_2m()
        syn = qr.syndromes(qr.add_errors(qr.encode(list(range(28))), {40: 0x55, 3: 0xA7}))
        F = cpp_field(qr.field)
        for name in NAMES:
            res = kesmodel.make_kes(name, F, 8).run(syn)
            assert qr.pack(res.lam, ascending=True) == 0x6674BB
        res = kesmodel.make_kes("ribm_systolic", F, 8).run(syn)
        assert qr.pack(res.omega, ascending=True) == 0xB837
