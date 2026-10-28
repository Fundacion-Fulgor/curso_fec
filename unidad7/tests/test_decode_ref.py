"""!
@file test_decode_ref.py
@brief Unit tests for decode_ref.py.
"""

import pytest

from galois import GaloisField, GFPoly
from decode_ref import chien_forney, decode, forney_exponent
from kes_ref import ARCHITECTURES
from rs_code import RSCode
from stimuli import StimulusGenerator

NAMES = sorted(ARCHITECTURES)

# GF(2^3), x^3 + x + 1: alpha^0..alpha^6 = 1, 2, 4, 3, 6, 7, 5
A1, A2, A3, A4, A5, A6 = 2, 4, 3, 6, 7, 5


@pytest.fixture(scope="module", params=["qr", "small"])
def code(request):
    return RSCode.qr_2m() if request.param == "qr" else RSCode.tp3_reduced()


@pytest.fixture(scope="module")
def zhang():
    return RSCode(GaloisField(3, 0b011), 7, 3, first_root=1)


class TestForneyExponent:
    """!
    @brief Exponent of X in eqs. (7) and (12).
    """

    @pytest.mark.parametrize("form, b, t, expected", [
        ("classic", 0, 8, 1), ("classic", 1, 8, 0), ("high", 0, 8, -15),
        ("high", 1, 2, -4), ("high", -16, 8, 1),
    ])
    def test_values(self, form, b, t, expected):
        assert forney_exponent(form, b, t) == expected

    def test_unknown_form(self):
        with pytest.raises(ValueError):
            forney_exponent("low", 0, 8)


class TestChienForney:
    """!
    @brief Chien search and error values on known locators.
    """

    def test_zhang_example_31_high(self, zhang):
        cf = chien_forney(zhang.field, 7, [A2, A3, A5], [A6, A4], "high", first_root=1)
        assert (cf.positions, cf.values, cf.fail) == ([0, 3], [1, A4], False)
        assert cf.errors == zhang.error_word({0: 1, 3: A4})

    def test_zhang_example_30_classic(self, zhang):
        # normalized Lambda = 1 + a x + a^3 x^2 and Omega = a x (Example 30)
        cf = chien_forney(zhang.field, 7, [1, A1, A3], [0, A1], "classic", first_root=1)
        assert (cf.positions, cf.values, cf.fail) == ([0, 3], [1, A4], False)

    def test_scaling_does_not_matter(self, zhang):
        F = zhang.field
        for s in range(1, 8):
            lam = [F.mul(s, c) for c in [1, A1, A3]]
            omega = [F.mul(s, c) for c in [0, A1]]
            cf = chien_forney(F, 7, lam, omega, "classic", first_root=1)
            assert cf.values == [1, A4]

    def test_repeated_root_fails(self, code):
        F = code.field
        X = F.exp(3)
        lam = (GFPoly(F, [X, 1]) ** 2).to_list(code.t + 1, ascending=True)
        cf = chien_forney(F, code.n, lam, [1] + [0] * (code.t - 1), "high")
        assert cf.positions == [3] and cf.fail

    def test_root_outside_a_shortened_code_fails(self):
        qr = RSCode.qr_2m()
        F = qr.field
        lam = [1, F.exp(50)] + [0] * 7    # root at position 50 > 43
        cf = chien_forney(F, 44, lam, [1] + [0] * 7, "high")
        assert cf.positions == [] and cf.fail
        assert not any(cf.errors)

    def test_zero_locator_fails(self, code):
        cf = chien_forney(code.field, code.n, [0] * (code.t + 1), [0] * code.t, "high")
        assert cf.fail and len(cf.positions) == code.n

    def test_constant_locator(self, code):
        cf = chien_forney(code.field, code.n, [5] + [0] * code.t, [0] * code.t, "classic")
        assert cf.positions == [] and not cf.fail and not any(cf.errors)

    @pytest.mark.parametrize("lam, omega, n, form", [
        ([1], [], 7, "high"), ([1, 0, 0], [0], 7, "high"), ([1, 0, 0], [0, 0], 8, "high"),
        ([1, 0, 0], [0, 0], 0, "high"), ([1, 0, 0], [0, 0], True, "high"),
        ([1, 0, 0], [0, 0], 7, "low"),
    ])
    def test_rejects_bad_arguments(self, zhang, lam, omega, n, form):
        with pytest.raises(ValueError):
            chien_forney(zhang.field, n, lam, omega, form)


class TestDecode:
    """!
    @brief Complete decoder with each architecture.
    """

    @pytest.mark.parametrize("name", NAMES)
    def test_correctable_cases(self, code, name):
        gen = StimulusGenerator(code, seed=21)
        cases = [c for c in gen.suite() if c.correctable]
        cases += list(gen.stream(60, {e: 1 for e in range(code.t + 1)}))
        for case in cases:
            res = decode(code, case.received, name)
            assert not res.uncorrectable
            assert res.codeword == case.codeword
            assert res.message == case.message
            assert res.errors == code.error_word(case.errors)
            assert res.positions == sorted(case.errors)
            assert res.values == [case.errors[p] for p in res.positions]
            assert res.corrected == (case.n_errors > 0)
            assert res.syndromes == case.syndromes
            assert res.consistent
            assert res.kes.degree == res.kes.length == case.n_errors

    @pytest.mark.parametrize("name", NAMES)
    def test_beyond_t(self, code, name):
        gen = StimulusGenerator(code, seed=22)
        cases = [c for c in gen.suite() if not c.correctable]
        cases += [gen.random_errors(e) for e in range(code.t + 1, 2 * code.t + 1) for _ in range(20)]
        flagged = 0
        for case in cases:
            res = decode(code, case.received, name)
            if res.uncorrectable:
                flagged += 1
                assert not res.corrected
            elif res.consistent:
                # miscorrection: a different codeword within distance t
                assert res.codeword != case.codeword
                assert len(res.positions) <= code.t
                assert res.kes.degree == res.kes.length
            else:
                # passes the root count without giving a codeword; deg Lambda != L
                assert res.kes.degree != res.kes.length
        # RS(15, 9) miscorrects about 9% of random patterns beyond t; RS(44, 28) almost never
        assert flagged >= (0.75 if code.field.m == 4 else 1.0) * len(cases)

    @pytest.mark.parametrize("name", NAMES)
    def test_tp3_exercise_3(self, name):
        small = RSCode.tp3_reduced()
        for word, n_errors in ((0x52D0AF17CB1EE64, 2), (0x56D09F171B1E06F, 3)):
            res = decode(small, small.unpack(word, 15), name)
            assert small.pack(res.message) == 0x52D09F17C
            assert len(res.positions) == n_errors
            assert res.corrected and not res.uncorrectable

    @pytest.mark.parametrize("name", NAMES)
    def test_zhang_example(self, zhang, name):
        res = decode(zhang, [0, A1, A4, A2, A2, A4, 1], name)
        assert (res.positions, res.values) == ([0, 3], [1, A4])
        assert zhang.is_codeword(res.codeword)

    def test_odd_redundancy(self):
        code = RSCode(GaloisField(4, 0x3), 15, 10)
        gen = StimulusGenerator(code, seed=23)
        for e in range(code.t + 1):
            case = gen.random_errors(e)
            assert decode(code, case.received).codeword == case.codeword

    def test_unknown_architecture(self, code):
        with pytest.raises(ValueError):
            decode(code, [0] * code.n, "euclid")

    def test_rejects_bad_words(self, code):
        with pytest.raises(ValueError):
            decode(code, [0] * (code.n - 1))
