"""!
@file test_rs_code.py
@brief Unit tests for rs_code.py.

External reference values: the generator of ISO/IEC 18004 Annex A for 16 EC
codewords, the parity of TP3 exercise 6, the control value of the reduced
generator and Examples 29-30 of Zhang, "VLSI Architectures for Modern Error
Correcting Codes" (RS(7, 3) over GF(2^3), first root alpha).
"""

import random

import pytest

from galois import GaloisField, GFPoly
from rs_code import RSCode

ANNEX_A_16 = [0, 120, 104, 107, 109, 102, 161, 76, 3, 91, 191, 147, 169, 182, 194, 225, 120]
QR_PARITY = [0xE8, 0x9B, 0xBB, 0x6A, 0x6F, 0x00, 0x3E, 0x3E,
             0xC9, 0x7B, 0xC7, 0x1E, 0x8E, 0x75, 0x65, 0x38]

# GF(2^3), x^3 + x + 1: alpha^0..alpha^6 = 1, 2, 4, 3, 6, 7, 5
A1, A2, A3, A4, A6 = 2, 4, 3, 6, 5


@pytest.fixture(scope="module")
def qr():
    return RSCode.qr_2m()


@pytest.fixture(scope="module")
def small():
    return RSCode.tp3_reduced()


@pytest.fixture(scope="module")
def zhang():
    """!
    @brief RS(7, 3) of Zhang, Examples 29-31: roots alpha^1 .. alpha^4.
    """
    return RSCode(GaloisField(3, 0b011), 7, 3, first_root=1)


@pytest.fixture
def rng():
    return random.Random(7)


class TestParameters:
    """!
    @brief Construction and derived parameters.
    """

    def test_qr(self, qr):
        assert (qr.n, qr.k, qr.nsym, qr.t, qr.d_min) == (44, 28, 16, 8, 17)
        assert qr.first_root == 0
        assert qr.field == GaloisField(8, 0x1D)
        assert qr.rate == pytest.approx(28 / 44)

    def test_reduced(self, small):
        assert (small.n, small.k, small.t, small.d_min) == (15, 9, 3, 7)
        assert small.field == GaloisField(4, 0x3)

    def test_odd_redundancy(self):
        code = RSCode(GaloisField(4, 0x3), 15, 10)
        assert code.nsym == 5 and code.t == 2

    @pytest.mark.parametrize("n, k, exc", [
        (1, 1, ValueError), (16, 9, ValueError), (15, 0, ValueError),
        (15, 15, ValueError), (15.0, 9, TypeError), (15, True, TypeError),
    ])
    def test_invalid_lengths(self, n, k, exc):
        with pytest.raises(exc):
            RSCode(GaloisField(4, 0x3), n, k)

    def test_invalid_field_and_root(self):
        with pytest.raises(TypeError):
            RSCode(16, 15, 9)
        with pytest.raises(TypeError):
            RSCode(GaloisField(4, 0x3), 15, 9, first_root="0")

    def test_text(self, qr):
        assert str(qr) == "RS(44, 28) over GF(2^8), t = 8"
        assert repr(qr) == "RSCode(GaloisField(m=8, poly=0x1d), n=44, k=28, first_root=0)"


class TestGenerator:
    """!
    @brief Generator polynomials against published values.
    """

    def test_qr_matches_annex_a(self, qr):
        assert [qr.field.log(c) for c in qr.generator.coeffs] == ANNEX_A_16

    def test_reduced_control_value(self, small):
        assert list(small.generator.coeffs) == [1, 0xA, 0xF, 2, 4, 3, 1]

    def test_zhang_example_29(self, zhang):
        # g(x) = x^4 + a^3 x^3 + x^2 + a x + a^3
        assert list(zhang.generator.coeffs) == [1, A3, 1, A1, A3]

    @pytest.mark.parametrize("name", ["qr", "small", "zhang"])
    def test_roots(self, name, request):
        code = request.getfixturevalue(name)
        g = code.generator
        assert g.degree == code.nsym and g.lead == 1
        assert code.roots == [code.field.exp(code.first_root + i) for i in range(code.nsym)]
        assert all(g(r) == 0 for r in code.roots)


class TestEncoding:
    """!
    @brief Systematic encoding.
    """

    def test_reference_vector(self, qr):
        msg = list(range(28))
        assert qr.parity(msg) == QR_PARITY
        assert qr.encode(msg) == msg + QR_PARITY

    @pytest.mark.parametrize("name", ["qr", "small", "zhang"])
    def test_codewords(self, name, request, rng):
        code = request.getfixturevalue(name)
        for _ in range(50):
            msg = [rng.randrange(code.field.order) for _ in range(code.k)]
            cw = code.encode(msg)
            assert len(cw) == code.n
            assert code.message(cw) == msg
            assert code.is_codeword(cw)
            assert (GFPoly(code.field, cw) % code.generator).is_zero

    def test_linearity(self, qr, rng):
        for _ in range(20):
            a = [rng.randrange(256) for _ in range(28)]
            b = [rng.randrange(256) for _ in range(28)]
            ab = [x ^ y for x, y in zip(a, b)]
            assert qr.encode(ab) == [x ^ y for x, y in zip(qr.encode(a), qr.encode(b))]

    def test_leading_zero_parity_keeps_its_length(self, qr, rng):
        found = 0
        for _ in range(2000):
            msg = [rng.randrange(256) for _ in range(28)]
            par = qr.parity(msg)
            assert len(par) == 16
            found += par[0] == 0
        assert found > 0

    def test_reduced_exercise_2(self, small):
        cw = small.encode(small.unpack(0x52D09F17C, 9))
        assert small.message(cw) == small.unpack(0x52D09F17C, 9)
        for j in (0, 5):
            assert GFPoly(small.field, cw)(small.field.exp(j)) == 0

    @pytest.mark.parametrize("msg, exc", [
        ([0] * 27, ValueError), ([0] * 29, ValueError),
        ([0] * 27 + [256], ValueError), ([0] * 27 + [1.0], TypeError),
    ])
    def test_rejects_bad_messages(self, qr, msg, exc):
        with pytest.raises(exc):
            qr.encode(msg)


class TestSyndromes:
    """!
    @brief Syndrome computation.
    """

    def test_single_error(self, qr, rng):
        F = qr.field
        cw = qr.encode([rng.randrange(256) for _ in range(28)])
        for p in (0, 1, 15, 16, 43):
            e = rng.randrange(1, 256)
            syn = qr.syndromes(qr.add_errors(cw, {p: e}))
            assert syn == [F.mul(e, F.exp(j * p)) for j in range(16)]

    def test_only_depend_on_the_errors(self, small, rng):
        cw = small.encode([rng.randrange(16) for _ in range(9)])
        errors = {2: 5, 9: 1, 14: 12}
        assert small.syndromes(small.add_errors(cw, errors)) == small.syndromes(small.error_word(errors))

    def test_zhang_example_30(self, zhang):
        # r(x) = 1 + a^4 x + a^2 x^2 + a^2 x^3 + a^4 x^4 + a x^5
        received = [0, A1, A4, A2, A2, A4, 1]
        assert zhang.syndromes(received) == [0, A1, A2, A6]
        assert not zhang.is_codeword(received)
        assert zhang.is_codeword(zhang.add_errors(received, {0: 1, 3: A4}))

    @pytest.mark.parametrize("word, exc", [([0] * 43, ValueError), ([0] * 43 + [-1], ValueError)])
    def test_rejects_bad_words(self, qr, word, exc):
        with pytest.raises(exc):
            qr.syndromes(word)


class TestErrors:
    """!
    @brief Error patterns and positions.
    """

    def test_position_is_the_degree(self, small):
        assert small.error_word({0: 7}) == [0] * 14 + [7]
        assert small.error_word({14: 3, 1: 2}) == [3] + [0] * 12 + [2, 0]
        assert small.error_word({}) == [0] * 15

    def test_add_errors(self, small):
        cw = small.encode([1] * 9)
        rx = small.add_errors(cw, {14: 1})
        assert rx[0] == 0 and rx[1:] == cw[1:]

    @pytest.mark.parametrize("errors, exc", [
        ({15: 1}, ValueError), ({-1: 1}, ValueError), ({True: 1}, ValueError),
        ({"3": 1}, ValueError), ({3: 16}, ValueError), ({3: 1.0}, TypeError),
    ])
    def test_rejects_bad_errors(self, small, errors, exc):
        with pytest.raises(exc):
            small.error_word(errors)


class TestPacking:
    """!
    @brief Bus packing: the degree-i symbol in bits [m*i + m - 1 : m*i].
    """

    def test_qr_buses(self, qr):
        msg = list(range(28))
        i_msg = qr.pack(msg)
        o_cw = qr.pack(qr.encode(msg))
        assert i_msg >> 216 == 0x00 and i_msg & 0xFF == 0x1B
        assert o_cw >> 128 == i_msg
        assert o_cw & ((1 << 128) - 1) == qr.pack(QR_PARITY)

    def test_ascending(self, qr):
        assert qr.pack([0x01, 0x02, 0x03], ascending=True) == 0x030201
        assert qr.pack([0x01, 0x02, 0x03]) == 0x010203
        assert qr.unpack(0x030201, 3, ascending=True) == [1, 2, 3]
        assert qr.unpack(0x030201, 4) == [0, 3, 2, 1]

    def test_round_trip(self, small, rng):
        for _ in range(100):
            word = [rng.randrange(16) for _ in range(15)]
            assert small.unpack(small.pack(word), 15) == word
            assert small.unpack(small.pack(word, ascending=True), 15, ascending=True) == word

    def test_errors(self, qr):
        with pytest.raises(ValueError):
            qr.pack([])
        with pytest.raises(ValueError):
            qr.unpack(1 << 24, 3)
