#!/usr/bin/env python3
"""!
@file test_galois_unit.py
@brief Unit tests for galois.py (pytest).

Run as a script to get line and branch coverage (needs pytest-cov):
    python3 test_galois_unit.py
    python3 test_galois_unit.py -k Poly      # extra arguments go to pytest

Or call pytest directly:
    pytest test_galois_unit.py --cov=galois --cov-branch --cov-report=term-missing
"""

import os
import random
import subprocess
import sys

import pytest

from galois import GFPoly, GaloisField

PRIMITIVE = [(1, 0b1), (2, 0b11), (3, 0b011), (4, 0b0011), (4, 0b1001),
             (5, 0b00101), (6, 0b000011), (8, 0x1D), (10, 0x009)]

QR_PARITY = [0xE8, 0x9B, 0xBB, 0x6A, 0x6F, 0x00, 0x3E, 0x3E,
             0xC9, 0x7B, 0xC7, 0x1E, 0x8E, 0x75, 0x65, 0x38]


@pytest.fixture(params=[(3, 0b011), (4, 0b0011), (4, 0b1001)],
                ids=lambda p: f"GF(2^{p[0]})-{p[1]:#x}")
def small_field(request):
    """!
    @brief Fields small enough for exhaustive sweeps.
    """
    return GaloisField(*request.param)


@pytest.fixture(scope="module")
def gf16():
    """!
    @brief GF(2^4) with x^4 + x + 1 (reduced code of TP3).
    """
    return GaloisField(4, 0b0011)


@pytest.fixture(scope="module")
def gf256():
    """!
    @brief GF(2^8) with x^8 + x^4 + x^3 + x^2 + 1 (QR code).
    """
    return GaloisField(8, 0x1D)


@pytest.fixture
def rng():
    return random.Random(2026)


def rand_poly(F, rng, max_degree=6):
    """!
    @brief Random polynomial; the zero polynomial is included.
    """
    d = rng.randrange(-1, max_degree + 1)
    if d < 0:
        return GFPoly.zero(F)
    return GFPoly(F, [rng.randrange(1, F.order)] + [rng.randrange(F.order) for _ in range(d)])


def mult_order(F, a):
    k, acc = 1, a
    while acc != 1:
        acc = F.mul(acc, a)
        k += 1
    return k


def mat_vec(M, a):
    """!
    @brief Apply a binary matrix to the bits of a; bit i is row i dotted with a.
    """
    out = 0
    for i, row in enumerate(M):
        bit = 0
        for j, mij in enumerate(row):
            bit ^= mij & (a >> j) & 1
        out |= bit << i
    return out


def odd_part(p):
    asc = p.to_list(ascending=True)
    return GFPoly.from_ascending(p.field, [c if i % 2 else 0 for i, c in enumerate(asc)])


class TestFieldConstruction:
    """!
    @brief GaloisField constructor and table generation.
    """

    @pytest.mark.parametrize("m, poly", PRIMITIVE)
    def test_primitive_polynomials(self, m, poly):
        F = GaloisField(m, poly)
        assert F.is_primitive
        assert F.alpha == F.generator == (2 if m > 1 else 1)
        assert len(F) == F.order == 2 ** m
        assert F.n == 2 ** m - 1
        assert sorted(F.exp(i) for i in range(F.n)) == list(range(1, F.order))
        assert list(F.elements) == list(range(F.order))

    @pytest.mark.parametrize("m, poly, exc", [
        (0, 1, ValueError),
        (-3, 1, ValueError),
        (2.5, 1, ValueError),
        (4, "0011", TypeError),
        (4, 0b10011, ValueError),
        (4, -1, ValueError),
        (4, 0b0010, ValueError),
    ], ids=["m=0", "m<0", "m-float", "poly-str", "poly-too-wide", "poly<0", "poly-even"])
    def test_invalid_arguments(self, m, poly, exc):
        with pytest.raises(exc):
            GaloisField(m, poly)

    def test_non_primitive_rejected_by_default(self):
        with pytest.raises(ValueError, match="not primitive"):
            GaloisField(4, 0b1111)

    def test_non_primitive_accepted_on_request(self):
        F = GaloisField(4, 0b1111, require_primitive=False)
        assert not F.is_primitive
        assert mult_order(F, 2) == 5
        assert F.generator != 2
        assert mult_order(F, F.generator) == F.n
        assert all(F.mul(a, F.inverse(a)) == 1 for a in range(1, F.order))

    def test_aes_field(self):
        F = GaloisField(8, 0x1B, require_primitive=False)
        assert F.generator == 3
        assert F.mul(0x57, 0x83) == F.mul_ref(0x57, 0x83) == 0xC1

    @pytest.mark.parametrize("m, poly", [(2, 0b01), (4, 0b0101), (4, 0b0111), (6, 0b000101)])
    def test_reducible_rejected(self, m, poly):
        with pytest.raises(ValueError, match="not irreducible"):
            GaloisField(m, poly, require_primitive=False)


class TestFieldValidation:
    """!
    @brief check() and membership.
    """

    @pytest.mark.parametrize("bad, exc", [
        (16, ValueError), (-1, ValueError), (3.0, TypeError),
        ("3", TypeError), (True, TypeError), (None, TypeError),
    ])
    def test_check_rejects(self, gf16, bad, exc):
        with pytest.raises(exc):
            gf16.check(bad)

    def test_check_returns_the_element(self, gf16):
        assert all(gf16.check(a) == a for a in gf16.elements)

    def test_contains(self, gf16):
        assert all(a in gf16 for a in range(16))
        for bad in (16, -1, True, 1.0, "1"):
            assert bad not in gf16


class TestFieldArithmetic:
    """!
    @brief Field operations, checked exhaustively on small fields.
    """

    def test_addition(self, small_field):
        F = small_field
        for a in F.elements:
            assert F.add(a, 0) == a
            assert F.add(a, a) == 0
            assert F.neg(a) == a
            for b in F.elements:
                assert F.add(a, b) == F.add(b, a) == F.sub(a, b) == a ^ b

    def test_product_axioms(self, small_field):
        F = small_field
        for a in F.elements:
            assert F.mul(a, 1) == a
            assert F.mul(a, 0) == F.mul(0, a) == 0
            for b in F.elements:
                assert F.mul(a, b) == F.mul(b, a)
                for c in F.elements:
                    assert F.mul(F.mul(a, b), c) == F.mul(a, F.mul(b, c))
                    assert F.mul(a, b ^ c) == F.mul(a, b) ^ F.mul(a, c)

    def test_no_zero_divisors(self, small_field):
        F = small_field
        nz = range(1, F.order)
        assert all(F.mul(a, b) != 0 for a in nz for b in nz)

    def test_mul_matches_reference(self, gf256):
        F = gf256
        assert all(F.mul(a, b) == F.mul_ref(a, b) for a in F.elements for b in F.elements)

    def test_inverse(self, small_field):
        F = small_field
        for a in range(1, F.order):
            assert F.mul(a, F.inverse(a)) == 1
            assert F.inverse(F.inverse(a)) == a
        with pytest.raises(ZeroDivisionError):
            F.inverse(0)

    def test_division(self, small_field):
        F = small_field
        for a in F.elements:
            for b in range(1, F.order):
                assert F.div(a, b) == F.mul(a, F.inverse(b))
                assert F.mul(F.div(a, b), b) == a
        for a in F.elements:
            with pytest.raises(ZeroDivisionError):
                F.div(a, 0)

    def test_power_matches_repeated_product(self, gf256, rng):
        F = gf256
        for _ in range(300):
            a, k = rng.randrange(F.order), rng.randrange(0, 40)
            acc = 1
            for _ in range(k):
                acc = F.mul(acc, a)
            assert F.power(a, k) == acc

    def test_power_identities(self, gf16):
        F = gf16
        nz = range(1, F.order)
        assert all(F.power(a, F.n) == 1 for a in nz)
        assert all(F.power(a, F.order) == a for a in F.elements)
        assert all(F.power(a, -3) == F.power(F.inverse(a), 3) for a in nz)
        assert all(F.power(a, 0) == 1 for a in F.elements)

    def test_power_of_zero(self, gf16):
        assert gf16.power(0, 0) == 1
        assert all(gf16.power(0, k) == 0 for k in range(1, 5))
        with pytest.raises(ZeroDivisionError):
            gf16.power(0, -1)

    @pytest.mark.parametrize("k", [1.0, "2", True])
    def test_power_rejects_bad_exponent(self, gf16, k):
        with pytest.raises(TypeError):
            gf16.power(3, k)

    def test_exp_and_log(self, gf256):
        F = gf256
        assert all(F.log(F.exp(i)) == i for i in range(F.n))
        assert all(F.exp(F.log(a)) == a for a in range(1, F.order))
        assert F.exp(F.n) == F.exp(0) == 1
        assert F.exp(-1) == F.inverse(F.alpha)
        with pytest.raises(ZeroDivisionError):
            F.log(0)

    def test_frobenius(self, gf256):
        F = gf256
        assert all(F.power(a ^ b, 2) == F.power(a, 2) ^ F.power(b, 2)
                   for a in F.elements for b in range(0, F.order, 5))

    def test_subfields(self, gf256):
        F = gf256
        for d, size in [(1, 2), (2, 4), (4, 16), (8, 256)]:
            assert sum(F.power(a, 2 ** d) == a for a in F.elements) == size

    def test_gf2(self):
        B = GaloisField(1, 0b1)
        assert [B.add(a, b) for a in (0, 1) for b in (0, 1)] == [0, 1, 1, 0]
        assert [B.mul(a, b) for a in (0, 1) for b in (0, 1)] == [0, 0, 0, 1]
        assert B.inverse(1) == 1


class TestHardwareHelpers:
    """!
    @brief inv_or_zero, mul_matrix and mul_masks.
    """

    def test_inv_or_zero(self, small_field):
        F = small_field
        assert F.inv_or_zero(0) == 0
        assert all(F.inv_or_zero(a) == F.inverse(a) for a in range(1, F.order))
        with pytest.raises(ValueError):
            F.inv_or_zero(F.order)

    def test_mul_matrix_computes_the_product(self, small_field):
        F = small_field
        for c in F.elements:
            M = F.mul_matrix(c)
            assert len(M) == F.m and all(len(row) == F.m for row in M)
            assert all(mat_vec(M, a) == F.mul(c, a) for a in F.elements)

    def test_mul_matrix_columns(self, gf256):
        F = gf256
        for c in (0x01, 0x02, 0x1D, 0x8E, 0xFF):
            M = F.mul_matrix(c)
            for j in range(F.m):
                assert sum(M[i][j] << i for i in range(F.m)) == F.mul(c, 1 << j)

    def test_mul_matrix_identity_and_zero(self, gf256):
        identity = tuple(tuple(int(i == j) for j in range(8)) for i in range(8))
        assert gf256.mul_matrix(1) == identity
        assert gf256.mul_matrix(0) == tuple((0,) * 8 for _ in range(8))

    def test_mul_matrix_is_a_ring_homomorphism(self, gf16, rng):
        F = gf16
        for _ in range(100):
            c1, c2 = rng.randrange(F.order), rng.randrange(F.order)
            M1, M2 = F.mul_matrix(c1), F.mul_matrix(c2)
            added = tuple(tuple(x ^ y for x, y in zip(r1, r2)) for r1, r2 in zip(M1, M2))
            product = tuple(
                tuple(sum(M1[i][k] & M2[k][j] for k in range(F.m)) & 1 for j in range(F.m))
                for i in range(F.m))
            assert added == F.mul_matrix(c1 ^ c2)
            assert product == F.mul_matrix(F.mul(c1, c2))

    def test_mul_masks(self, gf256):
        F = gf256
        for c in F.elements:
            masks = F.mul_masks(c)
            assert masks == tuple(sum(b << j for j, b in enumerate(r)) for r in F.mul_matrix(c))
            for a in F.elements:
                got = sum((bin(a & k).count("1") & 1) << i for i, k in enumerate(masks))
                assert got == F.mul(c, a)

    def test_mul_masks_known_values(self, gf256):
        assert gf256.mul_masks(0x01) == tuple(1 << i for i in range(8))
        assert gf256.mul_masks(0x02) == (0x80, 0x01, 0x82, 0x84, 0x88, 0x10, 0x20, 0x40)

    def test_mul_matrix_rejects_non_elements(self, gf16):
        with pytest.raises(ValueError):
            gf16.mul_matrix(16)
        with pytest.raises(TypeError):
            gf16.mul_masks(1.0)


class TestFieldRepresentation:
    """!
    @brief String forms, equality and hashing of fields.
    """

    def test_str_and_repr(self, gf256):
        assert str(gf256) == "GF(2^8) with P(x) = x^8 + x^4 + x^3 + x^2 + 1"
        assert repr(gf256) == "GaloisField(m=8, poly=0x1d)"
        assert repr(GaloisField(4, 3)) == "GaloisField(m=4, poly=0x3)"
        assert eval(repr(gf256), {"GaloisField": GaloisField}) == gf256

    def test_str_of_non_primitive_field(self):
        F = GaloisField(4, 0b1111, require_primitive=False)
        assert str(F).endswith(f"(not primitive, alpha = {F.generator})")

    def test_poly_str(self, gf16):
        assert gf16.poly_str(0) == "0"
        assert gf16.poly_str(1) == "1"
        assert gf16.poly_str(2) == "x"
        assert gf16.poly_str(0b1011) == "x^3 + x + 1"
        assert gf16.primitive_poly_str() == "x^4 + x + 1"
        with pytest.raises(ValueError):
            gf16.poly_str(16)

    def test_equality_and_hash(self):
        assert GaloisField(4, 3) == GaloisField(4, 3)
        assert GaloisField(4, 3) != GaloisField(4, 9)
        assert GaloisField(4, 3) != GaloisField(3, 3)
        assert GaloisField(4, 3) != "GF(16)"
        assert len({GaloisField(4, 3), GaloisField(4, 3), GaloisField(4, 9)}) == 2


class TestPolyConstruction:
    """!
    @brief GFPoly constructor, normalization and named constructors.
    """

    @pytest.mark.parametrize("given, stored, degree", [
        ([0, 0, 1, 5], (1, 5), 1),
        ([1, 0, 0], (1, 0, 0), 2),
        ([0, 0, 0], (0,), -1),
        ([], (0,), -1),
        ([0], (0,), -1),
        ([7], (7,), 0),
        ((3, 2, 1), (3, 2, 1), 2),
        (range(4), (1, 2, 3), 2),
    ])
    def test_normalization(self, gf16, given, stored, degree):
        p = GFPoly(gf16, given)
        assert p.coeffs == stored
        assert p.degree == degree
        assert len(p) == len(stored)
        assert p.is_zero == (degree == -1)

    def test_rejects_bad_coefficients(self, gf16):
        with pytest.raises(ValueError, match="position 1"):
            GFPoly(gf16, [1, 16])
        with pytest.raises(ValueError):
            GFPoly(gf16, [1, -1])
        with pytest.raises(TypeError, match="degree 0"):
            GFPoly(gf16, [1, 2.0])
        with pytest.raises(TypeError):
            GFPoly(gf16, [1, "a"])

    def test_rejects_bad_arguments(self, gf16):
        with pytest.raises(TypeError):
            GFPoly(16, [1, 2])
        with pytest.raises(TypeError):
            GFPoly(gf16, 5)
        with pytest.raises(TypeError):
            GFPoly(gf16, True)

    def test_named_constructors(self, gf16):
        assert GFPoly.zero(gf16).is_zero
        assert GFPoly.one(gf16).coeffs == (1,)
        assert GFPoly.monomial(gf16, 3).coeffs == (1, 0, 0, 0)
        assert GFPoly.monomial(gf16, 2, 7).coeffs == (7, 0, 0)
        assert GFPoly.monomial(gf16, 0, 5).coeffs == (5,)
        with pytest.raises(ValueError):
            GFPoly.monomial(gf16, -1)

    def test_from_roots(self, gf16, rng):
        F = gf16
        assert GFPoly.from_roots(F, []) == GFPoly.one(F)
        for _ in range(100):
            roots = rng.sample(range(F.order), rng.randint(1, 6))
            p = GFPoly.from_roots(F, roots)
            assert p.degree == len(roots) and p.lead == 1
            assert sorted(p.roots()) == sorted(roots)
            assert all((p % GFPoly(F, [1, r])).is_zero for r in roots)
        with pytest.raises(ValueError):
            GFPoly.from_roots(F, [16])

    def test_coefficient_access(self, gf16):
        p = GFPoly(gf16, [3, 0, 5, 1])
        assert [p.coeff(d) for d in range(4)] == [1, 5, 0, 3]
        assert p.coeff(4) == p.coeff(-1) == 0
        assert p.lead == 3
        assert GFPoly.zero(gf16).lead == 0
        assert GFPoly.zero(gf16).coeff(0) == 0


class TestPolyConversions:
    """!
    @brief Fixed-length lists, ascending order and packed integers.
    """

    def test_to_list(self, gf16):
        p = GFPoly(gf16, [3, 0, 5, 1])
        assert p.to_list() == [3, 0, 5, 1]
        assert p.to_list(ascending=True) == [1, 5, 0, 3]
        assert p.to_list(4) == [3, 0, 5, 1]
        assert p.to_list(6) == [0, 0, 3, 0, 5, 1]
        assert p.to_list(6, ascending=True) == [1, 5, 0, 3, 0, 0]
        assert GFPoly.zero(gf16).to_list() == [0]
        assert GFPoly.zero(gf16).to_list(3) == [0, 0, 0]

    @pytest.mark.parametrize("length", [3, 0, -1, True, 2.0])
    def test_to_list_rejects_bad_length(self, gf16, length):
        with pytest.raises(ValueError):
            GFPoly(gf16, [3, 0, 5, 1]).to_list(length)

    def test_from_ascending(self, gf16):
        assert GFPoly.from_ascending(gf16, [1, 5, 0, 3]) == GFPoly(gf16, [3, 0, 5, 1])
        assert GFPoly.from_ascending(gf16, (1, 0, 0)) == GFPoly.one(gf16)
        assert GFPoly.from_ascending(gf16, []).is_zero
        with pytest.raises(TypeError):
            GFPoly.from_ascending(gf16, 5)
        with pytest.raises(TypeError):
            GFPoly.from_ascending(gf16, False)

    def test_packing_layout(self, gf256):
        p = GFPoly(gf256, [0x00, 0x12, 0x34])
        assert len(p.coeffs) == 2
        assert p.to_list(3) == [0x00, 0x12, 0x34]
        assert p.to_int() == p.to_int(3) == 0x1234
        assert GFPoly.from_int(gf256, 0x1234, 3) == p
        assert GFPoly.from_int(gf256, 0, 5).is_zero
        assert GFPoly.from_int(gf256, 0xFFFF, 2).to_list(2) == [0xFF, 0xFF]
        assert GFPoly(GaloisField(4, 3), [0xA, 0xB, 0xC]).to_int() == 0xABC

    @pytest.mark.parametrize("value, length, exc", [
        (1 << 16, 2, ValueError),
        (-1, 2, ValueError),
        (1.0, 2, TypeError),
        (True, 2, TypeError),
        (0, 0, ValueError),
        (0, True, ValueError),
        (0, "2", ValueError),
    ])
    def test_from_int_rejects(self, gf256, value, length, exc):
        with pytest.raises(exc):
            GFPoly.from_int(gf256, value, length)

    def test_to_int_rejects_short_length(self, gf256):
        with pytest.raises(ValueError):
            GFPoly(gf256, [1, 0, 0]).to_int(2)

    def test_round_trips(self, gf16, gf256, rng):
        for _ in range(500):
            F = rng.choice((gf16, gf256))
            length = rng.randint(1, 44)
            word = [rng.randrange(F.order) for _ in range(length)]
            if rng.random() < 0.5:
                word[0] = 0
            p = GFPoly(F, word)
            assert p.to_list(length) == word
            assert GFPoly.from_ascending(F, word[::-1]) == p
            value = p.to_int(length)
            assert value < 1 << (F.m * length)
            assert GFPoly.from_int(F, value, length) == p
            assert all((value >> (F.m * i)) & (F.order - 1) == p.coeff(i) for i in range(length))


class TestPolyArithmetic:
    """!
    @brief Ring operations, scaling and powers.
    """

    def test_ring_axioms(self, gf16, rng):
        F = gf16
        zero, one = GFPoly.zero(F), GFPoly.one(F)
        for _ in range(300):
            a, b, c = (rand_poly(F, rng) for _ in range(3))
            assert a + b == b + a
            assert (a + b) + c == a + (b + c)
            assert a + zero == a and a + a == zero
            assert a - b == a + b and -a == a
            assert a * b == b * a
            assert (a * b) * c == a * (b * c)
            assert a * one == a and (a * zero).is_zero
            assert a * (b + c) == a * b + a * c
            if not (a.is_zero or b.is_zero):
                assert (a * b).degree == a.degree + b.degree

    def test_product_with_inner_zeros(self, gf16):
        F = gf16
        p = GFPoly(F, [1, 0, 1])      # x^2 + 1
        q = GFPoly(F, [1, 0, 0, 1])   # x^3 + 1
        assert p * q == GFPoly(F, [1, 0, 1, 1, 0, 1])

    def test_scalar_product(self, gf16, rng):
        F = gf16
        for _ in range(100):
            a = rand_poly(F, rng)
            s = rng.randrange(F.order)
            assert a.scale(s) == a * s == s * a
            assert all(a.scale(s).eval(x) == F.mul(a.eval(x), s) for x in F.elements)
        p = GFPoly(F, [3, 1])
        assert p.scale(0).is_zero
        assert p.scale(1) == p
        with pytest.raises(ValueError):
            p.scale(16)

    @pytest.mark.parametrize("other", [1, "x", 1.5, None])
    def test_addition_with_other_types(self, gf16, other):
        p = GFPoly(gf16, [1, 2])
        with pytest.raises(TypeError):
            p + other
        with pytest.raises(TypeError):
            p - other

    @pytest.mark.parametrize("other", [True, "x", 1.5, None])
    def test_product_with_other_types(self, gf16, other):
        with pytest.raises(TypeError):
            GFPoly(gf16, [1, 2]) * other

    def test_mixing_fields(self, gf16):
        G = GaloisField(4, 0b1001)
        p, q = GFPoly(gf16, [1, 2]), GFPoly(G, [1, 2])
        assert p != q
        for op in (lambda: p + q, lambda: p * q, lambda: divmod(p, q)):
            with pytest.raises(ValueError, match="different fields"):
                op()

    def test_power(self, gf16):
        p = GFPoly(gf16, [1, 3, 0, 7])
        assert p ** 0 == GFPoly.one(gf16)
        assert p ** 1 == p
        assert p ** 3 == p * p * p
        assert p ** 5 == p * p * p * p * p
        assert all((p ** 2).coeff(d) == 0 for d in range(1, (p ** 2).degree + 1, 2))
        assert GFPoly.zero(gf16) ** 0 == GFPoly.one(gf16)
        assert (GFPoly.zero(gf16) ** 3).is_zero

    @pytest.mark.parametrize("exponent", [-1, 2.0, True])
    def test_power_rejects_bad_exponent(self, gf16, exponent):
        with pytest.raises(ValueError):
            GFPoly(gf16, [1, 1]) ** exponent


class TestPolyDivision:
    """!
    @brief Long division and the //, % operators.
    """

    def test_division_identity(self, gf16, rng):
        F = gf16
        for _ in range(400):
            a = rand_poly(F, rng, 8)
            b = rand_poly(F, rng, 5)
            if b.is_zero:
                continue
            q, r = divmod(a, b)
            assert q * b + r == a
            assert r.degree < b.degree
            assert a // b == q and a % b == r
            assert a.divmod(b) == (q, r)
            assert (a * b) // b == a and ((a * b) % b).is_zero

    def test_skips_zero_steps(self, gf16):
        F = gf16
        q, r = divmod(GFPoly(F, [1, 0, 0, 0]), GFPoly(F, [1, 0, 1]))
        assert q == GFPoly(F, [1, 0]) and r == GFPoly(F, [1, 0])

    def test_edge_cases(self, gf16):
        F = gf16
        a = GFPoly(F, [1, 2, 3])
        q, r = divmod(a, GFPoly.monomial(F, 3))
        assert q.is_zero and r == a
        assert a // a == GFPoly.one(F) and (a % a).is_zero
        assert (GFPoly.zero(F) // a).is_zero
        assert a // GFPoly(F, [5]) == a.scale(F.inverse(5))

    def test_errors(self, gf16):
        a = GFPoly(gf16, [1, 2, 3])
        with pytest.raises(ZeroDivisionError):
            a // GFPoly.zero(gf16)
        with pytest.raises(ZeroDivisionError):
            a % GFPoly.zero(gf16)
        with pytest.raises(TypeError):
            a.divmod(3)


class TestPolyDerivative:
    """!
    @brief Formal derivative in characteristic 2.
    """

    def test_monomials(self, gf16):
        F = gf16
        for k in range(1, 12):
            expected = GFPoly.monomial(F, k - 1, 7) if k % 2 else GFPoly.zero(F)
            assert GFPoly.monomial(F, k, 7).derivative() == expected

    def test_constants(self, gf16):
        assert GFPoly(gf16, [9]).derivative().is_zero
        assert GFPoly.zero(gf16).derivative().is_zero

    def test_example(self, gf16):
        p = GFPoly(gf16, [7, 3, 0, 5, 1])
        assert p.derivative() == GFPoly(gf16, [3, 0, 5])

    def test_properties(self, gf16, rng):
        F = gf16
        x = GFPoly.monomial(F, 1)
        for _ in range(300):
            a, b = rand_poly(F, rng, 8), rand_poly(F, rng, 8)
            s = rng.randrange(F.order)
            assert (a + b).derivative() == a.derivative() + b.derivative()
            assert a.scale(s).derivative() == a.derivative().scale(s)
            assert (a * b).derivative() == a.derivative() * b + a * b.derivative()
            assert a.derivative().derivative().is_zero
            assert x * a.derivative() == odd_part(a)

    def test_repeated_roots(self, gf16, rng):
        F = gf16
        for _ in range(100):
            roots = rng.sample(range(F.order), rng.randint(1, 5))
            p = GFPoly.from_roots(F, roots)
            assert all(p.derivative().eval(r) != 0 for r in roots)
            double = p * GFPoly(F, [1, roots[0]])
            assert double.derivative().eval(roots[0]) == 0


class TestPolyEvaluation:
    """!
    @brief Horner evaluation and root search.
    """

    def test_horner_matches_definition(self, gf16, rng):
        F = gf16
        for _ in range(100):
            a, b = rand_poly(F, rng), rand_poly(F, rng)
            for x in F.elements:
                direct = 0
                for i, c in enumerate(reversed(a.coeffs)):
                    direct ^= F.mul(c, F.power(x, i))
                assert a.eval(x) == a(x) == direct
                assert (a + b)(x) == a(x) ^ b(x)
                assert (a * b)(x) == F.mul(a(x), b(x))

    def test_zero_polynomial(self, gf16):
        z = GFPoly.zero(gf16)
        assert all(z(x) == 0 for x in gf16.elements)
        assert z.roots() == list(gf16.elements)

    def test_every_element_is_a_root_of_x_q_plus_x(self, gf16):
        p = GFPoly(gf16, [1] + [0] * 14 + [1, 0])
        assert p.roots() == list(gf16.elements)

    def test_no_roots(self, gf16):
        assert GFPoly(gf16, [1, 1, 8]).roots() == []    # irreducible quadratic
        assert GFPoly(gf16, [5]).roots() == []

    def test_rejects_points_outside_the_field(self, gf16):
        with pytest.raises(ValueError):
            GFPoly(gf16, [1, 1]).eval(16)


class TestPolyRepresentation:
    """!
    @brief str, to_alpha_str, repr, equality and hashing.
    """

    def test_str(self, gf16):
        F = gf16
        assert str(GFPoly.zero(F)) == "0"
        assert str(GFPoly.one(F)) == "1"
        assert str(GFPoly(F, [7])) == "7"
        assert str(GFPoly(F, [1, 0])) == "x"
        assert str(GFPoly(F, [3, 0, 1, 5])) == "3·x^3 + x + 5"

    def test_to_alpha_str(self, gf16):
        F = gf16
        assert GFPoly(F, [2, 3, 0, 1]).to_alpha_str() == "α·x^3 + α^4·x^2 + 1"
        assert GFPoly(F, [2]).to_alpha_str() == "α"
        assert GFPoly.zero(F).to_alpha_str() == "0"

    def test_repr_round_trip(self, gf16):
        p = GFPoly(gf16, [3, 0, 5, 1])
        assert repr(p) == "GFPoly(GaloisField(m=4, poly=0x3), [3, 0, 5, 1])"
        assert eval(repr(p), {"GFPoly": GFPoly, "GaloisField": GaloisField}) == p

    def test_equality_and_hash(self, gf16):
        assert GFPoly(gf16, [1, 2]) == GFPoly(gf16, [0, 1, 2])
        assert GFPoly(gf16, [1, 2]) != GFPoly(gf16, [1, 3])
        assert GFPoly(gf16, [1]) != 1
        assert len({GFPoly(gf16, [1, 2]), GFPoly(gf16, [0, 1, 2]), GFPoly(gf16, [1, 3])}) == 2


class TestReedSolomon:
    """!
    @brief Integration checks against the TP3 reference values.
    """

    @staticmethod
    def encode(F, g, msg):
        shifted = GFPoly(F, msg) * GFPoly.monomial(F, g.degree)
        return shifted + shifted % g

    def test_qr_reference_vector(self, gf256):
        F = gf256
        g = GFPoly.from_roots(F, [F.exp(i) for i in range(16)])
        assert g.degree == 16 and g.lead == 1
        msg = list(range(28))
        cw = self.encode(F, g, msg)
        word = cw.to_list(44)
        assert word[:28] == msg
        assert word[28:] == QR_PARITY
        assert all(cw(F.exp(j)) == 0 for j in range(16))

    def test_qr_bus_layout(self, gf256):
        F = gf256
        g = GFPoly.from_roots(F, [F.exp(i) for i in range(16)])
        msg = list(range(28))
        i_msg = GFPoly(F, msg).to_int(28)
        o_cw = self.encode(F, g, msg).to_int(44)
        assert i_msg >> 216 == 0x00 and i_msg & 0xFF == 0x1B
        assert o_cw >> 128 == i_msg
        assert o_cw & ((1 << 128) - 1) == GFPoly(F, QR_PARITY).to_int(16)

    def test_reduced_code(self, gf16):
        F = gf16
        g = GFPoly.from_roots(F, [F.exp(i) for i in range(6)])
        assert list(g.coeffs) == [1, 0xA, 0xF, 2, 4, 3, 1]
        msg = GFPoly.from_int(F, 0x52D09F17C, 9)
        cw = self.encode(F, g, msg.to_list(9))
        assert cw.to_list(15)[:9] == msg.to_list(9)
        assert all(cw(F.exp(j)) == 0 for j in range(6))

    def test_single_error_correction(self, gf16, rng):
        F = gf16
        g = GFPoly.from_roots(F, [F.exp(i) for i in range(1, 5)])
        for _ in range(50):
            cw = self.encode(F, g, [rng.randrange(16) for _ in range(11)])
            pos, val = rng.randrange(15), rng.randrange(1, 16)
            rx = cw + GFPoly.monomial(F, pos, val)
            s1, s2 = rx(F.exp(1)), rx(F.exp(2))
            X = F.div(s2, s1)
            assert F.log(X) == pos and F.div(s1, X) == val


def main() -> int:
    """!
    @brief Run this file through pytest, with coverage when pytest-cov is available.
    @return pytest exit code.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    cmd = [sys.executable, "-m", "pytest", os.path.basename(__file__), "-q"]
    try:
        import pytest_cov  # noqa: F401
        cmd += ["--cov=galois", "--cov-branch", "--cov-report=term-missing",
                "--cov-fail-under=100"]
    except ImportError:
        print("pytest-cov is not installed; running without coverage")
    # a subprocess, so galois is imported after coverage starts
    return subprocess.call(cmd + sys.argv[1:], cwd=here)


if __name__ == "__main__":
    sys.exit(main())