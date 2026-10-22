#!/usr/bin/env python3
"""!
@file test_galois.py
@brief Functional test of galois.py (GaloisField and GFPoly).

Checks the algebraic properties of both classes and prints the field tables
along the way. The exit code is 0 when every check passes and 1 otherwise.

Usage:
    python3 test_galois.py              # everything
    python3 test_galois.py --quick      # skip the heaviest sweeps
    python3 test_galois.py --no-tables  # do not print tables
"""

import argparse
import random
import sys
import time
from math import gcd

from galois import GFPoly, GaloisField

WIDTH = 78


class Runner:
    """!
    @brief Check counter that collects failures without stopping the run.
    """

    def __init__(self) -> None:
        self.passed = 0
        self.failures = []
        self.current = "?"
        self.mark = 0

    def section(self, title: str) -> None:
        """!
        @brief Start a numbered section of the report.
        """
        self.current = title
        self.mark = len(self.failures)
        print()
        print("=" * WIDTH)
        print(f"  {title}")
        print("=" * WIDTH)

    def sub(self, title: str) -> None:
        """!
        @brief Start a block inside the current section.
        """
        self.mark = len(self.failures)
        print(f"\n--- {title}")

    def check(self, desc: str, cond: bool, detail: str = "") -> bool:
        """!
        @brief Record one check.
        @param desc Short description.
        @param cond Result of the check.
        @param detail Extra text printed on failure.
        @return cond as a bool.
        """
        if cond:
            self.passed += 1
        else:
            self.failures.append((self.current, desc, detail))
            print(f"  [FAIL] {desc}" + (f"\n         {detail}" if detail else ""))
        return bool(cond)

    def raises(self, desc: str, exc, fn, *args, **kwargs) -> bool:
        """!
        @brief Check that fn(*args, **kwargs) raises exc.
        """
        try:
            result = fn(*args, **kwargs)
        except exc:
            self.passed += 1
            return True
        except Exception as e:  # noqa: BLE001
            return self.check(desc, False, f"raised {type(e).__name__} instead of {exc.__name__}")
        return self.check(desc, False, f"did not raise, returned {result!r}")

    def error(self, part: str, e: BaseException) -> None:
        """!
        @brief Record an unexpected exception and keep going.
        """
        self.failures.append((part, f"exception {type(e).__name__}", str(e)))
        print(f"\n  [ERROR] {part} stopped with {type(e).__name__}: {e}")
        print("          continuing with the next sections")

    def summary(self, seconds: float) -> int:
        """!
        @brief Print the final report.
        @return Process exit code.
        """
        print()
        print("=" * WIDTH)
        if not self.failures:
            print(f"  ALL OK - {self.passed} checks in {seconds:.2f} s")
            print("=" * WIDTH)
            return 0
        print(f"  {len(self.failures)} FAILURES out of {self.passed + len(self.failures)} checks")
        for sec, desc, detail in self.failures:
            print(f"    - [{sec}] {desc} {detail}")
        print("=" * WIDTH)
        return 1


T = Runner()


def done(desc: str) -> None:
    """!
    @brief Close a block; prints ok only if nothing failed since the last sub().
    """
    if len(T.failures) == T.mark:
        print(f"  ok  {desc}")
    else:
        print(f"  ->  {desc}: NO (see failures above)")


def hexs(values, width: int = 2) -> str:
    """!
    @brief Format a sequence of field elements as hex digits.
    """
    return " ".join(f"{v:0{width}X}" for v in values)


def field_table(F: GaloisField) -> None:
    """!
    @brief Print every field element as power of alpha, binary and polynomial.
    """
    print(f"\n  Elements of {F}")
    print(f"  {'i':>3} | {'alpha^i':>7} | {'bin':>{F.m}} | {'hex':>4} | polynomial")
    print(f"  {'-'*3}-+-{'-'*7}-+-{'-'*F.m}-+-{'-'*4}-+-{'-'*28}")
    print(f"  {'-':>3} | {0:>7} | {0:0{F.m}b} | {0:>4} | 0")
    for i in range(F.n):
        a = F.exp(i)
        print(f"  {i:>3} | {a:>7} | {a:0{F.m}b} | {a:#04x} | {F.poly_str(a)}")


def op_table(F: GaloisField, op, name: str, skip_zero: bool = False) -> None:
    """!
    @brief Print the table of a binary operation in hex digits.
    """
    elems = [e for e in F.elements if not (skip_zero and e == 0)]
    print(f"\n  {name} table in GF(2^{F.m})  (hex digits)")
    header = "     | " + " ".join(f"{e:X}" for e in elems)
    print(header)
    print("   " + "-" * (len(header) - 3))
    for a in elems:
        row = " ".join(f"{op(a, b):X}" for b in elems)
        print(f"   {a:X} | {row}")


def inverse_table(F: GaloisField) -> None:
    """!
    @brief Print the multiplicative inverse of every non-zero element.
    """
    print(f"\n  Multiplicative inverses in GF(2^{F.m})")
    print("     a  | " + " ".join(f"{a:>2}" for a in range(1, F.order)))
    print("    a^-1| " + " ".join(f"{F.inverse(a):>2}" for a in range(1, F.order)))


def mult_order(F: GaloisField, a: int) -> int:
    """!
    @brief Multiplicative order of a non-zero element.
    """
    k, acc = 1, a
    while acc != 1:
        acc = F.mul(acc, a)
        k += 1
    return k


def totient(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def random_poly(F: GaloisField, rng: random.Random, max_degree: int = 6) -> GFPoly:
    """!
    @brief Random polynomial; sometimes the zero polynomial.
    """
    d = rng.randrange(-1, max_degree + 1)
    if d < 0:
        return GFPoly.zero(F)
    return GFPoly(F, [rng.randrange(1, F.order)] + [rng.randrange(F.order) for _ in range(d)])


def mat_vec(M, a: int) -> int:
    """!
    @brief Apply a binary matrix to the bits of a.
    @return Packed result; bit i is row i dotted with a over GF(2).
    """
    out = 0
    for i, row in enumerate(M):
        bit = 0
        for j, mij in enumerate(row):
            bit ^= mij & (a >> j) & 1
        out |= bit << i
    return out


def mask_vec(masks, a: int) -> int:
    """!
    @brief Bit i of the result is the parity of a & masks[i].
    """
    return sum((bin(a & k).count("1") & 1) << i for i, k in enumerate(masks))


def mat_add(A, B):
    return tuple(tuple(x ^ y for x, y in zip(ra, rb)) for ra, rb in zip(A, B))


def mat_mul(A, B):
    n = len(A)
    return tuple(
        tuple(sum(A[i][k] & B[k][j] for k in range(n)) & 1 for j in range(n))
        for i in range(n)
    )


def xor_count(masks) -> int:
    """!
    @brief Two-input XOR gates of a constant multiplier, without sharing.
    """
    return sum(max(bin(k).count("1") - 1, 0) for k in masks)


def part1_construction(args) -> None:
    """!
    @brief Field construction and validation of P(x).
    """
    T.section("1. Field construction and validation of P(x)")

    T.sub("Known primitive polynomials")
    print(f"  {'m':>2} | {'P(x) (hex)':>10} | {'2^m':>6} | polynomial")
    print(f"  {'-'*2}-+-{'-'*10}-+-{'-'*6}-+-{'-'*40}")
    for m, poly in [(1, 0b1), (2, 0b11), (3, 0b011), (4, 0b0011), (4, 0b1001),
                    (5, 0b00101), (6, 0b000011), (8, 0x1D), (10, 0x009), (12, 0x053)]:
        F = GaloisField(m, poly)
        T.check(f"GF(2^{m}) poly={poly:#x} is primitive", F.is_primitive)
        T.check(f"GF(2^{m}) size", len(F) == 2 ** m)
        print(f"  {m:>2} | {poly:>#10x} | {2**m:>6} | {F.primitive_poly_str()}")

    T.sub("Invalid inputs are rejected")
    T.raises("m = 0", ValueError, GaloisField, 0, 1)
    T.raises("negative m", ValueError, GaloisField, -3, 1)
    T.raises("non-integer m", ValueError, GaloisField, 2.5, 1)
    T.raises("non-integer poly", TypeError, GaloisField, 4, "0011")
    T.raises("poly out of range (>= 2^m)", ValueError, GaloisField, 4, 0b10011)
    T.raises("poly without constant term", ValueError, GaloisField, 4, 0b0010)
    done("invalid m and P(x) are rejected")

    T.sub("Irreducible but not primitive: P(x) = x^4 + x^3 + x^2 + x + 1")
    T.raises("rejected by default", ValueError, GaloisField, 4, 0b1111)
    G = GaloisField(4, 0b1111, require_primitive=False)
    T.check("accepted with require_primitive=False", not G.is_primitive)
    T.check("another generator is found", mult_order(G, G.generator) == G.n)
    T.check("still a field (every a != 0 is invertible)",
            all(G.mul(a, G.inverse(a)) == 1 for a in range(1, G.order)))
    print(f"  x has order {mult_order(G, 2)}, expected {G.n}")
    print(f"  chosen generator: {G.generator}, order {mult_order(G, G.generator)}")
    done("a non-primitive P(x) is detected and handled")

    T.sub("AES field: P(x) = x^8 + x^4 + x^3 + x + 1 (irreducible, not primitive)")
    T.raises("rejected by default", ValueError, GaloisField, 8, 0x1B)
    AES = GaloisField(8, 0x1B, require_primitive=False)
    T.check("0x57 * 0x83 = 0xC1 (FIPS-197 example)", AES.mul(0x57, 0x83) == 0xC1,
            f"got {AES.mul(0x57, 0x83):#x}")
    T.check("matches the reference implementation", AES.mul_ref(0x57, 0x83) == 0xC1)
    print(f"  0x57 * 0x83 = {AES.mul(0x57, 0x83):#04x}   (order of x: {mult_order(AES, 2)})")
    done("the AES field works with require_primitive=False")

    T.sub("Reducible polynomials do not define a field")
    for m, bad, why in [(4, 0b0101, "x^4 + x^2 + 1 = (x^2 + x + 1)^2"),
                        (4, 0b0111, "x^4 + x^2 + x + 1 = (x + 1)(x^3 + x^2 + 1)"),
                        (6, 0b000101, "x^6 + x^2 + 1 = (x^3 + x + 1)^2")]:
        T.raises(f"GF(2^{m}) poly={bad:#x} ({why})", ValueError,
                 GaloisField, m, bad, require_primitive=False)
    done("reducible polynomials are rejected")

    T.sub("Field equality and hash")
    T.check("same m and poly => equal", GaloisField(4, 3) == GaloisField(4, 3))
    T.check("different poly => not equal", GaloisField(4, 3) != GaloisField(4, 9))
    T.check("consistent hash", len({GaloisField(4, 3), GaloisField(4, 3), GaloisField(4, 9)}) == 2)
    done("fields compare by value")


def part2_tables(F: GaloisField, args) -> None:
    """!
    @brief Field, addition, product, division and inverse tables of GF(2^4).
    """
    T.section("2. Tables of GF(2^4) with P(x) = x^4 + x + 1")

    if args.no_tables:
        print("  (skipped with --no-tables)")
        return

    field_table(F)
    op_table(F, F.add, "Addition (XOR)")
    op_table(F, F.mul, "Product")
    op_table(F, F.div, "Division", skip_zero=True)
    inverse_table(F)

    print("\n  Multiplicative order of each element")
    print("     a  | " + " ".join(f"{a:>2}" for a in range(1, F.order)))
    print("   ord  | " + " ".join(f"{mult_order(F, a):>2}" for a in range(1, F.order)))

    prims = [a for a in range(1, F.order) if mult_order(F, a) == F.n]
    print(f"\n  Primitive elements: {prims}")
    print(f"  Count: {len(prims)}, expected phi({F.n}) = {totient(F.n)}")
    T.check("number of primitive elements = phi(2^m - 1)", len(prims) == totient(F.n))


def part3_axioms(args) -> None:
    """!
    @brief Exhaustive check of the field axioms in small fields.
    """
    T.section("3. Field axioms")

    fields = [GaloisField(2, 0b11), GaloisField(3, 0b011),
              GaloisField(4, 0b0011), GaloisField(4, 0b1001)]
    if not args.quick:
        fields.append(GaloisField(5, 0b00101))

    for F in fields:
        T.sub(f"Exhaustive sweep in {F}  ({F.order}^3 = {F.order**3} triples)")

        closure = add_assoc = add_comm = add_zero = add_inv = True
        mul_assoc = mul_comm = mul_one = distrib = True

        for a in F.elements:
            closure &= F.add(a, a) in F
            add_zero &= F.add(a, 0) == a
            add_inv &= F.add(a, a) == 0    # characteristic 2
            mul_one &= F.mul(a, 1) == a
            closure &= F.mul(a, a) in F
            for b in F.elements:
                add_comm &= F.add(a, b) == F.add(b, a)
                mul_comm &= F.mul(a, b) == F.mul(b, a)
                for c in F.elements:
                    add_assoc &= F.add(F.add(a, b), c) == F.add(a, F.add(b, c))
                    mul_assoc &= F.mul(F.mul(a, b), c) == F.mul(a, F.mul(b, c))
                    distrib &= F.mul(a, F.add(b, c)) == F.add(F.mul(a, b), F.mul(a, c))

        T.check("closure", closure)
        T.check("addition is associative", add_assoc)
        T.check("addition is commutative", add_comm)
        T.check("additive identity (0)", add_zero)
        T.check("additive inverse: a + a = 0", add_inv)
        T.check("product is associative", mul_assoc)
        T.check("product is commutative", mul_comm)
        T.check("multiplicative identity (1)", mul_one)
        T.check("distributivity", distrib)

        nz = range(1, F.order)
        T.check("every a != 0 has an inverse", all(F.mul(a, F.inverse(a)) == 1 for a in nz))
        T.check("a * 0 = 0", all(F.mul(a, 0) == 0 for a in F.elements))
        T.check("no zero divisors", all(F.mul(a, b) != 0 for a in nz for b in nz))
        T.check("mul matches mul_ref",
                all(F.mul(a, b) == F.mul_ref(a, b) for a in F.elements for b in F.elements))
        T.check("div(a, b) = a * b^-1",
                all(F.div(a, b) == F.mul(a, F.inverse(b)) for a in F.elements for b in nz))
        T.check("(a / b) * b = a",
                all(F.mul(F.div(a, b), b) == a for a in F.elements for b in nz))
        T.check("a / 1 = a", all(F.div(a, 1) == a for a in F.elements))
        T.check("a / a = 1", all(F.div(a, a) == 1 for a in nz))
        T.check("1 / a = a^-1", all(F.div(1, a) == F.inverse(a) for a in nz))
        T.check("inverse(inverse(a)) = a", all(F.inverse(F.inverse(a)) == a for a in nz))
        done(f"{F.order**3} triples checked")


def part4_structure(args) -> None:
    """!
    @brief Logarithms, powers, Frobenius map and subfields.
    """
    T.section("4. Logarithms, powers and field structure")

    fields = [GaloisField(4, 0b0011), GaloisField(8, 0x1D)]
    if not args.quick:
        fields.append(GaloisField(10, 0x009))

    for F in fields:
        T.sub(f"{F}")

        T.check("exp(log(a)) = a", all(F.exp(F.log(a)) == a for a in range(1, F.order)))
        T.check("log(exp(i)) = i", all(F.log(F.exp(i)) == i for i in range(F.n)))
        T.check("exp covers the whole group", len({F.exp(i) for i in range(F.n)}) == F.n)
        T.check("alpha^(2^m - 1) = 1", F.power(F.alpha, F.n) == 1)
        T.check("Lagrange: a^(2^m - 1) = 1 for every a != 0",
                all(F.power(a, F.n) == 1 for a in range(1, F.order)))
        T.check("a^(2^m) = a for every a", all(F.power(a, F.order) == a for a in F.elements))

        rng = random.Random(1234 + F.m)
        pow_ok = True
        for _ in range(200):
            a = rng.randrange(F.order)
            k = rng.randrange(0, 40)
            acc = 1
            for _ in range(k):
                acc = F.mul(acc, a)
            pow_ok &= F.power(a, k) == acc
        T.check("power(a, k) = repeated product", pow_ok)
        T.check("negative exponents",
                all(F.power(a, -3) == F.power(F.inverse(a), 3) for a in range(1, F.order)))
        T.check("power(a, 0) = 1", all(F.power(a, 0) == 1 for a in F.elements))

        step = max(1, F.order // 64)
        sample_a = range(0, F.order, step)
        sample_b = range(1, F.order, step)
        T.check("div(a, b) = a * b^-1",
                all(F.div(a, b) == F.mul(a, F.inverse(b)) for a in sample_a for b in sample_b))
        T.check("(a / b) * b = a",
                all(F.mul(F.div(a, b), b) == a for a in sample_a for b in sample_b))
        T.check("div(a, b) = a * power(b, -1)",
                all(F.div(a, b) == F.mul(a, F.power(b, -1)) for a in sample_a for b in sample_b))

        # squaring is an automorphism in characteristic 2
        T.check("Frobenius: (a + b)^2 = a^2 + b^2",
                all(F.power(a ^ b, 2) == F.power(a, 2) ^ F.power(b, 2)
                    for a in F.elements for b in sample_a))

        # {a : a^(2^d) = a} is the subfield GF(2^d) when d divides m
        for d in range(1, F.m + 1):
            if F.m % d:
                continue
            sub = [a for a in F.elements if F.power(a, 2 ** d) == a]
            T.check(f"subfield GF(2^{d}) has 2^{d} elements", len(sub) == 2 ** d,
                    f"found {len(sub)}")
        done("multiplicative structure is consistent")

    F = GaloisField(4, 0b0011)
    if not args.no_tables:
        print("\n  Subfields of GF(2^4)")
        for d in [1, 2, 4]:
            print(f"    GF(2^{d}) = {[a for a in F.elements if F.power(a, 2 ** d) == a]}")


def part5_degenerate(F: GaloisField, args) -> None:
    """!
    @brief Behaviour with 0, invalid elements and GF(2).
    """
    T.section("5. Degenerate cases")

    T.sub("Element 0")
    T.raises("inverse(0) raises ZeroDivisionError", ZeroDivisionError, F.inverse, 0)
    T.raises("div(a, 0) raises ZeroDivisionError", ZeroDivisionError, F.div, 7, 0)
    T.raises("log(0) raises ZeroDivisionError", ZeroDivisionError, F.log, 0)
    T.raises("0^-1 raises ZeroDivisionError", ZeroDivisionError, F.power, 0, -1)
    T.check("power(0, 0) = 1 (empty product)", F.power(0, 0) == 1)
    T.check("power(0, n > 0) = 0", all(F.power(0, k) == 0 for k in range(1, 10)))
    T.check("div(0, b) = 0 for b != 0", all(F.div(0, b) == 0 for b in range(1, F.order)))
    done("degenerate cases raise instead of returning garbage")

    T.sub("Field membership")
    T.raises("element above 2^m - 1", ValueError, F.check, F.order)
    T.raises("negative element", ValueError, F.check, -1)
    T.raises("float", TypeError, F.check, 3.0)
    T.raises("string", TypeError, F.check, "3")
    T.raises("bool", TypeError, F.check, True)
    T.check("__contains__ accepts valid elements", all(a in F for a in F.elements))
    T.check("__contains__ rejects invalid ones", (F.order not in F) and (-1 not in F))
    done("check() filters anything that is not a field element")

    T.sub("Edge case m = 1: GF(2)")
    B = GaloisField(1, 0b1)
    T.check("GF(2) has 2 elements", len(B) == 2)
    T.check("addition in GF(2) is XOR", [B.add(a, b) for a in (0, 1) for b in (0, 1)] == [0, 1, 1, 0])
    T.check("product in GF(2) is AND", [B.mul(a, b) for a in (0, 1) for b in (0, 1)] == [0, 0, 0, 1])
    T.check("inverse(1) = 1", B.inverse(1) == 1)
    done("GF(2) behaves as Boolean algebra")


def part6_poly_construction(F: GaloisField, args) -> None:
    """!
    @brief GFPoly construction, normalization and validation.
    """
    T.section("6. GFPoly: construction, normalization and validation")

    T.sub("Normalization (leading zeros)")
    cases = [
        ([0, 0, 1, 5], (1, 5), 1),
        ([1, 0, 0], (1, 0, 0), 2),
        ([0, 0, 0], (0,), -1),
        ([], (0,), -1),
        ([0], (0,), -1),
        ([7], (7,), 0),
    ]
    print(f"  {'input':>14} | {'stored':>12} | {'degree':>6} | printed")
    print(f"  {'-'*14}-+-{'-'*12}-+-{'-'*6}-+-{'-'*24}")
    for given, expected, degree in cases:
        p = GFPoly(F, given)
        T.check(f"normalizes {given}", p.coeffs == expected, f"got {p.coeffs}")
        T.check(f"degree of {given}", p.degree == degree, f"got {p.degree}")
        print(f"  {str(given):>14} | {str(list(p.coeffs)):>12} | {p.degree:>6} | {p}")
    T.check("only the zero polynomial keeps a leading 0",
            GFPoly.zero(F).coeffs == (0,) and GFPoly.zero(F).is_zero)
    done("no leading zeros except in the zero polynomial")

    T.sub("Coefficient validation")
    T.raises("coefficient out of the field", ValueError, GFPoly, F, [1, F.order])
    T.raises("negative coefficient", ValueError, GFPoly, F, [1, -1])
    T.raises("float coefficient", TypeError, GFPoly, F, [1, 2.0])
    T.raises("string coefficient", TypeError, GFPoly, F, [1, "a"])
    T.raises("field is not a GaloisField", TypeError, GFPoly, 16, [1, 2])
    T.raises("scalar coeffs", TypeError, GFPoly, F, 5)
    try:
        GFPoly(F, [1, 0, 99])
    except ValueError as e:
        print(f"  error message: {e}")
    done("invalid coefficients are rejected with context")

    T.sub("Alternative constructors")
    T.check("zero", GFPoly.zero(F).is_zero)
    T.check("one", GFPoly.one(F).coeffs == (1,) and GFPoly.one(F).degree == 0)
    T.check("monomial(3)", GFPoly.monomial(F, 3).coeffs == (1, 0, 0, 0))
    T.check("monomial(2, coeff=7)", GFPoly.monomial(F, 2, 7).coeffs == (7, 0, 0))
    T.raises("monomial with negative degree", ValueError, GFPoly.monomial, F, -1)
    done("zero / one / monomial")

    T.sub("Coefficient access")
    p = GFPoly(F, [3, 0, 5, 1])    # 3x^3 + 5x + 1
    T.check("coeff(3)", p.coeff(3) == 3)
    T.check("coeff(2)", p.coeff(2) == 0)
    T.check("coeff(1)", p.coeff(1) == 5)
    T.check("coeff(0)", p.coeff(0) == 1)
    T.check("coeff out of range", p.coeff(9) == 0 and p.coeff(-1) == 0)
    T.check("lead", p.lead == 3)
    print(f"  p = {p}   degree {p.degree}, lead {p.lead}")
    done("coeff() indexes by degree")

    T.sub("Equality, hash and repr")
    T.check("equal", GFPoly(F, [1, 2]) == GFPoly(F, [0, 1, 2]))
    T.check("not equal", GFPoly(F, [1, 2]) != GFPoly(F, [1, 3]))
    T.check("consistent hash",
            len({GFPoly(F, [1, 2]), GFPoly(F, [0, 1, 2]), GFPoly(F, [1, 3])}) == 2)
    T.check("repr round trip", eval(repr(p), {"GFPoly": GFPoly, "GaloisField": GaloisField}) == p)
    G = GaloisField(4, 0b1001)
    T.raises("operations across different fields", ValueError,
             lambda: GFPoly(F, [1, 2]) + GFPoly(G, [1, 2]))
    done("polynomials compare by value and repr can rebuild them")


def part7_poly_ring(F: GaloisField, args) -> None:
    """!
    @brief Ring axioms, powers and scaling of GFPoly.
    """
    T.section("7. GFPoly: ring structure")

    rng = random.Random(20260805)
    n = 100 if args.quick else 400
    T.sub(f"{n} random triples of polynomials over {F}")

    zero, one = GFPoly.zero(F), GFPoly.one(F)
    add_comm = add_assoc = add_zero = add_inv = True
    mul_comm = mul_assoc = mul_one = mul_zero = distrib = degree = True

    for _ in range(n):
        a, b, c = (random_poly(F, rng) for _ in range(3))
        add_comm &= (a + b) == (b + a)
        add_assoc &= ((a + b) + c) == (a + (b + c))
        add_zero &= (a + zero) == a
        add_inv &= (a + a) == zero
        mul_comm &= (a * b) == (b * a)
        mul_assoc &= ((a * b) * c) == (a * (b * c))
        mul_one &= (a * one) == a
        mul_zero &= (a * zero) == zero
        distrib &= (a * (b + c)) == (a * b + a * c)
        if not a.is_zero and not b.is_zero:
            degree &= (a * b).degree == a.degree + b.degree
        else:
            degree &= (a * b).is_zero

    T.check("addition is commutative", add_comm)
    T.check("addition is associative", add_assoc)
    T.check("additive identity", add_zero)
    T.check("a + a = 0 (characteristic 2)", add_inv)
    T.check("product is commutative", mul_comm)
    T.check("product is associative", mul_assoc)
    T.check("multiplicative identity", mul_one)
    T.check("a * 0 = 0", mul_zero)
    T.check("distributivity", distrib)
    T.check("deg(a * b) = deg(a) + deg(b)", degree)
    T.check("subtraction equals addition",
            all((a - b) == (a + b)
                for a, b in ((random_poly(F, rng), random_poly(F, rng)) for _ in range(50))))
    done(f"{n} triples checked")

    T.sub("Polynomial powers")
    p = GFPoly(F, [1, 3, 0, 7])
    T.check("p^0 = 1", p ** 0 == one)
    T.check("p^1 = p", p ** 1 == p)
    T.check("p^3 = p * p * p", p ** 3 == p * p * p)
    T.check("p^5 = p * p * p * p * p", p ** 5 == p * p * p * p * p)
    T.check("p(x)^2 only has even-degree terms",
            all((p ** 2).coeff(d) == 0 for d in range(1, (p ** 2).degree + 1, 2)))
    T.raises("negative exponent", ValueError, lambda: p ** -1)
    print(f"  p   = {p}")
    print(f"  p^2 = {p**2}")
    done("power by repeated squaring")

    T.sub("Scaling")
    scale_ok = True
    for _ in range(200):
        a = random_poly(F, rng)
        s = rng.randrange(F.order)
        sa = a.scale(s)
        scale_ok &= sa == a * s == s * a
        scale_ok &= all(sa.eval(x) == F.mul(a.eval(x), s) for x in F.elements)
    T.check("scale(s) matches multiplying by the scalar", scale_ok)
    T.check("scale(0) gives zero", random_poly(F, rng).scale(0).is_zero)
    q = random_poly(F, rng)
    T.check("scale(1) is the identity", q.scale(1) == q)
    T.raises("scalar out of the field", ValueError, p.scale, F.order)
    done("scaling is consistent with the product")


def part8_division(F: GaloisField, args) -> None:
    """!
    @brief Polynomial long division.
    """
    T.section("8. GFPoly: long division")

    rng = random.Random(31415)
    n = 150 if args.quick else 600
    T.sub(f"Division identity on {n} pairs: a = q * b + r, deg(r) < deg(b)")

    identity = rem_degree = operators = exact = True
    for _ in range(n):
        a = random_poly(F, rng, 8)
        b = random_poly(F, rng, 5)
        if b.is_zero:
            continue
        q, r = divmod(a, b)
        identity &= (q * b + r) == a
        rem_degree &= r.degree < b.degree
        operators &= (a // b) == q and (a % b) == r
        exact &= ((a * b) // b) == a and ((a * b) % b).is_zero

    T.check("a = q * b + r", identity)
    T.check("deg(r) < deg(b)", rem_degree)
    T.check("// and % match divmod", operators)
    T.check("(a * b) / b = a with zero remainder", exact)
    done(f"{n} pairs checked")

    T.sub("Edge cases")
    a = GFPoly(F, [1, 2, 3])
    b = GFPoly(F, [1, 0, 0, 0])
    q, r = divmod(a, b)
    T.check("lower-degree dividend: zero quotient", q.is_zero)
    T.check("lower-degree dividend: remainder = dividend", r == a)
    T.check("a / a = 1", (a // a) == GFPoly.one(F) and (a % a).is_zero)
    T.check("0 / a = 0", (GFPoly.zero(F) // a).is_zero)
    T.check("constant divisor", (a // GFPoly(F, [5])) == a.scale(F.inverse(5)))
    T.raises("division by the zero polynomial", ZeroDivisionError, lambda: a // GFPoly.zero(F))
    T.raises("divisor of the wrong type", TypeError, a.divmod, 3)
    done("division edge cases")

    T.sub("Worked example")
    a = GFPoly(F, [1, 0, 5, 3, 12])
    b = GFPoly(F, [1, 7, 2])
    q, r = divmod(a, b)
    print(f"  a(x) = {a}")
    print(f"  b(x) = {b}")
    print(f"  q(x) = {q}")
    print(f"  r(x) = {r}")
    print(f"  q * b + r = {q * b + r}")
    T.check("worked example", q * b + r == a)


def part9_evaluation(F: GaloisField, args) -> None:
    """!
    @brief Horner evaluation, from_roots and roots.
    """
    T.section("9. GFPoly: evaluation and construction from roots")

    rng = random.Random(2718)

    T.sub("Horner against the direct sum of terms")
    horner = hom = True
    for _ in range(200):
        a = random_poly(F, rng)
        b = random_poly(F, rng)
        for x in F.elements:
            direct = 0
            for i, c in enumerate(reversed(a.coeffs)):
                direct ^= F.mul(c, F.power(x, i))
            horner &= a.eval(x) == direct == a(x)
            hom &= (a + b).eval(x) == (a.eval(x) ^ b.eval(x))
            hom &= (a * b).eval(x) == F.mul(a.eval(x), b.eval(x))
    T.check("eval matches the definition", horner)
    T.check("eval is a ring homomorphism", hom)
    T.check("the zero polynomial evaluates to 0", all(GFPoly.zero(F).eval(x) == 0 for x in F.elements))
    T.raises("evaluation outside the field", ValueError, GFPoly(F, [1, 1]).eval, F.order)
    done("evaluation is consistent")

    T.sub("from_roots")
    roots_ok = degree_ok = monic = factor = True
    for _ in range(200):
        k = rng.randrange(0, 6)
        roots = [rng.randrange(F.order) for _ in range(k)]
        p = GFPoly.from_roots(F, roots)
        degree_ok &= p.degree == k
        monic &= p.lead == 1
        roots_ok &= all(p.eval(r) == 0 for r in roots)
        for r in set(roots):
            factor &= (p % GFPoly(F, [1, r])).is_zero
        if len(set(roots)) == k:
            roots_ok &= sorted(p.roots()) == sorted(set(roots))
    T.check("degree = number of roots", degree_ok)
    T.check("the result is monic", monic)
    T.check("vanishes at every root", roots_ok)
    T.check("(x - r) divides the polynomial", factor)
    T.check("from_roots([]) = 1 (empty product)", GFPoly.from_roots(F, []) == GFPoly.one(F))
    T.raises("root out of the field", ValueError, GFPoly.from_roots, F, [F.order])
    done("construction from roots")

    T.sub("Example: roots alpha^0, alpha^1, alpha^2, alpha^3")
    roots = [F.exp(i) for i in range(4)]
    g = GFPoly.from_roots(F, roots)
    print(f"  roots = {roots}")
    print(f"  g(x)  = {g}")
    print(f"  g(x)  = {g.to_alpha_str()}")
    print(f"  {'x':>4} | {'g(x)':>5}")
    print(f"  {'-'*4}-+-{'-'*5}")
    for x in F.elements:
        tag = "  <- root" if g.eval(x) == 0 else ""
        print(f"  {x:>4} | {g.eval(x):>5}{tag}")
    T.check("g vanishes exactly at its 4 roots", sorted(g.roots()) == sorted(roots))

    T.sub("Every element is a root of x^(2^m) - x")
    p = GFPoly(F, [1] + [0] * (F.order - 2) + [1, 0])    # x^(2^m) + x
    T.check("x^16 + x vanishes on all of GF(2^4)", all(p.eval(x) == 0 for x in F.elements))
    T.check("x^16 + x has 16 roots", len(p.roots()) == F.order)
    done("consistent with finite field theory")


def part10_reed_solomon(F: GaloisField, args) -> None:
    """!
    @brief Integration demo: RS(15, 11) over GF(2^4) with single-error correction.
    """
    T.section("10. Integration demo: Reed-Solomon (15, 11) over GF(2^4)")

    n, k = 15, 11
    npar = n - k    # 4 parity symbols, t = 2

    roots = [F.exp(i) for i in range(1, npar + 1)]
    g = GFPoly.from_roots(F, roots)
    print(f"  n = {n}, k = {k}, t = {npar // 2}")
    print(f"  roots of g: alpha^1 .. alpha^{npar} = {roots}")
    print(f"  g(x) = {g}")
    print(f"  g(x) = {g.to_alpha_str()}")
    T.check("deg(g) = n - k", g.degree == npar)
    T.check("g is monic", g.lead == 1)

    def encode(message):
        shifted = GFPoly(F, list(message)) * GFPoly.monomial(F, npar)
        return shifted + (shifted % g)

    def syndromes(c):
        return [c.eval(r) for r in roots]

    T.sub("Encoding")
    rng = random.Random(1)
    message = [rng.randrange(F.order) for _ in range(k)]
    c = encode(message)
    print(f"  message  = {message}")
    print(f"  parity   = {c.to_list(n)[k:]}")
    print(f"  codeword = {c.to_list(n)}")
    T.check("the codeword has at most n symbols", c.degree < n)
    T.check("systematic: the k high symbols are the message", c.to_list(n)[:k] == message)
    T.check("c is a multiple of g", (c % g).is_zero)
    T.check("syndromes are zero", syndromes(c) == [0] * npar, f"got {syndromes(c)}")
    done("systematic encoding")

    T.sub("Every codeword vanishes at the roots of g")
    all_ok = True
    for _ in range(200):
        msg = [rng.randrange(F.order) for _ in range(k)]
        cw = encode(msg)
        all_ok &= syndromes(cw) == [0] * npar
        all_ok &= cw.to_list(n)[:k] == msg
    T.check("200 codewords checked", all_ok)
    done("the code is the expected kernel")

    T.sub("Single-error correction with the TP classes only")
    print(f"  {'pos':>4} | {'error':>5} | {'S1':>3} {'S2':>3} {'S3':>3} {'S4':>3} | "
          f"{'X=S2/S1':>8} | {'found pos':>9} | {'found val':>9} | ok")
    print(f"  {'-'*4}-+-{'-'*5}-+-{'-'*15}-+-{'-'*8}-+-{'-'*9}-+-{'-'*9}-+----")

    fixed = True
    for pos in range(0, n, 2):
        for val in (1, 7, 13):
            received = c + GFPoly.monomial(F, pos, val)
            S = syndromes(received)

            # single error e at position j: S_i = e * (alpha^j)^i, so X = S2/S1 and e = S1/X
            X = F.div(S[1], S[0])
            found_pos = F.log(X)
            found_val = F.div(S[0], X)

            corrected = received + GFPoly.monomial(F, found_pos, found_val)
            good = found_pos == pos and found_val == val and corrected == c
            fixed &= good
            if val == 7:
                print(f"  {pos:>4} | {val:>5} | {S[0]:>3} {S[1]:>3} {S[2]:>3} {S[3]:>3} | "
                      f"{X:>8} | {found_pos:>9} | {found_val:>9} | {'yes' if good else 'NO'}")

    T.check("any single error is located and corrected", fixed)
    done("the algebra works end to end")

    T.sub("Double errors are detected")
    detected = True
    for p1 in range(0, n, 3):
        for p2 in range(p1 + 1, n, 4):
            received = c + GFPoly.monomial(F, p1, 5) + GFPoly.monomial(F, p2, 9)
            detected &= syndromes(received) != [0] * npar
    T.check("every double error gives a non-zero syndrome", detected)
    done("consistent with d_min = n - k + 1 = 5")


def part11_hardware_helpers(args) -> None:
    """!
    @brief inv_or_zero, mul_matrix and mul_masks.
    """
    T.section("11. Hardware helpers: inv_or_zero, mul_matrix, mul_masks")

    fields = [GaloisField(3, 0b011), GaloisField(4, 0b0011), GaloisField(8, 0x1D)]

    T.sub("inv_or_zero")
    for F in fields:
        T.check(f"GF(2^{F.m}): inv_or_zero(0) = 0", F.inv_or_zero(0) == 0)
        T.check(f"GF(2^{F.m}): inv_or_zero(a) = inverse(a) for a != 0",
                all(F.inv_or_zero(a) == F.inverse(a) for a in range(1, F.order)))
    T.raises("element out of the field", ValueError, fields[0].inv_or_zero, 8)
    done("inv_or_zero matches inverse and maps 0 to 0")

    rng = random.Random(8)
    for F in fields:
        T.sub(f"mul_matrix and mul_masks in GF(2^{F.m}), all {F.order}^2 pairs")
        shape = cols = packed = applied = True
        for c in F.elements:
            M = F.mul_matrix(c)
            K = F.mul_masks(c)
            shape &= len(M) == F.m and all(len(r) == F.m for r in M) and len(K) == F.m
            cols &= all(sum(M[i][j] << i for i in range(F.m)) == F.mul(c, 1 << j)
                        for j in range(F.m))
            packed &= K == tuple(sum(b << j for j, b in enumerate(r)) for r in M)
            for a in F.elements:
                applied &= mat_vec(M, a) == F.mul(c, a) == mask_vec(K, a)
        T.check("m x m matrix and one mask per output bit", shape)
        T.check("column j is c * x^j", cols)
        T.check("masks are the packed rows", packed)
        T.check("M(c) applied to a gives c * a", applied)

        identity = tuple(tuple(int(i == j) for j in range(F.m)) for i in range(F.m))
        T.check("M(1) is the identity", F.mul_matrix(1) == identity)
        T.check("M(0) is all zeros", not any(map(any, F.mul_matrix(0))))

        additive = multiplicative = True
        for _ in range(200):
            c1, c2 = rng.randrange(F.order), rng.randrange(F.order)
            M1, M2 = F.mul_matrix(c1), F.mul_matrix(c2)
            additive &= mat_add(M1, M2) == F.mul_matrix(c1 ^ c2)
            multiplicative &= mat_mul(M1, M2) == F.mul_matrix(F.mul(c1, c2))
        T.check("M(c1 + c2) = M(c1) + M(c2)", additive)
        T.check("M(c1 * c2) = M(c1) * M(c2)", multiplicative)
        done("the matrices represent the field over GF(2)")

    T.raises("constant out of the field", ValueError, fields[1].mul_matrix, 16)

    G = fields[2]
    T.sub("XOR count of the constant multipliers alpha^0 .. alpha^15 in GF(2^8)")
    T.check("alpha^0 = 1 needs no gates", xor_count(G.mul_masks(1)) == 0)
    T.check("alpha^1 needs 3 gates (taps of x^8 + x^4 + x^3 + x^2 + 1)",
            xor_count(G.mul_masks(2)) == 3, f"got {xor_count(G.mul_masks(2))}")
    if not args.no_tables:
        print(f"  {'j':>2} | {'alpha^j':>7} | {'masks (bit 0 .. bit 7)':<23} | XORs")
        print(f"  {'-'*2}-+-{'-'*7}-+-{'-'*23}-+-----")
        for j in range(16):
            c = G.exp(j)
            masks = G.mul_masks(c)
            print(f"  {j:>2} | {c:>#7x} | {hexs(masks):<23} | {xor_count(masks):>4}")
    done("XOR counts")

    F = fields[1]
    if not args.no_tables:
        print(f"\n  Matrix of the product by alpha = {F.alpha} in GF(2^4) (row i = output bit i)")
        for i, row in enumerate(F.mul_matrix(F.alpha)):
            print(f"    bit {i}: {' '.join(str(b) for b in row)}")


def part12_conversions(args) -> None:
    """!
    @brief Fixed-length lists, ascending order and packed integers.
    """
    T.section("12. GFPoly: fixed-length lists and packed integers")

    F4 = GaloisField(4, 0b0011)
    F8 = GaloisField(8, 0x1D)

    T.sub("to_list and from_ascending")
    p = GFPoly(F4, [3, 0, 5, 1])    # 3x^3 + 5x + 1
    T.check("to_list() is the normalized form", p.to_list() == [3, 0, 5, 1])
    T.check("to_list(ascending=True)", p.to_list(ascending=True) == [1, 5, 0, 3])
    T.check("to_list(6) pads the high degrees", p.to_list(6) == [0, 0, 3, 0, 5, 1])
    T.check("to_list(6, ascending=True)", p.to_list(6, ascending=True) == [1, 5, 0, 3, 0, 0])
    T.check("to_list(4) is exact", p.to_list(4) == [3, 0, 5, 1])
    T.check("zero polynomial", GFPoly.zero(F4).to_list() == [0]
            and GFPoly.zero(F4).to_list(3) == [0, 0, 0])
    T.check("from_ascending reverses the order", GFPoly.from_ascending(F4, [1, 5, 0, 3]) == p)
    T.check("from_ascending([]) is zero", GFPoly.from_ascending(F4, []).is_zero)
    T.raises("degree does not fit in length", ValueError, p.to_list, 3)
    T.raises("length = 0", ValueError, p.to_list, 0)
    T.raises("bool length", ValueError, p.to_list, True)
    T.raises("from_ascending with a scalar", TypeError, GFPoly.from_ascending, F4, 5)
    done("list conversions")

    rng = random.Random(4242)
    n = 200 if args.quick else 1000
    T.sub(f"Round trips with {n} random fixed-length words")
    lists_ok = ints_ok = layout_ok = True
    for _ in range(n):
        F = rng.choice((F4, F8))
        length = rng.randint(1, 44)
        word = [rng.randrange(F.order) for _ in range(length)]
        if rng.random() < 0.5:
            word[0] = 0
        p = GFPoly(F, word)
        lists_ok &= p.to_list(length) == word
        lists_ok &= GFPoly.from_ascending(F, word[::-1]) == p
        value = p.to_int(length)
        ints_ok &= GFPoly.from_int(F, value, length) == p and p.to_int() == value
        layout_ok &= value < 1 << (F.m * length)
        layout_ok &= all((value >> (F.m * i)) & (F.order - 1) == p.coeff(i) for i in range(length))
    T.check("to_list(length) keeps the leading zeros", lists_ok)
    T.check("from_int(to_int(p)) = p", ints_ok)
    T.check("coefficient i sits in bits [m*i + m - 1 : m*i]", layout_ok)
    done("round trips")

    T.sub("Invalid input")
    T.raises("value wider than m * length bits", ValueError, GFPoly.from_int, F8, 1 << 16, 2)
    T.raises("negative value", ValueError, GFPoly.from_int, F8, -1, 2)
    T.raises("non-integer value", TypeError, GFPoly.from_int, F8, 1.0, 2)
    T.raises("length = 0", ValueError, GFPoly.from_int, F8, 0, 0)
    T.raises("to_int with a short length", ValueError, GFPoly(F8, [1, 0, 0]).to_int, 2)
    T.check("from_int(0, n) is zero", GFPoly.from_int(F8, 0, 5).is_zero)
    T.check("largest value fits", GFPoly.from_int(F8, 0xFFFF, 2).to_list(2) == [0xFF, 0xFF])
    done("invalid input is rejected")

    T.sub("Example: a word that starts with 0x00")
    word = [0x00, 0x12, 0x34]
    p = GFPoly(F8, word)
    print(f"  word       = {hexs(word)}")
    print(f"  coeffs     = {hexs(p.coeffs)}   <- normalized, the 00 is gone")
    print(f"  to_list(3) = {hexs(p.to_list(3))}")
    print(f"  to_int(3)  = {p.to_int(3):#08x}")
    T.check("coeffs drops the leading 0x00", len(p.coeffs) == 2)
    T.check("to_list(3) restores it", p.to_list(3) == word)
    T.check("packed value is 0x001234", p.to_int(3) == 0x001234)


def part13_derivative(F: GaloisField, args) -> None:
    """!
    @brief Formal derivative in characteristic 2.
    """
    T.section("13. GFPoly: formal derivative")

    T.sub("Monomials and constants")
    T.check("d/dx (c x^k) is c x^(k-1) for odd k and 0 for even k",
            all(GFPoly.monomial(F, k, 7).derivative()
                == (GFPoly.monomial(F, k - 1, 7) if k % 2 else GFPoly.zero(F))
                for k in range(1, 12)))
    T.check("constants have zero derivative", GFPoly(F, [9]).derivative().is_zero)
    T.check("zero has zero derivative", GFPoly.zero(F).derivative().is_zero)
    done("monomials")

    rng = random.Random(1618)
    n = 100 if args.quick else 400
    x = GFPoly.monomial(F, 1)
    T.sub(f"{n} random polynomials over {F}")
    linear = product = second = odd = True
    for _ in range(n):
        a, b = random_poly(F, rng, 8), random_poly(F, rng, 8)
        s = rng.randrange(F.order)
        linear &= (a + b).derivative() == a.derivative() + b.derivative()
        linear &= a.scale(s).derivative() == a.derivative().scale(s)
        product &= (a * b).derivative() == a.derivative() * b + a * b.derivative()
        second &= a.derivative().derivative().is_zero
        odd_part = GFPoly.from_ascending(
            F, [c if i % 2 else 0 for i, c in enumerate(a.to_list(ascending=True))])
        odd &= x * a.derivative() == odd_part
    T.check("linear", linear)
    T.check("product rule", product)
    T.check("second derivative vanishes in characteristic 2", second)
    T.check("x * p'(x) is the odd part of p(x)", odd)
    done(f"{n} polynomials checked")

    T.sub("Simple and repeated roots")
    simple = repeated = True
    for _ in range(n):
        roots = rng.sample(range(F.order), rng.randint(1, 5))
        p = GFPoly.from_roots(F, roots)
        simple &= all(p.derivative().eval(r) != 0 for r in roots)
        q = p * GFPoly(F, [1, roots[0]])
        repeated &= q.derivative().eval(roots[0]) == 0
    T.check("p'(r) != 0 at simple roots", simple)
    T.check("p'(r) = 0 at double roots", repeated)
    done("derivative detects repeated roots")

    p = GFPoly(F, [7, 3, 0, 5, 1])
    print(f"\n  p(x)      = {p}")
    print(f"  p'(x)     = {p.derivative()}")
    print(f"  x * p'(x) = {x * p.derivative()}   <- odd-degree terms of p")
    T.check("worked example", x * p.derivative() == GFPoly(F, [3, 0, 5, 0]))


def part14_qr_code(args) -> None:
    """!
    @brief RS(44, 28) of QR version 2-M (TP3) and the reduced RS(15, 9).
    """
    T.section("14. Reed-Solomon (44, 28) of QR version 2-M over GF(2^8)")

    F = GaloisField(8, 0x1D)
    n, k = 44, 28
    npar = n - k

    g = GFPoly.from_roots(F, [F.exp(i) for i in range(npar)])
    print(f"  g(x) = product of (x + alpha^i), i = 0 .. {npar - 1}")
    print(f"  g    = {hexs(g.coeffs)}")
    T.check("deg(g) = 16 and g is monic", g.degree == npar and g.lead == 1)

    def encode(message):
        shifted = GFPoly(F, message) * GFPoly.monomial(F, npar)
        return shifted + shifted % g

    T.sub("Reference vector: message 00, 01, ..., 1B")
    msg = list(range(k))
    cw = encode(msg)
    word = cw.to_list(n)
    print(f"  parity = {hexs(word[k:])}")
    T.check("parity matches TP3 exercise 6",
            word[k:] == [0xE8, 0x9B, 0xBB, 0x6A, 0x6F, 0x00, 0x3E, 0x3E,
                         0xC9, 0x7B, 0xC7, 0x1E, 0x8E, 0x75, 0x65, 0x38])
    T.check("the leading 0x00 survives with to_list(44)", word[:k] == msg)
    T.check("the word vanishes at alpha^0 .. alpha^15",
            all(cw.eval(F.exp(j)) == 0 for j in range(npar)))

    i_msg = GFPoly(F, msg).to_int(k)
    o_cw = cw.to_int(n)
    T.check("i_msg[223:216] = 0x00 and i_msg[7:0] = 0x1B",
            i_msg >> 216 == 0x00 and i_msg & 0xFF == 0x1B)
    T.check("o_cw[351:128] = i_msg", o_cw >> 128 == i_msg)
    T.check("o_cw[127:0] = parity", o_cw & ((1 << 128) - 1) == GFPoly(F, word[k:]).to_int(npar))
    print(f"  o_cw   = {o_cw:088x}")
    done("reference vector and bus layout")

    rng = random.Random(44)
    trials = 200 if args.quick else 1000
    T.sub(f"{trials} random messages")
    lengths = systematic = codewords = True
    zero_lead = 0
    for _ in range(trials):
        msg = [rng.randrange(256) for _ in range(k)]
        if rng.random() < 0.25:
            msg[0] = 0
        cw = encode(msg)
        word = cw.to_list(n)
        lengths &= len(word) == n
        systematic &= word[:k] == msg
        codewords &= (cw % g).is_zero
        zero_lead += word[k] == 0
    T.check("every codeword has 44 bytes", lengths)
    T.check("systematic", systematic)
    T.check("every codeword is a multiple of g", codewords)
    print(f"  parity starting with 0x00: {zero_lead} of {trials}")
    done("random messages")

    T.sub("Reduced code RS(15, 9) over GF(2^4)")
    F16 = GaloisField(4, 0b0011)
    g16 = GFPoly.from_roots(F16, [F16.exp(i) for i in range(6)])
    print(f"  g = {hexs(g16.coeffs, 1)}")
    T.check("g matches the TP3 control value 1 A F 2 4 3 1",
            list(g16.coeffs) == [1, 0xA, 0xF, 2, 4, 3, 1])
    m16 = GFPoly.from_int(F16, 0x52D09F17C, 9)
    shifted = m16 * GFPoly.monomial(F16, 6)
    c16 = shifted + shifted % g16
    T.check("encoded 0x52D09F17C keeps the message", c16.to_list(15)[:9] == m16.to_list(9))
    T.check("encoded 0x52D09F17C vanishes at alpha^0 .. alpha^5",
            all(c16.eval(F16.exp(j)) == 0 for j in range(6)))
    done("reduced code")


def main() -> int:
    parser = argparse.ArgumentParser(description="Functional test of galois.py")
    parser.add_argument("--quick", action="store_true",
                        help="skip the heaviest exhaustive sweeps")
    parser.add_argument("--no-tables", action="store_true", dest="no_tables",
                        help="do not print field tables")
    args = parser.parse_args()

    start = time.perf_counter()

    print("=" * WIDTH)
    print("  FEC TP1 - functional test of galois.py")
    print("=" * WIDTH)

    try:
        F = GaloisField(4, 0b0011)    # working field: x^4 + x + 1
    except Exception as e:  # noqa: BLE001
        print(f"\n  [FATAL] could not build GF(2^4): {type(e).__name__}: {e}")
        print("  GaloisField is broken at construction; stopping here.")
        return 1

    parts = [
        ("1. construction",    part1_construction,      ()),
        ("2. tables",          part2_tables,            (F,)),
        ("3. axioms",          part3_axioms,            ()),
        ("4. structure",       part4_structure,         ()),
        ("5. degenerate",      part5_degenerate,        (F,)),
        ("6. poly: build",     part6_poly_construction, (F,)),
        ("7. poly: ring",      part7_poly_ring,         (F,)),
        ("8. poly: division",  part8_division,          (F,)),
        ("9. poly: eval",      part9_evaluation,        (F,)),
        ("10. Reed-Solomon",   part10_reed_solomon,     (F,)),
        ("11. hardware",       part11_hardware_helpers, ()),
        ("12. conversions",    part12_conversions,      ()),
        ("13. derivative",     part13_derivative,       (F,)),
        ("14. QR RS(44, 28)",  part14_qr_code,          ()),
    ]

    # an exception in one part must not hide the results of the others
    for name, fn, extra in parts:
        try:
            fn(*extra, args)
        except Exception as e:  # noqa: BLE001
            T.error(name, e)

    return T.summary(time.perf_counter() - start)


if __name__ == "__main__":
    sys.exit(main())