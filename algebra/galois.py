"""!
@file galois.py
@brief GF(2^m) arithmetic and polynomials over GF(2^m).

FEC Course 2026, Fundación Fulgor - TP1.

Field elements are integers in [0, 2^m - 1], where bit i is the coefficient
of x^i. The primitive polynomial is passed without the x^m term, e.g.
x^4 + x + 1 -> 0x3.

GFPoly takes coefficients in descending order [a_n, ..., a_0] and strips
leading zeros. Fixed-length words (codewords, error patterns, syndromes) must
be converted with to_list(length) or to_int(length), which keep them. Packed
integers hold the degree-i coefficient in bits [m*i + m - 1 : m*i].
"""

from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

__all__ = ["GaloisField", "GFPoly"]


class GaloisField:
    """!
    @brief Galois field GF(2^m) defined by a primitive polynomial.

    The exp/log tables use alpha = x as base, which also checks that P(x) is
    primitive. With require_primitive=False any irreducible P(x) is accepted
    and another generator is used for the tables.

    inverse(0) and div(a, 0) raise ZeroDivisionError, inv_or_zero(0) returns 0
    (hardware convention) and power(0, 0) returns 1.
    """

    __slots__ = ("m", "poly", "order", "n", "generator", "is_primitive", "_exp", "_log")

    def __init__(self, m: int, poly: int, require_primitive: bool = True) -> None:
        """!
        @brief Build GF(2^m) from P(x).
        @param m Field degree.
        @param poly P(x) without the x^m term.
        @param require_primitive If True, reject P(x) when x is not a generator.
        @throws TypeError If poly is not an int.
        @throws ValueError If m or poly are out of range, or P(x) is not
                irreducible (or not primitive when required).
        """
        if not isinstance(m, int) or m < 1:
            raise ValueError(f"m must be an integer >= 1, got {m!r}")
        if not isinstance(poly, int):
            raise TypeError(f"poly must be an int, got {type(poly).__name__}")

        self.m = m
        self.order = 1 << m        # number of elements
        self.n = self.order - 1    # order of the multiplicative group

        if not 0 <= poly < self.order:
            raise ValueError(
                f"poly = {poly:#x} does not fit in {m} bits "
                f"(the x^{m} term is implicit)"
            )
        if m > 1 and not poly & 1:
            raise ValueError(f"P(x) = {self._poly_str(poly, extra_degree=m)} is divisible by x")

        self.poly = poly
        self._build_tables(require_primitive)

    def _build_tables(self, require_primitive: bool) -> None:
        """!
        @brief Build the exp/log tables, trying x first as generator.
        """
        x = self._xtime(1)    # equals 1 only when m == 1
        tables = self._try_generator(x)
        if tables is not None:
            self.generator, self.is_primitive = x, True
            self._exp, self._log = tables
            return

        if require_primitive:
            raise ValueError(
                f"P(x) = {self.primitive_poly_str()} is not primitive in GF(2^{self.m}); "
                f"use require_primitive=False to search for another generator"
            )

        for g in range(2, self.order):
            tables = self._try_generator(g)
            if tables is not None:
                self.generator, self.is_primitive = g, False
                self._exp, self._log = tables
                return

        raise ValueError(f"P(x) = {self.primitive_poly_str()} is not irreducible")

    def _try_generator(self, g: int) -> Tuple[List[int], List[int]] | None:
        """!
        @brief Build exp/log tables with base g.
        @return (exp, log), or None if g does not have order 2^m - 1.
        """
        exp = [0] * self.n
        log = [-1] * self.order

        elem = 1
        for i in range(self.n):
            if log[elem] != -1:
                return None    # short cycle
            exp[i] = elem
            log[elem] = i
            elem = self.mul_ref(elem, g)

        return (exp, log) if elem == 1 else None

    def _xtime(self, a: int) -> int:
        """!
        @brief Multiply by x modulo P(x).
        """
        a <<= 1
        if a & self.order:
            a ^= self.order | self.poly
        return a

    def check(self, a: int) -> int:
        """!
        @brief Validate a field element.
        @param a Value to check.
        @return a, unchanged.
        @throws TypeError If a is not an int (bool is rejected).
        @throws ValueError If a is outside [0, 2^m - 1].
        """
        if not isinstance(a, int) or isinstance(a, bool):
            raise TypeError(f"field elements are ints, got {type(a).__name__}: {a!r}")
        if not 0 <= a < self.order:
            raise ValueError(f"{a} is not in GF(2^{self.m}), expected [0, {self.order - 1}]")
        return a

    def __contains__(self, a: object) -> bool:
        return isinstance(a, int) and not isinstance(a, bool) and 0 <= a < self.order

    @property
    def elements(self) -> range:
        """!
        @brief All field elements.
        """
        return range(self.order)

    @property
    def alpha(self) -> int:
        """!
        @brief Base of the exp/log tables (x when P(x) is primitive).
        """
        return self.generator

    def add(self, a: int, b: int) -> int:
        """!
        @brief Field addition (bitwise XOR).
        @param a Field element.
        @param b Field element.
        @return a + b.
        """
        return self.check(a) ^ self.check(b)

    sub = add    # characteristic 2

    def neg(self, a: int) -> int:
        """!
        @brief Additive inverse; -a = a in characteristic 2.
        """
        return self.check(a)

    def mul(self, a: int, b: int) -> int:
        """!
        @brief Field product using the log/antilog tables.
        @param a Field element.
        @param b Field element.
        @return a * b.
        """
        self.check(a)
        self.check(b)
        if a == 0 or b == 0:
            return 0
        return self._exp[(self._log[a] + self._log[b]) % self.n]

    def inverse(self, a: int) -> int:
        """!
        @brief Multiplicative inverse.
        @param a Non-zero field element.
        @return a^-1.
        @throws ZeroDivisionError If a == 0.
        """
        self.check(a)
        if a == 0:
            raise ZeroDivisionError("0 has no multiplicative inverse")
        return self._exp[(self.n - self._log[a]) % self.n]

    def inv_or_zero(self, a: int) -> int:
        """!
        @brief Inverse with the hardware convention inv(0) = 0.

        Matches an inverter implemented as a 2^m-entry lookup table
        (e.g. gf256_inv). Use inverse() for algebra.
        @param a Field element.
        @return a^-1, or 0 if a == 0.
        """
        self.check(a)
        if a == 0:
            return 0
        return self._exp[(self.n - self._log[a]) % self.n]

    def div(self, a: int, b: int) -> int:
        """!
        @brief Field division.
        @param a Dividend.
        @param b Non-zero divisor.
        @return a * b^-1.
        @throws ZeroDivisionError If b == 0.
        """
        self.check(a)
        self.check(b)
        if b == 0:
            raise ZeroDivisionError(f"division by zero in GF(2^{self.m}): {a} / 0")
        if a == 0:
            return 0
        return self._exp[(self._log[a] - self._log[b]) % self.n]

    def power(self, a: int, k: int) -> int:
        """!
        @brief Power a^k. Negative exponents are allowed for a != 0.
        @param a Field element.
        @param k Integer exponent.
        @return a^k, with 0^0 = 1.
        @throws TypeError If k is not an int.
        @throws ZeroDivisionError If a == 0 and k < 0.
        """
        self.check(a)
        if not isinstance(k, int) or isinstance(k, bool):
            raise TypeError(f"exponent must be an int, got {k!r}")
        if a == 0:
            if k == 0:
                return 1
            if k < 0:
                raise ZeroDivisionError("0 raised to a negative power")
            return 0
        return self._exp[(self._log[a] * k) % self.n]

    def exp(self, k: int) -> int:
        """!
        @brief Antilog.
        @param k Integer exponent (any sign).
        @return alpha^k.
        """
        return self._exp[k % self.n]

    def log(self, a: int) -> int:
        """!
        @brief Discrete logarithm in base alpha.
        @param a Non-zero field element.
        @return k in [0, 2^m - 2] such that alpha^k = a.
        @throws ZeroDivisionError If a == 0.
        """
        self.check(a)
        if a == 0:
            raise ZeroDivisionError("log(0) is undefined")
        return self._log[a]

    def mul_ref(self, a: int, b: int) -> int:
        """!
        @brief Reference product: carry-less multiply and reduce modulo P(x).

        Does not use the tables, so it can be used to cross-check mul().
        @param a Field element.
        @param b Field element.
        @return a * b.
        """
        self.check(a)
        self.check(b)
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a = self._xtime(a)
        return result

    def mul_matrix(self, c: int) -> Tuple[Tuple[int, ...], ...]:
        """!
        @brief Binary m x m matrix of the product by a constant.

        M[i][j] is bit i of c * x^j, so that
        bit i of (c * a) = XOR over j of (M[i][j] AND a_j).
        Rows and columns are indexed from bit 0.
        @param c Constant field element.
        @return Tuple of m rows with m bits each.
        """
        self.check(c)
        cols = [self.mul(c, 1 << j) for j in range(self.m)]
        return tuple(
            tuple((cols[j] >> i) & 1 for j in range(self.m))
            for i in range(self.m)
        )

    def mul_masks(self, c: int) -> Tuple[int, ...]:
        """!
        @brief Rows of mul_matrix(c) packed as integer masks.

        bit i of (c * a) = parity(a & masks[i]), which maps directly to
        RTL as assign o_p[i] = ^(i_a & MASK_i). Output bit i needs
        popcount(masks[i]) - 1 two-input XOR gates without sharing.
        @param c Constant field element.
        @return Tuple of m masks, one per output bit.
        """
        return tuple(
            sum(bit << j for j, bit in enumerate(row))
            for row in self.mul_matrix(c)
        )

    @staticmethod
    def _poly_str(value: int, var: str = "x", extra_degree: int | None = None) -> str:
        """!
        @brief Format a binary polynomial stored as an int.
        """
        terms: List[str] = []
        if extra_degree is not None:
            terms.append(f"{var}^{extra_degree}")
        for i in reversed(range(value.bit_length())):
            if value >> i & 1:
                terms.append("1" if i == 0 else (var if i == 1 else f"{var}^{i}"))
        return " + ".join(terms) if terms else "0"

    def poly_str(self, a: int) -> str:
        """!
        @brief Polynomial form of a field element, e.g. "x^3 + 1".
        """
        return self._poly_str(self.check(a))

    def primitive_poly_str(self) -> str:
        """!
        @brief P(x) as a string, including the x^m term.
        """
        return self._poly_str(self.poly, extra_degree=self.m)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, GaloisField):
            return NotImplemented
        return self.m == other.m and self.poly == other.poly

    def __hash__(self) -> int:
        return hash((GaloisField, self.m, self.poly))

    def __len__(self) -> int:
        return self.order

    def __repr__(self) -> str:
        return f"GaloisField(m={self.m}, poly={self.poly:#0{(self.m + 3) // 4 + 2}x})"

    def __str__(self) -> str:
        note = "" if self.is_primitive else f" (not primitive, alpha = {self.generator})"
        return f"GF(2^{self.m}) with P(x) = {self.primitive_poly_str()}{note}"


class GFPoly:
    """!
    @brief Polynomial with coefficients in GF(2^m).

    Coefficients are given in descending order: [a_n, ..., a_1, a_0].
    Leading zeros are stripped and the zero polynomial is stored as (0,),
    with degree -1. Because of this, coeffs cannot rebuild a fixed-length
    word: a message starting with 0x00 loses that byte. Use to_list(length)
    or to_int(length) for that.
    """

    __slots__ = ("field", "coeffs")

    def __init__(self, field: GaloisField, coeffs: Sequence[int]) -> None:
        """!
        @brief Build a polynomial from its coefficients.
        @param field Coefficient field.
        @param coeffs Coefficients in descending order of degree.
        @throws TypeError If field is not a GaloisField or coeffs is a scalar.
        @throws ValueError If a coefficient is not a field element.
        """
        if not isinstance(field, GaloisField):
            raise TypeError(f"field must be a GaloisField, got {type(field).__name__}")
        if isinstance(coeffs, (int, bool)):
            raise TypeError("coeffs must be a list or tuple")

        values = list(coeffs)
        for pos, c in enumerate(values):
            try:
                field.check(c)
            except (TypeError, ValueError) as exc:
                raise type(exc)(
                    f"invalid coefficient at position {pos} "
                    f"(degree {len(values) - 1 - pos}): {exc}"
                ) from None

        start = 0
        while start < len(values) and values[start] == 0:
            start += 1
        values = values[start:]

        self.field = field
        self.coeffs: Tuple[int, ...] = tuple(values) if values else (0,)

    @classmethod
    def zero(cls, field: GaloisField) -> GFPoly:
        """!
        @brief Zero polynomial.
        """
        return cls(field, [0])

    @classmethod
    def one(cls, field: GaloisField) -> GFPoly:
        """!
        @brief Constant polynomial 1.
        """
        return cls(field, [1])

    @classmethod
    def monomial(cls, field: GaloisField, degree: int, coeff: int = 1) -> GFPoly:
        """!
        @brief Monomial coeff * x^degree.
        @throws ValueError If degree < 0.
        """
        if degree < 0:
            raise ValueError(f"degree must be >= 0, got {degree}")
        return cls(field, [coeff] + [0] * degree)

    @classmethod
    def from_roots(cls, field: GaloisField, roots: Iterable[int]) -> GFPoly:
        """!
        @brief Product of (x + r) over the given roots.
        @param field Coefficient field.
        @param roots Field elements; an empty list gives 1.
        @return The monic polynomial with those roots.
        """
        result = cls.one(field)
        for r in roots:
            field.check(r)
            result = result * cls(field, [1, r])
        return result

    @classmethod
    def from_ascending(cls, field: GaloisField, coeffs: Sequence[int]) -> GFPoly:
        """!
        @brief Build a polynomial from coefficients in ascending order.
        @param field Coefficient field.
        @param coeffs [a_0, a_1, ..., a_n]; index equals degree.
        @throws TypeError If coeffs is a scalar.
        """
        if isinstance(coeffs, (int, bool)):
            raise TypeError("coeffs must be a list or tuple")
        return cls(field, list(coeffs)[::-1])

    @classmethod
    def from_int(cls, field: GaloisField, value: int, length: int) -> GFPoly:
        """!
        @brief Unpack length m-bit coefficients from an integer.

        The degree-i coefficient is read from bits [m*i + m - 1 : m*i], the
        bus convention of TP3 (e.g. o_cw[7:0] is the x^0 coefficient).
        @param field Coefficient field.
        @param value Packed value.
        @param length Number of coefficients.
        @throws TypeError If value is not an int.
        @throws ValueError If value is negative, length < 1, or value has bits
                above m * length.
        """
        if not isinstance(value, int) or isinstance(value, bool):
            raise TypeError(f"value must be an int, got {type(value).__name__}")
        if value < 0:
            raise ValueError(f"value must be >= 0, got {value}")
        if not isinstance(length, int) or isinstance(length, bool) or length < 1:
            raise ValueError(f"length must be an int >= 1, got {length!r}")
        m = field.m
        if value >> (m * length):
            raise ValueError(f"{value:#x} does not fit in {length} coefficients of {m} bits")
        mask = field.order - 1
        return cls.from_ascending(field, [(value >> (m * i)) & mask for i in range(length)])

    @property
    def degree(self) -> int:
        """!
        @brief Degree; -1 for the zero polynomial.
        """
        return -1 if self.is_zero else len(self.coeffs) - 1

    @property
    def is_zero(self) -> bool:
        return self.coeffs == (0,)

    @property
    def lead(self) -> int:
        """!
        @brief Leading coefficient (0 for the zero polynomial).
        """
        return self.coeffs[0]

    def coeff(self, degree: int) -> int:
        """!
        @brief Coefficient of x^degree, or 0 if out of range.
        """
        if degree < 0 or degree > self.degree:
            return 0
        return self.coeffs[self.degree - degree]

    def __len__(self) -> int:
        return len(self.coeffs)

    def _same_field(self, other: GFPoly) -> None:
        if self.field != other.field:
            raise ValueError(f"polynomials over different fields: {self.field} and {other.field}")

    def to_list(self, length: int | None = None, ascending: bool = False) -> List[int]:
        """!
        @brief Coefficients as a list, optionally zero-padded to a fixed length.
        @param length Number of coefficients (degrees 0 to length - 1). If None,
               the normalized form is returned.
        @param ascending If True, element i is the coefficient of x^i.
        @return List of coefficients.
        @throws ValueError If length < 1 or the degree does not fit in length.
        """
        if length is None:
            values = list(self.coeffs)
        else:
            if not isinstance(length, int) or isinstance(length, bool) or length < 1:
                raise ValueError(f"length must be an int >= 1, got {length!r}")
            if self.degree >= length:
                raise ValueError(f"degree {self.degree} does not fit in {length} coefficients")
            values = [0] * (length - len(self.coeffs)) + list(self.coeffs)
        return values[::-1] if ascending else values

    def to_int(self, length: int | None = None) -> int:
        """!
        @brief Pack the coefficients into an integer.

        The degree-i coefficient goes to bits [m*i + m - 1 : m*i].
        @param length If given, only checks that the polynomial fits.
        @return Packed value.
        @throws ValueError If the degree does not fit in length.
        """
        m = self.field.m
        value = 0
        for i, c in enumerate(self.to_list(length, ascending=True)):
            value |= c << (m * i)
        return value

    def __add__(self, other: GFPoly) -> GFPoly:
        if not isinstance(other, GFPoly):
            return NotImplemented
        self._same_field(other)
        size = max(len(self.coeffs), len(other.coeffs))
        a = (0,) * (size - len(self.coeffs)) + self.coeffs
        b = (0,) * (size - len(other.coeffs)) + other.coeffs
        return GFPoly(self.field, [x ^ y for x, y in zip(a, b)])

    __sub__ = __add__    # characteristic 2

    def __neg__(self) -> GFPoly:
        return self

    def scale(self, scalar: int) -> GFPoly:
        """!
        @brief Multiply every coefficient by a field element.
        """
        self.field.check(scalar)
        if scalar == 0:
            return GFPoly.zero(self.field)
        return GFPoly(self.field, [self.field.mul(c, scalar) for c in self.coeffs])

    def __mul__(self, other: GFPoly | int) -> GFPoly:
        if isinstance(other, int) and not isinstance(other, bool):
            return self.scale(other)
        if not isinstance(other, GFPoly):
            return NotImplemented
        self._same_field(other)
        if self.is_zero or other.is_zero:
            return GFPoly.zero(self.field)

        field = self.field
        result = [0] * (len(self.coeffs) + len(other.coeffs) - 1)
        for i, a in enumerate(self.coeffs):
            if a == 0:
                continue
            for j, b in enumerate(other.coeffs):
                if b:
                    result[i + j] ^= field.mul(a, b)
        return GFPoly(field, result)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> GFPoly:
        """!
        @brief Power by repeated squaring.
        @throws ValueError If exponent is not an int >= 0.
        """
        if not isinstance(exponent, int) or isinstance(exponent, bool) or exponent < 0:
            raise ValueError(f"exponent must be an int >= 0, got {exponent!r}")
        result = GFPoly.one(self.field)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent >>= 1
        return result

    def divmod(self, other: GFPoly) -> Tuple[GFPoly, GFPoly]:
        """!
        @brief Polynomial long division.
        @param other Non-zero divisor.
        @return (q, r) with self = q * other + r and deg(r) < deg(other).
        @throws TypeError If other is not a GFPoly.
        @throws ZeroDivisionError If other is the zero polynomial.
        """
        if not isinstance(other, GFPoly):
            raise TypeError(f"divisor must be a GFPoly, got {type(other).__name__}")
        self._same_field(other)
        if other.is_zero:
            raise ZeroDivisionError("polynomial division by zero")

        field = self.field
        if self.degree < other.degree:
            return GFPoly.zero(field), self

        rem = list(self.coeffs)
        divisor = other.coeffs
        inv_lead = field.inverse(other.lead)
        steps = len(rem) - len(divisor) + 1
        quotient = [0] * steps

        for i in range(steps):
            factor = rem[i]
            if factor == 0:
                continue
            factor = field.mul(factor, inv_lead)
            quotient[i] = factor
            for j, d in enumerate(divisor):
                if d:
                    rem[i + j] ^= field.mul(factor, d)

        return GFPoly(field, quotient), GFPoly(field, rem[steps:])

    def __divmod__(self, other: GFPoly) -> Tuple[GFPoly, GFPoly]:
        return self.divmod(other)

    def __floordiv__(self, other: GFPoly) -> GFPoly:
        return self.divmod(other)[0]

    def __mod__(self, other: GFPoly) -> GFPoly:
        return self.divmod(other)[1]

    def derivative(self) -> GFPoly:
        """!
        @brief Formal derivative.

        In characteristic 2, i * a_i is a_i for odd i and 0 for even i, so
        only odd-degree terms survive and x * p'(x) is the odd part of p(x).
        @return p'(x); zero for constants.
        """
        if self.degree < 1:
            return GFPoly.zero(self.field)
        asc = self.to_list(ascending=True)
        return GFPoly.from_ascending(
            self.field, [asc[i] if i % 2 else 0 for i in range(1, len(asc))]
        )

    def eval(self, x: int) -> int:
        """!
        @brief Evaluate at x with Horner's rule.
        @param x Field element.
        @return p(x).
        """
        field = self.field
        field.check(x)
        acc = 0
        for c in self.coeffs:
            acc = field.mul(acc, x) ^ c
        return acc

    def __call__(self, x: int) -> int:
        return self.eval(x)

    def roots(self) -> List[int]:
        """!
        @brief Roots in the field, by exhaustive search.
        """
        return [x for x in self.field.elements if self.eval(x) == 0]

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, GFPoly):
            return NotImplemented
        return self.field == other.field and self.coeffs == other.coeffs

    def __hash__(self) -> int:
        return hash((GFPoly, self.field, self.coeffs))

    def _format(self, coeff_str) -> str:
        if self.is_zero:
            return "0"
        terms: List[str] = []
        for i, c in enumerate(self.coeffs):
            if c == 0:
                continue
            deg = self.degree - i
            var = "" if deg == 0 else ("x" if deg == 1 else f"x^{deg}")
            if c == 1 and deg != 0:
                terms.append(var)
            elif deg == 0:
                terms.append(coeff_str(c))
            else:
                terms.append(f"{coeff_str(c)}·{var}")
        return " + ".join(terms)

    def __str__(self) -> str:
        return self._format(str)

    def to_alpha_str(self) -> str:
        """!
        @brief Like str(), with coefficients written as powers of alpha.
        """
        def as_alpha(c: int) -> str:
            k = self.field.log(c)
            if k == 0:
                return "1"
            return "α" if k == 1 else f"α^{k}"
        return self._format(as_alpha)

    def __repr__(self) -> str:
        return f"GFPoly({self.field!r}, {list(self.coeffs)})"