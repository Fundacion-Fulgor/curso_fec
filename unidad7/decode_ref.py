"""!
@file decode_ref.py
@brief Chien search, Forney formula and a complete reference decoder.

For syndromes S_j = r(alpha^(b + j)) and an error at X = alpha^p, the Forney
formula depends on the evaluator delivered by the key equation solver:

  classic, Omega = Lambda S mod x^t      Y = X^(1 - b)        Omega(X^-1)   / Lambda'(X^-1)   (eq. 7)
  high, Omega_h from riBM / RiBM         Y = X^-(b + 2t - 1)  Omega_h(X^-1) / Lambda'(X^-1)   (eq. 12)

Both are invariant to scaling Lambda and Omega by the same constant, so the
scaled outputs of the inversionless solvers can be used directly.

The flags follow rs_decoder of TP3: fail when the number of roots in the n
positions differs from the degree of Lambda. With more than t errors a word
can pass that test and still not be corrected. DecodeResult.consistent and
KesResult.length expose those cases: they all have deg Lambda != L.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

from galois import GaloisField, GFPoly
from kes_ref import KesResult, make_kes
from rs_code import RSCode

__all__ = ["ChienForneyResult", "DecodeResult", "forney_exponent", "chien_forney", "decode"]


@dataclass
class ChienForneyResult:
    """!
    @brief Output of the Chien search and error evaluation.

    errors is the error pattern as a word (highest degree first), positions
    are the roots found in ascending order and values their error values.
    fail is set when the number of roots differs from the degree of Lambda or
    a root is repeated.
    """

    errors: List[int]
    positions: List[int]
    values: List[int]
    fail: bool


@dataclass
class DecodeResult:
    """!
    @brief Output of the complete decoder, with every intermediate result.

    codeword and message are only meaningful when uncorrectable is False.
    consistent tells whether the corrected word has zero syndromes; the
    hardware flags do not check it.
    """

    syndromes: List[int]
    kes: KesResult
    errors: List[int]
    positions: List[int]
    values: List[int]
    codeword: List[int]
    message: List[int]
    corrected: bool
    uncorrectable: bool
    consistent: bool


def forney_exponent(form: str, first_root: int, t: int) -> int:
    """!
    @brief Exponent e such that Y = X^e * Omega(X^-1) / Lambda'(X^-1).
    @param form "classic" or "high".
    @param first_root Exponent b of the first root of the generator.
    @param t Error correcting capability.
    @throws ValueError If form is unknown.
    """
    if form == "classic":
        return 1 - first_root
    if form == "high":
        return -(first_root + 2 * t - 1)
    raise ValueError(f"form must be 'classic' or 'high', got {form!r}")


def chien_forney(field: GaloisField, n: int, lam: Sequence[int], omega: Sequence[int],
                 form: str, first_root: int = 0) -> ChienForneyResult:
    """!
    @brief Chien search over the n positions of the word and Forney formula.

    Only positions 0 .. n - 1 are searched, so a root outside a shortened
    code makes the root count smaller than the degree and sets fail.
    Lambda'(x) only keeps the odd-degree terms of Lambda.
    @param field Symbol field.
    @param n Code length.
    @param lam t + 1 locator coefficients, ascending.
    @param omega t evaluator coefficients, ascending.
    @param form Evaluator form, "classic" or "high".
    @param first_root Exponent b of the first root of the generator.
    @throws ValueError If the lengths do not match, n is out of range or form is unknown.
    """
    lam, omega = list(lam), list(omega)
    t = len(lam) - 1
    if t < 1 or len(omega) != t:
        raise ValueError(f"expected t + 1 locator and t evaluator coefficients, "
                         f"got {len(lam)} and {len(omega)}")
    if not isinstance(n, int) or isinstance(n, bool) or not 1 <= n <= field.n:
        raise ValueError(f"n must be in [1, {field.n}], got {n!r}")
    exponent = forney_exponent(form, first_root, t)

    locator = GFPoly.from_ascending(field, lam)
    derivative = locator.derivative()
    evaluator = GFPoly.from_ascending(field, omega)

    errors = [0] * n
    positions: List[int] = []
    values: List[int] = []
    fail = False
    for p in range(n):
        z = field.exp(-p)
        if locator.eval(z):
            continue
        den = derivative.eval(z)
        if den == 0:
            fail = True    # repeated root
            value = 0
        else:
            value = field.mul(field.exp(exponent * p), field.div(evaluator.eval(z), den))
        positions.append(p)
        values.append(value)
        errors[n - 1 - p] = value

    if len(positions) != locator.degree:
        fail = True
    return ChienForneyResult(errors, positions, values, fail)


def decode(code: RSCode, received: Sequence[int], architecture: str = "ribm_systolic") -> DecodeResult:
    """!
    @brief Syndromes, key equation solver, Chien search and Forney formula.

    The flags follow rs_decoder: corrected when at least one error was fixed,
    uncorrectable when the Chien search fails.
    @param code Reed-Solomon code with t >= 1.
    @param received n symbols, highest degree first.
    @param architecture Solver name, see kes_ref.ARCHITECTURES.
    """
    word = list(received)
    syndromes = code.syndromes(word)
    kes = make_kes(architecture, code.field, code.t).run(syndromes[: 2 * code.t])
    cf = chien_forney(code.field, code.n, kes.lam, kes.omega, kes.form, code.first_root)
    codeword = [a ^ e for a, e in zip(word, cf.errors)]
    return DecodeResult(
        syndromes=syndromes,
        kes=kes,
        errors=cf.errors,
        positions=cf.positions,
        values=cf.values,
        codeword=codeword,
        message=codeword[: code.k],
        corrected=any(cf.errors) and not cf.fail,
        uncorrectable=cf.fail,
        consistent=code.is_codeword(codeword),
    )
