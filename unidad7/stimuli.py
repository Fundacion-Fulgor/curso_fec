"""!
@file stimuli.py
@brief Reproducible test case generators for the Reed-Solomon blocks.

A Case carries the message, the codeword, the injected errors, the received
word and its syndromes, so the same case can drive the encoder, the syndrome
calculator, the key equation solver or the whole decoder, and the result can
be checked against the known error pattern.

Positions follow rs_code: position p is the coefficient of x^p, data symbols
are at positions n - k .. n - 1 and parity symbols at 0 .. n - k - 1.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Sequence

from rs_code import RSCode

__all__ = ["Case", "StimulusGenerator"]

MESSAGE_PATTERNS = ("random", "zeros", "ones", "ramp")


@dataclass
class Case:
    """!
    @brief One test case.

    errors maps position -> non-zero value. correctable is True when the
    number of errors is at most t.
    """

    kind: str
    message: List[int]
    codeword: List[int]
    errors: Dict[int, int]
    received: List[int]
    syndromes: List[int]
    correctable: bool

    @property
    def n_errors(self) -> int:
        return len(self.errors)


class StimulusGenerator:
    """!
    @brief Generator of test cases for one code, driven by a seeded RNG.
    """

    def __init__(self, code: RSCode, seed: int = 0) -> None:
        """!
        @brief Bind the generator to a code.
        @param code Reed-Solomon code.
        @param seed Seed of the random number generator.
        """
        self.code = code
        self.rng = random.Random(seed)

    def message(self, pattern: str = "random") -> List[int]:
        """!
        @brief Message of k symbols.
        @param pattern "random", "zeros", "ones" (all symbols 2^m - 1) or
               "ramp" (0, 1, 2, ...).
        @throws ValueError If the pattern is unknown.
        """
        k, order = self.code.k, self.code.field.order
        if pattern == "random":
            return [self.rng.randrange(order) for _ in range(k)]
        if pattern == "zeros":
            return [0] * k
        if pattern == "ones":
            return [order - 1] * k
        if pattern == "ramp":
            return [i % order for i in range(k)]
        raise ValueError(f"unknown message pattern {pattern!r}, expected one of {MESSAGE_PATTERNS}")

    def case(self, errors: Dict[int, int], message: Optional[Sequence[int]] = None,
             kind: str = "custom") -> Case:
        """!
        @brief Build a case from an explicit error pattern.
        @param errors Mapping position -> value; zero values are dropped.
        @param message k symbols; random if None.
        @param kind Label of the case.
        """
        code = self.code
        code.error_word(errors)    # validates positions and values
        errs = {p: v for p, v in errors.items() if v}
        msg = self.message() if message is None else list(message)
        codeword = code.encode(msg)
        received = code.add_errors(codeword, errs)
        return Case(kind, msg, codeword, errs, received, code.syndromes(received),
                    len(errs) <= code.t)

    def clean(self, pattern: str = "random") -> Case:
        """!
        @brief Codeword without errors.
        """
        return self.case({}, self.message(pattern), f"clean_{pattern}")

    def random_errors(self, count: int, positions: Optional[Sequence[int]] = None,
                      values: Optional[Sequence[int]] = None, kind: Optional[str] = None) -> Case:
        """!
        @brief count errors at distinct random positions.
        @param count Number of errors.
        @param positions Allowed positions; all of them if None.
        @param values Allowed non-zero values; all of them if None.
        @param kind Label; "errors_<count>" if None.
        @throws ValueError If count does not fit in the positions or a value is 0.
        """
        pool = list(range(self.code.n)) if positions is None else list(positions)
        if not 0 <= count <= len(pool):
            raise ValueError(f"count must be in [0, {len(pool)}], got {count}")
        vals = list(range(1, self.code.field.order)) if values is None else list(values)
        if not vals or 0 in vals:
            raise ValueError("values must be a non-empty list of non-zero symbols")
        errors = {p: self.rng.choice(vals) for p in self.rng.sample(pool, count)}
        return self.case(errors, kind=kind or f"errors_{count}")

    def data_errors(self, count: int) -> Case:
        """!
        @brief Errors only in the data symbols.
        """
        code = self.code
        return self.random_errors(count, range(code.nsym, code.n), kind=f"data_{count}")

    def parity_errors(self, count: int) -> Case:
        """!
        @brief Errors only in the parity symbols.
        """
        return self.random_errors(count, range(self.code.nsym), kind=f"parity_{count}")

    def extreme_values(self, count: int) -> Case:
        """!
        @brief Error values limited to 1 and 2^m - 1.
        """
        top = self.code.field.order - 1
        return self.random_errors(count, values=[1, top], kind=f"extreme_{count}")

    def edge_errors(self) -> Case:
        """!
        @brief Errors at the first and the last position of the word.
        """
        top = self.code.field.order - 1
        return self.case({0: self.rng.randrange(1, top + 1),
                          self.code.n - 1: self.rng.randrange(1, top + 1)}, kind="edges")

    def burst(self, length: int, start: Optional[int] = None) -> Case:
        """!
        @brief Errors at consecutive positions start .. start + length - 1.
        @param length Burst length, from 1 to n.
        @param start First position; random if None.
        @throws ValueError If the burst does not fit in the word.
        """
        n = self.code.n
        if not 1 <= length <= n:
            raise ValueError(f"length must be in [1, {n}], got {length}")
        if start is None:
            start = self.rng.randrange(n - length + 1)
        if not 0 <= start <= n - length:
            raise ValueError(f"start must be in [0, {n - length}], got {start}")
        order = self.code.field.order
        errors = {p: self.rng.randrange(1, order) for p in range(start, start + length)}
        return self.case(errors, kind=f"burst_{length}")

    def random_syndromes(self) -> List[int]:
        """!
        @brief 2t random syndromes. Every vector is the syndrome of some
               error pattern, so this stresses the solver beyond t errors.
        """
        order = self.code.field.order
        return [self.rng.randrange(order) for _ in range(2 * self.code.t)]

    def suite(self, per_count: int = 4) -> List[Case]:
        """!
        @brief Deterministic set covering the main categories.

        Clean words with every message pattern; for each number of errors
        from 1 to t, random, data-only, parity-only and extreme-value cases;
        errors at both edges; a burst of t; and patterns with t + 1, t + 2,
        2t and n errors.
        @param per_count Random cases for each number of errors up to t.
        """
        code, t = self.code, self.code.t
        cases = [self.clean(p) for p in MESSAGE_PATTERNS]
        for e in range(1, t + 1):
            cases += [self.random_errors(e) for _ in range(per_count)]
            cases.append(self.data_errors(min(e, code.k)))
            cases.append(self.parity_errors(min(e, code.nsym)))
            cases.append(self.extreme_values(e))
        cases.append(self.edge_errors())
        cases.append(self.burst(t))
        for e in sorted({t + 1, t + 2, 2 * t, code.n}):
            cases.append(self.random_errors(min(e, code.n), kind=f"beyond_{e}"))
        return cases

    def stream(self, count: int, weights: Optional[Dict[int, float]] = None) -> Iterator[Case]:
        """!
        @brief Random mixture of cases.
        @param count Number of cases.
        @param weights Relative weight of each number of errors. By default
               no errors has weight t, each count from 1 to t weight 1, and
               t + 1 and 2t weight 1 each.
        @throws ValueError If a number of errors does not fit in the word or
                no weight is positive.
        """
        t, n = self.code.t, self.code.n
        if weights is None:
            weights = {0: t, **{e: 1 for e in range(1, t + 1)}, t + 1: 1, 2 * t: 1}
        counts = list(weights)
        if any(not 0 <= e <= n for e in counts):
            raise ValueError(f"numbers of errors must be in [0, {n}]")
        if not any(w > 0 for w in weights.values()):
            raise ValueError("at least one weight must be positive")
        for _ in range(count):
            e = self.rng.choices(counts, weights=[weights[c] for c in counts])[0]
            yield self.random_errors(e)
