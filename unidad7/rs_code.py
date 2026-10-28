"""!
@file rs_code.py
@brief Reed-Solomon codes over GF(2^m): generator, systematic encoding and syndromes.

Words are lists of n symbols with the highest degree first: word[0] is the
coefficient of x^(n-1) and word[n-1] the one of x^0, as in clause 7.5.2 of
ISO/IEC 18004. Position p of a word is the coefficient of x^p, so it sits at
index n - 1 - p. Buses hold the degree-i coefficient in bits [m*i + m - 1 : m*i].
"""

from __future__ import annotations

from typing import Dict, List, Sequence

from galois import GaloisField, GFPoly

__all__ = ["RSCode"]


class RSCode:
    """!
    @brief Reed-Solomon code, possibly shortened (n < 2^m - 1).

    The generator is g(x) = prod (x + alpha^(b + i)) for i = 0 .. n - k - 1,
    with b = first_root, and the syndromes are S_j = r(alpha^(b + j)).
    """

    def __init__(self, field: GaloisField, n: int, k: int, first_root: int = 0) -> None:
        """!
        @brief Build the code.
        @param field Symbol field GF(2^m).
        @param n Code length, from 2 to 2^m - 1.
        @param k Message length, from 1 to n - 1.
        @param first_root Exponent b of the first root of g(x).
        @throws TypeError If field is not a GaloisField or n, k, first_root are not ints.
        @throws ValueError If n or k are out of range.
        """
        if not isinstance(field, GaloisField):
            raise TypeError(f"field must be a GaloisField, got {type(field).__name__}")
        for name, value in (("n", n), ("k", k), ("first_root", first_root)):
            if not isinstance(value, int) or isinstance(value, bool):
                raise TypeError(f"{name} must be an int, got {value!r}")
        if not 2 <= n <= field.n:
            raise ValueError(f"n must be in [2, {field.n}], got {n}")
        if not 1 <= k < n:
            raise ValueError(f"k must be in [1, {n - 1}], got {k}")

        self.field = field
        self.n = n
        self.k = k
        self.nsym = n - k
        self.t = self.nsym // 2
        self.first_root = first_root
        self.roots = [field.exp(first_root + i) for i in range(self.nsym)]
        self.generator = GFPoly.from_roots(field, self.roots)

    @classmethod
    def qr_2m(cls) -> RSCode:
        """!
        @brief RS(44, 28) of a version 2-M QR symbol (ISO/IEC 18004, Table 9).
        """
        return cls(GaloisField(8, 0x1D), 44, 28)

    @classmethod
    def tp3_reduced(cls) -> RSCode:
        """!
        @brief Reduced RS(15, 9) over GF(2^4) of TP3, part A.
        """
        return cls(GaloisField(4, 0x3), 15, 9)

    @property
    def d_min(self) -> int:
        """!
        @brief Minimum distance n - k + 1 (the code is MDS).
        """
        return self.nsym + 1

    @property
    def rate(self) -> float:
        return self.k / self.n

    def _symbols(self, values: Sequence[int], length: int, name: str) -> List[int]:
        symbols = list(values)
        if len(symbols) != length:
            raise ValueError(f"{name} must have {length} symbols, got {len(symbols)}")
        for s in symbols:
            self.field.check(s)
        return symbols

    def parity(self, message: Sequence[int]) -> List[int]:
        """!
        @brief Parity symbols: remainder of m(x) * x^(n-k) divided by g(x).
        @param message k symbols, highest degree first.
        @return n - k parity symbols, highest degree first.
        @throws ValueError If the message has the wrong length or invalid symbols.
        """
        msg = self._symbols(message, self.k, "message")
        shifted = GFPoly(self.field, msg + [0] * self.nsym)
        return (shifted % self.generator).to_list(self.nsym)

    def encode(self, message: Sequence[int]) -> List[int]:
        """!
        @brief Systematic encoding.
        @param message k symbols, highest degree first.
        @return Codeword of n symbols: the message followed by the parity.
        """
        msg = self._symbols(message, self.k, "message")
        return msg + self.parity(msg)

    def syndromes(self, word: Sequence[int]) -> List[int]:
        """!
        @brief Syndromes S_j = r(alpha^(b + j)), j = 0 .. n - k - 1.
        @param word n symbols, highest degree first.
        @return List of n - k syndromes.
        """
        r = GFPoly(self.field, self._symbols(word, self.n, "word"))
        return [r.eval(x) for x in self.roots]

    def is_codeword(self, word: Sequence[int]) -> bool:
        return not any(self.syndromes(word))

    def message(self, word: Sequence[int]) -> List[int]:
        """!
        @brief Message part of a codeword (its first k symbols).
        """
        return self._symbols(word, self.n, "word")[: self.k]

    def error_word(self, errors: Dict[int, int]) -> List[int]:
        """!
        @brief Error pattern as a word.
        @param errors Mapping position -> value, where position p is the
               coefficient of x^p.
        @return n symbols with the errors and zeros elsewhere.
        @throws ValueError If a position is outside [0, n - 1] or a value is
                not a field element.
        """
        word = [0] * self.n
        for pos, value in errors.items():
            if not isinstance(pos, int) or isinstance(pos, bool) or not 0 <= pos < self.n:
                raise ValueError(f"position must be in [0, {self.n - 1}], got {pos!r}")
            word[self.n - 1 - pos] = self.field.check(value)
        return word

    def add_errors(self, word: Sequence[int], errors: Dict[int, int]) -> List[int]:
        """!
        @brief Add an error pattern to a word.
        @param word n symbols, highest degree first.
        @param errors Mapping position -> value.
        @return word + e(x).
        """
        w = self._symbols(word, self.n, "word")
        return [a ^ e for a, e in zip(w, self.error_word(errors))]

    def pack(self, symbols: Sequence[int], ascending: bool = False) -> int:
        """!
        @brief Pack symbols into a bus value.

        The degree-i symbol goes to bits [m*i + m - 1 : m*i]. Words and
        messages are given highest degree first; syndromes and solver outputs
        are usually given in ascending order.
        @param symbols Symbols to pack; at least one.
        @param ascending True if symbols[i] is the coefficient of x^i.
        @return Packed value.
        """
        values = list(symbols)
        poly = GFPoly.from_ascending(self.field, values) if ascending else GFPoly(self.field, values)
        return poly.to_int(len(values))

    def unpack(self, value: int, count: int, ascending: bool = False) -> List[int]:
        """!
        @brief Inverse of pack().
        @param value Bus value.
        @param count Number of symbols.
        @param ascending Order of the returned list.
        """
        return GFPoly.from_int(self.field, value, count).to_list(count, ascending=ascending)

    def __repr__(self) -> str:
        return (f"RSCode({self.field!r}, n={self.n}, k={self.k}, "
                f"first_root={self.first_root})")

    def __str__(self) -> str:
        return f"RS({self.n}, {self.k}) over GF(2^{self.field.m}), t = {self.t}"
