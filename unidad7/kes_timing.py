"""!
@file kes_timing.py
@brief Cycle-accurate port model of key_equation_solver, following
       key_equation_solver/docs/kes_protocol.md.

Each call to clock() is one rising edge of i_clk with the inputs sampled on
that edge, and returns the registered outputs after it. The datapath is the
register-level model of kes_ref, so the internal state can be compared with
the RTL cycle by cycle while a word is being processed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

from galois import GaloisField, GFPoly
from kes_ref import make_kes

__all__ = ["KesPorts", "KesTimed", "simulate"]


@dataclass
class KesPorts:
    """!
    @brief Outputs after one clock edge.

    o_lambda and o_omega are ascending coefficient lists; busy and state are
    internal and only meant for debugging (state is None when idle).
    """

    o_valid: bool
    o_lambda: List[int]
    o_omega: List[int]
    busy: bool
    state: Optional[int]


class KesTimed:
    """!
    @brief key_equation_solver with the clocked interface of TP3.

    On every edge, in priority order:
      - i_rst clears everything; o_valid, o_lambda and o_omega go to 0.
      - o_valid drops, unless the step of this edge is the last one.
      - When busy, the solver advances one step. After the last step the
        output registers take the result and o_valid is set for one cycle.
      - i_valid loads the syndromes (state 0) and starts a new word. If a
        word was still in progress it is dropped. Loading on the same edge as
        the last step of the previous word gives back-to-back operation.
    """

    def __init__(self, architecture: str, field: GaloisField, t: int) -> None:
        """!
        @brief Build the model.
        @param architecture "ibm", "ribm" or "ribm_systolic".
        @param field Symbol field.
        @param t Error correcting capability.
        """
        self.kes = make_kes(architecture, field, t)
        self.field = field
        self.t = t
        self._clear()

    @property
    def latency(self) -> int:
        """!
        @brief Edges from the edge that loads a word to the edge that sets o_valid.
        """
        return self.kes.iterations

    @property
    def interval(self) -> int:
        """!
        @brief Minimum number of edges between two loads without dropping a word.
        """
        return self.kes.iterations

    def _clear(self) -> None:
        self.kes.reset()
        self.busy = False
        self.o_valid = False
        self.o_lambda = [0] * (self.t + 1)
        self.o_omega = [0] * self.t

    def ports(self) -> KesPorts:
        return KesPorts(self.o_valid, list(self.o_lambda), list(self.o_omega), self.busy,
                        self.kes.r if self.busy else None)

    def clock(self, rst: bool = False, valid: bool = False,
              syndromes: Optional[Sequence[int]] = None) -> KesPorts:
        """!
        @brief One rising edge.
        @param rst Synchronous reset, active high.
        @param valid i_valid.
        @param syndromes S_0 .. S_(2t-1); required when valid is set.
        @return Outputs after the edge.
        @throws ValueError If valid is set without syndromes, or they are invalid.
        """
        if rst:
            self._clear()
            return self.ports()
        if valid and syndromes is None:
            raise ValueError("i_valid is set but no syndromes were given")

        self.o_valid = False
        if self.busy:
            self.kes.step()
            if self.kes.done:
                res = self.kes.result()
                self.o_lambda, self.o_omega, self.o_valid = res.lam, res.omega, True
                self.busy = False
        if valid:
            self.kes.load(syndromes)
            self.busy = True
        return self.ports()

    def lambda_bus(self) -> int:
        """!
        @brief o_lambda packed: lambda_i in bits [m*i + m - 1 : m*i].
        """
        return GFPoly.from_ascending(self.field, self.o_lambda).to_int(self.t + 1)

    def omega_bus(self) -> int:
        """!
        @brief o_omega packed.
        """
        return GFPoly.from_ascending(self.field, self.o_omega).to_int(self.t)


def simulate(architecture: str, field: GaloisField, t: int,
             stimulus: Iterable[Tuple[bool, bool, Optional[Sequence[int]]]]) -> List[KesPorts]:
    """!
    @brief Run a whole stimulus through a fresh KesTimed.
    @param stimulus One (rst, valid, syndromes) tuple per clock edge.
    @return Outputs after every edge.
    """
    dut = KesTimed(architecture, field, t)
    return [dut.clock(rst, valid, syn) for rst, valid, syn in stimulus]
