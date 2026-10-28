"""!
@file kes_ref.py
@brief Register-level reference models of the key equation solvers of
       Sarwate and Shanbhag, "High-speed architectures for Reed-Solomon
       decoders", IEEE TVLSI 9(5), 2001.

  IBM           inversionless BM, Blahut version (section III-B, Figs. 1-3), 3t steps
  RiBMSplit     riBM, rDC block + ELU block (section IV-B, Fig. 5), 2t steps
  RiBMSystolic  RiBM, one array of 3t + 1 PE1 (section IV-B3, Fig. 6), 2t steps

Inputs are the 2t syndromes S_0 .. S_(2t-1). Outputs are lists in ascending
order: t + 1 coefficients of the error locator and t of the evaluator. All the
solvers return Lambda scaled by a non-zero constant (lambda_0 is the product
of the gammas, not 1).

Each step() is one clock edge. State r is what the registers hold after r
steps from load(), so state 0 is the loaded state and the outputs are read
from the last one.

With up to t errors the three architectures give the same Lambda. With more,
Lambda may need more than t + 1 coefficients and each architecture truncates
it differently, so their outputs can differ. In the RiBM, Lambda entering from
the right can also overwrite delta, since the bound that prevents it in the
paper assumes at most t errors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

from galois import GaloisField, GFPoly

__all__ = ["KesResult", "KesBase", "IBM", "RiBMSplit", "RiBMSystolic",
           "ARCHITECTURES", "make_kes"]


@dataclass
class KesResult:
    """!
    @brief Output of a key equation solver.

    form is "classic" when omega holds Omega(x) = Lambda(x) S(x) mod x^t and
    "high" when it holds the coefficients of x^(2t) .. x^(3t-1) of
    Lambda(x) S(x), as the riBM and RiBM deliver it. length is the LFSR
    length L = (2t - k) / 2 kept by the control unit; after a successful
    decoding it equals the degree of Lambda.
    """

    field: GaloisField
    lam: List[int]
    omega: List[int]
    form: str
    iterations: int
    length: int
    trace: Optional[List[dict]] = None

    @property
    def degree(self) -> int:
        """!
        @brief Degree of the error locator (-1 if it is zero).
        """
        return max((i for i, c in enumerate(self.lam) if c), default=-1)

    def lambda_bus(self) -> int:
        """!
        @brief Error locator packed as o_lambda: lambda_i in bits [m*i + m - 1 : m*i].
        """
        return GFPoly.from_ascending(self.field, self.lam).to_int(len(self.lam))

    def omega_bus(self) -> int:
        """!
        @brief Evaluator packed as o_omega.
        """
        return GFPoly.from_ascending(self.field, self.omega).to_int(len(self.omega))


class KesBase:
    """!
    @brief State handling shared by the three architectures.

    Subclasses define the registers (_clear, _init, _registers), one clock
    edge (_step), the discrepancy seen by the control unit and the outputs.
    """

    name = "base"
    form = "classic"

    def __init__(self, field: GaloisField, t: int) -> None:
        """!
        @brief Build a solver for a t-error-correcting code.
        @param field Symbol field.
        @param t Error correcting capability; the solver takes 2t syndromes.
        @throws TypeError If field is not a GaloisField.
        @throws ValueError If t is not an int >= 1.
        """
        if not isinstance(field, GaloisField):
            raise TypeError(f"field must be a GaloisField, got {type(field).__name__}")
        if not isinstance(t, int) or isinstance(t, bool) or t < 1:
            raise ValueError(f"t must be an int >= 1, got {t!r}")
        self.field = field
        self.t = t
        self.reset()

    @property
    def iterations(self) -> int:
        """!
        @brief Clock edges from load() to the result.
        """
        return 2 * self.t

    @property
    def done(self) -> bool:
        return self.loaded and self.r == self.iterations

    @property
    def mc(self) -> bool:
        """!
        @brief Control signal MC(r): discrepancy != 0 and k(r) >= 0, only
               during the 2t Berlekamp-Massey steps.
        """
        return self.loaded and self.r < 2 * self.t and self.discrepancy != 0 and self.k >= 0

    def reset(self) -> None:
        """!
        @brief Clear every register, as a synchronous reset does.
        """
        self.r = 0
        self.loaded = False
        self.gamma = 0
        self.k = 0
        self._clear()

    def load(self, syndromes: Sequence[int]) -> None:
        """!
        @brief Load the syndromes and set the initial state (state 0).
        @param syndromes S_0 .. S_(2t-1).
        @throws ValueError If there are not 2t syndromes or one is not a field element.
        """
        syn = list(syndromes)
        if len(syn) != 2 * self.t:
            raise ValueError(f"expected {2 * self.t} syndromes, got {len(syn)}")
        for s in syn:
            self.field.check(s)
        self.r = 0
        self.loaded = True
        self.gamma = 1
        self.k = 0
        self._init(syn)

    def step(self) -> None:
        """!
        @brief One clock edge.
        @throws RuntimeError If nothing is loaded or every step is already done.
        """
        if not self.loaded:
            raise RuntimeError("load() must be called before step()")
        if self.r >= self.iterations:
            raise RuntimeError(f"all {self.iterations} steps are done")
        self._step()
        self.r += 1

    def result(self, trace: Optional[List[dict]] = None) -> KesResult:
        """!
        @brief Outputs after the last step.
        @throws RuntimeError If the computation is not finished.
        """
        if not self.done:
            raise RuntimeError("the computation is not finished")
        lam, omega = self._outputs()
        length = (2 * self.t - self.k) // 2
        return KesResult(self.field, list(lam), list(omega), self.form, self.iterations,
                         length, trace)

    def run(self, syndromes: Sequence[int], trace: bool = False) -> KesResult:
        """!
        @brief load() followed by every step.
        @param syndromes S_0 .. S_(2t-1).
        @param trace If True, keep a snapshot of every state from 0 to the last.
        """
        self.load(syndromes)
        states = [self.snapshot()] if trace else None
        while self.r < self.iterations:
            self.step()
            if states is not None:
                states.append(self.snapshot())
        return self.result(states)

    def snapshot(self) -> Dict[str, object]:
        """!
        @brief Copy of the registers and of the control signals of the current state.
        """
        state: Dict[str, object] = {"r": self.r, "gamma": self.gamma, "k": self.k,
                                    "mc": self.mc, "discrepancy": self.discrepancy}
        state.update({name: list(value) for name, value in self._registers().items()})
        return state

    def _control(self, discrepancy: int) -> bool:
        """!
        @brief Control unit update of gamma and k.
        @return MC for this step.
        """
        mc = discrepancy != 0 and self.k >= 0
        if mc:
            self.gamma = discrepancy
            self.k = -self.k - 1
        else:
            self.k += 1
        return mc

    def _pe_update(self, upper: List[int], lower: List[int], d0: int) -> List[int]:
        """!
        @brief Datapath shared by PE0 and PE1: gamma * upper + d0 * lower.
        """
        mul, g = self.field.mul, self.gamma
        return [mul(g, a) ^ mul(d0, b) for a, b in zip(upper, lower)]


class IBM(KesBase):
    """!
    @brief Inversionless BM, Blahut version (section III-B).

    The DC block holds the syndromes in a circular register DS and multiplies
    t + 1 of its latches by lambda_0 .. lambda_t; the ELU block has t + 1 PE0.
    Steps 0 .. 2t-1 compute Lambda. At the last of them the latches
    DS_t .. DS_(2t-1) are cleared, and steps 2t .. 3t-1 reuse the DC block with
    Lambda frozen to produce omega_0 .. omega_(t-1) (step iBM.4).
    """

    name = "ibm"
    form = "classic"

    @property
    def iterations(self) -> int:
        return 3 * self.t

    @property
    def discrepancy(self) -> int:
        """!
        @brief DC block output: sum of DS_((2t - i) mod 2t) * lambda_i.
        """
        mul, n = self.field.mul, 2 * self.t
        acc = 0
        for i, lam in enumerate(self.lam):
            acc ^= mul(self.ds[-i % n], lam)
        return acc

    def _clear(self) -> None:
        t = self.t
        self.lam = [0] * (t + 1)
        self.b = [0] * (t + 1)
        self.ds = [0] * (2 * t)
        self.omega: List[int] = []

    def _init(self, syn: List[int]) -> None:
        t = self.t
        self.lam = [1] + [0] * t
        self.b = [1] + [0] * t
        self.ds = syn
        self.omega = []

    def _step(self) -> None:
        n = 2 * self.t
        delta = self.discrepancy
        rotated = self.ds[1:] + self.ds[:1]
        if self.r < n:
            lam = self._pe_update(self.lam, [0] + self.b[:-1], delta)
            self.b = self.lam if self._control(delta) else [0] + self.b[:-1]
            self.lam = lam
            if self.r == n - 1:
                rotated[self.t:] = [0] * self.t
        else:
            self.omega.append(delta)
        self.ds = rotated

    def _registers(self) -> Dict[str, List[int]]:
        return {"lam": self.lam, "b": self.b, "ds": self.ds, "omega": self.omega}

    def _outputs(self):
        return self.lam, self.omega


class RiBMSplit(KesBase):
    """!
    @brief riBM (section IV-B, Fig. 5): an rDC block of 2t PE1 holding
           delta-hat and theta-hat, and an ELU block of t + 1 PE0 holding
           Lambda and B. omega is the high-order evaluator.
    """

    name = "ribm"
    form = "high"

    @property
    def discrepancy(self) -> int:
        return self.delta[0]

    def _clear(self) -> None:
        t = self.t
        self.lam = [0] * (t + 1)
        self.b = [0] * (t + 1)
        self.delta = [0] * (2 * t)
        self.theta = [0] * (2 * t)

    def _init(self, syn: List[int]) -> None:
        t = self.t
        self.lam = [1] + [0] * t
        self.b = [1] + [0] * t
        self.delta = syn
        self.theta = list(syn)

    def _step(self) -> None:
        d0 = self.delta[0]
        shifted = self.delta[1:] + [0]
        lam = self._pe_update(self.lam, [0] + self.b[:-1], d0)
        delta = self._pe_update(shifted, self.theta, d0)
        if self._control(d0):
            self.b, self.theta = self.lam, shifted
        else:
            self.b = [0] + self.b[:-1]
        self.lam, self.delta = lam, delta

    def _registers(self) -> Dict[str, List[int]]:
        return {"lam": self.lam, "b": self.b, "delta": self.delta, "theta": self.theta}

    def _outputs(self):
        return self.lam, self.delta[: self.t]


class RiBMSystolic(KesBase):
    """!
    @brief RiBM (section IV-B3, Fig. 6): a single array of 3t + 1 PE1.

    PE1_0 .. PE1_(2t-1) start with the syndromes, PE1_2t .. PE1_(3t-1) with
    zeros and PE1_3t with 1. Lambda and B enter from the right and end in
    PE1_t .. PE1_2t; omega ends in PE1_0 .. PE1_(t-1).
    """

    name = "ribm_systolic"
    form = "high"

    @property
    def discrepancy(self) -> int:
        return self.delta[0]

    def _clear(self) -> None:
        size = 3 * self.t + 1
        self.delta = [0] * size
        self.theta = [0] * size

    def _init(self, syn: List[int]) -> None:
        self.delta = syn + [0] * self.t + [1]
        self.theta = list(self.delta)

    def _step(self) -> None:
        d0 = self.delta[0]
        shifted = self.delta[1:] + [0]
        self.delta = self._pe_update(shifted, self.theta, d0)
        if self._control(d0):
            self.theta = shifted

    def _registers(self) -> Dict[str, List[int]]:
        return {"delta": self.delta, "theta": self.theta}

    def _outputs(self):
        t = self.t
        return self.delta[t: 2 * t + 1], self.delta[:t]


ARCHITECTURES = {cls.name: cls for cls in (IBM, RiBMSplit, RiBMSystolic)}


def make_kes(name: str, field: GaloisField, t: int) -> KesBase:
    """!
    @brief Build a solver by name.
    @param name "ibm", "ribm" or "ribm_systolic".
    @param field Symbol field.
    @param t Error correcting capability.
    @throws ValueError If the name is unknown.
    """
    if name not in ARCHITECTURES:
        raise ValueError(f"unknown architecture {name!r}, expected one of {sorted(ARCHITECTURES)}")
    return ARCHITECTURES[name](field, t)
