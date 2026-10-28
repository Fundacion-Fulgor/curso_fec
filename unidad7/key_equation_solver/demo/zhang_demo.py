"""!
@file zhang_demo.py
@brief Example 31 of Zhang, "VLSI Architectures for Modern Error-Correcting Codes", on the RTL.

RS(7, 3) over GF(2^3) with P(x) = x^3 + x + 1, t = 2 and syndromes
S = (0, a, a^2, a^6). The test loads the syndromes, reads the registers of the
riBM after every clock edge, prints them next to Table 4.3 of the book and
fails if any value differs. Then it applies the Forney formula with the
convention of the book (first root alpha^1).
"""

import cocotb
from cocotb.clock import Clock
from cocotb.triggers import ReadOnly, RisingEdge

from decode_ref import chien_forney
from galois import GaloisField

F = GaloisField(3, 0b011)
M, T = 3, 2

# alpha^0 .. alpha^6 = 1, 2, 4, 3, 6, 7, 5
A1, A2, A3, A4, A5, A6 = (F.exp(i) for i in range(1, 7))
SYNDROMES = [0, A1, A2, A6]

# Table 4.3, ascending coefficients: (Lambda, B, k, gamma, delta-hat, theta-hat) at state r
BOOK = [
    ([1, 0, 0], [1, 0, 0], 0, 1, [0, A1, A2, A6], [0, A1, A2, A6]),
    ([1, 0, 0], [0, 1, 0], 1, 1, [A1, A2, A6, 0], [0, A1, A2, A6]),
    ([1, 0, A1], [1, 0, 0], -2, A1, [A2, 1, A3, 1], [A2, A6, 0, 0]),
    ([A1, A2, A2], [0, 1, 0], -1, A1, [A2, A2, A1, 0], [A2, A6, 0, 0]),
]
BOOK_LAMBDA = [A2, A3, A5]
BOOK_OMEGA = [A6, A4]


def unpack(value, count):
    return [(value >> (M * i)) & ((1 << M) - 1) for i in range(count)]


SUPERSCRIPT = str.maketrans("0123456789-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b")


def power(base, exponent):
    return base if exponent == 1 else base + str(exponent).translate(SUPERSCRIPT)


def poly(coeffs):
    """!
    @brief Ascending coefficients written as in the book: a x + a^2 x^2 + a^6 x^3.
    """
    terms = []
    for degree, c in enumerate(coeffs):
        if c == 0:
            continue
        coeff = "" if c == 1 and degree else ("1" if c == 1 else power("\u03b1", F.log(c)))
        var = "" if degree == 0 else power("x", degree)
        terms.append(coeff + var)
    return " + ".join(terms) if terms else "0"


def mark(ok):
    return "ok" if ok else "DIFFERS"


@cocotb.test()
async def zhang_example_31(dut):
    """!
    @brief Run the example and compare every state with the book.
    """
    Clock(dut.i_clk, 10, unit="ns").start()
    dut.i_rst.value = 1
    dut.i_valid.value = 0
    dut.i_syn.value = 0
    await RisingEdge(dut.i_clk)
    dut.i_rst.value = 0
    dut.i_valid.value = 1
    dut.i_syn.value = sum(s << (M * j) for j, s in enumerate(SYNDROMES))
    await RisingEdge(dut.i_clk)
    dut.i_valid.value = 0

    kes = dut.u_kes
    differences = 0
    print()
    print(f"Syndromes S(x) = {poly(SYNDROMES)}")
    print(f"{'r':>2} | {'Lambda':<16} | {'B':<4} | {'k':>2} | {'gamma':<5} | "
          f"{'delta-hat':<18} | {'theta-hat':<18} | MC    | book")
    for r, (lam, b, k, gamma, delta, theta) in enumerate(BOOK):
        await ReadOnly()
        rtl = (unpack(int(kes.lam_q.value), T + 1), unpack(int(kes.b_q.value), T + 1),
               kes.u_ctrl.k_q.value.to_signed(), int(kes.u_ctrl.gamma_q.value),
               unpack(int(kes.delta_q.value), 2 * T), unpack(int(kes.theta_q.value), 2 * T))
        ok = int(kes.u_ctrl.count_q.value) == r and rtl == (lam, b, k, gamma, delta, theta)
        differences += not ok
        print(f"{r:>2} | {poly(rtl[0]):<16} | {poly(rtl[1]):<4} | {rtl[2]:>2} | "
              f"{poly([rtl[3]]):<5} | {poly(rtl[4]):<18} | {poly(rtl[5]):<18} | "
              f"{str(bool(int(kes.mc.value))):<5} | {mark(ok)}")
        await RisingEdge(dut.i_clk)

    await ReadOnly()
    o_valid = int(dut.o_valid.value)
    lam = unpack(int(dut.o_lambda.value), T + 1)
    omega = unpack(int(dut.o_omega.value), T)
    for _ in range(3):    # a few more clocks, so the outputs show in the waveform
        await RisingEdge(dut.i_clk)
    assert o_valid == 1, "o_valid must rise after 2t = 4 iterations"
    differences += lam != BOOK_LAMBDA or omega != BOOK_OMEGA
    print(f"output: Lambda(x) = {poly(lam)}, Omega_h(x) = {poly(omega)}  "
          f"{mark(lam == BOOK_LAMBDA and omega == BOOK_OMEGA)}")

    scale = F.inverse(lam[0])
    print(f"lambda_0 = {poly([lam[0]])} = product of the gammas; "
          f"normalized Lambda(x) = {poly([F.mul(scale, c) for c in lam])}")

    # Forney with the first root of the book: S_j = r(alpha^(j+1)), so b = 1
    cf = chien_forney(F, 7, lam, omega, "high", first_root=1)
    errors = ", ".join(f"position {p}: {poly([v])}" for p, v in zip(cf.positions, cf.values))
    print(f"Chien + Forney (first root alpha^1): {errors}")
    assert (cf.positions, cf.values, cf.fail) == ([0, 3], [1, A4], False)

    assert differences == 0, f"{differences} states differ from Table 4.3"
