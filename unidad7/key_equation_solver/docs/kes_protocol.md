# key_equation_solver: interface, timing and implementation

`rtl/key_equation_solver.v` implements the reformulated inversionless Berlekamp-Massey architecture (riBM) of Sarwate and Shanbhag, "High-speed architectures for Reed-Solomon decoders", IEEE TVLSI 9(5), 2001. The timing below is also implemented by `unidad7/kes_timing.py` and by the SystemC model `kessc` (`model/`), and the three are checked against each other clock by clock.

## Ports and parameters

| Port | Width (default) | Description |
|---|---|---|
| `i_clk` | 1 | clock, rising edge |
| `i_rst` | 1 | synchronous reset, active high |
| `i_valid` | 1 | `i_syn` holds a new set of syndromes |
| `i_syn` | 2TM = 128 | S_j in bits `[M*j+M-1 : M*j]`, j = 0 .. 2T-1, with S_j = r(alpha^j) |
| `o_valid` | 1 | one-cycle pulse: `o_lambda` and `o_omega` hold a new result |
| `o_lambda` | (T+1)M = 72 | lambda_i in bits `[M*i+M-1 : M*i]`, i = 0 .. T |
| `o_omega` | TM = 64 | omega_i in bits `[M*i+M-1 : M*i]`, i = 0 .. T-1 |

Parameters `M = 8`, `POLY = 'h1D` (P(x) without the x^M term) and `T = 8` give the module of TP3. `M = 4`, `POLY = 'h3`, `T = 3` give the reduced code of part A.

Lambda is returned scaled by a non-zero constant: lambda_0 is the product of the gammas, not 1. `o_omega` holds the high-order evaluator Omega_h: the coefficients of x^2T .. x^(3T-1) of Lambda(x) S(x).

## Behaviour on each rising edge

In priority order:

1. If `i_rst` is high, the solver becomes idle and `o_valid`, `o_lambda` and `o_omega` go to 0. Nothing else happens on that edge, even with `i_valid` high.
2. `o_valid` goes to 0, unless rule 3 sets it.
3. If a word is in progress, the solver performs one iteration. After the last one (iteration 2T, counting from 1) the output registers load the result, `o_valid` goes to 1 and the solver becomes idle.
4. If `i_valid` is high, the solver loads `i_syn` (state 0) and starts a new word. A word still in progress is dropped without producing an output.

`o_lambda` and `o_omega` keep their value until the next result or reset. Loading on the same edge as the last iteration of the previous word gives back-to-back operation: latency and interval are both 2T = 16 clocks, as the TP3 decoder needs.

```
edge       E0   E1   E2  ...  E15  E16  E17  ...  E31  E32
i_valid     1    0    0  ...    0    1    0  ...    0    1
state       0    1    2  ...   15    0    1  ...   15    0
o_valid     0    0    0  ...    0    1    0  ...    0    1
o_lambda    .    .    .  ...    .    A    A  ...    A    B
```

## Forney formula

With $S_j = r(\alpha^j)$ and an error at $X = \alpha^p$, equation (12) of the paper with $m_0 = 0$ gives

$$
Y = X^{-(2T-1)} * \Omega_h(X^{-1}) / \Lambda'(X^{-1})
$$

and Lambda'(x) only keeps the odd-degree terms of Lambda. The scale factor of Lambda and Omega_h
cancels out.

## Implementation and the paper

| RTL | Paper |
|---|---|
| `key_equation_solver` | riBM architecture, Fig. 5: rDC block + ELU block |
| `g_rdc[i].u_pe1` (`kes_pe1`), i = 0 .. 2T-1 | PE1 processor, Fig. 4(a); rDC block, Fig. 4(b) |
| `g_elu[i].u_pe0` (`kes_pe0`), i = 0 .. T | PE0 processor, Fig. 3(a); ELU block, Fig. 3(b) |
| `u_ctrl` (`kes_ctrl`) | control unit, Fig. 2 (counter, k, gamma, MC) |
| `kes_gf_mul` | the m-bit parallel multipliers assumed in section III-B |
| `delta_q`, `theta_q` | delta-hat_i(r), theta-hat_i(r) of the riBM pseudocode |
| `lam_q`, `b_q` | lambda_i(r), b_i(r) |
| `u_ctrl.gamma_q`, `u_ctrl.k_q` | gamma(r), k(r) |
| `u_ctrl.count_q` | r |
| `delta_d`, `lam_d` (and the PE outputs) | step riBM.1 |
| `theta_d`, `b_d` multiplexers, `kes_ctrl` update | step riBM.2 |

The PE modules are combinational; their latches are the registers above, kept in the top module as packed vectors so the testbench can read them in one access. Each clock edge while busy is one iteration r of the pseudocode, and the discrepancy is always `delta_q[M-1:0]` (delta-hat_0).

What the paper leaves open and is fixed here:

- Loading: `delta_q` and `theta_q` take the syndromes, `lam_q` and `b_q` take 1, gamma takes 1 and k takes 0 (the initialization of the pseudocode).
- k is a 6-bit two's complement register: k >= 0 is its sign bit, and -k-1 is `~k`.
- Only the control (`busy_q`, `count_q`) and the outputs are reset. The datapath is loaded with every word, so it has no reset.
- The datapath holds its value while idle, so it does not toggle without words. The protocol does not require it; the testbench checks it as a design intent.
- The top coefficient of `b_q` is never read (b_T is shifted out); synthesis removes it. It is kept so that the ELU matches Fig. 3.

Resources for T = 8, as in Table I of the paper: 6T+2 = 50 GF multipliers, 3T+1 = 25 byte XORs, 3T+1 = 25 multiplexers and 6T+2 = 50 byte registers, plus 19 control flip-flops and 137 output flip-flops. Critical path: one multiplier, one XOR and the load multiplexer (eq. (13)).

## Internal state

While a word is in progress, the registers after the r-th edge following the load (r = 0 .. 2T-1) hold state r of the riBM model, `kes_ref.RiBMSplit(field, T).run(syn, trace=True).trace[r]`: `lam`, `b`, `delta`, `theta`, `gamma`, `k`, and `mc` for the combinational `mc` wire. State 2T only exists in the output registers.

## Beyond t errors

With more than t errors, Lambda may need more than t + 1 coefficients and is truncated, and the Chien search usually flags the failure. The models in `kes_ref` also include the iBM and the RiBM, useful in class to compare the three architectures; with more than t errors they can return different results from each other, which is expected.
