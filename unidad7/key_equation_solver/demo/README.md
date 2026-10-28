# Demo: Example 31 of Zhang on the riBM

Zhang, "VLSI Architectures for Modern Error-Correcting Codes", Example 31 and Table 4.3: an RS(7, 3) code over GF(2^3) with P(x) = x^3 + x + 1, t = 2 and syndromes S = (0, a, a^2, a^6), solved by the riBM by hand. `zhang_demo.v` instantiates `key_equation_solver` with `M = 3, POLY = 'h3, T = 2`, so the registers of the RTL show the values of the book.

## Run

From `unidad7`:

```
../docker/run.sh python3 key_equation_solver/demo/run_demo.py
```

It prints every iteration next to Table 4.3, the output, the scale factor of Lambda and the error values given by the Forney formula with the convention of the book (first root alpha^1), and writes `sim_build/dump.vcd`.

To see the waveform, outside the container (GTKWave is a graphical tool):

```
key_equation_solver/demo/view.sh
```

`view.sh` runs `gtkwave sim_build/dump.vcd zhang_demo.gtkw` from this directory, which is needed because GTKWave looks for the filter `alpha_gf8.txt` in the current directory. With the GTKWave of Windows, open a terminal in this folder and run the same `gtkwave` command.

## What the view shows

Signals are grouped as in the paper: control (r, k, gamma, MC; Fig. 2), rDC (delta-hat and theta-hat; Fig. 4) and ELU (lambda and b; Fig. 3), plus the outputs. Every coefficient has its own signal and is shown as a power of alpha.

| Time | State | What to look at |
|---|---|---|
| 10 ns | load | delta-hat and theta-hat take S(x) = a x + a^2 x^2 + a^6 x^3; Lambda = B = 1 |
| 20 ns | r = 1 | the discrepancy delta-hat_0 = a and k = 1 >= 0, so MC is high |
| 30 ns | r = 2 | MC acted: B took Lambda, theta-hat took delta-hat / x, gamma = a, k = -2 |
| 40 ns | r = 3 | k = -1 < 0: no MC, B is shifted, theta-hat holds |
| 50 ns | output | o_valid; Lambda = a^2 + a^3 x + a^5 x^2, Omega_h = a^6 + a^4 x |

MC is combinational: after the word it may be high again, but it only acts while busy.

`test_run_demo.py` runs the demo with the rest of the unit 7 tests, so a change in the RTL that breaks the example is noticed.
