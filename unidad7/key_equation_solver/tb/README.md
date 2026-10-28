# Testbench of key_equation_solver

pyuvm on Verilator through the cocotb runner. From `unidad7`:

```
../docker/run.sh pytest key_equation_solver/tb            # both configurations
../docker/run.sh pytest key_equation_solver/tb -k qr      # only the TP3 module
KES_WAVES=1 ...    dump.vcd   in tb/sim_build/<code>
KES_SAIF=1  ...    dump.saif  in tb/sim_build/<code>, for power estimation
KES_SEED=7  ...    another seed for the random sequence
```

`test_kes_rtl.py` builds the RTL once per configuration (`qr`: M = 8, T = 8; `tp3`: M = 4,
T = 3) and runs the three tests of `kes_tb.py`:

| Test | Stimulus |
|---|---|
| `KesReferenceTest` | known words: the exercise 6 vector with two errors (`qr`) or r1 and r2 of exercise 3 (`tp3`) |
| `KesSuiteTest` | every category of `stimuli.py`, back to back |
| `KesRandomTest` | 300 random words: up to t errors, beyond t, random syndromes, gaps, early loads and resets |

Components:

- `KesDriver` drives one item per call: idle edges, then a load or a reset.
- `KesMonitor` samples the inputs on each rising edge and then the outputs and the registers `lam_q`, `b_q`, `delta_q`, `theta_q`, `gamma_q`, `k_q`, `count_q`, `busy_q` and `mc`.
- `KesScoreboard` runs the SystemC model `kessc` in lockstep and compares outputs and registers on every edge, checks that the datapath holds while idle, and decodes every result with `decode_ref.chien_forney` against the injected errors.
- `KesCoverage` counts the degree of every result, back-to-back words, early loads, resets while busy and MC decisions; each test lists the bins it requires.
