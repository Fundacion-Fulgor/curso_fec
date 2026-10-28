# Reference models of key_equation_solver

Two Python modules built with CMake and pybind11 from header-only C++:

| Module | Sources | What it is |
|---|---|---|
| `kesmodel` | `gf.hpp`, `kes_core.hpp` | the three solvers at register level, step by step (C++ version of `kes_ref.py`) |
| `kessc` | `kes_sc.hpp` + the above | SystemC modules with the RTL ports and the protocol of `docs/kes_protocol.md` |

`kessc` is only built when CMake finds SystemC (`SystemCLanguage`), as in the course container.

## Build and test

From `unidad7`, `pytest` builds both modules into `model/build` (only changed sources are compiled
again) and runs their tests with the rest of the unit:

```
../docker/run.sh pytest
```

`tests/test_model.py` compares `kesmodel` with `galois.py` and with every state of `kes_ref.py`.
`tests/test_systemc.py` compares `kessc` with `kes_timing.py` on every clock edge.

## kesmodel

```python
import kesmodel
F = kesmodel.GaloisField(8, 0x1D)
kes = kesmodel.make_kes("ribm_systolic", F, 8)
res = kes.run(syndromes, trace=True)    # res.lam, res.omega, res.form, res.length, res.trace
kes.load(syndromes); kes.step(); kes.snapshot()    # same dict keys as kes_ref
```

## kessc

SystemC elaborates once per process, so one bench holds every instance on a single clock
(10 ns, first rising edge at 5 ns): codes `"qr"` (m = 8, t = 8) and `"tp3"` (m = 4, t = 3), each
with `"ibm"`, `"ribm"` and `"ribm_systolic"`.

```python
import kessc
h = kessc.harness()
p = h.clock("qr", "ribm_systolic", valid=True, syndromes=syn)    # one rising edge
p = h.clock("qr", "ribm_systolic")    # p: o_valid, o_lambda, o_omega, lambda_bus, omega_bus, busy, state
h.snapshot("qr", "ribm_systolic")     # internal registers, as in kes_ref
```

Each `clock()` drives one instance and gives an idle edge to the others. The module defines and
exports an empty `sc_main`: `libsystemc` refers to it, and pybind11 builds with hidden visibility.
