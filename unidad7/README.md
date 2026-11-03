# FEC course: Reed-Solomon key equation solver

Run every command from the root of the repository.

## Requirements

- Docker. `docker/run.sh <command>` runs the command inside the container and builds the image the first time.
- GTKWave, only to view waveforms (outside the container).

```
docker/run.sh                       # open a shell inside the container
docker/run.sh --rebuild <command>   # rebuild the image, then run the command
```

## Run everything

```
docker/run.sh pytest unidad7                       # 298 tests
docker/run.sh pytest algebra/test_galois_unit.py   # 141 tests
```

## Finite field library (algebra/)

```
docker/run.sh pytest algebra/test_galois_unit.py
docker/run.sh python3 algebra/test_galois.py
docker/run.sh pytest algebra/test_galois_unit.py --cov=galois --cov-branch --cov-report=term-missing
```

## Python reference models (unidad7/)

```
docker/run.sh pytest unidad7/tests
```

## C++ and SystemC models

```
docker/run.sh pytest unidad7/key_equation_solver/model/tests
```

## RTL testbench (pyuvm + cocotb)

```
docker/run.sh pytest unidad7/key_equation_solver/tb          # both codes
docker/run.sh pytest unidad7/key_equation_solver/tb -k qr    # RS(44,28) only
docker/run.sh pytest unidad7/key_equation_solver/tb -k tp3   # RS(15,9) only
```

Save waveforms or switching activity (written to `unidad7/key_equation_solver/tb/sim_build/qr/`):

```
docker/run.sh env KES_WAVES=1 pytest unidad7/key_equation_solver/tb -k qr   # dump.vcd
docker/run.sh env KES_SAIF=1 pytest unidad7/key_equation_solver/tb -k qr    # dump.saif
```

## Zhang demo

```
docker/run.sh pytest unidad7/key_equation_solver/demo
docker/run.sh python3 unidad7/key_equation_solver/demo/run_demo.py
unidad7/key_equation_solver/demo/view.sh    # outside the container, needs GTKWave
```

## Useful pytest options

| Option | Effect |
|---|---|
| `-q` | short output |
| `-v` | one line per test |
| `-x` | stop at the first failure |
| `-s` | show what the tests print |
| `-k NAME` | run only the tests whose name contains NAME |