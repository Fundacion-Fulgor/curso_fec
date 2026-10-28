# Course container

Ubuntu 24.04 with everything needed for the reference models and the RTL verification:

| Tool | Version | Notes |
|---|---|---|
| SystemC | 3.0.2 | shared library in `/opt/systemc`, C++17 |
| Verilator | 5.052 | in `/opt/verilator`; cocotb 2.x needs 5.036+, SAIF output needs 5.042+ |
| cocotb | 2.1.0 | |
| pyuvm | 5.0.0 | |
| cocotb-coverage | 2.0 | |
| pybind11 | 3.1.0 | CMake config in `/opt/venv/share/cmake` |
| pytest, pytest-cov | 9.1.1, 7.1.0 | |

The image builds natively on x86_64 and on arm64 (Apple Silicon); nothing is emulated.

## Requirements

- macOS: Docker Desktop.
- Windows: Docker Desktop with the WSL 2 backend, and the repository cloned inside the WSL file system (`~/...`, not `/mnt/c/...`), which is much faster.
- Linux: Docker Engine, with your user in the `docker` group.

## Usage

From anywhere inside the repository:

```
docker/run.sh                          # shell in the container
docker/run.sh pytest -q                # a command, run in the current directory
docker/run.sh docker/smoke/smoke_test.sh
docker/run.sh --rebuild                # rebuild the image
```

The repository is mounted at `/work` and the current directory is kept, so `cd unidad7` followed by `../docker/run.sh pytest` runs the unit 7 tests. The container user has the same ids as the host user, so generated files are not owned by root.

The first run builds the image, which compiles SystemC and Verilator: about 10 to 20 minutes on a laptop. With little memory assigned to Docker Desktop, limit the parallel jobs with `FEC_JOBS=2 docker/run.sh --rebuild`.

## Smoke test

`docker/smoke/smoke_test.sh` prints the versions and then checks, in order:

1. A SystemC model compiled as a Python module with pybind11, elaborated once and advanced with several `sc_start()` calls from Python.
2. A small RTL module built with Verilator and driven by a pyuvm test through the cocotb runner, including the SAIF activity file used for power estimation.
3. The `algebra` and `unidad7` test suites.

## Waveforms

Waveforms (`dump.vcd`, `dump.fst`) are written inside the repository and opened on the host with GTKWave or Surfer; the container has no graphical tools.

## Changing versions

Versions are build arguments (`SYSTEMC_TAG`, `VERILATOR_TAG`) in the `Dockerfile` and pinned packages in `requirements.txt`. After changing them, run `docker/run.sh --rebuild`.
