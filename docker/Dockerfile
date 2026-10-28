# Development environment for the FEC course: SystemC + pybind11 reference models,
# Verilator + cocotb + pyuvm for RTL verification.
# Builds natively on x86_64 (Linux, WSL) and arm64 (Apple Silicon).

FROM ubuntu:24.04

COPY install/ /tmp/install/
RUN bash /tmp/install/apt_packages.sh

# Slow builds first, so changing the Python packages does not rebuild them.
ARG JOBS=4
ARG SYSTEMC_TAG=3.0.2
RUN bash /tmp/install/systemc.sh "$SYSTEMC_TAG" /opt/systemc "$JOBS"

ARG VERILATOR_TAG=v5.052
RUN bash /tmp/install/verilator.sh "$VERILATOR_TAG" /opt/verilator "$JOBS"

COPY requirements.txt /tmp/requirements.txt
RUN python3 -m venv /opt/venv \
 && /opt/venv/bin/pip install --no-cache-dir --upgrade pip \
 && /opt/venv/bin/pip install --no-cache-dir -r /tmp/requirements.txt \
 && echo /opt/systemc/lib > /etc/ld.so.conf.d/systemc.conf \
 && ldconfig \
 && rm -rf /tmp/install /tmp/requirements.txt

ENV SYSTEMC_HOME=/opt/systemc \
    CMAKE_PREFIX_PATH=/opt/systemc:/opt/venv \
    VIRTUAL_ENV=/opt/venv \
    PATH=/opt/venv/bin:/opt/verilator/bin:$PATH \
    PYTHONDONTWRITEBYTECODE=1

# Regular user with the host ids, so files written to the mounted repository keep their owner.
ARG USER_UID=1000
ARG USER_GID=1000
RUN userdel --remove ubuntu 2>/dev/null || true \
 && groupadd --non-unique --gid "$USER_GID" dev \
 && useradd --non-unique --uid "$USER_UID" --gid "$USER_GID" --create-home --shell /bin/bash dev

USER dev
WORKDIR /work
CMD ["bash"]
