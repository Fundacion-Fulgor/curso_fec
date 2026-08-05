TOPLEVEL_LANG ?= verilog
SIM           ?= icarus

VERILOG_SOURCES := $(PWD)/hamming_7_4.sv \
                   $(PWD)/../hamming_7_4_encoder/hamming_encoder.sv \
                   $(PWD)/../hamming_7_4_decoder/hamming_decoder.sv \
                   $(PWD)/../hamming_7_4_decoder/syndrome_calculator.sv \
                   $(PWD)/../hamming_7_4_decoder/error_corrector.sv \
                   $(PWD)/dump.sv
TOPLEVEL        := hamming_7_4
MODULE          := hamming_7_4_tb

# para que el testbench encuentre galois.py
export PYTHONPATH := $(PWD):$(PWD)/../../algebra:$(PYTHONPATH)

COCOTB_HDL_TIMEUNIT      = 1ns
COCOTB_HDL_TIMEPRECISION = 1ps

ifeq ($(SIM),icarus)
    COMPILE_ARGS += -g2012 -s dump
endif

include $(shell cocotb-config --makefiles)/Makefile.sim

.PHONY: ondas
ondas: sim
	gtkwave hamming_7_4.vcd &