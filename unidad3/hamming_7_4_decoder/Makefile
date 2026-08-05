TOPLEVEL_LANG ?= verilog
SIM           ?= icarus

VERILOG_SOURCES := $(PWD)/hamming_decoder.sv \
                   $(PWD)/syndrome_calculator.sv \
                   $(PWD)/error_corrector.sv \
                   $(PWD)/dump.sv
TOPLEVEL        := hamming_decoder
MODULE          := hamming_decoder_tb

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
	gtkwave hamming_decoder.vcd &